#!/usr/bin/env python3
"""
bundle_pipeline_results.py

Collects all key outputs needed to (1) evaluate scientific validity of the pipeline run and
(2) assess biological interpretation derived from variant calling.

It creates an "analysis report bundle" folder containing:
  - QC (fastp + FastQC)
  - Alignment/Mapping BAMs + indices + (optional) stats files
  - PreProcess (MarkDuplicates metrics, BQSR tables, BQSR BAMs)
  - Variants (VCF/VCF.GZ + indices + any filtered/annotated outputs)
  - Coverage (mosdepth/bedtools/collect* metrics if present)
  - Logs (pipeline logs, stdout/stderr logs if present)
  - manifest.json (what was found/copied, what is missing)

Works with your repo structure like:
  runs/HCC1395_Germline/Bwa/QC
  runs/HCC1395_Germline/Bwa/Mapping
  runs/HCC1395_Germline/Bwa/PreProcess
  runs/HCC1395_Tumor/Bwa/...

Usage examples:
  python bundle_pipeline_results.py --repo-root . --basedir runs --sample HCC1395 --maptype Bwa
  python bundle_pipeline_results.py --repo-root /path/to/repo --basedir runs --sample HCC1395 --maptype Bwa --outdir analysis_report_bundle

Optional:
  --tumor-name Tumor --germline-name Germline (defaults)
  --copy (default is symlinks; use --copy to physically copy files)

Note:
  This script does NOT run any pipeline steps; it only bundles outputs.
"""

import argparse
import json
import os
import re
import shutil
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple

DEFAULT_TUMOR_NAME = "Tumor"
DEFAULT_GERMLINE_NAME = "Germline"

# What we *expect* to find (patterns). We will bundle whatever exists.
QC_PATTERNS = [
    "**/*fastqc.html",
    "**/*fastqc.zip",
    "**/*.json",          # fastp json often
    "**/*.html",          # fastp html often (QC dir)
]

MAPPING_PATTERNS = [
    "**/*.bam",
    "**/*.bai",
    "**/*.crai",
    "**/*.sam",
    "**/*.cram",
    "**/*.txt",           # flagstat/stats/metrics might be txt
    "**/*.log",
]

PREPROCESS_PATTERNS = [
    "**/*.table",         # GATK BQSR tables
    "**/*RECAL*",
    "**/*BQSR*.bam",
    "**/*BQSR*.bai",
    "**/*metrics*.txt",
    "**/*marked_dup*.txt",
    "**/*.log",
    "**/*.txt",
]

VARIANT_PATTERNS = [
    "**/*.vcf",
    "**/*.vcf.gz",
    "**/*.tbi",
    "**/*.idx",
    "**/*.stats",
    "**/*.log",
    "**/*.txt",
    "**/*.tsv",
    "**/*.csv",
]

COVERAGE_PATTERNS = [
    "**/*mosdepth*",
    "**/*coverage*",
    "**/*depth*",
    "**/*Collect*Metrics*",
    "**/*.bed",
    "**/*.tsv",
    "**/*.txt",
]

LOG_PATTERNS = [
    "**/*.log",
    "**/*command*",
    "**/*stderr*",
    "**/*stdout*",
    "**/*.out",
    "**/*.err",
]


def safe_mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def relpath_or_name(src: Path, base: Path) -> str:
    try:
        return str(src.relative_to(base))
    except Exception:
        return src.name


def find_files(root: Path, patterns: List[str]) -> List[Path]:
    files: List[Path] = []
    for pat in patterns:
        files.extend(root.glob(pat))
    # unique + only files
    uniq = []
    seen = set()
    for f in files:
        if f.is_file():
            key = str(f.resolve())
            if key not in seen:
                uniq.append(f)
                seen.add(key)
    return sorted(uniq)


def copy_or_link(src: Path, dst: Path, mode: str) -> None:
    """
    mode: 'symlink' (default) or 'copy'
    """
    safe_mkdir(dst.parent)
    if dst.exists() or dst.is_symlink():
        return
    if mode == "copy":
        shutil.copy2(src, dst)
    else:
        # symlink relative when possible (more portable inside repo bundle)
        try:
            rel = os.path.relpath(src.resolve(), dst.parent.resolve())
            dst.symlink_to(rel)
        except Exception:
            dst.symlink_to(src.resolve())


def stage_bundle(
    label: str,
    search_root: Path,
    patterns: List[str],
    out_root: Path,
    repo_root: Path,
    mode: str,
    manifest: Dict,
) -> None:
    found = find_files(search_root, patterns)
    dest_dir = out_root / label
    safe_mkdir(dest_dir)

    copied = []
    for src in found:
        # preserve some directory structure under the stage folder
        rel = relpath_or_name(src, search_root)
        # if rel includes parent dirs, keep them; otherwise just name
        dst = dest_dir / rel
        copy_or_link(src, dst, mode)
        copied.append({
            "src": str(src),
            "dst": str(dst),
            "rel_from_repo": relpath_or_name(src, repo_root),
        })

    manifest["stages"][label] = {
        "search_root": str(search_root),
        "patterns": patterns,
        "found_count": len(found),
        "bundled_count": len(copied),
        "files": copied,
    }


def guess_variant_dirs(sample_dir: Path, maptype: str) -> List[Path]:
    """
    Try to locate variant caller output dirs under:
      runs/<sample_type>/<maptype>/<CallerName>/*
    We search for common callers folders; if none, we just return maptype folder for pattern match.
    """
    map_dir = sample_dir / maptype
    if not map_dir.exists():
        return []

    # likely variant caller directories: anything under map_dir that isn't QC/Mapping/PreProcess
    candidates = []
    for p in map_dir.iterdir():
        if p.is_dir():
            if p.name.lower() in {"qc", "mapping", "preprocess"}:
                continue
            candidates.append(p)

    # fallback to map_dir itself
    return candidates if candidates else [map_dir]


def write_text(dst: Path, text: str) -> None:
    safe_mkdir(dst.parent)
    dst.write_text(text, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default=".", help="Path to repo root (default: current directory)")
    ap.add_argument("--basedir", default="runs", help="Base runs directory relative to repo root (default: runs)")
    ap.add_argument("--sample", required=True, help="Sample base name, e.g. HCC1395")
    ap.add_argument("--maptype", default="Bwa", choices=["Bwa", "Bowtie2", "NovoAlign"], help="Mapping type")
    ap.add_argument("--tumor-name", default=DEFAULT_TUMOR_NAME, help="Tumor label folder suffix (default: Tumor)")
    ap.add_argument("--germline-name", default=DEFAULT_GERMLINE_NAME, help="Germline label folder suffix (default: Germline)")
    ap.add_argument("--outdir", default=None, help="Output bundle directory (default: <repo>/analysis_report_bundle_<sample>_<timestamp>)")
    ap.add_argument("--copy", action="store_true", help="Copy files instead of symlinking (bundle becomes standalone but larger)")
    args = ap.parse_args()

    repo_root = Path(args.repo_root).resolve()
    runs_root = (repo_root / args.basedir).resolve()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir).resolve() if args.outdir else (repo_root / f"analysis_report_bundle_{args.sample}_{timestamp}").resolve()
    mode = "copy" if args.copy else "symlink"

    safe_mkdir(outdir)

    manifest: Dict = {
        "created_at": datetime.now().isoformat(),
        "repo_root": str(repo_root),
        "runs_root": str(runs_root),
        "sample": args.sample,
        "maptype": args.maptype,
        "mode": mode,
        "stages": {},
        "missing_roots": [],
        "notes": [],
    }

    # Define sample-specific paths
    tumor_dir = runs_root / f"{args.sample}_{args.tumor_name}"
    germ_dir = runs_root / f"{args.sample}_{args.germline_name}"

    # Stage roots
    def path_or_missing(p: Path, label: str) -> Path:
        if not p.exists():
            manifest["missing_roots"].append({"label": label, "path": str(p)})
        return p

    # Bundle each (Tumor, Germline) separately under top folders
    for sample_label, sample_dir in [("Tumor", tumor_dir), ("Germline", germ_dir)]:
        sample_out = outdir / sample_label
        safe_mkdir(sample_out)

        qc_root = path_or_missing(sample_dir / args.maptype / "QC", f"{sample_label}.QC")
        mapping_root = path_or_missing(sample_dir / args.maptype / "Mapping", f"{sample_label}.Mapping")
        preprocess_root = path_or_missing(sample_dir / args.maptype / "PreProcess", f"{sample_label}.PreProcess")

        # 1) QC
        if qc_root.exists():
            stage_bundle(
                label=f"{sample_label}/QC",
                search_root=qc_root,
                patterns=QC_PATTERNS,
                out_root=outdir,
                repo_root=repo_root,
                mode=mode,
                manifest=manifest,
            )

        # 2) Mapping
        if mapping_root.exists():
            stage_bundle(
                label=f"{sample_label}/Mapping",
                search_root=mapping_root,
                patterns=MAPPING_PATTERNS,
                out_root=outdir,
                repo_root=repo_root,
                mode=mode,
                manifest=manifest,
            )

        # 3) PreProcess (includes GATK4 folder, etc.)
        if preprocess_root.exists():
            stage_bundle(
                label=f"{sample_label}/PreProcess",
                search_root=preprocess_root,
                patterns=PREPROCESS_PATTERNS,
                out_root=outdir,
                repo_root=repo_root,
                mode=mode,
                manifest=manifest,
            )

        # 4) Variants: try to locate caller folders
        variant_dirs = guess_variant_dirs(sample_dir, args.maptype)
        if variant_dirs:
            # bundle from each candidate into its own subfolder
            for vdir in variant_dirs:
                # Use folder name (caller) if possible
                caller_name = vdir.name
                stage_bundle(
                    label=f"{sample_label}/Variants/{caller_name}",
                    search_root=vdir,
                    patterns=VARIANT_PATTERNS,
                    out_root=outdir,
                    repo_root=repo_root,
                    mode=mode,
                    manifest=manifest,
                )
        else:
            manifest["notes"].append(f"{sample_label}: could not find variant directories under {sample_dir / args.maptype}")

        # 5) Coverage: search broadly under maptype dir
        map_root = path_or_missing(sample_dir / args.maptype, f"{sample_label}.MapRoot")
        if map_root.exists():
            stage_bundle(
                label=f"{sample_label}/Coverage",
                search_root=map_root,
                patterns=COVERAGE_PATTERNS,
                out_root=outdir,
                repo_root=repo_root,
                mode=mode,
                manifest=manifest,
            )

        # 6) Logs: search under the sample run directory (not the entire repo)
        if sample_dir.exists():
            stage_bundle(
                label=f"{sample_label}/Logs",
                search_root=sample_dir,
                patterns=LOG_PATTERNS,
                out_root=outdir,
                repo_root=repo_root,
                mode=mode,
                manifest=manifest,
            )

    # Write manifest
    manifest_path = outdir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    # Write a quick guide
    guide = f"""\
Analysis Report Bundle: {args.sample} ({args.maptype})
Created: {manifest["created_at"]}
Mode: {mode}  (use --copy if you want a fully standalone bundle)

What this bundle is for:
  1) Scientific validity audit:
     - QC integrity (fastp/FastQC)
     - Mapping quality and BAM integrity
     - Duplicate metrics + BQSR evidence
     - Coverage/depth artifacts that affect confidence
     - Provenance/logs for reproducibility

  2) Biological interpretation:
     - Variant calls (VCFs) + caller logs
     - Any filtered/annotated tables if present
     - Coverage evidence for key loci

Important:
  - If you expected VCFs and none appear, open manifest.json and search for "Variants"
    to see what was found/missing. That usually means the variant calling step wrote
    outputs somewhere unexpected (or didn't run/filter/write due to an upstream error).

Bundle layout:
  Tumor/
    QC/
    Mapping/
    PreProcess/
    Variants/
    Coverage/
    Logs/
  Germline/
    QC/
    Mapping/
    PreProcess/
    Variants/
    Coverage/
    Logs/
  manifest.json

Next best step:
  - Send me the manifest.json plus the Tumor/Variants and Germline/Variants folders
    (or confirm they're empty). That’s enough to diagnose output-path issues quickly.
"""
    write_text(outdir / "BUNDLE_README.txt", guide)

    print(f"\n✅ Bundle created at:\n{outdir}")
    print(f"   - manifest: {manifest_path}")
    print("   - open BUNDLE_README.txt for how to use this bundle\n")


if __name__ == "__main__":
    main()