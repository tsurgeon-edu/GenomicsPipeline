import argparse
import glob
import os
import shutil
import sys
import variant_annotation
import variant_calling

def call_variant_caller(working_directory, tumor_bam, germline_bam, var_maptype, var_variantcaller, threads_p,
                        s_name, tumor_only, tumor_interval=None, germline_interval=None, to_temp=[False, None]):

    if to_temp[0]:
        os.makedirs(os.path.dirname(to_temp[1]), exist_ok=True)
        if os.path.isdir(working_directory):
            shutil.copytree(working_directory, to_temp[1], dirs_exist_ok=True)
        else:
            shutil.copy2(working_directory, to_temp[1])
        wd = to_temp[1]
    else:
        wd = working_directory

    if wd[-1] == "/" or wd[-1] == "\\":
        wd = wd[:-1]

    gm = germline_bam

    if gm[-1] == "/" or gm[-1] == "\\":
        gm = gm[:-1]

    gm_interval = germline_interval

    if gm_interval and (gm_interval[-1] == "/" or gm_interval[-1] == "\\"):
        gm_interval = gm_interval[:-1]

    if var_variantcaller == "Mutect2" or var_variantcaller == "Mutect2_gatk3":
        pipeline2 = variant_calling.VariantCall(variant_caller=var_variantcaller, thrds=threads_p, map_type=var_maptype,
                                                germline_bam=gm, germline_interval=gm_interval, wd=wd,
                                                tumor_bam=tumor_bam, tumor_interval=tumor_interval,
                                                sample_name=s_name, tumor_only=tumor_only)
    else:
        pipeline2 = variant_calling.VariantCall(variant_caller=var_variantcaller, thrds=threads_p, map_type=var_maptype,
                                                germline_bam=gm, germline_interval=None, wd=wd, tumor_bam=tumor_bam,
                                                tumor_interval=None, sample_name=s_name,  tumor_only=tumor_only)
    pipeline2_success = pipeline2.run_pipeline()

    # annotate = variant_annotation.VariantAnnotation(variant_annotater="Annovar", thread_v=threads_p, wd=pipeline2_success, sample_name=s_name, will_annotate=[], annotate_all=True)
    #
    # annotate.run_annotation()

    return pipeline2_success

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run somatic variant calling for a Tumor/Normal(Germline) pair."
    )

    parser.add_argument(
        "--basedir",
        default=os.getcwd(),
        help="Base directory containing sample folders (default: current directory).",
    )

    parser.add_argument(
        "--sample",
        required=True,
        help="Base sample ID (e.g., HCC1395). Folders inferred as '<sample>_Tumor' and '<sample>_Germline'.",
    )

    parser.add_argument(
        "--maptype",
        default="Bwa",
        choices=["Bwa", "Bowtie2", "NovoAlign"],
        help="Mapper used to generate BAMs (default: Bwa).",
    )

    parser.add_argument(
        "--variantcaller",
        default="Mutect2",
        choices=["Mutect2", "Mutect2_gatk3", "Strelka", "SomaticSniper", "Varscan"],
        help="Variant caller to run (default: Mutect2).",
    )

    parser.add_argument(
        "--threads",
        default=str(os.cpu_count() or 4),
        help="Number of threads (default: detected CPU count).",
    )

    parser.add_argument(
        "--tumor-only",
        action="store_true",
        help="Run in tumor-only mode (no germline required).",
    )

    parser.add_argument(
        "--sample-name",
        default=None,
        help="Sample label used in output filenames (default: same as --sample).",
    )

    # Optional explicit overrides (if you want to point to exact BAMs)
    parser.add_argument("--tumor-bam", default=None, help="Explicit path to tumor BAM.")
    parser.add_argument("--germline-bam", default=None, help="Explicit path to germline BAM.")

    # Optional intervals (only used by Mutect2 / Mutect2_gatk3 in your function)
    parser.add_argument("--tumor-interval", default=None, help="Tumor interval file (optional).")
    parser.add_argument("--germline-interval", default=None, help="Germline interval file (optional).")

    args = parser.parse_args()

    basedir = os.path.abspath(os.path.expanduser(args.basedir))
    sample_name = args.sample_name or args.sample

    tumor_dir = os.path.join(basedir, f"{args.sample}_Tumor")
    germline_dir = os.path.join(basedir, f"{args.sample}_Germline")

    if not os.path.isdir(tumor_dir):
        print(f"ERROR: Tumor folder does not exist: {tumor_dir}", file=sys.stderr)
        sys.exit(2)

    if not args.tumor_only and not os.path.isdir(germline_dir):
        print(f"ERROR: Germline folder does not exist: {germline_dir}", file=sys.stderr)
        sys.exit(2)

    def pick_bam(search_root: str) -> str | None:
        """
        Heuristic BAM picker. Does NOT change pipeline behavior; only chooses input BAM path.
        Preference order:
          1) PreProcess/*GATK4*.bam
          2) PreProcess/*.bam
          3) Any *.bam under sample folder
        """
        candidates = []
        candidates += glob.glob(os.path.join(search_root, "**", "PreProcess", "*GATK4*.bam"), recursive=True)
        candidates += glob.glob(os.path.join(search_root, "**", "PreProcess", "*.bam"), recursive=True)
        candidates += glob.glob(os.path.join(search_root, "**", "*.bam"), recursive=True)

        # Drop obvious non-final BAMs if present (best-effort)
        filtered = [p for p in candidates if not p.endswith(".bai")]
        if not filtered:
            return None

        # Prefer larger BAMs (often final merged BAM) as a simple heuristic
        filtered.sort(key=lambda p: os.path.getsize(p), reverse=True)
        return filtered[0]

    tumor_bam = args.tumor_bam or pick_bam(tumor_dir)
    if not tumor_bam or not os.path.isfile(tumor_bam):
        print(
            "ERROR: Could not find a tumor BAM automatically.\n"
            f"Looked under: {tumor_dir}\n"
            "Fix: provide --tumor-bam /full/path/to/tumor.bam",
            file=sys.stderr,
        )
        sys.exit(2)

    if args.tumor_only:
        germline_bam = args.germline_bam or tumor_bam  # placeholder; VariantCall may ignore in tumor-only mode
    else:
        germline_bam = args.germline_bam or pick_bam(germline_dir)
        if not germline_bam or not os.path.isfile(germline_bam):
            print(
                "ERROR: Could not find a germline BAM automatically.\n"
                f"Looked under: {germline_dir}\n"
                "Fix: provide --germline-bam /full/path/to/germline.bam",
                file=sys.stderr,
            )
            sys.exit(2)

    # Use the tumor folder as the working directory for outputs (common pattern)
    working_directory = tumor_dir

    call_variant_caller(
        working_directory=working_directory,
        tumor_bam=tumor_bam,
        germline_bam=germline_bam,
        var_maptype=args.maptype,
        var_variantcaller=args.variantcaller,
        threads_p=args.threads,
        s_name=sample_name,
        tumor_only="Yes" if args.tumor_only else "No",
        tumor_interval=args.tumor_interval,
        germline_interval=args.germline_interval,
    )

    