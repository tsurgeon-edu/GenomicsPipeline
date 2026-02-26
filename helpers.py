import glob
import os
import shutil
import subprocess
import tempfile
from log_command import log_command
from paths import GetPaths

def picard_cmd(subcommand_and_args: str) -> str:
    """
    Build a working Picard command regardless of whether picard_path is:
      - a wrapper executable (preferred), or
      - a picard.jar path.
    """
    picard = GetPaths().picard_path
    if not picard:
        raise RuntimeError("picard not found (set on PATH or ensure picard.jar exists).")

    picard = str(picard)
    if picard.endswith(".jar"):
        return f"java -jar {picard} {subcommand_and_args}"
    else:
        # wrapper script/binary
        return f"{picard} {subcommand_and_args}"

def get_fastq():
    """
    Get fastq files names with their extension

    Returns
    -------
    list
        A list of fastq files inside of given working directory.
    """

    all_fastq_files = glob.glob("*fastq.gz")
    split_names_v = [os.path.splitext(os.path.splitext(i)[0])[0] for i in all_fastq_files]
    return split_names_v


def get_info(sample_type, fastq_list, trimmed=False):

    """
    Prepare set of information in order to used in next steps especially creating read group in mapping function.

    Returns
    -------
    dict
        list of unique information inside dictionary
    """
    sample_id, germline_dna, index_seq, lanes, pairs_r, n_of_seq = (set() for i in range(6))
    if sample_type == "Tumor":
        for i in fastq_list:
            if trimmed:
                sample_id.add(i.split("_")[1])
                index_seq.add(i.split("_")[2])
                lanes.add(i.split("_")[3])
                pairs_r.add(i.split("_")[4])
                n_of_seq.add(i.split("_")[5])
            else:
                sample_id.add(i.split("_")[0])
                index_seq.add(i.split("_")[1])
                lanes.add(i.split("_")[2])
                pairs_r.add(i.split("_")[3])
                n_of_seq.add(i.split("_")[4])

        list_with_info = {"Sample_ID": list(sample_id), "Index": list(index_seq), "Lanes": list(lanes),
                          "Pairs": list(pairs_r), "Number_of_seq": list(n_of_seq)}
        return list_with_info
    elif sample_type == "Germline" or sample_type == "Normal":

        for i in fastq_list:
            if trimmed:
                sample_id.add(i.split("_")[1])
                germline_dna.add(i.split("_")[2])
                index_seq.add(i.split("_")[3])
                lanes.add(i.split("_")[4])
                pairs_r.add(i.split("_")[5])
                n_of_seq.add(i.split("_")[6])
            else:
                sample_id.add(i.split("_")[0])
                germline_dna.add(i.split("_")[1])
                index_seq.add(i.split("_")[2])
                lanes.add(i.split("_")[3])
                pairs_r.add(i.split("_")[4])
                n_of_seq.add(i.split("_")[5])

        list_with_info = {"Sample_ID": list(sample_id), "Germline": list(germline_dna), "Index": list(index_seq),
                          "Lanes": list(lanes), "Pairs": list(pairs_r), "Number_of_seq": list(n_of_seq)}
        return list_with_info
    else:
        print("raise error and ask again for a valid sample type")


def create_folder(working_directory, all_files, map_type=None, step="Other", folder_directory=None):
    # Always include the log if present
    files = set(all_files)
    files.add("log_file.txt")

    # QC stays inside the working directory (because Mapping switches into QC/)
    if step == "QC":
        mk_dir = os.path.join(working_directory, "QC")
        os.makedirs(mk_dir, exist_ok=True)
        for f in files:
            src = os.path.join(working_directory, f)
            if os.path.exists(src):
                shutil.move(src, os.path.join(mk_dir, f))
        return

    # Everything else goes under: folder_directory/<map_type>/<step>
    if folder_directory is None:
        folder_directory = working_directory

    if map_type is None:
        # fallback: old behavior
        mk_dir = os.path.join(folder_directory, step)
    else:
        mk_dir = os.path.join(folder_directory, map_type, step)

    os.makedirs(mk_dir, exist_ok=True)

    for f in files:
        # keep fastqs in place (matches your existing intent)
        if f.endswith(".gz"):
            continue
        src = os.path.join(working_directory, f)
        if os.path.exists(src):
            shutil.move(src, os.path.join(mk_dir, f))

def create_index(lastbam, function, threads, step):
    indexcol = picard_cmd(f"BuildBamIndex I={lastbam}")
    log_command(indexcol, function, threads, step)
    bai = os.path.splitext(lastbam)[0] + ".bai"
    return bai

def delete_files_from_folder(base_dir: str, map_type: str, step: str, keep_files):
    """
    Delete intermediate files inside: <base_dir>/<map_type>/<step>
    while preserving the files listed in keep_files (plus their .bai indexes, and log_file.txt).

    This matches the call signature used in run_pipeline_mapping.py:
        helpers.delete_files_from_folder(wd, mt, "Mapping", mapping_files)
    """
    folder = os.path.join(base_dir, map_type, step)
    if not os.path.isdir(folder):
        return []

    keep = set(keep_files or [])
    keep.add("log_file.txt")

    # If keep contains BAMs, also keep corresponding BAIs
    keep_bais = set()
    for f in keep:
        if f.endswith(".bam"):
            keep_bais.add(os.path.splitext(f)[0] + ".bai")
    keep |= keep_bais

    deleted = []
    for name in os.listdir(folder):
        path = os.path.join(folder, name)

        # never touch dirs here
        if os.path.isdir(path):
            continue

        if name in keep:
            continue

        try:
            os.remove(path)
            deleted.append(name)
        except FileNotFoundError:
            pass

    return deleted

def get_sample_name(bamfile):
    cmd = ["samtools", "view", "-H", bamfile]
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT).decode("utf-8", errors="replace")
        for line in out.splitlines():
            if line.startswith("@RG"):
                for field in line.split("\t"):
                    if field.startswith("SM:"):
                        return field[3:]
    except Exception:
        return False

def delete_file_custom(file):
    return True
