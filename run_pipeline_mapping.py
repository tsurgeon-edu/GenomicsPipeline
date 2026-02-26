import argparse
import gatk_pre_processing
import glob
import helpers
import mapping
import os
import pre_processing
import qc_trim
import sys


def callmapping(var_maptype, var_sampletype, working_directory, library, threads, var_gatk_tools, issplitchr, trim,
                middle_files="Yes"):
    mt = var_maptype
    if middle_files == "Yes":
        mdf_keep = True
    else:
        mdf_keep = False
    st = var_sampletype
    wd = working_directory
    if wd[-1] == "/" or wd[-1] == "\\":
        wd = wd[:-1]
    lb = library
    th = threads
    gt = var_gatk_tools
    sc = issplitchr
    tr = trim
    os.chdir(wd)

    fastq_list = helpers.get_fastq()
    info_dict = helpers.get_info(st, fastq_list)

    if tr == "Yes":
        if not os.path.exists(wd + "/QC"):
            qc = qc_trim.QC(wd, st, th, fastq_list, info_dict, mt)
            qc.run_qc()
    else:
        if os.path.exists(wd+"/QC"):
            tr = "Yes"


    mapping_step = mapping.Mapping(working_directory=wd, map_type=mt, sample_type=st, library_matching_id=lb,
                                   thrds=th, trim=tr)

    mapping_files = mapping_step.mapping()
    #mapping_files = ["SortedBAM_Bwa_NOB01_AACGTGA_L001_001.bam"]

    if not mdf_keep:
        helpers.delete_files_from_folder(wd, mt, "Mapping", mapping_files)


    print("---------------------------")
    print(mapping_files)
    pre_processing_step = pre_processing.PreProcessing(working_directory=wd, map_type=mt, sample_type=st,
                                                       library_matching_id=lb, thrds=th, issplitchr=sc)

    print("---------------------------")
    print(fastq_list)
    print(info_dict)
    gatk_file_list = []
    if gt == "Yes":
        if issplitchr != "No":
            mark_duplicate_file = pre_processing_step.pre_process(info_dict, mapping_files)
            for file in mark_duplicate_file:
                gatk_pre_processing_step = gatk_pre_processing.GatkPreProcessing(working_directory=wd, map_type=mt,
                                                                                 sample_type=st, library_matching_id=lb,
                                                                                 thrds=th)
                return_files = gatk_pre_processing_step.run_gatks4(file)
                print(return_files)
                gatk_file_list.append(return_files)
                print(gatk_file_list)

        else:
            mark_duplicate_file = pre_processing_step.pre_process(info_dict, mapping_files)
            gatk_pre_processing_step = gatk_pre_processing.GatkPreProcessing(working_directory=wd, map_type=mt,
                                                                             sample_type=st, library_matching_id=lb,
                                                                             thrds=th)
            gatk_files = gatk_pre_processing_step.run_gatks4(mark_duplicate_file)

            if not mdf_keep:
                helpers.delete_files_from_folder(wd, mt, "PreProcess", gatk_files)
    else:
        mark_duplicate_file = pre_processing_step.pre_process(info_dict, mapping_files)
        if not mdf_keep:
            helpers.delete_files_from_folder(wd, mt, "PreProcess", mark_duplicate_file)

    return True

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run mapping + preprocessing pipeline for one sample-type folder."
    )

    parser.add_argument(
        "--basedir",
        default=os.getcwd(),
        help="Base directory containing sample folders (default: current directory).",
    )

    parser.add_argument(
        "--sample",
        required=True,
        help="Base sample ID (e.g., HCC1395). Folder is inferred as '<sample>_<sampletype>'.",
    )

    parser.add_argument(
        "--sampletype",
        required=True,
        choices=["Tumor", "Germline", "Normal"],
        help="Sample type. Used for both parsing FASTQ names and folder inference.",
    )

    parser.add_argument(
        "--maptype",
        default="Bwa",
        choices=["Bwa", "Bowtie2", "NovoAlign"],
        help="Mapper to use (default: Bwa).",
    )

    parser.add_argument(
        "--threads",
        default=str(os.cpu_count() or 4),
        help="Number of threads (default: detected CPU count). Keep as a number; will be converted to string.",
    )

    parser.add_argument(
        "--gatk",
        action="store_true",
        help="Enable GATK preprocessing (sets var_gatk_tools='Yes').",
    )

    parser.add_argument(
        "--split-chr",
        action="store_true",
        help="Enable chromosome splitting behavior (sets issplitchr != 'No').",
    )

    parser.add_argument(
        "--no-trim",
        action="store_true",
        help="Disable trimming/QC step (sets trim='No').",
    )

    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Keep intermediate files (sets middle_files='Yes').",
    )

    parser.add_argument(
        "--library",
        default="1",
        help="Library matching id used by the pipeline (default: '1').",
    )

    args = parser.parse_args()

    sample_folder = f"{args.sample}_{args.sampletype}"
    basedir = os.path.abspath(os.path.expanduser(args.basedir))
    wd = os.path.join(basedir, sample_folder)

    if not os.path.isdir(wd):
        print(f"ERROR: Working directory does not exist: {wd}", file=sys.stderr)
        sys.exit(2)

    # The pipeline's helpers.get_fastq() uses glob('*fastq.gz') in the working directory
    fastqs = glob.glob(os.path.join(wd, "*fastq.gz"))
    if not fastqs:
        print(
            f"ERROR: No FASTQ files found in {wd} matching '*fastq.gz'.\n"
            f"Expected files like: {args.sample}_S1_L001_R1_001.fastq.gz (Tumor)\n"
            f"or {args.sample}_Germline_S1_L001_R1_001.fastq.gz (Germline).",
            file=sys.stderr,
        )
        sys.exit(2)

    callmapping(
        working_directory=wd,
        var_maptype=args.maptype,
        var_sampletype=args.sampletype,
        library=args.library,
        threads=str(args.threads),
        var_gatk_tools="Yes" if args.gatk else "No",
        issplitchr="Yes" if args.split_chr else "No",  # anything other than "No" triggers the split path
        trim="No" if args.no_trim else "Yes",
        middle_files="Yes" if args.keep_intermediate else "No",
    )
