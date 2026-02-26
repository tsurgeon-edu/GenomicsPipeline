import os
import glob
import helpers
from log_command import log_command
from paths import GetPaths
from pathlib import Path


class GatkPreProcessing(object):

    def __init__(self, working_directory, map_type, sample_type, library_matching_id, thrds):
        self.get_paths = GetPaths()
        self.main_directory = working_directory
        self.folder_directory = working_directory + "/" + map_type
        self.working_directory = working_directory + "/" + map_type + "/PreProcess"
        self.map_type = map_type
        self.sample_type = sample_type
        self.library_matching_id = library_matching_id
        self.threads = thrds
        self.bundle_dir = self.get_paths.ref_dir
        # self.bundle_dir = self.get_paths.ref_dir + "hg19_bundle"
        self.file_list = []
        os.makedirs(self.working_directory, exist_ok=True)
        os.chdir(self.working_directory)

    def gatk3_realign_target_creator(self, lastbam):
        realign_target = str(lastbam).split(".")[0] + "_realign_target.intervals"
        bcal = "java -jar " + self.get_paths.gatk_path + " -T RealignerTargetCreator -nt " + \
               self.threads + " -R " + self.bundle_dir + "/ucsc.hg19.fasta -known " + \
               self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf -I " + lastbam + \
               " -o " + realign_target
        print(bcal)
        log_command(bcal, "Realign Target Creator", self.threads, "GatkPreProcessing")
        self.file_list.append(realign_target)
        return realign_target

    def gatk3_indel_realigner(self, lastbam, realign_target):

        realigned_last_bam = "IR_" + lastbam
        bcal = "java -jar " + self.get_paths.gatk_path + " -T IndelRealigner -R " + self.bundle_dir + \
               "/ucsc.hg19.fasta -known " + self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + \
               " -targetIntervals " + realign_target + " --noOriginalAlignmentTags -I " + lastbam + " -o " + \
               realigned_last_bam

        log_command(bcal, "Indel Realigner", self.threads, "GatkPreProcessing")
        self.file_list.append(realigned_last_bam)
        return realigned_last_bam

    def gatk3_base_recalibrator(self, lastbam):
        basequalityscore = str(lastbam).split(".")[0] + "_bqsr.grp"
        nct = " -nct " + str(self.threads)
        bcal = "java -jar " + self.get_paths.gatk_path + nct + " -T BaseRecalibrator -R " + self.bundle_dir +\
               "/ucsc.hg19.fasta -I " + lastbam + " -knownSites " + self.bundle_dir +\
               "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + " -o " + basequalityscore
        log_command(bcal, "Base Recalibrator", self.threads, "GatkPreProcessing")
        self.file_list.append(basequalityscore)
        return basequalityscore


    def gatk3_print_reads(self, lastbam, bqsr):
        nct = " -nct " + str(self.threads)

        aftercalibratorBam = "GATK_PR" + lastbam
        bcal = "java -jar " + self.get_paths.gatk_path + nct + " -T PrintReads -R " + self.bundle_dir + \
               "/ucsc.hg19.fasta -I " + lastbam + " --BQSR " + bqsr + " -o " + aftercalibratorBam
        log_command(bcal, "Print Reads", self.threads, "GatkPreProcessing")
        self.file_list.append(aftercalibratorBam)
        indexed = helpers.create_index(aftercalibratorBam, "Create Index by GATK_PrintReads", self.threads,
                                   "GatkPreProcess")
        self.file_list.append(indexed)


    def gatk4_base_recalibrator(self, lastbam):
        gatk4_dir = os.path.join(self.working_directory, "GATK4")
        os.makedirs(gatk4_dir, exist_ok=True)

        lastbam_name = os.path.basename(lastbam)
        stem = os.path.splitext(lastbam_name)[0]
        recal_table = os.path.join(gatk4_dir, f"{stem}_RECAL.table")

        ref_fa = os.path.join(self.bundle_dir, "Homo_sapiens_assembly38.fasta")
        bcal = self.get_paths.gatk4_path + " BaseRecalibrator -R " + ref_fa + \
               " -I " + lastbam + " --known-sites " + self.get_paths.mills_indel + \
               " --known-sites " + self.get_paths.dbsnp + " --known-sites " + self.get_paths.one_thousand_g + " -O " + \
               recal_table
        log_command(bcal, "Base Recalibrator", self.threads, "Gatk4PreProcessing")
        self.file_list.append(recal_table)
        return recal_table

    def gatk4_applybsqr(self, lastbam, recaltable):
        # Make sure output folder exists
        gatk4_dir = os.path.join(self.working_directory, "GATK4")
        os.makedirs(gatk4_dir, exist_ok=True)

        # Use ONLY the filename, never the full path
        lastbam_name = os.path.basename(lastbam)                 # e.g. MDUP_Bwa_HCC1395_MergedBAM.bam
        stem = os.path.splitext(lastbam_name)[0]                 # e.g. MDUP_Bwa_HCC1395_MergedBAM

        afterbqsrbam = os.path.join(gatk4_dir, f"{stem}.BQSR.bam")

        ref_fa = os.path.join(self.bundle_dir, "Homo_sapiens_assembly38.fasta")
        apply_command = (
            self.get_paths.gatk4_path
            + " ApplyBQSR -R " + ref_fa
            + " -I " + lastbam
            + " --bqsr-recal-file " + recaltable
            + " -O " + afterbqsrbam
        )

        log_command(apply_command, "ApplyBQSR", self.threads, "Gatk4PreProcessing")
        self.file_list.append(afterbqsrbam)

        indexed = helpers.create_index(afterbqsrbam, "Create Index by GATK_ApplyBSQR", self.threads, "GatkPreProcess")
        self.file_list.append(indexed)

    def run_gatks3(self, after_markdpl):

        realign_target = self.gatk3_realign_target_creator(after_markdpl)
        realigned_bam = self.gatk3_indel_realigner(after_markdpl, realign_target)
        basequality = self.gatk3_base_recalibrator(realigned_bam)
        self.gatk3_print_reads(realigned_bam, basequality)
        gatk_files = glob.glob("GATK_PR*.bam")
        return gatk_files

    def run_gatks4(self, after_markdpl):
        basequality = self.gatk4_base_recalibrator(after_markdpl)
        self.gatk4_applybsqr(after_markdpl, basequality)
        gatk_files = glob.glob(os.path.join("GATK4", "*.BQSR.bam"))
        return gatk_files


# if __name__ == "__main__":
#     os.chdir("/home/bioinformaticslab/Desktop/GitHub_Repos/test_files/NOB01/Bwa/PreProcess")
#     print(os.getcwd())
#     after_markdpl_file = glob.glob("MDUP_*.bam")
#     print(after_markdpl_file)
#     gatk_file_list = []
#     for file in after_markdpl_file:
#         gatk_pre_processing_step = GatkPreProcessing(
#             working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/test_files/NOB01",
#             map_type="Bwa", sample_type="Germline", library_matching_id="203", thrds="2")
#
#         return_files = gatk_pre_processing_step.run_gatks4(file)
#         print(return_files)
#         gatk_file_list.append(return_files)
#         print(gatk_file_list)