import os
import subprocess
from paths import GetPaths


class PonCreation(object):
    def __init__(self, normal_bam, normal_interval=None, ref="hg38"):
        self.get_paths = GetPaths(ref=ref)
        self.normal_bam = normal_bam
        self.normal_interval = normal_interval

        # output name beside the BAM (or change to your preferred output folder)
        bam_base = os.path.basename(normal_bam).replace(".bam", "")
        self.output_vcf = f"{bam_base}.pon.vcf.gz"

    def create_normal_for_pon(self):
        """
        Create a single-sample VCF to later combine into a Panel of Normals (PoN).

        NOTE: Modern GATK4 PoN workflows typically use Mutect2:
          gatk Mutect2 --max-mnp-distance 0 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
            ... then CreateSomaticPanelOfNormals

        This provides a sane starting point; adjust if your pipeline expects a different PoN format.
        """
        ref_fasta = os.path.join(self.get_paths.ref_dir, "Homo_sapiens_assembly38.fasta")

        cmd = [
            self.get_paths.gatk4_path,
            "Mutect2",
            "-R", ref_fasta,
            "-I", self.normal_bam,
            "-O", self.output_vcf,
        ]

        # Optional: targets/intervals
        if self.normal_interval:
            cmd += ["-L", self.normal_interval]

        # Optional: known sites (dbSNP)
        if self.get_paths.dbsnp:
            cmd += ["--dbsnp", self.get_paths.dbsnp]

        # Optional: COSMIC
        if self.get_paths.cosmic:
            cmd += ["--cosmic", self.get_paths.cosmic]

        print("Running:", " ".join(cmd))
        subprocess.check_call(cmd)

    def combine_pon(self, input_vcfs, output_pon="Mutect2_PON.vcf.gz"):
        """
        Combine per-normal VCFs into a PoN.
        GATK4 tool: CreateSomaticPanelOfNormals
        """
        cmd = [self.get_paths.gatk4_path, "CreateSomaticPanelOfNormals", "-O", output_pon]
        for v in input_vcfs:
            cmd += ["-V", v]

        print("Running:", " ".join(cmd))
        subprocess.check_call(cmd)