import os
import shutil


class GetPaths(object):
    def _resolve_picard(self):
        """
        Returns either:
          - a picard wrapper/executable path (preferred), OR
          - a picard.jar path if wrapper not available.
        """
        # Prefer a runnable wrapper on PATH (Homebrew often installs this)
        p = shutil.which("picard")
        if p:
            return p

        # Fall back to your explicit default wrapper location
        wrapper_default = "/usr/local/bin/picard"
        if os.path.exists(wrapper_default):
            return wrapper_default

        # Otherwise try common Picard jar locations
        jar_candidates = [
            "/usr/local/share/picard/picard.jar",
            "/opt/homebrew/share/picard/picard.jar",
        ]
        for jar in jar_candidates:
            if os.path.exists(jar):
                return jar

        return None

    def __init__(self, ref="hg38"):

        # Reference bundle roots (set these to where YOUR reference bundles live)
        HG38_BUNDLE = os.environ["HG38_BUNDLE"]
        HG19_BUNDLE = os.environ["HG19_BUNDLE"]
        ANNOVAR_DB  = os.environ["ANNOVAR_DB"]

        # normalize to always end with "/"
        HG38_BUNDLE = os.path.join(HG38_BUNDLE, "")
        HG19_BUNDLE = os.path.join(HG19_BUNDLE, "")
        ANNOVAR_DB  = os.path.join(ANNOVAR_DB, "")

        # Helper: find executables on PATH, fall back to a default, or keep None
        def which(cmd, default=None):
            return shutil.which(cmd) or default

        # Tools (prefer command on PATH)
        self.picard_path   = self._resolve_picard()
        self.gatk4_path    = which("gatk",     "/usr/local/bin/gatk")

        # This pipeline also mentions GATK 3 (GenomeAnalysisTK.jar) in the template.
        # If the code really calls GATK3, you'll need that JAR. Otherwise set None.
        self.gatk_path     = os.environ.get("GATK3_JAR", None)  # e.g. /path/to/GenomeAnalysisTK.jar

        self.varscan_path  = os.environ.get("VARSCAN_JAR", None)  # e.g. /path/to/VarScan.v2.3.9.jar
        self.fastqc        = which("fastqc")
        self.fastp         = which("fastp")

        # If your pipeline uses these, set them via env vars or absolute paths
        self.strelka = "/usr/local/bin/configureStrelkaSomaticWorkflow.py"
        self.novoalign = os.path.join("/usr/local/share/novocraft", "")
        self.somaticsniper = "/usr/local/bin/bam-somaticsniper"

        # ANNOVAR (if you created wrappers like table_annovar / annotate_variation)
        self.annovar       = os.environ.get("ANNOVAR_HOME", "/usr/local/share/annovar")
        self.table_annovar = which("table_annovar")
        self.convert2annovar = which("convert2annovar")
        self.annotate_variation = which("annotate_variation")

        if ref == "hg38":
            self.ref_dir = HG38_BUNDLE
            self.dbsnp = os.path.join(HG38_BUNDLE, "dbsnp_146.hg38.vcf.gz")
            self.mills_indel = os.path.join(HG38_BUNDLE, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
            self.one_thousand_g = os.path.join(HG38_BUNDLE, "1000G_phase1.snps.high_confidence.hg38.vcf.gz")
            self.cosmic = os.environ.get("COSMIC_HG38", None)
            self.annovar_db = ANNOVAR_DB
        else:
            self.ref_dir = HG19_BUNDLE
            self.dbsnp = os.path.join(HG19_BUNDLE, "dbsnp_138.hg19.vcf.gz")
            self.mills_indel = os.path.join(HG19_BUNDLE, "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz")
            self.one_thousand_g = os.path.join(HG19_BUNDLE, "1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz")
            self.cosmic = os.environ.get("COSMIC_HG19", None)
            self.annovar_db = ANNOVAR_DB
            