from latch.resources.launch_plan import LaunchPlan
from latch.types.file import LatchFile
from .entrypoint import nf_cs_genetics_simplecell_pipeline, ReferenceGenome
from wf.entrypoint import Sample

LaunchPlan(
    nf_cs_genetics_simplecell_pipeline,
    "4 sample PBMC demo",
    {
        "input_csv": [
            Sample(
                sample="sample1",
                fastq_1=LatchFile("s3://latch-public/test-data/36794/241023_Dorothy_HIVE1802_1400pM_Sample273_S273_L001_R1_001.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/36794/241023_Dorothy_HIVE1802_1400pM_Sample273_S273_L001_R2_001.fastq.gz"),
                manual_cellcaller_threshold = None,
                hsap_manual_cell_caller_threshold = None,
                mmus_manual_cell_caller_threshold = None
            ),
            Sample(
                sample="sample2",
                fastq_1=LatchFile("s3://latch-public/test-data/36794/241023_Dorothy_HIVE1802_1400pM_Sample275_S275_L001_R1_001.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/36794/241023_Dorothy_HIVE1802_1400pM_Sample275_S275_L001_R2_001.fastq.gz"),
                manual_cellcaller_threshold = None,
                hsap_manual_cell_caller_threshold = None,
                mmus_manual_cell_caller_threshold = None
            ),
            Sample(
                sample="sample3",
                fastq_1=LatchFile("s3://latch-public/test-data/36794/241023_Dorothy_HIVE1802_1400pM_Sample277_S277_L001_R1_001.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/36794/241023_Dorothy_HIVE1802_1400pM_Sample277_S277_L001_R2_001.fastq.gz"),
                manual_cellcaller_threshold = None,
                hsap_manual_cell_caller_threshold = None,
                mmus_manual_cell_caller_threshold = None
            ),
            Sample(
                sample="sample4",
                fastq_1=LatchFile("s3://latch-public/test-data/36794/241023_Dorothy_HIVE1802_1400pM_Sample280_S280_L001_R1_001.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/36794/241023_Dorothy_HIVE1802_1400pM_Sample280_S280_L001_R2_001.fastq.gz"),
                manual_cellcaller_threshold = None,
                hsap_manual_cell_caller_threshold = None,
                mmus_manual_cell_caller_threshold = None
            )
        ],
        "genome": ReferenceGenome.grch38
    }
)