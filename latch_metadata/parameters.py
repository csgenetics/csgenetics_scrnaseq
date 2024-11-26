
from dataclasses import dataclass
import typing
import typing_extensions

from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import NextflowParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

from pathlib import Path
from enum import Enum

class ReferenceGenome(Enum):
    no_selection = ""
    grch38 = 'GRCh38'
    grcm39 = 'GRCm39'
    mouse_human_mix = 'mouse_human_mix'

@dataclass(frozen=True)
class Sample:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile
    manual_cellcaller_threshold: typing.Optional[float] = None
    hsap_manual_cell_caller_threshold: typing.Optional[float] = None
    mmus_manual_cell_caller_threshold: typing.Optional[float] = None

generated_parameters = {
    'input_csv': NextflowParameter(
        display_name="Input CSV",
        batch_table_column=False,
        type=typing.List[Sample],
        samplesheet=True,
        samplesheet_type='csv',
        section_title=None,
        description=
    """
    This input defines the sample names and their corresponding fastq files. There should be one row for every fastq.gz file pair.
    Rows with the same sample name will have their R1 and R2 fastq.gz files merged together (E.g. for multiple lanes of sequencing data per sample).
    """
    ),
    'outdir': NextflowParameter(
        display_name="Output directory",
        batch_table_column=False,
        type=LatchDir,
        description='The output directory where the results will be exported.',
        results_paths=[
            Path("/"), Path("/report/"), Path("/pipeline_info/"), Path("/RSeQC/"), Path("/multiqc/")
        ]
    ),
    'genome': NextflowParameter(
        display_name="Genomic reference",
        batch_table_column=False,
        type=ReferenceGenome,
        default=ReferenceGenome.no_selection,
        description="A string representing the preconfigured genomic resources to use. If using a preconfigured reference, do not set the star_index, gtf_dir or mitochondra_chromosome variables, they will automatically be set to the correct values for you. If using a custom set of references, set this to 'custom'.",
    ),
    'star_index': NextflowParameter(
        display_name="STAR index dir path",
        batch_table_column=False,
        type=typing.Optional[LatchDir],
        description="The path to the STAR index directory."
    ),
    'gtf': NextflowParameter(
        display_name="GTF file path",
        batch_table_column=False,
        type=typing.Optional[LatchFile],
        description="The path to the GTF file used to build the STAR index."
    ),
    'mitochondria_chromosome': NextflowParameter(
        display_name="Mitochondria chromosome identifier",
        batch_table_column=False,
        type=str,
        default="",
        description="The name of the mitochondria chromosome in the GTF file e.g. 'MT' or 'chrM'."
    ),
    'minimum_count_threshold': NextflowParameter(
        display_name="Minimum count threshold for cell calling",
        batch_table_column=False,
        type=float,
        default=100.0,
        description="The minimum total counts for a given barcode to classify it as a cell."
    ),
}

