import sys
from dataclasses import dataclass
from enum import Enum
import os
import subprocess
import requests
import shutil
from pathlib import Path
import typing
import typing_extensions

from latch.resources.workflow import workflow
from latch.resources.tasks import nextflow_runtime_task, custom_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.ldata.path import LPath
from latch.executions import report_nextflow_used_storage
from latch_cli.nextflow.workflow import get_flag
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.utils import urljoins
from latch.types import metadata
from flytekit.core.annotation import FlyteAnnotation
from latch import message
from functools import partial

from latch_cli.services.register.utils import import_module_by_path

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)

@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


class ReferenceGenome(Enum):
    no_selection = ''
    grch38 = 'GRCh38'
    grcm39 = 'GRCm39'
    mouse_human_mix = 'mouse_human_mix'




@dataclass
class Sample:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile
    manual_cellcaller_threshold: typing.Optional[float]
    hsap_manual_cell_caller_threshold: typing.Optional[float]
    mmus_manual_cell_caller_threshold: typing.Optional[float]

# Depending on whether the workflow is single or mixed species, the input CSV will have different columns for setting the cell caller thresholds.
# The `SingleSpeciesSample` and `MixedSpeciesSample` dataclasses are used to define the input CSV schema for the two different workflows.
# They will be generated by the curate_samplesheet task and passed into the nextflow_runtime task.
@dataclass
class SingleSpeciesSample:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile
    manual_cellcaller_threshold: typing.Optional[float]

@dataclass
class MixedSpeciesSample:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile
    hsap_manual_cell_caller_threshold: typing.Optional[float]
    mmus_manual_cell_caller_threshold: typing.Optional[float]

input_csv_construct_samplesheet = metadata._nextflow_metadata.parameters['input_csv'].samplesheet_constructor

# Function will check that a valid selection has been made by the user.
# I.e. either the genome as been selected from one of our curated genomes or the user has provided a STAR index and GTF file.
# NOT both. The function will return the genome, star_index, gtf and mitochondria_chromosome variables after adjustment (if necessary).
@nextflow_runtime_task(cpu=1, memory=4, storage_gib=10)
def validate_genome_selection(genome: ReferenceGenome, star_index: typing.Optional[LatchDir], gtf: typing.Optional[LatchFile], mitochondria_chromosome: str) -> typing.Tuple[ReferenceGenome, typing.Optional[LatchDir], typing.Optional[LatchFile], str, bool]:
    single_species = True
    if genome == ReferenceGenome.no_selection:
        # No curated genome has been selected so we are expecting a user-supplied STAR index and GTF file.
        if star_index is None or gtf is None:
            raise ValueError("If you are not using a curated genome, you must provide a STAR index and GTF file.")
        message(typ='info', data={'title': "Genome configuration", "body":"Custom genome configured"})
    else:
        # A curated genome has been selected so the user should not provide a STAR index or GTF file.
        if star_index is not None or gtf is not None:
            raise ValueError("If you are using a preconfigured genome, you cannot provide a STAR index or GTF file.")
        message(typ='info', data={'title': "Genome configuration", "body":"Curated genome configured"})
        
    if genome == ReferenceGenome.mouse_human_mix:
        single_species = False

    return genome, star_index, gtf, mitochondria_chromosome, single_species

@nextflow_runtime_task(cpu=1, memory=4, storage_gib=10)
def curate_samplesheet(input_csv: typing.List[Sample], single_species: bool) -> typing.Union[typing.List[SingleSpeciesSample], typing.List[MixedSpeciesSample]]:
    if single_species:
        curated_samplesheet = [SingleSpeciesSample(sample=sample.sample, fastq_1=sample.fastq_1, fastq_2=sample.fastq_2, manual_cellcaller_threshold=sample.manual_cellcaller_threshold) for sample in input_csv]
        message(typ='info', data={'title': "Species configuration", "body":f"Single species configured: {str(curated_samplesheet)}"})
        return curated_samplesheet
    else:
        curated_samplesheet = [MixedSpeciesSample(sample=sample.sample, fastq_1=sample.fastq_1, fastq_2=sample.fastq_2, hsap_manual_cell_caller_threshold=sample.hsap_manual_cell_caller_threshold, mmus_manual_cell_caller_threshold=sample.mmus_manual_cell_caller_threshold) for sample in input_csv]
        message(typ='info', data={'title': "Species configuration", "body":f"Mixed species configured: {str(curated_samplesheet)}"})
        return curated_samplesheet

@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(pvc_name: str, input_csv: typing.Union[typing.List[SingleSpeciesSample], typing.List[MixedSpeciesSample]], outdir: LatchDir, star_index: typing.Optional[LatchDir], gtf: typing.Optional[LatchFile], genome: ReferenceGenome, mitochondria_chromosome: str, minimum_count_threshold: float, ) -> None:
    shared_dir = Path("/nf-workdir")

    exec_name = _get_execution_name()
    if exec_name is None:
        print("Failed to get execution name.")
        exec_name = "unknown"

    latch_log_dir = urljoins("latch:///your_log_dir/nf_cs_genetics_simplecell_pipeline", exec_name)
    print(f"Log directory: {latch_log_dir}")

    from latch.executions import add_execution_results

    results = []
    results.append(os.path.join(outdir.remote_path, ''))
    results.append(os.path.join(outdir.remote_path, 'report'))
    results.append(os.path.join(outdir.remote_path, 'pipeline_info'))
    results.append(os.path.join(outdir.remote_path, 'RSeQC'))
    results.append(os.path.join(outdir.remote_path, 'multiqc'))
    add_execution_results(results)

    # The below call to input_csv_contruct_samplesheet will write out the csv
    # file that will be supplied to the NextflowPipeline to the input_csv flag.
    # The input_csv_contruct_samplesheet function was originally created
    # with the t (type) set to <Sample>. We need to update the type to 
    # SingleSpeciesSample or MixedSpeciesSample depending on the workflow
    # so that the correct cell caller count threshold columns are written out. 
    if isinstance(input_csv[0], SingleSpeciesSample):
        t = SingleSpeciesSample
    else:
        t = MixedSpeciesSample

    input_csv_samplesheet = input_csv_construct_samplesheet(input_csv, t=t)

    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        "work",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared_dir,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    profile_list = ['docker']
    if False:
        profile_list.extend([p.value for p in execution_profiles])

    if len(profile_list) == 0:
        profile_list.append("standard")

    profiles = ','.join(profile_list)

    if star_index is not None:
        message(typ='info', data={'title': "Pipeline Launch", "body":" Launching Nextflow with custom genome configuration."})
        cmd = [
            "/root/nextflow",
            "run",
            str(shared_dir / "main.nf"),
            "-work-dir",
            str(shared_dir),
            "-profile",
            profiles,
            "-c",
            "latch.config",
            "-resume",
            "--genome", "custom",
            *get_flag('input_csv', input_csv_samplesheet),
            *get_flag('outdir', outdir),
            *get_flag('star_index', star_index),
            *get_flag('gtf', gtf),
            *get_flag('mitochondria_chromosome', mitochondria_chromosome),
            *get_flag('minimum_count_threshold', minimum_count_threshold)
        ]
    else:
        message(typ='info', data={'title': "Pipeline Launch", "body":" Launching Nextflow with curated genome configuration."})
        cmd = [
            "/root/nextflow",
            "run",
            str(shared_dir / "main.nf"),
            "-work-dir",
            str(shared_dir),
            "-profile",
            profiles,
            "-c",
            "latch.config",
            "-resume",
            "--genome", "custom",
            *get_flag('input_csv', input_csv_samplesheet),
            *get_flag('outdir', outdir),
            *get_flag('genome', genome),
            *get_flag('minimum_count_threshold', minimum_count_threshold)
        ]

    print("Launching Nextflow Runtime")
    print(' '.join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
            "NXF_ENABLE_FS_SYNC": "true",
        }

        if False:
            env["LATCH_LOG_DIR"] = latch_log_dir

        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            remote = LPath(urljoins(latch_log_dir, "nextflow.log"))
            print(f"Uploading .nextflow.log to {remote.path}")
            remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ['du', '-sb', str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print("Failed to compute storage size: Operation timed out after 5 minutes.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)


@workflow(metadata._nextflow_metadata)
def nf_cs_genetics_simplecell_pipeline(input_csv: typing.List[Sample], outdir: LatchDir, star_index: typing.Optional[LatchDir], gtf: typing.Optional[LatchFile], genome: ReferenceGenome = ReferenceGenome.no_selection, mitochondria_chromosome: str = '', minimum_count_threshold: float = 100.0) -> None:
    """
    CS Genetics SimpleCell pipeline

    Sample Description
    """
    local_var_message_str = "\n".join([f"{name}: {value}" for name, value in locals().items()])
    message(typ='info', data={'title': "Local variables", "body":local_var_message_str})

    global_var_message_str = "\n".join([f"{name}: {value}" for name, value in globals().items()])
    message(typ='info', data={'title': "Global variables", "body":global_var_message_str})

    pvc_name: str = initialize()
    # Run validate_genome_selection to check that the user has provided a valid selection of genomic references.
    genome, star_index, gtf, mitochondria_chromosome, single_species = validate_genome_selection(genome=genome, star_index=star_index, gtf=gtf, mitochondria_chromosome=mitochondria_chromosome)
    # Run curate_samplesheet to adjust the input CSV schema based on whether the workflow is single or mixed species.
    input_csv_curated = curate_samplesheet(input_csv=input_csv, single_species=single_species)
    nextflow_runtime(pvc_name=pvc_name, input_csv=input_csv_curated, outdir=outdir, genome=genome, star_index=star_index, gtf=gtf, mitochondria_chromosome=mitochondria_chromosome, minimum_count_threshold=minimum_count_threshold, single_species=single_species)
