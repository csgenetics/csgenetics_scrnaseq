# 06 - Test data inventory & local run feasibility (EXTERNAL vs INTERNAL)

## How EXTERNAL test profiles are wired

`nextflow.config:93-146` defines five test-related profiles. Each `includeConfig`s
`conf/base.config` + a per-profile `conf/test*.config` + `conf/images.config`, and
sets the container engine.

| Profile | Container engine | Config file | nextflow.config lines |
|---|---|---|---|
| `test` | docker | `conf/test.config` | 93-103 |
| `test_pbmc_4_sample_full` | docker | `conf/test_pbmc_4_sample_full.config` | 104-114 |
| `test_hsap_mmus_2_sample_full` | docker | `conf/test_hsap_mmus_2_sample_full.config` | 115-125 |
| `test_singularity` | singularity | `conf/test.config` (+ genomes.config) | 126-136 |
| `test_conda` | conda | `conf/test.config` (+ genomes.config, conda_envs.config) | 137-146 |

`test_singularity` and `test_conda` reuse the SAME tiny `test` dataset; they only swap
the container backend. Run command pattern: `nextflow run main.nf -profile <name>`.

Barcode list for the `test` profile: `s3://csgx.public.readonly/resources/barcode_lists/IDT_IO_kit_v2.csv` (1.3 MB). The two full profiles inherit barcode_list_path from `nextflow.config` defaults (not overridden in their config).

## Per-profile detail (FASTQ source, reference, size, local runnability)

All S3 sources below are in `s3://csgx.public.readonly/` (public-read bucket) and STAR
indices in `s3://csgx.public.readonly/resources/references/`. AWS access for the
`scrnaseq-agent` IAM user is confirmed working (see "AWS access" below); both buckets read fine.

### 1. `test` (the smoke test) — `conf/test.config`
- **input_csv:** `s3://.../test_profile_resources/test_input.csv` (2 samples, Sample1/Sample2). A LOCAL copy exists at `input_csv/test_input.csv` (identical paths).
- **FASTQs:** `s3://.../test_profile_resources/reduced/Sample{1,2}_L001_R{1,2}_001.sub.chr1.16th.1M.fastq.gz`. ~1M reads/sample, chr1-subset. Total ~99 MB (4 files: 27/21/29/22 MB).
- **Reference (genome `GRCh38_test`):** reduced STAR index `.../test_profile_resources/reduced/STAR_2_7_11b/GRCh38.Ensembl109.GENCODEv44.test_profile/star/` + matching reduced GTF (`...gene_subset_lof.reduced.gtf`). STAR index = **~4.5 GB** (26 objects, 4,460,841,460 bytes).
- **Expected output:** ~309 + 563 cells, ~60 genes/cell (from config description).
- **mixed_species:** false; mito chr = MT.
- **Runnable locally NOW:** YES. Smallest by far. Only needs ~4.6 GB download + docker. Best smoke test.

### 2. `test_pbmc_4_sample_full` — `conf/test_pbmc_4_sample_full.config`
- **input_csv:** `s3://.../test_profile_resources/HIVE1802_PBMC_4_sample_full/input_csv/HIVE1802_PBMC_4_sample_full.input.csv` (4 PBMC samples: S275,273,277,280).
- **FASTQs:** `.../HIVE1802_PBMC_4_sample_full/fastq/241023_Dorothy_HIVE1802_1400pM_Sample*_L001_R{1,2}_001.fastq.gz`. 54-66M reads/sample. Total fastq = **~12.9 GB** (9 objects).
- **Reference (genome `GRCh38_pbmc_test`):** FULL production GRCh38 STAR index `.../references/Ensembl_Gencode_resources/STAR_2_7_11b/GRCh38.Ensembl109.GENCODEv44_with_pcLoF_without_readthrough/star/` = **~29.7 GB** (16 objects, 29,703,629,625 bytes) + matching full GTF.
- **mixed_species:** false.
- **Runnable locally:** Heavy. ~13 GB FASTQ + ~30 GB STAR index + 4-sample compute. Functional but slow; not a smoke test.

### 3. `test_hsap_mmus_2_sample_full` — `conf/test_hsap_mmus_2_sample_full.config`
- **input_csv:** `s3://.../test_profile_resources/HIVE1809_mixed_2_sample_full/input_csv/HIVE1809_Hsap_Mmus_2_sample_full.input.csv` (2 mixed-species samples).
- **FASTQs:** `.../HIVE1809_mixed_2_sample_full/fastq/`. 99M & 133M reads. Total fastq = **~12.0 GB** (5 objects).
- **Reference (genome `mouse_human_mix_test`):** FULL combined human+mouse STAR index `.../STAR_2_7_11b/GRCh38.Ensembl109.GENCODEv44_GRCm39.Ensembl110.GENCODEvM33_with_pcLoF_without_readthrough/star/` = **~56.4 GB** (16 objects, 56,362,404,131 bytes) + species-tagged GTF.
- **mixed_species:** TRUE. Exercises the mixed-species/barnyard code path (hsap/mmus prefixes GRCh38_/GRCm39_, separate mito chrs GRCh38_MT/GRCm39_MT). This is the ONLY external test profile that covers mixed-species.
- **Runnable locally:** Heaviest. ~12 GB FASTQ + **~56 GB** STAR index. Functional but the slowest/largest.

## Other test-data artifacts in EXTERNAL repo (not wired to a profile)

- `input_csv/` (local helper CSVs):
  - `test_input.csv` — local copy of the `test` profile CSV (S3 paths).
  - `template.csv`, `manual_threshold_template.csv`, `manual_threshold_mixed_species_template.csv` — user templates with placeholder `/home/example_user/...` paths (NOT runnable).
  - `badsamples.csv` — intentionally-bad samples (likely for failure-mode testing).
- `mixed_testing.csv` (repo root): 2 mixed samples (Sample289_S289, Sample129_S129) with FASTQs at **local `/nssd2/humebc/v4/mixed_v2_testing_reads/...`** (existence not verified; ad-hoc dev CSV, not a profile).
- `MOR027_test_data/MOR027.brian.external.lite.input.metadata.csv`: 3 samples (Sample1/2/3_S*) pointing at **`s3://csg-nextflow/251212_UCSD_MOR027_1_1400pM/...`**. These are LARGE raw FASTQs (~1.2 GB R1 + ~1.0 GB R2 per sample observed) on the private `csg-nextflow` bucket (read access confirmed via scrnaseq-agent). Ad-hoc dataset, not a profile; "lite" in name but files are full-size. Folder appears recently added (Jan 2026, untracked dir per git status).
- `images/`: NOT test data — these are Docker container build contexts (aws_download, Dockerized, gtf2bed, html_build, io_extract, multiqc, Report, samtools_star, tidyverse, umi-tools-csgx). Referenced by `conf/images.config`.

## AWS access (confirmed)
- External creds: `/nssd2/humebc/agent-secrets/scrnaseq/aws/{credentials,config}`.
  Export `AWS_SHARED_CREDENTIALS_FILE` + `AWS_CONFIG_FILE` to these.
  `sts get-caller-identity` => `arn:aws:iam::837012689225:user/scrnaseq-agent`.
- Reads OK for `s3://csgx.public.readonly/...` (all 3 profiles' data) and `s3://csg-nextflow/...` (MOR027). Same account 837012689225 as internal rnaseq-agent.

## What INTERNAL (reference) pipeline uses for test data

Internal uses different data layout AND different CSV schema (long-form metadata columns:
`sample_id,cell_type,num_input_cell,...,s3_path_R1,s3_path_R2`), NOT external's 3-column
`sample,fastq_1,fastq_2`. So internal CSVs are NOT directly reusable by external.

- **Skill `run-pipeline`** (`.claude/skills/run-pipeline/references/run-types.yml`): 4 curated medium-weight run types, each 4 samples on a 2x2 design, run with `-profile csgxserver`:
  - `mixed-species-kol0054`, `human-cellline-kol0054`, `mouse-cellline-kol0054` (same KOL0054 FASTQs vs 3 references), `pbmc-mor034`.
  - FASTQs live LOCALLY: `/nssd2/humebc/pipeline-runs/fastqs/{KOL0054,MOR034}/`.
  - References live LOCALLY: STAR/GTF at `/10TB_B/resources/...` via `conf/species_profiles_local.config` (profiles: GRCh38, GRCh38_v110, GRCm39, GRCm39_with_readthrough_no_LoF, mouse_human_mix, mouse_human_mix_v110). `csgxserver` profile caps STAR maxForks=3 for the 256 GB box.
- **`input_csvs/` lite/CI datasets** (the small ones, S3-backed on `s3://csg-reference/circleci-reference-fastqs_v3/...`):
  - `human.100000.test.input_csv.csv` (2 samples, 100k-read human-reduced)
  - `mixed.100000.test.input_csv.csv` (2 samples, 100k-read mixed; reuses the same Sample129/131 HIVE1326 reads family as external's `mixed_testing.csv`)
  - `mouse.100000.test.input_csv.csv`
  - `test_medium_{human_cellline,mixed_species,mouse_cellline,pbmc}.csv` (the 4 run-types above)
  - `circle_ci_input.data_provenance.single_lane.csv`
- Internal nf-test fixtures: `s3://csg-reference/internal_nf_tests_data/*` (agent has read/write).

### Reusability for EXTERNAL
- **Raw FASTQs are reusable conceptually**, but the CSV schema differs — would need conversion to external's 3-column format. The `circleci-reference-fastqs_v3` 100k-read reduced FASTQs (human/mixed/mouse) on `s3://csg-reference/` are the smallest real-data fixtures available and would make excellent fast external smoke tests IF rewrapped in a 3-column external CSV and paired with a small external-compatible reference.
- **References are NOT directly reusable:** internal's `/10TB_B/resources/` local indices and the full `s3://csgx.public.readonly/.../references/` indices are full-size (~30-56 GB). External's own reduced `test_profile` STAR index (~4.5 GB) is the only small reference and is human-only chr1-subset.
- The mixed 100k FASTQs (Sample129/131, HIVE1326) overlap with external's `mixed_testing.csv` Sample129 — same underlying source library, so behaviourally comparable.

## RECOMMENDATION

1. **Fastest smoke test to run locally first: `-profile test`.** It is by far the lightest:
   ~99 MB FASTQ + ~4.5 GB reduced STAR index, 2 samples, ~1M reads each, docker.
   Everything is on the public `csgx.public.readonly` bucket and AWS access is confirmed.
   Use this to validate the nf26 migration end-to-end before touching the heavy profiles.
   (`test_singularity`/`test_conda` reuse the same data if backend coverage is also wanted.)

2. **Keep `test_pbmc_4_sample_full` and `test_hsap_mmus_2_sample_full` as the heavier
   functional/regression profiles** (only the mixed one exercises the barnyard path).
   They are runnable on beast but are 12-13 GB FASTQ + 30-56 GB references each.

3. **No new tiny dataset is strictly required** — `-profile test` already serves as the
   fast smoke test. However, two gaps worth closing during the overhaul:
   - There is no *small/fast* MIXED-species smoke test (the only mixed profile pulls a 56 GB
     index). The internal `mixed.100000` reduced FASTQs + a reduced mixed STAR index could be
     adapted into a fast external mixed-species smoke profile. This would need (a) rewrapping
     the FASTQs in a 3-column external CSV and (b) building/locating a reduced mixed STAR index
     (external currently has no reduced mixed reference). Recommend generating one (mirror how
     the existing human `test_profile` reduced index was built) if mixed-species coverage in CI
     is desired.
   - The `MOR027_test_data/` and `mixed_testing.csv` are ad-hoc, untracked/dev-only; do not
     rely on them as canonical test profiles.
