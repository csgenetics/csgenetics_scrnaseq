"""Static fail-loud contracts for published coordinate-sorted BAM producers."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_initial_featurecounts_bam_is_sorted_after_annotation():
    source = (ROOT / "modules/local/initial_feature_count/main.nf").read_text(
        encoding="utf-8"
    )

    featurecounts = source.index("featureCounts -a")
    coordinate_sort = source.index("samtools sort -@", featurecounts)
    final_move = source.index("mv ${sample_id}.featureCounts.coordinate.bam")

    assert featurecounts < coordinate_sort < final_move
    assert "${sample_id}_Aligned.sortedByCoord.out.bam.featureCounts.bam" in source[
        coordinate_sort:final_move
    ]


def test_final_annotated_bam_is_sorted_while_merging():
    source = (
        ROOT
        / "modules/local/merge_annotated_UMRs_with_annotated_multimappers/main.nf"
    ).read_text(encoding="utf-8")

    merge = source.index("samtools merge -u")
    coordinate_sort = source.index("| samtools sort", merge)
    output = source.index(
        "-o ${sample_id}.mapped.sorted.filtered.annotated.bam", coordinate_sort
    )

    assert merge < coordinate_sort < output
