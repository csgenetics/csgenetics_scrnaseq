# Count-statistics arithmetic

The pipeline treats every entry in a raw count matrix as a non-negative integer
observation. H5AD stores the matrix as `float32` for compatibility with the R
and Seurat readers supported by the pipeline, but that storage choice does not
make counts continuous measurements.

Before calculating cell statistics or categorising counts as inside/outside
cells, the pipeline therefore checks that all stored values are real, finite,
non-negative integers in the exact-integer range of their floating data type.
It also checks that the total can be reduced without overflowing a signed
64-bit integer. Invalid matrices fail rather than producing rounded or wrapped
metrics. Sparse duplicate entries and explicit zeros are canonicalised on a
copy; the H5AD matrix itself is not changed.

All count totals, detection counts, and integer medians are then calculated
with sparse `int64` arithmetic. The calculation scales with stored entries plus
the number of cells and genes; it does not create a dense cell-by-gene matrix.

For mixed-species samples, per-cell metrics continue to use only the genes for
the cell's assigned species. The two combined across-sample detection metrics
continue to count every detected gene in every single cell, including an
off-species detection. If one species has no cells, its metrics are numeric
zero.

This correction can change the low-order decimal spelling of a raw
`*.metrics.csv` mean or percentage that was previously accumulated in
`float32` (for example, a shortened float approximation becomes the exact
integer-ratio result). Reports remain formatted to two decimal places. For the
representative validation fixtures the correction was below that display
precision; a valid sufficiently high-count matrix can, however, change the
last displayed decimal digits as well as the raw CSV spelling. That is an
intentional arithmetic correction rather than a change to a metric's
definition. Integer count totals and the H5AD/tripartite matrix schemas are
unchanged.
