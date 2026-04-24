# Test Data

This folder contains small, stable inputs for testing `CMplot-FastAnno`.

## Files

- `pig60K_example.tsv.gz`: main CMplot input table. Columns are `SNP`, `Chromosome`, `Position`, `trait1`, `trait2`, and `trait3`.
- `pig60K_trait1_annotation_targets.tsv`: annotation table for `trait1`. Columns are `SNP`, `Label`, and `Trait`.
- `pig60K_trait1_top_snps.csv`: selected top SNPs used by the example tests.

## Run

```r
Rscript test/test_pig60k_annotation.R
Rscript test/test_cmplot_fastanno_compatibility.R
```

To recreate these files from `data("pig60K", package = "CMplot")`:

```r
Rscript test/create_test_data.R
```
