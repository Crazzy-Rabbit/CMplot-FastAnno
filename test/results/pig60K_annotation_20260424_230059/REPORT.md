# pig60K CMplot_modify Annotation Test Report

- Run id: 20260424_230059
- Data source: F:/Shi/AI_Rewrite/CMplot/CMplot_modify/test/data/pig60K_example.tsv.gz
- Dataset dimensions: 44580 x 6
- Tested trait: trait1
- Highlighted top SNPs: 28
- Best P value: 5.26e-09
- Max -log10(P): 8.279
- Manhattan threshold: 5e-08

## Test Results

- original_rejects_new_argument: PASS (unused argument (highlight.text.mode = "top"))
- modified_top_annotation_png: PASS (file size 59588 bytes; elapsed 1.44 sec)
- modified_line_mode_auto_png: PASS (file size 59588 bytes; elapsed 2.89 sec)
- modified_line_mode_straight_png: PASS (file size 63915 bytes; elapsed 0.95 sec)
- modified_line_mode_elbow_png: PASS (file size 59588 bytes; elapsed 0.93 sec)
- modified_line_mode_none_png: PASS (file size 49286 bytes; elapsed 0.99 sec)
- modified_threshold_default_size_png: PASS (file size 37909 bytes; elapsed 0.97 sec)
- modified_scatter_compatibility_png: PASS (file size 56084 bytes; elapsed 1.33 sec)
- modified_nearby_annotation_png: PASS (file size 56187 bytes; elapsed 0.91 sec)
- modified_multitrack_top_annotation_png: PASS (file size 145106 bytes; elapsed 3.14 sec)

## Top SNP Preview

```
              SNP Chromosome Position       trait1   log10P       label
41953 MARC0066784          8 53910480 5.260000e-09 8.279014 MARC0066784
40400 MARC0040492          8 80539938 6.870000e-08 7.163043 MARC0040492
26779 DRGA0003569          2 67385863 1.080000e-07 6.966576 DRGA0003569
18698 ASGA0039359          8 45703730 3.565820e-07 6.447841 ASGA0039359
36009 INRA0056207         18 34446892 4.483550e-07 6.348378 INRA0056207
13460 ALGA0106996          3 21472898 5.343700e-06 5.272158 ALGA0106996
29178 H3GA0001667          1 43201866 5.805990e-06 5.236124 H3GA0001667
20194 ASGA0054364         12 18258585 6.752110e-06 5.170560 ASGA0054364
14309 ASGA0002752          1 43358326 8.013220e-05 4.096193 ASGA0002752
11681 ALGA0088449         15 95748847 8.536440e-05 4.068723 ALGA0088449
6719  ALGA0049156          8 54271131 9.465000e-04 3.023879 ALGA0049156
30000 H3GA0009305          3 12593770 1.006106e-03 2.997356 H3GA0009305
```

## Output Files

- top annotation default red highlight: Rect_Manhtn.pig60K_trait1_top_annotation.png
- line mode auto: Rect_Manhtn.pig60K_trait1_line_auto.png
- line mode straight: Rect_Manhtn.pig60K_trait1_line_straight.png
- line mode elbow: Rect_Manhtn.pig60K_trait1_line_elbow.png
- line mode none: Rect_Manhtn.pig60K_trait1_line_none.png
- threshold default-size plot: Rect_Manhtn.pig60K_trait1_threshold_default_size.png
- scatter compatibility: Rect_Manhtn.pig60K_trait1_scatter_annotation.png
- nearby annotation: Rect_Manhtn.pig60K_trait1_nearby_annotation.png
- multitrack top annotation: Multi-tracks_Manhtn.pig60K_multitrack_top_annotation.png
- selected top SNPs: pig60K_trait1_top_snps.csv
- annotation target file: pig60K_trait1_annotation_targets.tsv
- raw test results: test_results.csv
- summary: summary.csv
- checksums: checksums.md5
- log: run.log
