# CMplot-FastAnno Compatibility Report

- Run id: 20260424_230126
- Data source: F:/Shi/AI_Rewrite/CMplot/CMplot_modify/test/data/pig60K_example.tsv.gz
- Dataset dimensions: 44580 x 6
- Scope: smoke/regression coverage for density, Manhattan, circular Manhattan, QQ, highlight text, multi-track, and multi-trait calls.
- `BASELINE_FAIL` means the original local CMplot failed that call; it is reported for context and does not block CMplot-FastAnno validation.
- Visual differences are expected only for explicitly requested CMplot-FastAnno behavior, such as GWASLab-style colors, chromosome x-axis title/labels, tighter y-axis padding, and threshold amplification default.

## Results

- original / density: PASS (0.15 sec; Marker_Density.density.png)
- modified / density: PASS (0.24 sec; Marker_Density.density.png)
- original / manhattan_single: PASS (2.27 sec; Rect_Manhtn.manhattan_single.png)
- modified / manhattan_single: PASS (2.45 sec; Rect_Manhtn.manhattan_single.png)
- original / qq_single: PASS (0.39 sec; QQplot.qq_single.png)
- modified / qq_single: PASS (0.31 sec; QQplot.qq_single.png)
- original / circular_single: PASS (0.58 sec; Cir_Manhtn.circular_single.png)
- modified / circular_single: PASS (0.52 sec; Cir_Manhtn.circular_single.png)
- original / legacy_highlight_text: PASS (0.69 sec; Rect_Manhtn.legacy_highlight_text.png)
- modified / legacy_highlight_text: PASS (0.75 sec; Rect_Manhtn.legacy_highlight_text.png)
- original / multitrack_manhattan: PASS (1.72 sec; Multi-tracks_Manhtn.multitrack_manhattan.png)
- modified / multitrack_manhattan: PASS (1.74 sec; Multi-tracks_Manhtn.multitrack_manhattan.png)
- original / multitrait_manhattan: BASELINE_FAIL (NA sec; )
- modified / multitrait_manhattan: PASS (3.83 sec; Multi-traits_Manhtn.multitrait_manhattan.png)
- original / multitrack_qq: PASS (0.81 sec; Multi-tracks_QQplot.multitrack_qq.png)
- modified / multitrack_qq: PASS (0.78 sec; Multi-tracks_QQplot.multitrack_qq.png)

## Runtime Summary

- modified: 10.62 sec total
- original: 6.61 sec total

## Files

- detailed CSV: compatibility_results.csv
- summary CSV: compatibility_summary.csv
- original output folder: original
- modified output folder: modified
- log: run.log
