BTM (Blood Transcription Modules)
=================================

Gene modules to interpret blood transcriptomics data.

Used as alternative to conventional pathways, offering granular immunology and often better sensitivity.
The modules can also be used as gene sets for GSEA analysis.
This is the BTM modules described in
https://www.nature.com/articles/ni.2789

Li, S., Rouphael, N., et al, (2014).
Molecular signatures of antibody responses derived from a systems biology study of five human vaccines.
Nature immunology, 15(2), p.195.

The btm_tool.py is to illustrate

1. Converting gene level data to BTM activity table. (Can also convert Affy probeset level data to gene level data.)
2. Enrichment test of an input gene list.
3. Testing antibody correlation to gene expression at module level.

Installation:
This program requires Python 2.x, Numpy and Scipy (verion 0.10+, http://scipy.org/install.html).
The optional plotting function depends on Python plot library matplotlib.

## Example use

To convert gene level data to BTM module activity scores:
```
from btm.btm_tool import genetable_to_activityscores
genetable_to_activityscores(infile, outfile)
```

Download tutorial package at
https://media.nature.com/original/nature-assets/ni/journal/v15/n2/extref/ni.2789-S5.zip

This "BTM_tutorial_package" download package should contain -
btm_tool.py, btm_example_data.py, MCV4_D3v0_probesets.txt,
gene_ab_correlation.rnk, BTM_for_GSEA_20131008.gmt, monocytes_vs_bcells.txt.
