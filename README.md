# BTM (Blood Transcription Modules)

Gene modules to interpret blood transcriptomics data.

Used as alternative to conventional pathways, offering granular immunology and often better sensitivity.
The modules can also be used as gene sets for GSEA analysis.
This is the BTM modules described in
https://www.nature.com/articles/ni.2789

Li, S., Rouphael, N., et al, (2014).
Molecular signatures of antibody responses derived from a systems biology study of five human vaccines.
Nature immunology, 15(2), p.195.

## Installation

···
pip install BTM
···

## BTM_Plus at 2021

An update with the help of Amnah Siddiqa.

This set of modules excludes all TBA modules, and replaces the cell specific markers by a new set of markers
compiled from these three papers:

- Aran, D., Hu, Z. and Butte, A.J., 2017. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome biology, 18(1), pp.1-14.

- Zhang, X., Lan, Y., Xu, J., Quan, F., Zhao, E., Deng, C., Luo, T., Xu, L., Liao, G., Yan, M. and Ping, Y., 2019. CellMarker: a manually curated resource of cell markers in human and mouse. Nucleic acids research, 47(D1), pp.D721-D728.

- Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y.Y., Carre, C., Burdin, N., Visan, L., Ceccarelli, M., Poidinger, M. and Zippelius, A., 2019. RNA-Seq signatures normalized by mRNA abundance allow absolute deconvolution of human immune cell types. Cell reports, 26(6), pp.1627-1640.

The cell specific markers were merged by keeping genes in at least half of the source sets. 
I.e., if 5 genesets exist for the same cell population, a gene is kept if it appears in at least 3 of the genesets.
All genesets/modules larger than 100 genes are excluded from BTM_Plus.

```
>>> from BTM_Plus import BTM_Plus as B2
>>> len(B2)
276
>>> B2[88]
{'id': 'P089', 
'name': 'mismatch repair (I)', 
'src': ['Li-Pulendran, Li_M22.0_mismatch repair (I)'], 
'genes': ['SMC1A', 'POLA1', 'NCAPG2', 'RFC5', 'RFC4', 'MSH2', 'TMPO', 'MSH6', 'RFC2', 'GMNN', 'BUB1', 'RMI1', 'RACGAP1', 'EXO1', 'POLD3', 'PRIM1', 'ZWINT', 'CHEK1', 'CENPK', 'FIGNL1', 'MCM6', 'RFC3', 'SSBP1', 'TOPBP1', 'RPA3', 'SMC2']}
```

# To generate permutations of BTM_Plus 

In same directory, change Number_Permutation to your need -
```
python permutate.py
```

This writes two files,
'resampled_BTM_Plus.txt', 
'random_genesets.txt',
each line of a geneset of size 50 (default).

## Example use of the 2013 version (Companion to original paper)

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

The btm_tool.py is to illustrate

1. Converting gene level data to BTM activity table. (Can also convert Affy probeset level data to gene level data.)
2. Enrichment test of an input gene list.
3. Testing antibody correlation to gene expression at module level.

