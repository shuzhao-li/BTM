"""
To generate permuation lists of genesets of fixed sizes.

Tier 1: resembling the real BTM_Plus, but in size of ~50 genes.
Tier 2: shuffling real BTM_Plus in slices.
Tier 3: random draw of genes not in BTM_Plus

>>> B2[0]
{'src': ['Xcell_naive B-cells_BLUEPRINT_2', 'Xcell_naive B-cells_HPCA_1', 'Xcell_naive B-cells_HPCA_2', 'Xcell_naive B-cells_HPCA_3', 'Xcell_naive B-cells_NOVERSHTERN_1', 'Xcell_naive B-cells_NOVERSHTERN_2', 'Xcell_naive B-cells_NOVERSHTERN_3', 'Monaco_B Naive'], 
'genes': ['PRDM4', 'COL19A1', 'CD22', 'P2RY10', 'MS4A1', 'MBD4', 'KHDRBS2', 'TCL1A', 'FCER2', 'SMC6', 'CAPN3', 'TSPAN13', 'FCRL2', 'CD19', 'TCL6', 'STAG3', 'CD72', 'GGA2', 'CXCR5', 'RRAS2', 'SIPA1L3', 'TCL1B', 'DSP', 'CSNK1G3', 'MMP17'], 
'id': 'P001', 
'name': 'naive B-cells'}
"""

SIZE = 50                   # geneset size
BinSize = 5                 # bin of genes, for efficiency
Number_Permutation = 1000   # counts for tier 3 random sets 

from random import sample
from BTM_Plus import BTM_Plus as B2
from btm_example_data import affy_probeset_dict

B2Genes = []
for m in B2:
    B2Genes += m['genes']

B2Genes = set(B2Genes)
allGenes = set(affy_probeset_dict.values())
otherGenes = allGenes - B2Genes
print(len(B2Genes), len(otherGenes), len(allGenes))

def slice_geneset(geneset, BinSize = BinSize):
    '''
    return slices of geneset in bin of 5 genes. The last bin may be smaller.
    Bins are sets, as many operations are on sets.
    '''
    steps = int( len(geneset)/BinSize )
    slices = []
    for ii in range(steps):
        slices.append( set(geneset[5*ii: 5*(ii+1)]) )

    return slices

# all sliced bins of BTM_Plus
geneBins = []
for m in B2: geneBins += slice_geneset(m['genes'])
print("geneBins: ", len(geneBins))


############################################################
# 
############################################################

def draw_random_otherGenes(otherGenes = otherGenes, SIZE = SIZE, Number_Permutation = 1000):
    '''
    Tier3.
    Draw random genes in group of SIZE.
    otherGenes needs to be a set, otherwise Python3 returns groups in lexicographic order.

    Different from Tier1 and Tier2, this returns List of lists (not sets).
    '''
    N = 0
    tmp = []
    while N < Number_Permutation:
        tmp.append(sample( otherGenes, SIZE ))
        N += 1

    return tmp

def combine_fit_bins(seed, geneBins, SIZE=SIZE):
    '''
    Input
    =====
    seed: a set of genes. seed must be smaller than SIZE.
    geneBins: a list of sets.

    Output size can be trimmed to precision but only unique genes count.

    Return
    ======
    Target size of geneset is SIZE (default 50)
    '''
    # 2 x number of genes that are needed, in case too many repeats - there'd always be overlap btw sets
    needed = SIZE - len(seed)
    missing = 2 * SIZE/BinSize
    missing = sample(geneBins, missing)
    seed2 = set(seed)
    while len(seed2) < SIZE:
        for ss in missing:
            seed2 = seed2.union(ss)

    # trim
    added_genes = list(seed2 - seed)[:needed]
    return seed.union(set(added_genes))

def seed_tier1_genesets(BTMs=B2, geneBins=geneBins, SIZE=SIZE):
    '''
    To create a list of genesets most resembling the real BTM_Plus, 
    but in size of SIZE (default 50) genes (Tier1).
    If an input geneset is larger than 50, split and fit the 2nd part.

    Return
    ======
    List of sets, each is a geneset of SIZE 
    '''
    tmp = []
    for m in BTMs:
        print(m['name'])
        if len(m['genes']) > SIZE:
            tmp.append( set(m['genes'][:SIZE]) )
            tmp.append( combine_fit_bins(set(m['genes'][SIZE:]), geneBins, SIZE=SIZE) )
        else:
            tmp.append( combine_fit_bins(set(m['genes']), geneBins, SIZE=SIZE) )

    return tmp

def merge_sets(ListOfSets):
    L = set()
    for ss in ListOfSets:
        L.union(ss)
    return L

def make_tier2_genesets(geneBins=geneBins, SIZE=SIZE, Number_Permutation = 10000):
    '''
    Return
    ======
    List of sets, each is a geneset of SIZE 
    '''
    permutationList = []
    min_num = SIZE/BinSize
    for ii in xrange(Number_Permutation):
        permutationList.append( combine_fit_bins(
                                    merge_sets(sample(geneBins, min_num)), 
                                    geneBins, 
                                    SIZE=SIZE) 
                                )
    return permutationList


# -------------------------------------------------------------------
# main - change Number_Permutation to your need
#

if __name__ == '__main__':
    tier1 = seed_tier1_genesets()
    print("Done tier1.")
    tier2 = make_tier2_genesets(Number_Permutation = 1000)
    print("Done tier2.")

    # export each line a geneset
    s1 = ''
    for x in tier1: s1 += ','.join(x) + '\n'
    s2 = ''
    for x in tier2: s2 += ','.join(x) + '\n'
    with open('resampled_BTM_Plus.txt', 'w') as O:
        O.write(s1+s2)


    tier3 = draw_random_otherGenes(Number_Permutation = 1000)
    print("Done tier3.")
    s3 = ''
    for x in tier3: s3 += ','.join(x) + '\n'
    with open('random_genesets.txt', 'w') as O:
        O.write(s3)

