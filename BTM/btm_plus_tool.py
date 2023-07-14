‘’‘
## BTM_Plus at 2021

An update with the help of Amnah Siddiqa.

This set of modules excludes all TBA modules, and replaces the cell specific markers by a new set of markers
compiled from these three papers:

- Aran, D., Hu, Z. and Butte, A.J., 2017. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome biology, 18(1), pp.1-14.

- Zhang, X., Lan, Y., Xu, J., Quan, F., Zhao, E., Deng, C., Luo, T., Xu, L., Liao, G., Yan, M. and Ping, Y., 2019. CellMarker: a manually curated resource of cell markers in human and mouse. Nucleic acids research, 47(D1), pp.D721-D728.

- Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y.Y., Carre, C., Burdin, N., Visan, L., Ceccarelli, M., Poidinger, M. and Zippelius, A., 2019. RNA-Seq signatures normalized by mRNA abundance allow absolute deconvolution of human immune cell types. Cell reports, 26(6), pp.1627-1640.
’‘’

import numpy as np
from scipy import stats

def dequote(s): return s.replace('"', '').replace("'", "")

def read_gene_data(infile, start_col=1, sep='\t'):
    '''
    returns gene_data dictionary and header
    '''
    gene_data = {}
    w = open(infile).readlines()
    header = w[0].rstrip().split(sep)
    for line in w[1:]:
        a = line.rstrip().split(sep)
        if len(a) > 1:
            data = np.array([float(x) for x in a[1:]])
            # using dequote in case user mistakes file format - 
            # some spreadsheet program may add quotes automatically
            gene_data[dequote(a[0]).upper()] = data

    print("Got %d lines of gene data: " %len(gene_data))

    return gene_data, header

def gene_data_to_activityscores(gene_data,
                                header, 
                                module_list,
                                zscore=True,
                                outfile='btm_plus_converted_activities.tsv'):
    '''
    Convert gene level data to BTM activity scores.
    
    Note: this works without zscore when difference of expression is taken between time points,
    which self-normalize to subjects. 
    Standardization is needed if the activity score is computed for original gene expression,
    e.g. z-scores are calculated across samples and module scores are then computed on z-scores. 
    Default for gene_data_to_activityscores.
    '''
    if zscore:
        for g,v in gene_data.items():
            v = stats.zscore(v)
        
    s = '\t'.join(header) + '\n'
    for x in module_list:
        ascores = compute_activity_score(x['genes'], gene_data)
        s += x + '\t' + '\t'.join([str(d.round(3)) for d in ascores]) + '\n'
        
    out = open(outfile, 'w')
    out.write(s)
    out.close()

def compute_activity_score(M, datadict):
    '''
    compute module M activity score as mean expression value of member genes, 
    using numpy.array function.
    datadict is formatted as {gene: expression_vector}.
    '''
    data = []
    for x in M:
        try:
            data.append(datadict[x])
        except KeyError:
            pass
            
    N = len(data)
    return np.array(data).sum(0)/N
    