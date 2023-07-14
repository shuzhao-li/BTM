# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# Copyright (c) 2013 Shuzhao Li.
# This is intended to demonstrate the application of BTM modules.

'''
Features:
1) Converting gene level data to BTM activity table. 
   (Can also convert Affy probeset level data to gene level data.)
2) Enchriment test of an input gene list.
3) Testing antibody correlation to gene expression at module level.

Installation:
This program requires Python 2.x, Numpy and Scipy (verion 0.10+, http://scipy.org/install.html).
The optional plotting function depends on Python plot library matplotlib.

This "BTM_tutorial_package" download package should contain -
btm_tool.py, btm_example_data.py, MCV4_D3v0_probesets.txt,
gene_ab_correlation.rnk, BTM_for_GSEA_20131008.gmt, monocytes_vs_bcells.txt.

Usage Example:

1. Evoke Python in the same directory where you unpacked the downloaded files, which should include
btm_tool.py and others. In Python interpreter environment, do
>>> from btm_tool import *
This makes BTM data and a few functions available for your use.

2. The supplied file, MCV4_D3v0_probesets.txt, is Affymetrix probeset level data from the MCV4
study. We can convert it to gene level data by
>>> probeset_to_genetable('MCV4_D3v0_probesets.txt', affy_probeset_dict, 'MCV4_D3v0_genes.txt')
The output file is MCV4_D3v0_genes.txt. This program supports only Affymetrix platform. You will
need to prepare your own gene level data for the next step.

3. To convert gene level data to BTM module activity scores:
>>> genetable_to_activityscores('MCV4_D3v0_genes.txt', 'MCV4_D3v0_BTMactivity.txt')
The output file is MCV4_D3v0_BTMactivity.txt. The module activity scores are computed as the mean
value of member genes. You can use these activity scores to perform further statistical test of your
choice.

4. A common bioinformatics task is to test the over-representation in a list of genes. To do this with
BTM modules, we first need to get a gene list of interest.
>>> genelist = [x.split('\t')[0] for x in open('gene_ab_correlation.rnk').readlines()[20: 220]]
This is a trick to pull 200 genes from the gene_ab_correlation.rnk file we used earlier. You may
want to get a more interesting genelist from your own data.

5. To do enrichment test on this gene list using BTM:
>>> enrichment_test(genelist, 'my_enrichment_test.txt')
This performs Fisher Exact Test on this gene list and each BTM module. The output file,
my_enrichment_test.txt, contains enrichment p-values of each module in tab-delimited format. You
can import it into a spreadsheet program for further formatting and editing.

6. To do antibody correlation analysis using BTM framework:
>>> do_antibody_correlation('MCV4_D3v0_genes.txt', mcv4_log2_antibody, 'my_correlation_test.txt')
This takes some time because the correlation significance is estimated by permutations of both samplelabels and gene memberships. The input data are gene level expression file MCV4_D3v0_genes.txt and
mcv4_log2_antibody, which is pre-loaded as an example. The antibody data, in matched sample
order, will of course have to be supplied by users in a real analysis. The output file,
my_correlation_test.txt, contains p-values of each module in tab-delimited format. You can import
it into a spreadsheet program for further formatting and editing. If the plot library matplotlib is
installed, a probability distribution figure will also be produced by this step.
This btm_tool program is provided as demonstration code. Users should feel free to modify and
incorporate it into their own analysis.
'''

import random as rd
from scipy import stats, std
from numpy import array, random

# check availability of plotting library matplotlib
HAS_matplotlib = True
try:
    from matplotlib import pyplot
except ImportError:
    HAS_matplotlib = False


#
# btm_example_data.py
# includes affy probeset -> gene dictionary; 
# BTM modules ModuleIndex, ModuleDict; 
# example antibody data (log2 MCV4 anti-PS IgG D30/0)
#
from btm_example_data import *


#
# a few common functions
#

def cmpr_data(x1, x2):
    if abs(sum(x1)) > abs(sum(x2)):
        return x1
    else:
        return x2

def dequote(s): return s.replace('"', '').replace("'", "")

def get_url(module_id):
    htmlpath = 'http://www.interactivefigures.com/meni/btm416_annotation/btmdata/'
    return htmlpath + module_id + '.htm'

def compute_activity_score(M, datadict):
    '''
    compute module M activity score as mean expression value of member genes, 
    using numpy.array function.
    datadict is formatted as {gene: expression_vector}.

    Note: this works when difference of expression is taken between time points,
    which self-normalize to subjects. 
    Standardization is needed if the activity score is computed for original gene expression,
    e.g. z-scores are calculated across samples and module scores are then computed on z-scores. 
    '''
    data = []
    for x in M:
        try:
            data.append(datadict[x])
        except KeyError:
            pass
            
    N = len(data)
    return array(data).sum(0)/N

def set_flatten_list(L):
    new = []
    for x in L: new += x
    return set(new)

#
# data conversion functions
#

def probeset_to_genetable(infile, probeset_dict, outfile=''):
    '''
    Convert probeset level data to gene level data.
    Input file has 1st row as header, 1st col probeset ID.
    User can specify output file name.
    Unmatched probeset IDs are ignored.
    '''
    genedict = {}
    if not outfile:
        outfile = 'genetable_' + infile
    w = open(infile).readlines()
    for line in w[1:]:
        a = line.rstrip().split('\t')
        if len(a) > 1:
            data = [float(x) for x in a[1:]]
            g = probeset_dict.get(dequote(a[0]), '')
            if g:
                if genedict.has_key(g):
                    genedict[g] = cmpr_data( data, genedict[g] )
                else:
                    genedict[g] = data
                
    genes = genedict.keys()
    genes.sort()
    s = w[0]
    for g in genes: s += g + '\t' + '\t'.join( [str(x) for x in genedict[g]] ) + '\n'
    
    out = open(outfile, 'w')
    out.write(s)
    out.close()


def genetable_to_activityscores(infile, outfile=''):
    '''
    Convert gene level data to BTM activity scores.
    Input file has 1st row as header, 1st col gene symbol.
    User can specify output file name.
    '''
    global ModuleIndex, ModuleDict, compute_activity_score
    genedict = {}
    if not outfile:
        outfile = 'activityscore_' + infile
    w = open(infile).readlines()
    for line in w[1:]:
        a = line.rstrip().split('\t')
        if len(a) > 1:
            data = [float(x) for x in a[1:]]
            # using dequote in case user mistakes file format - 
            # some spreadsheet program may add quotes automatically
            genedict[dequote(a[0]).upper()] = data
        
    s = w[0]
    for x in ModuleIndex:
        ascores = compute_activity_score(ModuleDict[x], genedict)
        s += x + '\t' + '\t'.join([str(d.round(3)) for d in ascores]) + '\n'
        
    out = open(outfile, 'w')
    out.write(s)
    out.close()


#
# enrichment test of a gene list
#

def enrichment_test(genelist, outfile=''):
    '''
    Enrichment test (Fisher Exact Test) on an input gene list, 
    using BTM modules similarly to canonical pathways.
    User can specify output file name.
    Fisher exact test is using scipy.stats.fisher_exact
    for right-tail p-value:
    >>> stats.fisher_exact([[12, 5], [29, 2]], 'greater')[1]
    0.99452520602188932
    
    return foramt:  (Module, enrich_pvalue, module_size, overlap_size, overlap_features) 
    '''
    global ModuleDict
    if not outfile:
        outfile = 'btm_enrich_' + infile
    
    allgenes = set_flatten_list(ModuleDict.values())
    total_feature_num = len(allgenes)
    qset = set([x.upper() for x in genelist]).intersection( allgenes )
    query_set_size = len(qset)
    print ("Input list contains %d genes out of total %d genes in BTMs." 
            %(query_set_size, total_feature_num))
    
    result = []
    for modulename, members in ModuleDict.items():
        module_size = len(members)
        overlap_features = qset.intersection(members)
        overlap_size = len(overlap_features)
        
        negneg = total_feature_num + overlap_size - module_size - query_set_size
        # Fisher's exact test
        p_FET = stats.fisher_exact([[overlap_size, query_set_size - overlap_size],
                               [module_size - overlap_size, negneg]], 'greater')[1]
        result.append( (p_FET, modulename, module_size, overlap_size, overlap_features) )
    
    result.sort()
    s = 'Module\tEnrichment_p_value\tModule_size\tOverlap_size\tOverlap_genes\n'
    for line in result:
        s += '\t'.join(
            [line[1], str(line[0]), str(line[2]), str(line[3]), ','.join(line[4])]
        ) + '\n'
    
    out = open(outfile, 'w')
    out.write(s)
    out.close()


#
# antibody correlation analysis
#

def make_moduledict(mdict, datadict):
    '''
    Compute activity scores for all modules in mdict. 
    datadict is formatted as {gene: expression_vector}.
    Return {module: activity_score}.
    '''
    outdict = {}
    for k,v in mdict.items():
        outdict[k] = compute_activity_score(v, datadict)

    return outdict

def activity_score_scrambled(M, datadict):
    '''
    Compute scrambled activity score of a single module M by shuffling samples.
    datadict is formatted as {gene: expression_vector}.
    '''
    data = []
    for x in M:
        try:
            t = datadict[x]
            rd.shuffle(t)
            data.append(t)
        except KeyError:
            pass
            
    N = len(data)
    return array(data).sum(0)/N

def make_scrambled_modules():
    '''
    Make random modules by shuffling gene members.
    '''
    global ModuleDict, ModuleIndex
    scrambled_ModuleDict = {}
    allgenes = set_flatten_list(ModuleDict.values())
    module_sizes = [len(ModuleDict[x]) for x in ModuleIndex]
    for ii in range(len(module_sizes)):
        m = 'random_' + str(ii)
        #scrambled_ModuleIndex.append( m )
        scrambled_ModuleDict[m] = rd.sample(allgenes, module_sizes[ii])

    return scrambled_ModuleDict

def make_random_moduledict(datadict):
    '''
    Generate permutation data by 100x num of random modules.
    datadict is formatted as {gene: expression_vector}.
    '''
    outdict = {}
    for ii in range(100):
        scrambled_ModuleDict = make_scrambled_modules()
        for k,v in scrambled_ModuleDict.items():
            outdict[rd.random()] = activity_score_scrambled(v, datadict)

    return outdict


def get_Pearson_rlist(datalist, ab):
    '''
    return a list of Pearson r values, between each row in datalist and ab.
    '''
    return [stats.pearsonr(x, ab)[0] for x in datalist]



# main function of antibody correlation analysis

def do_antibody_correlation(infile, antibody, outfile=''):
    '''
    Antibody correlation analysis using BTM modules.
    This function computes the Pearson correlation between module activity and antibody data,
    then estimate the p-values based on 100x permutation of sample labels and module memberships.
    
    Input:    
    infile is gene table data, 1st row as headers and 1st col gene symbols.
    antibody is a vector, which must be matched to sample order in gene table.
    
    Output:
    A report of Pearson correlation coefficients and p-values is written in tab delimited format.
    If matplotlib library is installed, a probability distribution figure is also produced.
    '''
    genedict = {}
    if not outfile:
        outfile = 'btm_abcorr_' + infile
    w = open(infile).readlines()
    for line in w[1:]:
        a = line.rstrip().split('\t')
        if len(a) > 1:
            data = [float(x) for x in a[1:]]
            # using dequote in case user mistakes file format - 
            # some spreadsheet program may add quotes automatically
            genedict[dequote(a[0]).upper()] = data

    real_scoredict = make_moduledict(ModuleDict, genedict)
    random_scoredict = make_random_moduledict(genedict)
    
    m = get_Pearson_rlist([ real_scoredict[x] for x in ModuleIndex ], antibody)
    msr = get_Pearson_rlist(random_scoredict.values(), antibody)
    
    # rank significance
    m2 = [abs(x) for x in m]
    total = m2 + [abs(x) for x in msr]
    total.sort(reverse=True)
    D = float(len(total))
    p_values = [(total.index(x)+1)/D for x in m2]
    
    # write report
    s = 'Module\tPearson_r\tp-value\n'
    for ii in range(len(ModuleIndex)):
        s += '\t'.join([str(x) for x in [ModuleIndex[ii], m[ii], p_values[ii]]]) + '\n'
    out = open(outfile, 'w')
    out.write(s)
    out.close()
    
    # plot
    if HAS_matplotlib:
        # one can change output_format to 'pdf' if so desired.
        output_format = 'png'
        output_file = outfile.replace('.txt', '') + '.' + output_format
        X = [0.001*x for x in range(-1000, 1000)]
        pyplot.figure()
        pyplot.fill(X, [stats.norm.pdf(x, 0, std(msr)) for x in X], 'grey', alpha=0.2)
        n, bins, patches = pyplot.hist(m, 30, normed=1, facecolor='r', alpha=0.5)
        pyplot.xlabel('Pearson r')
        pyplot.ylabel('PDF')
        pyplot.title(output_file)
        # one can change dpi should higher resolution be desired.
        pyplot.savefig(output_file, dpi=100, format=output_format)
        print("A probability distribution plot is saved to %s" %output_file)
