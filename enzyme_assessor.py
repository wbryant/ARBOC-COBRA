'''
Created on 20 Sep 2015

@author: wbryant
'''

import re
from itertools import combinations, product, chain
from local_utils import dict_append
from random import randrange
from types import StringTypes
from copy import deepcopy

string_file = '/Users/wbryant/work/BTH/data/stringdb/226186.protein.links.v10.txt'
taxonomy_id = '226186'

def load_string_data(string_file=None, taxonomy_id=None):
    """Import String DB data from a single file; it must be of a single organism or risks getting confused.""" 
    string_file=string_file or '/Users/wbryant/work/BTH/data/stringdb/226186.protein.links.v10.txt'
    taxonomy_id=taxonomy_id or '226186'
    
    f_in = open(string_file, 'r')
    f_in.readline()
    relatedness_dict = {}
    gene_scores_dict = {}
    for line in f_in:
        entry = line.strip().split(" ")
        gene1 = re.sub('^.*\.','',entry[0])
        gene2 = re.sub('^.*\.','',entry[1])
        score = float(entry[2])
        relatedness_dict[frozenset([gene1,gene2])] = score
        dict_append(gene_scores_dict, gene1, score)
        dict_append(gene_scores_dict, gene2, score)
    f_in.close()
    
    return relatedness_dict, gene_scores_dict

def load_rel_list(rel_list_file):
    f_in = open(rel_list_file, 'r')
    
    #enzyme_list = []
    rel_list = []
    for line in f_in:
        entry = line.strip().split("\t")
        rel_list.append(tuple(entry))
        #enzyme_list.append(entry[1].split(","))
    f_in.close()
    return rel_list

def get_gene_pair_score(gene_pair, relatedness_dict):
    try:
        score = relatedness_dict[frozenset(gene_pair)]
    except:
        score = 0
    return score

def get_enzyme_scores(enzyme, relatedness_dict):
    """Get scores for all pairs of genes in enzyme; must be comma separated locus_tag list."""
    
    if isinstance(enzyme, StringTypes):
        gene_list = enzyme.split(",")
    else:
        gene_list = enzyme
    if len(gene_list) == 1:
        return None
    
    gene_pairs = list(combinations(gene_list, 2))
    
    gene_pair_scores = []
    for gene_pair in gene_pairs:
        gene_pair_scores.append(get_gene_pair_score(gene_pair, relatedness_dict))
    
    return gene_pair_scores

def analyze_rel_list(rel_list, relatedness_dict):
    """Calculate various scores for all members of an imported REL list."""
    
    enzyme_data_dict = {}
    
    for rel in rel_list:
        enzyme = rel[1]
        gene_set = frozenset(enzyme.split(","))
        if gene_set not in enzyme_data_dict:
            gene_pair_scores = get_enzyme_scores(enzyme, relatedness_dict)
            if gene_pair_scores:
                enzyme_data = []
                enzyme_data.append(len(gene_pair_scores))
                enzyme_data.append(enzyme)
                enzyme_data.append(min(gene_pair_scores))
                enzyme_data.append(sum(gene_pair_scores)/len(gene_pair_scores))
                enzyme_data_dict[gene_set] = enzyme_data
    
    return enzyme_data_dict

def get_benchmark_relatedness(rel_list, relatedness_dict, N=1000, num_genes_max=10):
    """Get relatedness for randomly picked sets of genes within the REL list (assumed metabolic genes)."""
    
    ## Get list of metabolic genes
    metabolic_genes = []
    for entry in rel_list:
        genes = entry[1].split(",")
        metabolic_genes.extend(genes)
    metabolic_genes = list(set(metabolic_genes))
    num_metabolic_genes = len(metabolic_genes)
    
    for num_genes in range(2,num_genes_max+1):
        mean_sum = 0
        min_sum = 0
        for rand_enzyme_idx in range(N):
            gene_indices = []
            while len(gene_indices) < num_genes:
                gene_indices.append(randrange(num_metabolic_genes))
                gene_indices = list(set(gene_indices))
            gene_list = [metabolic_genes[gene_index] for gene_index in gene_indices]
            gene_pair_scores = get_enzyme_scores(gene_list, relatedness_dict)
            mean_sum += sum(gene_pair_scores)/len(gene_pair_scores)
            min_sum += min(gene_pair_scores)
        
        rand_mean = mean_sum/N
        rand_min = min_sum/N
        
        print("For {}-gene enzymes, randomly picked (N={}), mean = {}, min = {}".format(
            num_genes,
            N,
            rand_mean,
            rand_min))

def infer_all_gene_combinations_for_enzyme(cog_list, cog_gene_dict, gene_cog_dict):
    gene_lists = []
    for cog in cog_list:
        gene_lists.append(cog_gene_dict[cog])
    all_gene_combinations = [element for element in product(*gene_lists)]
    all_gene_combinations = [list(set(enzyme)) for enzyme in all_gene_combinations]
    return all_gene_combinations

def powerset_min(iterable, min_length=1):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(min_length,len(s)+1))

def remove_genes_keeping_enzyme(gene_combination, cog_set, gene_cog_dict):
    for gene_subset in powerset_min(gene_combination, min_length=2):
        gene_subset_cog_list = []
        for gene in gene_subset:
            gene_subset_cog_list.extend(gene_cog_dict[gene])
        gene_subset_cog_set = set(gene_subset_cog_list)

def get_enzyme_cogs(gene_list, gene_cog_dict):
    cog_list = []
    for gene in gene_list:
        cog_list.extend(gene_cog_dict[gene])
    return set(cog_list)

def get_gene_by_gene_cogs(gene_list, gene_cog_dict):
    cogs_list = []
    for gene in gene_list:
        cogs_list.append(set(gene_cog_dict[gene]))
    return cogs_list

def add_gene(enzyme_cogs, gene_cog_dict, cog_gene_dict, enzyme_genes = None, fulfilled_cogs = None):
    print "\nIn 'add_gene2()'"
    base_enzyme_genes = enzyme_genes or []
    fulfilled_cogs = fulfilled_cogs or set([])
    
    print "Pre-append enzyme:",enzyme_genes
    next_cog_to_fill = sorted(list(enzyme_cogs-fulfilled_cogs))[0]
    print "Next COG to fill:", next_cog_to_fill
    
    for gene in cog_gene_dict[next_cog_to_fill]:
        
        print "\nIn gene loop, COG:", next_cog_to_fill
        
        ## Check whether any current gene's COGs are a subset of the proposed gene, if so pass
        subset_found = False
        new_gene_cogs = set(gene_cog_dict[gene])
        for gene_cogs_set in get_gene_by_gene_cogs(base_enzyme_genes, gene_cog_dict):
            if gene_cogs_set.issubset(new_gene_cogs):
                subset_found = True
                break
        if subset_found:
            print "Subset found ..."
            continue
        
        print "Appending gene: ", gene
        enzyme_genes = deepcopy(base_enzyme_genes)
        enzyme_genes.append(gene)
        
        fulfilled_cogs = get_enzyme_cogs(enzyme_genes, gene_cog_dict)
        print "Current enzyme: ", enzyme_genes
        print "Fulfilled COGs: ", fulfilled_cogs
        
        if (fulfilled_cogs & enzyme_cogs) == enzyme_cogs:
            print "Yielding: ", enzyme_genes
            enzyme = deepcopy(enzyme_genes)
            enzyme_genes.remove(gene)
            yield enzyme
        else:
            print "Not yielding, looping"
            for enzyme in add_gene(enzyme_cogs, gene_cog_dict, cog_gene_dict, enzyme_genes, fulfilled_cogs):
                yield enzyme


def add_gene_verbose(enzyme_cogs, gene_cog_dict, cog_gene_dict, enzyme_genes = None, fulfilled_cogs = None):
    print "\nIn 'add_gene2()'"
    base_enzyme_genes = enzyme_genes or []
    fulfilled_cogs = fulfilled_cogs or set([])
    
    print "Pre-append enzyme:",enzyme_genes
    next_cog_to_fill = sorted(list(enzyme_cogs-fulfilled_cogs))[0]
    print "Next COG to fill:", next_cog_to_fill
    
    for gene in cog_gene_dict[next_cog_to_fill]:
        
        print "\nIn gene loop, COG:", next_cog_to_fill
        
        ## Check whether any current gene's COGs are a subset of the proposed gene, if so pass
        subset_found = False
        new_gene_cogs = set(gene_cog_dict[gene])
        for gene_cogs_set in get_gene_by_gene_cogs(base_enzyme_genes, gene_cog_dict):
            if gene_cogs_set.issubset(new_gene_cogs):
                subset_found = True
                break
        if subset_found:
            print "Subset found ..."
            continue
        
        print "Appending gene: ", gene
        enzyme_genes = deepcopy(base_enzyme_genes)
        enzyme_genes.append(gene)
        
        fulfilled_cogs = get_enzyme_cogs(enzyme_genes, gene_cog_dict)
        print "Current enzyme: ", enzyme_genes
        print "Fulfilled COGs: ", fulfilled_cogs
        
        if (fulfilled_cogs & enzyme_cogs) == enzyme_cogs:
            print "Yielding: ", enzyme_genes
            enzyme = deepcopy(enzyme_genes)
            enzyme_genes.remove(gene)
            yield enzyme
        else:
            print "Not yielding, looping"
            for enzyme in add_gene_verbose(enzyme_cogs, gene_cog_dict, cog_gene_dict, enzyme_genes, fulfilled_cogs):
                yield enzyme
        
        
        
def benchmark_cogzymes(rel_list, relatedness_dict, cog_membership_file=None):           
    """Enumerate and test every COGzyme-produced enzyme prediction against curated enzymes."""
    
    cog_membership_file=cog_membership_file or "/Users/wbryant/work/BTH/data/COG/BT_cog_associations.csv"
    
    f_in = open(cog_membership_file,'r')
    cog_gene_dict = {}
    gene_cog_dict = {}
    
    for line in f_in:
        entry = line.strip().split("\t")
        cog = entry[0]
        gene_entry = entry[1]
        gene_locus = re.sub('^.*\.','',gene_entry)
        dict_append(cog_gene_dict,cog,gene_locus)
        dict_append(gene_cog_dict,gene_locus,cog)
    f_in.close()
    
    for rel in rel_list:
        gene_list = rel[1].split(",")
        
        ## Establish COGzyme
        cog_list = []
        for gene in gene_list:
            cog_list.extend(gene_cog_dict[gene])
        cog_list = list(set(cog_list))
        
        cog_members_count = {}
        for cog in cog_list:
            cog_members_count[cog] = 0
        
        
        for gene in cog_gene_dict[cog_list[0]]:
            working_cog_set = deepcopy(set(cog_list))
            gene_cogs = gene_cog_dict[gene]
            for cog in gene_cogs:
                working_cog_set.remove(cog)
            
            
    
    
        
        
        
        
    
     
    
    