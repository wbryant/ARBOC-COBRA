'''
Created on 20 Sep 2015

@author: wbryant
'''

import re
from itertools import combinations, product, chain
from local_utils import dict_append, preview_dict, loop_counter
from random import randrange
from types import StringTypes
from copy import deepcopy

from abcsmc import create_extended_model
from local_gene_parser import gene_parser
from annotation.models import Reaction

def mean(obj_list):
    return sum(obj_list)/float(len(obj_list))

def infer_rel_list(model_file):
    """Infer the list of RELs from an SBML file, preferring MNX IDs where possible."""
    
    model = create_extended_model(model_file, require_solver=False)
    
    rel_list = []
    counter = loop_counter(len(model.reactions), "Running through model reactions")
    for reaction in model.reactions:
        model_rxn_id = re.sub('(\_enz\d*)*$','',reaction.id)
        db_reactions = Reaction.objects.filter(reaction_synonym__synonym=model_rxn_id).distinct()
        if len(db_reactions) == 1:
            reaction_id = db_reactions[0].name
        else:
            reaction_id = model_rxn_id   
        for enzyme in gene_parser(reaction.gene_reaction_rule):
            if ((len(enzyme) > 0) & (enzyme != '0000000.0.peg.0')):
                rel_list.append((reaction_id, enzyme))
        counter.step()
    counter.stop()
    return rel_list

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

def get_enzyme_cogs(proposed_enzyme, target_enzyme_cogs, gene_cog_dict):
    """Get all COGs of which there is at least one member gene in the enzyme, limited to COGs in the target enzyme."""
    cog_list = []
    for gene in proposed_enzyme:
        cog_list.extend(set(gene_cog_dict[gene]) & target_enzyme_cogs)
    return set(cog_list)

def get_gene_by_gene_cogs(proposed_enzyme, target_enzyme_cogs, gene_cog_dict):
    """Get COG memberships for each gene in enzyme, limited to COGs in the target enzyme."""
    cogs_list = []
    for gene in proposed_enzyme:
        cogs_list.append(set(gene_cog_dict[gene]) & target_enzyme_cogs)
    return cogs_list

def add_gene(target_enzyme_cogs, gene_cog_dict, cog_gene_dict, proposed_enzyme = None, fulfilled_cogs = None):
    """Generator for all enzymes with membership of all target_enzyme_cogs, without duplicate enzymes or redundant genes."""
    
    base_enzyme_genes = proposed_enzyme or []
    fulfilled_cogs = get_enzyme_cogs(base_enzyme_genes, target_enzyme_cogs, gene_cog_dict)
    
    ## Which COG will we try to find a member of?
    next_cog_to_fill = sorted(list(target_enzyme_cogs-fulfilled_cogs))[0]
    gene_members_of_cog = cog_gene_dict[next_cog_to_fill] 
    
    for gene in gene_members_of_cog:
        
        ## Check whether any already-present gene's COG set is a subset of the proposed gene's COG set, if so skip addition
        subset_found = False
        proposed_gene_cogs = set(gene_cog_dict[gene]) & target_enzyme_cogs
        for gene_cogs_set in get_gene_by_gene_cogs(base_enzyme_genes, target_enzyme_cogs, gene_cog_dict):
            if gene_cogs_set.issubset(proposed_gene_cogs):
                subset_found = True
                break
        if subset_found:
            continue
        
        ## Add gene to proposed enzyme
        proposed_enzyme = deepcopy(base_enzyme_genes)
        proposed_enzyme.append(gene)
        
        ## Determine which COG memberships are fulfilled by the genes in the proposed enzyme
        fulfilled_cogs = get_enzyme_cogs(proposed_enzyme, target_enzyme_cogs, gene_cog_dict)
        
        if (fulfilled_cogs & target_enzyme_cogs) == target_enzyme_cogs:
            ## Proposed enzyme has members of every required COG, so yield 
            enzyme = deepcopy(proposed_enzyme)
            proposed_enzyme.remove(gene)
            yield enzyme
        else:
            ## Proposed enzyme is still missing some COG members
            for enzyme in add_gene(target_enzyme_cogs, gene_cog_dict, cog_gene_dict, proposed_enzyme, fulfilled_cogs):
                yield enzyme
            
def load_cog_data(taxonomy_id, rel_list, cog_member_file=None):
    """Get all relevant gene-COG relationships from member file (eggNOG 4.1)."""
    
    cog_member_file=cog_member_file or "/Users/wbryant/work/data/NOG/NOG.members.tsv"
    cog_gene_dict = {}
    gene_cog_dict = {}
    
    f_in = open(cog_member_file, 'r')
    for line in f_in:
        cog_data = line.strip().split("\t")
        cog = cog_data[1]
        #print "COG:", cog
        cog_members = cog_data[5]
        organism_query = "(" + taxonomy_id + "\.)" + "([^,\s]+)"
        #print "Query:", organism_query
        for entry in re.findall(organism_query, cog_members):
            
            gene = entry[1]
            #print "Gene:", gene
            dict_append(gene_cog_dict, gene, cog)
            dict_append(cog_gene_dict, cog, gene)
    f_in.close()
    
    ## For those genes without COGs ensure they are in the dictionary
    all_genes_in_rels = []
    for entry in rel_list:
        gene_list = entry[1].split(",")
        all_genes_in_rels.extend(gene_list)
    all_genes_in_rels = list(set(all_genes_in_rels))
    for gene in all_genes_in_rels:
        if gene not in gene_cog_dict:
            gene_cog_dict[gene] = []
    
    return gene_cog_dict, cog_gene_dict


def benchmark_cogzymes(rel_list, relatedness_dict, taxonomy_id, cog_membership_file=None):           
    """Enumerate and test every COGzyme-produced enzyme prediction against curated enzymes."""
    
    cog_membership_file=cog_membership_file or "/Users/wbryant/work/data/NOG/NOG.members.tsv"
    
    gene_cog_dict, cog_gene_dict = load_cog_data(taxonomy_id, rel_list, cog_membership_file)
    
    ## Establish which enzymes from each COGzyme are in the model
    cogzyme_enzyme_dict = {}
    valid_rel_list = []
    for rel in rel_list:
        gene_list = rel[1].split(",")
        if len(gene_list) == 1:
            continue
        cog_list = []
        non_cog_gene=False
        for gene in gene_list:
            gene_cogs = gene_cog_dict[gene]
            if len(gene_cogs) > 0:    
                cog_list.extend(gene_cogs)
            else:
                non_cog_gene=True
                break
        if non_cog_gene:
#             print "REL had gene without COGs:", gene, ";", rel
            continue
        valid_rel_list.append(rel)
        enzyme_cogs = frozenset(cog_list)
        gene_set = set(gene_list)
        dict_append(cogzyme_enzyme_dict, enzyme_cogs, gene_set)
        
    
    all_model_mean_scores = {}
    all_pred_mean_scores = {}
    all_model_min_scores = {}
    all_pred_min_scores = {}
    mean_ratios = []
    min_ratios = []
    
    counter = loop_counter(len(valid_rel_list), "Inferring predicted enzymes from REL list")
    for rel in valid_rel_list:
#         print "\nREL:", rel
        
        gene_list = rel[1].split(",")
        
        ## Establish COGzyme
        cog_list = []
        for gene in gene_list:
            cog_list.extend(gene_cog_dict[gene])
        target_enzyme_cogs = set(cog_list)
        
        ## Get list of enzymes for this COGzyme already in the model
        enzymes_in_model = cogzyme_enzyme_dict[frozenset(target_enzyme_cogs)]
        
        model_mean_scores = []
        pred_mean_scores = []
        model_min_scores = []
        pred_min_scores = []
        
        enzyme_list = []
        for enzyme in add_gene(target_enzyme_cogs, gene_cog_dict, cog_gene_dict):
            enzyme_list.append(set(enzyme))
            
        rel_enzyme = set(gene_list)
        
        if rel_enzyme not in enzyme_list:
            enzyme_list.append(rel_enzyme)
        
        for enzyme in enzyme_list:
            if len(enzyme) > 1:
                gene_pair_scores = get_enzyme_scores(enzyme, relatedness_dict)
                if set(enzyme) in enzymes_in_model:
                    dict_append(all_model_mean_scores, len(enzyme), sum(gene_pair_scores)/float(len(gene_pair_scores)))
                    dict_append(all_model_min_scores, len(enzyme), min(gene_pair_scores))
                    model_mean_scores.append(sum(gene_pair_scores)/float(len(gene_pair_scores)))
                    model_min_scores.append(min(gene_pair_scores))
                else:
                    dict_append(all_pred_mean_scores, len(enzyme), sum(gene_pair_scores)/float(len(gene_pair_scores)))
                    dict_append(all_pred_min_scores, len(enzyme), min(gene_pair_scores))
                    pred_mean_scores.append(sum(gene_pair_scores)/float(len(gene_pair_scores)))
                    pred_min_scores.append(min(gene_pair_scores))
        counter.step()
        
    counter.stop()
    output_answers = {
        'all_model_min_scores': all_model_min_scores,
        'all_pred_min_scores': all_pred_min_scores,
        'all_model_mean_scores': all_model_mean_scores,
        'all_pred_mean_scores': all_pred_mean_scores,
    }
     
    return output_answers


def get_enzyme_scores_for_model(model_file=None, taxonomy_id=None, string_file=None, cog_file=None):
    model_file=model_file or "/Users/wbryant/git/incognito/static/models/BTH_iAH991_w_gprs.xml"
    cog_file=cog_file or "/Users/wbryant/work/data/NOG/NOG.members.tsv"
    string_file=string_file or '/Users/wbryant/work/BTH/data/stringdb/226186.protein.links.v10.txt'
    taxonomy_id=taxonomy_id or '226186'
        
    print "Inferring REL list ..."
    model_rel_list = infer_rel_list(model_file)
    
    print "Inferring relatedness dictionary ..."
    relatedness_dict, gene_scores_dict = load_string_data(string_file)
    
    print "Inferring COGzyme enzymes and scoring ..."
    results = benchmark_cogzymes(model_rel_list, relatedness_dict, taxonomy_id)
    
#     for dict in results:
#         print ""
#         for num_genes, scores in dict.iteritems():
#             print num_genes,len(scores), sum(scores)/len(scores)
        
    return results
    
def report_on_results(results):
    
    if 'mean_ratios' in results:
        mean_ratios = results['mean_ratios']
        print "Average ratio of model/prediction mean enzyme scores is", mean(mean_ratios)
    
    if 'min_ratios' in results:
        min_ratios = results['min_ratios']
        print "Average ratio of model/prediction minimum gene-gene enzyme scores is", mean(min_ratios)
    
    return None

def collate_results_sets(results_set_list):
    return None

    
    
    