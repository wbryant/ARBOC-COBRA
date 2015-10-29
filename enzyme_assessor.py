'''
Created on 20 Sep 2015

@author: wbryant
'''

import re, sys
from itertools import combinations, product, chain
from local_utils import dict_append, preview_dict, loop_counter
from random import randrange
from types import StringTypes
from copy import deepcopy

from abcsmc import create_extended_model
from local_gene_parser import gene_parser
from annotation.models import Reaction

from statsmodels.api import Logit
from scipy.misc import comb
from pandas import qcut

def mean(obj_list):
    return sum(obj_list)/float(len(obj_list))

def infer_rel_list(model_file):
    """Infer the list of RELs from an SBML file, preferring MNX IDs where possible."""
    
    model = create_extended_model(model_file, require_solver=False)
    
    rel_list = []
    counter = loop_counter(len(model.reactions), "Running through model reactions")
    gene_list = []
    model_rxn_ids = []
    for reaction in model.reactions:
        model_rxn_id = re.sub('(\_enz\d*)*$','',reaction.id)
        model_rxn_ids.append(model_rxn_id)
    db_reaction_data = Reaction.objects.filter(reaction_synonym__synonym__in=model_rxn_ids).values_list('reaction_synonym__synonym','name')
    db_reaction_dict = {}
    for item in db_reaction_data:
        dict_append(db_reaction_dict, item[0], item[1], ignore_duplicates=True)    
    
    for reaction in model.reactions:
        model_rxn_id = re.sub('(\_enz\d*)*$','',reaction.id)
        try:
            db_reactions_specific = db_reaction_dict[model_rxn_id]
        except:
            db_reactions_specific = ''
        if len(db_reactions_specific) == 1:
            reaction_id = db_reactions_specific[0]
        else:
            reaction_id = model_rxn_id   
        for enzyme in gene_parser(reaction.gene_reaction_rule):
            if ((len(enzyme) > 0) & (enzyme != '0000000.0.peg.0')):
                rel_list.append((reaction_id, enzyme, reaction.source))
                genes = enzyme.split(",")
                for gene in genes:
                    gene_list.append(gene)
        counter.step()
    counter.stop()
    gene_list = list(set(gene_list))
    return rel_list, gene_list

def infer_rel_list_slow(model_file):
    """Infer the list of RELs from an SBML file, preferring MNX IDs where possible."""
    
    model = create_extended_model(model_file, require_solver=False)
    
    rel_list = []
    counter = loop_counter(len(model.reactions), "Running through model reactions")
    gene_list = []
    for reaction in model.reactions:
        model_rxn_id = re.sub('(\_enz\d*)*$','',reaction.id)
        db_reactions = Reaction.objects.filter(reaction_synonym__synonym=model_rxn_id).distinct()
        if len(db_reactions) == 1:
            reaction_id = db_reactions[0].name
        else:
            reaction_id = model_rxn_id   
        for enzyme in gene_parser(reaction.gene_reaction_rule):
            if ((len(enzyme) > 0) & (enzyme != '0000000.0.peg.0')):
                rel_list.append((reaction_id, enzyme, reaction.source))
                genes = enzyme.split(",")
                for gene in genes:
                    gene_list.append(gene)
        counter.step()
    counter.stop()
    gene_list = list(set(gene_list))
    return rel_list, gene_list

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
    
    cog_member_file=cog_member_file or "/Users/wbryant/work/data/NOG/NOG.members_adapted.tsv"
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
            #print "Gene:", gene, "COG:", cog
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
            print("Gene not found in NOG file: '{}'".format(gene))
            gene_cog_dict[gene] = []
            print(" -> dict entry: '{}'".format(gene_cog_dict[gene]))
    
    return gene_cog_dict, cog_gene_dict


def benchmark_cogzymes(rel_list, relatedness_dict, taxonomy_id, operon_dict, gene_cog_dict_in, cog_gene_dict_in, max_preds_per_rel=None):           
    """Enumerate and test every COGzyme-produced enzyme prediction against curated enzymes."""

#     print("Gene-cog dict: {}".format(len(gene_cog_dict_in)))
#     preview_dict(gene_cog_dict_in, 20)
#     
#     print("Cog-gene dict: {}".format(len(cog_gene_dict_in)))
#     preview_dict(cog_gene_dict_in, 20)

    ## Check all genes for presence in Operon data (currently used as the standard)
    gene_cog_dict = {} 
    cog_gene_dict = {}
    for gene, coglist in gene_cog_dict_in.iteritems():
        if gene in operon_dict:
            gene_cog_dict[gene] = coglist
    for cog, genelist in cog_gene_dict_in.iteritems():
        for gene in genelist:
            if gene not in operon_dict:
                genelist.remove(gene)
        if len(genelist) > 0:
            cog_gene_dict[cog] = genelist
    
#     print("Gene-cog dict: {}".format(len(gene_cog_dict)))
#     preview_dict(gene_cog_dict, 20)
#     
#     print("Cog-gene dict: {}".format(len(cog_gene_dict)))
#     preview_dict(cog_gene_dict, 20)
    
    ## Establish which enzymes from each COGzyme are in the model
    cogzyme_enzyme_dict = {}
    valid_rel_list = []
    max_preds_per_rel=max_preds_per_rel or 100
    
    
    print("There are {} RELs".format(len(rel_list)))
    num_rels_failed = 0
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
                num_rels_failed += 1
                #print "REL:", rel
                #print(" - gene '{}' not found".format(gene))
                non_cog_gene=True
                break
        if non_cog_gene:
#             print "REL had gene without COGs:", gene, ";", rel
            continue
        valid_rel_list.append(rel)
        enzyme_cogs = frozenset(cog_list)
        gene_set = set(gene_list)
        dict_append(cogzyme_enzyme_dict, enzyme_cogs, gene_set)
    
    print("{} RELs failed".format(num_rels_failed))
    
    all_model_mean_scores = {}
    all_pred_mean_scores = {}
    all_model_min_scores = {}
    all_pred_min_scores = {}
    mean_ratios = []
    min_ratios = []
    
    all_min_scores = {}
    all_mean_scores = {}
    dataset = []
    all_enzymes_found = []
    
    counter = loop_counter(len(valid_rel_list), "Inferring predicted enzymes from REL list ({})".format(len(valid_rel_list)))
    for rel in valid_rel_list:
#         print "\nREL:", rel
        
        rel_source=rel[2]
        
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
        
        ## Cutoff for huge COGzymes - nothing useful from them as far too many predictions
        enzyme_list = []
        num_preds = 0
        too_many_preds=False
        for enzyme in add_gene(target_enzyme_cogs, gene_cog_dict, cog_gene_dict):
            enzyme_list.append(set(enzyme))
            num_preds += 1
            if num_preds > max_preds_per_rel:
                too_many_preds=True
                break
        if too_many_preds:
            continue
        
        ## Check that original enzyme is in list
        rel_enzyme = set(gene_list)
        if rel_enzyme not in enzyme_list:
            enzyme_list.append(rel_enzyme)
        
        ## Do scoring
        for enzyme in enzyme_list:
            enzyme_length = len(enzyme)
            enzyme_orig = deepcopy(enzyme)
            if enzyme_length > 1:
                    
                gene_pair_scores = get_enzyme_scores(enzyme, relatedness_dict)
                try:
                    mean_score = sum(gene_pair_scores)/float(len(gene_pair_scores))
                except:
                    print enzyme
                    print enzyme_length
                    print gene_pair_scores
                    mean_score = sum(gene_pair_scores)/float(len(gene_pair_scores))
                    
                min_score = min(gene_pair_scores)
                
                if set(enzyme) in enzymes_in_model:
                    present_in_model = 1
                    dict_append(all_model_mean_scores, enzyme_length, mean_score)
                    dict_append(all_model_min_scores, enzyme_length, min_score)
                    model_mean_scores.append(mean_score)
                    model_min_scores.append(min_score)
                    dict_append(all_min_scores, enzyme_length, (min_score, 1))
                    dict_append(all_mean_scores, enzyme_length, (mean_score, 1))
                else:
                    present_in_model = 0
                    dict_append(all_pred_mean_scores, enzyme_length, mean_score)
                    dict_append(all_pred_min_scores, enzyme_length, min_score)
                    pred_mean_scores.append(mean_score)
                    pred_min_scores.append(min_score)
                    dict_append(all_min_scores, enzyme_length, (min_score, 0))
                    dict_append(all_mean_scores, enzyme_length, (mean_score, 0))
                
                if enzyme not in all_enzymes_found: 
                    all_enzymes_found.append(enzyme)
                    operon_memberships = []
                    for gene in list(enzyme):
                        try:
                            operon_memberships.append(operon_dict[gene])
                        except:
                            print rel, enzyme, enzyme_orig, gene, "\n"
                            operon_memberships.append(operon_dict[gene])
                    operon_score = 0 
                    for member1 in range(enzyme_length):
                        for member2 in range(member1+1, enzyme_length):
                            if operon_memberships[member1] == operon_memberships[member2]:
                                operon_score += 1
                    dataset.append([mean_score, min_score, enzyme_length, operon_score/comb(enzyme_length,2), present_in_model,rel_source])
                
        counter.step()
        
    counter.stop()
    output_answers = {
        'all_model_min_scores': all_model_min_scores,
        'all_pred_min_scores': all_pred_min_scores,
        'all_model_mean_scores': all_model_mean_scores,
        'all_pred_mean_scores': all_pred_mean_scores,
        'all_min_scores': all_min_scores,
        'all_mean_scores': all_mean_scores,
        'dataset': dataset
    }
     
    return output_answers

def load_operon_data(operon_data_file):
    
    operon_membership_dict = {}
    
    f_in = open(operon_data_file, 'r')
    for line in f_in:
        try:
            cols = line.strip().split("\t")
            operon_id = int(cols[0])
            gene_locus_tag = cols[2]
            operon_membership_dict[gene_locus_tag] = operon_id
        except:
            pass
    f_in.close()
    return operon_membership_dict
    

def get_enzyme_scores_for_model(model_file=None, taxonomy_id=None, string_file=None, cog_file=None, operon_file=None, model_rel_list=None, relatedness_dict=None):
    model_file=model_file or "/Users/wbryant/git/incognito/static/models/BTH_iAH991_w_gprs.xml"
    cog_file=cog_file or "/Users/wbryant/work/data/NOG/NOG.members_adapted.tsv"
    string_file=string_file or '/Users/wbryant/work/BTH/data/stringdb/226186.protein.links.v10.txt'
    operon_file=operon_file or "/Users/wbryant/work/BTH/data/proopdb/bth_operon_predictions.txt"
    taxonomy_id=taxonomy_id or '226186'
    model_rel_list=model_rel_list or None
    relatedness_dict=relatedness_dict or None
    
    if not model_rel_list:
        print "Inferring REL list ..."
        model_rel_list, _ = infer_rel_list(model_file)
    else:
        print("REL list given ...")
    
    if not relatedness_dict:
        print "Inferring relatedness dictionary ..."
        relatedness_dict, gene_scores_dict = load_string_data(string_file)
    else:
        print("Relatedness dictionary given ...")

    print("Loading COG data ...")
    gene_cog_dict, cog_gene_dict = load_cog_data(taxonomy_id, model_rel_list, cog_file)
    
    print("Importing operon data ...")
    operon_dict = load_operon_data(operon_file)
    
    print "Inferring COGzyme enzymes and scoring ..."
    results = benchmark_cogzymes(model_rel_list, relatedness_dict, taxonomy_id, operon_dict, gene_cog_dict, cog_gene_dict)
    
#     for dict in results:
#         print ""
#         for num_genes, scores in dict.iteritems():
#             print num_genes,len(scores), sum(scores)/len(scores)
        
    return results, model_rel_list, relatedness_dict
    
def report_on_results(results):
    
    if 'mean_ratios' in results:
        mean_ratios = results['mean_ratios']
        print "Average ratio of model/prediction mean enzyme scores is", mean(mean_ratios)
    
    if 'min_ratios' in results:
        min_ratios = results['min_ratios']
        print "Average ratio of model/prediction minimum gene-gene enzyme scores is", mean(min_ratios)
    
    return None

def solve_logit_fixed_var(dataset, fixed_var, value, target_var=None):
    """Take 'dataset' and for a particular 'value' of 'fixed_col' do logistic regression."""
    
    target_var=target_var or 'Enzyme'
    training_cols= [col for col in dataset.columns[:]]
    try:
        training_cols.remove(target_var)
        training_cols.remove(fixed_var)
    except:
        print("Column names could not be found")
        sys.exit(1)
        
        
    selected_dataset=dataset[dataset[fixed_var] == value]
    logit = Logit(selected_dataset[target_var],selected_dataset[training_cols])
    result = logit.fit(method="bfgs")
    return result, selected_dataset
    
def get_qcut_bin_centres(column, num_bins):
    """Use a quantile binning, and return the mid-points of each bin as a list."""
    
    _, edges = qcut(column, num_bins, retbins=True)
    bin_centres = []
    for idx, _ in enumerate(edges[1:]):
        centre = (edges[idx+1] + edges[idx])/float(2)
        bin_centres.append(centre)
    
    return bin_centres