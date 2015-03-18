'''
Created on 10 Mar 2014

@author: wbryant
'''
import re, sys
from cobra.core.Model import Model
from cobra.core import ArrayBasedModel
from cobra.io import write_sbml_model
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.delete import delete_model_genes, undelete_model_genes, find_gene_knockout_reactions
from cobra import Reaction, Metabolite
from local_gene_parser import gene_parser
from local_utils import *
import numpy as np
from numpy import log as ln
from copy import deepcopy
from random import random
from math import sqrt, log10, ceil, exp
import shelve
from collections import Counter

def import_prior_dict(
        prior_files,
        rem_default = 0.9,
        add_default = 0.1):
    """Create prior dictionary from prior files, giving default prior values 
    where not specified."""
    prior_dict = {}
    if not isinstance(prior_files, list):
        prior_files = [prior_files]
    for filename in prior_files:
        f_in = open(filename, 'r')
        for line in f_in:
            entry = line.strip().split("\t")
            if entry[1] == "rem":
                prior_dict[entry[0]] = rem_default
            elif entry[1] == "add":
                prior_dict[entry[0]] = add_default
            else:
                try:
                    prior_val = float(entry[1])
                    prior_dict[entry[0]] = prior_val
                except:
                    print -1
        f_in.close()
    return prior_dict      

def import_media(media_shelf_file=None):
    """Import media values from a media shelf file for input into experiments."""
    media_shelf_file = media_shelf_file or '/Users/wbryant/work/BTH/analysis/fba/media.data'
    media_shelf = shelve.open(media_shelf_file)
    media = {}
    media['min'] = media_shelf['min']
    media['tyg'] = media_shelf['tyg']
    media_shelf.close()
    return media


    
def import_expt_data(model, media = None, objective_id="Biomass_BT_v2", data_file = None):
    """From a table of genotype/growth information, create the experiments 
    for a particular organism.
    
    If there are any lines at the top of the file starting with the word 
    "MEDIUM" in the first column the next column will be the medium name and 
    the following cells will contain medium components for specifying media 
    used in the experiments. 
    """
    
    objective = model.reactions.get_by_id(objective_id)
    
    data_file = data_file or '/Users/wbryant/work/BTH/analysis/gene_essentiality/essentiality_data_complete.csv'
    experiments = []
    expts_in = open(data_file,'r')
    if not media:
        media = {}
    line_number = 0
    for line in expts_in:
#         line_number += 1
        if line[0] != "#":
            details = line.strip().split("\t")
            if details[0] == 'MEDIUM':
                medium_name = details[1]
                medium_components = details[2:]
                medium_dict = {}
                for component in medium_components:
                    medium_dict[component.strip()] = -1
                media[medium_name] = medium_dict
                try:
                    model.set_medium(medium_dict)
                except:
                    print("Medium '{}' cannot be set on the model - incorrect component ID?".format(medium_name))
                    for component in medium_components:
                        print (" - {}".format(component))        
            else:
                
                medium_base = details[1]
                c_source = details[2]
                genotype = details[3].split(" and ")
                result = details[4]
                
#                 print line_number
#                 print details
#                 print medium_base
                
                c_sources = []
                
                ## Add carbon source to medium
                if c_source == "-":
                    c_sources = []
                elif 'and' in c_source:
                    c_sources = c_source.split(' and ')
                elif c_source:
                    c_sources = [c_source]
                else:
                    c_sources = None
                expt_medium = deepcopy(media[medium_base])
                if c_sources:
                    for source in c_sources:
                        c_source_exchange = "EX_" + source + "(e)"
                        expt_medium[c_source_exchange] = -100
                
                
                experiment = Experiment(details[0], expt_medium, result, genotype, objective)
                
                experiments.append(experiment)
    expts_in.close()
    return experiments
 

def import_rxn_ids_from_file(file_in, criterion = None):
    """Import rxn IDs from the first column in file_in.
    Filter by criterion value in the second column, if specified.
    """
    
    f_in = open(file_in, 'r')
    
    rxn_ids = []
    
    line_num = 0
    
    for line in f_in:
        line_num += 1
        rxn_id = line.split("\t")[0].strip()
        if criterion is not None:
            try:
                if line.split("\t")[1].strip() == criterion:
                    rxn_ids.append(rxn_id)
            except:
                print line_num, ": '", line, "'"      
        else:
            rxn_ids.append(rxn_id)
    f_in.close()
    
    return rxn_ids
        
# def get_p_transition(t,T,p_0 = 0.1, type='linear'):
#     if type == 'linear':
#         step_size = p_0 / T
#         p = p_0 - step_size * t
#         return p
#     else:
#         return None
#    
# 
# def get_epsilon(t,T, epsilon_0 = 0.5,type='linear'):
#     if type == 'linear':
#         step_size = epsilon_0 / T
#         epsilon = epsilon_0 - step_size * t
#         return epsilon
#     else:
#         return None
      

class AbcProblem():
    
    def __init__(self,
            particles_per_population,
            num_populations_max,
            model,
            abc_reactions = None,
            prior_dict = None,
            rxn_enz_priors = None,
            default_prior_value = 0.95,
            enzyme_limit = 10,
            objective_name = 'Biomass',
            epsilon_0 = 0.3,
            epsilon_T = 0.285,
            p_0 = 0.2,
            experiments=None,
            alpha = None,
            include_all=False,
            prior_file=None):

        """Set up ABC for given model
        
        num_populations = number of iterations to run through;
        model = ExtendedCobraModel to be tested;
        get_epsilon = a function of t to calculate epsilon;
        get_p_transition = a function of t to calculate p_transition;
        abc_reactions: list of reaction IDs to be included;
        prior_dict: dictionary of predefined prior values for specific reactions (rxn ID is key);
        rxn_enz_priors: a list with entries so - ['base_rxn_id',[gene_list],prior];
        default_prior: default prior;
        include_all: if True, include all non-exchange reactions in the ABC
        
        """
        self.t = 0
        self.N = particles_per_population
        self.T = num_populations_max
        self.model = deepcopy(model)
        if experiments:
            self.experiments = experiments
        else:
            self.experiments = import_experiments(self.model, objective_name)
        
        ## Give total numbers of positive and negative results for calculating test cutoffs
        self.num_expts_pos = 0
        self.num_expts_neg = 0
        for expt in self.experiments:
            if expt.result == 1:
                self.num_expts_pos += 1
            else:
                self.num_expts_neg += 1
        
        
        self.num_rxns = len(self.model.reactions)
        
#         ## Only a few enzymes will be completely incorrect, and be replaceable 
#         ## by an alternative.  Set prior for non-enzyme enzrxns to 1/num_rxns 
#         ## in model.
#         self.non_enz_prior = 1/float(self.num_rxns)
        
        
        self.p_0 = p_0
        self.epsilon_0 = epsilon_0
        self.epsilon_T = epsilon_T
        self.alpha = alpha or 0.01
        self.get_p_transition()
        self.get_epsilon()
        self.theta_set_prev = None
        self.w_set_prev_unnorm = None
        self.w_set_prev = None
        self.theta_set = None
        self.w_set = None
        self.include_all=include_all
        self.default_prior_value=default_prior_value
        
        if self.include_all:
            abc_reactions = []
            for reaction in self.model.reactions:
                if reaction.id[0:2] != 'EX_':
                    abc_reactions.append(reaction.id) 
        
        ## Set particle ID format according to how many particles there are
        id_length = str(int(ceil(log10(self.N))))
        self.id_format = "{:0" + id_length + "d}"
        
        ## Initialise results variables, for storing run results
        self.results_theta = []
        self.results_w = []
        self.results_d = []
        
        ## This is for storing the theta sets at each timepoint
        self.intermediate_theta_dict = {}
        
        if not prior_dict:
            prior_dict = {}
        
        if prior_file:
            f_in = open(prior_file, 'r')
            for line in f_in:
                details = line.strip().split("\t")
                prior_dict[details[0]] = float(details[1])
            f_in.close()
        
        if abc_reactions:
            ## Split all included reactions into individual enzyme/reaction pairs
            counter = loop_counter(len(abc_reactions),'Splitting ABC reactions')
            
            for rxn_id in abc_reactions:        
                counter.step()
                enzrxn_ids, non_enz_rxn_id = self.model.split_rxn_by_enzymes(rxn_id, enzyme_limit)
                num_enzrxns = len(enzrxn_ids)
                if rxn_id in prior_dict:
                    prior_value = prior_dict[rxn_id]/sqrt(float(num_enzrxns))
                else:
                    prior_value = default_prior_value/sqrt(float(num_enzrxns))
                for enzrxn_id in enzrxn_ids:
                    prior_dict[enzrxn_id] = prior_value
                if non_enz_rxn_id:
                    prior_dict[non_enz_rxn_id] = self.default_prior_value
                
            counter.stop()
        
        ## Apply belief about rxn/gene/enzyme relationships
        prior_dict = self.set_rxn_enz_priors(rxn_enz_priors, prior_dict)
        
        ## All beliefs about reactions included in the model are now in prior_dict.
#         ## Any reaction not in prior_dict is not in the ABC and should always be included.
        self.prior_set = np.ones(len(prior_dict))
        self.model.theta_rxn_map = {}
        idx = -1
        for rxn_id, prior_value in prior_dict.iteritems():
            idx += 1
            self.prior_set[idx] = prior_value
            self.model.theta_rxn_map[idx] = self.model.reactions.get_by_id(rxn_id)
        
        ## Create original lb/ub vars for reference
        for rxn in self.model.reactions:
            rxn.lb_orig = deepcopy(rxn.lower_bound)
            rxn.ub_orig = deepcopy(rxn.upper_bound)
        
        ## Calculate ln_pi_max for calculation of pi_rel in weight calculations
        ln_pi_max = 0
        for p in self.prior_set:
            if p < 0.5:
                ln_pi_max += ln(1-p)
            else:
                ln_pi_max += ln(p) 
        self.ln_pi_max = ln_pi_max
        
        ## Finding media with >10 experiments, for ensuring testing does not waste time
        media_list = []
        for expt in self.experiments:
            media_list.append(expt.medium_components)
        media_count = Counter(media_list)
        precalc_media_lists = []
        for medium_count in media_count.most_common():
            if medium_count[1] <= 10:
                break
            else:
                precalc_media_lists.append(medium_count[0])        
        precalc_media = []
        for medium_list in precalc_media_lists:
            precalc_medium = {}
            for component in medium_list:
                precalc_medium[component] = -10
            precalc_media.append(precalc_medium)
        self.precalc_media = precalc_media
        self.precalc_media_frozensets = precalc_media_lists
        
        print("ABC initialisation complete.")
                
    
    def set_rxn_enz_priors(self, rxn_enz_priors, prior_dict, include_all=False):
        """Where there is belief about certain enzyme/reaction pairs, 
        apply this to the relevant reaction priors.
        """
        
        if rxn_enz_priors:
            for rxn_enz_prior in rxn_enz_priors:
                base_rxn_id = rxn_enz_prior[0]
                gene_set = set(rxn_enz_prior[1])
                prior_value = rxn_enz_prior[2]
                id_string = base_rxn_id + "(_enz[0-9]+)*$"
                for rxn in self.model.reactions:
                    if re.match(id_string, rxn.id):
                        enzrxn_gene_set = set(rxn.genes)
                        if (enzrxn_gene_set | gene_set) == (enzrxn_gene_set & gene_set):
                            prior_dict[rxn.id] = prior_value
        return prior_dict


    def initialise_particle(self, particle_id):
        """Create a new particle from the current ABC timepoint."""
        particle = Particle(
            self.theta_set_prev,
            self.w_set_prev,
            self.model,
            self.experiments,
            self.epsilon_t,
            self.prior_set,
            self.p_t,
            self.t,
            self.N,
            particle_id,
            self.num_expts_pos,
            self.num_expts_neg,
            self.ln_pi_max,
            self.precalc_media,
            self.precalc_media_frozensets)
        
        return particle
    
    def create_population_t(self):
        self.population_t = []
        for particle_id in range(self.N):
            particle_id_string = self.id_format.format(particle_id)
            self.population_t.append(self.initialise_particle(particle_id_string))
    
        
    def record_previous_results(self, theta_accepted_set, ln_w_accepted_set, distance_set):
        """Add the results from the previous run to the full lists of results."""

        self.results_theta.append(theta_accepted_set)
        self.results_d.append(distance_set)
        self.results_w.append(self.w_set_prev)
        self.intermediate_theta_dict[self.t] = theta_accepted_set
        self.theta_set_prev = deepcopy(theta_accepted_set)        
        
        
    def step_forwards(self, ln_w_accepted_set):
        """Increment time and calculate new problem parameters."""
        
        self.update_weights(ln_w_accepted_set)
        if self.t < self.T-1:            
            self.t += 1
        else:
            print("This was the final run.")
            return None
        self.get_p_transition()
        self.get_epsilon()               
           
        
    def update_weights(self, ln_weights):
        """Normalise outputed weights and update w_set_prev."""
        max_ln_w = max(ln_weights)
        ln_ws_full = []
        for ln_w in ln_weights:
            if not ln_w:
                ln_ws_full.append(0)
            else:
                ln_ws_full.append(ln_w-max_ln_w)
        weights_unnorm = [exp(ln_w) for ln_w in ln_ws_full]
        sum_weights = float(sum(weights_unnorm))
        weights_norm = [w/sum_weights for w in weights_unnorm]
        self.w_set_prev_unnorm = weights_unnorm 
        self.w_set_prev = weights_norm
    
    def get_epsilon_linear(self):
        step_size = (self.epsilon_0 - self.epsilon_T) / (float(self.T) - 1)
        epsilon_t = self.epsilon_0 - step_size * self.t
        self.epsilon_t = epsilon_t
    
    def get_epsilon(self):
        ## Taking into account distance scores for the previous population, 
        ## propose a new epsilon, or default to linear epsilon selection
        
        if self.t==0:
            self.epsilon_t = self.epsilon_0
            return None
        distance_set = self.results_d[-1]
        distance_set_sorted = sorted(distance_set)
        quantile_idx = int(self.alpha*len(distance_set))
        self.epsilon_t = distance_set_sorted[quantile_idx]
        
    def get_p_transition(self):
        step_size = self.p_0 / float(self.T)
        p_t = self.p_0 - step_size * self.t
        self.p_t = p_t
    
class Particle():
    """A self-contained ABC particle that will return an accepted parameter set and its weight."""
     
    def __init__(self,
            theta_set_prev,
            w_set_prev,
            model,
            experiments,
            epsilon_curr,
            prior_set,
            p_transition_curr,
            t,
            N,
            particle_id,
            num_expts_pos,
            num_expts_neg,
            ln_pi_max,
            precalc_media,
            precalc_media_frozensets):
         
        self.theta_set_prev = theta_set_prev
        self.w_set_prev = w_set_prev
        self.model = model
        self.experiments = experiments
        self.prior_set = prior_set
        self.epsilon = epsilon_curr
        self.p_transition = p_transition_curr
        self.t = t
        self.N = N
        self.id = particle_id
        self.num_expts_pos = num_expts_pos
        self.num_expts_neg = num_expts_neg
        self.ln_pi_max = ln_pi_max
        self.precalc_media = precalc_media
        self.precalc_media_frozensets = precalc_media_frozensets
        
        self.num_params = len(prior_set)
         
        self.theta_sampled = None
        self.theta_accepted = None
        self.theta_proposed = None
        self.weight = None
        self.full_results = None
        self.result = None
         
        return None

         
    def propose_theta(self):
        """Sample theta from previous iteration and perturb according to K_t."""
         
        if self.t > 0: 
            ## Get theta by sampling
             
            theta_sampled_idx = np.random.choice(self.N,p=self.w_set_prev)
            self.theta_sampled = self.theta_set_prev[theta_sampled_idx]
             
            ## Perturb theta using perturbation kernel
            self.theta_proposed = self.K_t()
        else:
            ## Sample theta from prior
            self.theta_proposed = self.pi()
         
    def apply_proposed_theta(self):
        """Set constraints on model according to proposed theta."""
        for idx, param_value in enumerate(self.theta_proposed):
            reaction = self.model.theta_rxn_map[idx]
            if param_value == 1:
                reaction.lower_bound = reaction.lb_orig
                reaction.upper_bound = reaction.ub_orig
            else:
                reaction.lower_bound = 0   
                reaction.upper_bound = 0          
     
#     def calculate_weight(self):
#         """Calculate weight from accepted theta value."""
#         if self.t == 0:
#             self.weight=1
#         else:
#             prior_prob = self.pi()
#             w_perturbation_sum = 0
#             sum_data = zip(self.w_set_prev, self.theta_set_prev)
#             for w_prev, theta_prev in sum_data:
#                 contribution = w_prev * self.K_t(theta_prev)
#                 w_perturbation_sum += contribution
#             w_out = prior_prob / w_perturbation_sum
#             self.weight = w_out
#     
#     def calculate_weight_ln(self):
#         """To avoid underflow use log/exp to calulate weight."""
#         if self.t == 0:
#             self.weight=1
#         else:
#             a_j_set = []
#             for idx, theta_prev in enumerate(self.theta_set_prev):
#                 w_j = self.w_set_prev[idx]
#                 a_j = 0
#                 added_thetas = [sum(entry) for entry in zip(self.theta_accepted, theta_prev)]
#                 for param_idx in range(self.num_params):
#                     if added_thetas[param_idx] == 1:
#                         ##Different
#                         a_j += ln(w_j**(1.0/self.num_params)*0.2)
#                     else:
#                         a_j += ln(w_j**(1.0/self.num_params)*0.8)
#                 a_j_set.append(a_j)
#             m = max(a_j_set)
#             a_j_minus_m_set = [a_j - m for a_j in a_j_set]
#             exp_ajmm_sum = sum([exp(ajmm) for ajmm in a_j_minus_m_set])
#             ln_w_i = ln(self.pi()) - m + ln(exp_ajmm_sum)
#             w_i = exp(ln_w_i)
#             self.weight = w_i
    
    def ln_pi(self):
        """find logarithm of the prior for calculating logarithm of weight."""
        
        ln_pi = 0
        for idx, param in enumerate(self.theta_accepted):
            if (param == 1):
                ln_pi += ln(self.prior_set[idx])
            else:
                ln_pi += ln(1-self.prior_set[idx])        
        self.ln_pi = ln_pi
    
    def n_diff_min_and_residuals(self):
        """From the previous accepted theta set, calculate the minimum number 
        of changes from any of them to the accepted theta."""
        
        n_diff_min = self.num_params
        n_diff_list = []
        for theta_prev in self.theta_set_prev:
            n_diff = sum([abs(diff[0]-diff[1]) for diff in zip(theta_prev, self.theta_accepted)])
            n_diff_list.append(n_diff)
            if n_diff < n_diff_min:
                n_diff_min = n_diff
        self.n_diff_res = [n_diff - n_diff_min for n_diff in n_diff_list]
        self.n_diff_min = n_diff_min
    
    def calculate_ln_weight_denominator(self):
        """Calculate the logarithm of the weighted sum of the K_ts for the 
        calculation of weight."""
        
        p = self.p_transition
        
        weighted_sum = 0
        for w_j, n_diff_res in zip(self.w_set_prev, self.n_diff_res):
            weighted_residual = w_j * p**n_diff_res * (1-p)**n_diff_res
            weighted_sum += weighted_residual
        
        try:
            self.ln_weight_denominator = n_diff_res*(ln(p) + ln(1-p)) + ln(weighted_sum)
        except:
            print("This particle is so distant from every other particle that it has underflown, it will take the maximum value of the weights of the other particles in the population")
            self.ln_weight_denominator = None
            
    def calculate_ln_w(self):
        """Calculate ln(weight) for this particle."""
        
        if self.t == 0:
            self.ln_w = 0
        else:
            self.ln_pi()
            self.n_diff_min_and_residuals()
            self.calculate_ln_weight_denominator()
            
            try:
                self.ln_w = self.ln_pi - self.ln_weight_denominator
            except:
                self.ln_w = None
        
  
    
    
#     def pi_rel(self):
#         """Find pi_rel = pi/pi_max using logs to avoid underflow."""
#         ## Return pi_i / pi_max (== 0 if out of float bounds)
#         ln_pi_rel = -self.ln_pi_max
#         for idx, param in enumerate(self.theta_accepted):
#             if (param == 1):
#                 ln_pi_rel += ln(self.prior_set[idx])
#             else:
#                 ln_pi_rel += ln(1-self.prior_set[idx])
#         self.pi_rel = exp(ln_pi_rel)    
#     
#     def calculate_weight_simple_K(self):
#         """Weight calculations are still very tricky - use simplified K_t."""
#         
#         if self.t == 0:
#             self.weight = 1
#         else:
#             self.pi_rel()
#             if self.pi_rel != 0:
#                 denominator = 0
#                 for j, theta_prev in enumerate(self.theta_set_prev):
#                     w_j = self.w_set_prev[j]
#                     sum_diff = 0
#                     theta_diff = [entry[0] - entry [1] for entry in zip(self.theta_accepted, theta_prev)]
#                     for param_idx in range(self.num_params):
#                         if theta_diff[param_idx] != 0:
#                             ## Different
#                             sum_diff += 1
#                     K_t_j = 1.0/(sum_diff + 1)**2
#                     denominator += w_j*K_t_j
#                 w_i = self.pi_rel/denominator
#             else:
#                 w_i = 0        
#             self.weight = w_i
        
    def find_accepted_theta(self, debug = True, use_ln = True, use_simple_K=False):
        """Find accepted theta and return it with weight."""
         
        count_attempts = 0
        print("\nFinding particle:".format(self.t, self.id))
        
        while True:
            count_attempts += 1
            if debug:
                sys.stdout.write("Particle {} (epsilon = {}) ... ".format(count_attempts, self.epsilon))    
                sys.stdout.flush()
            self.propose_theta()
            self.apply_proposed_theta()
            self.conduct_experiments() 
            if debug:
                sys.stdout.write("\rresult = {:.3f}".format(self.result))
                if self.result != 2:
                    sys.stdout.write(" after {}/{} tests".format(self.num_tests_checked, self.num_tests_total))
                sys.stdout.flush()
             
            ## If model is close enough to experiments, accept and return theta and calculated weight
            if self.result < self.epsilon:
#                 print(" - Attempt successful, returning accepted theta ...")
                self.theta_accepted = self.theta_proposed
                self.calculate_ln_w()
                return self.theta_accepted, self.ln_w, self.result
 
    def K_t(self, theta_previous = None):
        """Perturbation function!
         
        Return a perturbed theta, theta_perturbed, from theta_sampled according to a simple distribution."""
         
        if theta_previous is None:
            ## Perturb the parameter set self.theta_sampled
            theta_perturbed = np.zeros(len(self.theta_sampled))
            for idx, param_sampled in enumerate(self.theta_sampled):
                param_perturbed = param_sampled
                if random() < self.p_transition:
                    if param_sampled == 1:
                        param_perturbed = 0
                    if param_sampled == 0:
                        param_perturbed = 1
                theta_perturbed[idx] = param_perturbed 
            return theta_perturbed
        else:
            ## Output the probability of a transition to self.theta_sampled from theta_previous    
            added_thetas = self.theta_accepted + theta_previous
            probability = 1.0
            ln_p = 0
            for x in added_thetas:
                if x == 1:
                    ln_p += ln(self.p_transition)
                    probability *= self.p_transition
                else:
                    ln_p += ln(1-self.p_transition)
                    probability *= 1 - self.p_transition
            return probability, ln_p
     
    def pi(self):
        """Return probability of getting a particular theta by sampling at random from the prior.
        If theta is not specified, sample randomly from the prior."""
         
        if self.theta_accepted is not None:
            ## Calculate the chances of picking the accepted parameter set at random from the prior
            probability = 1.0 
            param_prior_pairs = zip(self.theta_accepted, self.prior_set)
            for param, prior in param_prior_pairs:
                if param == 1:
                    probability *= prior
                else:
                    probability *= 1 - prior
            return probability
        else:
            ## Sample from the prior
            theta_proposed = np.zeros(self.prior_set.size)
#             print theta_proposed, "\n"
            prior_iter = np.nditer(self.prior_set)
            for idx, param_prior in enumerate(prior_iter):
#                 print param_prior
                param_value = np.random.binomial(1,param_prior)
#                 print(" - {} ({})".format(param_value, type(param_value)))
                theta_proposed[idx] = param_value
            return theta_proposed
 
    def conduct_experiments(self, epsilon=None):
        """Conduct experiments for model in current state and return 1 - (balanced accuracy)."""
         
#         ## Check that the model can run without changes, else return 1
#         if self.model.opt() <= 0:
#             print("Failed on original media")
#             self.result = 1
#             return None
        
        self.num_tests_checked = 0
        self.num_tests_total = 0
        
        ## Must run on all common media before experimental testing
        for precalc_medium in self.precalc_media:
            self.model.set_medium(precalc_medium)
            if self.model.opt() <= 0:
                self.result = 2
#                 print("Failed on precalc media")
                return None
         
#         print"No failure, continuing with test"
         
        ## Check all genotypes for presence in model and create a list of valid experiments
        valid_experiments = []
        num_pos_remaining = 0
        num_neg_remaining = 0
        gene_ids_in_model = self.model.get_relevant_gene_ids()
        for idx, expt in enumerate(self.experiments):
            all_genes_present = True
            for gene in expt.genotype:
                if gene not in gene_ids_in_model:
                    print("Expt {}: '{}' not present in model ...".format(idx+1, gene))
                    all_genes_present = False
                    break
            if all_genes_present:
                valid_experiments.append(expt)
                if expt.result == 1:
                    num_pos_remaining += 1
                else:
                    num_neg_remaining += 1
         
        self.num_tests_total = len(valid_experiments)
        print("{} valid experiments ...".format(len(valid_experiments)))
        
        num_failed_tests = 0
        num_succeeded_tests = 0 
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        
        running_results = ResultSet(0,0,0,0)
        min_dist = 'unk'
        
        for idx, experiment in enumerate(valid_experiments):
            
            
#             print("Testing experiment {}".format(idx))
            expt_result, tp_add, tn_add, fp_add, fn_add = experiment.test(self.model, self.precalc_media_frozensets)
            running_results.tp += tp_add
            running_results.tn += tn_add
            running_results.fp += fp_add
            running_results.fn += fn_add
            
            if experiment.result == 1:
                num_pos_remaining -= 1
            else:
                num_neg_remaining -= 1
            
            
            if expt_result == 1:
                num_succeeded_tests += 1
            if (expt_result != 1) and (expt_result != -2):
                num_failed_tests += 1

            ## If it is impossible to get below epsilon with all following 
            ## tests correct, stop early
            if fn_add + fp_add > 0:
                check_results = deepcopy(running_results)
                check_results.tp += num_pos_remaining 
                check_results.tn += num_neg_remaining
                min_dist = 1 - check_results.balanced_accuracy()
                if min_dist > self.epsilon:
#                     print("\nMin dist too high ({}), aborting".format(min_dist))
                    break
            
            print("{}: current bal.acc. = {}".format(idx, running_results.balanced_accuracy()))
            
        self.num_tests_checked = num_succeeded_tests + num_failed_tests

        self.full_results = deepcopy(running_results)
        distance = 1.0 - running_results.balanced_accuracy() 
        self.result = distance
        return None
    



       
def get_reaction(model, reaction_id):
    """Get reaction object from model using reaction ID."""
    
    try:
        return model.reactions.get_by_id(reaction_id)
    except:
        return None  
    
def create_medium(medium_file):
    """
    Import file containing a model medium and return a relevant dictionary.
    """
    
    f_in = open(medium_file, 'r')
    
    medium_dict = {}
    
    for line in f_in:
        if line[0] != "#":
            try:
                exch_id, stoichiometry = line.split("\t")
            except:
                print line
                continue
            
            medium_dict[exch_id] = float(stoichiometry)
    
    return medium_dict

def import_experiments(model, objective_name = "Biomass_BT_v2"):

    objective = model.reactions.get_by_id(objective_name)

    ## Import media
    
    media_shelf_file = '/Users/wbryant/work/BTH/analysis/fba/media.data'
    media_shelf = shelve.open(media_shelf_file)
    media = {}
    media['min'] = media_shelf['min']
    media['tyg'] = media_shelf['tyg']
    media_shelf.close()

    ## Create list of experiments for doing model assessment
    
    essentiality_data_file = '/Users/wbryant/work/BTH/analysis/gene_essentiality/essentiality_data_complete.csv'
    experiments = []
    expts_in = open(essentiality_data_file,'r')
    
    for line in expts_in:
        if line[0] != "#":
            details = line.split("\t")
            medium_base = details[1]
            c_source = details[2]
            genotype = details[3].split(" and ")
            result = details[4]
            
            c_sources = []
            
            ## Add carbon source to medium
            if c_source == "-":
                c_sources = []
            elif 'and' in c_source:
                c_sources = c_source.split(' and ')
            else:
                c_sources = [c_source]
            #expt_medium = media[medium_base]
            expt_medium = deepcopy(media[medium_base])
            for source in c_sources:
                c_source_exchange = "EX_" + source + "(e)"
                expt_medium[c_source_exchange] = -100
            
            
            experiment = Experiment(details[0], expt_medium, result, genotype, objective)
            
            experiments.append(experiment)
    
    return experiments


class Experiment():
    """
    The details of a real experiment formatted for testing a constraint-based model.
    """
    
    def __init__(self, expt_id, medium, result, genotype = None, objective = None, old_expt = None):
        if old_expt:
            self.id = old_expt.id
            self.medium = old_expt.medium
            self.result = old_expt.result
            self.genotype = old_expt.genotype
            self.objective = old_expt.objective
        else:
            self.id = expt_id
            self.medium = medium
            try: 
                float(result)
                self.result = result
            except:
                if result[0] == 'y':
                    self.result = 1
                else:
                    self.result = 0
            self.genotype = genotype or None
            self.objective = objective or None
        
        self.medium_components = frozenset([component for component in self.medium])
    
    def test(self, ec_model, precalc_frozensets = None):
        """
        Calculate whether ec_model is consistent with this experiment.
        """  
        
        precalc_frozensets = precalc_frozensets or []
        
        ec_model.set_medium(self.medium)
        if self.objective:
            ec_model.change_objective(self.objective)
        
        change_to_model = ec_model.set_genotype(self.genotype)
        if (change_to_model == 0) & (self.medium_components in precalc_frozensets):
#             sys.stdout.write(" -> Precalculated medium, no change in model")
#             sys.stdout.flush()
            return 1, 1, 0, 0, 0 
        elif change_to_model == -1:
            ## Gene is not in model.
            return -2,0,0,0,0
        
        model_growth = ec_model.opt()
        ec_model.unset_genotype()
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        
        if (model_growth > 0) and (self.result > 0):
            tp = 1
            expt_result = 1
        elif (model_growth == 0) and (self.result == 0):
            expt_result = 1
            tn = 1
#         elif (model_growth < 0):
#             expt_result = -1            
        else:
            expt_result = 0
            if self.result > 0:
                fn = 1
            else:
                fp = 1
            
            
#         print("%3s (%d): %d" % (self.id, self.result, answer))
        return expt_result, tp, tn, fp, fn
        
def create_extended_model(model_file, objective_id = 'Biomass_BT_v2'):
    """Take an ArrayBasedModel and convert to an Extended_Cobra_Model."""
    
    print("Creating extended COBRA model")
    model = create_cobra_model_from_sbml_file(model_file)
    ecm_model = ExtendedCobraModel(model)
    ecm_model.set_solver()
    ecm_model.set_medium()
    print("done.")
    try:
        ecm_model.change_objective(objective_id)
    except:
        print("Objective could not be set ...")
    return ecm_model    

class ExtendedCobraModel(ArrayBasedModel):
    """Additional functionality for helping to use COBRA model."""
    
    def set_biomass_id(self, biomass_id):
        """Make a note of the biomass equation for setting the objective."""
        
        self.biomass_id = biomass_id
    
    def set_objective_to_biomass(self):
        try:
            self.change_objective(self.biomass_id)
        except:
            print("Objective could not be set to biomass: bimoass_id not set?")
            
    def get_relevant_gene_ids(self):
        """Get a list of gene IDs for all genes implicated in non-zero flux reactions."""
        
        relevant_genes_list = []
        for reaction in self.reactions:
            if (reaction.lower_bound < 0) &\
               (reaction.upper_bound > 0):
                
                gene_ids = [gene.id for gene in reaction.genes]
                relevant_genes_list.extend(gene_ids)
        relevant_genes_list = list(set(relevant_genes_list))
        return relevant_genes_list
    
    def split_rxn_by_enzymes(self, reaction_id, enzyme_limit = 10):
        """Take a reaction from the model and add new reactions for each potential enzyme.
        
        Do not split the reaction if the number of predicted enzymes is > enzyme_limit.
        
        Original reaction is retained as a 'hypothetical' (unidentified) enzyme.
        
        Return list of IDs created for these enzrxns.
        """
        
        rxn = self.reactions.get_by_id(reaction_id)
        gpr_list = self.enzymes_as_gprs(rxn)
        num_rxns = len(gpr_list)
        
        
        if len(gpr_list) > enzyme_limit:
#             print("Too many enzymes ({}) for reaction {}.".format(len(gpr_list),reaction_id))
            return [rxn.id], None
#         elif len(gpr_list) == 1:
# #             print("Single enzyme for reaction {}.".format(rxn.id))
#             return[rxn.id], rxn.id
        elif len(gpr_list) == 0:
#             print("Spontaneous reaction: {}.".format(rxn.id))
            return[rxn.id], None
        
        enzyme_index = 0
        enzrxn_id_list = []
        
#         print("{} enzymes for reaction {}.".format(len(gpr_list), rxn.id))
        for gpr in gpr_list:
            enzyme_index += 1
            rxn_dup = deepcopy(rxn)
            rxn_dup.gene_reaction_rule = gpr
            rxn_dup.id = rxn_dup.id + "_enz" + str(enzyme_index)
            enzrxn_id_list.append(rxn_dup.id)
            self.add_reaction(rxn_dup)

        rxn.gene_reaction_rule = ''
        
        return enzrxn_id_list, rxn.id
        
        
        
    def enzyme_to_gpr(self, enzyme_string, warn_if_absent = False, reject_if_absent = False):
        """Convert enzyme string to GPR string.
        Enzyme string should be a comma separated list of valid 
        gene loci.
        """
        
        gene_ids = enzyme_string.split(",")
        if reject_if_absent or warn_if_absent:
            genes_model_list = []
            for gene in self.genes:
                genes_model_list.append(gene.name)
        reject = False 
        for gene_id in gene_ids:
            if reject_if_absent or warn_if_absent:
                if gene_id not in genes_model_list:
                    if warn_if_absent:
                        print("Gene {} not currently in model".format(gene_id))
                    if reject_if_absent:
                        reject = True
        if reject == True:
            print(" - Gene(s) not in current model, conversion aborted.")
            return None
        else:
            gpr_string =\
                "( " +\
                " and ".join(gene_ids) +\
                " )"
            return gpr_string
    
    
    def enzymes_as_gprs(self, reaction):
        """Expand the GPR for the given reaction and list all enzymes as separate GPRs.
        """
        enzymes = gene_parser(reaction.gene_reaction_rule)
        gpr_list = []
        for enzyme in enzymes:
            gpr_list.append(self.enzyme_to_gpr(enzyme))
        return gpr_list
    
  
    def halt_reactions(self, reaction_list):
        """Remove each reaction in list (of IDs), ignoring reactions that break the model.
        
        USE WITH CARE! Will be dependent on current medium and objective.
        """
        ## Check that model runs
        if self.opt() <= 0:
            print("Model does not currently run ...")
            return None
        self.halted_reactions = []
        print("Reactions halted ...")
        for reaction_id in reaction_list:
            reaction = self.reactions.get_by_id(reaction_id)
            old_lb = 0
            old_ub = 0
            if reaction.lower_bound < 0:
                old_lb = reaction.lower_bound
                reaction.lower_bound = 0
                self.halted_reactions.append((reaction, 'lb', old_lb))
            if reaction.upper_bound > 0:
                old_ub = reaction.upper_bound
                reaction.upper_bound = 0
                self.halted_reactions.append((reaction, 'ub', old_ub))
            if self.opt() <= 0:
                reaction.lower_bound = old_lb
                reaction.upper_bound = old_ub
                print("(Not {} - {}, {})".format(reaction.id, old_lb, old_ub))
            else:
                print(" - {} ({})".format(reaction.id, self.opt()))
    
    def unhalt_reactions(self):
        """Reinstate reactions halted by 'halt_reactions_limited'."""
        try:
            halted_list = self.halted_reactions
        except:
            print("No halted reactions ...")
            return None
        for halted_rxn in halted_list:
            if halted_rxn[1] == 'ub':
                halted_rxn[0].upper_bound = halted_rxn[2]
            else:
                halted_rxn[0].lower_bound = halted_rxn[2]   
        self.halted_reactions = []
    
    def set_genotype(self, genotype):
        """
        Run delete_model_genes on model.  Return False if no change is made to the model
        """
        
        try:
            ## CODE FROM COBRAPY: to access reaction deletions
            # Allow a single gene to be fed in as a string instead of a list.
            if not hasattr(genotype, '__iter__') or \
                    hasattr(genotype, 'id'):  # cobra.Gene has __iter__
                genotype = [genotype]
        
            if not hasattr(genotype[0], 'id'):
                if genotype[0] in self.genes:
                        tmp_gene_dict = dict([(x.id, x) for x in self.genes])
                else:
                    # assume we're dealing with names if no match to an id
                    tmp_gene_dict = dict([(x.name, x) for x in self.genes])
                genotype = [tmp_gene_dict[x] for x in genotype]    
            # Make the genes non-functional
            for x in genotype:
                x.functional = False
        except:
            return -1
       
        
        knocked_out_reactions = find_gene_knockout_reactions(self, genotype)
        if len(knocked_out_reactions) == 0:
            return 0
        else:
            delete_model_genes(self, genotype)
            return 1
        
    def unset_genotype(self):
        """
        Run undelete_model_genes on model.
        """    
        undelete_model_genes(self)
        
    
    def set_solver(self, solver_string = 'glpk'):
        self.solver = solver_string
        
    def opt(self, new_objective = None):
        if new_objective:
            self.change_objective(new_objective)
        
        self.optimize(solver=self.solver)
        
        if self.solution.f:
            if self.solution.f < 0:
                return 0
            else:
                return self.solution.f
        else:
            return 0
    
    def set_medium(self,medium_dict=None):
        """
        Set exchanges of all "EX_" reactions to zero, or what is in medium_dict.
        medium_dict: key = reaction ID, value = new lower bound.
        """
        exchange_rxn_dict = {}
        for reaction in self.reactions:
            if reaction.id[0:3] == "EX_":
                ## Exchange reaction, so if lb < 0, print
                reaction.lower_bound = 0
                exchange_rxn_dict[reaction.id] = reaction
        if medium_dict: 
            for component in medium_dict:
                if component in exchange_rxn_dict:
                    exchange_rxn_dict[component].lower_bound = medium_dict[component]
                elif 'EX_' + component + '(e)' in exchange_rxn_dict:
                    exchange_rxn_dict['EX_' + component + '(e)'].lower_bound = medium_dict[component]
#                 else:
#                     print("Reaction ID '%s' not found, ignoring ..." % component)
        self.medium_dict = medium_dict
        
        
    def show_medium(self):
        """
        Find all exchange lower bounds < 0 and print to screen.
        """
        for reaction in self.reactions:
            if reaction.id[0:3] == "EX_":
                ## Exchange reaction, so if lb < 0, print
                
                if reaction.lower_bound < 0:
                    reaction_identifier = reaction.name + " (" + reaction.id + ")"
                    print("%45s %5.0f %5.0f" % (reaction_identifier, reaction.lower_bound, reaction.upper_bound))
    
    def add_source_exchange_reaction(self, compound):
        """Create a source reaction for a metabolite to allow the take-up of 
        that compound without a transport reaction."""
        rxn_name = "{} exchange".format(compound.name) 
        rxn_id = "EX_{}".format(compound.id)
        new_rxn = Reaction(rxn_id)
        new_rxn.name = rxn_name
        new_rxn.lower_bound = -1
        new_rxn.upper_bound = 1000
        new_rxn.add_metabolites({compound: -1.0})
        self.add_reaction(new_rxn)
        print("Reaction '{}' added to model".format(new_rxn.name))
    
    def find_blocked_reaction_components(self, show_not_blocked=False, tol=1e-10):
        """Find those components of the objective reaction that cannot be created in the given conditions."""
        try:
            self.blocked_metabolites
        except:
            self.blocked_metabolites = []
        objective_reactions = []
        for reaction in self.reactions:
            if reaction.objective_coefficient != 0:
                objective_reactions.append(reaction)
        if len(objective_reactions) == 0:
            print("No objective set")
            return None
        elif len(objective_reactions) > 1:
            print("Too many reactions in objective")
            return None
        else:
            subject_reaction = objective_reactions[0]
        for metabolite in subject_reaction.reactants:
            test_rxn = Reaction('test')
            test_rxn.add_metabolites({metabolite:-1})
            self.add_reaction(test_rxn)
            self.change_objective('test')
#             if (self.opt() < tol) or show_not_blocked:
#                 print("{:50s}: {}".format(metabolite.name + " (" + metabolite.id + ")", self.opt()))
            if self.opt() < tol:
                self.blocked_metabolites.append(metabolite)
            test_rxn.remove_from_model()
        self.change_objective(subject_reaction.id)
    
    def reset_blocked_metabolites(self):
        self.blocked_metabolites = []
        
    def source_blocked_biomass_components(self, blocked_met_file, biomass_id = None, general_metabolite_ids = None, xml_output_file = '\Users\wbryant\work\MTU\MTU_SEED_unblocked.xml'):
        """Add sources for all specific metabolites blocked in the biomass.
        Adds exchange and uncatalysed transport reactions.
        
        N.B. THIS DEPENDS ON THE CURRENT MEDIUM
        
        general_metabolites is a list of 'general' metabolites 
        included in the biomass, such as 'protein'.  The reactions producing these metabolites
        will be individually tested and sourced if blocked.""" 
        
        general_metabolite_ids = general_metabolite_ids or None
        if (not hasattr(general_metabolite_ids, '__iter__')) and (general_metabolite_ids is not None):
            general_metabolite_ids = [general_metabolite_ids]
        
        if biomass_id:
            self.set_biomass_id(biomass_id)
        
        self.set_objective_to_biomass()
        self.find_blocked_reaction_components()
        
        blocked_reaction_ids = []
        for metabolite_id in general_metabolite_ids:
            try:
                metabolite = self.metabolites.get_by_id(metabolite_id)
                if metabolite not in self.blocked_metabolites:
                    print("Metabolite '{}' is either not blocked, or is not in the biomass equation, skipping ...".format(metabolite.id))
                    continue
                else:
                    self.blocked_metabolites.remove(metabolite)
            except:
                print("Metabolite ID '{}' could not be found, skipping ...".format(metabolite_id))
                print("Some typical metabolite IDs are:\n")
                met_ids = [metabolite.id for metabolite in self.metabolites]
                met_ids = met_ids[0:4]
                for met_id in met_ids:
                    print met_id
                continue
            
            ## Check that there is only one producing reaction
            if len(metabolite.reactions) != 2:
                print("Metabolite '{}' did not have a single generating reaction, ignoring ...".format(metabolite.id))
            else:
                ## Check it is in biomass
                this_met_reaction_ids = [reaction.id for reaction in metabolite.reactions]
                if self.biomass_id not in this_met_reaction_ids:
                    print("Metabolite '{}' is not a biomass component, ignoring ...".format(metabolite.id)) 
                else:
                    this_met_reaction_ids.remove(self.biomass_id)
                    blocked_reaction_ids.append(this_met_reaction_ids[0])
                                   
        for reaction_id in blocked_reaction_ids:
            reaction = self.reactions.get_by_id(reaction_id)
            self.change_objective(reaction)
            self.find_blocked_reaction_components()
        
        print("Adding source and transporter for the following metabolites:")
        for metabolite in self.blocked_metabolites:
            print("- {} ({})".format(metabolite.name, metabolite.id))
        
        ## For now don't bother with boundary metabolites, just exchange of 
        ## extracellular metabolites and transport
        
        ## N.B. if the EXC metabolite differs in ID by more than just the final 
        ## character (i.e. 'e'/'c'), this will not find it. 
        
        f_out = open(blocked_met_file, 'w')
        ex_reaction_ids = []
        for metabolite_c in self.blocked_metabolites:
            
            ### Check for extracellular metabolite and create if not existent
            metabolite_c_id = metabolite_c.id
            metabolite_e_id = metabolite_c_id[:-1] + 'e'
            met_e_exists = False
            try:
                metabolite_e = self.metabolites.get_by_id(metabolite_e_id)
                met_e_exists = True
            except:
                metabolite_e = Metabolite(metabolite_e_id,
                    formula=metabolite_c.formula,
                    name=metabolite_c.name,
                    compartment='e')
            
            ### Create transport reaction if one does not exist
            tr_rxn_exists = False
            if met_e_exists:
                for reaction in self.reactions:
                    if ((metabolite_c in reaction.reactants) and (metabolite_e in reaction.products))\
                    or ((metabolite_e in reaction.reactants) and (metabolite_c in reaction.products)):
                        print("Transport reaction for '{}' already exists ('{}') ...".format(
                                        metabolite_c.name, reaction.id))
                        tr_rxn_exists = True
                        continue    
            if not tr_rxn_exists:
                tr_reaction_id = "TR_{}".format(metabolite_c_id[:-2])
                tr_reaction_name = "{} Transport".format(metabolite_c.name)
                transport_reaction = Reaction(tr_reaction_id)
                transport_reaction.add_metabolites({
                    metabolite_e: -1.0,
                    metabolite_c: 1.0
                    })
                self.add_reaction(transport_reaction)
                
            ### Create exchange reaction if one does not exist
            ex_reaction_id = "EX_{}".format(metabolite_c_id[:-2])
            try:
                ex_reaction = self.reactions.get_by_id(ex_reaction_id)
                print("Exchange reaction '{}' already exists ...".format(ex_reaction_id))
            except:
                ex_reaction_name = "{} Exchange".format(metabolite_c.name)
                exchange_reaction = Reaction(ex_reaction_id)
                exchange_reaction.add_metabolites({
                    metabolite_e: -1.0
                    })
                self.add_reaction(exchange_reaction)
            
            ### Log to 'blocked_met_file' exchange reactions for the 
            ### metabolites that have been sourced this way
            ex_reaction_ids.append(ex_reaction_id)
            f_out.write("{}\n".format(ex_reaction_id))
        f_out.close()
        
        ### Test addition of exchange reactions for blocked metabolites
        self.set_objective_to_biomass()
        blocked_test = self.opt()
        if blocked_test != 0:
            print("Model runs without the addition of sources for blocked metabolites - there are no blocked metabolites")
        else:
            unblocked_medium = self.medium_dict
            for reaction_id in ex_reaction_ids:
                unblocked_medium[reaction_id] = -1
            self.change_objective(unblocked_medium)
            unblocked_opt = self.opt()
            print("With unblocked metabolites, biomass growth rate is {}".format(unblocked_opt))
            
        ## Output model to XML file
        write_sbml_model(self, xml_output_file)
             
if __name__ == '__main__':
    pass