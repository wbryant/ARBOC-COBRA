'''
ABC-SMC for COBRA models v0.2
Copyright 2017 William A. Bryant

ABSTRACT
========

This module undertakes Approximate Bayesian Computation using a Sequential
Monte Carlo approach (ABC-SMC) to assess proposed additions and uncertain 
enzyme/reaction relationships in a COBRA constraint-based metabolic model, and
find a combination of additions and removals of reactions and enzyme/reaction 
relationships to and from an original model that produces a model with 
maximum parsimony with a particular set of experimental observations.   

SUMMARY
=======

A specific reaction can be catalysed by one or more (or no) enzymes, which are
either known or unknown.  The assertion that a particular reaction is catalysed
by a particular enzyme, a Reaction/Enzyme Link (or REL), is tested in this
algorithm.  The algorithm requires as input a full model including all known and
proposed RELs, a set of experiments to determine the congruence of model and 
experiment, and a confidence estimation for the existence of each REL.

The algorithm requires a set of parameters to estimate.  In this case the
existence of a particular REL is cast as a boolean parameter and prior confidence
for the existence of that REL is taken from the confidence estimation provided.
The existence of all RELs with a confidence less than 1 are included as 
parameters in the ABC-SMC and will be either included in or excluded from each
proposed model created by the algorithm.

The algorithm will propose a large number of alternative models, calculating
for them the congruence with the experimental data provided.  Models will be 
accepted if their distance measure form the experimental data falls below a
specified value (epsilon).  Once N models have been found that are accepted,
epsilon is reduced and new models are proposed by selecting a random accepted 
model from the previous step and perturbing it.  At each step both epsilon and
the perturbation strength are reduced until a set of N models is settled upon,
having a low epsilon.

The output is therefore a set of N models, each with a high congruence with the
observed experimental data.  The confidence for the existence of each REL in 
reality can then be calculated as the proportion of these models containing 
that REL.      

For further information on the algorithm and formulation of inputs, see the 
accompanying paper.

LICENSE
=======

This program is free software, distributed under the terms of the GNU GPL:
you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

INSTALLATION
============

This module does not need to be installed; if the directory containing this
file is in the PYTHONPATH, it can be imported using:

`import abcsmc`

PREREQUISITES
=============

This module has been tested with Python 2.7

Required non-standard Python packages (tested versions specified in brackets):

- COBRApy (0.4.0b4)
- NumPy (1.10.4)
- LibSBML (5.12.1) for SBML import and export
- Pyparsing (2.0.3) for GPR <--> list of enzymes translation

This package also requires the installation of a compatible linear optimisation
program.  Testing has shown Gurobi to be up to 30 times faster than the free
GLPK program, and is available for free with an academic license when used 
solely for the purposes of research.

EXAMPLE USAGE
=============
Most C{AbcProblem} arguments are optional (or conditionally optional), the
following are required:

- An SBML file containing the model including predicted RELs
- An experiment file [1]_ containing results of gene essentiality analyses
- if the 'include_all' variable is set to False (default) then a list 
  of candidate RELs for inclusion/exclusion must be supplied (e.g.
  'example/candidate_rels.tsv')

Default values for variables are set for the purposes of testing, rather than
full analysis.  This example uses non-default values for several optional
variables.    

Create problem object:

.. code-block:: python

    abc_problem = AbcProblem(
        model_file='example/model.xml',
        experiments_file='example/experiments.xml',
        prior_file='example/candidate_rels.tsv',
        particles_per_population=20,
        num_populations_max=6
    )
    


.. [1] The experiment file has a specific format, as set out in the file
       'example/experiments.tsv'

@author: wbryant
'''
from cobra.flux_analysis import find_blocked_reactions
from cobra.core import ArrayBasedModel
from cobra.io.sbml import create_cobra_model_from_sbml_file 
from cobra.manipulation.delete import delete_model_genes, undelete_model_genes, find_gene_knockout_reactions
from local_gene_parser import gene_parser
import numpy as np
from numpy import log as ln
from copy import deepcopy
from random import random
from math import sqrt, log10, ceil, exp
import shelve
from collections import Counter
import sys, re
 
class AbcProblem():
    
    def __init__(self,
            model=None,
            model_file=None,
            prior_dict = None,
            prior_file=None,
            experiments=None,
            experiments_file=None,
            particles_per_population=None,
            num_populations_max=None,
            rel_priors = None,
            default_prior_value = 0.99,
            enzyme_limit = 10,
            objective_name = 'Biomass',
            epsilon_0 = 0.5,
            epsilon_T = 0.285,
            p_0 = 0.2,
            alpha = None,
            include_all=False,
            biomass_id=None,
            original_model_rxn_ids=None,
            solver=None,
            prior_multiplier=None,
            test_type=None,
            test_model_pre_abc=None,
            include_fn_reactions=None):

        """Set up a COBRA model for ABC-SMC improvement
        
        num_populations = number of iterations to run through;
        model = ExtendedCobraModel to be tested;
        get_epsilon = a function of t to calculate epsilon;
        get_p_transition = a function of t to calculate p_transition;
        abc_reactions: list of reaction IDs to be included;
        prior_dict: dictionary of predefined prior values for specific reactions (rxn ID is key);
        rel_priors: a list with entries so - ['base_rxn_id',[gene_list],prior];
        default_prior: default prior;
        include_all: if True, include all non-exchange reactions in the ABC
        
        """
        self.t = 0
        self.N = particles_per_population or 10
        self.T = num_populations_max or 5
        solver = solver or 'glpk'
        biomass_id = biomass_id or 'Biomass0'
        experiments_file = experiments_file or None
        self.test_type = test_type or None
        test_model_pre_abc = test_model_pre_abc or False
        if include_fn_reactions is None:
            include_fn_reactions = True
        
        print("Loading model ...")
        if model_file:
            self.model = create_extended_model(model_file,biomass_id,solver=solver)
        elif model:
            self.model = deepcopy(model)
        else:
            print("Model not specified, exiting ...")
            sys.exit(1)    
        
        print("Loading experiments ...")
        if experiments_file:
            self.experiments = import_expt_data(self.model, objective_id = biomass_id, data_file=experiments_file)
        elif experiments is not None:
            self.experiments = experiments
        else:
            print("Experiments not specified, exiting ...")
            sys.exit(1)
        
        
        if original_model_rxn_ids:
            ## Test original model
            print("Testing original model ...")
            model_orig = deepcopy(self.model)
            for reaction in model_orig.reactions:
                if reaction.id not in original_model_rxn_ids:
                    reaction.remove_from_model()
            model.repair()     
            original_results, _ = conduct_experiments(model_orig, self.experiments)
            print(" => Results:\n")
            original_results.stats()    
                 
        if test_model_pre_abc:
            ## Test full model
            print("Testing full model ...")
            initial_results, _ = conduct_experiments(self.model, self.experiments)
            print(" => Results:\n")
            initial_results.stats()        
        
        
        
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
        abc_reactions = []
        
        ## Separate reactions specified in rel_priors, to allow specification 
        ## of individual REL priors
        if rel_priors:
            for rxn_id, _, _ in rel_priors:
                abc_reactions.append(rxn_id)
        
        ## Set particle ID format according to how many particles there are
        id_length = str(int(ceil(log10(self.N))))
        self.id_format = "{:0" + id_length + "d}"
        
        ## Initialise results variables, for storing run results
        self.results_theta = []
        self.results_w = []
        self.results_d = []
        self.thetas_selected = []
        self.intermediate_theta_dict = {}
        
        if self.include_all:
            print("Including all non-exchange reactions")
            for reaction in self.model.reactions:
                if not reaction.id.startswith('EX_'):
                    abc_reactions.append(reaction.id)
            print("\t{} reactions in total".format(len(abc_reactions)))
        print("Pre-priors: {} reactions included".format(len(abc_reactions)))
        print("Loading prior values ...")
        model_reaction_ids = [reaction.id for reaction in self.model.reactions]
        if not prior_dict:
            prior_dict = {}
            if prior_file:
                f_in = open(prior_file, 'r')
                for line in f_in:
                    details = line.strip().split("\t")
                    if details[0] not in model_reaction_ids:
                        print("At least one reaction ID is not recognised ('{}')".format(details[0]))
                        for reaction_id in model_reaction_ids[0:9]:
                            print reaction_id
                        sys.exit(1)
                    try:
                        prior_dict[details[0]] = float(details[1])
                    except:
                        prior_dict[details[0]] = self.default_prior_value
                    #print("\t'{}':\t'{}'".format(details[0],prior_dict[details[0]]))
                f_in.close()
        
        print("{} in prior_dict".format(len(prior_dict)))
        
        for rxn_id, _ in prior_dict.iteritems():
            abc_reactions.append(rxn_id)
        abc_reactions = list(set(abc_reactions))
        
        print("Post-priors: {} reactions included".format(len(abc_reactions)))
        
        ## Split all included reactions into individual enzyme/reaction pairs 
        ## and assign prior values according to prior_dict
        counter = loop_counter(len(abc_reactions),'Splitting ABC reactions')
        for rxn_id in abc_reactions:        
            enzrxn_ids, non_enz_rxn_id = self.model.split_rxn_by_enzymes(rxn_id, enzyme_limit)
            num_enzrxns = len(enzrxn_ids)
            if num_enzrxns == 0:
                num_enzrxns = 1
            if rxn_id in prior_dict:
                prior_value = prior_dict[rxn_id]/sqrt(float(num_enzrxns))
            else:
                #prior_value = default_prior_value/sqrt(float(num_enzrxns))
                prior_value = default_prior_value
            if enzrxn_ids:
                for enzrxn_id in enzrxn_ids:
                    prior_dict[enzrxn_id] = prior_value
                if non_enz_rxn_id:
                    prior_dict[non_enz_rxn_id] = prior_value/100.0
            else:        
                if non_enz_rxn_id:
                    prior_dict[non_enz_rxn_id] = self.default_prior_value
            counter.step()
        counter.stop()


        if rel_priors:
            ## Apply belief about rxn/gene/enzyme relationships
            print("Adding REL priors ...")
            prior_dict = self.set_rel_priors(rel_priors, prior_dict)
            print("REL priors added.")
        
        abc_running_total = len(prior_dict)
        print(" => {} RELs included".format(abc_running_total))

        ## Create original lb/ub vars for reference
        for rxn in self.model.reactions:
            rxn.lb_orig = deepcopy(rxn.lower_bound)
            rxn.ub_orig = deepcopy(rxn.upper_bound)
        
        if include_fn_reactions:
            ## FIND RELS PRODUCING FALSE NEGATIVE PREDICTIONS AND ADD TO ABC-SMC        
            ## For each gene in experiments, determine set of reactions that are stopped by gene deletion.
            prior_fn_value = 0.95  
            counter = loop_counter(len(self.experiments),"Finding inessential genes incorrectly predicted as essential")
            expt_num = 0
            for expt in self.experiments:
                expt_num += 1
                expt_no_genotype = deepcopy(expt)
                expt_no_genotype.genotype = []
                if len(expt.genotype) == 1:
                    a, b, c, d, fn = expt.test(self.model)
                    if fn:
                        ## Which reaction is implicated?
                        ko_rxns = self.model.set_genotype(expt.genotype, return_ko_rxns=True)
                        self.model.unset_genotype()
                        if len(ko_rxns) == 0:
                            pass
                        elif len(ko_rxns) == 1:
                            ko_rxns_essential = [ko_rxns[0].id]
                        else:
                            ## Find model-breaking reaction(s)
                            ko_rxns_essential = []
                            for rxn in ko_rxns:
                                rxn.lower_bound = 0
                                rxn.upper_bound = 0
                                _, _, _, _, fn_rxn = expt_no_genotype.test(self.model)
                                if fn_rxn:
                                    ko_rxns_essential.append(rxn.id)
                                rxn.lower_bound = rxn.lb_orig
                                rxn.upper_bound = rxn.ub_orig
                        for rxn_id in ko_rxns_essential:
                            try:
                                ## Already in prior dict? Amend if new prior lower
                                prior_val = prior_dict[rxn_id]
                                if prior_val > prior_fn_value:
                                    prior_dict[rxn_id] = prior_fn_value
                            except:
                                ## Add to prior dict
                                prior_dict[rxn_id] = prior_fn_value   
                counter.step()
            counter.stop()
             
            print(" => {} RELs added".format(len(prior_dict)-abc_running_total))
            
            
        ## All beliefs about reactions included in the model are now in prior_dict.
        ## Any reaction not in prior_dict is not in the ABC and should always be included.
        self.prior_set = np.ones(len(prior_dict))
        self.model.theta_rxn_map = {}
        rxn_theta_map = {}
        idx = -1
        for rxn_id, prior_value in prior_dict.iteritems():
            idx += 1
            self.prior_set[idx] = prior_value
            self.model.theta_rxn_map[idx] = self.model.reactions.get_by_id(rxn_id)
            rxn_theta_map[rxn_id] = idx
        self.prior_set_original = deepcopy(self.prior_set)       
        ## Use prior_multiplier to adjust entire prior enabling tuning of particle finders
        if prior_multiplier or (prior_multiplier == 0):
            print("Using prior multiplier ...")
            for idx, prior_val in enumerate(self.prior_set):
                prior_val_new = prior_val + (1.0-prior_val)*prior_multiplier
                self.prior_set[idx] = prior_val_new
        self.prior_estimate = deepcopy(self.prior_set)
        
        
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
        
        ## FIND ESSENTIAL REACTIONS HERE for precalc media
        ## For each reaction, find all RELs
        abc_reaction_lists = []
        for rxn_id in abc_reactions:
            rxn_list = []
            for reaction in self.model.reactions:
                if reaction.id.startswith(rxn_id):
                    rxn_list.append(reaction)
            abc_reaction_lists.append(rxn_list)
        
        if test_model_pre_abc:
            print("Testing full model with precalculated media ...")
            for medium in self.precalc_media:
                self.model.set_medium(medium)
                precalc_fail = False
                if self.model.opt() < 1e-5:
                    precalc_fail = True
                    print("Initial model did not run on precalc medium:")
                    for met, _ in medium.iteritems():
                        print(" - {}".format(met))
            if precalc_fail:
                print("Failed on precalc media, exiting ...")
                sys.exit(1)
            else:
                print("... done.")

        essential_abc_reaction_sets = []
        count_rxn_lists = loop_counter(
            len(abc_reaction_lists),
            'Testing ABC reactions for essentiality'
        )
        for rxn_list in abc_reaction_lists:
            essential_rxn_set = False
            for rxn in rxn_list:
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            for medium in self.precalc_media:
                self.model.set_medium(medium)
                if self.model.opt() < 1e-5:
                    essential_rxn_set = True
                    break
            if essential_rxn_set:
                essential_abc_reaction_sets.append(rxn_list)
            for rxn in rxn_list:
                rxn.lower_bound = rxn.lb_orig
                rxn.upper_bound = rxn.ub_orig
            count_rxn_lists.step()
        count_rxn_lists.stop()
        print(" => {} ABC reactions essential.".format(len(essential_abc_reaction_sets)))
        
        ## CONVERT ESSENTIAL REACTION SETS TO THETA SUBSETS - IF ONLY A SINGLE
        ## THETA, SET PRIOR TO ONE!
        essential_theta_sets = []
        for rxn_list in essential_abc_reaction_sets:
            if len(rxn_list) == 1:
                self.prior_set[rxn_theta_map[rxn.id]] = 1
                self.prior_estimate[rxn_theta_map[rxn.id]] = 1
            else:
                theta_indices = []
                for rxn in rxn_list:
                    theta_indices.append(rxn_theta_map[rxn.id])
                essential_theta_sets.append(set(theta_indices))
        
        self.essential_theta_sets = essential_theta_sets

        ## FIND BLOCKED REACTIONS AND SET BOUNDS TO 0
        print("Setting bounds of blocked reactions to 0 ...")
        blocked_rxn_ids = find_blocked_reactions(self.model, open_exchanges=True)
        num_restricted = 0
        for rxn_id in blocked_rxn_ids:
            try:
                blocked_reaction = self.model.reactions.get_by_id(rxn_id)
                blocked_reaction.upper_bound = 0
                blocked_reaction.lower_bound = 0
                num_restricted += 1
            except:
                print(" - 'Blocked' reaction '{}' could not be found".format(
                    rxn_id
                ))
        
        print(" => {} reactions had bounds restricted to 0.".format(num_restricted))
        
        
        num_essential_expts = None
        if not(not(self.test_type)) & (self.test_type == 'essential_and_orphan_mets'):        
            # get genes in model
            gene_set_all = [gene.id for gene in self.model.genes]
            num_essential_expts = 0
            for expt in self.experiments:
                if expt.result == 0:
                    gene_set = set(expt.genotype)
                    if gene_set.issubset(gene_set_all):
                        num_essential_expts += 1
            print("There are {} experiments showing essentiality".format(num_essential_expts))
        self.num_essential_expts = num_essential_expts
        
        
#         ## Test full model balanced accuracy
#         print("Testing full ABC model ...")
#         particle_full = self.initialise_particle(0)
#         particle_full.conduct_experiments()
#         print(" => balanced accuracy = {}".format(1.0-particle_full.result))
# #         print(" => Results:\n")
# #         initial_results.stats() 

        if test_model_pre_abc:
            ## Test full model
            print("\nTesting full ABC model ...")
            initial_results, _ = conduct_experiments(self.model, self.experiments)
            print(" => Results:\n")
            initial_results.stats()         
       
        print("ABC initialisation complete.")
                
    
    def set_rel_priors(self, rel_priors, prior_dict, include_all=False):
        """Where there is belief about certain enzyme/reaction pairs, 
        apply this to the relevant reaction priors.
        """
        if rel_priors:
            for rel_prior in rel_priors:
                base_rxn_id = rel_prior[0]
                gene_set = set(rel_prior[1])
                prior_value = rel_prior[2]
                id_string = base_rxn_id + "(_enz[0-9]+)*$"
                for rxn in self.model.reactions:
                    if re.match(id_string, rxn.id):
                        enzrxn_gene_set = set(gene.id for gene in rxn.genes)
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
            self.precalc_media_frozensets,
            self.prior_estimate,
            self.essential_theta_sets,
            self.test_type,
            self.num_essential_expts,
            self.prior_set_original)
        
        return particle
    
    def create_population_t(self):
        self.population_t = []
        for particle_id in range(self.N):
            particle_id_string = self.id_format.format(particle_id)
            self.population_t.append(self.initialise_particle(particle_id_string))
    
        
    def record_previous_results(self, theta_accepted_set, ln_w_accepted_set, distance_set, proposed_theta_idx_set):
        """Add the results from the previous run to the full lists of results."""
        
        self.results_theta.append(theta_accepted_set)
        self.results_d.append(distance_set)
        self.results_w.append(self.w_set_prev)
        self.thetas_selected.append(proposed_theta_idx_set)
        self.intermediate_theta_dict[self.t] = theta_accepted_set
        self.theta_set_prev = deepcopy(theta_accepted_set)   
        self.ln_w_accepted_set = ln_w_accepted_set
        self.create_prior_estimate()     
     
    def create_prior_estimate(self):
        """Use theta_accepted_set to estimate new prior for determining transition probabilities."""
        
        self.prior_estimate = []
        rxn_presence = zip(*self.results_theta[-1])
        for rxn in rxn_presence:
            rxn_prior = sum(rxn)/float(len(rxn))
            self.prior_estimate.append(rxn_prior)
        
    def step_forwards(self):
        """Increment time and calculate new problem parameters."""
        
        self.update_weights()
        if self.t < self.T-1:            
            self.t += 1
        else:
            print("This was the final run.")
            return None
        self.get_p_transition()
        self.get_epsilon()               
           
        
    def update_weights(self, ln_weights = None):
        """Normalise outputed weights and update w_set_prev."""
        if ln_weights is None:
            ln_weights = self.ln_w_accepted_set
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
    
    def test_theta(self, thetaTest):
        """Apply a manually defined theta and test against the experiments in 
        the ABC problem, returning a ResultSet instance with the results in it."""
        
        particleTest = self.initialise_particle(0)
        particleTest.theta_proposed = thetaTest    
        particleTest.epsilon = 1
        particleTest.apply_proposed_theta()
        particleTest.conduct_experiments()
        return particleTest.full_results
    
    def export_particle_rels(self, timepoint, particle_index, rel_output_file):
        """Export RELs from any particle from any of the timepoints."""
        f_out = open(rel_output_file,"w")
        ### For reaction in model,
        ###    if theta[t, index] == 1: export RELs 
        ###    elif not in index: export RELs.
        f_out.close()
    
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
            precalc_media_frozensets,
            prior_estimate,
            essential_theta_sets,
            test_type,
            num_essential_expts,
            prior_set_original):
         
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
        self.prior_estimate = prior_estimate
        self.essential_theta_sets = essential_theta_sets
        self.test_type = test_type
        self.num_essential_expts = num_essential_expts
        self.prior_set_original = prior_set_original
        
        self.num_params = len(prior_set)
         
        self.theta_sampled = None
        self.theta_sampled_idx = None
        self.theta_accepted = None
        self.theta_proposed = None
        self.weight = None
        self.full_results = None
        self.result = None
        
        self.verbose = False
        
        return None
    
    def vprint(self, string, verbose=None):
        """Print if verbose is true - override class verbosity with specified verbosity
        """
        if verbose is None:
            verbose = self.verbose
        if verbose:
            print(string)
        return None
    
    def reset_prior(self, prior_multiplier=None):
        """Set prior to original values, and apply prior_multiplier if specified 
        """
        self.prior_set = deepcopy(self.prior_set_original)
        if prior_multiplier or (prior_multiplier == 0):
            for idx, prior_val in enumerate(self.prior_set):
                prior_val_new = prior_val + (1.0-prior_val)*prior_multiplier
                self.prior_set[idx] = prior_val_new
        self.prior_estimate = deepcopy(self.prior_set)
        self.theta_accepted = None
        return None
         
    def propose_theta(self):
        """Sample theta from previous iteration and perturb according to K_t."""
         
        if self.t > 0: 
            ## Get theta by sampling
             
            self.theta_sampled_idx = np.random.choice(self.N,p=self.w_set_prev)
            self.theta_sampled = self.theta_set_prev[self.theta_sampled_idx]
             
            ## Perturb theta using perturbation kernel
            self.theta_proposed = self.K_t()
        else:
            ## Sample theta from prior
            self.theta_proposed = self.pi()
            self.theta_sampled_idx = -1
         
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
    
    def calculate_ln_weight_denominator_modified(self):
        """Calculate the logarithm of the weighted sum of the K_ts for the 
        calculation of weight.
        
        For the modified K_t, the shortcut leveraging the fact that all transition 
        probabilities are the same does not work so all must be calculated sequentially."""
        
        p = self.p_transition
        
        weighted_sum = 0
        for w_j, theta_j in zip(self.w_set_prev, self.theta_set_prev):
            weighted_residual = w_j * self.K_t_modified(theta_j)
            weighted_sum += weighted_residual
        
        try:
            pass
#            self.ln_weight_denominator = n_diff_res*(ln(p) + ln(1-p)) + ln(weighted_sum)
        except:
            print("This particle is so distant from every other particle that \
                    it has underflown, it will take the maximum value of the \
                    weights of the other particles in the population")
            self.ln_weight_denominator = None
    
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
            print("This particle is so distant from every other particle that \
                    it has underflown, it will take the maximum value of the \
                    weights of the other particles in the population")
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
              
    def find_accepted_theta(self, debug = False, use_ln = True, use_simple_K=False, max_tests=100000, verbose=None):
        """Find accepted theta and return it with weight."""
         
        count_attempts = 0
        failed_precalc = 0
        failed_theta_sets = 0
        print("\nFinding particle:".format(self.t, self.id))
        
        while True:
            if max_tests:
                if count_attempts >= max_tests:
                    print("\nToo many failed, exiting ...")
                    print("Attempts\tPrecalc\tThetasets")
                    print("{}\t{}\t{}".format(count_attempts,failed_precalc,failed_theta_sets))
                    return -1,-1,-1,-1
            count_attempts += 1            

            self.propose_theta()            
            excluded_theta_set = set([
                idx for idx, t_value
                in enumerate(self.theta_proposed)
                if t_value == 0])
            running_model=True
            for essential_theta_set in self.essential_theta_sets:
                if essential_theta_set.issubset(excluded_theta_set):
                    running_model=False
           
            if running_model:
                self.apply_proposed_theta()
                if not(not(self.test_type)) & (self.test_type == 'essential_and_orphan_mets'):
                    self.conduct_experiments_essential_and_orphan_mets(debug=debug, verbose=False)
                else:
                    self.conduct_experiments(debug=debug)
                if debug:
                    sys.stdout.write("\rresult = {:.3f}                                    ".format(self.result))
                    if self.result != 2:
                        sys.stdout.write(" after {}/{} tests".format(self.num_tests_checked, self.num_tests_total))
                    sys.stdout.flush() 
            else:
                self.result = 3
                
            ## If model is close enough to experiments, accept and return theta and calculated weight
            if self.result == 2:
                failed_precalc += 1
            elif self.result == 3:
                failed_theta_sets += 1
            elif self.result < self.epsilon:
                self.theta_accepted = self.theta_proposed
                self.calculate_ln_w()
                print("\nAccepted particle (result = {})\n".format(self.result))
                self.vprint(" => {} total attempts".format(count_attempts),verbose)
                self.vprint(" => {} failed on precalc".format(failed_precalc),verbose)
                self.vprint(" => {} failed on theta sets".format(failed_theta_sets),verbose)
                #sys.stdout.write(".")
                sys.stdout.flush()
                return self.theta_accepted, self.ln_w, self.result, self.theta_sampled_idx
            else:
                print("result = {}".format(self.result))
            sys.stdout.write("\r{}\t{}\t{}\t{}                            "
                .format(
                    count_attempts,
                    failed_precalc,
                    failed_theta_sets,
                    count_attempts-failed_precalc-failed_theta_sets))
            sys.stdout.flush()

        return None
        
    def perturb_param_using_prior(self, idx):
        """Return p_transition based on current prior estimate for REL idx."""
        
        rel_prior = self.prior_estimate[idx]
        if rel_prior > 0.95:
            rel_prior = 0.95
        if rel_prior < 0.05:
            rel_prior = 0.05
        
        param_perturbed = np.random.binomial(1,rel_prior)
        return param_perturbed

    def param_p_sourced(self,idx,param,param_previous):
        """Return the probability of a parameter having come from a previous parameter given a pertubation."""

        r_p = self.prior_estimate[idx]
        p_t = self.p_transition
        
        if r_p > 0.95:
            r_p = 0.95
        if r_p < 0.05:
            r_p = 0.05        
 
        if (param == 1) & (param_previous == 1):
            return p_t*r_p + 1-p_t
        if (param == 1) & (param_previous == 0):
            return p_t*r_p
        if (param == 0) & (param_previous == 1):
            return p_t*(1-r_p)
        if (param == 0) & (param_previous == 0):
            return p_t*(1-r_p) + 1-p_t
            
    def K_t_modified(self, theta_previous = None):
        """Perturbation function based on current estimated priors for reaction presence parameters. 
         
        Return a perturbed theta, theta_perturbed, from theta_sampled according to a modified distribution."""
        
        if theta_previous is None:
            ## Perturb the parameter set self.theta_sampled
            theta_perturbed = np.zeros(len(self.theta_sampled))
            for idx, param_sampled in enumerate(self.theta_sampled):
                if random() < self.p_transition:
                    param_perturbed = self.perturb_param_using_prior(idx)
                else:
                    param_perturbed = param_sampled
                theta_perturbed[idx] = param_perturbed 
            return theta_perturbed
        else:
            ## Output the probability of a transition to self.theta_sampled from theta_previous    
            #added_thetas = self.theta_accepted + theta_previous
            probability = 1.0
            ln_p = 0
            for idx, param in enumerate(self.theta_accepted):
                p_transition = self.param_p_sourced(idx, param,theta_previous[idx])
                ln_p += ln(p_transition)
                probability *= p_transition
            return probability, ln_p
   
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
            prior_iter = np.nditer(self.prior_set)
            for idx, param_prior in enumerate(prior_iter):
                param_value = np.random.binomial(1,param_prior)
                theta_proposed[idx] = param_value
            return theta_proposed
 
    def conduct_experiments(self, debug=False, epsilon=None, verbose=False):
        """Conduct experiments for model in current state and return 1 - (balanced accuracy)."""
         
#         ## Check that the model can run without changes, else return 1
#         if self.model.opt() <= 0:
#             print("Failed on original media")
#             self.result = 1
#             return None
        
        opt_cutoff = 1e-5
        
        self.num_tests_checked = 0
        self.num_tests_total = 0
        
        if debug:
            sys.stdout.write("\rChecking precalculated media ...                                     ")    
            sys.stdout.flush()
        ## Must run on all common media before experimental testing
        for precalc_medium in self.precalc_media:
            self.model.set_medium(precalc_medium)
            if self.model.opt() <= opt_cutoff:
                self.result = 2
                if debug:
                    print("Failed on precalc media")
                return None
         
#         print"No failure, continuing with test"
        if debug:
            sys.stdout.write("\rCreating list of valid experiments ...                           ")    
            sys.stdout.flush()
        ## Check all genotypes for presence in model and create a list of valid experiments
        valid_experiments = []
        num_pos_remaining = 0
        num_neg_remaining = 0
        gene_ids_in_model = self.model.get_relevant_gene_ids()
        for idx, expt in enumerate(self.experiments):
            all_genes_present = True
            for gene in expt.genotype:
                if gene not in gene_ids_in_model:
#                     print("Expt {}: '{}' not present in model ...".format(idx+1, gene))
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
        
        if verbose:
            print("Beginning tests ...")
        
        counter = loop_counter(len(valid_experiments), "Testing experiments")
        
        for idx, experiment in enumerate(valid_experiments):
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
                    break
            counter.step()        
        counter.stop()
        self.num_tests_checked = num_succeeded_tests + num_failed_tests

        self.full_results = deepcopy(running_results)
        distance = 1.0 - running_results.balanced_accuracy() 
        self.result = distance
        return None
    
    def conduct_experiments_essential_and_orphan_mets(self, debug=False, epsilon=None, verbose=None):
        """Conduct experiments for model in current state and return 1 - (balanced accuracy)."""
         
        opt_cutoff = 1e-5
        self.num_tests_checked = 0
        self.num_tests_total = 0
        self.vprint("\rChecking precalculated media ...                                     ",verbose)    
        ## Must run on all common media before experimental testing
        for precalc_medium in self.precalc_media:
            self.model.set_medium(precalc_medium)
            self.vprint(" => precalc: {}".format(self.model.opt()),verbose)
            if self.model.opt() <= opt_cutoff:
                self.result = 2
                self.vprint("Failed on precalc media",verbose)
                return None 

        ## Check all genotypes for presence in model and create a list of valid experiments
        valid_experiments = []
        num_pos_remaining = 0
        num_neg_remaining = 0        
        gene_ids_in_model = self.model.get_relevant_gene_ids()
        num_expts_essential = 0
        for expt in self.experiments:
            if expt.result == 0:
                num_expts_essential += 1
                all_genes_present = True
                for gene in expt.genotype:
                    if gene not in gene_ids_in_model:
                        all_genes_present = False
                        break
                if all_genes_present:
                    valid_experiments.append(expt)
                    if expt.result == 1:
                        num_pos_remaining += 1
                    else:
                        num_neg_remaining += 1
         
        self.num_tests_total = len(valid_experiments)
        self.vprint("\n{} valid experiments ...".format(len(valid_experiments)),verbose)
        self.vprint("Beginning tests ...",verbose)
        num_tn = 0
        num_fp = 0
        num_tests_remaining = len(valid_experiments)
        self.num_tests_checked = 0
        for experiment in valid_experiments:
            _, _, tn_add, fp_add, _ = experiment.test(self.model, self.precalc_media_frozensets)
            self.num_tests_checked += 1
            num_tn += tn_add
            num_fp += fp_add
            num_tests_remaining -= 1
                
            max_distance = 1.0 - (1.0 * num_tn) / self.num_essential_expts
            min_distance = 1.0 - (1.0 * num_tn + num_tests_remaining)/self.num_essential_expts
            
            self.vprint("{}\t{}\t{}".format(num_tests_remaining,min_distance, max_distance),verbose)
            
            if min_distance > self.epsilon:
                self.vprint("Minimum distance > epsilon, aborting ...",verbose)
                self.result = min_distance
                return None
            if max_distance < self.epsilon:
                self.vprint("Maximum distance < epsilon, finishing ...",verbose)
                self.result = min_distance
                return None
            
        distance = 1.0 - num_tn / self.num_essential_expts
        self.result = distance
        return None
        
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
        
        opt_cutoff = 1e-5
        
        precalc_frozensets = precalc_frozensets or []
        
        ec_model.set_medium(self.medium)
        if self.objective:
            ec_model.change_objective(self.objective)
        
        change_to_model = ec_model.set_genotype(self.genotype)
        if (change_to_model == 0) & (self.medium_components in precalc_frozensets):
            return 1, 1, 0, 0, 0 
        elif change_to_model == -1:
            ## Gene is not in model.
            return -2,0,0,0,0
        model_growth = timeout(ec_model.opt, default=0)
        ec_model.unset_genotype()
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        if (model_growth > opt_cutoff) and (self.result > 0):
            tp = 1
            expt_result = 1
        elif (model_growth <= opt_cutoff) and (self.result == 0):
            expt_result = 1
            tn = 1
        else:
            expt_result = 0
            if self.result > 0:
                fn = 1
            else:
                fp = 1
        return expt_result, tp, tn, fp, fn
        
def create_extended_model(model_file, objective_id = 'Biomass_BT_v2', require_solver=True, solver='cglpk',tidy_ids=True):
    """Take an ArrayBasedModel and convert to an Extended_Cobra_Model."""
    
    print("Creating extended COBRA model")
    model = create_cobra_model_from_sbml_file(model_file)
    ecm_model = ExtendedCobraModel(model)
    ecm_model.set_solver(solver, require_solver=require_solver)
    ecm_model.set_medium()
    print("done.")
    try:
        ecm_model.change_objective(objective_id)
    except:
        print("Objective could not be set ...")
    if tidy_ids:
        for rxn in ecm_model.reactions:
            newID = re.sub('_LPAREN_','_',deepcopy(rxn.id))
            newID = re.sub('_RPAREN_','',newID)
            newID = re.sub('\(','_',newID)
            newID = re.sub('\)','',newID)
            newID = re.sub('_LSQBKT_','_',newID)
            newID = re.sub('_RSQBKT_','',newID)
            if newID.endswith('_er'):
                newID = newID[:-1]           
            rxn.id = newID       
        ecm_model.repair() 
        for gene in ecm_model.genes:
            gene_id = re.sub('^(G_)+','',gene.id)
            gene.id = gene_id
        ecm_model.repair()
        for reaction in ecm_model.reactions:
            grr = re.sub('(G_)','',reaction.gene_reaction_rule)
            reaction.gene_reaction_rule = grr        
        ecm_model.repair()
    return ecm_model    

class ExtendedCobraModel(ArrayBasedModel):
    """Additional functionality for helping to use COBRA model."""
            
    def get_relevant_gene_ids(self):
        """Get a list of gene IDs for all genes implicated in non-zero flux reactions."""
        
        relevant_genes_list = []
        for reaction in self.reactions:
            if (reaction.lower_bound != reaction.upper_bound) or (reaction.lower_bound != 0):
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
        if len(gpr_list) > enzyme_limit:
            return [rxn.id], None
        elif len(gpr_list) == 0:
            return [], rxn.id
        enzyme_index = 0
        enzrxn_id_list = []
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
        
    def set_genotype(self, genotype, return_ko_rxns=False):
        """
        Run delete_model_genes on model.  Return False if no change is made to the model
        """
        if len(genotype) == 0:
            return -2
        try:
            ## CODE FROM COBRAPY: to access reaction deletions
            # Allow a single gene to be fed in as a string instead of a list.
            if not hasattr(genotype, '__iter__') or hasattr(genotype, 'id'):
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
        if return_ko_rxns:
            return knocked_out_reactions
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
        
    
    def set_solver(self, solver_string = 'cglpk', require_solver=True):
        try:
            self.optimize(solver=solver_string)
            self.solver = solver_string
        except:
            try:
                self.optimize(solver='glpk')
                self.solver = 'glpk'
            except:
                if not require_solver:
                    self.solver = None
                else:
                    self.optimize(solver='glpk')
                    
        
        print("Solver is '{}'.".format(self.solver))
        
        
    def opt(self, new_objective = None, time_limit=1):
        if new_objective:
            self.change_objective(new_objective)
        
#         sys.stdout.write("\rrunning optimize() ...                         ")    
#         sys.stdout.flush()
        try:
            self.optimize(solver=self.solver, time_limit=time_limit)
        except:
            
#             print("Optimize failed.")
#             print("Solver: '{}'".format(self.solver))
#             print("Time limit: {}".format(time_limit))
            return 0

#         sys.stdout.write("\rfinished running optimize() ...                ")    
#         sys.stdout.flush()
        
        if self.solution.f:
            if self.solution.f < 0:
#                 print("Solution less than 0.")
                return 0
            else:
                return self.solution.f
        else:
#             print("No solution.")
            return 0
    
    def set_medium(self,medium_dict=None,debug=False):
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
                else:
                    if debug:
                        print("Reaction ID '{}' not found, ignoring ...".format(component))
        self.medium_dict = medium_dict
               
    def show_medium(self):
        """
        Find all exchange lower bounds < 0 and print to screen.
        """
        for reaction in self.reactions:
            if reaction.id.startswith("EX_"):
                ## Exchange reaction, so if lb < 0, print
                if reaction.lower_bound < 0:
                    reaction_identifier = reaction.name + " (" + reaction.id + ")"
                    print("%45s %5.0f %5.0f" % (reaction_identifier, reaction.lower_bound, reaction.upper_bound))
    
def conduct_experiments(model, experiments, debug = False, epsilon=None, verbose=False):
    """Conduct experiments for model in current state and return 1 - (balanced accuracy)."""
          
    if debug:
        sys.stdout.write("\rCreating list of valid experiments ...                           ")    
        sys.stdout.flush()
    ## Check all genotypes for presence in model and create a list of valid experiments
    valid_experiments = []
    num_pos_remaining = 0
    num_neg_remaining = 0
    gene_ids_in_model = model.get_relevant_gene_ids()
    for expt in experiments:
        all_genes_present = True
        for gene in expt.genotype:
            if gene not in gene_ids_in_model:
                all_genes_present = False
                break
        if all_genes_present:
            valid_experiments.append(expt)
            if expt.result == 1:
                num_pos_remaining += 1
            else:
                num_neg_remaining += 1
     
    print("\n{} valid experiments ...".format(len(valid_experiments)))
    
    num_failed_tests = 0
    num_succeeded_tests = 0 
    running_results = ResultSet(0,0,0,0)
        
    if verbose:
        print("Beginning tests ...")
    
    counter = loop_counter(len(valid_experiments), "Testing experiments")
    
    for experiment in valid_experiments:
        expt_result, tp_add, tn_add, fp_add, fn_add = experiment.test(model)
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
        counter.step()
    
    counter.stop()
    num_tests_checked = num_succeeded_tests + num_failed_tests
    full_results = deepcopy(running_results) 
    return full_results, num_tests_checked        

def timeout(func, args=(), kwargs={}, timeout_duration=60, default=None):
    """Run func, but exit after timeout_duration if func hasn't completed.
    """
    import signal
    class TimeoutError(Exception):
        pass
    def handler(signum, frame):
        raise TimeoutError()
    # set the timeout handler
    signal.signal(signal.SIGALRM, handler) 
    signal.alarm(timeout_duration)
    try:
        result = func(*args, **kwargs)
    except TimeoutError as exc:
        print("\ntimed out ...\n")
        result = default
    finally:
        signal.alarm(0)
    return result

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
    
def import_expt_data(model, objective_id="Biomass_BT_v2", media=None, data_file=None, debug=False):
    """From a table of genotype/growth information, create the experiments 
    for a particular organism.
    
    If there are any lines at the top of the file starting with the word 
    "MEDIUM" in the first column the next column will be the medium name and 
    the following cells will contain medium components for specifying media 
    used in the experiments. 
    """
    print("Importing data ...")
    objective = model.reactions.get_by_id(objective_id)
    data_file = data_file or '/Users/wbryant/work/BTH/analysis/gene_essentiality/essentiality_data_complete.csv'
    experiments = []
    expts_in = open(data_file,'r')
    if not media:
        media = {}
    for line in expts_in:
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
                        c_source_exchange = "EX_" + source + "_e"
                        expt_medium[c_source_exchange] = -100
                experiment = Experiment(details[0], expt_medium, result, genotype, objective)
                experiments.append(experiment)
    expts_in.close()
    return experiments

from time import time, sleep
from math import floor

class loop_counter:
    """Use to track progress of a loop of known length."""
    
    def __init__(self, length, message = 'Entering loop', timed=False, add_delay=False):
        self.stopped = False
        self.length = length
        self.num_done = 0
        self.next_progress = 1
        self.percent_done = 0
        self.timed = timed
        self.add_delay = add_delay
        if self.timed:
            self.time_0 = time()
        print("{}:".format(message))
        try:
            if self.add_delay:
                sleep(0.001)
            sys.stdout.write("\r - 0 %%")
            sys.stdout.flush()
        except:
            pass
    
    def step(self):
        self.num_done += 1
        if not self.stopped:
            if self.num_done >= self.length:
                self.stop()
            else:
                percent_done = floor(1000.0*self.num_done/self.length)/10
                try:
                    if self.add_delay:
                        sleep(0.001)
                    sys.stdout.write("\r - {} % ({} / {})".format(
                        percent_done,
                        self.num_done,
                        self.length
                    ))
                    sys.stdout.flush()
                except:
                    pass
                
    def stop(self):
        if not self.stopped:
            try:
                if self.add_delay:
                    sleep(0.001)
                sys.stdout.write("\r - 100 % ({})                        \n"
                    .format(self.length))
            except:
                pass
            self.stopped = True
            sys.stdout.flush()

class ResultSet:
    """A container for a set of results that can return stats"""
    
    def __init__(self, tp=0, tn=0, fp=0, fn=0):
        
        self.tp = tp
        self.tn = tn
        self.fp = fp
        self.fn = fn
    
    def precision(self):
        try:
            return self.tp / float(self.tp+self.fp)
        except:
            return 0
    
    def recall(self):
        try:
            return self.tp / float(self.tp+self.fn)
        except:
            return 0
            
    def f_measure(self):
        precision = self.precision()
        recall = self.recall()
        try:
            return 2*precision*recall/(precision + recall)
        except:
            return 0
    
    def accuracy(self):
        try:
            return float(self.tp + self.tn) / (self.tp + self.tn + self.fp + self.fn)
        except:
            return 0
    
    def balanced_accuracy(self):
        try:
            return 0.5 * (float(self.tp)/(self.tp + self.fn) + float(self.tn)/(self.tn+self.fp))
        except:
            return 0
    
    def max_balanced_accuracy(self, num_tests_remaining):
        """Given a number of tests to be added to the results, what is the 
        maximum balanced accuracy that could be achieved?  
        """
        
        pass
    
    def calc_balanced_tfpn_addition(self, num_tests_remaining):
        denominator = self.fn + self.fp
        numerator = self.fn*self.tn - self.tp*self.fp + self.fn*num_tests_remaining
        
        x = round(numerator/float(denominator))
        y = num_tests_remaining - x
        
        return x, y
         
        
        
    def stats(self):
        if self.tp > 0:
            precision = self.precision()
            recall = self.recall()
            f_measure = self.f_measure()
            accuracy = self.accuracy()
            balanced_accuracy = self.balanced_accuracy()
            
            print("tp\ttn\tfp\tfn")
            print("{}\t{}\t{}\t{}\n".format(
                self.tp,
                self.tn,
                self.fp,
                self.fn
            ))
            
            print("%12s = %1.3f" % ('Precision', precision))
            print("%12s = %1.3f" % ('Recall', recall))
            print("%12s = %1.3f" % ('F-measure', f_measure))
            print("%12s = %1.3f" % ('Accuracy', accuracy))
            print("%12s = %1.3f" % ('Bal. Accuracy', balanced_accuracy))
        else:
            print("There are no true positives.")
            print("tp\ttn\tfp\tfn")
            print("{}\t{}\t{}\t{}\n".format(
                self.tp,
                self.tn,
                self.fp,
                self.fn
            ))

if __name__ == '__main__':
    pass