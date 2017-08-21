# ARBOC-SMC
An implementation of AppRoximate Bayesian COmputation by Sequential Monte Carlo for genome scale metabolic models (using COBRApy).

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
