'''
Created on 10 Dec 2014

@author: wbryant
'''

"""
The contents of this file should enable ABC-SMC to be run for a particular 
set of model/experiment data.
"""

import paramiko
# import subprocess32 as subprocess
import sys
import pickle
from time import sleep as wait
from time import time
from local_utils import loop_counter
from copy import deepcopy
from os import remove as remove_file
import os

base_dir = "/bmm/home/wbryant/BTH/analysis/ABC-SMC/single_particle_test/"

def run_remote_abc(abc_problem, particles_per_job, queue_length):
    """Use the PopulationSubmitter to run an entire ABC-SMC on the farm."""
    
    for timepoint in range(abc_problem.T):
        pop_sub = PopulationSubmitter(abc_problem, particles_per_job, queue_length)
        pop_sub.submit_all(abc_problem)
    
    abc_problem.get_posterior()
            
def get_job_hours_elapsed():
    """Use qstat to determine how many hours the current jobs have been running."""
    
    ## This command gets elapsed time from the qstat command and selects the 
    ## hours column, returning the maximum value
    command = "qstat -u wbryant | tr -s ' ' | rev | cut -d' ' -f 1 | rev | cut -c2"
    
    stdout, _ = run_q_on_bamboon(command, ret_stdout = True)
    entries = stdout.readlines()
    max_hours_elapsed = 0
    for entry in entries:
        try:
            hours_elapsed = int(entry.strip())
            if hours_elapsed > max_hours_elapsed:
                max_hours_elapsed = hours_elapsed
        except:
            continue
    return max_hours_elapsed

def run_q_on_bamboon(command, ret_stdout = False):
    """Execute a command on Bamboon and return the stripped stdout first line."""
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        ssh.connect('bamboon.bc.ic.ac.uk', username='wbryant', password='h73B,xh9R')
    except:
        wait(10)
        ssh.connect('bamboon.bc.ic.ac.uk', username='wbryant', password='h73B,xh9R')
#     print("\n\nSubmitting '{}' to Bamboon.".format(command))
    command = command.strip() 
    stdin, stdout, stderr = ssh.exec_command(command)
    if ret_stdout:
        return stdout, stderr
    text = stdout.readlines()
    stderr_text = stderr.readlines()
#     stdin_text = stdin.readlines()
    try:
        return float(text[0].split(".")[0].strip())
    except:
#         print("\nCOMMAND:")
#         print("'" + command + "'")
#         print("\nSTDOUT:")
#         print("'" + "\n".join(text) + "'")
#         print("STDERR:")
#         print("'" + "\n".join(stderr_text) + "'")
#         sys.stdout.flush()
        return None

def bins(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

class PopulationSubmitter():
    """Takes a list of particles and bins them into a number of jobs to be
    submitted to the farm.  It then has submission and monitoring capabilities! 
    """
    
    def __init__(self, abc_problem, particles_per_job, queue_length, wall_time = 10):
        
        self.particles_per_job = particles_per_job
        self.queue_length = queue_length
        self.abc_problem = abc_problem
        self.wall_time = wall_time
        
        print("Using wall time = {}".format(self.wall_time))
        
        sys.stdout.write("Clearing running directories ...")
        run_dirs = ['logs', 'output', 'stdout', 'errors', 'jobs']
        rm_command_template = "rm -f {}{}/* "
        rm_command_chain = "; ".join([rm_command_template.format(base_dir, run_dir) for run_dir in run_dirs])
        #print rm_command_chain.strip()
        run_q_on_bamboon(rm_command_chain.strip())
        sys.stdout.write(" done.\n")
        
    def update_submitter(self):
        """Update submitter when ABC problem moves to next population."""
        
#         print("Moving submitter to time = {} ...".format(self.abc_problem.t))
        print("Moving submitter to time = {} (epsilon = {}) ...".format(self.abc_problem.t, self.abc_problem.epsilon_t))
        
        self.queue = []
        self.job_submit_id = -1
        
        ## Create particles with job data
        print("Creating particles ...")
        self.particle_ds = []
        self.abc_problem.create_population_t()
#         counter = loop_counter(self.particles_per_job, "Creating particles ...")
        for particle in self.abc_problem.population_t:
#             counter.step()
            particle_d = ParticleData(particle)
            self.particle_ds.append(particle_d)
#         counter.stop()
        
        ## Create jobs
        print("Creating jobs ...")
        self.jobs = []
        counter = loop_counter(len(list(bins(self.particle_ds, self.particles_per_job))), "Creating jobs")
        for idx, particle_d_set in enumerate(bins(self.particle_ds, self.particles_per_job)):
            self.jobs.append(Job(particle_d_set, self.abc_problem.t, idx))
            counter.step()
        counter.stop()
        print(" - done.")
        
    def submit_all_direct(self):
        """Submit all jobs directly.  Do not submit at all if the number of 
        jobs is greater than the maximum queue size."""
        
        num_jobs = len(self.jobs)
        if num_jobs > self.queue_length:
            print("Too many jobs - please adjust parameters!")
            sys.exit(1)
        else:
            time_0 = time()
            print("Submitting all jobs ...")
            submitted_ids = []
            for job in self.jobs:
                if os.path.isfile(pickle_file):
                    print("Submitting job {}, file present.".format(job.name))
                job.submit()
                submitted_ids.append(job.job_id)
                wait(20)
            
            print("Monitoring jobs ...")
            ticker = 0
            interval = 5
            while True:
                wait(30)
                num_jobs_running = int(run_q_on_bamboon("qstat -u wbryant | wc -l")) - 5
                ticker += 1
                if num_jobs_running <= 0:
                    print("\nAll jobs finished after {:.1f} hours.".format((time()-time_0)/3600))
                    sys.stdout.flush()
                    break                    
                hours_elapsed = get_job_hours_elapsed()
                if ticker/interval == ticker/float(interval):
                    sys.stdout.write("\r{} jobs running after {:.1f} minutes"
                        .format(int(num_jobs_running), (time()-time_0)/60))
                    print("\n - {} hours elapsed, {} hours wall time".format(hours_elapsed, self.wall_time))
                    sys.stdout.flush()
                if hours_elapsed >= self.wall_time:
                    print("Reached wall time, cancelling jobs ...")
                    qdel_commands = ["qdel {}".format(int(qid)) for qid in submitted_ids]
                    for command in qdel_commands:
                        wait(10)
                        print("Running: '{}'".format(command))
                        sys.stdout.flush()
                        run_q_on_bamboon(command)
                    return 1


                    
        ## When all jobs are complete, update self.abc_problem
        theta_accepted_set = []
        ln_w_accepted_set = []
        distance_set = []
        print("Accumulating data from accepted particles ...")
        sys.stdout.flush()
        for particle_d in self.particle_ds:
            f_in = open(particle_d.output_file, 'r')
            data = f_in.readline().strip().split("\t")
            f_in.close()
            try:
                theta_accepted = [int(theta) for theta in data[0].split(",")]
                ln_w_accepted = float(data[1])
                theta_accepted_set.append(theta_accepted)
                ln_w_accepted_set.append(ln_w_accepted)
                distance_set.append(float(data[2]))
            except:
                print("Results could not be read in for {} with data '{}'"
                    .format(particle_d.name, data))
                print("File '{}' did not conform to correct specification"
                    .format(particle_d.output_file))
                return 1

        print("Updating ABC problem ...")
        sys.stdout.flush()
        try:
            self.abc_problem.step_forwards(theta_accepted_set, ln_w_accepted_set, distance_set)
        except:
            print("step_forwards failed to execute")
            print theta_accepted_set
            print ln_w_accepted_set
            print distance_set
            print self.abc_problem.t
            print(" -> Pickling abc_problem ...")
            abc_results_file = base_dir + 'abc_problem_final.pickle'
            f_out = open(abc_results_file, 'w')
            pickle.dump(self.abc_problem, f_out)
            f_out.close()
        return None


class Job():
    """A single job comprised of 1 or more particles for serial calculation."""
    
    def __init__(self, particle_data, timepoint, parallel_id):
        
        self.parallel_id = parallel_id
        self.name = "job_" + str(timepoint) + "." + str(parallel_id)
        self.job_id = None
        ## Job status: 0 = not submitted, 1 = submitted and running, 2 = finished, -1 = failed
        self.job_status = 0
        self.particle_list = []
        
        ## Make text for insertion into job command
        if isinstance(particle_data, basestring):
            self.particle_data = [particle_data]
        else:
            self.particle_data = particle_data
        
        ## Set the same input file, and create that file from the first particle_d
        self.pickle_file = base_dir + "particles/" + self.name + ".pickle"
        open(self.pickle_file, 'w').close()
        particle_data[0].make_pickle(self.pickle_file)
        self.run_command_text = ""
        for particle_d in particle_data:
            self.particle_list.append(particle_d.name)
            particle_text = particle_d.job_text.format(self.pickle_file)
            self.run_command_text += particle_text
        
        ## Create job file from template
        self.job_file = base_dir + "jobs/" + self.name + ".com"
        self.log_file = base_dir + "logs/" + self.name + ".log"
        open(self.log_file, 'w').close()
        self.job_text = job_file_template.format(
            base_dir,
            self.log_file,
            self.run_command_text)
        f_job = open(self.job_file, 'w')
        f_job.write(self.job_text)
        f_job.close()
        
        ## Create qsub command - requires ID, stdout file and error file
        self.stdout_file = base_dir + "stdout/" + self.name + ".o"
        open(self.stdout_file, 'w').close()
        self.error_file = base_dir + "errors/" + self.name + ".e"
        open(self.error_file, 'w').close() 
        self.qsub_command = qsub_template.format(self.stdout_file, self.error_file, self.job_file)       
        

    def submit(self):
        """Run qsub command and update with relevant."""
        qsub_status = run_q_on_bamboon(self.qsub_command)
        if qsub_status > 0:
            self.job_status = 1
            self.job_id = qsub_status
            return 0
        else:
            self.job_status = -1
            return -1
          

class ParticleData():
    """Creates all data needed for adding a particle to a job. 
    """
    
    def __init__(self, particle):
        
        self.particle = particle
        
        ## Job status: 0 = not submitted, 1 = submitted and running, 2 = finished, -1 = failed
#         self.job_status = 0
#         self.job_id = None
        self.name = "particle_{}.{}".format(self.particle.t, self.particle.id)
        
#         ## Pickle particle into particles directory
#         self.pickle_file = base_dir + "particles/" + self.name + ".pickle"
#         f_pickle = open(self.pickle_file, 'w')
#         pickle.dump(self.particle, f_pickle)
#         f_pickle.close()
        
        self.output_file = base_dir + "output/" + self.name + ".out"
        self.job_text = particle_command_template.format("{}", self.output_file)

    def make_pickle(self, pickle_file):
        """Output a pickle file of the particle for job input."""
        
        f_pickle = open(pickle_file, 'w')
        pickle.dump(self.particle, f_pickle)
        f_pickle.close()

qsub_template="""qsub -lnodes=limbo -q long -o {} -e {} {}"""

particle_command_template = """/bmm/home/wbryant/envs/cobra/bin/python run_particle.py -i {} -o {} >> $MY_TMP_FILE\n"""

job_file_template = """#!/bin/sh

# A job is started as
#     qsub -lnodes=limbo -t 1-4 job.com
#
# This submits job.com 4 times and set the variable
# $PBS_ARRAYID values, 1,2,3 and 4
#
# A long job is started as
#     qsub -lnodes=limbo -t 1-10 -q long job.com 

# --- Set job number and name (if array id is blank set to 1)

# JOB PARAMETERS

export MY_RUN_DIRECTORY="{}"
export MY_TMP_FILE="{}"

# JOB SUBMISSION

cd /bmm/home/wbryant/myutils/abcsmc
echo $HOSTNAME >> $MY_TMP_FILE 
echo "Date started:" `date` >> $MY_TMP_FILE
{}
echo "Date finished: " `date` >> $MY_TMP_FILE
"""

if __name__=="__main__":
    import myutils.abcsmc.abcsmc as abcsmc
    
    model_file = "/bmm/home/wbryant/BTH/analysis/ABC-SMC/full_run/datafiles/bth_model.xml"
    ess_file = "/bmm/home/wbryant/BTH/analysis/ABC-SMC/full_run/datafiles/essentiality.csv"
    media_file = "/bmm/home/wbryant/BTH/analysis/ABC-SMC/full_run/datafiles/media.data"
    priors_file = "/bmm/home/wbryant/BTH/analysis/ABC-SMC/full_run/datafiles/priors.txt"
    
    print("Importing model ...")
    sys.stdout.flush()
    bth_model = abcsmc.create_extended_model(model_file)
    print("Importing media ...")
    sys.stdout.flush()
    media = abcsmc.import_media(media_file)
    print("Importing priors ...")
    sys.stdout.flush()
    prior_dict = abcsmc.import_prior_dict(priors_file)
    print("Importing experimental data ...")
    sys.stdout.flush()
    expts = abcsmc.import_expt_data(bth_model, media, data_file = ess_file)
    
    
    print("Setting up ABC problem ...")

    
    abc_options = {}
    abc_options['default_prior_value'] = 0.99
    abc_options['epsilon_0'] = 0.45
    abc_options['epsilon_T'] = 0.4
    abc_options['alpha'] = 0.3
    abc_options['particles_per_population'] = 1050
    abc_options['num_populations_max'] = 5
    abc_options['model'] = bth_model
    abc_options['prior_dict'] = prior_dict
    abc_options['experiments'] = expts
    
    abc_problem = abcsmc.AbcProblem(**abc_options)
    
    print("Initialising ABC job submitter ...")
    sys.stdout.flush()
    submitter_options = {}
    submitter_options['wall_time'] = 10 # hours
    submitter_options['abc_problem'] = abc_problem
    submitter_options['particles_per_job'] = 35
    submitter_options['queue_length'] = 30
    pop_sub = PopulationSubmitter(**submitter_options)
    
    time_0 = time()
    time_t = time()
    
    
    
    for t in range(abc_problem.T):
        print("\nSubmitting population {} to Bamboon ...".format(t))
        pop_sub.update_submitter()    
        no_result = pop_sub.submit_all_direct()
        time_end =  time()
        print("... population {} took {} seconds.".format(t, time_end-time_t))
        time_t = time_end
        if no_result:
            abc_problem = pop_sub.abc_problem
            if abc_problem.t == 0:
                print("Not enough acceptable particles could be found in the first population, exiting ...")
                sys.exit(1)
            abc_problem.t -= 1
            print("Total time for ABC-SMC = {}.".format(time()-time_0))
            break

        
    
    ## Get data from final run and return results
    ## N.B. Final population is accessible in pop_sub.particle_ds
    
    model_out = abc_problem.model
    theta_rxn_map = model_out.theta_rxn_map
    results_dict = {}
    mean_dict = {}
    mean_w_dict = {}
    weight_dict = {}
    rxn_ids = [reaction.id for reaction in model_out.reactions]
    try:
        theta_accepted_set = abc_problem.intermediate_theta_dict[abc_problem.t]
    except:
        print("'t' was too high")
        theta_accepted_set = abc_problem.intermediate_theta_dict[abc_problem.t-1]
    weights = abc_problem.w_set_prev
    
    theta_weighted_accepted_set = []
    
    for idx, w in enumerate(weights):
        theta_weighted_accepted_set.append([w*param for param in theta_accepted_set[idx]])
    
    results_per_param = zip(*theta_accepted_set)
    results_w_per_param = zip(*theta_weighted_accepted_set)
    
#     print("Rxn_ids (model_out):")
#     print("\n".join(rxn_ids))
    
#     print("Rxn_ids (theta_rxn_map):")
    for idx, rxn in theta_rxn_map.iteritems():
        rxn_id = rxn.id
#         print "'" + rxn_id + "'", sum(results_per_param[idx])/float(len(results_per_param[idx])), sum(results_w_per_param[idx])
        results_dict[rxn_id] = results_per_param[idx]
        mean_dict[rxn_id] = sum(results_per_param[idx])/float(len(results_per_param[idx]))
        mean_w_dict[rxn_id] = sum(results_w_per_param[idx])
#         weight_dict[rxn_id] = weights[idx]
    
    print("Pickling abc_problem ...")
    abc_results_file = base_dir + 'abc_problem_final.pickle'
    f_out = open(abc_results_file, 'w')
    pickle.dump(abc_problem, f_out)
    f_out.close()
    
#     print("\nResults:\n")
#     for rxn_id in rxn_ids:
#         if rxn_id not in results_dict:
# #             print "'" + rxn_id + "'"
#             mean_w_dict[rxn_id] = 1
#             mean_dict[rxn_id] = 1
#         
#         print("{}:\t{}\t({})".format(rxn_id, mean_dict[rxn_id], mean_w_dict[rxn_id]))
       
    print("All finished.")
        