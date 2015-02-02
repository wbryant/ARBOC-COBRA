import sys, argparse, os
import cPickle as pickle


def main(argv=None): # IGNORE:C0111
    '''Command line options.'''
#     print("Running run_particle!")
    sys.stdout.flush()
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)
        
    parser = argparse.ArgumentParser(description="Run single particle until result accepted.")
    parser.add_argument("-i", "--input", action="store", dest="pickle_file", required=True, help="Path and name of pickle file containing particle.")
    parser.add_argument("-o", "--output", action="store", dest="result_file", required=True, help="Path and name of file for output.")
    args = parser.parse_args()
    
    pickle_file = args.pickle_file
    result_file = args.result_file
    
    print("\nLoading particle (dev run_particle) ...")
    if os.path.isfile(pickle_file):
        print("Reading pickle file ...")
        f_in = open(pickle_file, 'r')
        particle = pickle.load(f_in)
        f_in.close()
    else:
        print("'{}' does not exist".format(pickle_file)) 
        sys.exit()
    
    print("Finding acceptable theta ...")
    theta_accepted, ln_w, distance = particle.find_accepted_theta()
    print("... accepted theta.")
    
    ## Write results to output file
    f_out = open(result_file, 'w')
    f_out.write("{}\t{}\t{}".format(",".join([str(int(theta)) for theta in theta_accepted]), str(ln_w), str(distance)))
    f_out.close()
    print("Printed results to file.")
    
    

if __name__ == "__main__":
    main()
    
    