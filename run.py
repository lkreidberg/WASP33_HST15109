import sys
sys.path.insert(0, './fit_funcs')
import numpy as np
import glob
import os
from astropy.io import ascii
import getopt
import time as pythontime
from read_data import Data
from model import Model
from least_squares import lsq_fit
                    
def usage():
    cmd = sys.argv[0]
    sys.stderr.write('Usage: python %s OPTION\n\n' % os.path.basename(cmd))
    sys.stderr.write(
        'Allowed OPTION flags are:\n'
        '  --show-plot       displays fitted light curve plots\n'
        '  --run-mcmc        runs MCMC starting from least-squares parameters\n'
        '  --plot-raw-data   plots raw light curve separated by visit\n'   
        '  --path PATH	     specifies PATH to light curves to be fit\n' 
        '  --fit-white FILE  fits the white light curve stored in FILE\n' 
        '  -v                prints fit diagnostic information\n'
        '  -o                saves fit output to file\n'
        '  -h                lists instructions for usage\n'
        '\n')
    sys.exit(1)


def main():
    myfuncs = ['constant']
    #model = ['constant', 'polynomial1', 'eclipse', 'model-ramp']
    #other options: polynomial2, transit, sine, divide-white

    #significance above which to mask outliers
    #outlier_cut = 10.

    #parses command line input
    try: opts, args = \
            getopt.getopt(sys.argv[1:], 
                "hov", ["help", "show-plot", "run-mcmc", "plot-raw-data", 
                "plot-sys", "path=", "fit-white=", "divide-white"]
            )
    except getopt.GetoptError: usage()

    #defaults for command line flags
    verbose         = False
    output          = False
    show_plot       = False
    run_mcmc        = False
    run_lsq         = True
    plot_raw_data   = False
    path            = "spec_lc"
    fit_white       = False
    divide_white    = False

    for o, a in opts:
        if o in ("-h", "--help"): usage()
        elif o == "-o": output = True
        elif o == "-v": verbose = True
        elif o == "--show-plot": show_plot = True
        elif o == "--run-mcmc": run_mcmc, run_lsq = True, False
        elif o == "--run-lsq": run_lsq = True
        elif o == "--plot-raw-data": plot_raw_data = True
        elif o == "--path": path = a
        elif o == "--fit-white": fit_white, white_file = True, a
        elif o == "--divide-white": divide_white = True
        else: assert False, "unhandled option"

    flags = {'verbose': verbose, 'show-plot': show_plot, 
            'plot-raw-data': plot_raw_data, 'output': output, 
            'out-name': 'none.txt', 'run-lsq': run_lsq, 
            'run-mcmc': run_mcmc, 'divide-white': divide_white, 
            'fit-white': fit_white}		

    #reads in observation and fit parameters
    obs_par = {x['parameter']: x['value'] for x in 
                ascii.read("config/obs_par.txt", Reader=ascii.CommentedHeader)
              }
    fit_par =   ascii.read("config/fit_par.txt", Reader=ascii.CommentedHeader)

    files = glob.glob(os.path.join(path, "*"))		
    if fit_white: files = glob.glob(white_file)

    flags['out-name'] = "fit_" + pythontime.strftime("%Y_%m_%d_%H:%M") + ".txt"

    for f in files:
        data = Data(f, obs_par, fit_par)
        if flags['divide-white']:
            sys_vector = np.genfromtxt("white_systematics.txt")
            data.all_sys = sys_vector
        model = Model(data, myfuncs)
        data, model = lsq_fit(f, obs_par, fit_par, data, flags, model, myfuncs)

        #outfile = open("white_systematics.txt", "w")
        #for i in range(len(m.all_sys)): print>>outfile, m.all_sys[i]
        #outfile.close()

        if flags['run-mcmc']:
            output = mcmc_fit(f, obs_par, fit_par, flags)
            for i in range(len(output)): 
                    temp = "mcmc_out_"+'{0:0.2f}'.format(output[i][0])
                    np.save(temp, output[i][3])
                            

if __name__ == '__main__':
    main()
