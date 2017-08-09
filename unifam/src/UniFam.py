#!/usr/bin/python

'''
UniFam.py

pipeline

Created by JJ Chai on 02/24/2014
Last modified Mon Apr 20 16:58:56 EDT 2015
Copyright (c) 2014 JJ Chai (ORNL). All rights reserved.

'''
# Import Python modules
import ConfigParser
import argparse
import sys
from datetime import datetime

# Import local modules
import UniFam_lib # in this directory
import UniFam_lib_batch

## Version
version_str = "1.1.0"
''' 0.0.1.  first version of pipeline, including prodigal, hmmsearch, and annotation
            all configuration options are in a file. User can(and must) change the options in the config file to customize.
    0.0.2.  Added more information for pathologic module of the analysis, and rRNA, tRNA            
            analysis.
    1.0.0   First Stable version for release, UniFam 1.0.0
    1.1.0   Added README file to describe output files; zip pathway inference output results
'''

parser = argparse.ArgumentParser(description="Annotation of contigs/proteins using UniFam",
                                 prog = 'UniFam', #program name
                                 prefix_chars='-', # prefix for options
                                 fromfile_prefix_chars='@', # if options are read from file, '@args.txt'
                                 conflict_handler='resolve', # for handling conflict options
                                 add_help=True, # include help in the options
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter # print default values for options in help message
                                 )
## version control
parser.add_argument("--version", action="version",version='%(prog)s {}'.format(version_str))

## --verbose mode, default: false (quiet mode)
parser.add_argument("-v", "--verbose", action="store_true",help="verbose mode, more output")
#parser.add_argument("-n", "--dryrun", action="store_true",help="dryrun, only print commands, do not execute")

parser.add_argument("-b", "--batch", action="store_true",help="batch mode to construct pathway", dest="batch")

## input files and directories
## configuration file, required
parser.add_argument("-c",help="configuration file",dest='configFile',required=True)
## input file, required
parser.add_argument("-i",help="input fasta file (contig or protein fasta/faa file)",dest='inputFile',required=True)


## output file, now this is removed, determined from the prefix argument instead
## parser.add_argument("-o",help="output annotation file",dest='outputFile',required=True)

## =================================================================
## main function
## =================================================================
def main(argv=None):
    
    if argv is None:
        args = parser.parse_args()
    
    ## print some information
    if args.verbose:
        sys.stdout.write('running verbosely\n')
        sys.stdout.write('configuration file is: {0}\n'.format(args.configFile))
        sys.stdout.write('input fasta file is: {0}\n'.format(args.inputFile))
        #sys.stdout.write('output file is: {}\n'.format(args.outputFile))
    else:
        sys.stdout.write('\n')

    # display work start, and time record
    start_time = datetime.now()
    sys.stderr.write("\n===============================================================================\n")
    sys.stderr.write("Welcome to UniFam v{0}: \n".format(version_str))

    # read configuration file
    config = ConfigParser.ConfigParser()
    config.read(args.configFile)


    # Annotating with UniFam
    if args.batch:
        UniFam_lib_batch.UniFam(args.inputFile,config,args.verbose)
    else:
        UniFam_lib.UniFam(args.inputFile,config,args.verbose)

    # write the configuration file to standard output for checking
    # config.write(sys.stdout)

    ## display work end, and time record, elapsed time
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write("\nTotal Elapsed Time = [%s] \n" % duration)
    sys.stderr.write("===============================================================================\n")

##==============================================================
## call from command line (instead of interactively)
##==============================================================

if __name__ == '__main__':
    sys.exit(main())
