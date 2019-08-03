#! /usr/bin/env python

"""
Front-end script for skeletonizing a p-map stack represented as an .h5 file.
"""

import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description='Skeletonize a p-map stack represented as an .h5 file')
parser.add_argument('sample_date', help='Date of the sample to use, e.g. "2019-05-27"')
parser.add_argument('--size_threshold', nargs='?', help='Size threshold to use', type=int, default=300)
parser.add_argument('--probability_threshold', nargs='?', help='Probability threshold, as unsigned 8-bit integer', type=int, default=100)
parser.add_argument('--fullh', nargs='?', help='not sure what this does', type=int, default=15)
args = parser.parse_args()

matlab_command_line_template = 'try; modpath; submit_bsub_jobs_for_skeletonization(\'%s\', %d, %d, %d); catch err; fprintf(2, \'%%s\\n\', err.getReport()); quit(1); end; quit(0);'
print(matlab_command_line_template)
matlab_command_line = (matlab_command_line_template % (args.sample_date, args.size_threshold, args.probability_threshold, args.fullh))
print(matlab_command_line)

script_file_path = os.path.abspath(__file__)
script_folder_path = os.path.dirname(script_file_path)
os.chdir(script_folder_path)
child = subprocess.Popen(['/misc/local/matlab-2018b/bin/matlab', '-nodisplay', '-r', matlab_command_line])
child.communicate()
rc = child.returncode
sys.exit(rc)
