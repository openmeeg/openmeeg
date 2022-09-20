import os
import sys
import time
import argparse
import pathlib
import openmeeg as om

def parse_arguments(argv):
    parser = argparse.ArgumentParser(description='Compute the forward problem given a leadfield (gain matrix) and sources.')
    parser.add_argument('leadfield', metavar='gain_matrix', type=open, nargs=1, help='file name of the gain matrix')
    parser.add_argument('sources', metavar='source_matrix', type=open, nargs=1, help='file name of the source matrix')
    parser.add_argument('dest_file', metavar='simulated_data', type=pathlib.Path, nargs=1, help='output: simulated data file name')
    parser.add_argument('noise_level', type=float, help='noise level to add on simulated data')

    args = parser.parse_args(argv[1:])
    args.processing_function = 'forward'
    return args

def forward(args):
    leadfield = om.Matrix(args.leadfield)
    sources = om.Matrix(args.sources)
    noise = args.noise_level

    simulated_data = om.Forward(leadfield, sources, noise)
    simulated_data.save(args.dest_file)
