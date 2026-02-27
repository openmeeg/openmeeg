import os
import sys
import time
import argparse
import pathlib
import openmeeg as om

def parse_arguments(argv):
    parser = argparse.ArgumentParser(description='Compute the inverse of an head matrix.')
    parser.add_argument('head_matrix', type=open, nargs=1, help='head matrix file name')
    parser.add_argument('dest_file', type=pathlib.Path, nargs=1, help='output: inverse head matrix file name')

    args = parser.parse_args(argv[1:])
    args.processing_function = 'invert'

    return args

def invert(args):
    head_matrix = om.SymMatrix(args.head_matrix)
    head_matrix.invert()
    head_matrix.save(args.dest_file)
