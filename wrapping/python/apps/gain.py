import os
import sys
import time
import argparse
import pathlib
import openmeeg as om

def parse_arguments(argv):
    parser = argparse.ArgumentParser(description='Compute various leadfields (gain matrices).')
    subparsers = parser.add_subparsers(dest='subparser_name', help='sub-command help')

    # create the parser for the EEG case
    EEGparser = subparsers.add_parser('EEG', help='Compute the EEG gain matrix')
    EEGparser.add_argument('HeadMatInv', type=open, nargs=1, help='file name containing the inverse of the head matrix')
    EEGparser.add_argument('SourceMat', type=open, nargs=1, help='file name containing the source matrix')
    EEGparser.add_argument('Head2EEGMat', type=open, nargs=1, help='file name containing the head to EEG conversion matrix')
    EEGparser.add_argument('dest_file', metavar='EEGGainMatrix', type=str, nargs=1, help='output: EEG gain matrix file name')

    # create the parser for the MEG case
    MEGparser = subparsers.add_parser('MEG', help='Compute the MEG gain matrix')
    MEGparser.add_argument('HeadMatInv', type=open, nargs=1, help='file name containing the inverse of the head matrix')
    MEGparser.add_argument('SourceMat', type=open, nargs=1, help='file name containing the source matrix')
    MEGparser.add_argument('Head2MEGMat', type=open, nargs=1, help='file name containing the head to MEG conversion matrix')
    MEGparser.add_argument('Source2MEGMat', type=open, nargs=1, help='file name containing the source to MEG conversion matrix')
    MEGparser.add_argument('dest_file', metavar='MEGGainMatrix', type=str, nargs=1, help='output: MEG gain matrix file name')

    # create the parser for the internal potential case
    IPparser = subparsers.add_parser('InternalPotential', aliases=['IP'], help='Compute the internal potential gain matrix')
    IPparser.add_argument('HeadMatInv', type=open, nargs=1, help='file name containing the inverse of the head matrix')
    IPparser.add_argument('SourceMat', type=open, nargs=1, help='file name containing the source matrix')
    IPparser.add_argument('Head2IPMat', type=open, nargs=1, help='file name containing the head to internal potential conversion matrix')
    IPparser.add_argument('Source2IPMat', type=open, nargs=1, help='file name containing the source to internal potential conversion matrix')
    IPparser.add_argument('dest_file', metavar='IPGainMatrix', type=str, nargs=1, help='output: EEG gain matrix file name')

    # create the parser for the EIT internal potential case
    EITIPparser = subparsers.add_parser('EITInternalPotential', aliases=['EITIP'], help='Compute the gain matrix for internal potentials using boundary normal current as input')
    EITIPparser.add_argument('HeadMatInv', type=open, nargs=1, help='file name containing the inverse of the head matrix')
    EITIPparser.add_argument('SourceMat', type=open, nargs=1, help='file name containing the source matrix')
    EITIPparser.add_argument('Head2IPMat', type=open, nargs=1, help='file name containing the head to internal potential conversion matrix')
    EITIPparser.add_argument('dest_file', metavar='EEGGainMatrix', type=str, nargs=1, help='output: EEG gain matrix file name')

    # create the parser for the EEG adjoint case
    EEGAdjointparser = subparsers.add_parser('EEGadjoint', help='Compute the EEG gain matrix')
    EEGAdjointparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    EEGAdjointparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    EEGAdjointparser.add_argument('sources', type=open, nargs=1, help='file name containing the description of the dipoles (positions, orientations)')
    EEGAdjointparser.add_argument('HeadMat', type=open, nargs=1, help='file name containing the head matrix')
    EEGAdjointparser.add_argument('Head2EEGMat', type=open, nargs=1, help='output: file name containing the head to EEG conversion matrix')
    EEGAdjointparser.add_argument('dest_file', metavar='EEGGainMatrix', type=str, nargs=1, help='output: EEG gain matrix file name')

    # create the parser for the MEG adjoint case
    MEGAdjointparser = subparsers.add_parser('MEGadjoint', help='Compute the MEG gain matrix')
    MEGAdjointparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    MEGAdjointparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    MEGAdjointparser.add_argument('sources', type=open, nargs=1, help='file name containing the description of the dipoles (positions, orientations)')
    MEGAdjointparser.add_argument('HeadMat', type=open, nargs=1, help='file name containing the head matrix')
    MEGAdjointparser.add_argument('Head2MEGMat', type=open, nargs=1, help='file name containing the head to MEG conversion matrix')
    MEGAdjointparser.add_argument('Source2MEGMat', type=open, nargs=1, help='file name containing the source to MEG conversion matrix')
    MEGAdjointparser.add_argument('dest_file', metavar='MEGGainMatrix', type=str, nargs=1, help='output: MEG gain matrix file name')

    # create the parser for the joint EEG MEG adjoint case
    EEGMEGAdjointparser = subparsers.add_parser('EEGMEGadjoint', help='Compute both the EEG and MEG gain matrices')
    EEGMEGAdjointparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    EEGMEGAdjointparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    EEGMEGAdjointparser.add_argument('sources', type=open, nargs=1, help='file name containing the description of the dipoles (positions, orientations)')
    EEGMEGAdjointparser.add_argument('HeadMat', type=open, nargs=1, help='file name containing the head matrix')
    EEGMEGAdjointparser.add_argument('Head2EEGMat', type=open, nargs=1, help='file name containing the head to EEG conversion matrix')
    EEGMEGAdjointparser.add_argument('Head2MEGMat', type=open, nargs=1, help='file name containing the head to MEG conversion matrix')
    EEGMEGAdjointparser.add_argument('Source2MEGMat', type=open, nargs=1, help='file name containing the source to MEG conversion matrix')
    EEGMEGAdjointparser.add_argument('EEG_dest_file', metavar='EEGGainMatrix', type=str, nargs=1, help='output: EEG gain matrix file name')
    EEGMEGAdjointparser.add_argument('MEG_dest_file', metavar='MEGGainMatrix', type=str, nargs=1, help='output: MEG gain matrix file name')

    args = parser.parse_args(argv[1:])
    args.processing_function = args.subparser_name

    return args

def EEG(args):
    HeadMatInv  = om.SymMatrix(args.HeadMatInv)
    Head2EEGMat = om.SparseMatrix(args.Head2EEGMat)
    tmp         = Head2EEGMat*HeadMatInv
    SourceMat   = om.Matrix(args.SourceMat)
    EEGGainMat  = tmp*SourceMat
    EEGGainMat.save(args.dest_file)

def EEGadjoint(args):
    geometry    = om.Geometry(args.geometry,args.conductivities)
    dipoles     = om.Matrix(args.sources)
    HeadMat     = om.SymMatrix(args.HeadMat)
    Head2EEGMat = om.SparseMatrix(args.Head2EEGMat)
    EEGGainMat  = om.GainEEGadjoint(geometry, dipoles, HeadMat, Head2EEGMat)
    EEGGainMat.save(args.dest_file)

def MEG(args):
    HeadMatInv    = om.SymMatrix(args.HeadMatInv)
    Head2MEGMat   = om.Matrix(args.Head2MEGMat)
    tmp1          = Head2MEGMat*HeadMatInv
    SourceMat     = om.Matrix(args.SourceMat)
    tmp2          = tmp1*SourceMat
    Source2MEGMat = om.Matrix(args.Source2MEGMat)
    MEGGainMat    = Source2MEGMat+tmp2
    MEGGainMat.save(args.dest_file)

def MEGadjoint(args):
    geometry      = om.Geometry(args.geometry,args.conductivities)
    dipoles       = om.Matrix(args.sources)
    HeadMat       = om.SymMatrix(args.HeadMat)
    Head2MEGMat   = om.Matrix(args.Head2MEGMat)
    Source2MEGMat = om.Matrix(args.Source2MEGMat)
    MEGGainMat    = om.GainMEGadjoint(geometry, dipoles, HeadMat, Head2MEGMat, Source2MEGMat)
    MEGGainMat.save(args.dest_file)

def EEGMEGadjoint(args):
    geometry      = om.Geometry(args.geometry,args.conductivities)
    dipoles       = om.Matrix(args.sources)
    HeadMat       = om.SymMatrix(args.HeadMat)
    Head2EEGMat   = om.SparseMatrix(args.Head2EEGMat)
    Head2MEGMat   = om.Matrix(args.Head2MEGMat)
    Source2MEGMat = om.Matrix(args.Source2MEGMat)
    EEGMEGGainMat = om.EEGMEGGainMat(geometry, dipoles, HeadMat, Head2EEGMat, Head2MEGMat, Source2MEGMat)
    EEGMEGGainMat.saveEEG(args.EEG_dest_file)
    EEGMEGGainMat.saveMEG(args.MEG_dest_file)

def InternalPotential(args):
    HeadMatInv         = om.SymMatrix(args.HeadMatInv)
    Head2IPMat         = om.SparseMatrix(args.Head2IPMat)
    tmp1               = Head2IPMat*HeadMatInv
    SourceMat          = om.Matrix(args.SourceMat)
    tmp2               = tmp1*SourceMat
    Source2IPMat       = om.Matrix(args.Source2IPMat)
    InternalPotGainMat = Source2IPMat+tmp2
    InternalPotGainMat.save(args.dest_file)

def EITInternalPotential(args):
    HeadMatInv         = om.SymMatrix(args.HeadMatInv)
    SourceMat          = om.Matrix(args.SourceMat)
    Head2IPMat         = om.SparseMatrix(args.Head2IPMat)
    InternalPotGainMat = (Head2IPMat*HeadMatInv)*SourceMat
    InternalPotGainMat.save(args.dest_file)
