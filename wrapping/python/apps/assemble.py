import sys
import argparse
#import openmeeg as om

def parse_arguments(argv):
    parser = argparse.ArgumentParser(prog=argv[0],description='Compute various head matrices).')
    parser.add_argument('--old-ordering', type=bool, default=False, help='Using old ordering i.e using (V1, p1, V2, p2, V3) instead of (V1, V2, V3, p1, p2)')
    subparsers = parser.add_subparsers(dest='subparser_name', help='sub-command help')

    # create the parser for the HeadMat case
    HMparser = subparsers.add_parser('HeadMat', aliases=['HM', 'hm'], help='Compute Head Matrix for Symmetric BEM (left-hand side of the linear system).')
    HMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    HMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    HMparser.add_argument('dest_file', metavar='HeadMatrix', type=str, nargs=1, help='output: head matrix file name')

    # create the parser for the CorticalMat case
    CMparser = subparsers.add_parser('CorticalMat', aliases=['CM', 'cm'], help='Compute Cortical Matrix for Symmetric BEM (left-hand side of the linear system).')
    CMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    CMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    CMparser.add_argument('electrodes', type=open, nargs=1, help='EEG electrode positions file name (.patches)')
    CMparser.add_argument('domain', type=open, nargs=1, help='domain name (containing the sources)')
    CMparser.add_argument('dest_file', metavar='CorticalMatrix', type=str, nargs=1, help='output: matrix file name')
    CMparser.add_argument('options', type=str, nargs='*', help='alpha, gamma, filename')

    # create the parser for the SurfSourceMat case
    SSMparser = subparsers.add_parser('SurfSourceMat', aliases=['SSM', 'ssm'], help='Compute Surfacic Source Matrix for Symmetric BEM (right-hand side of the linear system).')
    SSMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    SSMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    SSMparser.add_argument('sources', type=open, nargs=1, help='sources mesh file name (.tri .vtk .mesh .bnd)')
    SSMparser.add_argument('dest_file', metavar='SSMMatrix', type=str, nargs=1, help='output: matrix file name')

    # create the parser for the DipSourceMat case
    DSMparser = subparsers.add_parser('DipSourceMat', aliases=['DSM', 'dsm'], help='Compute Dipolar Source Matrix for Symmetric BEM (right-hand side of the linear system).')
    DSMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    DSMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    DSMparser.add_argument('dipoles', type=open, nargs=1, help='dipoles file name (.dip)')
    DSMparser.add_argument('dest_file', metavar='DSMMatrix', type=str, nargs=1, help='output: matrix file name')
    DSMparser.add_argument('domain', type=str, nargs='?', help='domain name where lie all dipoles')

    # create the parser for the DipSourceMatNoAdapt case
    DSMNAparser = subparsers.add_parser('DipSourceMatNoAdapt', aliases=['DSMNA', 'dsmna'], help='Compute Dipolar Source Matrix for Symmetric BEM (right-hand side of the linear system).')
    DSMNAparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    DSMNAparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    DSMNAparser.add_argument('dipoles', type=open, nargs=1, help='dipoles file name (.dip)')
    DSMNAparser.add_argument('dest_file', metavar='DSMMatrix', type=str, nargs=1, help='output: matrix file name')
    DSMNAparser.add_argument('domain', type=str, nargs='?', help='domain name where lie all dipoles')

    # create the parser for the EITSourceMat case
    EITSMparser = subparsers.add_parser('EITSourceMat', aliases=['EITSM', 'eitsm'], help='Compute the EIT Source Matrix from an injected current (right-hand side of the linear system).')
    EITSMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    EITSMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    EITSMparser.add_argument('electrodes', type=open, nargs=1, help='EIT electrode positions file name (.patches)')
    EITSMparser.add_argument('dest_file', metavar='EITSMMatrix', type=str, nargs=1, help='output: matrix file name')

    # create the parser for the Head2EEGMat case
    H2EMparser = subparsers.add_parser('Head2EEGMat', aliases=['H2EM', 'h2em'], help='Compute the linear application which maps the potential on the scalp to the EEG electrodes.')
    H2EMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    H2EMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    H2EMparser.add_argument('electrodes', type=open, nargs=1, help='EEG electrode positions file name (.patches)')
    H2EMparser.add_argument('dest_file', metavar='Head2EEGMat', type=str, nargs=1, help='output: matrix file name')

    # create the parser for the Head2ECoGMat case
    H2ECOGMparser = subparsers.add_parser('Head2ECoGMat', aliases=['H2ECogM', 'h2ecogm', 'H2ECOGM'], help='Compute the linear application which maps the potential on the scalp to the ECoG electrodes.')
    H2ECOGMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    H2ECOGMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    H2ECOGMparser.add_argument('electrodes', type=open, nargs=1, help='ECoG electrode positions file name (.patches)')
    H2ECOGMparser.add_argument('interface', type=open, nargs=1, help='name of the interface on which to project the electrodes (\"name\")')
    H2ECOGMparser.add_argument('dest_file', metavar='Head2ECoGMat', type=str, nargs=1, help='output: matrix file name')

    # create the parser for the Head2MEGMat case
    H2MMparser = subparsers.add_parser('Head2MEGMat', aliases=['H2MM', 'h2mm'], help='Compute the linear application which maps the potential on the scalp to the MEG sensors.')
    H2MMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    H2MMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    H2MMparser.add_argument('sensors', type=open, nargs=1, help='MEG sensor positions and orientations file name (.squids)')
    H2MMparser.add_argument('dest_file', metavar='Head2MEGMat', type=str, nargs=1, help='output: matrix file name')

    # create the parser for the SurfSource2MEGMat case
    SS2MMparser = subparsers.add_parser('SurfSource2MEGMat', aliases=['SS2MM', 'ss2mm'], help='Compute the linear mapping from distributed source to MEG sensors.')
    SS2MMparser.add_argument('sources', type=open, nargs=1, help='sources mesh file name (.tri .vtk .mesh .bnd)')
    SS2MMparser.add_argument('sensors', type=open, nargs=1, help='MEG sensor positions and orientations file name (.squids)')
    SS2MMparser.add_argument('dest_file', metavar='SurfSource2MEGMat', type=str, nargs=1, help='output: matrix file name')

    # create the parser for the DipSource2MEGMat case
    DS2MMparser = subparsers.add_parser('DipSource2MEGMat', aliases=['DS2MM', 'ds2mm'], help='Compute the linear mapping from current dipoles to MEG sensors.')
    DS2MMparser.add_argument('dipoles', type=open, nargs=1, help='dipoles file name (.dip)')
    DS2MMparser.add_argument('sensors', type=open, nargs=1, help='MEG sensor positions and orientations file name (.squids)')
    DS2MMparser.add_argument('dest_file', metavar='DipSource2MEGMat', type=str, nargs=1, help='output: matrix file name')

    # create the parser for the Head2InternalPotMat case
    H2IPMparser = subparsers.add_parser('Head2InternalPotMat', aliases=['H2IPM', 'h2ipm'], help='Compute the linear mapping from surface potential and normal current to internal potential at a set of points.')
    H2IPMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    H2IPMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    H2IPMparser.add_argument('points', type=open, nargs=1, help='mesh file name or point positions file name')
    H2IPMparser.add_argument('dest_file', metavar='Head2InternalPotMat', type=str, nargs=1, help='output: matrix file name')

    # create the parser for the DipSource2InternalPotMat case
    DS2IPMparser = subparsers.add_parser('DipSource2InternalPotMat', aliases=['DS2IPM', 'ds2ipm'], help='Compute the linear mapping from current dipoles to infinite potential at a set of points.')
    DS2IPMparser.add_argument('geometry', type=open, nargs=1, help='head geometry file name (.geom)')
    DS2IPMparser.add_argument('conductivities', type=open, nargs=1, help='head conductivities file name (.cond)')
    DS2MMparser.add_argument('dipoles', type=open, nargs=1, help='dipoles file name (.dip)')
    DS2IPMparser.add_argument('points', type=open, nargs=1, help='mesh file name or point positions file name')
    DS2IPMparser.add_argument('dest_file', metavar='SSMMatrix', type=str, nargs=1, help='output: matrix file name')
    DS2IPMparser.add_argument('domain', type=open, nargs='?', help='domain name where lie all dipoles')

    args = parser.parse_args(argv[1:])
    args.processing_function = args.subparser_name
    if args.subparser_name=='DipSourceMatNoAdapt':
        args.processing_function = 'DipSourceMat'

    return args

def HeadMat(args):
    geometry = om.Geometry(args.geometry,args.conductivities, args.old_ordering)
    if not geometry.selfCheck():
        sys.exit(1);

    HeadMat = om.HeadMat(geometry)
    HeadMat.save(args.dest_file)

def CorticalMat(args):
    numopt = len(args.options)
    filename = ''
    alpha = -1.0
    beta  = -1.0
    gamma = -1.0
    if numopt==1:
        try:
            gamma = float(args.options[0])
        except ValueError:
            filename = args.options[0]

    elif numopt==2:
        alpha = float(args.options[0])
        try:
            beta = float(args.options[1])
        except ValueError:
            filename = args.options[1]
            gamma = alpha
            
    elif numopt==3:
        alpha = float(args.options[0])
        beta  = float(args.options[1])
        filename = args.options[2]
        
    geometry = om.Geometry(args.geometry,args.conductivities, args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    electrodes = on.Sensors(args.electrodes)
    M          = om.Head2EEGMat(geometry, electrodes)

    if gamma>0:
        CM = om.CorticalMat2(geometry, M, args.domain, gamma, filename)
    else:
        CM = om.CorticalMat2(geometry, M, args.domain, gamma, alpha, beta, filename)

    CM.save(args.dest_file)

def SurfSourceMat(args):

    geometry = om.Geometry(args.geometry,args.conductivities,args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    mesh_sources = om.Mesh(args.sources)
    ssm          = SurfSourceMat(geometry,mesh_sources);

    ssm.save(args.dest_file)

def DipSourceMat(args):

    domain = ''
    if args.domain:
        domain = args.domain
    print(domain)

    geometry = om.Geometry(args.geometry,args.conductivities,args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    dipoles = om.Matrix(args.dipoles)
    if dipoles.ncol()!=6:
        print('Dipoles File Format Error')
        sys.exit(1)

    integration_levels = 10
    if args.subparser_name=='DipSourceMatNoAdapt':
        integration_levels = 0

    dsm = om.DipSourceMat(geometry, dipoles, om.Integrator(3,integration_levels,0.001))
    dsm.save(args.dest_file)

def EITSourceMat(args):
    
    geometry = om.Geometry(args.geometry,args.conductivities,args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    electrodes = om.Sensors(args.electrodes)
    EITsource = om.EITSourceMat(geometry,electrodes)
    EITsource.save(args.dest_file)

def Head2EEGMat(args):
    
    geometry = om.Geometry(args.geometry,args.conductivities,args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    electrodes = om.Sensors(args.electrodes)
    EEGmat = om.Head2EEGMat(geometry,electrodes)
    EEGmat.save(args.dest_file)

def Head2ECoGMat(args):
    
    # TODO: The C++ interface supports that the interface argument is omitted, which
    # is only meaningful for nested geometries (in which case the inner layer is used).
    # Do we want to continue to support that ?

    geometry = om.Geometry(args.geometry,args.conductivities,args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    electrodes = om.Sensors(args.electrodes)
    ECoG_layer = geometry.interface(args.interface)
    ECoGmat = om.Head2ECoGMat(geometry,electrodes,ECoG_layer)
    ECoGmat.save(args.dest_file)

def Head2MEGMat(args):
    
    geometry = om.Geometry(args.geometry,args.conductivities,args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    sensors = om.Sensors(args.sensors)
    Head2MEGmat = om.Head2MEGMat(geo,sensors)
    Head2MEGmat.save(args.dest_file)

def SurfSource2MEGMat(args):
    
    dipoles = om.Matrix(args.dipoles)
    sensors = om.Sensors(args.sensors)
    
    Source2MEGmat = om.SurfSource2MEGMat(geo,sensors)
    Source2MEGmat.save(args.dest_file)

def DipSource2MEGMat(args):
    
    dipoles = om.Matrix(args.dipoles)
    sensors = om.Sensors(args.sensors)
    
    DipSource2MEGMat = om.DipSource2MEGMat(geo,sensors)
    DipSource2MEGMat.save(args.dest_file)

def Head2InternalPotMat(args):
    
    geometry = om.Geometry(args.geometry,args.conductivities,args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    points = om.Matrix(args.points)
    Head2InternalPotMat = om.Surf2VolMat(geometry, points)
    Head2InternalPotMat.save(args.dest_file)

def DipSource2InternalPotMat(args):
    
    geometry = om.Geometry(args.geometry,args.conductivities,args.old_ordering)

    if not geometry.selfCheck():
        sys.exit(1);

    dipoles = om.Matrix(args.dipoles)
    points  = om.Matrix(args.points)

    domain_name = ''
    if args.domain!=None:
        domain_name = args.domain
        print("Dipoles are considered to be in \""+domain_name+"\" domain.")

    DipSource2InternalPotMat = om.DipSource2InternalPotMat(geometry, dipoles, points, domain_name)
    DipSource2InternalPotMat.save(args.dest_file)
