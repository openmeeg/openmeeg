#! /usr/bin/python

from __future__ import print_function

import sys
import io
import os
import argparse
import subprocess
from shutil import rmtree

from lxml import etree
import StringIO

if sys.platform == 'win32':
    CMAKE_COMMAND = "C:\Program Files\CMake 2.8\\bin\cmake.exe"
    CTEST_COMMAND = "C:\Program Files\CMake 2.8\\bin\ctest.exe"
else:
    CMAKE_COMMAND = "cmake"
    CTEST_COMMAND = "ctest"

def CTest2Unit(xsl_file,build_dir,result_file):
    TestDir     = build_dir+'/Testing/'
    TAGFileName = TestDir+'TAG'
    if os.path.isfile(TAGFileName):
        TAGfile = open(TAGFileName,'r')
        dirname = TAGfile.readline().strip()

        xmlfile = open(TestDir+dirname+"/Test.xml",'r')
        xslfile = open(xsl_file,'r')

        xmlcontent = xmlfile.read()
        xslcontent = xslfile.read()

        xmldoc = etree.parse(StringIO.StringIO(xmlcontent))
        xslt_root = etree.XML(xslcontent)
        transform = etree.XSLT(xslt_root)

        result_tree = transform(xmldoc)
        with open(result_file,"w") as text_file:
            print(result_tree,file=text_file)

#   Call an external command and log the result.

def CallAndLog(command,logfile,debug=False):
    if debug:
        print(command)
    with open(logfile,'w') as f:
        process = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        for line in iter(process.stdout.readline, ''):
            sys.stdout.write(line)
            f.write(line)
        process.wait()

#   Build cmake command line and execute it

def add_cmake_parameter(boolean,base,cmd):
    str = '-D'+base+'='
    if boolean:
        cmd.append(str+"ON")
    else:
        cmd.append(str+"OFF")

def cmake_configuration(args):

    cmake_command_line = [ CMAKE_COMMAND ]
    if sys.platform == 'win32':
        cmake_command_line.extend(['-G', 'Visual Studio 11', '-DCMAKE_BUILD_TYPE=RelWithDebInfo'])

    add_cmake_parameter(args.python,'ENABLE_PYTHON',cmake_command_line)
    add_cmake_parameter(args.testing,'BUILD_TESTING',cmake_command_line)
    add_cmake_parameter(args.packaging,'ENABLE_PACKAGING',cmake_command_line)
    add_cmake_parameter(args.matlab,'MATLAB_TESTING',cmake_command_line)
    add_cmake_parameter(args.omp,'USE_OMP',cmake_command_line)

    cmake_command_line.append('..')
    CallAndLog(cmake_command_line,'configure.log',args.debug)

def cmake_build(args):

    cmake_command_line = [ CMAKE_COMMAND, '--build', '.' ]
    if sys.platform == 'win32':
        cmake_command_line.extend(['--config', 'RelWithDebInfo'])

    if args.incremental_build:
        cmcomline = list(cmake_command_line).extend(['--target','update'])
        call(cmcomline,args.debug)

    CallAndLog(cmake_command_line,'build.log',args.debug)

def cmake_test(test_type,args):
    cmake_command_line = [ CTEST_COMMAND, '-D', test_type]
    if sys.platform == 'win32':
        cmake_command_line.extend(['-C', 'RelWithDebInfo'])
    if test_type=='ExperimentalTest':
        cmake_command_line.append('--no-compress-output')

    CallAndLog(cmake_command_line,'test.log',args.debug)

basedir = os.path.dirname(os.path.realpath(__file__))

#   Parse arguments

parser = argparse.ArgumentParser(description='Options for building OpenMEEG.')
parser.add_argument('--debug',dest='debug',action='store_true',help='Debugging this script (print actions instead of executing them')
parser.add_argument('--incremental',dest='incremental_build',action='store_true',help='incremental build')
parser.add_argument('--disable-python',dest='python',action='store_false',help='disable python support')
parser.add_argument('--disable-testing',dest='testing',action='store_false',help='disable testing')
parser.add_argument('--disable-packaging',dest='packaging',action='store_false',help='disable packaging')
parser.add_argument('--enable-matlab-testing',dest='matlab',action='store_true',help='enable matio matlab comparison (requires matio build)')
parser.add_argument('--use-atlas',dest='atlas',action='store_true',help='use atlas library')
parser.add_argument('--use-mkl',dest='mkl',action='store_true',help='use mkl library')
parser.add_argument('--use-lapack',dest='lapack',action='store_true',help='use lapack library')
parser.add_argument('--use-OpenMP',dest='omp',action='store_true',help='use OpenMP acceleration')
args = parser.parse_args()

if args.lapack+args.mkl+args.atlas>1:
    sys.exit("ERROR: you can only use one of --use-atlas --use-mkl --use-lapack")

#   Create the proper build directory.

if not args.incremental_build and os.path.exists('build'):
    rmtree('build')

if os.path.isfile('build'):
    os.remove('build')

if not os.path.isdir('build'):
    os.mkdir('build')

os.chdir('build')

#   Configure the project with cmake if ut is a new one.

cmake_configuration(args)

#   Build the project

cmake_build(args)

#   Generate test reports and dashboards.

os.chdir('OpenMEEG/build')
if os.path.exists('JUnitTestResults.xml'):
    os.remove('JUnitTestResults.xml')

cmake_test('ExperimentalConfigure',args)
cmake_test('ExperimentalBuild',args)
cmake_test('ExperimentalTest',args)

CTest2Unit(basedir+'/CTest2JUnit.xsl','.','JUnitTestResults.xml')

# cdash backward compatibility.

cmake_test('ExperimentalSubmit',args)
