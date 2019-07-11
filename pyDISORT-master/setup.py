#!/usr/bin/env python
# ------------------------------------------------------
# Adapted from setup.py from CliMT project
# ------------------------------------------------------

import os
import sys
import glob
import string
from numpy.distutils.core import setup
from numpy.distutils import fcompiler
from distutils.dep_util import newer

from numpy import f2py
f2py_path = f2py.__path__[0]
direct = f2py_path
while 1:
    base = os.path.basename(direct)
    if base=='':
        print '*f2py* not found in your python installation'
        raise SystemExit(0)
    direct = os.path.dirname(direct)
    if base=='lib' or base=='lib64':
        root = direct
        if 'bin' in os.listdir(root):
            f2py_exec_path = os.path.join(root,'bin')
            if 'f2py' in os.listdir(f2py_exec_path):
                f2py_bin = os.path.join(f2py_exec_path,'f2py')
                break

Extensions = [
    {'name':'disort',
     'dir':'src/disort'},
    ]

# figure out which compiler we're goint to use
compiler = fcompiler.get_default_fcompiler(requiref90=True)
compiler = 'gnu95'  # ensure gfortran
for i in range(len(sys.argv)):
    if '--fcompiler' in sys.argv[i]:
        compiler = sys.argv.pop(i)
        compiler = compiler[compiler.index('=')+1:]
print 'Using %s compiler' % compiler

# set some fortran compiler-dependent flags
if compiler == 'gnu95' or compiler == 'gnu':
    # f77flags='-ffixed-line-length-132 -fdefault-real-8'
    # f90flags='-fdefault-real-8'
    f77flags='-O3'
    f90flags='-O3'
elif compiler == 'intel' or compiler == 'intelem':
    f77flags='-132 -r8 -w95 -w90 -mp'
    f90flags='-r8 -w95 -mp'
elif compiler == 'ibm':
    f77flags='-qautodbl=dbl4 -qsuffix=f=f:cpp=F -qfixed=132'
    f90flags='-qautodbl=dbl4 -qsuffix=f=f90:cpp=F90 -qfree=f90'
else:
    print 'Sorry, compiler %s not supported' % compiler
    sys.exit()

for i in range(len(Extensions)):
    Defaults = {  # 'cppflags':'-DIM=%i -DJM=%i -DKM=%i' % (IM,JM,KM),
                'cppflags':'',
                'f77flags':f77flags,
                'f90flags':f90flags}
    Defaults.update(Extensions[i])
    Extensions[i] = Defaults

if compiler == 'ibm':
    for ext in Extensions:
        ext['cppflags']='-WF,'+string.join(ext['cppflags'].split(),',')

def getSources(dir):
    # Gets list of source files for extensions
    SrcFile = os.path.join(dir,'sources_in_order_of_compilation')
    if os.path.exists(SrcFile):
        Sources = open(SrcFile).readlines()
        Sources = [os.path.join(dir,s[:-1]) for s in Sources]
    else:
        Sources = []
        for pattern in ['*.f','*.F','*.f90','*.F90']:
            Sources += glob.glob(os.path.join(dir,'src',pattern))  # sgg: i changed to order
            Sources += glob.glob(os.path.join(dir,pattern))        #
    return Sources

def buildNeeded(target,src):
    # Checks if source code is newer than extension, so extension needs to be rebuilt
    target = os.path.join('lib/disort',target)
    if not os.path.exists(target):
        return True
    for file in src:
        if newer(file,target):
            return True
    print 'Extension %s is up to date' % os.path.basename(target)
    return False

def build_ext(name=None, dir=None, cppflags='', f77flags='', f90flags='',
              lib='', libdir='', incdir=''):
    # Builds an extension
    src = getSources(dir)
    target = '_%s.so' % name
    driver = glob.glob(os.path.join(dir,'Driver.f'))[0]
    f77flags = '"%s %s"' % (cppflags,f77flags)
    f90flags = '"%s %s"' % (cppflags,f90flags)
    if buildNeeded(target,src):
        print '\n Building %s ... \n' % os.path.basename(target)
        # generate signature file
        #os.system('f2py --overwrite-signature %s -m _%s -h _%s.pyf'%(driver,name,name))
        ff = '%s '*len(src)
        #sformat = '/usr/bin/f2py --overwrite-signature '+ff+' -m _%s -h _%s.pyf'
        sformat = f2py_bin+' --overwrite-signature '+ff+' -m _%s -h _%s.pyf'
        args = src + [name,name]
        print sformat % tuple(args)
        os.system(sformat % tuple(args))
        # compile extension
        F2pyCommand = []
#        F2pyCommand.append('/usr/bin/f2py -c -m _%s' % name)
        F2pyCommand.append('%s -c ' % f2py_bin)
        F2pyCommand.append('--fcompiler=%s' % compiler)
        F2pyCommand.append('-I%s' % dir)
        F2pyCommand.append('-I%s' % os.path.join(dir,'include'))
        F2pyCommand.append('-I%s' % os.path.join(dir,'src'))
        F2pyCommand.append('-I%s' % os.path.join(dir,'src','include'))
        if incdir is not '':
            for i in incdir:
                F2pyCommand.append('-I%s' % i)
        if libdir is not '':
            for i in libdir:
                F2pyCommand.append('-L%s' % i)
        if lib is not '':
            for i in lib:
                F2pyCommand.append('-l%s' % i)
        F2pyCommand.append('--f77flags=%s' % f77flags)
        F2pyCommand.append('--f90flags=%s' % f90flags)
        F2pyCommand.append('_%s.pyf' % name)
        F2pyCommand.append('%s' % string.join(src))
        F2pyCommand = string.join(F2pyCommand)
        print F2pyCommand
        if os.system(F2pyCommand) > 0:
            print '+++ Compilation failed'
            sys.exit()
        os.system('mv -f _%s.so lib/disort' % name)
        # os.system('rm -f _%s.pyf' % name)

# Build all extensions
for ext in Extensions:
    build_ext(**ext)

# Finish the setup
# note: setup() cannot copy directories, and falls over
# trying to copy the CVS directory in climt/lib/data
# workaround: make data list which specifically excludes CVS
os.chdir('lib/')
DataFiles = []
for File in glob.glob('test/*.py'):
    if 'CVS' not in File:
        DataFiles.append('../'+File)
print DataFiles
os.chdir('..')

setup(name         = "disort",
      version      = open('Version').read()[:-1],
      description  = "Python wrapper to the DISORT library",
      author       = "Sebati\'an Gimeno Garc\'\ia",
      author_email = "sebastian.gimenogarcia@gmail.com",
      packages     = ['disort'],
      package_dir  = {'':'lib'},
      package_data = {'disort':['*.so']+DataFiles})
