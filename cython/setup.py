## cython build file

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import glob
import numpy
import os

def findfiles(directory):
    files = glob.glob(directory+"/*.pyx")
    print "checking directory: %s ..." %directory
    filename = []
    for eachfile in files:
        eachfile.strip()
        item = eachfile.split("/")
        filename.append(item[-1][:-4])
        print item[-1]
    return filename, files

direct = "./package"
name, files = findfiles(direct)
numpy_inc = [numpy.get_include()]

ext_modules = []

for i in range(len(name)):
    ext_modules.append( Extension(name[i], [files[i]], include_dirs=numpy_inc) )

setup(name="postprocess", cmdclass={"build_ext":build_ext}, ext_modules=ext_modules)

## remove generated c file
cfile = glob.glob(direct + "/*.c")
if len(cfile)>0:
    for eachfile in cfile:
        os.remove(eachfile)
