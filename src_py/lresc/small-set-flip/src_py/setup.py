#from distutils.core import setup
#from setuptools import Extension
#from Cython.Compiler.Options import get_directive_defaults
#directive_defaults = get_directive_defaults()
#
#directive_defaults['linetrace'] = True
#directive_defaults['binding'] = True
#
#from Cython.Build import cythonize
#
#extensions = [
#    Extension("decoder", ["decoder.pyx"], define_macros=[('CYTHON_TRACE', '1')])
#]
#
#setup(
#    ext_modules = cythonize("decoder.pyx")
#)

from distutils.extension import Extension
from distutils.core import setup
from Cython.Build import cythonize

extensions = [
    Extension("decoder", ["decoder.pyx"], define_macros=[('CYTHON_TRACE', '1')]),
]

setup(
    ext_modules = cythonize(extensions)
)


# setup(ext_modules=cythonize('decoder.pyx'))

