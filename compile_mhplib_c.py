from distutils.core import setup
from Cython.Build import cythonize

setup(
        ext_modules = cythonize('mhplib_c.pyx')
)

setup(
        ext_modules = cythonize('sas.pyx')
)
