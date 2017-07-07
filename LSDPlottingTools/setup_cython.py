from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules = [
    Extension(
        "fast_hillshade",
        ["fast_hillshade.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
]

setup(
    name='fast_hillshade',
    ext_modules=cythonize(ext_modules),
)
# build with python setup_cython.py build_ext --inplace
