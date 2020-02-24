import platform
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extra_compile_args = []
system = platform.system()
if system == "Linux":
    extra_compile_args.append('-std=c++11')
# TODO: macOS

sources = ["Ligand.cpp", "Settings.cpp", "SimulationState.cpp",
           "forces.cpp", "helpers.cpp", "interpolated.cpp", "velocities.cpp"]
sources = [f"../engine/sources/{s}" for s in sources]

setup(ext_modules=cythonize(Extension(
    "engine",
    sources=["engine.pyx"] + sources,
    include_dirs=[np.get_include()],
    extra_compile_args=extra_compile_args,
    language="c++"
), language_level='3', annotate=False))
