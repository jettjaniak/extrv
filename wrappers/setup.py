import platform
from distutils.core import setup, Extension
from Cython.Build import cythonize

extra_compile_args = []
system = platform.system()
if system == "Linux":
    extra_compile_args.append('-std=c++11')
# TODO: macOS

sources = ["Ligand.cpp", "SimulationState.cpp", "forces.cpp",
           "helpers.cpp", "interpolated.cpp", "velocities.cpp"]
sources = [f"../engine/{s}" for s in sources]

setup(ext_modules=cythonize(Extension(
    "engine",
    sources=["engine.pyx"] + sources,
    extra_compile_args=extra_compile_args,
    language="c++"
), language_level='3', annotate=False))
