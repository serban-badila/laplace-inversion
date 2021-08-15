import os
import subprocess
import shlex
import shutil

from pathlib import Path

from setuptools import setup  # a problem with the setup override from skbuild is that it wants to install
# everything it can find in the CMakeLists.txt (nested also) which is not a good idea in this case
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import CMaker
from skbuild.constants import CMAKE_BUILD_DIR
from skbuild.cmaker import pop_arg


BUILD_TARGET = "_py_laplace"
PYTHON_MODULE = "py_laplace_inversion"


class PatchedCMaker(CMaker):
    """Relax the CMaker generator target specification.
    
    TODO: try get rid of this once https://github.com/scikit-build/scikit-build/pull/477 is merged
    """

    def make(self, clargs=(), config="Release", source_dir=".", env=None):
        """Calls the system-specific make program to compile code.
        """
        clargs, config = pop_arg('--config', clargs, config)
        if not os.path.exists(CMAKE_BUILD_DIR()):
            raise SKBuildError(("CMake build folder ({}) does not exist. "
                                "Did you forget to run configure before "
                                "make?").format(CMAKE_BUILD_DIR()))

        cmd = [self.cmake_executable, "--build", source_dir,
               "--target", BUILD_TARGET, "--config", config, "--"]
        cmd.extend(clargs)
        cmd.extend(
            filter(bool,
                   shlex.split(os.environ.get("SKBUILD_BUILD_OPTIONS", "")))
        )

        rtn = subprocess.call(cmd, cwd=CMAKE_BUILD_DIR(), env=env)
        if rtn != 0:
            raise SKBuildError(
                "An error occurred while building with CMake.\n"
                "  Command:\n"
                "    {}\n"
                "  Source directory:\n"
                "    {}\n"
                "  Working directory:\n"
                "    {}\n"
                "Please see CMake's output for more informationpy_laplace_inversion"
                )
readme_path = Path("README.md")
long_description = readme_path.read_text(encoding="utf-8")
version_path = Path("VERSION")
version = version_path.read_text(encoding="utf-8").strip()
license_path = Path("LICENSE")
license = license_path.read_text(encoding="utf-8").strip()

cmaker = PatchedCMaker()
cmaker.configure()
cmaker.make()

pp = Path()
compiled_extension = list(pp.glob(f"**/{BUILD_TARGET}.*.*")).pop()
shutil.copyfile(compiled_extension.absolute(), pp.cwd() / PYTHON_MODULE / compiled_extension.name)

setup(
    name="py-laplace-inversion",
    version=version,
    description="Numerical inversion of Laplace transforms using Den Iseger's algorithm.",
    long_description=long_description,
    url="https://github.com/serban-badila/laplace-inversion",
    author="Serban Badila",
    author_email="badilaserban@gmail.com",
    license=license,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD 3-Clause",
        "Programming Language :: Python :: 3.7",
    ],
    packages=['py_laplace_inversion'],
    package_dir={'py_laplace_inversion': 'py_laplace_inversion'},
    package_data={'': [compiled_extension.name]},
)
