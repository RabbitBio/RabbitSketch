import os
import pathlib
import platform
import re
import subprocess
import sys

import pybind11
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


ROOT = pathlib.Path(__file__).resolve().parent
VERSION = "2.0.0"


class CMakeExtension(Extension):
    def __init__(self, name, source_directory=""):
        super().__init__(name, sources=[])
        self.source_directory = str(pathlib.Path(source_directory).resolve())


class CMakeBuild(build_ext):
    def run(self):
        try:
            output = subprocess.check_output(["cmake", "--version"], text=True)
        except (OSError, subprocess.CalledProcessError) as error:
            raise RuntimeError("CMake >= 3.16 is required to build rabbitsketch") from error
        match = re.search(r"version\s+(\d+)\.(\d+)(?:\.(\d+))?", output)
        if not match or tuple(int(part or 0) for part in match.groups()) < (3, 16, 0):
            raise RuntimeError("CMake >= 3.16 is required to build rabbitsketch")
        super().run()

    def build_extension(self, extension):
        output_directory = pathlib.Path(
            self.get_ext_fullpath(extension.name)
        ).resolve().parent
        configuration = "Debug" if self.debug else "Release"
        cmake_arguments = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={output_directory}{os.sep}",
            f"-DPython3_EXECUTABLE={sys.executable}",
            f"-DCMAKE_PREFIX_PATH={sys.prefix}",
            f"-Dpybind11_DIR={pybind11.get_cmake_dir()}",
            "-DCXXAPI=OFF",
            "-DRABBITSKETCH_BUILD_TESTS=OFF",
            "-DRABBITSKETCH_NATIVE_ARCH=OFF",
        ]
        build_arguments = ["--config", configuration]
        if platform.system() == "Windows":
            cmake_arguments.append(
                f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{configuration.upper()}="
                f"{output_directory}{os.sep}"
            )
            if sys.maxsize > 2**32:
                cmake_arguments.extend(["-A", "x64"])
        else:
            cmake_arguments.append(f"-DCMAKE_BUILD_TYPE={configuration}")
        parallel = self.parallel or min(2, os.cpu_count() or 1)
        build_arguments.extend(["--parallel", str(parallel)])

        build_directory = pathlib.Path(self.build_temp) / extension.name
        build_directory.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(
            ["cmake", "-S", extension.source_directory, "-B", str(build_directory),
             *cmake_arguments],
        )
        subprocess.check_call(
            ["cmake", "--build", str(build_directory), *build_arguments],
        )


setup(
    name="rabbitsketch",
    version=VERSION,
    author="Zekun Yin",
    author_email="zekun.yin@mail.sdu.edu.cn",
    description="High-performance genomic sketching and unified queries",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.8",
    ext_modules=[CMakeExtension("rabbitsketch", ROOT)],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
