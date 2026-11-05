#!/usr/bin/env python3
"""Build an isolated accuracy adapter against the frozen corrected native source."""
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = Path('/home/bigssd/zt/rabbitsketch2_benchmark')
SOURCE = ROOT / 'source/RabbitSketch-validated'
BUILD = ROOT / 'build_comprehensive_20261030'


def main():
    BUILD.mkdir(exist_ok=True)
    original = (HERE.parent / 'regimes/main_regimes.cpp').read_text()
    needle = '}else if(c.algorithm=="probminhash"){'
    assert original.count(needle) == 1
    replacement = '''}else if(c.algorithm=="probminhash_unweighted"){
  measured(c,[&]{return std::make_unique<Sketch::ProbMinHash4>(c.p,c.k,c.seed,0);},up,finish,[](auto& s){return s.getM()*8;},notempty,[](auto&){return NAN;},query,[](auto&,auto&){return NAN;});
 }else if(c.algorithm=="probminhash"){'''
    adapted = original.replace(needle, replacement)
    adapted = adapted.replace('../screening/real_input.h', 'real_input.h')
    # dr=0 is legal for KSSD; exact uses parameter=1.
    adapted = adapted.replace('c.p<1||c.passes<1', '(c.p<1&&c.algorithm!="kssd")||c.passes<1')
    (BUILD / 'real_adapter.cpp').write_text(adapted)
    (BUILD / 'real_input.h').write_bytes((HERE.parent / 'screening/real_input.h').read_bytes())
    (BUILD / 'CMakeLists.txt').write_text(f'''cmake_minimum_required(VERSION 3.16)
project(comprehensive_accuracy LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
set(CXXAPI ON CACHE BOOL "" FORCE)
set(RABBITSKETCH_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory("{SOURCE}" rabbitsketch EXCLUDE_FROM_ALL)
find_package(OpenMP REQUIRED)
find_library(DEFLATE_LIBRARY deflate REQUIRED)
add_executable(rabbitsketch_accuracy "{SOURCE}/examples/main.cpp" "{SOURCE}/examples/sub_command.cpp" "{HERE}/accuracy_adapter.cpp" real_adapter.cpp)
target_include_directories(rabbitsketch_accuracy PRIVATE "{SOURCE}/examples" "{BUILD}")
target_link_libraries(rabbitsketch_accuracy PRIVATE rabbitsketch_static OpenMP::OpenMP_CXX ${{DEFLATE_LIBRARY}})
''')
    subprocess.run(['cmake', '-S', str(BUILD), '-B', str(BUILD / 'cmake')], check=True)
    subprocess.run(['cmake', '--build', str(BUILD / 'cmake'), '-j', '16', '--target', 'rabbitsketch_accuracy'], check=True)
    paths = [HERE / 'accuracy_adapter.cpp', HERE / 'build.py', BUILD / 'real_adapter.cpp',
             HERE.parent / 'regimes/main_regimes.cpp', HERE.parent / 'screening/synthetic_pilot.cpp',
             BUILD / 'cmake/rabbitsketch_accuracy']
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (BUILD / 'build_manifest.json').write_text(json.dumps(hashes, indent=2) + '\n')


if __name__ == '__main__':
    main()
