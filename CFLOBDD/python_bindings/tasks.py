import invoke
import pathlib
import sys
import os
import shutil
import re
import glob

def print_banner(msg):
    print("==================================================")
    print("= {} ".format(msg))

@invoke.task()
def build_cflobdd(c):
    """Build the shared library for the sample C++ code"""
    print_banner("Building C++ Library")
    invoke.run(
        "g++ -O3 -Wall -Werror -shared -std=c++17 -w -fPIC -I../ -I../../../boost_1_80_0/ -I../Solver/uwr/bit_vector -I../Solver/uwr/matrix -I../Solver/uwr/parsing -lm ../*.cpp ../Solver/uwr/bit_vector/*.cpp ../Solver/uwr/parsing/*.cpp "
        "-o libcflobdd.so "
    )
    print("* Complete")

def compile_python_module(cpp_name, extension_name):
    invoke.run(
        "g++ -O3 -Wall -Werror -shared -std=c++17 -fPIC -I../../../boost_1_80_0/ "
        "`python3.6 -m pybind11 --includes` "
        "-I /usr/include/python3.6 -I .  "
        "{0} "
        "-o {1}`python3.6-config --extension-suffix` "
        "-L. -lcflobdd -Wl,-rpath,.".format(cpp_name, extension_name)
    )

#@invoke.task(build_cflobdd)
@invoke.task()
def build_pybind11(c):
    """Build the pybind11 wrapper library"""
    print_banner("Building PyBind11 Module")
    compile_python_module("python_wrapper.cpp", "python_wrapper_cflobdd")
    print("* Complete")

@invoke.task()
def test_pybind11(c):
    """Run the script to test PyBind11"""
    print_banner("Testing PyBind11 Module")
    invoke.run("python3 pybind11_test.py", pty=True)
