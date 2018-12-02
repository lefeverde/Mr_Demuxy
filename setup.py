
import os
import site
import subprocess
import sys


from distutils.core import setup

# Functions to help w/ install #
def file_path_joiner(file_dir):
    files_list = os.listdir(file_dir)
    corrected_file_path = []
    for cur_file in files_list:
        cur_fp = os.path.join(file_dir, cur_file)
        corrected_file_path.append(cur_fp)
    return corrected_file_path

# def bash_checker(in_file, line_want):
#     append_file = True
#     for line in in_file:
#         if line_want in line:
#             append_file = False
#             break
#
#         if line_want not in line:
#             append_file = True
#     return append_file

# vars for making program work out of box #

script_files = file_path_joiner('bin')



setup(
    name = 'Mr_Demuxy',
    version = '1.2.2',
    description = 'demultiplexs combinatorially tagged reads',

    author = 'Daniel E. Lefever',
    author_email = 'lefeverde@gmail.com',

    #packages = ['Mr_Demuxy']
    packages = ['mr_demuxy'],
    package_dir = {
    'Mr_Demuxy': 'mr_demuxy',
    },
    scripts = script_files,

    package_data = {
    'Mr_Demuxy': ['example_data/example*'],
    'Mr_Demuxy': ['example_data/*.txt']
    },
    data_files = [
    ('Mr_Demuxy', ['MANIFEST.in']),
    ('', ['LICENSE/biopython_license.txt'])],
    url = 'https://pypi.python.org/pypi/Mr_Demuxy/',

    license = 'MIT',
    python_requires='>=2.7, >=3.6, <4',
    platforms = 'any',
    classifiers = [
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ]

    )
