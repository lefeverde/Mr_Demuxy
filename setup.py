
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

def bash_checker(in_file, line_want):
    append_file = True
    for line in in_file:
        if line_want in line:
            append_file = False
            break
            
        if line_want not in line:
            append_file = True
    return append_file

# vars for making program work out of box #

script_files = file_path_joiner('bin')
user_home = os.path.expanduser('~')
user_path = site.USER_BASE
export_line = 'export PATH=' + user_path + '/bin:$PATH'
bash_file = os.path.join(user_home, '.bash_profile')


setup(
    name = 'Mr_Demuxy',
    version = '1.2.0',
    description = 'demultiplexs combinatorially tagged reads',

    author = 'Daniel E. Lefever',
    author_email = 'lefeverde@gmail.com',
    
    packages = ['Mr_Demuxy'],
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

    platforms = 'any',
    classifiers = [
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ]
    
    )

sys.stdout.write('\nensuring local dir is in path\n')
sys.stdout.flush() 
try:
    sys.stdout.write('bash_profile found\nadding local dir to path\n')
    sys.stdout.flush()
    
    pe_alias = 'alias pe_demuxer'
    merg_alias = 'alias merged_demuxer'
 
    open_bf = open(bash_file, 'rU')
    bash_lines = open_bf.readlines()
    open_bf.close()

    add_local_path = bash_checker(bash_lines, export_line)

    with open(bash_file, 'w+') as f:
        if add_local_path is True:
            f.write(
                '\n#Added by Mr. Demuxy\n'
                + export_line + '\n\n'
                )
        for line in bash_lines:
            if pe_alias not in line:
                if merg_alias not in line:
                    f.write(line)
                
            
        

except IOError:
    sys.stdout.write('creating new bash_profile\nadding local dir to path\n')
    sys.stdout.flush()
    with open(bash_file, 'w+') as f:
        f.write(
            '\n#Added by Mr. Demuxy\n'
             + export_line + '\n'
             )      
