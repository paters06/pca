#######################################################################
# DO NOT REMOVE THIS SEGMENT
import os
import sys

# Get current working directory
# Command found in https://note.nkmk.me/en/python-script-file-path/
dir1 = os.path.dirname(os.path.abspath(__file__))
# Insert .. command to go to the upper directory
dir2 = dir1 + '/..'
# Setting the package directory path for the modules execution
sys.path.append(dir1)

os.chdir(dir1)
#######################################################################

def read_input_files(file_str: str):
    with open(file_str) as file:
        while True:
            line = file.readline()
            if not line:
                break
            print(line.strip())

if __name__ == '__main__':
    read_input_files('input_file_cantilever.txt')
