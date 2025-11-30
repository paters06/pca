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

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

def read_input_files(file_str: str):
    array_list = []
    with open(file_str) as file:
        while True:
            line = file.readline()
            if not line:
                break
            array_list.append(line.split())

    num_rows = len(array_list)
    num_cols = len(array_list[0])
    imported_mat = np.zeros((num_rows,num_cols))

    for i in range(0,num_rows):
        imported_mat[i,:] = np.array(array_list[i], dtype=float)
    
    return imported_mat

def plotTemperatureField(temp_mat: np.ndarray) -> None:
    cx = temp_mat[:,0]
    cy = temp_mat[:,1]
    T = temp_mat[:,3]

    fig, ax = plt.subplots()

    field = ax.tripcolor(cx,cy,T,vmin=T.min(),vmax=T.max())
    ax.set_title('Temperature Field')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    # ax1.set_aspect('equal')
    cb1 = fig.colorbar(field,label='[°C]')

    # plt.tight_layout()
    plt.show()


def main() -> None:
    file_name = sys.argv[1]
    output_name = 'output_' + file_name
    temp_mat = read_input_files(output_name)
    plotTemperatureField(temp_mat)


if __name__ == "__main__":
    main()