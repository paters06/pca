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


def plot_geometry(contour, P_mat, dirichlet_points) -> None:
    """
    cbpts: contour boundary points
    P_mat: control points
    dirichlet_points: control points where Dirichlet conditions are enforced
    load_coor: points with Neumann boundary conditions
    """
    fig = plt.figure()
    ax = plt.axes()
    plt.axis('equal')
    ax.use_sticky_edges = False
    titlestring = "Geometry with boundary conditions"
    ax.axis("off")

    # Boundary Geometry
    fieldplot = ax.fill(contour[:,0],contour[:,1],facecolor='none',edgecolor='black',linewidth=1.5)

    # Control Points
    controlplot = ax.scatter(P_mat[:,0],P_mat[:,1])

    dirichletplot = plt.scatter(dirichlet_points[:,0],dirichlet_points[:,1],c = "g",marker = 'o')

    first_string = "Control points"
    second_string = "Temperature BC"
    third_string = "Flux BC"

    # if len(loadcoor) != 0:
    #     # neumannplot = ax.scatter(loadcoor[:,0],loadcoor[:,1],c = "b",marker = "s")
    #     neumannplot = ax.quiver(loadcoor[:,0],loadcoor[:,1],loadfield[:,0],loadfield[:,1],color=['b'])
    #     # Uncomment for example2
    #     plt.legend((dirichletplot,neumannplot),(first_string,second_string),loc='upper right',bbox_to_anchor=(1.2,1.0))
    #     # Uncomment for example3
    #     # plt.legend((dirichletplot,neumannplot),('Displacement restrictions','Load conditions'),loc='lower right',bbox_to_anchor=(1.2,0.0))
    # else:
    #     # Uncomment for example2
    #     plt.legend([dirichletplot],[first_string],loc='upper right',bbox_to_anchor=(1.2,1.0))
    #     # Uncomment for example3
    #     # plt.legend((dirichletplot),('Displacement restrictions'),loc='lower right',bbox_to_anchor=(1.2,0.0))

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(titlestring)

    plt.legend([controlplot, dirichletplot],[first_string, second_string],loc='upper right',bbox_to_anchor=(1.2,1.0))
    plt.tight_layout()
    plt.show()


def main() -> None:
    file_name = sys.argv[1]
    contour_name = 'boundary_points_' + file_name
    P_name = 'control_points_' + file_name
    dirichlet_name = 'dirichlet_points_' + file_name

    contour_mat = read_input_files(contour_name)
    P_mat = read_input_files(P_name)
    dirichlet_points = read_input_files(dirichlet_name)
    plot_geometry(contour_mat, P_mat, dirichlet_points)


if __name__ == "__main__":
    main()