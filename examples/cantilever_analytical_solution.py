import numpy as np
import matplotlib.pyplot as plt


def create_UV_evaluation_points(x_value:float, y_range:list[float], num_points:int):
    # UV_eval_pts = np.zeros((num_points, 2))
    start_pt = np.array((x_value, y_range[0]))
    end_pt = np.array((x_value, y_range[1]))
    points = np.linspace(start_pt, end_pt, num_points)
    print(points)


def create_evaluation_points(x_value, y_range, num_points: int):
    eval_pts = np.zeros((num_points, 2))
    y_list = np.linspace(y_range[0], y_range[1], num_points)
    for i in range(num_points):
        eval_pts[i, 0] = x_value
        eval_pts[i, 1] = y_list[i]
    
    return eval_pts


def evaluate_fields(eval_pts, D: float, L: float, P: float):
    num_points = eval_pts.shape[0]
    eval_fields = np.zeros((num_points, 2))
    for i in range(num_points):
        y = eval_pts[i, 0]
        x = eval_pts[i, 1]
        # print("x: {} m, y: {} m".format(x, y))
        ux, uy = displacements(x, y, D, L, P)
        sx, sy, txy = stresses(x, y, D, L, P)
        eval_fields[i, 0] = uy
        eval_fields[i, 1] = sx

    return eval_fields    


def displacements(x: float, y: float, D: float, L: float, P: float, E:float = 210e9, nu:float = 0.3):
    I = (D**3)/12
    u1 = -(P/(6*E*I))*(y - (D/2))*((3*x*(2*L - x)) + (2 + nu)*y*(y-D))
    u2 = (P/(6*E*I))*((x**2)*(3*L - x) + 3*nu*(L-x)*(y-(D/2))**2 + ((4+5*nu)/4)*(D**2)*x)

    return u1, u2

def stresses(x: float, y: float, D: float, L: float, P: float):
    I = (D**3)/12
    sigma_11 = -(P/I)*(L-x)*(y-(D/2))
    sigma_22 = 0.0
    sigma_33 = -((P*y)/(2*I))*(y-D)

    return sigma_11, sigma_22, sigma_33


def test1():
    L = 1.0
    D = 0.2
    x1 = L
    y1 = 0.5*D
    x2 = 0.0
    y2 = D
    P = -1
    E = 2e5
    nu = 0.31
    ux, uy = displacements(x1, y1, D, L, P)
    sx, sy, txy = stresses(x2, y2, D, L, P)

    print("Ux: {} m, Uy: {} m".format(ux, uy))
    print("Sx: {} Pa".format(sx))

def test2():
    L = 1.0
    D = 0.2
    x_value = 1.0
    y_range = [0.0, 1.0]
    num_points = 20
    y_list = np.linspace(y_range[0], y_range[1], num_points)
    P = -1
    E = 2e5
    nu = 0.31
    eval_pts = create_evaluation_points(x_value, y_range, num_points)
    UV_eval_pts = create_UV_evaluation_points(x_value, y_range, num_points)
    eval_field = evaluate_fields(eval_pts, D, L, P)

    # print(eval_field)
    plt.plot(y_list, eval_field[:, 1], '.', label='Analytical solution. Sx')
    plt.legend()
    # plt.show()

if __name__ == "__main__":
    # test1()
    test2()