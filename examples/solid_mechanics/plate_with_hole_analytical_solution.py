import numpy as np
import matplotlib.pyplot as plt

def create_UV_evaluation_points(x_range:list[float], y_range:list[float], num_points:int):
    start_pt = np.array((x_range[0], x_range[1]))
    end_pt = np.array((y_range[0], y_range[1]))
    points = np.linspace(start_pt, end_pt, num_points)
    return points

def create_evaluation_points(x_range:list[float], y_range:list[float], num_points:int):
    start_pt = np.array((x_range[0], x_range[1]))
    end_pt = np.array((y_range[0], y_range[1]))
    points = np.linspace(start_pt, end_pt, num_points)
    return points

def stresses(pt: np.ndarray, a: float, sigma_o: float) -> np.ndarray:
    """
    Params
    ------
    a is the radius of the hole
    sigma_o is the applied tensions

    Notes
    ------

    sigma[0,0] -> s_xx
    sigma[0,1] -> s_yy
    sigma[0,2] -> t_xy
    """
    sigma = np.zeros((pt.shape[0],3))

    radius = np.sqrt(pt[:,0]**2 + pt[:,1]**2)
    theta = np.arctan2(pt[:,1],pt[:,0])

    sigma[:,0] = sigma_o*(1.0 - ((a/radius)**2)*(1.5*np.cos(2*theta) + np.cos(4*theta)) + (1.5*(a/radius)**4)*np.cos(4*theta))
    sigma[:,1] = sigma_o*(-((a/radius)**2)*(0.5*np.cos(2*theta) - np.cos(4*theta)) - (1.5*(a/radius)**4)*np.cos(4*theta))
    sigma[:,2] = sigma_o*(-((a/radius)**2)*(0.5*np.sin(2*theta) + np.sin(4*theta)) + (1.5*(a/radius)**4)*np.sin(4*theta))
    return sigma

def test_01():
    L = 4.
    R = 1.
    nu = 0.3
    E = 1000.
    sigma_o = 1.0

    x_range = [0.0, R]
    y_range = [0.0, L]
    num_points = 20
    path_pt = create_evaluation_points(x_range, y_range, num_points)
    eval_field = stresses(path_pt, R, sigma_o)

    # print(eval_field)
    plt.plot(path_pt[:,1], eval_field[:,0], label='Analytical solution. Sx')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    test_01()
