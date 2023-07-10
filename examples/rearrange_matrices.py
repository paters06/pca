import numpy as np

def find_match_rows(A_mat: np.ndarray, B_mat: np.ndarray) -> np.ndarray:
    """
    A is the sample matrix
    B is the global matrix

    This functions return the indices of the rows of the B matrix
    that corresponds to each row of the A matrix
    """
    match_indices_list = []

    for i in range(0, A_mat.shape[0]):
        for j in range(0, B_mat.shape[0]):
            equal_row = np.where(A_mat[i,:] == B_mat[j,:])[0]
            if len(equal_row) == 2:
                match_indices_list.append(j)
    
    match_indices = np.array(match_indices_list, dtype=int)
    return match_indices

def main():
    A_mat = np.array([[0.,0.],
                      [0.2,0.],
                      [0.4,0.],
                      [0.,0.1],
                      [0.2,0.1],
                      [0.4,0.1],
                      [0.,0.2],
                      [0.2,0.2],
                      [0.4,0.2]])

    B_mat = np.array([[0.4,0.],
                      [0.7,0.],
                      [1.,0.],
                      [0.4,0.1],
                      [0.7,0.1],
                      [1.,0.1],
                      [0.4,0.2],
                      [0.7,0.2],
                      [1.,0.2]])

    C_mat = np.array([[0.,0.],
                      [0.,0.1],
                      [0.,0.2],
                      [0.2,0.],
                      [0.2,0.1],
                      [0.2,0.2],
                      [0.4,0.],
                      [0.4,0.1],
                      [0.4,0.2],
                      [0.7,0.],
                      [0.7,0.1],
                      [0.7,0.2],
                      [1.,0.],
                      [1.,0.1],
                      [1.,0.2]])
    
    A_indices = find_match_rows(A_mat, C_mat)
    print('A matrix')
    print(A_indices)

    B_indices = find_match_rows(B_mat, C_mat)
    print('B matrix')
    print(B_indices)
    

if __name__ == '__main__':
    main()
