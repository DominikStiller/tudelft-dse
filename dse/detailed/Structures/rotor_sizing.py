import numpy as np


def y_transformation(theta, arr, fix_points):
    T1 = np.array([[1, 0, fix_points[0][0]], [0, 1, fix_points[1][0]], [0, 0, 1]])
    T2 = np.array([[1, 0, -fix_points[0][0]], [0, 1, -fix_points[1][0]], [0, 0, 1]])
    if len(np.shape(arr)) == 2:
        a = np.vstack((arr, np.zeros(np.shape(arr)[1])))
        R = np.array([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        new = T1 @ R @ T2 @ a
        return new[:-1]
    elif len(np.shape(arr)) == 3:
        dummy = np.zeros((np.shape(arr)[0], 3, np.shape(arr)[2]))
        for i in range(np.shape(arr)[0]):
            a = np.vstack((arr[i], np.zeros(np.shape(arr[i])[1])))
            R = np.array([[np.cos(theta[i]), np.sin(theta[i]), 0], [-np.sin(theta[i]), np.cos(theta[i]), 0], [0, 0, 1]])
            dummy[i] = T1 @ R @ T2 @ a
        return dummy[:, :-1]
    else:
        raise TypeError("The array needs to be either 2D or 3D")


if __name__ == "__main__":
    ...
