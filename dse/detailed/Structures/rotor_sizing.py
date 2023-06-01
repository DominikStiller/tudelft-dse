from material_properties import  materials
from StructureClasses import Beam, Force
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def y_transformation(theta, arr):
    if len(np.shape(arr)) == 2:
        av = np.mean(arr, 1)
        T = np.matrix([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
        rotated = np.asarray(T * arr)
        rotated += np.reshape(av-np.mean(rotated, 1), (2, 1))
        return rotated
    elif len(np.shape(arr)) == 3:
        dummy = np.zeros(np.shape(arr))
        for i in range(np.shape(arr)[0]):
            av = np.mean(arr[i], 1)
            T = np.matrix([[np.cos(theta[i]), np.sin(theta[i])], [-np.sin(theta[i]), np.cos(theta[i])]])
            dummy[i] = T * arr[i]
            dummy[i] += np.reshape(av - np.mean(dummy[i], 1), (2, 1))
        return np.asarray(dummy)
    else:
        raise TypeError('The array needs to be either 2D or 3D')


if __name__ == '__main__':
    ...
