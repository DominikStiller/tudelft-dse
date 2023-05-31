import pandas as pd
import numpy as np
import csv


def read_material_properties(filename):
    material_properties = pd.read_csv(filename, delimiter=";", index_col=0)
    return material_properties


a = read_material_properties("materials.csv")
print(read_material_properties("materials.csv"))
