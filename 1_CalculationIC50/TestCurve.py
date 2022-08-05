###SCRIPT TO TEST SOME CURVES TO FIT
##In inhibition curve the controls (these determine the plateau) can either be or not be included in the curve. 

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
from sklearn.metrics import r2_score
from re import S
from pandas.core.frame import DataFrame
import curves

##hard coded controls
min = 52
max = 334826

csvdata = pd.read_csv("../Z_Data/data_testcurve.csv", sep = ";", decimal=",")

s1 = csvdata[["Sample 1a", "Sample 1b", "Sample 1c"]]
s2 = csvdata[["Sample 2a", "Sample 2b", "Sample 2c"]]
s3 = csvdata[["Sample 3a", "Sample 3b", "Sample 3c"]]


logDil = []
for i in csvdata["dilfactor"]:
    k = math.log(int(i), 10)
    logDil.append(k)

##change active function below to try other calculation method:
#curves.IC50_Norm_Abs(logDil, [s1, s2, s3], min, max) ##Absolute IC50, ctr not included in dataset
#curves.IC50_Norm_Rel(logDil, [s1, s2, s3]) ##Relative IC50, ctr not included in dataset
curves.IC50_Norm_Abs_plateaus(logDil, [s1, s2, s3], min, max) ##Absolute IC50, ctr included in dataset
#curves.IC50_Norm_Rel_plateaus(logDil, [s1, s2, s3]) ##Relative IC50, ctr included in dataset
#curves.IC50_allpoints(logDil, [s1, s2, s3], min, max) ##uses all points instead of average only for eacht x-value