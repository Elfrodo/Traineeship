### Load Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
from sklearn.metrics import r2_score
from re import S
from pandas.core.frame import DataFrame
from scipy import stats


###read CSV file
csvdata = pd.read_csv("../Z_Data/CSVread3.csv", sep = ";", decimal=",")

###Set log(DilutionF) column as index and remove 'min' and 'max' text so only the real x-values remain
logspre = csvdata["log(dil)"].values.tolist()
logspre.remove("min")
logspre.remove("max")
logs = []
for i in logspre:
    k = i.replace(",",".")
    k = float(k)
    logs.append(k)

alldata = csvdata.set_index("log(dil)")
alldata.index.names = [None]

###Make Subsets of data for each sample
##generate list of names of samples:
rep = int(len(alldata.columns)/3)
dataNames = []
for i in range(rep):
    naam = "sample_" + str(i+1)
    dataNames.append(naam)

##Add Data of sample (value) to name of sample (key) in dictionary
tel = 0
frames = []
for k in range(rep):
    frame = alldata.iloc[:,[tel, tel+1, tel+2]]
    tel = tel+3
    frames.append(frame)

joined = dict(zip(dataNames, frames))

##IC50 curve, initial
def IC50(logs, dict):
    identity = []
    IC50 = []
    rsq = []
    hill = []
    for names, frame in dict.items():

        identity.append(names)

        mins = frame.loc["min"][0]
        maxs = frame.loc["max"][0]
        dataOnly = frame.drop(["min","max"])

        xvalues = []
        for i in range(len(logs)):
            extra = [logs[i]] * 3
            xvalues.extend(extra)
        
        lijstdat = dataOnly.stack().tolist()
        lijstdata = [100*(1-(m - mins) / ( maxs - mins )) for m in lijstdat]
        
        ###correct values outside of 0-100 range to fall in this range:
        #for i in range(len(lijstdata)):
        #    if lijstdata[i]<0:
        #        lijstdata[i] = 0
        #    elif lijstdata[i]>100:
        #        lijstdata[i] = 100

        def logistic_curve(logx, logIC50, Hill):
            return 100 / (1 + 10**((logIC50-np.asarray(logx))*Hill))
            ###Or other formula:
            # return ((-100) / (1 + ((np.asarray(logx)/logIC50)**Hill))) + 100

        line = np.arange(0, 8, 0.1)
        para, cov = curve_fit(logistic_curve, xvalues, lijstdata, maxfev=5000)


        ###restrict Hill Slope to be between two values: (Not correct way to do so, should be set before fitting, not after)
        #if para[1]>-0.8:
        #    para[1]=-0.8
        #elif para[1]<-1.2:
        #    para[1]=-1.2

        fit1 = logistic_curve(xvalues, para[0], para[1])
        fit2 = logistic_curve(line, para[0], para[1])

        ##Fill lists with results:
        IC50.append(int(10**para[0]))
        hill.append("%.2f" %para[1])
        score = r2_score(lijstdata, fit1)
        rsq.append("%.2f" %score)

        plt.plot(line, fit2)
        plt.ylim(0, 100)
        plt.xlabel("Log(Dilution Factor)")
        plt.ylabel("% Neutralization")
        plt.xlim(1.5, 5.1)
        # plot datapoints
        #plt.scatter(xvalues, lijstdata)

        plt.scatter(para[0], 50, c='red', marker="X")
    results_dict = {"Sample":identity, "IC50":IC50, "HillSlope":hill, "r²": rsq}
    df_results = pd.DataFrame(results_dict)
    print(df_results)
    plt.show()

IC50(logs,joined)