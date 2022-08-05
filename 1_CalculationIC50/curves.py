import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
from sklearn.metrics import r2_score
from re import S
from pandas.core.frame import DataFrame

####DONE####
#IC50, Normalized data, min max not included, absolute
def IC50_Norm_Abs(xval, datas, min, max):
    
    for dats in datas:
        #Average and Normalize
        data = 100*(1-(dats.mean(axis=1) - min) / ( max - min ))

        #Curve Function
        def logistic_curve(x, logIC50, Hill):
            #return 100 / (1 + 10**((logIC50-np.asarray(x))*Hill))
            return ((-100) / (1 + ((np.asarray(x)/logIC50)**Hill))) + 100

        #Curve fit to data, Print parameters
        para, cov = curve_fit(logistic_curve, xval, data, maxfev=5000)
        print("The value of IC50 is %.0f and the value of the Hill Slope is %.2f." %(10**para[0], para[1]))

        #Create the model curve with these parameters, Compare with data
        fit1 = logistic_curve(xval, para[0], para[1])

        score = r2_score(data, fit1)
        print("The R² score (for average values) is %.3f" %(score))

        lijstdat = dats.stack().tolist()
        lijstdats = [100*(1-(m - min) / ( max - min )) for m in lijstdat]

        fitss = []
        for i in range(len(fit1)):
            extra = [fit1.tolist()[i]] * 3
            fitss.extend(extra)
        
        score2 = r2_score(lijstdats, fitss)
        print("The R² score (for all values) is %.3f \n" %(score2))

        #Plot Generation
        plt.plot(xval, fit1)
        plt.ylim(0, 100)
        plt.xlabel("Log(Dilution Factor)")
        plt.ylabel("% Neutralization")
        plt.xlim(1.5, 5.1)
        plt.scatter(xval, data)
        plt.scatter(para[0], 50, c='red', marker="X")
    plt.show()

####DONE####
#IC50, Normalized data, min max not included, relative
def IC50_Norm_Rel(xval, datas):
    
    for dats in datas:
        #Choose min and max repsonse
        low = min(dats.min())
        high = max(dats.max())

        #Averaging and Normalization
        data = 100*(1-(dats.mean(axis=1) - low) / ( high - low ))

        #Curve
        def logistic_curve(x, logIC50, Hill):
            return 100 / (1 + 10**((logIC50-np.asarray(x))*Hill))

        #Fit to curve
        para, cov = curve_fit(logistic_curve, xval, data, maxfev=5000)
        print("The value of IC50 is %.0f and the value of the Hill Slope is %.2f." %(10**para[0], para[1]))

        fit1 = logistic_curve(xval, para[0], para[1])

        score = r2_score(data, fit1)
        print("The R² score (for average values) is %.3f" %(score))

        #Calculate r² for all data points (not only averages)
        lijstdat = dats.stack().tolist()
        lijstdats = [100*(1-(m - low) / ( high - low )) for m in lijstdat]

        fitss = []
        for i in range(len(fit1)):
            extra = [fit1.tolist()[i]] * 3
            fitss.extend(extra)
        
        score2 = r2_score(lijstdats, fitss)
        print("The R² score (for all values) is %.3f \n" %(score2))

        plt.plot(xval, fit1)
        plt.ylim(0, 100)
        plt.xlabel("Log(Dilution Factor)")
        plt.ylabel("% Neutralization")
        plt.xlim(1.5, 5.1)
        plt.scatter(xval, data)
        plt.scatter(para[0], 50, c='red', marker="X")
    plt.show()

####DONE####
#IC50, Normalized data, min/max included in sample data, absolute
def IC50_Norm_Abs_plateaus(xval, datas, min, max):
    #Add min and max x times to dataset to set plateaus
    rep = 6
    
    for i in range(rep):
        xval.insert(0, 0)
        xval.append(8)

    #Average and Normalize
    for dats in datas:
        data = 100*(1-(dats.mean(axis=1) - min) / ( max - min ))

    #Add Min and Max
        for i in range(rep):
            data.loc[-1] = 100
            data.index = data.index + 1
            data.sort_index(inplace=True) 
            data[len(data)+i] = 0

        def logistic_curve(x, logIC50, Hill):
            return 100 / (1 + 10**((logIC50-np.asarray(x))*Hill))

        para, cov = curve_fit(logistic_curve, xval, data, maxfev=8000)
        print("The value of IC50 is %.0f and the value of the Hill Slope is %.2f." %(10**para[0], para[1]))

        fit1 = logistic_curve(xval, para[0], para[1])

        score = r2_score(data, fit1)
        print("The R² score (for average values) is %.3f" %(score))

        lijstdat = dats.stack().tolist()
        lijstdats = [100*(1-(m - min) / ( max - min )) for m in lijstdat]

        fitss = []
        for i in range(len(fit1)):
            extra = [fit1.tolist()[i]] * 3
            fitss.extend(extra)
        
        #score2 = r2_score(lijstdats, fitss)
        #print("The R² score (for all values) is %.3f \n" %(score2))

        plt.plot(xval, fit1)
        plt.ylim(0, 100)
        plt.xlabel("Log(Dilution Factor)")
        plt.ylabel("% Neutralization")
        plt.xlim(1.5, 5.1)
        plt.scatter(xval, data)
        plt.scatter(para[0], 50, c='red', marker="X")
    plt.show()

####TO ADJUST####
#IC50, Normalized data, min max repeat, relative
def IC50_Norm_Rel_plateaus(xval, datas):
    for dats in datas:
        #Choose min and max repsonse
        low = min(dats.min())
        high = max(dats.max())

        rep = 6
    
        for i in range(rep):
            xval.insert(0, 0)
            xval.append(8)
        #Averaging and Normalization
        data = 100*(1-(dats.mean(axis=1) - low) / ( high - low ))

        for i in range(rep):
            data.loc[-1] = 100
            data.index = data.index + 1
            data.sort_index(inplace=True) 
            data[len(data)+i] = 0
        
        #Curve
        def logistic_curve(x, logIC50, Hill):
            return 100 / (1 + 10**((logIC50-np.asarray(x))*Hill))

        #Fit to curve
        para, cov = curve_fit(logistic_curve, xval, data, maxfev=5000)
        print("The value of IC50 is %.0f and the value of the Hill Slope is %.2f." %(10**para[0], para[1]))

        fit1 = logistic_curve(xval, para[0], para[1])

        score = r2_score(data, fit1)
        print("The R² score (for average values) is %.3f" %(score))

        #Calculate r² for all data points (not only averages)
        lijstdat = dats.stack().tolist()
        lijstdats = [100*(1-(m - low) / ( high - low )) for m in lijstdat]

        fitss = []
        for i in range(len(fit1)):
            extra = [fit1.tolist()[i]] * 3
            fitss.extend(extra)
        
        #score2 = r2_score(lijstdats, fitss)
        #print("The R² score (for all values) is %.3f \n" %(score2))

        plt.plot(xval, fit1)
        plt.ylim(0, 100)
        plt.xlabel("Log(Dilution Factor)")
        plt.ylabel("% Neutralization")
        plt.xlim(1.5, 5.1)
        plt.scatter(xval, data)
        plt.scatter(para[0], 50, c='red', marker="X")
    plt.show()

def IC50_allpoints(xval, datas, min, max):
    
        for dats in datas:
            #Keep all data, Normalize
            data = 100*(1-(dats - min) / ( max - min ))
                        
            #Change data format from df to list
            ListRLU = data.stack().tolist()
            ListData = []
            for i in range(len(xval)):
                three = [xval[i]] * 3
                ListData.extend(three)

            AllData = pd.DataFrame(np.column_stack([ListData, ListRLU]), columns = ["Log(dil)", "RLU"])
            
            def logistic_curve(x, logIC50, Hill):
                return 100 / (1 + 10**((logIC50-np.asarray(x))*Hill))

            para, cov = curve_fit(logistic_curve, ListData, ListRLU, maxfev=5000)
            print("The value of IC50 is %.0f and the value of the Hill Slope is %.2f." %(10**para[0], para[1]))

            fit1 = logistic_curve(ListData, para[0], para[1])

            Comp = pd.DataFrame(np.column_stack([AllData["RLU"].tolist(), fit1.tolist()]), columns = ["True", "Predict"])

            score2 = r2_score(Comp["True"], Comp["Predict"])
            print("The R² score (for all values) is %.3f \n" %(score2))

            plt.plot(ListData, fit1)
            plt.ylim(0, 100)
            plt.xlabel("Log(Dilution Factor)")
            plt.ylabel("% Neutralization")
            plt.xlim(1.5, 5.1)
            plt.scatter(ListData, ListRLU)
            plt.scatter(para[0], 50, c='red', marker="X")
        plt.show()

