
import numpy as np
import pandas as pd
import csv
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from datetime import date

###Open Files ###For testing script
#file = "E:\\temporarySkanitExport May Tue 31 2022 15-38-47-2687998.csv"
#d_csv = "DEMO.csv"
#s_csv = "DEMO_setup.csv"

#print(pd_data)
letters = ["A","B","C","D","E","F","G","H"]

###Find first row of plate(based on A)
def findA(df):
    ##determine number of rows and columns in df
    vorm = df.shape
    ##for columns in range( number of columns )
    for col in range(vorm[1]):
        ##for row number in range( length of column ): iteration over columns above, for this column length is determined en each item is checked for "A"
        for rownr in range(len(df[col])):
            if df[col][rownr] == 'A':
                ## return the coordinates of the letter "A"
                return (rownr, col)

###From raw data file to list with values   
def plate_to_list(file_rawdata):
    ##Open the raw data file, read as CSV, transform into dataframe
    file = open(file_rawdata)
    if file_rawdata.endswith(".csv"):
        csvdata = csv.reader(file, delimiter = ";")
    #elif file_rawdata.endswith(".xlsx") or file_rawdata.endswith(".xls"):
    #    csvdata = pd.read_excel(file)
    df = pd.DataFrame(csvdata)
    ## For this file, determine the position of "A". Assign Rownumber and Columnnumber to "r" and "c" respectively
    Aloc = findA(df)
    r = Aloc[0]
    c = Aloc[1]
    ##List with wells, list with the values. Iterate to fill both lists.
    well = []
    values = []
    for rows in range(8): ##A 96-well plate has 8 rows (A-H)
        for col in range(12): ## A 96-well plate has 12 columns (1-12)
            well_id = str(letters[rows]) + str(col +1) ## Combine both the letter and the column number to find well-ID on plate
            well.append(well_id) #Fill Well list with Well ID

            rplus = r + rows ##rows is number relative to startingvalue r (letter A)
            cplus = c + col + 1 ##col is number relative to startingvalue c (letter A), +1 Since a right shift is needed to jump from "A" to first column containg data
            val = df[cplus][rplus] ##find the value in the dataframe corresponding to the well-ID
            if val not in ("X", ""):
                val2 = val.replace(",",".") ##Change the decimal point: python uses ".", not ","
                value = float(val2) ##Make sure number is a float
            else:
                value = "" ##Empty value is "X"
            values.append(value) ##Fill value list with the found value
    return (well, values)

###Read Plate Setup/dilutions from CSV
def readfile(platesetup):
    ##Open plate, read as csv, store as dataframe
    plate = open(platesetup)
    csvplate = csv.reader(plate, delimiter = ";")
    df_plate = pd.DataFrame(csvplate)

    ##Make two lists to store the sampleID/dilution and the Well-ID
    samples = []
    well = []
    ##Iterate through 8 rows and 12 columns of plate
    for rows in range(8):
        for col in range(12):
            ##Make list with well_IDs
            well_id = str(letters[rows]) + str(col +1)
            well.append(well_id)
            
            ##List with value: location of well relative to first cell in csv; add 1 on row and column to start on first well of plate in file
            a = rows + 1
            b = col + 1
            naam = df_plate[b][a]
            samples.append(naam)
    return(well, samples)
  
###Compare and make df from both CSV files
def dataFramed(df, platesetup, dilutions):
    data = plate_to_list(df)
    names = readfile(platesetup)
    dilutions = readfile(dilutions)
    ##Make empty Lists
    wellid = []
    namesample = []
    dilute = []
    reads = []
    ##Put all data from 3 files into one dataframe
    ##Iterate for all 96 wells of the plate
    for i in range(96):
        wellid.append(data[0][i]) ##Every item of well-IDs
        namesample.append(names[1][i]) ##Every item in column with sampleID
        dilute.append(dilutions[1][i]) ##Every item in column with diltion factor
        reads.append(data[1][i]) ##Every item from RLU data column
    framed = pd.DataFrame(zip(wellid, namesample, dilute, reads), columns = ['well', 'sampleID', 'Dilution Factor', 'RLU value']) #Combine 4 lists to Dataframe
    return framed

###Function to normalize data
def normalize(value, min, max):
    x = 100* (1-(value - min)/(max-min))
    return float(x)

###Calculate average
def gemiddeld(lijst):
    som = 0
    for i in lijst:
        som += i
    return som/len(lijst)

###Function to calculate and plot curve/IC50/hill slope of the dataframe of the complete plate
def curvefit_calculate_IC50(df_plate, naming, l1_check):
    ##Selecting the control values
    global neutr_100
    neutr_100 = df_plate[(df_plate == 'ctrcell').any(axis=1)]['RLU value'].tolist()
    global minRLU
    minRLU = gemiddeld(neutr_100)
    global neutr_0
    neutr_0 = df_plate[(df_plate == 'ctrpp').any(axis=1)]['RLU value'].tolist()
    global maxRLU
    maxRLU = gemiddeld(neutr_0)
    global IC50_list
    global Hill_list
    IC50_list = []
    Hill_list = []    
    ##Calculations on each subset, including controls
    for samplename in naming:
        #Select data for one sampleID
        sub_data = df_plate[(df_plate.isin([samplename, "ctrpp", "ctrcell"])).any(axis=1)]
        pd.to_numeric(sub_data["Dilution Factor"])
        #Normalize every value by the controls
        normalized_RLU_value = []
        for RLU in sub_data['RLU value']:
            normalized_RLU_value.append(normalize(RLU, minRLU, maxRLU))
        #Make sure data type of dilutions is integer
        dilution_integer = []
        for dil in sub_data["Dilution Factor"]:
            dilution_integer.append(int(dil))
        #Log transform xvalues
        dilution_log = np.log10(dilution_integer)
        #Make new datastructure containing only log(dil) and normalized RLU (x and y values)
        sub_data_normalized = pd.DataFrame(list(zip(dilution_log, normalized_RLU_value)), columns = ["x_val", "y_val"])
        sub_data_sorted = sub_data_normalized.sort_values("x_val")
        #Make fit and print IC50
        curvefitted(sub_data_sorted, samplename, l1_check)
    plt.show()
    return IC50_list, Hill_list

###Make curvefit of data
def curvefitted(data_to_fit, name_of_sample, l1_check):
    ##Define the curve to fit
    def logistic_curve(logx, logIC50, Hill):
        return 100 / (1 + 10**((logIC50-np.asarray(logx))*Hill))
    ###Determine distance of real value to expected value
    def residuals(param, xs, ys):
        return (100 / (1 + 10**((param[0]-np.asarray(xs))*param[1])))-ys    
    ##Make fit on data, parameter will contain (LogIC50, Hill)
    if l1_check.get() == 0:
        parameter, covariance = curve_fit(logistic_curve, data_to_fit.x_val, data_to_fit.y_val, maxfev=5000)
    elif l1_check.get() == 1:
        soft_l1 = least_squares(residuals, [3,-1], loss='soft_l1', f_scale=15, args=(data_to_fit.x_val, data_to_fit.y_val))    
        parameter = [soft_l1.x[0], soft_l1.x[1]]
    else:
        print("SOMETHING'S WRONG")   
    ##Print IC50 and Hill slope
    IC50 = 10**(parameter[0])    ##Change log(IC50) to IC50
    IC50_list.append(IC50)
    Hill_list.append(parameter[1])
    print("The IC50 value for ",name_of_sample, "is %.0f" %IC50)
    print("The Hill Slope value for ",name_of_sample, "is %.2f" %parameter[1])
    #create y-values based on x_val and parameters (this is the fitted curve)
    fit1 = logistic_curve(data_to_fit["x_val"], parameter[0], parameter[1])
    #with same parameters create a smooth line curve
    line_range=np.arange(0, 8, 0.1)
    line = logistic_curve(line_range, parameter[0], parameter[1])
    ##Make a plot
    plt.ylim(0, 100)
    plt.xlabel("Log(Dilution Factor)")
    plt.ylabel("% Neutralization")
    plt.xlim(1,6.5)
    plt.plot(line_range, line)
    plt.scatter(data_to_fit["x_val"], data_to_fit["y_val"])

###SAVE DATA TO MASTER FILE
def master_write(file, names, IC50s, Hills, experi, l1):
    with open(file, 'a', newline='') as file_open:
        write_to_file = csv.writer(file_open, delimiter=";")
        for i in range(len(names)):
            IDsample = names[i]
            IC50 = IC50s[i]
            hill = Hills[i]
            fields = [date.today(), IDsample, experi, str(IC50).replace(".", ","), str(hill).replace(".", ","), l1]
            write_to_file.writerow(fields)
        file_open.close()
        print("DONE!")

#print(curvefit_calculate_IC50(dataFramed(file, s_csv, d_csv), ["s1", "s2"], 0)) ##TEST SCRIPT WITH THIS