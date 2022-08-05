import tkinter as tk
from tkinter import Toplevel, ttk
from tkinter.filedialog import askopenfilename
import mymodule as mm
import pandas as pd
import numpy as np

###Set up Window "root" and set title and icon. Define canvas size
root = tk.Tk()
root.title("Calculation of IC50 Value")
#root.iconbitmap('etn.ico')
canvas = tk.Canvas(root, width = 1100, height = 500)
canvas.grid(columnspan=6, rowspan=9)
###Set an image
img = tk.PhotoImage(file="testfile.gif")
canvas.create_image(0, 0, anchor='nw',  image=img)

##instructions TEXT
info = tk.Label(root, text="Tool written by Elfriede. Picture too!")
info.grid(column=0, row=9)
instruct = tk.Label(root, text="Select CSV file to calculate IC50 of samples.")
instruct.grid(column=2, row=0)
samp_format = tk.Label(root, text="Choose the number of dilutions and replications /sample")
samp_format.grid(column=2, row=1)
dil_serie = tk.Label(root, text="Choose start dilution and dilution factor")
dil_serie.grid(column=2, row=2)
check_l1loss = tk.Label(root, text="Check to reduce influence of outliers")
check_l1loss.grid(column=2, row=3)
names_fill = tk.Label(root, text="Fill in names of samples on plate.\nEach name on new line, no need for comma's")
names_fill.grid(column=2, row=4)
Project = tk.Label(root, text="Choose Project!")
Project.grid(column=2, row=5)
Experiment = tk.Label(root, text="Fill in Experiment code (INF0xx/NTR0xx)")
Experiment.grid(column=2, row=6)

##Look for file FUNCTION
def open_data():
    browse.set("loading...")
    file_path = askopenfilename(title="Choose csv file(s)", filetypes=[('CSV files', '*.csv'),("All files","*.*"), ("Excel files",'*.xlsx')])
    global file_loc
    file_loc = file_path
    if file_path:
        browse.set("Browse")
        ##Display last 25 characters of choosen filename
        fileChoosen = tk.Label(root, text= ("..." + file_path[-25: -1: 1] + "v"))
        fileChoosen.grid(column=4, row=0, columnspan=2)

##Browse Button; set text to BROWSE, use open_data function to retrieve file location when clicked
browse = tk.StringVar()
browse.set("Browse")
browse_b = ttk.Button(root, textvariable = browse, command = open_data)
browse_b.grid(column = 3, row = 0)

##Checkbox, choose to use L1-or not
l1_checked = tk.IntVar()
l1regul = ttk.Checkbutton(root, text= "Robust Regression (L1-loss)", variable=l1_checked)
l1regul.grid(column = 3, row = 3)

##Dropdown menu with choices for plate setup/ dropdown menu dilutions and choise of project
plate_Setup_list = {"Choose Sample Format":"NA", "8 dilutions, 3 replicates":"platesetups/PlateSetup8x3_ctr12.csv","DEMO": "platesetups/DEMO_setup.csv"}
platesetup = tk.StringVar(root)
platesetup.set(plate_Setup_list)
dropd_Setup = ttk.OptionMenu(root, platesetup, *plate_Setup_list)
dropd_Setup.grid(column=3, row=1)
dropd_Setup.configure(width = 35)

plate_Dilutions_list = {"Choose dilution start and factor":"NA","Start 100, dilute 3 (8x3)":"dilutions/dil100_3x.csv", "DEMO": "dilutions/DEMO.csv"}
dilutions = tk.StringVar(root)
dilutions.set(plate_Dilutions_list)
dropd_Dilutions = ttk.OptionMenu(root, dilutions, *plate_Dilutions_list)
dropd_Dilutions.grid(column=3, row=2)
dropd_Dilutions.configure(width = 35)

project_list = {"Choose Project":"NA","DEMO":"MasterFiles/MASTER_demo.csv", "PseudoParticles":"MasterFiles/MASTER_pseudoparticles.csv"}
project = tk.StringVar(root)
project.set(project_list)
drop_project = ttk.OptionMenu(root, project, *project_list)
drop_project.grid(column=3, row=5)
drop_project.configure(width = 35)

##Fill in names of samples, each line one name
txt_names = tk.Text(width=20, height=6, font =("Helvetica", 10))
txt_names.grid(column = 3, row = 4)
names = txt_names.get("1.0", tk.END)
print(names)

##Fill in experiment number
txt=tk.StringVar()
experiment = tk.Entry(root, textvariable=txt)
experiment.grid(column = 3, row = 6)

###Make calculate button. run function from mymodule to generate one dataframe
##The functions
def process():
    #run the dataframed function on the choosen file, store in global variable so other functions can use output
    global df_dataCombined
    df_dataCombined = mm.dataFramed(file_loc, plate_Setup_list[platesetup.get()], plate_Dilutions_list[dilutions.get()])
    #store seperate sample names in file to global variable "samples"
    global samples
    samples = pd.unique(df_dataCombined[~df_dataCombined["sampleID"].isin(["ctrcell", "ctrpp", ""])]["sampleID"])
    
def sample_naming():
    ##Retreive filled in names, store in global variable, split names by 'enter'
    global names
    names = txt_names.get("1.0", tk.END+"-1c").split("\n")
    ##Make dictionary containing names from csv file and given names, these will be assigned to each other
    toekenning = dict(zip(samples, names))
    print(toekenning)
    ##Replace original s1, s2, s3 by the choosen names, skip original 'ctrcell' and 'ctrpp' as to not use these as 'samples'.
    global df_renamed
    df_renamed = df_dataCombined.replace({"sampleID": toekenning})

#The button to check the data by showing the plots and values, before saving the results to the master file   
button = ttk.Button(root, text='Process data', command = lambda:(process(), sample_naming(), check_graphs()))#
button.grid(column= 4, row= 8)

##Generate Popup window to check graphs before saving IC50 values to master file
def check_graphs():
    popup= Toplevel()
    IC50s, Hills = mm.curvefit_calculate_IC50(df_renamed, names, l1_checked)
    text_results = ""
    for i in range(len(names)):
        result_zin = ("The IC50 value of " + names[i] + " is " + str(np.round(IC50s)[i]) + ". Hill slope is " + str(np.round(Hills,2)[i]) + ". \n")
        text_results = (text_results+ result_zin)
    label_graphs= tk.Label(popup, text=text_results)
    label_graphs.grid(column = 1, row = 0)

    ##Get status l1-checkbox  
    global l1_status
    if l1_checked.get() == 1:
        l1_status = "L1 applied"
    else:
        l1_status = "Not applied"
    ##write input experiment nr to variable
    experiment_id = experiment.get()

    ##Make button to save results to masterfile (and report)
    #popup= Toplevel()
    save_status = tk.StringVar()
    status_text = "Save to Master File"
    save_status.set(status_text)   
    save_files = ttk.Button(popup, text=status_text, command = lambda:(mm.master_write(project_list[project.get()], names, IC50s, Hills, experiment_id, l1_status)))
    save_files.grid(column= 0, row= 2)

    exit = ttk.Button(popup, text = "Exit",command = popup.destroy)
    exit.grid(column = 2, row = 2) 

root.mainloop()