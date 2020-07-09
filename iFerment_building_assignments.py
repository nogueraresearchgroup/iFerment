#Building assignments	Abel Trinidad Ingle

#Date: June 25, 2020

#This file is intended to be used in conjunction with the iFerment python model and iFerment_assignments.
#The execution of this file is intended for terminal with the above-mentioned documents, which are to be located
#in local machine directories.

import cobra
import logging
import pandas as pd
import copy

#Below, I am using a SBML formatted version of the model, which is generated when the iFerment model is ran (Line 5235), 
#effectively copying the model as "iFerment".

iFerment = cobra.io.read_sbml_model("iFerment186_Plus.xml")

#A copy of iFerment is created for downstream purposes
iFerment2 = copy.deepcopy(iFerment)

#microbe_select is a string of corresponding to the headers in the iFerment_assignments file
microbe_select = ["Melsdenii", "Pacidpropionici", "Wbifida", "Ckluyveri", "Lacidophilus", "Dsuccinatiphilus", "Lbrevis"]

#NOTE: Again, the below path is specific to my personal Mac
Assignment = pd.read_excel("iFerment_assignments.xlsx")

#Counter below
n = 0

#Some of the lines in this for-loop are present to observe script runs in Jupyter Notebook

#For each microbe in the microbe_select string. . . 
for microbe in microbe_select:
	#. . . for each reaction
    for index, row in Assignment.iterrows():
    	#. . . if there's a 1
        if row[microbe_select[n]] == 1:
        	#Effectively, nothing happens. The master iFerment still possesses the reaction of interest.
        	#This line is intended for Jupyter Notebook.
            print(iFerment.reactions.get_by_id(row['id']).upper_bound)
        
        #And if there's not a 1, if there's a 0 for that reaction. . . 
        else:
            print("I don't have this reaction" ,iFerment.reactions.get_by_id(row['id']))
            #Then knock the reaction out. The master iFerment has been changed to mimmic the organism in microbe_select
            iFerment.reactions.get_by_id(row['id']).knock_out()
            
        #If iterated reaction is currently the last reaction. . . 

#This is for Jupyter Notebook
print("You're on the last one")
#Then perform a pfba using the iFerment model given the knock outs, and. . . 
pfba_solution = cobra.flux_analysis.pfba(iFerment)
#Generate an excel sheet with pfba solutions.
writer = pd.ExcelWriter("iFermentAs" + microbe_select[n] + ".xlsx")
pfba_solution.fluxes.to_excel(writer,'Sheet1')
writer.save()
            
#The current version of iFerment possesses knockouts, so. . . 
            
#Revert it to its master version with no knockouts.
iFerment = copy.deepcopy(iFerment2)
            
#Next microbe in microbe_select.
n=n+1