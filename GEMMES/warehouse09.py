# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 20:24:05 2023

@author: Achilleas
"""


import os
import win32com.client as win32
import subprocess
import pandas as pd

import subprocess
import os
import win32com.client as win32
import pandas as pd
import logging

# Set the working directory
working_directory = os.getcwd()
os.chdir(working_directory)


# Configure logging
logging.basicConfig(filename='leap_gemmes.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define constants
AREA_NAME = "coupling_test_2"
SCENARIO_NAME = "Net zéro"
FAVORITE_CHART_NAME = "inv_costs"
#LEAP_CSV_PATH = r"C:\Users\Achilleas\Documents\Gemmes\foo111.csv"
#CLEANED_CSV_PATH = r"C:\Users\Achilleas\Documents\Gemmes\cleaned_foo2.csv"

leap = win32.Dispatch('LEAP.LEAPApplication')
#leap.Visible = True
#leap.ActiveArea = "coupling_test_2"
#leap.ActiveScenario = "Net zéro"
#leap.ActiveView = "Results"
#leap.ActivateScenario = "Reference"



def run_gemmes_init():
    logging.info("Running Gemmes initialization...")
    result = subprocess.call(["C:/Program Files/R/R-4.2.1/bin/Rscript", "C:/Users/Achilleas/Documents/Gemmes/RunGemmesInit.R"], shell=True)
    logging.info(f"Gemmes initialization result: {result}")

def run_gemmes():
    logging.info("Running Gemmes...")
    subprocess.call(["C:/Program Files/R/R-4.2.1/bin/Rscript", "C:/Users/Achilleas/Documents/Gemmes/RunGemmes.R"], shell=True)
    logging.info("Gemmes execution completed.")
    

def run_leap_init(area, scenario, favorite_name):
    logging.info("Running LEAP...")

    try:
        # Initialize LEAP
        leap = win32.Dispatch('LEAP.LEAPApplication')
        leap.Visible = True
        leap.ActiveArea = area
        leap.ActiveScenario = scenario
        leap.ActiveView = "Results"

        # Activate the specified favorite chart
        leap.Favorites(favorite_name).Activate()

        # Define the CSV file name
        csv_file_name = f"{favorite_name}_raw.csv"
        csv_file_path = os.path.join(r"C:\Users\Achilleas\Documents\Gemmes", csv_file_name)

        # Export data from the favorite chart to a CSV file
        leap.ExportresultsCSV(csv_file_path)
        logging.info("LEAP execution and data export completed.")
    except Exception as e:
        logging.error(f"An error occurred while running LEAP: {str(e)}")
        raise e        
        


def run_leap_analysis(area):
    logging.info("Running LEAP...")

    try:
        
        leap = win32.Dispatch('LEAP.LEAPApplication')
        leap.Visible = True
        leap.ActiveArea = area
        leap.ActiveView = "Analysis"

        logging.info("LEAP Analsyis update completed.")
    except Exception as e:
        logging.error(f"An error occurred while running LEAP: {str(e)}")
        raise e   


 


def clean_leap_output(leap_csv_path, cleaned_csv_path):
    logging.info("Cleaning LEAP output...")
    
    try:
        # Load the LEAP output CSV file
        df = pd.read_csv(leap_csv_path, encoding='ISO-8859-1', skiprows=4, header=0)


        # Replace "-" with 0 in the entire DataFrame
        df.replace('-', 0, inplace=True)

        # Remove the last column (axis=1) using iloc
        df = df.iloc[:, :-1]
        # Remove the last row
        df.drop(df.index[-1], inplace=True)
        
        df = df.transpose()

        # Set the first row as the header
        new_header = df.iloc[0]
        df = df[1:]
        df.columns = new_header
        
        # Use international language
        df.rename(columns={"Net zéro": "netzero"}, inplace=True)
        df.rename(columns={"Référence": "reference"}, inplace=True)
        
        # Export the cleaned data to a new CSV file
        df.to_csv(cleaned_csv_path,encoding='ISO-8859-1', index=False)
        logging.info("LEAP output cleaning completed.")
    except Exception as e:
        logging.error(f"An error occurred while cleaning LEAP output: {str(e)}")
        raise e
        



#Initialize the iteration
converged = False
iteration = 0
counter = 0
max_iterations = 3
inv_convergence = 0.001
prev_inv_sum = 0.1
# Start the iteration

run_gemmes_init()
run_leap_init("coupling_test",SCENARIO_NAME,"inv_costs")

while not converged and iteration < max_iterations:
    iteration += 1
    counter = iteration -1

    # Run Gemmes to update input data
    run_gemmes()

    run_leap_analysis("coupling_test")
    # Run LEAP
    run_leap_init("coupling_test",SCENARIO_NAME,"inv_costs")
    
    clean_leap_output(r"C:\Users\Achilleas\Documents\Gemmes\inv_costs_raw.csv", r"C:\Users\Achilleas\Documents\GitHub\Gemmes-Leap\GEMMES\inv_costs_clean.csv")



    # Check for convergence based on cleaned data
    # Load the cleaned data
    #cleaned_data = pd.read_csv(CLEANED_CSV_PATH)
    #df_conv = pd.read_csv(r"C:\Users\Achilleas\Documents\Gemmes\inv_costs_clean.csv")
    
    #if counter > 1:
    #    df_conv[f"netzero_{iteration}"] = df_conv["netzero"]
    #    df_conv[f"diff_{iteration}"] = df_conv[f"netzero_{iteration}"] - df_conv[f"netzero_{counter}"]
    # Calculate the sum of contents in the cleaned data
    #inv_sum = cleaned_data['column_name'].sum()  # Replace 'column_name' with the actual column name

    # Check for convergence
    #if abs(df_conv[f"diff_{iteration}"].sum() / df_conv[f"diff_{counter}"].sum() - 1) < inv_convergence:
    #    print(abs(df_conv[f"diff_{iteration}"].sum() / df_conv[f"diff_{counter}"].sum() - 1))
   #    converged = True
    #else:
    #    print(abs(df_conv[f"diff_{iteration}"].sum() / df_conv[f"diff_{counter}"].sum() - 1))
    #    converged = False

if converged:
    logging.info("Convergence achieved.")
else:
    logging.info("Maximum iterations reached without convergence.")