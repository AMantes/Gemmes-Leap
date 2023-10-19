# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:01:50 2023

@author: Achilleas
"""

import os
import win32com.client as win32
import subprocess
import pandas as pd

# Set the working directory
working_directory = os.getcwd()
os.chdir(working_directory)

import subprocess
import os
import win32com.client as win32
import pandas as pd
import logging

# Configure logging
logging.basicConfig(filename='leap_gemmes.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define constants
AREA_NAME = "morocco_1"
SCENARIO_NAME = "Net zéro"
FAVORITE_CHART_NAME = "foo"
LEAP_CSV_PATH = r"C:\Users\Achilleas\Documents\Gemmes\foo111.csv"
CLEANED_CSV_PATH = r"C:\Users\Achilleas\Documents\Gemmes\cleaned_foo2.csv"

def run_gemmes_init():
    logging.info("Running Gemmes initialization...")
    result = subprocess.call(["C:/Program Files/R/R-4.2.1/bin/Rscript", "C:/Users/Achilleas/Documents/Gemmes/RunGemmesInit.R"], shell=True)
    logging.info(f"Gemmes initialization result: {result}")

def run_gemmes():
    logging.info("Running Gemmes...")
    subprocess.call(["C:/Program Files/R/R-4.2.1/bin/Rscript", "C:/Users/Achilleas/Documents/Gemmes/RunGemmes.R"], shell=True)
    logging.info("Gemmes execution completed.")

def run_leap(area, scenario, favorite_name):
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
        csv_file_name = f"{favorite_name}.csv"
        csv_file_path = os.path.join(working_directory, csv_file_name)

        # Export data from the favorite chart to a CSV file
        leap.ExportresultsCSV(csv_file_path)
        logging.info("LEAP execution and data export completed.")
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

        # Select the last row (Total)
        total = df[df["Branch"] == "Total"]

        # Export the cleaned data to a new CSV file
        total.to_csv(cleaned_csv_path, index=False)
        logging.info("LEAP output cleaning completed.")
    except Exception as e:
        logging.error(f"An error occurred while cleaning LEAP output: {str(e)}")
        raise e

if __name__ == "__main__":
    try:
        run_gemmes_init()
        run_gemmes()
        run_leap(AREA_NAME, SCENARIO_NAME, FAVORITE_CHART_NAME)
        clean_leap_output(LEAP_CSV_PATH, CLEANED_CSV_PATH)
    except Exception as e:
        logging.error(f"An error occurred during execution: {str(e)}")

# Close LEAP when done (optional)
# leap.Quit()



#Initialize the iteration
converged = False
iteration = 0
max_iterations = 100
inv_convergence = 0.0001
# Start the iteration
while not converged and iteration < max_iterations:
    # Run GEMMES
    
    get_and_parse_convergence_data()
    gemmes_results =run_gemmes()
    gemmes_prod = gemmes_results[["GDP", "PROD","YD"]]
    gemmes_prd.to_csv("gemmes_prod.csv",index=False) # export GDP, Prod, YD from Gemmes' boundary conditions
    
    # Load the updated GEMMES results into LEAP (ATM leap.LoadData is not defined)
    leap.LoadData("Gemmes Aggregate Demand, Production and Disposable Income", gemmes_prod.csv) # check the LEAP software documentation to determine how to load data; Eric suggested SEI to built-in this inside LEAP
    # run leap
    leap

    
    # Recalculate the sum of energy investment in LEAP and check for convergence
    prev_inv_sum = inv_sum
    inv_sum = 0
    for i in range(1, proc_branch.Children.Count):
        inv_sum += proc_branch.Children(i).Value(baseyear)
        
    prev_energy_price = energy_price
    energy_price = leap.Branch("...") # energy prices branch to be added here

    # Check for convergence; parameterise convergence criteria
    if abs(inv_sum/prev_inv_sum - 1)< inv_convergence and abs(energy_price/prev_energy_price - 1) < energy_p_convergence:
       converged = True
    
    # Increment the iteration counter
    iteration += 1
# End the iteration