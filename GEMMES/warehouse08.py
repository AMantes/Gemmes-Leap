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

# Usage
area_name = "morocco_1" # set the appropriate LEAP model
scenario_name = "Net zéro" #set the scenarios you want to run
favorite_chart_name = "foo" # set the favorite chart you want to export
leap_csv_path = r"C:\Users\Achilleas\Documents\Gemmes\foo111.csv" # set the path and csv name for LEAP's output
cleaned_csv_path = r"C:\Users\Achilleas\Documents\Gemmes\cleaned_foo2.csv" # set the path and name for the cleaned csv; to be imported in Gemmes


# Initialise Gemmes
# Gemmes is compiled and runs with initial values.
# RunGemmesInit exports an annualised GDP growth rate in "combined_vector.csv" that feeds directly in LEAP    
def run_gemmes_init():
    resG= subprocess.call(["C:/Program Files/R/R-4.2.1/bin/Rscript","C:/Users/Achilleas/Documents/Gemmes/RunGemmesInit.R"], shell=True)  # change paths as needed; it runs from powershell 
    resG 
    
# Run Gemmes
# To be used after Gemmes is compiled
def run_gemmes():
    subprocess.call(["C:/Program Files/R/R-4.2.1/bin/Rscript", "C:/Users/Achilleas/Documents/Gemmes/RunGemmes.R"], shell=True)

# Run LEAP
def run_leap(area, scenario, favorite_name):
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
    
    # Close LEAP when done (optional)
    #leap.Quit()


def clean_leap_output(leap_csv_path, cleaned_csv_path):
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
















# Call the functions
run_leap(area_name, scenario_name, favorite_chart_name)
clean_leap_output(leap_csv_path, cleaned_csv_path)