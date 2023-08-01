# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 16:19:07 2023
This is to test github. the script opens LEAP, loads version 101, sets Zer Net as the active scenario and computes it. 
Then it closes Leap

@author: Achilleas
"""

import win32com.client as win32
import subprocess
import logging
import csv
import os


# Set the working directory
working_directory = r"C:\Users\Achilleas\Documents\Gemmes"
os.chdir(working_directory)


# Set up logging
logging.basicConfig(
    filename='script.log',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.DEBUG
)



# Verify the current working directory
current_directory = os.getcwd()
logging.info("Current working directory: %s", current_directory)



# Set up LEAP
try:
    leap = win32.Dispatch('LEAP.LEAPApplication')
    leap.Visible = True
    leap.ActiveArea = "lt-leds maroc v0.101_wresults"  # Morocco v101 
    baseyear = leap.BaseYear
    logging.info("LEAP initialized successfully.")
    # Log leap setup information
    logging.debug("LEAP setup: Visible=%s, ActiveArea=%s, BaseYear=%d", leap.Visible, leap.ActiveArea, baseyear)

except Exception as e:
    logging.error("Error occurred while initializing LEAP: %s", str(e))

#print(baseyear)



#### Load Scenarios

# Load the "Baseline" scenario
scenario_name = "Net zéro"  #"Current Accounts"  #"Référence"
scenario = leap.Scenarios(scenario_name).ResultsShown

if scenario:
    logging.info("Scenario '%s' loaded successfully.", scenario_name)
else:
    logging.error("Scenario '%s' could not be loaded.", scenario_name)



#Runs everything
# leap.Calculate()
leap.ActiveScenario = "Net zéro"
leap.ActiveView="Results"
# Set the active scenario
scenario_name = "Net zéro"
scenario = leap.Scenarios(scenario_name)

# Check if the active scenario was set successfully
if scenario:
    logging.info(f"Successfully set '{scenario_name}' as the active scenario.")
    
    # Run calculations for the active scenario
    try:
        scenario.Calculate()
        logging.info(f"Scenario '{scenario_name}' calculations completed successfully.")
    except Exception as e:
        logging.error(f"Error occurred while running scenario '{scenario_name}' calculations: {str(e)}")
else:
    logging.error(f"Failed to set '{scenario_name}' as the active scenario.")


# Close LEAP
try:
    leap.Quit()
    logging.info("LEAP closed successfully.")
except Exception as e:
    logging.error("Error occurred while closing LEAP: %s", str(e))