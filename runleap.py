# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 16:19:07 2023
This is to test github. the script opens LEAP, loads version 101, sets Zer Net as the active scenario and computes it. 
Then it closes Leap

@author: Achilleas
"""

import win32com.client as win32
import pandas as pd
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
    logging.info("LEAP opens successfully.")
except Exception as e:
    logging.error("Error occurred while opening LEAP: %s", str(e))

    
    leap.Visible = True
    
try:    
    leap.ActiveArea = "lt-leds maroc v0.101_wresults"  # Morocco v101 
    logging.info("Morocco area opens successfully.")
except Exception as e:
    logging.error("Error occurred while setting the active area in LEAP: %s", str(e))

        
   # baseyear = leap.BaseYear
    #logging.info("LEAP initialized successfully.")
    # Log leap setup information
    #logging.debug("LEAP setup: Visible=%s, ActiveArea=%s, BaseYear=%d", leap.Visible, leap.ActiveArea, baseyear)

#except Exception as e:
#    logging.error("Error occurred while initializing LEAP: %s", str(e))

#print(baseyear)


# =============================================================================
# 
# #### Load Scenarios
# 
# # Load the "Baseline" scenario
# #scenario_name = "Net zéro"  #"Current Accounts"  #"Référence"
# #scenario = leap.Scenarios(scenario_name).ResultsShown
# 
# #if scenario:
# #    logging.info("Scenario '%s' loaded successfully.", scenario_name)
# #else:
# #    logging.error("Scenario '%s' could not be loaded.", scenario_name)
# 
# 
# 
# #Runs everything
# # leap.Calculate()
# leap.ActiveScenario = "Net zéro"
# leap.ActiveView="Results"
# 
# 
# # Check if the active scenario was set successfully
# if scenario:
#     logging.info(f"Successfully set '{scenario_name}' as the active scenario.")
#     
#     # Run calculations for the active scenario
#     try:
#         scenario.Calculate()
#         logging.info(f"Scenario '{scenario_name}' calculations completed successfully.")
#     except Exception as e:
#         logging.error(f"Error occurred while running scenario '{scenario_name}' calculations: {str(e)}")
# else:
#     logging.error(f"Failed to set '{scenario_name}' as the active scenario.")
# 
# 
# # List of favorite charts
# favorite_charts = ['industry', 'residential', 'tertiaire', 'transport', 'transformation', 'resources', 'non_energy']
# 
# # Loop through the favorite views and export each one as a CSV
# for chart in favorite_charts:
#     try:
#         leap.Favorites(chart).Activate()
#         
#         # Set the file name for the CSV export
#         csv_file_name = f"{chart}.csv"
# 
#         # Concatenate the working directory and the file name to create the full path
#         csv_file_path = os.path.join(working_directory, csv_file_name)
# 
#         # Export the current favorite view as a CSV
#         leap.ExportresultsCSV(csv_file_path)
#         logging.info(f"Results for '{chart}' successfully exported to '{csv_file_path}'.")
#     except Exception as e:
#         logging.error(f"Error occurred while exporting results for '{chart}': {str(e)}")
# 
# 
# =============================================================================





# List of favorite charts
favorite_charts = ['industry', 'residential', 'tertiaire', 'transport', 'transformation', 'resources', 'non_energy']

# List of scenarios
scenarios = ['Référence', 'CDN Plus', 'Net zéro']

# Loop through the scenarios
for scenario in scenarios:
    try:
        # Set the active scenario
        leap.ActiveScenario = scenario
        leap.ActiveView = "Results"

        # Loop through the favorite charts and export each one as a CSV
        for chart in favorite_charts:
            leap.Favorites(chart).Activate()

            # Set the file name for the CSV export
            csv_file_name = f"{scenario}_{chart}.csv"

            # Concatenate the working directory and the file name to create the full path
            csv_file_path = os.path.join(working_directory, csv_file_name)

            # Export the current favorite chart as a CSV
            leap.ExportresultsCSV(csv_file_path)
            logging.info(f"Results for scenario '{scenario}', chart '{chart}' successfully exported to '{csv_file_path}'.")
    except Exception as e:
        logging.error(f"Error occurred while exporting results for scenario '{scenario}', chart '{chart}': {str(e)}")


# =============================================================================
# Close LEAP
try:
    leap.Quit()
    logging.info("LEAP closed successfully.")
except Exception as e:
    logging.error("Error occurred while closing LEAP: %s", str(e))