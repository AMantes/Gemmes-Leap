# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 11:43:49 2023

@author: Achilleas Mantes
"""

import os
import logging
import win32com.client as win32

# Set the working directory
working_directory = r"C:\Users\Achilleas\Documents\Gemmes"
os.chdir(working_directory)

# Set up logging
logging.basicConfig(
    filename='script.log',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.DEBUG
)

def open_leap():
    try:
        leap = win32.Dispatch('LEAP.LEAPApplication')
        leap.Visible = True
        leap.ActiveArea = "lt-leds maroc v0.101_wresults"  # Morocco v101 
        logging.info("LEAP initialized successfully.")
        return leap
    except Exception as e:
        logging.error("Error occurred while initializing LEAP: %s", str(e))
        return None

def export_favorite_chart(leap, chart, scenario, working_directory):
    try:
        leap.Favorites(chart).Activate()
        csv_file_name = f"{scenario}_{chart}.csv"
        csv_file_path = os.path.join(working_directory, csv_file_name)
        leap.ExportresultsCSV(csv_file_path)
        logging.info(f"Results for scenario '{scenario}', chart '{chart}' successfully exported to '{csv_file_path}'.")
    except Exception as e:
        logging.error(f"Error occurred while exporting results for scenario '{scenario}', chart '{chart}': {str(e)}")

def close_leap(leap):
    try:
        leap.Quit()
        logging.info("LEAP closed successfully.")
    except Exception as e:
        logging.error("Error occurred while closing LEAP: %s", str(e))

def main():
    # Open LEAP
    leap = open_leap()
    
    if leap:
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
                    export_favorite_chart(leap, chart, scenario, working_directory)
            except Exception as e:
                logging.error(f"Error occurred while processing scenario '{scenario}': {str(e)}")
        
        # Close LEAP
        close_leap(leap)

if __name__ == "__main__":
    main()