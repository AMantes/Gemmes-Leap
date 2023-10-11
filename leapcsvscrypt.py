# -*- coding: utf-8 -*-
"""
Script that opens leap and extracts a favorite csv.
foo111 at the moment is solar panels investment.
It exports back a cleaned version that cen be imported in Gemmes.
The cleaned version needs to be transposed and cut to 32 years (atm 57); this takes place inside the Gemmes script.


@author: Achilleas
"""
import pandas as pd
import os
import win32com.client as win32

# Set the working directory
working_directory = os.getcwd()
os.chdir(working_directory)

# Initialize LEAP
leap = win32.Dispatch('LEAP.LEAPApplication')
leap.Visible = True
leap.ActiveArea = "lt-leds maroc v0.105_wresults_am"  # Set your appropriate area
leap.ActiveScenario = "Net zéro"  # Set your scenario
leap.ActiveView = "Results"

# Activate the "foo" favorite chart
favorite_name = "foo"
leap.Favorites(favorite_name).Activate()

# Define the CSV file name
csv_file_name = f"{favorite_name}111.csv"
csv_file_path = os.path.join(working_directory, csv_file_name)

# Export data from the favorite chart to a CSV file
leap.ExportresultsCSV(csv_file_path)

# Close LEAP
leap.Quit()

print(f"Data from '{favorite_name}' successfully exported to '{csv_file_path}'.")



# Load the "foo111.csv" as a DataFrame
file_path = r"C:\Users\Achilleas\Documents\Gemmes\foo111.csv"  # Replace with the actual path to your "foo.csv" file
df = pd.read_csv(file_path, encoding='ISO-8859-1', skiprows=4,header =0)

# Replace "-" with 0 in the entire DataFrame
df.replace('-', 0, inplace=True)


# Remove the last column (axis=1) using iloc
df = df.iloc[:, :-1]



# Select the last row (Total)
total = df[df["Branch"] == "Total"]
#transposed_total = total.transpose()
#transposed_total = transposed_total.reset_index(drop=True)
#transposed_total = transposed_total.iloc[1:]


# Create a new DataFrame with just the last row
cleaned_df = pd.DataFrame([total], columns=df.columns)

# Export the cleaned data to a new CSV file
cleaned_csv_file_path = r"C:\Users\Achilleas\Documents\Gemmes\cleaned_foo2.csv"  # Replace with your desired path
total.to_csv(cleaned_csv_file_path,index=False)

cleaned_df.to_csv(cleaned_csv_file_path, index=False)

print(f"Cleaned data saved to '{cleaned_csv_file_path}'.")
