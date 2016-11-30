"""
This script is for converting files that were generated online
via the ewsf modeling website http://ewsf.engr.colostate.edu/

The website is currently accesible by invite only and offers a really nice
way to generate climate by hru data files for use with PRMS.

The online platform takes data in a different format so this script converts the
files for PRMS that is developed by the USGS and used locally.
"""

from fileinput import filename
import __main__
import pandas as pd
import numpy as np
from __main__ import __name__

class OnlineClimateFile(object):
    def __init__(self,filename,**kwargs):
        if type(filename) == str:
            self.filename = filename
        else:
            ValueError("\nExpected filename to be type str not {0}".format(type(filename)))
        
        try:
            self.df = pd.read_csv(self.filename,parse_dates=True)
        except:
            print "\nError: File Failed to have expected data."
            print "\nTry the following: "
            print " * Delete the header.\n * Delete excess date columns. \n * Make sure the only date column is the first column."
            raise SystemExit
        
        #Gather and edit the date to match the format needs for prms
    	for index, row in self.df.iterrows():

            orig_date = self.df.get_value(index,"date")

            try:
                if "/" in orig_date:
                    month, day, year = orig_date.split("/")
                    
                elif "-" in orig_date:
                    year, month, day = orig_date.split("-")                    
            
            except:
                
                if index == 0:
                    print "\nError: Unexpected line before data in csv file, try deleting the row with value types.\n CSVFILE line {0}".format(index+2)
                else:
                    print "\nError: Unknown issue with the date column on line {0} in csv file".format(index+2)
                raise SystemExit
            
            prms_date = {"year":year,"month":month,"day":day,"hour":"0","minute":"0", "second":"0"}
        
            for time_name, time_value in prms_date.items():
                self.df.set_value(index, time_name,time_value)

    def output_data_file(self,col_str,filename):
        """
        Searches dataframe work for columns having the col_str,
        creates a new data frame work with the date append at the front
        of each line.
        Writes the file in space delimited format to filename
        """
        frames = []
        df2 = pd.DataFrame(self.df[["year","month","day","hour","minute","second"]])         
        frames.append(df2)
        
        for name in list(self.df.columns.values):
            header_name = name.split("[")
            if col_str in header_name[0]:
                #If we have a name match then add the name plus the associated hru i.e. tmin[23]
                df2[name] = self.df[name]
        
        vals_len = len(list(df2.columns.values))-6
        #Provide feedback is usr provided string was not found.
        if vals_len <= 1 :
            print "\nWarning: Column name {0} was not found in the file.".format(col_str)
        
        else:
            print "\nWriting {0} columns that contained the string {1}". format(vals_len, col_str)
                    
        with open(filename,'w') as f:
            #Append the PRMS expected Header
            f.write("File Generated using Micahs hru_climate_converter.py\n")
            f.write("{0} {1}\n".format(col_str,vals_len))
            
            f.write("########################################\n")
            df2.to_csv(f," ", header = False, index = False)
            f.close()
        print "\tData file outputted to {0}".format(filename)

    def write_climate_data(self,prms_input_dir):
        """
        Writes all the data required for PRMS Climate by HRU
        data to their respective files in the appropriate format
        for local prms runs.
        """
        
        climate = ["tmin","tmax", "precip", "swe"]
        for data in climate:
            self.output_data_file(data,prms_input_dir + data + ".data")
   
    def write_runoff_data(self,prms_input_dir):
        """
        eWSF allows the user to collact station data to be used for prms with 
        ease. This converts the file to prms local executable to read.
        
        args:
            prms_input_dir    This is the location of the directory where PRMS will looks for data
        """
        
        data = "runoff"
        self.output_data_file(data,prms_input_dir + data + ".data")


if __name__=='__main__':

#     my_file = "/home/micahjohnson/projects/eWFS/RCEW/data/climate/data-hru.csv"
#     prms_input_dir = "/home/micahjohnson/projects/RCEW_Model/prms_input/"
    
    my_file = "/home/scotthavens/Documents/Projects/PRMS/BRB/data/stationData/runoff_data-sta_final.csv"
    prms_input_dir = "/home/scotthavens/Documents/Projects/PRMS/BRB/prms_input"
        
    ui_file = OnlineClimateFile(my_file)   
    #ui_file.write_climate_data(prms_input_dir)
    ui_file.write_runoff_data(prms_input_dir)


   
