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
                month, day, year = orig_date.split("/")

            except:
                if index == 0:
                    print "\nError: Unexpected line before data in csv file, try deleting the row with value types.\n CSVFILE line {0}".format(index+2)
                else:
                    print "\nError: Unknown issue with the date column on line {0} in csv file".format(index+2)
                raise SystemExit
            
            prms_date = ' '.join([year,month,day,'0','0', '0'])
            self.df.set_value(index,"date",prms_date)

    def output_data_file(self,col_str,filename):
        """
        Searches dataframe work for columns having the col_str,
        creates a new data frame work with the date append at the front
        of each line.
        Writes the file in space delimited format to filename
        """
        df2 = pd.DataFrame(self. df["date"])
        frames = []
        frames.append(df2)
        for name in list(self.df.columns.values):
            header_name = name.split("[")
            if col_str in header_name[0]:
                #If we have a name match then add the name plus the associated hru i.e. tmin[23]
                df2[name] = self.df[name]
        
        vals_len = len(list(df2.columns.values))-1
        #Provide feedback is usr provided string was not found.
        if vals_len <= 1 :
            print "\nWarning: Column name {0} was not found in the file.".format(col_str)
        
        else:
            print "\nWriting {0} columns that contained the string {1}". format(vals_len,col_str)
                    
        with open(filename,'w') as f:
            #Append the PRMS expected Header
            f.write("File Generated using Micahs hru_climate_converter.py\n")
            f.write("{0} {1}\n".format(col_str,vals_len))
            
            f.write("########################################\n")
            df2.to_csv(f,"\t", header = False, index = False)
            f.close()
        print "\tData file outputted to {0}".format(filename)

if __name__=='__main__':
    my_file = "/home/micahjohnson/projects/eWFS/RCEW/data/climate/data-hru.csv"
    ui_file = OnlineClimateFile(my_file)

    ui_file.output_data_file("tmin","/home/micahjohnson/projects/RCEW_Model/prms_input/RCEW_Tmin.data")
    ui_file.output_data_file("tmax","/home/micahjohnson/projects/RCEW_Model/prms_input/RCEW_Tmax.data")
    ui_file.output_data_file("precip","/home/micahjohnson/projects/RCEW_Model/prms_input/RCEW_precip.data")
    ui_file.output_data_file("swe","/home/micahjohnson/projects/RCEW_Model/prms_input/RCEW_swe.data")
    

   
