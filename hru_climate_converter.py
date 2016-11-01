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
        
        self.df = pd.read_csv(self.filename,parse_dates=True)
        
        #Gather and edit the date to match the format needs for prms
    	for index, row in self.df.iterrows():
            orig_date = self.df.get_value(index,"date")
            month, day, year = orig_date.split("/")
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
                
        #Provide feedback is usr provided string was not found.
        if len(list(df2.columns.values)) <= 1 :
            print "\nWarning: Column name {0} was not found in the file.".format(col_str)
        
        else:
            print "\nWriting {0} columns that contained the string {1}". format(len(list(df2.columns.values))-1,col_str)
                    
        with open(filename,'w') as f:
            #Append the PRMS expected Header
            f.write("########################################\n")
            df2.to_csv(f,"\t", header = False, index = False)
            f.close()
        print "\tData file outputted to {0}".format(filename)

if __name__=='__main__':
    my_file = "/home/micahjohnson/Documents/data-hru.csv"
    
    ui_file = OnlineClimateFile(my_file)
    
    ui_file.output_data_file("tmin","/home/micahjohnson/Documents/RCEW_Tmin.data")
    ui_file.output_data_file("tmax","/home/micahjohnson/Documents/RCEW_Tmax.data")
    ui_file.output_data_file("precip","/home/micahjohnson/Documents/RCEW_precip.data")
    ui_file.output_data_file("swe","/home/micahjohnson/Documents/RCEW_swe.data")
    

   
