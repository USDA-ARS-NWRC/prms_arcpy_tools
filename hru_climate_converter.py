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
        
        df = pd.read_csv(self.filename,parse_dates=True)
        
        #Gather and edit the date to match the format needs for prms
    	for index, row in df.iterrows():
            orig_date = df.get_value(index,"date")
            month, day, year = orig_date.split("/")
            prms_date = " ".join([year,month,day,"0","0", "0"])
            df.set_value(index,"date",prms_date)
            
        self.output_data_file(df,"tmin","RCEW_Tmin.data")
            
        def output_data_file(self,df,col_str,filename):
            """
            Searches dataframe work for columns having the col_str,
            creates a new data frame work with the date append at the front
            of each line.
            Writes the file in space delimited format to filename
            """
            df2 = pd.DataFrame(df["date"])
            frames = []
            frames.append(df2)
            for name in list(df.columns.values):
                if col_str in name:
                    frames.append(df[name])
            result = pd.concat(frames)
            
            result.to_csv(filename," ")


if __name__=='__main__':
    my_file = "/home/micahjohnson/Documents/data-hru.csv"
    
    ui_file = OnlineClimateFile(my_file)
   
