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


class OnlineClimateFile(object):
    def __init__(self,filename,**kwargs):
        if type(filename) == str:
            self.filename = filename
        else:
            ValueError("\nExpected filename to be type str not {0}".format(type(filename)))
        
        df = pd.read_csv(self.filename,parse_dates=True)
        df
#         #Is the param file from the online interface or generated locally?
# 
#         #First line @s,parameters
#         if param_lines[0][0]=="@":
#             print "\nWeb interface generated file detected."
#             self.names  = get_interface_params(param_lines,f_len)
# 
#         #First line text from user.
#         else:
#             print "\nLocally generated file detected."
#             self.names  = get_local_params(param_lines,f_len)
# 
#         
# #             if i>0:
# #                 if param_lines[i-1][0:3]=="####":
# #                     print param_lines[i]
# #  
def print_missing_names(their_lst, our_lst):
    for their_name in their_lst:
        for our_name in our_lst:
            if their_name==our_name:
                break
        else:
            print their_name
if __name__=='__main__':
    my_file = "C:/Users/micah.johnson/Downloads/data-hru.csv"
    
    ui_file = OnlineClimateFile(my_file)
   
