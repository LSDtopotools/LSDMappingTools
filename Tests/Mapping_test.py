"""
This is a tester function for some of the mapping scripts

Author: SMM

Date 05/01/2018
"""

import sys
import os

# We need to import the path since we don't install mapping tools directly
this_dir = os.getcwd()
print("This dir is: "+this_dir)

LSDMT_dir = this_dir.split(os.sep)

LSDMT_dir = LSDMT_dir[:-1]
LSDMT_path = (os.sep).join(LSDMT_dir)

print("The path is: "+LSDMT_path)
sys.path.append(LSDMT_path)

# import Maping tools
import LSDMapWrappers as LSDMW

# Plot a simple hillshade
DataDir = this_dir+os.sep
LSDMW.SimpleHillshade(DataDir,"WA")