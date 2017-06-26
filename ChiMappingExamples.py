"""
This function drives the Chi mapping examples that are used with our documentation.

Written by Simon Mudd and Fiona Clubb
June 2017

GPL3
"""


import LSDChiMappingExamples as CME

if __name__ == "__main__":

    
    import time
    tic = time.clock()

    # Change these filenames and paths to suit your own files
    DataDirectory = "T:\\analysis_for_papers\\Xian\\"
    #DataDirectory = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\LSDTT_chi_examples\\'
    Base_file = 'Xian'


    # This is for the first example. Uncomment to get a hillshade image
    CME.ExampleOne_PartOne_SimpleHillshade(DataDirectory,Base_file)
    CME.ExampleOne_PartTwo_PrintBasins(DataDirectory,Base_file)
    
    toc = time.clock()
    print("This took: "+str(toc - tic)+" units of time")
    
