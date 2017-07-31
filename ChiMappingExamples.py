"""
This function drives the Chi mapping examples that are used with our documentation.
The website is: https://lsdtopotools.github.io/LSDTopoTools_ChiMudd2014/

Written by Simon Mudd and Fiona Clubb
June 2017
git
GPL3
"""


import LSDChiMappingExamples as CME
import LSDMapWrappers as MW

if __name__ == "__main__":

    """
    This is just a few lines for keeping track of how long the program is taking.
    You can ignore it.
    """
    import time
    tic = time.clock()

    """
    These lines tell the example functions where your files are. If you are using
    the recommended file setup you won't need to modify these lines.

    If you have set up your own directory structure you will need to modify
    the directory names.
    """
    #DataDirectory = "S:\\movern_analysis\hokkaido\\"
    #Base_file = "hokkaido_points"
    #DataDirectory = "T:\\analysis_for_papers\\Xian\\"
    #DataDirectory = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\LSDTT_chi_examples\\'
    #Base_file = 'Xian2'
    DataDirectory = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Meghalaya\\Paper\\'
    Base_file = 'Meghalaya'

    """
    These lines are used to run the examples. Each line calls a function in the directory
    ../LSDChiMappingExamples/

    The individual examples are different python scripts which you can inspect at your
    lesiure. By playing with these scripts you can learn how to use our plotting tools.

    To run the examples simply comment or uncomment these lines by adding or
    removing the comment symbol, #, below.
    """
    
    MW.SimpleHillshade(DataDirectory,Base_file,cmap = "terrain", dpi = 250)
    MW.PrintBasins(DataDirectory,Base_file)
    #MW.PrintBasins_Complex(DataDirectory,Base_file,Remove_Basins = [6], Rename_Basins = {3:"A",4:"B",8:"C",10:"D",12:"E"})
    #MW.PrintChannels(DataDirectory,Base_file)
    MW.PrintChannelsAndBasins(DataDirectory,Base_file, dpi = 300)
    
    #CME.ExampleOne_PartOne_SimpleHillshade(DataDirectory,Base_file)
    #CME.ExampleOne_PartTwo_PrintBasins(DataDirectory,Base_file)
    #CME.ExampleOne_PartThree_PrintBasinsWithLabels(DataDirectory,Base_file)
    #CME.ExampleOne_PartFour_MaskBasins(DataDirectory,Base_file)
    #CME.ExampleOne_PartFive_MaskBasinsMF(DataDirectory,Base_file)
    
    toc = time.clock() 
    print("This took: "+str(toc - tic)+" units of time")
