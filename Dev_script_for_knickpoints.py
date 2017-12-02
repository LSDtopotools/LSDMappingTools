# Knickpoint plotting stuffs - B.G. if you need help
import matplotlib
matplotlib.use("Agg")
from LSDPlottingTools import LSDMap_KnickpointPlotting as KP


print("This is the example test script for the development version of the knickpoint plotting tools")

Smugglers = KP.KP_dev("/home/s1675537/PhD/LSDTopoData/knickpoint/Calibration/SC/Smugglers/", "smugglers", sources = [30,0,11,16,20,25,27,30,69], bd_KDE = 0.5)

Smugglers.plot_KDE_per_river()

print("done")