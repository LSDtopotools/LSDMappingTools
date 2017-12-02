# Knickpoint plotting stuffs - B.G. if you need help
import matplotlib
matplotlib.use("Agg")
from LSDPlottingTools import LSDMap_KnickpointPlotting as KP


print("This is the example test script for the development version of the knickpoint plotting tools")

Smugglers = KP.KP_dev("C:\Users\s1675537\Desktop\LSD\data\Smugglers\Smugglers/", "smugglers",binning = "source_key", sources = [30,0,11,16,20,25,27,30,69], bd_KDE = "scott", first_min_coeff = 0.000001)

Smugglers.plot_KDE( method = "deriv_delta_ksn", group = "source_key")
Smugglers.map_of_knickpoint( color_extent_Mx = [], method = "deriv_ksn")
Smugglers.plot_mchi_segments(method = "deriv_ksn", group = "source_key")
Smugglers.plot_chi_profiles( method = "deriv_ksn", group = "source_key")

print("done")