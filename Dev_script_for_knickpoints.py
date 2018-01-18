# Knickpoint plotting stuffs - B.G. if you need help
import matplotlib
matplotlib.use("Agg")
from LSDPlottingTools import LSDMap_KnickpointPlotting as KP
import pandas as pd
from matplotlib import pyplot as plt


# # print("This is the example test script for the development version of the knickpoint plotting tools")

# # Smugglers = KP.KP_dev("/home/s1675537/PhD/LSDTopoData/knickpoint/model/movern_0p65_n_is_one94/", "movern_0p65_n_is_one94",binning = "source_key",min_chi_dist_for_kp = 0.0001, sources = [], bd_KDE = "scott", first_min_coeff =0.001, max_coeff = 0.5) #[30,0,11,16,20,25,27,30,69]

# # Smugglers.plot_KDE( method = "deriv_ksn", group = "basin_key")
# # Smugglers.map_of_knickpoint(method = "deriv_ksn", color_extent_Mx = [],color_extent_kp = [], cut_off = 0)
# # Smugglers.plot_mchi_segments(method = "deriv_ksn", group = "basin_key")
# # Smugglers.plot_chi_profiles( method = "deriv_ksn", group = "basin_key")
# # Smugglers.plot_KDE( method = "deriv_ksn", group = "source_key")
# # Smugglers.map_of_knickpoint(method = "deriv_ksn", color_extent_Mx = [],color_extent_kp = [], cut_off = 0)
# # Smugglers.plot_mchi_segments(method = "deriv_ksn", group = "source_key")
# # Smugglers.plot_chi_profiles( method = "deriv_ksn", group = "source_key")

# # print("done")

# df = pd.read_csv("/home/s1675537/PhD/LSDTopoData/knickpoint/model/movern_0p65_n_is_one94/movern_0p65_n_is_one94_MChiSegmented_Ks.csv")
# df = df[df["chi"]>=0]
# for source in df["source_key"].unique():
# 	print("source %s out of %s"%(source,df["source_key"].unique().shape[0] ))
# 	this_df = df[df["source_key"] == source]
# 	plt.plot(this_df["chi"], this_df["elevation"], c = "b")
# 	plt.plot(this_df["chi"], this_df["segmented_elevation"], c = "r")
# 	plt.savefig("/home/s1675537/PhD/LSDTopoData/knickpoint/model/movern_0p65_n_is_one94/source_%s.png"%source)
# 	plt.clf()


path = "/home/s1675537/PhD/LSDTopoData/knickpoint/puerto_rico/"
prefix = "puerto_rico"

Jaguar = KP.KP_dev_v2(path,prefix, min_length = 20000)

# Jaguar.DEBUG_print_ksn_filters()
# Jaguar.DEBUG_print_KDE()
# Jaguar.DEBUG_print_ksn_outliers()
# Jaguar.DEBUG_print_ksn_dksndchi()

Jaguar.DEBUG_print_ksn_comparison()