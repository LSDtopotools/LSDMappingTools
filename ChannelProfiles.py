# quick script for plotting channel profiles

import LSDPlottingTools as LSDP

DataDirectory = '/home/fionaclubb/muddpile/movern_0p5/n_is_one/'
fname = 'movern_0p5_n_is_one200'

LSDP.LSDMap_ChiPlotting.ChannelProfilePlot(DataDirectory,fname,FigFormat='png')
