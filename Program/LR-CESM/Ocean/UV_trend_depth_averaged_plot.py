#Program plots the local velocity trend 

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap

#Making pathway to folder with all data
directory	= '../../../Data/LR-CESM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min	= 0
depth_max	= 100

year_start	= 2000
year_end	= 2100

month_start	= 1 	#1 = January, 2 = February, 3 = March, ..., 13 = January (+ 1), ...
month_end	= 12	#12 = December, 13 = January (+ 1), 14 = February (+ 1), ...

sig_level	= 95.0
#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/UV_trend_depth_'+str(depth_min)+'-'+str(depth_max)+'_m_year_'+str(year_start)+'-'+str(year_end)+'_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

lon		= HEAT_data.variables['lon'][:] 			
lat		= HEAT_data.variables['lat'][:] 			
vel_trend	= HEAT_data.variables['VEL_trend'][:] 
vel_trend_sig	= HEAT_data.variables['VEL_trend_sig'][:] 

HEAT_data.close()

#Set all the non-significant trends to masked elements
vel_trend	= ma.masked_where(vel_trend_sig < (sig_level / 100.0), vel_trend)
vel_trend_sig	= ma.masked_where(vel_trend_sig < (sig_level / 100.0), vel_trend_sig)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots()

m = Basemap(
projection = 'merc',
llcrnrlat=-5, urcrnrlat=25,
llcrnrlon=270, urcrnrlon=330,
resolution='l', area_thresh=0.01
) 

m.drawcoastlines(linewidth=0.2)
m.drawcountries()
m.fillcontinents(color='#cc9966',lake_color='#99ffff')
par = m.drawparallels(np.arange(-80,81,10),labels=[1,0,0,0])
mer = m.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1])

x, y	= m(lon, lat)

levels	= np.arange(-1.0, 1.01, 0.1)

CS	= m.contourf(x, y, vel_trend, levels, extend = 'both', cmap = 'RdBu_r')
cbar 	= m.colorbar(CS, ticks = np.arange(-1, 1.1, 0.5))
cbar.set_label('Velocity trend (cm s$^{-1}$ per decade)')

ax.set_title('LR-CESM')

show()