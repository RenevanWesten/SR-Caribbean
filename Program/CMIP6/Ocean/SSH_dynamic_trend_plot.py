#Programs plots the local DSL trend

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
directory_CMIP6		= '../../../Data/CMIP6/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

year_start		= 1
year_end		= 101

month_start		= 1 	#1 = January, 2 = February, 3 = March, ..., 13 = January (+ 1), ...
month_end		= 12	#12 = December, 13 = January (+ 1), 14 = February (+ 1), ...

sig_level		= 95.0
#-----------------------------------------------------------------------------------------

	
#Get the model names and path
models = glob.glob(directory_CMIP6+'*')
models.sort()

for model_i in range(len(models)):
	#Only retain the model names
	models[model_i]	= models[model_i][len(directory_CMIP6):]
	

for model_i in range(len(models)):
	#Get for each model the file

	HEAT_data 		= netcdf.Dataset(directory_CMIP6+models[model_i]+'/Ocean/SSH_dynamic_trend_year_'+str(year_start)+'-'+str(year_end)+'_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

	lon			= HEAT_data.variables['lon'][:] 			
	lat			= HEAT_data.variables['lat'][:] 			
	SSH_trend		= HEAT_data.variables['SSH_trend'][:] 
	SSH_trend_sig		= HEAT_data.variables['SSH_trend_sig'][:] 

	HEAT_data.close()

	#Set all the non-significant trends to masked elements
	SSH_trend_sig	      = ma.masked_where(SSH_trend_sig < (sig_level / 100.0), SSH_trend_sig)

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


	if models[model_i] == 'CESM2' or models[model_i] == 'CESM2-WACCM' or models[model_i] == 'SAM0-UNICON':
		x, y = m(lon, lat)
		
	else:
		x, y	= np.meshgrid(lon, lat)
		x, y	= m(x, y)
		
	levels	= np.arange(0, 4.01, 0.2)

	CS		= m.contourf(x, y, SSH_trend, levels, extend = 'both', cmap = 'Reds')
	cbar 	= m.colorbar(CS, ticks = np.arange(0, 4.01, 1))
	cbar.set_label('DSL trend (mm year$^{-1}$)')

	CS = m.contourf(x, y, SSH_trend_sig, levels, extend = 'both', colors='none', linewidth = 2.0, hatches= '\\')
	CS = m.contourf(x, y, SSH_trend_sig, levels, extend = 'both', colors='none', linewidth = 2.0, hatches= '//')

	ax.set_title(models[model_i])
	
	#Region 1
	x = [278, 283, 283, 278, 278]
	y = [10, 10, 14, 14, 10]

	x, y 	= m(x, y)
	m.plot(x, y, linestyle = '-', color = 'black', linewidth = 3)

        #Region 2
	x = [300, 305, 305, 300, 300]
	y = [15, 15, 19, 19, 15]

	x, y 	= m(x, y)
	m.plot(x, y, linestyle = '-', color = 'black', linewidth = 3)

	x, y = m(280.5, 12)
	plt.text(x,y,'1', color ='k',fontsize=20,ha='center',va='center')

	x, y = m(302.5, 17)
	plt.text(x,y,'2', color ='k',fontsize=20,ha='center',va='center')	

show()