#Program plots the AMOC strength

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats

#Making pathway to folder with all data
directory_CMIP6		= '../../../Data/CMIP6/'
directory_cesm 		= '../../../Data/HR-CESM/'
directory_cesm_low 	= '../../../Data/LR-CESM/'

def ReadinData(filename):
	"""Reads-in the data"""
	HEAT_data 		= netcdf.Dataset(filename, 'r')

	time			= HEAT_data.variables['time'][:] 
	transport		= HEAT_data.variables['Transport'][:] 

	HEAT_data.close()

	return time, transport

def YearlyConverter(time, data, month_start = 1, month_end = 12):
	"""Determines the averaged value over different months of choice,
	default is set to January - December"""

	#Take twice the amount of years for the month day (no leap year)
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31., 31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days[month_start - 1:month_end]
	month_days	= month_days / np.sum(month_days)

	if month_end <= 12:
		#Normal average over a single year, for example, January 100 - December 100
		time_year		= np.zeros(len(time) / 12)

	else:
		#If you take the average, for example, over November 100 - May 101
		#Take year 101 as the average over this period
		#There is one year less compared to the period analysed
		time_year		= np.zeros(len(time) / 12 - 1)

	#-----------------------------------------------------------------------------------------
	data_year	= ma.masked_all(len(time_year))

	for year_i in range(len(time_year)):
		#Determine the average over the selected months

		#The year is defined as the current year
		year			= int(str(datetime.date.fromordinal(int(time[year_i * 12])))[0:4])

		if month_end	>= 13:
			#If average is taken over, for example, November 100 - May 101, the year is defined as 101
			year = year + 1

		time_year[year_i] 	= datetime.datetime(year, 1, 1).toordinal()

		#Determine the time mean over the months of choice
		data_year[year_i]		= np.sum(data[year_i * 12 + month_start - 1: year_i * 12 + month_end] * month_days, axis = 0)

	return time_year, data_year

def TrendRemover(time, data, trend_type):
	"""Removes trend of choice"""
	
	rank = polyfit(time, data, trend_type)
	fitting = 0.0 
		
	for rank_i in range(len(rank)):
			
		fitting += rank[rank_i] * (time**(len(rank) - 1 - rank_i))

	data -= fitting
	
	return data

def SignificantTrend(time, data):
	"""Finds whether trend is significant
	Returns the trend and if it significant (= 1)"""

	#Set time similar to Santer et al. (2000), time array from 1 till N
	#Statistical significance of trends and trend differences in layer-average atmospheric temperature time series
	time		= np.arange(1, len(time) + 1)

	#Determine the detrended time series
	trend, base 	= polyfit(time, data, 1)
	data_res	= data - ((trend * time) + base)

	#Effective sample size, based on the lag-1 correlation
	corr_1		= np.corrcoef(data_res[:-1], data_res[1:])[0, 1]
	N_eff		= int(len(time) * (1.0 - corr_1) / (1.0 + corr_1))

	#Determine the variance of the anomalies
	data_var	= np.sum(data_res**2.0) / (N_eff - 2.0)

	#Determine the standard error
	standard_error	=  np.sqrt(data_var) / np.sqrt(np.sum((time - np.mean(time))**2.0))

	#Determine the Student-T value
	t_value		= trend / standard_error

	#Get the significance levels and the corresponding critical values (two-sided)
	sig_levels 	= np.arange(50, 100, 0.5) / 100.0
	t_crit 		= stats.t.ppf((1.0 + sig_levels) / 2.0, N_eff - 2)

	#Get the indices where the significance is exceeding the critical values
	sig_index	= np.where(fabs(t_value) > t_crit)[0]
	significant	= 0.0

	if len(sig_index) > 0:
		#If there are significance values, take the highest significant level
		significant = sig_levels[sig_index[-1]]

	return trend, np.sqrt(standard_error), significant

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0 
depth_max	= 1000

month_start	= 1 	#1 = January, 2 = February, 3 = March, ..., 13 = January (+ 1), ...
month_end	= 12	#12 = December, 13 = January (+ 1), 14 = February (+ 1), ...

#-----------------------------------------------------------------------------------------
	
#Get the model names and path
models = glob.glob(directory_CMIP6+'*')
models.sort()

for model_i in range(len(models)):
	#Only retain the model names
	models[model_i]	= models[model_i][len(directory_CMIP6):]

for model_i in range(len(models)):
	#Get for each model the corresponding file
	time, transport            = ReadinData(directory_CMIP6+models[model_i]+'/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'_m.nc')
	
	#Determine yearly averaged time series
	time_year, transport_year = YearlyConverter(time, transport, month_start, month_end)

	if model_i == 0:
		#For the CMIP6 mean
		transport_all	= ma.masked_all((len(models), len(time_year)))

        #Save all the CMIP6 output in an array
	transport_all[model_i]	= transport_year

	#Print the mean over the first 10 years
	print models[model_i],':', np.mean(transport_year[:10])	

print


time_cesm, transport_cesm			= ReadinData(directory_cesm+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'_m.nc')
time_year_cesm, transport_year_cesm		= YearlyConverter(time_cesm, transport_cesm, month_start, month_end)

print 'HR-CESM:', np.mean(transport_year_cesm[:10])
#-----------------------------------------------------------------------------------------

time_cesm_low, transport_cesm_low	        = ReadinData(directory_cesm_low+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'_m.nc')
time_year_cesm_low, transport_year_cesm_low	= YearlyConverter(time_cesm_low, transport_cesm_low, month_start, month_end)
print 'LR-CESM:', np.mean(transport_year_cesm_low[:10])
#-----------------------------------------------------------------------------------------

print
print 'HR-CESM:', SignificantTrend(time_year_cesm, transport_year_cesm)
print 'LR-CESM low:', SignificantTrend(time_year_cesm_low, transport_year_cesm_low)
print 'CMIP6:', SignificantTrend(time_year, np.mean(transport_all, axis = 0))
print

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

CMIP6_graph	= ax.plot_date(time_year, np.mean(transport_all, axis = 0), '-k', linewidth = 2.0, label = 'CMIP6')
ax.fill_between(time_year, np.min(transport_all, axis = 0), np.max(transport_all, axis = 0), color = 'k', alpha = 0.25)

HR_CESM_graph	= ax.plot_date(time_year, transport_year_cesm, '-r', linewidth = 2.0, label = 'HR-CESM')
LR_CESM_graph	= ax.plot_date(time_year, transport_year_cesm_low, '-b', linewidth = 2.0, label = 'LR-CESM')

ax.set_xlabel('Model year')
ax.set_ylabel('Transport (Sv)')
ax.set_ylim(0, 25)
ax.grid()

ax.set_xticks([	datetime.datetime(1, 1, 1).toordinal(),
		datetime.datetime(10, 1, 1).toordinal(),
		datetime.datetime(20, 1, 1).toordinal(),
		datetime.datetime(30, 1, 1).toordinal(), 
		datetime.datetime(40, 1, 1).toordinal(),
		datetime.datetime(50, 1, 1).toordinal(), 
		datetime.datetime(60, 1, 1).toordinal(),
		datetime.datetime(70, 1, 1).toordinal(),
		datetime.datetime(80, 1, 1).toordinal(),
		datetime.datetime(90, 1, 1).toordinal(),
		datetime.datetime(100, 1, 1).toordinal()])

graphs	      = HR_CESM_graph + LR_CESM_graph + CMIP6_graph 

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)

ax.set_title('b) AMOC strength')

ax2 = fig.add_axes([0.15, 0.10, 0.2, 0.2])

ax2.set_ylim(-0.1, 1.1)
ax2.set_xlim(0, 2.6)
ax2.axis('off')

x_legend	= np.arange(0, 1.51, 0.1)
ax2.fill_between(x_legend, 0.20, 0.80, color = 'k', alpha = 0.25)
ax2.plot(x_legend, 0.5 + np.zeros(len(x_legend)), linestyle = '-', color = 'k', linewidth = 3.0)


ax2.text(2.30, 0.20,'Min', color ='k',fontsize=15,ha='left',va='center')
ax2.plot([1.5, 2.28], [0.20, 0.20], '--k', linewidth = 0.5)

ax2.text(1.9, 0.5,'Mean', color ='k',fontsize=15,ha='left',va='center')
ax2.plot([1.5, 1.88], [0.5, 0.5], '--k', linewidth = 0.5)

ax2.text(2.30, 0.80,'Max', color ='k',fontsize=15,ha='left',va='center')
ax2.plot([1.5, 2.28], [0.80, 0.80], '--k', linewidth = 0.5)


show()
