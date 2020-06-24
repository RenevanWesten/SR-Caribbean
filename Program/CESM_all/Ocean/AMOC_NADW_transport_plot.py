#Program plots the AMOC and NADW strength

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

depth_min_1 	= 0 
depth_max_1	= 1000
depth_min_2 	= 1500
depth_max_2	= 4000

month_start	= 1 	#1 = January, 2 = February, 3 = March, ..., 13 = January (+ 1), ...
month_end	= 12	#12 = December, 13 = January (+ 1), 14 = February (+ 1), ...

#-----------------------------------------------------------------------------------------
#------------------------------------------AMOC-------------------------------------------
#-----------------------------------------------------------------------------------------

print
print '-----------------------------------------------------------------------------------------'

time_cesm, transport_cesm		= ReadinData(directory_cesm+'/Ocean/AMOC_transport_depth_'+str(depth_min_1)+'-'+str(depth_max_1)+'_m.nc')
time_year_cesm, AMOC_year_cesm		= YearlyConverter(time_cesm, transport_cesm, month_start, month_end)

print 'HR-CESM AMOC', np.mean(AMOC_year_cesm[:10])
#-----------------------------------------------------------------------------------------

time_cesm_low, transport_cesm_low	= ReadinData(directory_cesm_low+'/Ocean/AMOC_transport_depth_'+str(depth_min_1)+'-'+str(depth_max_1)+'_m.nc')
time_year_cesm_low, AMOC_year_cesm_low	= YearlyConverter(time_cesm_low, transport_cesm_low, month_start, month_end)

print 'LR-CESM AMOC', np.mean(AMOC_year_cesm_low[:10])
print '-----------------------------------------------------------------------------------------'
print
#-----------------------------------------------------------------------------------------
#------------------------------------------NADW-------------------------------------------
#-----------------------------------------------------------------------------------------

print '-----------------------------------------------------------------------------------------'
time_cesm, transport_cesm		= ReadinData(directory_cesm+'/Ocean/NADW_transport_depth_'+str(depth_min_2)+'-'+str(depth_max_2)+'_m.nc')
time_year_cesm, NADW_year_cesm		= YearlyConverter(time_cesm, transport_cesm, month_start, month_end)

print 'HR-CESM NADW', np.mean(NADW_year_cesm[:10])
#-----------------------------------------------------------------------------------------

time_cesm_low, transport_cesm_low	= ReadinData(directory_cesm_low+'/Ocean/NADW_transport_depth_'+str(depth_min_2)+'-'+str(depth_max_2)+'_m.nc')
time_year_cesm_low, NADW_year_cesm_low	= YearlyConverter(time_cesm_low, transport_cesm_low, month_start, month_end)

print 'LR-CESM NADW', np.mean(NADW_year_cesm_low[:10])
print '-----------------------------------------------------------------------------------------'
print
#-----------------------------------------------------------------------------------------

print '-----------------------------------------------------------------------------------------'
print 'HR-CESM AMOC:', SignificantTrend(time_year_cesm, AMOC_year_cesm)[0] * 100.0, '+/-', stats.linregress(np.arange(101), AMOC_year_cesm)[4] * 100.0, ', sig = ', SignificantTrend(time_year_cesm, AMOC_year_cesm)[2]
print 'LR-CESM AMOC:', SignificantTrend(time_year_cesm_low, AMOC_year_cesm_low)[0] * 100.0, '+/-', stats.linregress(np.arange(len(time_year_cesm_low)), AMOC_year_cesm_low)[4] * 100.0, ', sig = ',  SignificantTrend(time_year_cesm_low, AMOC_year_cesm_low)[2]
print '-----------------------------------------------------------------------------------------'
print
print '-----------------------------------------------------------------------------------------'
print 'HR-CESM NADW:', SignificantTrend(time_year_cesm, NADW_year_cesm)[0] * 100.0, '+/-', stats.linregress(np.arange(101), NADW_year_cesm)[4] * 100.0, ', sig = ', SignificantTrend(time_year_cesm, NADW_year_cesm)[2]  
print 'LR-CESM NADW:', SignificantTrend(time_year_cesm_low, NADW_year_cesm_low)[0] * 100.0, '+/-', stats.linregress(np.arange(len(time_year_cesm_low)), NADW_year_cesm_low)[4] * 100.0, ', sig = ',   SignificantTrend(time_year_cesm_low, NADW_year_cesm_low)[2]
print '-----------------------------------------------------------------------------------------'
#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

HR_CESM_graph	= ax.plot_date(time_year_cesm, AMOC_year_cesm, '-r', linewidth = 2.0, label = 'HR-CESM')
LR_CESM_graph	= ax.plot_date(time_year_cesm, AMOC_year_cesm_low, '-b', linewidth = 2.0, label = 'LR-CESM')

ax.plot_date(time_year_cesm, NADW_year_cesm, '--r', linewidth = 2.0, label = 'HR-CESM')
ax.plot_date(time_year_cesm, NADW_year_cesm_low, '--b', linewidth = 2.0, label = 'LR-CESM')

ax.axhline(y = 0, color = 'k', linestyle = ':', linewidth = 3.0)
ax.text(0.9, 0.5, 'AMOC', horizontalalignment='center', verticalalignment='bottom', size = 16, transform=ax.transAxes)
ax.text(0.9, 0.49, 'NADW', horizontalalignment='center', verticalalignment='top', size = 16, transform=ax.transAxes)

ax.set_xlabel('Model year')
ax.set_ylabel('Transport (Sv)')
ax.set_ylim(-26, 26)
ax.grid()

ax.set_xticks([	datetime.datetime(2000, 1, 1).toordinal(),
		datetime.datetime(2010, 1, 1).toordinal(),
		datetime.datetime(2020, 1, 1).toordinal(),
		datetime.datetime(2030, 1, 1).toordinal(), 
		datetime.datetime(2040, 1, 1).toordinal(),
		datetime.datetime(2050, 1, 1).toordinal(), 
		datetime.datetime(2060, 1, 1).toordinal(),
		datetime.datetime(2070, 1, 1).toordinal(),
		datetime.datetime(2080, 1, 1).toordinal(),
		datetime.datetime(2090, 1, 1).toordinal(),
		datetime.datetime(2100, 1, 1).toordinal()])
		
ax.set_yticks(np.arange(-25, 26, 5))
graphs	      = HR_CESM_graph + LR_CESM_graph

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='center left',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)

ax.set_title('b) Meridional volume transport')

show()