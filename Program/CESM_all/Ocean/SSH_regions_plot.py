#Program plots the dynamical SLR components for region 1 and 2

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
directory_cesm 		    = '../../../Data/HR-CESM/'
directory_cesm_low          = '../../../Data/LR-CESM/'
directory_cesm_control 	    = '../../../Data/HR-CESM_Control/'
directory_cesm_low_control  = '../../../Data/LR-CESM_Control/'

def ReadinData(filename):
	"""Reads-in the data"""
	HEAT_data 		= netcdf.Dataset(filename, 'r')

	time			= HEAT_data.variables['time'][:] 
	ssh_1			= HEAT_data.variables['SSH_region_1'][:] 						
	ssh_2			= HEAT_data.variables['SSH_region_2'][:] 

	HEAT_data.close()

	return time, ssh_1, ssh_2

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

month_start	= 1 	#1 = January, 2 = February, 3 = March, ..., 13 = January (+ 1), ...
month_end	= 12	#12 = December, 13 = January (+ 1), 14 = February (+ 1), ...

#-----------------------------------------------------------------------------------------
time_cesm, ssh_1_cesm, ssh_2_cesm	        = ReadinData(directory_cesm+'Ocean/SSH_regions.nc')
time_year_cesm, ssh_year_1_cesm			= YearlyConverter(time_cesm, ssh_1_cesm, month_start, month_end)
time_year_cesm, ssh_year_2_cesm			= YearlyConverter(time_cesm, ssh_2_cesm, month_start, month_end)

ssh_year_1_cesm					= ssh_year_1_cesm - ssh_year_1_cesm[0]
ssh_year_2_cesm					= ssh_year_2_cesm - ssh_year_2_cesm[0]	

print 'HR-CESM SSH 1:', SignificantTrend(time_year_cesm, ssh_year_1_cesm)
print 'HR-CESM SSH 2:', SignificantTrend(time_year_cesm, ssh_year_2_cesm)

#-----------------------------------------------------------------------------------------
time_cesm, steric_1_cesm, steric_2_cesm		= ReadinData(directory_cesm+'/Ocean/SSH_steric_regions.nc')
time_year_cesm, steric_year_1_cesm		= YearlyConverter(time_cesm, steric_1_cesm, month_start, month_end)
time_year_cesm, steric_year_2_cesm		= YearlyConverter(time_cesm, steric_2_cesm, month_start, month_end)

steric_year_1_cesm				= steric_year_1_cesm - steric_year_1_cesm[0]
steric_year_2_cesm				= steric_year_2_cesm - steric_year_2_cesm[0]	

print 'HR-CESM STERIC 1:', SignificantTrend(time_year_cesm, steric_year_1_cesm)
print 'HR-CESM STERIC 2:', SignificantTrend(time_year_cesm, steric_year_2_cesm)

#-----------------------------------------------------------------------------------------

time_cesm_low, ssh_1_cesm_low, ssh_2_cesm_low	= ReadinData(directory_cesm_low+'/Ocean/SSH_regions.nc')
time_year_cesm_low, ssh_year_1_cesm_low		= YearlyConverter(time_cesm_low, ssh_1_cesm_low, month_start, month_end)
time_year_cesm_low, ssh_year_2_cesm_low		= YearlyConverter(time_cesm_low, ssh_2_cesm_low, month_start, month_end)

ssh_year_1_cesm_low				= ssh_year_1_cesm_low - ssh_year_1_cesm_low[0]
ssh_year_2_cesm_low				= ssh_year_2_cesm_low - ssh_year_2_cesm_low[0]

print 'LR-CESM SSH 1:', SignificantTrend(time_year_cesm_low, ssh_year_1_cesm_low)
print 'LR-CESM SSH 2:', SignificantTrend(time_year_cesm_low, ssh_year_2_cesm_low)

#-----------------------------------------------------------------------------------------
time_cesm, steric_1_cesm_low, steric_2_cesm_low	= ReadinData(directory_cesm_low+'/Ocean/SSH_steric_regions.nc')
time_year_cesm_low, steric_year_1_cesm_low	= YearlyConverter(time_cesm_low, steric_1_cesm_low, month_start, month_end)
time_year_cesm_low, steric_year_2_cesm_low	= YearlyConverter(time_cesm_low, steric_2_cesm_low, month_start, month_end)

steric_year_1_cesm_low				= steric_year_1_cesm_low - steric_year_1_cesm_low[0]
steric_year_2_cesm_low				= steric_year_2_cesm_low - steric_year_2_cesm_low[0]	

print 'LR-CESM STERIC 1:', SignificantTrend(time_year_cesm, steric_year_1_cesm_low)
print 'LR-CESM STERIC 2:', SignificantTrend(time_year_cesm, steric_year_2_cesm_low)
print

#-----------------------------------------------------------------------------------------
print '-----------------------------------------------------------------------------------------'
print 'HR-CESM SSH 1:', SignificantTrend(time_year_cesm, ssh_year_1_cesm)[0] * 10.0, '+/-', stats.linregress(np.arange(101), ssh_year_1_cesm)[4] * 10.0, ', sig = ', SignificantTrend(time_year_cesm, ssh_year_1_cesm)[2]
print 'HR-CESM SSH 2:', SignificantTrend(time_year_cesm, ssh_year_2_cesm)[0] * 10.0, '+/-', stats.linregress(np.arange(101), ssh_year_2_cesm)[4] * 10.0, ', sig = ', SignificantTrend(time_year_cesm, ssh_year_2_cesm)[2]
'-----------------------------------------------------------------------------------------'
print
print '-----------------------------------------------------------------------------------------'
print 'LR-CESM SSH 1:', SignificantTrend(time_year_cesm, ssh_year_1_cesm_low)[0] * 10.0, '+/-', stats.linregress(np.arange(101), ssh_year_1_cesm_low)[4] * 10.0, ', sig = ', SignificantTrend(time_year_cesm, ssh_year_1_cesm_low)[2]
print 'LR-CESM SSH 2:', SignificantTrend(time_year_cesm, ssh_year_2_cesm_low)[0] * 10.0, '+/-', stats.linregress(np.arange(101), ssh_year_2_cesm_low)[4] * 10.0, ', sig = ', SignificantTrend(time_year_cesm, ssh_year_2_cesm_low)[2]
'-----------------------------------------------------------------------------------------'

#-----------------------------------------------------------------------------------------

HEAT_data 	                                = netcdf.Dataset(directory_cesm+'Ocean/SSH_global_steric.nc', 'r')

time		                                = HEAT_data.variables['time'][:] 
steric_global_cesm                              = HEAT_data.variables['SSH'][:] * 100.0

HEAT_data.close()

time_year_cesm, steric_global_year_cesm	        = YearlyConverter(time, steric_global_cesm, month_start, month_end)
steric_global_year_cesm			        = steric_global_year_cesm - steric_global_year_cesm[0]

#-----------------------------------------------------------------------------------------

HEAT_data 	                                = netcdf.Dataset(directory_cesm_control+'Ocean/SSH_global_steric.nc', 'r')

time		                                = HEAT_data.variables['time'][:] 
steric_global_cesm_control                      = HEAT_data.variables['SSH'][:] * 100.0

HEAT_data.close()

time_year_2, steric_global_year_cesm_control	= YearlyConverter(time, steric_global_cesm_control, month_start, month_end)
steric_global_year_cesm_control			= steric_global_year_cesm_control - steric_global_year_cesm_control[0]

#Adjust for drift
steric_global_year_cesm		                = steric_global_year_cesm - steric_global_year_cesm_control
ssh_total_year_1_cesm	                        = ssh_year_1_cesm + steric_global_year_cesm	
ssh_total_year_2_cesm	                        = ssh_year_2_cesm + steric_global_year_cesm

print 'HR-CESM:', SignificantTrend(time_year_cesm, steric_global_year_cesm)
#-----------------------------------------------------------------------------------------

HEAT_data 	                                = netcdf.Dataset(directory_cesm_low+'Ocean/SSH_global_steric.nc', 'r')

time_low		                        = HEAT_data.variables['time'][:] 
steric_global_cesm_low                          = HEAT_data.variables['SSH'][:] * 100.0

HEAT_data.close()

time_year_cesm_low, steric_global_year_cesm_low = YearlyConverter(time_low, steric_global_cesm_low, month_start, month_end)
steric_global_year_cesm_low			= steric_global_year_cesm_low - steric_global_year_cesm_low[0]

#-----------------------------------------------------------------------------------------

HEAT_data 	                                = netcdf.Dataset(directory_cesm_low_control+'Ocean/SSH_global_steric.nc', 'r')

time		                                = HEAT_data.variables['time'][:] 
steric_global_cesm_low_control                  = HEAT_data.variables['SSH'][:] * 100.0

HEAT_data.close()

time_year_2, steric_global_year_cesm_low_control= YearlyConverter(time, steric_global_cesm_low_control, month_start, month_end)
steric_global_year_cesm_low_control		= steric_global_year_cesm_low_control - steric_global_year_cesm_low_control[0]

#Adjust for drift
steric_global_year_cesm_low		        = steric_global_year_cesm_low - steric_global_year_cesm_low_control
ssh_total_year_1_cesm_low	                = ssh_year_1_cesm_low + steric_global_year_cesm_low	
ssh_total_year_2_cesm_low	                = ssh_year_2_cesm_low + steric_global_year_cesm_low

print 'LR-CESM:', SignificantTrend(time_year_cesm_low, steric_global_year_cesm_low)
#-----------------------------------------------------------------------------------------

dynamic_trend_cesm_1	        = SignificantTrend(time_year_cesm, ssh_total_year_1_cesm)[0] * 10.0 / 1.8
dynamic_trend_cesm_2	        = SignificantTrend(time_year_cesm, ssh_total_year_2_cesm)[0] * 10.0 / 1.8
dynamic_trend_cesm_low_1        = SignificantTrend(time_year_cesm, ssh_total_year_1_cesm_low)[0] * 10.0 / 1.8
dynamic_trend_cesm_low_2	= SignificantTrend(time_year_cesm, ssh_total_year_2_cesm_low)[0] * 10.0 / 1.8

#-----------------------------------------------------------------------------------------

OBP_1_cesm	= steric_global_year_cesm - steric_year_1_cesm
OBP_2_cesm	= steric_global_year_cesm - steric_year_2_cesm
OBP_1_cesm_low	= steric_global_year_cesm_low - steric_year_1_cesm_low
OBP_2_cesm_low	= steric_global_year_cesm_low - steric_year_2_cesm_low

print 'HR-CESM OBP 1:', SignificantTrend(time_year_cesm, OBP_1_cesm)
print 'HR-CESM OBP 2:', SignificantTrend(time_year_cesm, OBP_2_cesm)
print 'LR-CESM OBP 1:', SignificantTrend(time_year_cesm, OBP_1_cesm_low)
print 'LR-CESM OBP 2:', SignificantTrend(time_year_cesm, OBP_2_cesm_low)
#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

HR_CESM_SSH_1_graph	    = ax.plot_date(time_year_cesm, ssh_year_1_cesm, '-b', linewidth = 2.0, label = '$\eta_{M}^{1}$')
HR_CESM_SSH_2_graph	    = ax.plot_date(time_year_cesm, ssh_year_2_cesm, '--b', linewidth = 2.0, label = '$\eta_{M}^{2}$')
HR_CESM_SSH_total_1_graph   = ax.plot_date(time_year_cesm, ssh_total_year_1_cesm, '-k', linewidth = 2.0, label = '$\eta^{1}$')
HR_CESM_SSH_total_2_graph   = ax.plot_date(time_year_cesm, ssh_total_year_2_cesm, '--k', linewidth = 2.0, label = '$\eta^{2}$')
HR_CESM_SSH_global_graph    = ax.plot_date(time_year_cesm, steric_global_year_cesm, '-r', linewidth = 3.0, label = '$\eta_{S}^{g}$')

ax.set_xlabel('Model year')
ax.set_ylabel('Sea level anomaly (cm)')
ax.set_ylim(-10, 25)
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

ax.set_title('a) HR-CESM')

graphs	      = HR_CESM_SSH_1_graph + HR_CESM_SSH_2_graph + HR_CESM_SSH_global_graph + HR_CESM_SSH_total_1_graph + HR_CESM_SSH_total_2_graph

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

HR_CESM_STERIC_1_graph	= ax.plot_date(time_year_cesm, steric_year_1_cesm, '-r', linewidth = 2.0, label = '$\eta_{S}^{1}$')
HR_CESM_STERIC_2_graph	= ax.plot_date(time_year_cesm, steric_year_2_cesm, fmt = '-', color = 'darkred', linewidth = 2.0, label = '$\eta_{S}^{2}$')
HR_CESM_OBP_1_graph	= ax.plot_date(time_year_cesm, OBP_1_cesm, '-b', linewidth = 2.0, label = '$\eta_{S}^{M, 1}$')
HR_CESM_OBP_2_graph	= ax.plot_date(time_year_cesm, OBP_2_cesm, '-c', linewidth = 2.0, label = '$\eta_{S}^{M, 2}$')

ax.set_xlabel('Model year')
ax.set_ylabel('Sea level anomaly (cm)')
ax.set_ylim(-10, 25)
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

ax.set_title('c) HR-CESM')

graphs	      = HR_CESM_STERIC_1_graph + HR_CESM_STERIC_2_graph + HR_CESM_OBP_1_graph + HR_CESM_OBP_2_graph

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

LR_CESM_SSH_1_graph		= ax.plot_date(time_year_cesm_low, ssh_year_1_cesm_low, '-b', linewidth = 2.0, label = '$\eta_{M}^{1}$')
LR_CESM_SSH_2_graph		= ax.plot_date(time_year_cesm_low, ssh_year_2_cesm_low, '--b', linewidth = 2.0, label = '$\eta_{M}^{2}$')
LR_CESM_SSH_total_1_graph	= ax.plot_date(time_year_cesm_low, ssh_total_year_1_cesm_low, '-k', linewidth = 2.0, label = '$\eta^{1}$')
LR_CESM_SSH_total_2_graph	= ax.plot_date(time_year_cesm_low, ssh_total_year_2_cesm_low, '--k', linewidth = 2.0, label = '$\eta^{2}$')
LR_CESM_SSH_global_graph	= ax.plot_date(time_year_cesm_low, steric_global_year_cesm_low, '-r', linewidth = 3.0, label = '$\eta_{S}^{g}$')

ax.set_xlabel('Model year')
ax.set_ylabel('Sea level anomaly (cm)')
ax.set_ylim(-10, 25)
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

ax.set_title('b) LR-CESM')

graphs	      = LR_CESM_SSH_1_graph + LR_CESM_SSH_2_graph + LR_CESM_SSH_global_graph + LR_CESM_SSH_total_1_graph + LR_CESM_SSH_total_2_graph

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

LR_CESM_STERIC_1_graph	= ax.plot_date(time_year_cesm_low, steric_year_1_cesm_low, '-r', linewidth = 2.0, label = '$\eta_{S}^{1}$')
LR_CESM_STERIC_2_graph	= ax.plot_date(time_year_cesm_low, steric_year_2_cesm_low, fmt = '-', color = 'darkred', linewidth = 2.0, label = '$\eta_{S}^{2}$')
LR_CESM_OBP_1_graph	= ax.plot_date(time_year_cesm_low, OBP_1_cesm_low, '-b', linewidth = 2.0, label = '$\eta_{S}^{M, 1}$')
LR_CESM_OBP_2_graph	= ax.plot_date(time_year_cesm_low, OBP_2_cesm_low, '-c', linewidth = 2.0, label = '$\eta_{S}^{M, 2}$')

ax.set_xlabel('Model year')
ax.set_ylabel('Sea level anomaly (cm)')
ax.set_ylim(-10, 25)
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

ax.set_title('d) LR-CESM')

graphs	      = LR_CESM_STERIC_1_graph + LR_CESM_STERIC_2_graph + LR_CESM_OBP_1_graph + LR_CESM_OBP_2_graph

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)


show()
