#Program fits the GEV distributions

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy.stats import genextreme
from scipy import stats
from matplotlib.colors import LogNorm

#Making pathway to folder with all data
directory 	= '../../../Data/HR-CESM/'

def ReturnValue(prob, shape, loc = 0.0, scale = 1.0):
	"""Return the return value at a given probability"""
	
	return loc - (scale / shape)* (1.0 - (-np.log(1 - prob))**(-shape))

def ReturnTime(value, shape, loc = 0.0, scale = 1.0):
	"""Returns the return time of a given event"""

	prob	= 1.0 - np.exp(-(1.0 + shape * ( (value - loc) / scale))**(-1.0 / shape))

	return 1.0 / prob

def TrendRemover(time, data, trend_type):
	"""Removes trend of choice"""
	
	rank = polyfit(time, data, trend_type)
	fitting = 0.0 
		
	for rank_i in range(len(rank)):
			
		fitting += rank[rank_i] * (time**(len(rank) - 1 - rank_i))

	data -= fitting
	
	return data

def MonthConverter(X):
        """Converts number of months to number of years"""
	V = ((X)/12.0)

	return ["%.0f" % z for z in V]

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
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

trend_type	= 0
section		= 30	#Number of years for GEV fit
#-----------------------------------------------------------------------------------------
HEAT_data = netcdf.Dataset(directory+'Ocean/SSH_NBC_monthly_maximum.nc', 'r')

#Writing data to correct variable	
time_all	= HEAT_data.variables['time'][:]     	 
ssh		= HEAT_data.variables['SSH'][:] 

HEAT_data.close()

time_year, ssh_year	= YearlyConverter(time_all, ssh)
number_years		= len(time_all) / 12 + 1 - section

data_all		= ma.masked_all((number_years, section * 12))
shape_all		= ma.masked_all((number_years))
loc_all			= ma.masked_all((number_years))
scale_all		= ma.masked_all((number_years))
time_year		= ma.masked_all(101)
return_time_all		= ma.masked_all(number_years)
return_level_all	= ma.masked_all(number_years)

for year_i in range(len(time_year)):
	time_year[year_i]	= datetime.datetime(2000 + year_i, 7, 1).toordinal()

#-----------------------------------------------------------------------------------------

for year_i in range(number_years):
	#Take a sliding-window over each year

	#Get the data for each period of time
	data	= ssh[year_i * 12:(year_i + section) * 12]

	#Time-array same length
	time	= np.arange(len(data))

	if trend_type > 0:
	       #Remove the trend to get extreme values stationary before the fit
                data	= TrendRemover(time, data, trend_type)

	#Fit GEV distribution to stationary data
	shape, loc, scale	= genextreme.fit(data)

	#Reverse shape (Python reverses the shape fit...)
	shape			= -shape

	#Save the fitted parameters 
	data_all[year_i]	= data
	shape_all[year_i]	= shape
	loc_all[year_i]		= loc
	scale_all[year_i]	= scale

	#Save return time and level and convert to years
	return_level_all[year_i]	= ReturnValue(1/60.0, shape, loc, scale)
	return_time_all[year_i]		= ReturnTime(return_level_all[0], shape, loc, scale) / 12.0

#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(directory+'Ocean/SSH_NBC_monthly_maximum.nc', 'r')

#Writing data to correct variable	
time_all	= HEAT_data.variables['time'][:]     	 
ssh		= HEAT_data.variables['SSH'][:] 

HEAT_data.close()

print return_level_all[0], return_level_all[70]
#-----------------------------------------------------------------------------------------

ssh_level	= return_level_all[0]

fig, ax = subplots()

graph_ssh		= ax.plot_date(time_all, ssh, '-k', linewidth = 1.5, label = '$\eta_M^{Max}$')

graph_return_level	= plot_date([time_year[0], time_year[29]], [return_level_all[0], return_level_all[0]], '-r', linewidth = 3.0, label = '$\eta_M^{Max}$ 1:5 year')
graph_return_level	= plot_date([time_year[35], time_year[64]], [return_level_all[35], return_level_all[35]], '-r', linewidth = 3.0, label = '$\eta_M^{Max}$ 1:5 year')
graph_return_level	= plot_date([time_year[70], time_year[99]], [return_level_all[70], return_level_all[70]], '-r', linewidth = 3.0, label = '$\eta_M^{Max}$ 1:5 year')

ax.set_xlabel('Model year')
ax.set_ylabel('Sea surface height (cm)')
ax.set_ylim(0, 55)
ax.grid()

ax2 			= ax.twinx()
graph_return_time	= plot_date([time_year[0], time_year[29]], [return_time_all[0], return_time_all[0]], '-b', linewidth = 3.0, label = r'$\tau_{\mathrm{return}}$ '+str(round(ssh_level, 1))+' cm')
graph_return_time	= plot_date([time_year[35], time_year[64]], [return_time_all[35], return_time_all[35]], '-b', linewidth = 3.0, label = r'$\tau_{\mathrm{return}}$ '+str(round(ssh_level, 1))+' cm')
graph_return_time	= plot_date([time_year[70], time_year[99]], [return_time_all[70], return_time_all[70]], '-b', linewidth = 3.0, label = r'$\tau_{\mathrm{return}}$ '+str(round(ssh_level, 1))+' cm')

ax2.set_ylabel('Return time (years)', color = 'b')
ax2.set_yscale('log')
ax2.set_ylim(1, 1000)
#ax2.set_yticks([0, 20, 40, 60, 80, 100])

for tl in ax2.get_yticklabels():
    tl.set_color('b')

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

graphs	      = graph_ssh + graph_return_level + graph_return_time

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right',
 		  ncol=3, fancybox=True, shadow=False, numpoints = 1)

ax.set_title('c) HR-CESM')

#-----------------------------------------------------------------------------------------

prob = np.linspace(1.0 / 10000.0, 0.9999, 100000)

z_p	= ReturnValue(prob, shape_all[0], loc_all[0], scale_all[0])
z_p_50	= ReturnValue(prob, shape_all[35], loc_all[35], scale_all[35])
z_p_100	= ReturnValue(prob, shape_all[70], loc_all[70], scale_all[70])

freq = np.linspace(1.0 / len(time), 1, len(time))

fig, ax = subplots()

graph_2000	= ax.plot(1.0 / prob, z_p, '-k', linewidth = 2.0, label = 'Model year 2000 - 2029')
graph_2050	= ax.plot(1.0 / prob, z_p_50, '-b', linewidth = 2.0, label = 'Model year 2035 - 2064')
graph_2100	= ax.plot(1.0 / prob, z_p_100, '-r', linewidth = 2.0, label = 'Model year 2070 - 2099')

ax.plot(1.0/freq , sorted(data_all[0])[::-1], 'ok')
ax.plot(1.0/freq , sorted(data_all[35])[::-1], 'ob')
ax.plot(1.0/freq , sorted(data_all[70])[::-1], 'or')

	
ax.set_xlim(1, 1000)
ax.set_ylim(0, 55)
ax.grid(axis = 'y')
ax.set_xlabel('Return time (months)')
ax.set_ylabel('Sea surface height (cm)')
ax.set_xscale('log')

ax2 = ax.twiny()

ax2.set_xlim(ax.get_xlim())
ax2.set_xscale('log')

new_tick_locations = np.array([12, 60, 120, 600])
ax2.set_xlabel('Return time (years)')

ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(MonthConverter(new_tick_locations))
ax2.grid()


graphs	      = graph_2000 + graph_2050 + graph_2100

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower right',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)


ax.text(0.01, 0.95, 'HR-CESM', size = 16, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

show()

