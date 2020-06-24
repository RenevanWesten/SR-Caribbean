#Program determines the DSL trend over the two regions

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
directory_CMIP6		    = '../../../Data/CMIP6/'
directory_cesm 		    = '../../../Data/HR-CESM/'
directory_cesm_low          = '../../../Data/LR-CESM/'
directory_cesm_control 	    = '../../../Data/HR-CESM_Control/'
directory_cesm_low_control  = '../../../Data/LR-CESM_Control/'
directory_aviso 	    = '../../../Data/AVISO/'

def ReadinData(filename):
	"""Reads-in the data"""
	HEAT_data 		= netcdf.Dataset(filename, 'r')

	time			= HEAT_data.variables['time'][:] 
	ssh_1			= HEAT_data.variables['SSH_region_1'][:] 						
	ssh_2			= HEAT_data.variables['SSH_region_2'][:] 

	HEAT_data.close()

	return time, ssh_1, ssh_2

def ReadinDataGlobal(filename):
	"""Reads-in the global data"""
	HEAT_data 		= netcdf.Dataset(filename, 'r')

	time			= HEAT_data.variables['time'][:] 
	ssh			= HEAT_data.variables['SSH'][:] * 100.0	#Global steric (cm)				

	HEAT_data.close()

	return time, ssh

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

month_start	= 1 	#1 = January, 2 = February, 3 = March, ..., 13 = January (+ 1), ...
month_end	= 12	#12 = December, 13 = January (+ 1), 14 = February (+ 1), ...

#-----------------------------------------------------------------------------------------
	
#Get the model names and path
models = glob.glob(directory_CMIP6+'*')
models.sort()

for model_i in range(len(models)):
	#Only retain the model names
	models[model_i]	= models[model_i][len(directory_CMIP6):]

#Change the order for the final plots
models.remove('CanESM5')
models.insert(2, 'CanESM5')

for model_i in range(len(models)):
	#Get the model name for each file
	file_ssh		= directory_CMIP6+models[model_i]+'/Ocean/SSH_regions.nc'
	file_ssh_global		= directory_CMIP6+models[model_i]+'/Ocean/SSH_global_steric.nc'
	file_ssh_global_control	= directory_CMIP6+models[model_i]+'/Ocean/SSH_global_steric_Control.nc'

	time, ssh_1, ssh_2			= ReadinData(file_ssh)
	time, ssh_global			= ReadinDataGlobal(file_ssh_global)
	time, ssh_global_control		= ReadinDataGlobal(file_ssh_global_control)

	#-----------------------------------------------------------------------------------------
	time_year, ssh_year_1			= YearlyConverter(time, ssh_1, month_start, month_end)
	time_year, ssh_year_2			= YearlyConverter(time, ssh_2, month_start, month_end)
	time_year, ssh_global_year		= YearlyConverter(time, ssh_global, month_start, month_end)
	time_year, ssh_global_control_year	= YearlyConverter(time, ssh_global_control, month_start, month_end)

	#Set first year to 0 and remove drift
	ssh_year_1			= ssh_year_1 - ssh_year_1[0]
	ssh_year_2			= ssh_year_2 - ssh_year_2[0]
	ssh_global_control_year		= ssh_global_control_year - ssh_global_control_year[0]
	ssh_global_year			= ssh_global_year - ssh_global_year[0] - ssh_global_control_year

	if model_i == 0:
		#For the ensemble mean
		ssh_all_1		= ma.masked_all((len(models), len(time_year)))
		ssh_all_2		= ma.masked_all((len(models), len(time_year)))
		ssh_global_all		= ma.masked_all((len(models), len(time_year)))
		ssh_global_control_all	= ma.masked_all((len(models), len(time_year)))

	ssh_all_1[model_i]		= ssh_year_1
	ssh_all_2[model_i]		= ssh_year_2
	ssh_global_all[model_i]	        = ssh_global_year
	ssh_global_control_all[model_i]	= ssh_global_control_year

#Get the dynamical sea level (SSH + steric)
ssh_total_1		= ssh_all_1 + ssh_global_all
ssh_total_2		= ssh_all_2 + ssh_global_all	

#Now get the AVISO, HR-CESM and LR-CESM
#-----------------------------------------------------------------------------------------
time_aviso, ssh_1_aviso, ssh_2_aviso		= ReadinData(directory_aviso+'/Ocean/SSH_regions.nc')
time_year_aviso, ssh_year_1_aviso		= YearlyConverter(time_aviso, ssh_1_aviso, month_start, month_end)
time_year_aviso, ssh_year_2_aviso		= YearlyConverter(time_aviso, ssh_2_aviso, month_start, month_end)

trend_ssh_aviso_1, base_ssh_aviso_1		= polyfit(np.arange(len(time_year_aviso)), ssh_year_1_aviso, 1) * 10.0
trend_ssh_aviso_2, base_ssh_aviso_2		= polyfit(np.arange(len(time_year_aviso)), ssh_year_2_aviso, 1) * 10.0
#-----------------------------------------------------------------------------------------
time_cesm, ssh_1_cesm, ssh_2_cesm		= ReadinData(directory_cesm+'/Ocean/SSH_regions.nc')
time_year_cesm, ssh_year_1_cesm			= YearlyConverter(time_cesm, ssh_1_cesm, month_start, month_end)
time_year_cesm, ssh_year_2_cesm			= YearlyConverter(time_cesm, ssh_2_cesm, month_start, month_end)

ssh_year_1_cesm					= ssh_year_1_cesm - ssh_year_1_cesm[0]
ssh_year_2_cesm					= ssh_year_2_cesm - ssh_year_2_cesm[0]
#-----------------------------------------------------------------------------------------	
time_cesm_low, ssh_1_cesm_low, ssh_2_cesm_low	= ReadinData(directory_cesm_low+'/Ocean/SSH_regions.nc')
time_year_cesm_low, ssh_year_1_cesm_low		= YearlyConverter(time_cesm_low, ssh_1_cesm_low, month_start, month_end)
time_year_cesm_low, ssh_year_2_cesm_low		= YearlyConverter(time_cesm_low, ssh_2_cesm_low, month_start, month_end)

ssh_year_1_cesm_low				= ssh_year_1_cesm_low - ssh_year_1_cesm_low[0]
ssh_year_2_cesm_low				= ssh_year_2_cesm_low - ssh_year_2_cesm_low[0]

#Get the globally averaged steric sea level rise
#-----------------------------------------------------------------------------------------

time, ssh				                = ReadinDataGlobal(directory_cesm+'Ocean/SSH_global_steric.nc')
time_year_cesm, ssh_global_cesm_year	                = YearlyConverter(time, ssh, month_start, month_end)
ssh_global_cesm_year			                = ssh_global_cesm_year - ssh_global_cesm_year[0]
#-----------------------------------------------------------------------------------------
time, ssh_low				                = ReadinDataGlobal(directory_cesm_low+'Ocean/SSH_global_steric.nc')
time_year_low, ssh_global_cesm_low_year	                = YearlyConverter(time, ssh_low, month_start, month_end)
ssh_global_cesm_low_year		                = ssh_global_cesm_low_year - ssh_global_cesm_low_year[0]
#-----------------------------------------------------------------------------------------
time, ssh				                = ReadinDataGlobal(directory_cesm_control+'Ocean/SSH_global_steric.nc')
time_year_cesm_control, ssh_global_cesm_control_year	= YearlyConverter(time, ssh, month_start, month_end)
ssh_global_cesm_control_year			        = ssh_global_cesm_control_year - ssh_global_cesm_control_year[0]
#-----------------------------------------------------------------------------------------
time, ssh_low				                = ReadinDataGlobal(directory_cesm_low_control+'Ocean/SSH_global_steric.nc')
time_year_low, ssh_global_cesm_low_control_year	        = YearlyConverter(time, ssh_low, month_start, month_end)
ssh_global_cesm_low_control_year		        = ssh_global_cesm_low_control_year - ssh_global_cesm_low_control_year[0]
#-----------------------------------------------------------------------------------------


#Remove drift from models
ssh_global_cesm_year		                        = ssh_global_cesm_year - ssh_global_cesm_control_year
ssh_global_cesm_low_year	                        = ssh_global_cesm_low_year - ssh_global_cesm_low_control_year

#Retain the global trend (in mm / yr)
trend_ssh_global_cesm, base_ssh_global_cesm		= polyfit(np.arange(len(time_year)), ssh_global_cesm_year, 1) * 10.0
trend_ssh_global_cesm_low, base_ssh_global_cesm_low	= polyfit(np.arange(len(time_year)), ssh_global_cesm_low_year, 1) * 10.0

#Add the global steric effect
ssh_total_year_1_cesm		= ssh_year_1_cesm + ssh_global_cesm_year
ssh_total_year_2_cesm		= ssh_year_2_cesm + ssh_global_cesm_year
ssh_total_year_1_cesm_low	= ssh_year_1_cesm_low + ssh_global_cesm_low_year	
ssh_total_year_2_cesm_low	= ssh_year_2_cesm_low + ssh_global_cesm_low_year	

ssh_trend_all_1		= ma.masked_all(len(models))
ssh_trend_all_2		= ma.masked_all(len(models))
ssh_dynamic_trend_all_1	= ma.masked_all(len(models))
ssh_dynamic_trend_all_2	= ma.masked_all(len(models))
ssh_global_trend	= ma.masked_all(len(models))

for model_i in range(len(models)):
	#Determine the trend (mm / year) for the regions
	trend_ssh_1, base_ssh_1			= polyfit(np.arange(len(time_year)), ssh_all_1[model_i], 1) * 10.0
	trend_ssh_2, base_ssh_2			= polyfit(np.arange(len(time_year)), ssh_all_2[model_i], 1) * 10.0
	trend_ssh_dynamic_1, base_ssh_dynamic_1	= polyfit(np.arange(len(time_year)), ssh_total_1[model_i], 1) * 10.0
	trend_ssh_dynamic_2, base_ssh_dynamic_2	= polyfit(np.arange(len(time_year)), ssh_total_2[model_i], 1) * 10.0
	trend_ssh_global, base_ssh_global	= polyfit(np.arange(len(time_year)), ssh_global_all[model_i], 1) * 10.0

	ssh_trend_all_1[model_i]		= trend_ssh_1
	ssh_trend_all_2[model_i]		= trend_ssh_2
	ssh_dynamic_trend_all_1[model_i]	= trend_ssh_dynamic_1
	ssh_dynamic_trend_all_2[model_i]	= trend_ssh_dynamic_2
	ssh_global_trend[model_i]	        = trend_ssh_global

	print models[model_i],':', np.round(trend_ssh_global, 1), ssh_dynamic_trend_all_1[model_i] / trend_ssh_global, ssh_dynamic_trend_all_2[model_i] / trend_ssh_global

#Determine the trend for the HR-CESM and LR-CESM
trend_ssh_1_cesm, base_ssh_1_cesm				= polyfit(np.arange(len(time_year)), ssh_year_1_cesm, 1) * 10.0
trend_ssh_2_cesm, base_ssh_2_cesm				= polyfit(np.arange(len(time_year)), ssh_year_2_cesm, 1) * 10.0
trend_ssh_1_cesm_low, base_ssh_1_cesm_low			= polyfit(np.arange(len(time_year)), ssh_year_1_cesm_low, 1) * 10.0
trend_ssh_2_cesm_low, base_ssh_2_cesm_low			= polyfit(np.arange(len(time_year)), ssh_year_2_cesm_low, 1) * 10.0
trend_ssh_dynamic_1_cesm, base_ssh_dynamic_1_cesm		= polyfit(np.arange(len(time_year)), ssh_total_year_1_cesm, 1) * 10.0
trend_ssh_dynamic_2_cesm, base_ssh_dynamic_2_cesm		= polyfit(np.arange(len(time_year)), ssh_total_year_2_cesm, 1) * 10.0
trend_ssh_dynamic_1_cesm_low, base_ssh_dynamic_1_cesm_low	= polyfit(np.arange(len(time_year)), ssh_total_year_1_cesm_low, 1) * 10.0
trend_ssh_dynamic_2_cesm_low, base_ssh_dynamic_2_cesm_low	= polyfit(np.arange(len(time_year)), ssh_total_year_2_cesm_low, 1) * 10.0

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(np.arange(-5, 5.1, 1), np.arange(-5, 5.1, 1), '--k', linewidth = 2.0, alpha = 0.5)

ax.scatter(trend_ssh_1_cesm, trend_ssh_2_cesm, marker = 'o', color = 'r' , s = 75, label = 'HR-CESM')
ax.scatter(trend_ssh_1_cesm_low, trend_ssh_2_cesm_low, marker = 'o', color = 'b' , s = 75, label = 'LR-CESM')
ax.errorbar(np.mean(ssh_trend_all_1), np.mean(ssh_trend_all_2), yerr = np.std(ssh_trend_all_2), xerr = np.std(ssh_trend_all_1), color = 'k') 
ax.scatter(np.mean(ssh_trend_all_1), np.mean(ssh_trend_all_2), marker = 'o', color = 'k' , s = 75, label = 'CMIP6 mean + std')

ax.scatter(ssh_trend_all_1[0], ssh_trend_all_2[0], marker = 'p', color = 'k' , s = 75, label = models[0])
ax.scatter(ssh_trend_all_1[1], ssh_trend_all_2[1], marker = 'x', color = 'k' , s = 75, label = models[1])
ax.scatter(ssh_trend_all_1[2], ssh_trend_all_2[2], marker = '<', color = 'k' , s = 75, label = models[2])
ax.scatter(ssh_trend_all_1[3], ssh_trend_all_2[3], marker = 's', color = 'k' , s = 75, label = models[3])
ax.scatter(ssh_trend_all_1[4], ssh_trend_all_2[4], marker = 'D', color = 'k' , s = 75, label = models[4])
ax.scatter(ssh_trend_all_1[5], ssh_trend_all_2[5], marker = '*', color = 'k' , s = 75, label = models[5])
ax.scatter(ssh_trend_all_1[6], ssh_trend_all_2[6], marker = '>', color = 'k' , s = 75, label = models[6])
ax.scatter(ssh_trend_all_1[7], ssh_trend_all_2[7], marker = 'h', color = 'k' , s = 75, label = models[7])
ax.scatter(ssh_trend_all_1[8], ssh_trend_all_2[8], marker = 'v', color = 'k' , s = 75, label = models[8])
ax.scatter(ssh_trend_all_1[9], ssh_trend_all_2[9], marker = 'd', color = 'k' , s = 75, label = models[9])
ax.scatter(ssh_trend_all_1[10], ssh_trend_all_2[10], marker = '+', color = 'k' , s = 75, label = models[10])
ax.scatter(ssh_trend_all_1[11], ssh_trend_all_2[11], marker = '1', color = 'k' , s = 75, label = models[11])
ax.scatter(ssh_trend_all_1[12], ssh_trend_all_2[12], marker = '2', color = 'k' , s = 75, label = models[12])
ax.scatter(ssh_trend_all_1[13], ssh_trend_all_2[13], marker = '3', color = 'k' , s = 75, label = models[13])
ax.scatter(ssh_trend_all_1[14], ssh_trend_all_2[14], marker = '4', color = 'k' , s = 75, label = models[14])

ax.set_xlabel('$\eta_M$ trend region 1 (mm year$^{-1}$)')
ax.set_ylabel('$\eta_M$ trend region 2 (mm year$^{-1}$)')

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)

ax.grid()


left, width = .01, .9
bottom, height = 0.01, .9
right = left + width
top = bottom + height

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, frameon = False, prop={'size': 9.6})

ax.set_title('a) $\eta_M$ trends')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.axhline(y = 1.0, xmin = 0.275, xmax = 5, linestyle = '--', color = 'k', linewidth = 2.0, alpha = 0.5)
ax.axvline(x = 1.0, linestyle = '--', color = 'k', linewidth = 2.0, alpha = 0.5)

#Normalise the CMIP6 trend to the steric global trend
ssh_dynamic_trend_all_1 	= ssh_dynamic_trend_all_1 / ssh_global_trend
ssh_dynamic_trend_all_2 	= ssh_dynamic_trend_all_2 / ssh_global_trend

print 'CMIP6 region 1:', np.mean(ssh_dynamic_trend_all_1), np.std(ssh_dynamic_trend_all_1)
print 'CMIP6 region 2:', np.mean(ssh_dynamic_trend_all_2), np.std(ssh_dynamic_trend_all_2)

ax.scatter(trend_ssh_aviso_1 / 3.0, trend_ssh_aviso_2 / 3.0, marker = 'o', color = 'c' , s = 75, label = 'AVISO')
ax.scatter(trend_ssh_dynamic_1_cesm / trend_ssh_global_cesm, trend_ssh_dynamic_2_cesm / trend_ssh_global_cesm, marker = 'o', color = 'r' , s = 75, label = 'HR-CESM')
ax.scatter(trend_ssh_dynamic_1_cesm_low / trend_ssh_global_cesm_low, trend_ssh_dynamic_2_cesm_low / trend_ssh_global_cesm_low, marker = 'o', color = 'b' , s = 75, label = 'LR-CESM')
ax.errorbar(np.mean(ssh_dynamic_trend_all_1), np.mean(ssh_dynamic_trend_all_2), yerr = np.std(ssh_dynamic_trend_all_2), xerr = np.std(ssh_dynamic_trend_all_1), color = 'k') 
ax.scatter(np.mean(ssh_dynamic_trend_all_1), np.mean(ssh_dynamic_trend_all_2), marker = 'o', color = 'k' , s = 75, label = 'CMIP6 mean + std')

ax.scatter(ssh_dynamic_trend_all_1[0], ssh_dynamic_trend_all_2[0], marker = 'p', color = 'k' , s = 75, label = models[0])
ax.scatter(ssh_dynamic_trend_all_1[1], ssh_dynamic_trend_all_2[1], marker = 'x', color = 'k' , s = 75, label = models[1])
ax.scatter(ssh_dynamic_trend_all_1[2], ssh_dynamic_trend_all_2[2], marker = '<', color = 'k' , s = 75, label = models[2])
ax.scatter(ssh_dynamic_trend_all_1[3], ssh_dynamic_trend_all_2[3], marker = 's', color = 'k' , s = 75, label = models[3])
ax.scatter(ssh_dynamic_trend_all_1[4], ssh_dynamic_trend_all_2[4], marker = 'D', color = 'k' , s = 75, label = models[4])
ax.scatter(ssh_dynamic_trend_all_1[5], ssh_dynamic_trend_all_2[5], marker = '*', color = 'k' , s = 75, label = models[5])
ax.scatter(ssh_dynamic_trend_all_1[6], ssh_dynamic_trend_all_2[6], marker = '>', color = 'k' , s = 75, label = models[6])
ax.scatter(ssh_dynamic_trend_all_1[7], ssh_dynamic_trend_all_2[7], marker = 'h', color = 'k' , s = 75, label = models[7])
ax.scatter(ssh_dynamic_trend_all_1[8], ssh_dynamic_trend_all_2[8], marker = 'v', color = 'k' , s = 75, label = models[8])
ax.scatter(ssh_dynamic_trend_all_1[9], ssh_dynamic_trend_all_2[9], marker = 'd', color = 'k' , s = 75, label = models[9])
ax.scatter(ssh_dynamic_trend_all_1[10], ssh_dynamic_trend_all_2[10], marker = '+', color = 'k' , s = 75, label = models[10])
ax.scatter(ssh_dynamic_trend_all_1[11], ssh_dynamic_trend_all_2[11], marker = '1', color = 'k' , s = 75, label = models[11])
ax.scatter(ssh_dynamic_trend_all_1[12], ssh_dynamic_trend_all_2[12], marker = '2', color = 'k' , s = 75, label = models[12])
ax.scatter(ssh_dynamic_trend_all_1[13], ssh_dynamic_trend_all_2[13], marker = '3', color = 'k' , s = 75, label = models[13])
ax.scatter(ssh_dynamic_trend_all_1[14], ssh_dynamic_trend_all_2[14], marker = '4', color = 'k' , s = 75, label = models[14])

ax.set_xlabel('Normalised DSL trend region 1')
ax.set_ylabel('Normalised DSL trend region 2')

ax.set_xlim(0, 2.0)
ax.set_ylim(0, 2.0)

ax.grid()

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, frameon = False, prop={'size': 9.6})

x_square	= np.asarray([0.0, 0.2, 0.2, 0.0, 0.0])
y_square	= np.asarray([0.0, 0.0, 0.18, 0.18, 0.0])

ax.plot(x_square + 0.02, y_square + 0.02, '-k', linewidth = 2.0)
ax.fill(x_square + 0.02, y_square + 0.02, 'b', alpha = 0.3)
ax.text(0.1 + 0.02, 0.09 + 0.02, '1', horizontalalignment='center', verticalalignment='center', fontsize = 16)
ax.plot(x_square + 0.4, y_square + 0.20, '-k', linewidth = 2.0)
ax.fill(x_square + 0.4, y_square + 0.20, 'b', alpha = 0.3)
ax.text(0.1 + 0.4, 0.09 + 0.20, '2', horizontalalignment='center', verticalalignment='center', fontsize = 16)

ax.plot(x_square + 1.4, y_square + 0.02, '-k', linewidth = 2.0)
ax.fill(x_square + 1.4, y_square + 0.02, 'r', alpha = 0.3)
ax.text(0.1 + 1.4, 0.09 + 0.02, '1', horizontalalignment='center', verticalalignment='center', fontsize = 16)
ax.plot(x_square + 1.78, y_square + 0.20, '-k', linewidth = 2.0)
ax.fill(x_square + 1.78, y_square + 0.20, 'b', alpha = 0.3)
ax.text(0.1 + 1.78, 0.09 + 0.20, '2', horizontalalignment='center', verticalalignment='center', fontsize = 16)

ax.plot(x_square + 1.4, y_square + 1.62, '-k', linewidth = 2.0)
ax.fill(x_square + 1.4, y_square + 1.62, 'r', alpha = 0.3)
ax.text(0.1 + 1.4, 0.09 + 1.62, '1', horizontalalignment='center', verticalalignment='center', fontsize = 16)
ax.plot(x_square + 1.78, y_square + 1.8, '-k', linewidth = 2.0)
ax.fill(x_square + 1.78, y_square + 1.8, 'r', alpha = 0.3)
ax.text(0.1 + 1.78, 0.09 + 1.8, '2', horizontalalignment='center', verticalalignment='center', fontsize = 16)

ax.set_title('b) Normalised DSL trends')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CMIP6_graph	= ax.plot_date(time_year, np.mean(ssh_global_all, axis = 0), '-k', linewidth = 2.0, label = 'CMIP6')
ax.fill_between(time_year, np.min(ssh_global_all, axis = 0), np.max(ssh_global_all, axis = 0), color = 'k', alpha = 0.25)

HR_CESM_graph	= ax.plot_date(time_year, ssh_global_cesm_year, '-r', linewidth = 2.0, label = 'HR-CESM')
LR_CESM_graph	= ax.plot_date(time_year, ssh_global_cesm_low_year, '-b', linewidth = 2.0, label = 'LR-CESM')

ax.set_xlabel('Model year')
ax.set_ylabel('Sea level anomaly (cm)')
ax.set_ylim(-1, 25)
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
ax.legend(graphs, legend_labels, loc='upper left',
 		  ncol=1, fancybox=True, shadow=False, numpoints = 1)

ax.set_title('a) Global steric contribution')

ax2 = fig.add_axes([0.68, 0.10, 0.2, 0.2])

ax2.set_ylim(-0.1, 1.1)
ax2.set_xlim(0, 2.6)
ax2.axis('off')

x_legend	= np.arange(1, 2.51, 0.1)
ax2.fill_between(x_legend, 0.20, 0.80, color = 'k', alpha = 0.25)
ax2.plot(x_legend, 0.5 + np.zeros(len(x_legend)), linestyle = '-', color = 'k', linewidth = 3.0)


ax2.text(0.2, 0.20,'Min', color ='k',fontsize=15,ha='right',va='center')
ax2.plot([0.22, 1], [0.20, 0.20], '--k', linewidth = 0.5)

ax2.text(0.6, 0.5,'Mean', color ='k',fontsize=15,ha='right',va='center')
ax2.plot([0.62, 1], [0.5, 0.5], '--k', linewidth = 0.5)

ax2.text(0.2, 0.80,'Max', color ='k',fontsize=15,ha='right',va='center')
ax2.plot([0.22, 1], [0.80, 0.80], '--k', linewidth = 0.5)

show()