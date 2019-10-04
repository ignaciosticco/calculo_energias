'''
Grafica todos los graficos de energia del paper //doi.org/10.1016/j.physa.2007.06.033
INPUT: archivo con densidades y la energia de cada tipo para cada densidad.
La energia de cada tipo es la suma sobre todos los individuos que se encuentran en 
una region determinada. 
NOTA: no hay normalizacion de la energia  
'''


import sys
import os
import pandas as pd
import pylab
import numpy as np
import matplotlib.pyplot as plt
import math

# a dos columnas: 3+3/8 (ancho, in)
# a una columna : 6.5   (ancho  in)

golden_mean = (math.sqrt(5)-1.0)/2.0        # Aesthetic ratio
fig_width = 3+3/8 			    # width  in inches
fig_height = fig_width*golden_mean          # height in inches
fig_size =  [fig_width,fig_height]

params = {'backend': 'ps',
          'axes.titlesize': 8,
          'axes.labelsize': 9,
          'axes.linewidth': 0.5, 
          'axes.grid': True,
          'axes.labelweight': 'normal',  
          'font.family': 'serif',
          'font.size': 8.0,
          'font.weight': 'normal',
          'text.color': 'black',
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': True,
          'legend.fontsize': 8,
          'figure.dpi': 300,
          'figure.figsize': fig_size,
          'savefig.dpi': 300,
         }

pylab.rcParams.update(params)

def main():

	df = pd.read_csv('sum_wallenergies_knE5.txt',delimiter = '\t') 
	clase = 'knE5_onlywalls'

	list_dens_local=df['list_dens_local'].tolist()
	list_wfg_wall=df['list_wfg_wall'].tolist()
	list_wfs_wall=df['list_wfs_wall'].tolist()
	list_wfc_wall=df['list_wfc_wall'].tolist()
	plot_energies_all_walls(list_dens_local,list_wfg_wall,list_wfs_wall,list_wfc_wall,clase)
	#plot_sum_energies(list_dens_local,list_wfg,list_wfd,list_wfs,list_wfc,clase)
	#plot_energies_ignoring(list_dens_local,list_wfg,list_wfd,list_wfs,list_wfc,"list_wfg",clase)
	#plot_energies_ignoring(list_dens_local,list_wfg,list_wfd,list_wfs,list_wfc,"list_wfc",clase)
	#plot_energies_ignoring(list_dens_local,list_wfg,list_wfd,list_wfs,list_wfc,"list_wfs",clase)


def plot_energies_all_walls(list_dens_local,list_wfg_wall,list_wfs_wall,list_wfc_wall,clase):

	plt.cla()	
	plt.plot(list_dens_local,list_wfc_wall,'r-x',mew=0.7,linewidth = '1',markersize=4,label='Granular potential') 
	plt.plot(list_dens_local,list_wfg_wall,'c-d',mew=0.7,linewidth = '1',markerfacecolor='none',markersize=4,markeredgecolor='c',label='Granular Tang. Dissip.') 
	plt.plot(list_dens_local,list_wfs_wall,'v',mec='mediumvioletred',mew=0.7,markerfacecolor='none',markersize=4,zorder=3,label='Social') 
	plt.plot(list_dens_local,list_wfs_wall,'-',color='mediumvioletred',linewidth = '1') 
	pylab.grid(linewidth=0.3)
	pylab.xlabel('Density~(p m$^{-2}$)')
	pylab.ylabel('Work~(J)')
	#plt.axis([0, 9.5, -1e7, 1.2e7])
	#plt.ticklabel_format(style='sci', axis='y')
	plt.yscale('symlog')
	lgd=plt.legend(numpoints=1,handlelength=0.8) 
	plt.legend(frameon=False,loc='best',labelspacing=-0.1,borderpad=0.3,handletextpad=0.5,fontsize=6,numpoints=1) 
	pylab.savefig('all_works_{}.png'.format(clase), format='png', dpi=300, bbox_inches='tight')


def plot_energies_ignoring(list_dens_local,list_wfg,list_wfd,list_wfs,list_wfc,ignore,clase):
	'''
	Grafica la energia total ignorando la contribucion de un tipo de trabajo 
	Esta contribucion se pasa con la variable ignore. 
	'''
	
	plt.cla()
	total_work =  np.array(list_wfg)+np.array(list_wfd)+np.array(list_wfs)+np.array(list_wfc)

	if ignore=="list_wfg":
		work = total_work-np.array(list_wfg)
		simbol ='c-d'
		title = "Friction" 
		color = 'c'
	elif ignore=="list_wfc":
		work = total_work-np.array(list_wfc)
		simbol='r-x'
		title = "Compresion" 
		color = 'r'
	elif ignore=="list_wfs":
		work = total_work-np.array(list_wfs)
		simbol='-v'
		title = "Social" 
		color = 'mediumvioletred'
	else:
		print("ERROR: Revisar el campo 'ignore' ")

	plt.plot(list_dens_local,work,simbol,color=color,mfc='none',mew=0.7,linewidth = '1',markersize=5) 	
	pylab.grid(linewidth=0.3)
	pylab.xlabel('Density~(p m$^{-2}$)')
	pylab.ylabel('Work~(J)')
	pylab.title("Ignoring {}".format(title))
	plt.yscale('symlog')
	#plt.axis([0, 9.5, -1e7, 1.2e7])
	#plt.ticklabel_format(style='sci', axis='y')
	#lgd=plt.legend(numpoints=1,handlelength=0.8) 
	#plt.legend(frameon=False,loc='best',labelspacing=-0.1,borderpad=0.3,handletextpad=0.5,fontsize=6,numpoints=1) 
	pylab.savefig('work_ignore_{}_{}.png'.format(ignore,clase), format='png', dpi=300, bbox_inches='tight')


def plot_sum_energies(list_dens_local,list_wfg,list_wfd,list_wfs,list_wfc,clase):
	'''
	Suma todas las energias y la grafica vs la densidad
	'''

	plt.cla()
	total_work =  np.array(list_wfg)+np.array(list_wfd)+np.array(list_wfs)+np.array(list_wfc)

	plt.plot(list_dens_local,total_work,'k-o',mew=0.7,linewidth = '1',markersize=5) 	
	pylab.grid(linewidth=0.3)
	pylab.xlabel('Density~(p m$^{-2}$)')
	pylab.ylabel('Total Work~(J)')
	plt.yscale('symlog')
	#plt.axis([0, 9.5, -1e7, 1.2e7])
	#plt.ticklabel_format(style='sci', axis='y')
	pylab.savefig('total_work_{}.png'.format(clase), format='png', dpi=300, bbox_inches='tight')


if __name__=='__main__':
	main()