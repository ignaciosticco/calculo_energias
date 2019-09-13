'''
Loop sobre configuraciones para distintas vd y calcula las energias en 
cada caso. 
'''
import sys
import os
import pandas as pd
import pylab
import numpy as np
import matplotlib.pyplot as plt
import math
#import numarray

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

	# WARNING: elimina todos los archivos que comiencen con "energias_bulk"
	os.system("rm energias_bulk_* ")
	file_constantes = "constantes.dat"
	vdi = 0.5
	vdf = 6
	list_vd =  list(np.linspace(vdi,vdf,12))
	for v in list_vd:
		if v%1==0:
			list_vd[list_vd.index(v)]=int(v)
	list_ecin = []
	list_wfg = []
	list_wfd = []
	list_wfs = []
	list_wfc = []
	for vd in list_vd:
		name_input = "config_bottleneck_verific_vd{}_N225".format(vd)
		name_output = "energias_bulk_225_vd{}".format(vd)
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,vd,file_constantes))
		sum_energies =  suma_energias(name_output)
		list_ecin += [sum_energies[0]]
		list_wfg += [sum_energies[1]]
		list_wfd += [sum_energies[2]]
		list_wfs += [sum_energies[3]]
		list_wfc += [sum_energies[4]]
	plot_energies(list_vd,list_wfg,list_wfd,list_wfs,list_wfc)
	
	df = pd.DataFrame()
	df['list_ecin']  = list_ecin
	df['list_wfg']  = list_wfg
	df['list_wfd']  = list_wfd
	df['list_wfs']  = list_wfs
	df['list_wfc']  = list_wfc
	df.to_csv('sum_energies.txt', index=None, sep='\t')

def suma_energias(data_energias):
	df = pd.read_csv("{}".format(data_energias),delimiter='\t\t',skiprows=3,engine="python") 
	sum_avg_kinetic_energy=df['avg_kinetic_energy'].sum()
	sum_work_fgranular=df["work_fgranular"].sum()
	sum_work_fdesired=df["work_fdesired"].sum()
	sum_work_fsocial=df["work_fsocial"].sum()
	sum_work_fcompresion=df["work_fcompresion"].sum()
	return sum_avg_kinetic_energy,sum_work_fgranular,sum_work_fdesired,sum_work_fsocial,sum_work_fcompresion


def plot_energies(list_vd,list_wfg,list_wfd,list_wfs,list_wfc):
	
	plt.plot(list_vd,list_wfd,'b-+',mew=0.7,linewidth = '1',markersize=5,label='Desire') 
	plt.plot(list_vd,list_wfc,'r-x',mew=0.7,linewidth = '1',markersize=4,label='Granular potential') 
	plt.plot(list_vd,list_wfg,'c-d',mew=0.7,linewidth = '1',markerfacecolor='none',markersize=4,markeredgecolor='c',label='Granular Tang. Dissip.') 
	plt.plot(list_vd,list_wfs,'v',mec='mediumvioletred',mew=0.7,markerfacecolor='none',markersize=4,zorder=3,label='Social') 
	plt.plot(list_vd,list_wfs,'-',color='mediumvioletred',linewidth = '1') 


	pylab.grid(linewidth=0.3)
	pylab.xlabel('Desire velocity~(m s$^{-1}$)')
	plt.axis([0, 6, -3e6, 5e6])
	plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	lgd=plt.legend(numpoints=1,handlelength=0.8) 
	plt.legend(frameon=False,loc='best',labelspacing=-0.1,borderpad=0.3,handletextpad=0.5,fontsize=6,numpoints=1) 
	pylab.savefig('work_bulk.png', format='png', dpi=300, bbox_inches='tight')


if __name__=='__main__':
	main()