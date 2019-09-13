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

	df = pd.read_csv('sum_energies.txt',delimiter = '\t') 
	vdi = 0.5
	vdf = 6
	list_vd =  list(np.linspace(vdi,vdf,12))
	for v in list_vd:
		if v%1==0:
			list_vd[list_vd.index(v)]=int(v)

	list_wfg=df['list_wfg'].tolist()
	list_wfd=df['list_wfd'].tolist()
	list_wfs=df['list_wfs'].tolist()
	list_wfc=df['list_wfc'].tolist()
	plot_energies(list_vd,list_wfg,list_wfd,list_wfs,list_wfc)


def plot_energies(list_vd,list_wfg,list_wfd,list_wfs,list_wfc):
	
	plt.plot(list_vd,list_wfd,'b-+',mew=0.7,linewidth = '1',markersize=5,label='Desire') 
	plt.plot(list_vd,list_wfc,'r-x',mew=0.7,linewidth = '1',markersize=4,label='Granular potential') 
	plt.plot(list_vd,list_wfg,'c-d',mew=0.7,linewidth = '1',markerfacecolor='none',markersize=4,markeredgecolor='c',label='Granular Tang. Dissip.') 
	plt.plot(list_vd,list_wfs,'v',mec='mediumvioletred',mew=0.7,markerfacecolor='none',markersize=4,zorder=3,label='Social') 
	plt.plot(list_vd,list_wfs,'-',color='mediumvioletred',linewidth = '1') 


	pylab.grid(linewidth=0.3)
	pylab.xlabel('Desire velocity~(m s$^{-1}$)')
	pylab.ylabel('Work~(J)')
	plt.axis([0, 6, -3e6, 5e6])
	plt.ticklabel_format(style='sci', axis='y')
	lgd=plt.legend(numpoints=1,handlelength=0.8) 
	plt.legend(frameon=False,loc='best',labelspacing=-0.1,borderpad=0.3,handletextpad=0.5,fontsize=6,numpoints=1) 
	pylab.savefig('work_bulk.png', format='png', dpi=300, bbox_inches='tight')


if __name__=='__main__':
	main()