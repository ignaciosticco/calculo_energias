'''
Loop sobre configuraciones para distintas densidades gloables y calcula las energias
y densidades locales en cada caso 
'''
import sys
import os
import pandas as pd
import pylab
import numpy as np
import math


def main():

	# WARNING: elimina todos los archivos que comiencen con "energias_bulk"
	os.system("rm energias_* ")

	# Escribir lista de densidades gloabales correspondientes a las config. 
	list_dens_global =  [1,2,3,4,5,6,7,8,9]  
	vd = 1
	clase = "test"  
	file_constname = "constantes.dat"

	list_dens_local = []
	list_std_density = []
	list_ecin = []
	list_wfg = []
	list_wfd = []
	list_wfs = []
	list_wfc = []
	for density in list_dens_global:
		name_input = "config_density{}_width22".format(density)
		name_output = "energias_corridor_density{}_{}".format(density,clase)
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,vd,file_constname))
		sum_energies =  suma_energias(name_output)
		list_dens_local += [sum_energies[0]]
		list_std_density += [sum_energies[1]]
		list_ecin += [sum_energies[2]]
		list_wfg += [sum_energies[3]]
		list_wfd += [sum_energies[4]]
		list_wfs += [sum_energies[5]]
		list_wfc += [sum_energies[6]]
	
	df = pd.DataFrame()
	df['list_dens_local']  = np.round(list_dens_local,2) # Es la densidad local medida
	df['list_ecin']  = np.round(list_ecin,2)
	df['list_wfg']  = np.round(list_wfg,2)
	df['list_wfd']  = np.round(list_wfd,2)
	df['list_wfs']  = np.round(list_wfs,2)
	df['list_wfc']  = np.round(list_wfc,2)
	df.to_csv('sum_energies_{}.txt'.format(clase), index=None, sep='\t')

def suma_energias(data_energias):
	'''
	Suma las energias de todas las particulas. 
	'''
	#Primero extrae la densidad y su desviacion standard del header. 
	df_header = pd.read_csv("{}".format(data_energias), nrows=1,delimiter='\t\t',engine="python")
	mean_density = df_header.mean_density[0]
	std_density = df_header.std_density[0]

	df = pd.read_csv("{}".format(data_energias),skiprows=3,delimiter='\t\t',engine="python") 
	sum_avg_kinetic_energy=df['avg_kinetic_energy'].sum()
	sum_work_fgranular=df["work_fgranular"].sum()
	sum_work_fdesired=df["work_fdesired"].sum()
	sum_work_fsocial=df["work_fsocial"].sum()
	sum_work_fcompresion=df["work_fcompresion"].sum()
	return mean_density,std_density,sum_avg_kinetic_energy,sum_work_fgranular,sum_work_fdesired,sum_work_fsocial,sum_work_fcompresion


if __name__=='__main__':
	main()
