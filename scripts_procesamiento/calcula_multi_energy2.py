'''
Loop sobre configuraciones para distintas densidades gloables y calcula 
el promedio de las energias y densidades locales en cada caso.  
'''
import sys
import os
import pandas as pd
import numpy as np
import math


def main():

	# WARNING: elimina todos los archivos que comiencen con "energias_bulk"
	os.system("rm energias_* ")

	# Escribir lista de densidades gloabales correspondientes a las config. 
	list_dens_global =  [1,2]  
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
		name_input = "config_bottleneck_verific_vd{}_N225".format(density)
		name_output = "energias_corridor_density{}_{}".format(density,clase)
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,density,file_constname))
		promedia_energies =  promedia_energias(name_output)
		list_dens_local += [promedia_energies[0]]
		list_std_density += [promedia_energies[1]]
		list_ecin += [promedia_energies[2]]
		list_wfg += [promedia_energies[3]]
		list_wfd += [promedia_energies[4]]
		list_wfs += [promedia_energies[5]]
		list_wfc += [promedia_energies[6]]
	
	df = pd.DataFrame()
	df['list_dens_local']  = np.round(list_dens_local,2) # Es la densidad local medida
	df['list_ecin']  = np.round(list_ecin,2)
	df['list_wfg']  = np.round(list_wfg,2)
	df['list_wfd']  = np.round(list_wfd,2)
	df['list_wfs']  = np.round(list_wfs,2)
	df['list_wfc']  = np.round(list_wfc,2)
	df.to_csv('promedia_energies_{}.txt'.format(clase), index=None, sep='\t')

def promedia_energias(data_energias):
	'''
	Promedia las energias de todas las particulas. 
	'''
	#Primero extrae la densidad y su desviacion standard del header. 
	df_header = pd.read_csv("{}".format(data_energias), nrows=1,delimiter='\t\t',engine="python")
	mean_density = df_header.mean_density[0]
	std_density = df_header.std_density[0]

	df = pd.read_csv("{}".format(data_energias),skiprows=3,delimiter='\t\t',engine="python") 

	mean_avg_kinetic_energy=df['avg_kinetic_energy'].loc[df['avg_kinetic_energy']!=0.0].mean()
	mean_work_fgranular=df["work_fgranular"].loc[df["work_fgranular"]!=0.0].mean()
	mean_work_fdesired=df["work_fdesired"].loc[df["work_fdesired"]!=0.0].mean()
	mean_work_fsocial=df["work_fsocial"].loc[df["work_fsocial"]!=0.0].mean()
	mean_work_fcompresion=df["work_fcompresion"].loc[df["work_fcompresion"]!=0.0].mean()


	return mean_density,std_density,mean_avg_kinetic_energy,mean_work_fgranular,mean_work_fdesired,mean_work_fsocial,mean_work_fcompresion


if __name__=='__main__':
	main()
