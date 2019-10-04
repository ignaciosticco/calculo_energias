'''
Suma los trabajos de las fuerzas asociadas a la interaccion con las paredes. 
Loop sobre configuraciones para distintas densidades gloables y calcula las energias
y densidades locales en cada caso.
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
	#list_dens_global =  [0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9]
	vd = 1
	clase = "knE5"
	file_constname = "constantes.dat"

	list_dens_local = []
	list_std_density = []
	list_wall_wfg = []
	list_wall_wfs = []
	list_wall_wfc = []
	for density in list_dens_global:
		name_input = "config_density{}_width22_knE5".format(density)
		name_output = "energias_wall_corridor_density{}_{}".format(density,clase)
		os.system("./calcula_energias_rectangulo_wall.exe {} {} {} {}".format(name_input,name_output,vd,file_constname))

		sum_energies =  suma_energias(name_output)
		list_dens_local += [sum_energies[0]]
		list_std_density += [sum_energies[1]]
		list_wall_wfg += [sum_energies[2]]
		list_wall_wfs += [sum_energies[3]]
		list_wall_wfc += [sum_energies[4]]


	df = pd.DataFrame()
	df['list_dens_local']  = np.round(list_dens_local,2) # Es la densidad local medida
	df['list_wfg_wall']  = np.round(list_wall_wfg,2)
	df['list_wfs_wall']  = np.round(list_wall_wfs,2)
	df['list_wfc_wall']  = np.round(list_wall_wfc,2)
	df.to_csv('sum_wallenergies_{}.txt'.format(clase), index=None, sep='\t')

def suma_energias(data_energias):
	'''
	Suma las energias de todas las particulas. 
	'''
	#Primero extrae la densidad y su desviacion standard del header. 
	df_header = pd.read_csv("{}".format(data_energias), nrows=1,delimiter='\t\t',engine="python")
	mean_density = df_header.mean_density[0]
	std_density = df_header.std_density[0]

	df = pd.read_csv("{}".format(data_energias),skiprows=3,delimiter='\t\t',engine="python") 
	sum_work_fgranular=df["work_fgranular_wall"].sum()
	sum_work_fsocial=df["work_fsocial_wall"].sum()
	sum_work_fcompresion=df["work_fcompresion_wall"].sum()
	
	return mean_density,std_density,sum_work_fgranular,sum_work_fsocial,sum_work_fcompresion


if __name__=='__main__':
	main()
