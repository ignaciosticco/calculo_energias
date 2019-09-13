'''
WARNING: Este testing solo funciona con los parametros originales de SFM
Con un corridor como el de Jamaraat (vd = 1 XT = 28 YD = 0.0 YUP = 22). 
<<<Ver los parametros abajo de todo>>>
Testea el codigo que calcula las energias de todas las fuerzas. 
Toma como input configuraciones de prueba (ficticias), calcula las energias 
correspondientes y checkea si el calculo de las energias es correcto. 
'''
import sys
import os
import unittest
import pandas as pd

vd=1
param_filename="constantes_originalsfm_anillo.dat" 

class testeos(unittest.TestCase):

	def test_energia_cinetica(self):
		name_input = "config_test_ecin.txt"
		name_output = "out__test_ecin.txt"
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,vd,param_filename))
		df=pd.read_csv("{}".format(name_output),delimiter='\t\t',skiprows=3,engine='python')
		self.assertEqual(float(df['avg_kinetic_energy'][0]),35.00)
		self.assertEqual(float(df['avg_kinetic_energy'][1]),140.00)


	def test_work_fgranular(self):
		name_input = "config_test_wfg.txt"
		name_output = "out__test_wfg.txt"
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,vd,param_filename))
		df=pd.read_csv("{}".format(name_output),delimiter='\t\t',skiprows=3,engine='python')
		self.assertEqual(float(df['work_fgranular'][0]),-21600)
		self.assertEqual(float(df['work_fgranular'][1]),0.00)
	

	
	def test_work_fdesired1(self):
		
				
		name_input = "config_test_wfd1.txt"
		name_output = "out__test_wfd1.txt"
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,vd,param_filename))
		df=pd.read_csv("{}".format(name_output),delimiter='\t\t',skiprows=3,engine='python')		
		self.assertEqual(float(df['work_fdesired'][0]),-420.0)
		self.assertEqual(float(df['work_fdesired'][1]),-840.0)

	
	def test_work_fdesired2(self):
		
		#Testea las particulas que estan arriba de YTARGETUP
		# y abajo de YTARGETDOWN
		name_input = "config_test_wfd2.txt"
		name_output = "out__test_wfd2.txt"
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,vd,param_filename))
		df=pd.read_csv("{}".format(name_output),delimiter='\t\t',skiprows=3,engine='python')			
		self.assertEqual(int(df['work_fdesired'][0]),-472) # Test down 
		self.assertEqual(int(df['work_fdesired'][1]),-18) # Test up
	

	
	def test_work_fsocial(self):
		name_input = "config_test_wfs.txt"
		name_output = "out__test_wfs.txt"
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,vd,param_filename))
		df=pd.read_csv("{}".format(name_output),delimiter='\t\t',skiprows=3,engine='python')			

		self.assertEqual(int(df['work_fsocial'][0]),-4946)
		self.assertEqual(int(df['work_fsocial'][1]),-9892)

	def test_work_fcompresion(self):
		name_input = "config_test_wfc.txt"
		name_output = "out__test_wfc.txt"
		os.system("./calcula_energias_anillo.exe {} {} {} {}".format(name_input,name_output,vd,param_filename))
		df=pd.read_csv("{}".format(name_output),delimiter='\t\t',skiprows=3,engine='python')			

		self.assertEqual(int(df['work_fcompresion'][0]),-7200)
		self.assertEqual(int(df['work_fcompresion'][1]),-14400)


if __name__=='__main__':

	# Warning: primero elimina a todos los archivos que comiencen con "out__"
	if(os.system("rm out__*")!=0):
		print("Devolvio este error porque no encontro ningun archivo 'out__' ")
	unittest.main()


'''
Estos son los parametros para los cuales se penso este testing

CANTATOMS_MAX=5
TI=0.0
TF=1.5
XCENTER=0.0
YCENTER=0.0
L_12=0.0
L_23=100.0

KAPPA=240000.0
KCOMP=120000.0
XTARGET=28.0
YTARGETUP=22.00
YTARGETDOWN=0.0

INTEGRATION_STEP=0.0001
MASS=70.0
A=2000.0
B=0.08
MOT=140.0
DIAM=0.46
TIMESTEP=0.5
RCUT2=0.7744
'''