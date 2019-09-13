/*
Input: Lammps configuration dump with: id x y vx vy diameter
Output: For each pedestrians returns the energies along the trajectory
Calcula la energia de todas las particulas que estan en la region determinada
por la funcion esta_en_anillo. 
*/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <math.h>
#include <stdio.h>

using namespace std;

void calcula_energia_cinetica(vector<int> &vector_id,vector<double> &vector_vx,
	vector<double> &vector_vy,vector<double> &kinetic_energy,vector<double> & vector_x,vector<double> & vector_y);

void calcula_avg_energia_cinetica(vector<double> &avg_kinetic_energy,
	vector<double> &kinetic_energy,int iter);

void calcula_work_fgranular(vector<int> &vector_id,vector<double> &vector_x, 
	vector<double> &vector_y,vector<double> &vector_vx,vector<double> &vector_vy,vector<double> &work_fg);
void calcula_granular_force(vector<double> &vector_fg,vector<double> &vector_x,
 vector<double> &vector_y, vector<double> &vector_vx,vector<double> &vector_vy,int id);

void calcula_work_fdesired(vector<int> &vector_id,vector<double> &vector_x,
 vector<double> &vector_y,vector<double> &vector_vx,vector<double> &vector_vy,vector<double> &work_fd,double vd);
void calcula_desired_force(vector<double> &force_desired,
	double xi,double yi,double vxi,double vyi,double vd);

void calcula_work_fsocial(vector<int> &vector_id,vector<double> &vector_x,
 vector<double> &vector_y,vector<double> &vector_vx,vector<double> &vector_vy,vector<double> &work_fs);
void calcula_social_force(vector<double> &vector_fsocial,vector<double> &vector_x,
 vector<double> &vector_y, vector<double> &vector_vx,vector<double> &vector_vy,int i);

void calcula_work_fcompresion(vector<int> &vector_id,vector<double> &vector_x, 
	vector<double> &vector_y,vector<double> &vector_vx,vector<double> &vector_vy,vector<double> &work_fcompresion);
void calcula_compresion_force(vector<double> &vector_fcompresion,vector<double> &vector_x,
 vector<double> &vector_y, int id);

void escribir( vector<double> &observable1,vector<double> &observable2,vector<double> &observable3,
	vector<double> &observable4,vector<double> &observable5,string output_file);

bool esta_en_anillo(double x, double y);

void read_constants(string archivo_ctes);

int    CANTATOMS_MAX;
double TI,TF,XCENTER,YCENTER,L_12,L_23,KAPPA,KCOMP,XTARGET, YTARGETUP,YTARGETDOWN,INTEGRATION_STEP, MASS,A,B,MOT, DIAM,TIMESTEP,RCUT2;



int main(int argc, char const *argv[]){

	double x,y,vx,vy,id,time,diameter;
	int cantAtoms, n_time,i;
	string  line;
	int    iter=0;
	
	//////////// INPUTS ////////////////
	string archivoConfig=argv[1];
	string archivoOut=argv[2];
	string vds =  argv[3];
	double vd = atof(vds.c_str());
	////////////////////////////////////

	string archivo_ctes=argv[4];
	string coso;
	ifstream filectes(archivo_ctes.c_str());
	read_constants(archivo_ctes);

	vector<double> kinetic_energy(CANTATOMS_MAX+1, 0.0);
	vector<double> avg_kinetic_energy(CANTATOMS_MAX+1, 0.0);
	vector<double> work_fgranular(CANTATOMS_MAX+1, 0.0);
	vector<double> work_fdesired(CANTATOMS_MAX+1, 0.0);
	vector<double> work_fsocial(CANTATOMS_MAX+1, 0.0);
	vector<double> work_fcompresion(CANTATOMS_MAX+1, 0.0);

	ifstream fileIn(archivoConfig.c_str());
	while(fileIn.good()){
		getline(fileIn,line,'P');
		getline(fileIn,line,'I');
		n_time=atoi(line.c_str());
		time = n_time*INTEGRATION_STEP;
		getline(fileIn,line,'S');
		getline(fileIn,line,'I');
		cantAtoms=atoi(line.c_str());
		getline(fileIn,line,'y');
		getline(fileIn,line,'r');// ultima letra en tabla configuraciones es r (de diameter)
		getline(fileIn,line,' ');

		vector<int> vector_id;
		vector<double> vector_x;
		vector<double> vector_y;
		vector<double> vector_vx;
		vector<double> vector_vy;
		vector<double> vector_diameter;

		i = 0;
		while(i<cantAtoms){		
			getline(fileIn,line,' ');
			id=atoi(line.c_str());
			getline(fileIn,line,' ');
			x=atof(line.c_str());
			getline(fileIn,line,' ');
			y=atof(line.c_str());
			getline(fileIn,line,' ');
			vx=atof(line.c_str());
			getline(fileIn,line,' ');
			vy=atof(line.c_str());
			getline(fileIn,line,' ');
			diameter=atof(line.c_str());
			vector_id.push_back(id);
			vector_x.push_back(x);
			vector_y.push_back(y);			
			vector_vx.push_back(vx);
			vector_vy.push_back(vy);	
			vector_diameter.push_back(diameter);
			i++;		
		}

		//// Calculo de energias en Bulk ////////
		if (time>=TI && time<=TF && fileIn.good()){
			calcula_energia_cinetica(vector_id,vector_vx,vector_vy,kinetic_energy,vector_x, vector_y);
			calcula_work_fgranular(vector_id,vector_x, vector_y,vector_vx,vector_vy,work_fgranular);
			calcula_work_fdesired(vector_id,vector_x, vector_y,vector_vx,vector_vy,work_fdesired,vd);
			calcula_work_fsocial(vector_id,vector_x,vector_y,vector_vx,vector_vy,work_fsocial);
			calcula_work_fcompresion(vector_id,vector_x,vector_y,vector_vx,vector_vy,work_fcompresion);
			iter++;
		}
		//////////////////////////////////////////		
	}
	fileIn.close();
	calcula_avg_energia_cinetica(avg_kinetic_energy,kinetic_energy,iter);
	escribir(avg_kinetic_energy,work_fgranular,work_fdesired,work_fsocial,work_fcompresion,archivoOut);
}

void calcula_work_fcompresion(vector<int> &vector_id,vector<double> &vector_x,
 vector<double> &vector_y,vector<double> &vector_vx,vector<double> &vector_vy,
 vector<double> &work_fcompresion){

	int id;
	int n = vector_id.size();
	for (int i = 0; i < n; ++i){
		if (esta_en_anillo(vector_x[i],vector_y[i])){
			id = vector_id[i];
			vector<double> vector_fcompresion(2,0.0);
			calcula_compresion_force(vector_fcompresion,vector_x, vector_y,i);	
			work_fcompresion[id]+=(vector_vx[i]*vector_fcompresion[0]+vector_vy[i]*vector_fcompresion[1])*TIMESTEP;	
		}
	}
}

void calcula_compresion_force(vector<double> &vector_fcompresion,
	vector<double> &vector_x,vector<double> &vector_y,int i){
 
	double xi,yi,xj,yj,delx,dely,r,gpair,rsq,compresion_factor;
	double sum_rads = DIAM;

	xi = vector_x[i];
	yi = vector_y[i];
	int j = 0;
	while(j<(int)vector_x.size()){
		xj = vector_x[j];
		yj = vector_y[j];
		delx = xi-xj;
		dely = yi-yj;
	
		rsq = delx*delx+dely*dely;

		if (rsq<(sum_rads*sum_rads) && rsq>0.0 ){
			r = sqrt(rsq);	
			gpair = sum_rads - r;
			compresion_factor = KCOMP*gpair/r; 
	      	vector_fcompresion[0] += compresion_factor*delx; 
	      	vector_fcompresion[1] += compresion_factor*dely;
		}
		j++;
	}
}


void calcula_work_fsocial(vector<int> &vector_id,vector<double> &vector_x,
 vector<double> &vector_y,vector<double> &vector_vx,vector<double> &vector_vy,vector<double> &work_fsocial){

	int id;
	int n = vector_id.size();
	for (int i = 0; i < n; ++i){
		if (esta_en_anillo(vector_x[i],vector_y[i])){
			id = vector_id[i];
			vector<double> vector_fsocial(2,0.0);
			calcula_social_force(vector_fsocial,vector_x, vector_y, vector_vx,vector_vy, i);
			work_fsocial[id]+=(vector_vx[i]*vector_fsocial[0]+vector_vy[i]*vector_fsocial[1])*TIMESTEP;	
		}
	}
}


void calcula_social_force(vector<double> &vector_fsocial,vector<double> &vector_x,
 vector<double> &vector_y, vector<double> &vector_vx,vector<double> &vector_vy,int i){

	double xi,yi,xj,yj,delx,dely,r,rsq,fpair;
	double sum_rads = DIAM;

	xi = vector_x[i];
	yi = vector_y[i];
	int j = 0;
	while(j<(int)vector_x.size()){
		xj = vector_x[j];
		yj = vector_y[j];
		delx = xi-xj;
		dely = yi-yj;
	
		rsq = delx*delx+dely*dely;

		if (rsq<RCUT2 && rsq>0.0 ){
			r = sqrt(rsq);	
			fpair = A*exp((sum_rads-r)/B);
	      	vector_fsocial[0] += fpair*delx/r; 
	      	vector_fsocial[1] += fpair*dely/r;
		}
		j++;
	}
}

void calcula_work_fdesired(vector<int> &vector_id,vector<double> &vector_x,
 vector<double> &vector_y,vector<double> &vector_vx,vector<double> &vector_vy,
 vector<double> &work_fdesired,double vd){

	int id;
	double xi,yi,vxi,vyi;
	int n = vector_id.size();
	for (int i = 0; i < n; ++i){
		if (esta_en_anillo(vector_x[i],vector_y[i])){
			id = vector_id[i];
			xi = vector_x[i];
			yi = vector_y[i];;
			vxi = vector_vx[i];
			vyi = vector_vy[i];
			vector<double> force_desired(2,0.0); // vector Fd para el individuo i 
			calcula_desired_force(force_desired,xi,yi,vxi,vyi,vd);
			work_fdesired[id]+=(vector_vx[i]*force_desired[0]+vector_vy[i]*force_desired[1])*TIMESTEP;	
		}
	}
}


void calcula_desired_force(vector<double> &force_desired,double xi,
	double yi,double vxi,double vyi,double vd){
	/*
	Calcula la fuerza de deseo para un opening vertical
	Asume que los individuos van de Izquierda a derecha. 
	*/

	double dy,rsq,nx,ny,rinv;

  	force_desired[0] = 0.0;
  	force_desired[1] = 0.0;

	double dx = XTARGET - xi;
	if (dx> 0){       //significa que no paso por el opening
		if (yi >= YTARGETDOWN && yi <= YTARGETUP){  // Medio del opening
			dy = 0.0;
			rsq = dx*dx;
			nx = dx/sqrt(rsq);
			ny = 0.0;
			force_desired[0] = MOT*(vd*nx - vxi); 
			force_desired[1] = MOT*(vd*ny - vyi);  
			
		}
		else if (yi < YTARGETDOWN){ //Abajo del opening
			dy = YTARGETDOWN - yi;
			rsq = dx*dx + dy*dy;      
			rinv = 1.0/sqrt(rsq);
			nx = dx*rinv;
			ny = dy*rinv;
			force_desired[0] = MOT*(vd*nx-vxi);  
			force_desired[1] = MOT*(vd*ny-vyi);  
			
		}
		else {	// Arriba del opening
			dy = YTARGETUP - yi;
			rsq = dx*dx + dy*dy;
			rinv = 1.0/sqrt(rsq);
			nx = dx*rinv;
			ny = dy*rinv;
			force_desired[0] = MOT*(vd*nx-vxi);  
			force_desired[1] = MOT*(vd*ny-vyi); 
		}
	}

}


void calcula_work_fgranular(vector<int> &vector_id,vector<double> &vector_x, vector<double> &vector_y,vector<double> &vector_vx,vector<double> &vector_vy,vector<double> &work_fgranular){

	int id;
	int n = vector_id.size();
	for (int i = 0; i < n; ++i){
		if (esta_en_anillo(vector_x[i],vector_y[i])){
			id = vector_id[i];
			vector<double> vector_fg(2,0.0);
			calcula_granular_force(vector_fg,vector_x, vector_y, vector_vx,vector_vy, i);
			work_fgranular[id]+=(vector_vx[i]*vector_fg[0]+vector_vy[i]*vector_fg[1])*TIMESTEP;	
		}
	}
}


void calcula_granular_force(vector<double> &vector_fg,vector<double> &vector_x,
 vector<double> &vector_y, vector<double> &vector_vx,vector<double> &vector_vy,int i){
   
	double xi,yi,xj,yj,vxi,vyi,vxj,vyj,delx,dely,delvx,delvy,
	r,gpair,granular_factor,rsq;

	double sum_rads = DIAM;
	xi = vector_x[i];
	yi = vector_y[i];
	vxi = vector_vx[i];
	vyi = vector_vy[i];	
	int j = 0;
	while(j<(int)vector_x.size()){
		xj = vector_x[j];
		yj = vector_y[j];
		vxj = vector_vx[j];
		vyj = vector_vy[j];		
		delx = xi-xj;
		dely = yi-yj;
		delvx = vxi - vxj;
		delvy = vyi - vyj;		
		
		rsq = delx*delx+dely*dely;

		if (rsq<(sum_rads*sum_rads) && rsq>0.0 ){
			r = sqrt(rsq);	
			gpair = sum_rads - r;
			granular_factor = KAPPA*gpair*(dely*delvx-delx*delvy)/rsq; 
	      	vector_fg[0] +=- granular_factor*dely; 
	      	vector_fg[1] += granular_factor*delx;
		}
		j++;
	}
}

void calcula_energia_cinetica(vector<int> &vector_id,vector<double> &vector_vx,
	vector<double> &vector_vy,vector<double> &kinetic_energy,vector<double> & vector_x,vector<double> & vector_y){
	/*
	Acumula la energia cinetica de cada particula (suma Ecin timestep a timestep). 
	*/

	int id;
	int n = vector_id.size();
	for (int i = 0; i < n; ++i){
		if (esta_en_anillo(vector_x[i],vector_y[i])){
			id = vector_id[i];
			kinetic_energy[id]+=0.5*MASS*(vector_vx[i]*vector_vx[i]+vector_vy[i]*vector_vy[i]);
		}
	}
}

void calcula_avg_energia_cinetica(vector<double> &avg_kinetic_energy,vector<double> &kinetic_energy,int iter){
	/*
	Divide los valores de la energia cinetica por la cantidad de veces que se midio la 
	energia cinetica (devuelve el promedio de la energia cinetica). 
	*/

	for (int i = 0; i <  (int) avg_kinetic_energy.size(); ++i){
		avg_kinetic_energy[i]=kinetic_energy[i]/((double)iter);
	}
}


void escribir( vector<double> &observable1,vector<double> &observable2,vector<double> &observable3,
	vector<double> &observable4,vector<double> &observable5,string output_file){

	char char_output_file[output_file.size() + 1];
	strcpy(char_output_file, output_file.c_str());	
	FILE *fp;
	fp=fopen(char_output_file,"a");
	fprintf(fp, "avg_kinetic_energy\t\twork_fgranular\t\twork_fdesired\t\twork_fsocial\t\twork_fcompresion\n");
	for (int i = 1; i < (int)observable1.size(); ++i){
		fprintf(fp,"%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n",observable1[i],observable2[i],observable3[i],
			observable4[i],observable5[i]);
	}
	fclose(fp);
}

bool esta_en_anillo(double x, double y){
	/*
	Devuelve True si la particula esta adentro del anillo  L_12<r<L_23
	Segun la notacion del paper doi.org/10.1016/j.physa.2007.06.033
	*/
	double r = sqrt((x-XCENTER)*(x-XCENTER)+(y-YCENTER)*(y-YCENTER));
	return (r<=L_23 && r>=L_12);
}

void read_constants(string archivo_ctes){

	string tmp;
	ifstream filectes(archivo_ctes.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	CANTATOMS_MAX=atoi(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	TI=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	TF=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	XCENTER=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	YCENTER=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	L_12=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	L_23=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	KAPPA=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	KCOMP=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	XTARGET=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	YTARGETUP=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	YTARGETDOWN=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	INTEGRATION_STEP=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	MASS=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	A=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	B=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	MOT=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	DIAM=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	TIMESTEP=atof(tmp.c_str());
	getline(filectes,tmp,'=');
	getline(filectes,tmp,'\n');
	RCUT2=atof(tmp.c_str());

}