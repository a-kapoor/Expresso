#include <stdlib.h> 
#include <time.h>  
#include <random>
#include <iostream>

using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;

using Vec_t = const ROOT::RVec<float>&;
using RVecF = RVec<float>;
using RVecI = RVec<int>;


RVecF GenerateDecayTime(RVecI Particle_pid, RVecF Particle_px, RVecF Particle_py, RVecF Particle_pz, RVecF Particle_energy) {
    int count_mesons=0;
	int l=0;
	for (size_t k = 0; k < Particle_pid.size(); k++){
	if (Particle_pid[k]==4900111){
	    count_mesons=count_mesons+1;
		}
    }

    RVecF Particle_DecayTime(count_mesons);


    for (size_t i = 0; i < Particle_pid.size(); i++){
		if (Particle_pid[i]==4900111){
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0.0, 1.0);
        float random=dis(gen);
		float t=-log(random);
		Particle_DecayTime[l]=t;
		l=l+1;
		}
	}
    return Particle_DecayTime;
}

RVecF Particle_DecayLength(RVecI Particle_pid, RVecF Particle_px, RVecF Particle_py, RVecF Particle_pz, RVecF Particle_energy, RVecF Particle_DecayTime){
	int l=0;
	int count_mesons=0;
	for (size_t k = 0; k < Particle_pid.size(); k++){
	if (Particle_pid[k]==4900111){
	    count_mesons=count_mesons+1;
		}
    }

    RVecF Particle_Decaylength(count_mesons);

	for (size_t i=0;i<Particle_pid.size();i++){
		if (Particle_pid[i]==4900111){
		float total_p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
		float beta=total_p/Particle_energy[i];
		float gamma=1/(sqrt(1-beta*beta));
		float t=Particle_DecayTime[l];
		float decaylength=beta*gamma*t;
		Particle_Decaylength[l]=decaylength;
		l=l+1;
		}
	}
	return  Particle_Decaylength;
}

RVecF Particle_Theta(RVecI Particle_pid, RVecF Particle_px, RVecF Particle_py, RVecF Particle_pz){
	float pi=3.1415926535;
	RVecF Particle_theta(Particle_pid.size());
	for (size_t i=0;i<Particle_pid.size();i++){
		float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
		float pxy=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2));
		
		if(Particle_pz[i]>=0){
			float theta=asin(pxy/p);
			Particle_theta[i]=theta;
		}
		if (Particle_pz[i]<0){
			float theta=pi-asin(pxy/p);
			Particle_theta[i]=theta;
		}
	}
	return Particle_theta;
}

RVecF Particle_Phi(RVecI Particle_pid, RVecF Particle_px, RVecF Particle_py, RVecF Particle_pz){
	float pi=3.1415926535;
	RVecF Particle_phi(Particle_pid.size());
	for (size_t i=0;i<Particle_pid.size();i++){
		float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
		float pxy=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2));
		if (Particle_px[i]>0 && Particle_py[i]>0){
			float phi=acos(Particle_px[i]/pxy);
			Particle_phi[i]=phi;
		}
		if (Particle_px[i]<0 && Particle_py[i]>0){
			float phi=pi/2+acos(Particle_py[i]/pxy);
			Particle_phi[i]=phi;
		}
		if (Particle_px[i]<0 && Particle_py[i]<0){
			float phi=pi/2+acos(Particle_py[i]/pxy);
			Particle_phi[i]=phi;
		}

		if (Particle_px[i]>0 && Particle_py[i]<0){
			float phi=2*pi-acos(Particle_px[i]/pxy);
			Particle_phi[i]=phi;
		}

		

	}
	return Particle_phi;
}

RVecF Particle_X_Meson(RVecI Particle_pid, RVecF Particle_DecayTime, RVecF Particle_theta, RVecF Particle_phi, RVecF Particle_Decaylength) {
    int count_mesons=0;
    int l=0;
	int j=0;
    for (size_t k = 0; k < Particle_pid.size(); k++){
	if (Particle_pid[k]==4900111){
	    count_mesons=count_mesons+1;
		}
    }
    RVecF Particle_x_meson(count_mesons*2);

    for (size_t i = 0; i < Particle_pid.size(); i++){
	if (Particle_pid[i]==4900111){
	float theta=Particle_theta[i];
	float phi=Particle_phi[i];
		float decay_x=Particle_Decaylength[l]*sin(theta)*cos(phi)*1000;
		Particle_x_meson[j]=decay_x;
		Particle_x_meson[j+1]=decay_x;
		l=l+1;
		j=j+2;
	    }
		}
    return Particle_x_meson;
}

RVecF Particle_X(RVecI Particle_pid, RVecF Particle_x_meson, RVecF Particle_status, RVecF Particle_x, RVecF Particle_y) {
    int count_muons=0;

    int norm=1;
	int l=0;

    RVecF N_Particle_x(Particle_pid.size());

    for (size_t i = 0; i < Particle_pid.size(); i++){
		N_Particle_x[i]=Particle_x[i];
		if (fabs(Particle_pid[i])==13 && Particle_status[i]==1 && sqrt(pow(Particle_x[i],2)+pow(Particle_y[i],2))>10){
			count_muons=count_muons+1;
		float x = Particle_x_meson[l];
		float x1= x * norm;
		N_Particle_x[i]=x1;
		l=l+1;
	    }
	}
	
    return N_Particle_x;
}

RVecF Particle_Y_Meson(RVecI Particle_pid, RVecF Particle_DecayTime, RVecF Particle_theta, RVecF Particle_phi, RVecF Particle_Decaylength) {
    int count_mesons=0;

    int l=0;
	int j=0;
    for (size_t k = 0; k < Particle_pid.size(); k++){
	if (Particle_pid[k]==4900111){
	    count_mesons=count_mesons+1;
		}
    }
    RVecF Particle_y_meson(count_mesons*2);

    for (size_t i = 0; i < Particle_pid.size(); i++){
	if (Particle_pid[i]==4900111){
        float theta=Particle_theta[i];
		float phi=Particle_phi[i];
		float decay_y=Particle_Decaylength[l]*sin(theta)*sin(phi)*1000;
		Particle_y_meson[j]=decay_y;
		Particle_y_meson[j+1]=decay_y;
		l=l+1;
		j=j+2;
	    }
		}
    return Particle_y_meson;
}

RVecF Particle_Y(RVecI Particle_pid, RVecF Particle_y_meson, RVecF Particle_status, RVecF Particle_x, RVecF Particle_y) {
    int count_muons=0;

    int norm=1;
	int l=0;

    RVecF N_Particle_y(Particle_pid.size());

    for (size_t i = 0; i < Particle_pid.size(); i++){
		N_Particle_y[i]=Particle_y[i];
		if (fabs(Particle_pid[i])==13 && Particle_status[i]==1 && sqrt(pow(Particle_x[i],2)+pow(Particle_y[i],2))>10){
			count_muons=count_muons+1;
		float y = Particle_y_meson[l];
		float y1= y * norm;
		N_Particle_y[i]=y1;
		l=l+1;
	    }
	}
    return N_Particle_y;
}

RVecF Particle_Z_Meson(RVecI Particle_pid, RVecF Particle_DecayTime, RVecF Particle_theta, RVecF Particle_phi, RVecF Particle_Decaylength) {
    int count_mesons=0;

    int l=0;
	int j=0;
    for (size_t k = 0; k < Particle_pid.size(); k++){
	if (Particle_pid[k]==4900111){
	    count_mesons=count_mesons+1;
		}
    }
    RVecF Particle_z_meson(count_mesons*2);

    for (size_t i = 0; i < Particle_pid.size(); i++){
	if (Particle_pid[i]==4900111){
        float theta=Particle_theta[i];
		float phi=Particle_phi[i];
		float decay_z=Particle_Decaylength[l]*cos(theta)*1000;
		Particle_z_meson[j]=decay_z;
		Particle_z_meson[j+1]=decay_z;
		l=l+1;
		j=j+2;
	    }
		}
    return Particle_z_meson;
}

RVecF Particle_Z(RVecI Particle_pid, RVecF Particle_z_meson, RVecF Particle_status, RVecF Particle_x, RVecF Particle_y, RVecF Particle_z) {
    int count_muons=0;

    int norm=1;
	int l=0;
    RVecF N_Particle_z(Particle_pid.size());

    for (size_t i = 0; i < Particle_pid.size(); i++){
		N_Particle_z[i]=Particle_z[i];
		if (fabs(Particle_pid[i])==13 && Particle_status[i]==1 && sqrt(pow(Particle_x[i],2)+pow(Particle_y[i],2))>10){
		float z = Particle_z_meson[l];
		float z1= z * norm;
		N_Particle_z[i]=z1;
		l=l+1;
	    }
	}
    return N_Particle_z;
}

RVecF Particle_X_CMS(RVecI Particle_pid, RVecF N_Particle_x){
	float pi=3.1415926535897;
	float d=27000000/pi;

RVecF N_Particle_x_CMS(Particle_pid.size());

for(int i=0;i<Particle_pid.size();i++){
	N_Particle_x_CMS[i]=d-N_Particle_x[i];
}
return N_Particle_x_CMS;
}

RVecF Particle_Y_CMS(RVecI Particle_pid, RVecF N_Particle_y){

RVecF N_Particle_y_CMS(Particle_pid.size());

for(int i=0;i<Particle_pid.size();i++){
	N_Particle_y_CMS[i]=N_Particle_y[i];
}
return N_Particle_y_CMS;
}

RVecF Particle_Z_CMS(RVecI Particle_pid, RVecF N_Particle_z){

RVecF N_Particle_z_CMS(Particle_pid.size());

for(int i=0;i<Particle_pid.size();i++){
	N_Particle_z_CMS[i]=N_Particle_z[i];
}
return N_Particle_z_CMS;
}

RVecI Particle_Check(RVecI Particle_pid, RVecF N_Particle_x_CMS, RVecF N_Particle_y_CMS, RVecF N_Particle_z_CMS){
	float R=7500;//mm
	float L=21000;//mm
	RVecI Particle_check(Particle_pid.size());
	for(int i=0;i<Particle_pid.size();i++){
		float x_CMS=N_Particle_x_CMS[i];
		float y_CMS=N_Particle_y_CMS[i];
		float z_CMS=N_Particle_z_CMS[i];
		float xy_CMS=sqrt(pow(x_CMS,2)+pow(y_CMS,2));
		if (xy_CMS < R && fabs(z_CMS)<=L/2){
			int check=0;//inside
			Particle_check[i]=check;
		}
		if (xy_CMS>R){
			int check=1;//outside
			Particle_check[i]=check;
		}
		if (xy_CMS<R && z_CMS<-L/2){
			int check=2;//outside-inside-left
			Particle_check[i]=check;
		}
		if(xy_CMS<R && z_CMS>L/2){
			int check=3;//outside-inside-right
			Particle_check[i]=check;
		}	
	}
	return Particle_check;
}

RVecF Particle_TOF1(RVecI Particle_pid, RVecF Particle_px, RVecF Particle_py, RVecF Particle_pz, RVecF Particle_energy, RVecF N_Particle_x,  RVecF N_Particle_y,  RVecF N_Particle_z){
float pi=3.1415926535897;
float d=27000000/pi;//mm
float R=7500;//mm
float L=21000;//mm
int a=0;
int b=0;
int c=-1;
int x0=d;
int y0=0;
int z0=0;

float Particle_tof1=-999.;

RVecF Particle_Tof1(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
	float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
	float beta=p/Particle_energy[i];
	float gamma=1/sqrt(1-beta*beta);
	float beta_x=Particle_px[i]/Particle_energy[i];
	float beta_y=Particle_py[i]/Particle_energy[i];
	float beta_z=Particle_pz[i]/Particle_energy[i];
	float vx=beta_x*gamma;
	float vy=beta_y*gamma;
	float vz=beta_z*gamma;

	float x1=N_Particle_x[i];
	float y1=N_Particle_y[i];
	float z1=N_Particle_z[i];

	float a1i=c*c*vy*vy+b*b*vz*vz-2*b*c*vy*vz;
    float a1j=a*a*vz*vz+c*c*vx*vx-2*c*a*vz*vx;
    float a1k=b*b*vx*vx+a*a*vy*vy-2*a*b*vx*vy;
              
    float b1i=2*c*c*y1*vy-2*c*c*y0*vy+2*b*b*z1*vz-2*b*b*z0*vz+2*b*c*y0*vz-2*b*c*y1*vz+2*b*c*vy*z0-2*b*c*z1*vy;
    float b1j=2*a*a*z1*vz-2*a*a*z0*vz+2*c*c*x1*vx-2*c*c*x0*vx+2*c*a*z0*vx-2*c*a*z1*vx+2*c*a*vz*x0-2*c*a*x1*vz;
    float b1k=2*b*b*x1*vx-2*b*b*x0*vx+2*a*a*y1*vy-2*a*a*y0*vy+2*a*b*x0*vy-2*a*b*x1*vy+2*a*b*vx*y0-2*a*b*y1*vx;

              
    float c1i=c*c*y0*y0+c*c*y1*y1-2*c*c*y0*y1+b*b*z0*z0+b*b*z1*z1-2*b*b*z0*z1-2*b*c*y0*z0+2*b*c*y0*z1+2*b*c*y1*z0-2*b*c*y1*z1;     
    float c1j=a*a*z0*z0+a*a*z1*z1-2*a*a*z0*z1+c*c*x0*x0+c*c*x1*x1-2*c*c*x0*x1-2*c*a*z0*x0+2*c*a*z0*x1+2*c*a*z1*x0-2*c*a*z1*x1;
    float c1k=b*b*x0*x0+b*b*x1*x1-2*b*b*x0*x1+a*a*y0*y0+a*a*y1*y1-2*a*a*y0*y1-2*a*b*x0*y0+2*a*b*x0*y1+2*a*b*x1*y0-2*a*b*x1*y1;    
              
    float a1f=a1i+a1j+a1k;
    float b1f=b1i+b1j+b1k;
    float c1f=c1i+c1j+c1k-R*R;
    float del1=b1f*b1f-4*a1f*c1f;

    if (del1>0){
		Particle_tof1=(-b1f-sqrt(b1f*b1f-4*a1f*c1f))/(2*a1f);
		Particle_Tof1[i]=Particle_tof1;
	}
	Particle_Tof1[i]=Particle_tof1;
	}
return Particle_Tof1;
}


RVecF Particle_TOF2(RVecI Particle_pid, RVecF Particle_px, RVecF Particle_py, RVecF Particle_pz, RVecF Particle_energy, RVecF N_Particle_x,  RVecF N_Particle_y,  RVecF N_Particle_z){

float pi=3.1415926535897;
float d=27000000/pi;
float R=7500;
float L=21000;
int a=0;
int b=0;
int c=-1;
float ox_left=d;
float oy_left=0;
float oz_left=-L/2;

float ox_left_CMS=0;
float oy_left_CMS=0;
float oz_left_CMS=L/2;

RVecF Particle_Tof2(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
	float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
	float beta=p/Particle_energy[i];
	float gamma=1/sqrt(1-beta*beta);
	float beta_x=Particle_px[i]/Particle_energy[i];
	float beta_y=Particle_py[i]/Particle_energy[i];
	float beta_z=Particle_pz[i]/Particle_energy[i];
	float vx=beta_x*gamma;
	float vy=beta_y*gamma;
	float vz=beta_z*gamma;


	float x1=N_Particle_x[i];
	float y1=N_Particle_y[i];
	float z1=N_Particle_z[i];

	float Particle_tof2=(a*(ox_left-x1)+b*(oy_left-y1)+c*(oz_left-z1))/(a*vx+b*vy+c*vz);
	Particle_Tof2[i]=Particle_tof2;
}
return Particle_Tof2;
}

RVecF Particle_TOF3(RVecI Particle_pid, RVecF Particle_px, RVecF Particle_py, RVecF Particle_pz, RVecF Particle_energy, RVecF N_Particle_x,  RVecF N_Particle_y,  RVecF N_Particle_z){

float pi=3.1415926535897;
float d=27000000/pi;
float R=7500;
float L=21000;
int a=0;
int b=0;
int c=1;
float ox_right=d;
float oy_right=0;
float oz_right=L/2;

float ox_right_CMS=0;
float oy_right_CMS=0;
float oz_right_CMS=-L/2;

RVecF Particle_Tof3(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
	float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
	float beta=p/Particle_energy[i];
	float gamma=1/sqrt(1-beta*beta);
	float beta_x=Particle_px[i]/Particle_energy[i];
	float beta_y=Particle_py[i]/Particle_energy[i];
	float beta_z=Particle_pz[i]/Particle_energy[i];
	float vx=beta_x*gamma;
	float vy=beta_y*gamma;
	float vz=beta_z*gamma;


	float x1=N_Particle_x[i];
	float y1=N_Particle_y[i];
	float z1=N_Particle_z[i];

	float Particle_tof3=(a*(ox_right-x1)+b*(oy_right-y1)+c*(oz_right-z1))/(a*vx+b*vy+c*vz);
	Particle_Tof3[i]=Particle_tof3;
}
return Particle_Tof3;
}

RVecF N_Particle_X_Final1(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_x, RVecF Particle_Tof1, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_x_final1(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_x=Particle_px[i]/Particle_energy[i];
float vx=beta_x*gamma;
float x1i=N_Particle_x[i];
float time1=Particle_Tof1[i];
float x1f=x1i+vx*time1;
N_Particle_x_final1[i]=x1f;
}
return N_Particle_x_final1;
}

RVecF N_Particle_Y_Final1(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_y, RVecF Particle_Tof1, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_y_final1(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_y=Particle_py[i]/Particle_energy[i];
float vy=beta_y*gamma;
float y1i=N_Particle_y[i];
float time1=Particle_Tof1[i];
float y1f=y1i+vy*time1;
N_Particle_y_final1[i]=y1f;
}
return N_Particle_y_final1;
}

RVecF N_Particle_Z_Final1(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_z, RVecF Particle_Tof1, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_z_final1(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_z=Particle_pz[i]/Particle_energy[i];
float vz=beta_z*gamma;
float z1i=N_Particle_z[i];
float time1=Particle_Tof1[i];
float z1f=z1i+vz*time1;
N_Particle_z_final1[i]=z1f;
}
return N_Particle_z_final1;
}

RVecF N_Particle_X_Final2(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_x, RVecF Particle_Tof2, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_x_final2(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_x=Particle_px[i]/Particle_energy[i];
float vx=beta_x*gamma;
float x2i=N_Particle_x[i];
float time2=Particle_Tof2[i];
float x2f=x2i+vx*time2;
N_Particle_x_final2[i]=x2f;
}
return N_Particle_x_final2;
}

RVecF N_Particle_Y_Final2(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_y, RVecF Particle_Tof2, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_y_final2(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_y=Particle_py[i]/Particle_energy[i];
float vy=beta_y*gamma;
float y2i=N_Particle_y[i];
float time2=Particle_Tof2[i];
float y2f=y2i+vy*time2;
N_Particle_y_final2[i]=y2f;
}
return N_Particle_y_final2;
}

RVecF N_Particle_Z_Final2(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_z, RVecF Particle_Tof2, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_z_final2(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_z=Particle_pz[i]/Particle_energy[i];
float vz=beta_z*gamma;
float z2i=N_Particle_z[i];
float time2=Particle_Tof2[i];
float z2f=z2i+vz*time2;
N_Particle_z_final2[i]=z2f;
}
return N_Particle_z_final2;
}

RVecF N_Particle_X_Final3(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_x, RVecF Particle_Tof3, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_x_final3(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_x=Particle_px[i]/Particle_energy[i];
float vx=beta_x*gamma;
float x3i=N_Particle_x[i];
float time3=Particle_Tof3[i];
float x3f=x3i+vx*time3;
N_Particle_x_final3[i]=x3f;
}
return N_Particle_x_final3;
}

RVecF N_Particle_Y_Final3(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_y, RVecF Particle_Tof3, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_y_final3(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_y=Particle_py[i]/Particle_energy[i];
float vy=beta_y*gamma;
float y3i=N_Particle_y[i];
float time3=Particle_Tof3[i];
float y3f=y3i+vy*time3;
N_Particle_y_final3[i]=y3f;
}
return N_Particle_y_final3;
}

RVecF N_Particle_Z_Final3(RVecI Particle_pid, RVecF Particle_px, RVecF N_Particle_z, RVecF Particle_Tof3, RVecF Particle_energy, RVecF Particle_py, RVecF Particle_pz){

RVecF N_Particle_z_final3(Particle_pid.size());

for(size_t i=0; i<Particle_pid.size();i++){
float p=sqrt(pow(Particle_px[i],2)+pow(Particle_py[i],2)+pow(Particle_pz[i],2));
float beta=p/Particle_energy[i];
float gamma=1/sqrt(1-beta*beta);
float beta_z=Particle_pz[i]/Particle_energy[i];
float vz=beta_z*gamma;
float z3i=N_Particle_z[i];
float time3=Particle_Tof3[i];
float z3f=z3i+vz*time3;
N_Particle_z_final3[i]=z3f;
}
return N_Particle_z_final3;
}

RVecF Particle_Check_Final(RVecI Particle_check, RVecI Particle_pid, RVecF N_Particle_x_final1, RVecF N_Particle_x_final2, RVecF N_Particle_x_final3, RVecF N_Particle_y_final1, RVecF N_Particle_y_final2, RVecF N_Particle_y_final3, RVecF N_Particle_z_final1, RVecF N_Particle_z_final2, RVecF N_Particle_z_final3, RVecF Particle_Tof1, RVecF Particle_Tof2, RVecF Particle_Tof3, RVecF Particle_status){

float pi=3.1415926535897;
float d=27000000/pi;
float L=21000;
float R=7500;

RVecI Particle_check_final(Particle_pid.size());

for(int i=0;i<Particle_pid.size();i++){

	int check=Particle_check[i];
	int check_final=-999;

	float x1f=N_Particle_x_final1[i];
	float y1f=N_Particle_y_final1[i];
	float z1f=N_Particle_z_final1[i];

	float x2f=N_Particle_x_final2[i];
	float y2f=N_Particle_y_final2[i];
	float z2f=N_Particle_z_final2[i];

	float x3f=N_Particle_x_final3[i];
	float y3f=N_Particle_y_final3[i];
	float z3f=N_Particle_z_final3[i];
	
	float x1f_CMS=d-x1f;
	float y1f_CMS=y1f;
	float z1f_CMS=-z1f;

	float x2f_CMS=d-x2f;
	float y2f_CMS=y2f;
	float z2f_CMS=-z2f;

	float x3f_CMS=d-x3f;
	float y3f_CMS=y3f;
	float z3f_CMS=-z3f;

	float time1=Particle_Tof1[i];
	float time2=Particle_Tof2[i];
	float time3=Particle_Tof3[i];

	float ox_right_CMS=0;
	float oy_right_CMS=0;
	float oz_right_CMS=-L/2;

	float ox_left_CMS=0;
	float oy_left_CMS=0;
	float oz_left_CMS=L/2;

    float dis_m2=sqrt(pow(x2f_CMS-ox_left_CMS,2)+pow(y2f_CMS-oy_left_CMS,2)+pow(z2f_CMS-oz_left_CMS,2));
	float dis_m3=sqrt(pow(x3f_CMS-ox_right_CMS,2)+pow(y3f_CMS-oy_right_CMS,2)+pow(z3f_CMS-oz_right_CMS,2));

	if (check==0){
		check_final=0;
		Particle_check_final[i]=check_final;
	}
	if (check==1 && time1 > 0 && fabs(x1f_CMS)<L/2){
		check_final=1;//hit barrel
		Particle_check_final[i]=check_final;
	}
	if(check==1 && time1>0 && fabs(x1f_CMS)>L/2 && dis_m2 <= R && dis_m3 > R){
		check_final=2;//hit left endcap
		Particle_check_final[i]=check_final;
	}
	if(check==1 && time1>0 && fabs(x1f_CMS)>L/2 && dis_m2 > R && dis_m3 <= R){
		check_final=3;//hit right endcap
		Particle_check_final[i]=check_final;
	}
	if (check==2 && dis_m2<=R){
		check_final=2;
		Particle_check_final[i]=check_final;
	}
	if (check==3 && dis_m3<=R){
		check_final=3;
		Particle_check_final[i]=check_final;
	}
	Particle_check_final[i]=check_final;

}
return Particle_check_final;
}

RVecF Particle_Depth(RVecI Particle_pid, RVecI Particle_check_final, RVecF N_Particle_x_final1, RVecF N_Particle_x_final2, RVecF N_Particle_x_final3, RVecF N_Particle_y_final1, RVecF N_Particle_y_final2, RVecF N_Particle_y_final3, RVecF N_Particle_z_final1, RVecF N_Particle_z_final2, RVecF N_Particle_z_final3, RVecF N_Particle_x, RVecF N_Particle_y, RVecF N_Particle_z){
float depth=0.;
RVecI Particle_depth(Particle_pid.size());

for(int i=0;i<Particle_pid.size();i++){ 
	float x1f=N_Particle_x_final1[i];
	float y1f=N_Particle_y_final1[i];
	float z1f=N_Particle_z_final1[i];

	float x2f=N_Particle_x_final2[i];
	float y2f=N_Particle_y_final2[i];
	float z2f=N_Particle_z_final2[i];

	float x3f=N_Particle_x_final3[i];
	float y3f=N_Particle_y_final3[i];
	float z3f=N_Particle_z_final3[i];

	float xi=N_Particle_x[i];
	float yi=N_Particle_y[i];
	float zi=N_Particle_z[i];



	int check_final=Particle_check_final[i];

	if (check_final==1){
		 depth=sqrt(pow(x1f-xi,2)+pow(y1f-yi,2)+pow(z1f-zi,2));
		
	}

	if (check_final==2){
		 depth=sqrt(pow(x2f-xi,2)+pow(y2f-yi,2)+pow(z2f-zi,2));
		
	}

	if (check_final==3){
		 depth=sqrt(pow(x3f-xi,2)+pow(y3f-yi,2)+pow(z3f-zi,2));
		
	}
	Particle_depth[i]=depth;


}
return Particle_depth;
}





































