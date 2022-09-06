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

RVecF Get_Photon_PT(RVecI Particle_PID, RVecF Particle_PT){
    int count_photon=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==22){
            count_photon=count_photon+1;
        }
    }
    RVecF Photon_PT(count_photon);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==22){
            Photon_PT[l]=Particle_PT[i];
            l=l+1;
        }
    }
    return Photon_PT;
}

RVecF Get_Photon_Px(RVecI Particle_PID, RVecF Particle_Px){
    int count_photon=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==22){
            count_photon=count_photon+1;
        }
    }
    RVecF Photon_Px(count_photon);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==22){
            Photon_Px[l]=Particle_Px[i];
            l=l+1;
        }
    }
    return Photon_Px;
}

RVecF Get_Photon_Py(RVecI Particle_PID, RVecF Particle_Py){
    int count_photon=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==22){
            count_photon=count_photon+1;
        }
    }
    RVecF Photon_Py(count_photon);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==22){
            Photon_Py[l]=Particle_Py[i];
            l=l+1;
        }
    }
    return Photon_Py;
}


RVecF Get_Photon_Pz(RVecI Particle_PID, RVecF Particle_Pz){
    int count_photon=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==22){
            count_photon=count_photon+1;
        }
    }
    RVecF Photon_Pz(count_photon);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==22){
            Photon_Pz[l]=Particle_Pz[i];
            l=l+1;
        }
    }
    return Photon_Pz;
}

RVecF Get_Photon_energy(RVecI Particle_PID, RVecF Particle_energy){
    int count_photon=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==22){
            count_photon=count_photon+1;
        }
    }
    RVecF Photon_energy(count_photon);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==22){
            Photon_energy[l]=Particle_energy[i];
            l=l+1;
        }
    }
    return Photon_energy;
}

RVecF Get_Photon_mass(RVecI Particle_PID, RVecF Particle_mass){
    int count_photon=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==22){
            count_photon=count_photon+1;
        }
    }
    RVecF Photon_mass(count_photon);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==22){
            Photon_mass[l]=Particle_mass[i];
            l=l+1;
        }
    }
    return Photon_mass;
}

RVecF Get_Photon_eta(RVecI Particle_PID, RVecF Particle_eta){
    int count_photon=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==22){
            count_photon=count_photon+1;
        }
    }
    RVecF Photon_eta(count_photon);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==22){
            Photon_eta[l]=Particle_eta[i];
            l=l+1;
        }
    }
    return Photon_eta;
}

RVecF Get_Photon_phi(RVecI Particle_PID, RVecF Particle_phi){
    int count_photon=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==22){
            count_photon=count_photon+1;
        }
    }
    RVecF Photon_phi(count_photon);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==22){
            Photon_phi[l]=Particle_phi[i];
            l=l+1;
        }
    }
    return Photon_phi;
}









