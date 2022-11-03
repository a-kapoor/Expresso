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

RVecI Get_P_Position(RVecI Particle_PID, int PID){
    int count_p=0;
    int p_position=0;
    int l=0;
    for(size_t k=0; k<Particle_PID.size(); k++){
        if (Particle_PID[k]==PID){
            count_p=count_p+1;
        }
    }
    RVecI P_Position(count_p);
    for(size_t i=0;i<Particle_PID.size();i++){
        if (Particle_PID[i]==PID){
            P_Position[l]=i;
            l=l+1;
        }
    }
    return P_Position;
}

RVecF Get_P_Vector(RVecI P_Position, RVecF Particle_Vector){
    int l=0;
    RVecF New_Particle_Vector(P_Position.size());
    for(size_t i=0;i<P_Position.size();i++){
            New_Particle_Vector[l]=Particle_Vector[P_Position[i]];
            l=l+1;
    }
    return New_Particle_Vector;
}








