using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;

using Vec_t = const ROOT::RVec<float>&;
using RVecF = RVec<float>;
using RVecI = RVec<int>;

RVecI TagWWbTriplet(Vec_t btagDeepFlavB, Vec_t pT)
{
    RVecI triptag(3);
    int btagindex=0;
    if(btagDeepFlavB[0]>btagDeepFlavB[1] && btagDeepFlavB[0]>btagDeepFlavB[2]) btagindex=0;
    else{
	if(btagDeepFlavB[1]>btagDeepFlavB[2]) btagindex=1;
	else{btagindex=2;}
    }
    float pt=0;
    for(int i=0;i<3;i++){
	if(i!=btagindex){
	    if(pT[i]>pt){
		pt=pT[i];
		triptag[i]=9100;
	    }
	    else{
		triptag[i]=9010;
	    }
	}
	else{
	    triptag[i]=9001;
	}
    }
    cout<<"<-------------->"<<endl;
    cout<<pT[0]<<" "<<pT[1]<<" "<<pT[2]<<" "<<endl;
    cout<<btagDeepFlavB[0]<<" "<<btagDeepFlavB[1]<<" "<<btagDeepFlavB[2]<<" "<<endl;
    cout<<triptag[0]<<" "<<triptag[1]<<" "<<triptag[2]<<" "<<endl;
    return triptag;
}
	
	    
