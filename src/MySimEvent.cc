#include "MySimEvent.hh"


MySimEvent::MySimEvent(){
  ClearEvent();
}


void MySimEvent::ClearEvent() {
  ETotCal  = 0;
  for(int i=0;i<1000;i++) {
    ecry[i] = 0;
  }

  
  particles.clear();
  clusters.clear();
  
}
void MySimEvent::SetCrystalEnergies(double *en, int n){
  if(n>1000) {
    std::cout << "Event with more than 1000 cells: storing only 1000 energies" << std::endl;
    n=1000;
  }
  for(int i=0;i<n;i++) {
    ecry[i] = en[i];
  }
}
void MySimEvent::AddCluster(double x, double y,double e){
  MyCaloCluster clus;
  clus.pos[0] = x;
  clus.pos[1] = y;
  clus.energy = e;
  clusters.push_back(clus);
}
