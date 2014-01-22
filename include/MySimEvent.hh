#ifndef _MY_SIM_EVENT_
#define _MY_SIM_EVENT_

#include "MyParticle.hh"
#include <iostream>
#include <vector>
#include "MyCaloCluster.hh"


class MySimEvent{
  
public:
  MySimEvent();
  ~MySimEvent(){;};

  void ClearEvent();
  
  void SetCrystalEnergies(double *en, int n);
  double * GetCrystalEnergies(){return ecry;};
  
  void AddCluster(double x, double y,double e);
  std::vector<MyCaloCluster>* GetClusters(){return &clusters;};

  void AddParticle(MyParticle part){particles.push_back(part);};
  std::vector<MyParticle>* GetParticles(){return &particles;};

  double ecry[1000];
  std::vector<MyCaloCluster> clusters;
  double ETotCal;

  int pid;
  
  //Particles generated from the primary track
  std::vector<MyParticle> particles;
  

};
#endif
