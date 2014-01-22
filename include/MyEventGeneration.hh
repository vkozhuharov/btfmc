#ifndef _MY_EVENT_GENERATOR_
#define _MY_EVENT_GENERATOR_


#include "globals.hh"
#include "G4Decay.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include <iostream>
#include <vector>
#include "MyParticle.hh"



//using namespace std;

//class  MyParticle;

class MyEventGenerator {
  
public:
  MyEventGenerator();
  ~MyEventGenerator();


  void SetUbosonMass(double val){MU = val;}
  double GetUbosonMass(){return MU;}

  void SetUbosonEps(double val){eps = val;}
  double GetUbosonEps(){return eps;}


  void SetBeamEnergy(double val){positronE = val;}
  double GetBeamEnergy() {return positronE;}

  void AddPrimary(MyParticle &part) {primaries.push_back(part);};
  void AddProduct(MyParticle &part) {products.push_back(part);};

  void GenerateEvent();
  
  void CreateInitialState();
  void CreateIntermediateState();
  void CreateFinalState();





  void ClearEvent();


  std::vector<MyParticle>* getParticles(){return &products;};


  
  
private:
  int nPart;
  double MU;
  double eps;
  double positronE;
  double vtx[3];

  std::vector<MyParticle> primaries;
  MyParticle virtualState;
  std::vector<MyParticle> products;
};

#endif
