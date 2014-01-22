
#include "MyEventGeneration.hh"

MyEventGenerator::MyEventGenerator()
 :virtualState(0){
}


MyEventGenerator::~MyEventGenerator(){
}


void MyEventGenerator::GenerateEvent(){
  ClearEvent();
  CreateInitialState();
  CreateIntermediateState();
  CreateFinalState();
}


void MyEventGenerator::ClearEvent(){
  primaries.clear();
  virtualState.clear();
  products.clear();
}

void MyEventGenerator::CreateInitialState(){
  MyParticle positron(CLHEP::electron_mass_c2); //This sets the particles at rest
  MyParticle electron(CLHEP::electron_mass_c2);
    
  positron.setEnergy(positronE);
  positron.calcMomentum(0.,0.,1.);
  //

  primaries.push_back(positron);
  primaries.push_back(electron); 
}



void MyEventGenerator::CreateFinalState(){
  //Really doing the decay of the particles
  
  //Define the process: X-> gU -> g e+e-
  
  MyParticle Uboson;
  MyParticle gamma;
  
  //Do the calculation in the Center of Mass system 
  Uboson.setMass(MU);
  Uboson.setType("Uboson");
  
  


  gamma.setMass(0.);
  gamma.setType("gamma");

  //Do the decay positronium -> gU
  //In the CM:
  
  double eu= (SQR(Uboson.m()) + SQR(virtualState.m()))/(2.*virtualState.m());
  double eg= (-SQR(Uboson.m()) + SQR(virtualState.m()))/(2.*virtualState.m());
  double pu = sqrt(SQR(eu) - SQR(Uboson.m()));
  
  Uboson.setEnergy(eu);
  gamma.setEnergy(eg);
  
  //  std::cout << "Energy of the U: " << eu << "\t Energy of g: "<< eg << std::endl;
  //  std::cout << "Momentum of the U: " << pu << "\t Momentum of g: "<< eg << std::endl;
  
  //Generate a random vector
  G4ThreeVector dir = G4RandomDirection();
  //  std::cout << "Random direction: " << dir.x() << "\t" << dir.y() << "\t" << dir.z() << std::endl;
  
  Uboson.calcMomentum(dir.x(),dir.y(),dir.z());
  gamma.calcMomentum(-dir.x(),-dir.y(),-dir.z());
  
  gamma.boostLAB(virtualState.get4Momentum());
  //  std::cout << "Gamma energy :"<< gamma.E() << std::endl;
  Uboson.boostLAB(virtualState.get4Momentum());


  //The primary decay is generated, set some parameters for U and get kinematics variables
  Uboson.calcGamma();
  Uboson.setPars(&eps,1);
  Uboson.calcWidthU();
  Uboson.calcTau();
  Uboson.calcBeta();

  //Generate a decay point
  Uboson.genDecayVtx();  

  products.push_back(gamma);
  //Now perform the decay of the U-boson, if possible:
  if(Uboson.m() > 2.*CLHEP::electron_mass_c2) {
    
    //In the CM of the U-boson:
    MyParticle ee(CLHEP::electron_mass_c2);
    ee.setType("e-");
    MyParticle ep(CLHEP::electron_mass_c2);
    ep.setType("e+");

    dir = G4RandomDirection();
    
    ee.setEnergy(Uboson.m()/2);
    ep.setEnergy(Uboson.m()/2);
    
    ee.calcMomentum(dir.x(),dir.y(),dir.z());
    ep.calcMomentum(-dir.x(),-dir.y(),-dir.z());
    
    ee.boostLAB(Uboson.get4Momentum());
    ep.boostLAB(Uboson.get4Momentum());
    
    ee.setPVtx(Uboson.getDVtx());
    ep.setPVtx(Uboson.getDVtx());

    products.push_back(ee);
    products.push_back(ep);
  }


  
}

void MyEventGenerator::CreateIntermediateState(){
  //loop over the primaries
  //virtualState.clear();
  std::vector<MyParticle>::iterator it;
  it = primaries.begin();
  while (it != primaries.end()) {
    virtualState += *it;
    it++;
  }
  virtualState.calcMass();
}

