//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: ExN02PrimaryGeneratorAction.cc,v 1.2 2014/01/22 16:56:13 veni Exp $
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "MyEventGeneration.hh"
#include "Constants.hh"
#include "MyEvent.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02PrimaryGeneratorAction::ExN02PrimaryGeneratorAction(
                                               ExN02DetectorConstruction* myDC)
:myDetector(myDC)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  evt = MyEvent::GetInstance()->GetGenEvent();
  evt->SetUbosonMass(UMass*MeV);
  evt->SetUbosonEps(epsilon);
  evt->SetBeamEnergy(BeamEnergy*MeV);


  // default particle  

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e+");
  
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(BeamEnergy*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02PrimaryGeneratorAction::~ExN02PrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  std::cout << "============================ Generate Primaries" << std::endl;
  evt->ClearEvent();
  
  static int nev;
  nev++;
  
  if(nev%10000 == 0) 
    std::cout << "Generating event numner " << nev << std::endl;
  
  G4double position = -0.5*(myDetector->GetWorldFullLength());
  position = -1.*m;
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));  
  
  static G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();   
  static G4ParticleDefinition* particle ;
  
  double xf = G4UniformRand();
  //to get a fraction 1/N is the number of U you get
  //  if(xf > 1./10000.*epsilon) {

  if(IsCalibRun==1){
    std::cout << "==== Generating calibration event === " << std::endl;


    G4double ECalX      = ECalSizeX*cm;
    G4double ECalY      = ECalSizeY*cm;
    G4double ECalLength = ECalSizeZ*cm;
    G4double ECryLength = ECalLength;
    G4double ECryX      = ECalX/ECalNRow;
    G4double ECryY      = ECalY/ECalNCol;
    G4int i=((nev-1)/ECalNCol)%ECalNRow; //put nrow Nrow !!!
    G4int j=(nev-1)%ECalNCol;
    G4double PosXCry=-ECalX*0.5+0.5*ECryX+i*ECryX;
    G4double PosYCry=-ECalY*0.5+0.5*ECryY+j*ECryY;
    G4ThreeVector positionCry = G4ThreeVector(PosXCry,PosYCry,ECalPosiZ*cm-ECryLength/2.);
    G4int HoleFlag=0; //should be parametric in NRow NCol
    if( fabs(PosXCry)<ECalInnHole*cm && fabs(PosYCry)<ECalInnHole*cm ){
      HoleFlag=1;  //means is in the hole 
    }else{
      HoleFlag=0;
    }	
    if(HoleFlag==0){
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
      particleGun->SetParticleDefinition(particle);
      particleGun->SetNumberOfParticles(NPrimaries);  
      particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm)); 
      particleGun->SetParticleMomentumDirection(positionCry.unit()); //with Y=0 goes to the Gap!!!

      //Fill all the information about the primaries before generating the event:
      MyParticle part(particleGun->GetParticleDefinition ()->GetPDGMass ());
      part.setType(particleGun->GetParticleDefinition ()->GetParticleName ()); 
      part.setProcess("beam");

      part.setEnergy(particleGun->GetParticleEnergy());
      part.calcMomentum(
			particleGun->GetParticleMomentumDirection().x(),
			particleGun->GetParticleMomentumDirection().y(),
			particleGun->GetParticleMomentumDirection().z()  );
      part.setPVtx(particleGun->GetParticlePosition ().x(),
		   particleGun->GetParticlePosition ().y(),
		   particleGun->GetParticlePosition ().z());
      evt->AddProduct(part);

      particleGun->GeneratePrimaryVertex(anEvent);
    }
  } else {
    std::cout << "==== Generating general event === " << std::endl;
    if(xf>1E-14) {
      G4double position = -0.5*(myDetector->GetWorldFullLength());
      particleGun->SetNumberOfParticles(NPrimaries);  
      //for  energy resolution use
      //    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
      //   particleGun->SetParticleDefinition(particle);
      //  G4ThreeVector positionCry = G4ThreeVector(5.5*cm,5.5*cm,200.*cm);
      // particleGun->SetParticleMomentumDirection(positionCry.unit()); //with Y=0 goes to the Gap!!!
      position = -1.*m;
      
      if(BeamSpot==0){
	particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position)); 
      }else{
	G4double StartX = G4RandGauss::shoot(0.,SigmaBeamX);
	G4double StartY = G4RandGauss::shoot(0.,SigmaBeamY);
	//      G4cout<<"Sigma "<<SigmaBeamX<<" Startx "<<StartX << G4endl;
	// particleGun->SetParticlePosition(G4ThreeVector(StartX*cm,StartY*cm,position)); 
	particleGun->SetParticlePosition(G4ThreeVector(StartX,StartY,position)); 
	//      G4cout<< "XY of beam " <<StartX<<" StartY "<<StartY<<G4endl;
      }
      
      if(BeamESpread==1){
	G4double NewBeamEnergy = G4RandGauss::shoot(BeamEnergy,SigmaBeamE);
	particleGun->SetParticleEnergy(NewBeamEnergy*MeV);
      }
      
      //    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.026,0.03,0.95).unit()); //with Y=0 goes to the Gap!!!
      
      MyParticle part(particleGun->GetParticleDefinition ()->GetPDGMass ());
      part.setType(particleGun->GetParticleDefinition ()->GetParticleName ()); 
      part.setProcess("beam");
      part.setEnergy(particleGun->GetParticleEnergy());
      part.calcMomentum(particleGun->GetParticleMomentumDirection().x(),particleGun->GetParticleMomentumDirection().y(),particleGun->GetParticleMomentumDirection().z()  );
      part.setPVtx(particleGun->GetParticlePosition ().x(),particleGun->GetParticlePosition ().y(),particleGun->GetParticlePosition ().z());
      evt->AddProduct(part);
      
      particleGun->GeneratePrimaryVertex(anEvent);
    } else {  //U Boson MC 
      G4double NewBeamEnergy = BeamEnergy;
      if(BeamESpread==1){
	NewBeamEnergy = G4RandGauss::shoot(BeamEnergy,SigmaBeamE);
	//particleGun->SetParticleEnergy(NewBeamEnergy*MeV);
      }
      evt->SetBeamEnergy(NewBeamEnergy);
      evt->GenerateEvent();
      
      G4double StartX=0;
      G4double StartY=0;
      
      if(BeamSpot==1){
	StartX= G4RandGauss::shoot(0.,SigmaBeamX);
	StartY= G4RandGauss::shoot(0.,SigmaBeamY);
      }
      
      //Loop over the final state particles and trace them through the detector
      std::vector<MyParticle>* particles = evt->getParticles();
      std::vector<MyParticle>::iterator it;
      
      //    std::cout << "Number of primaries :" << particles->size() << std::endl;
      
      it = particles->begin();
      while (it != particles->end()) {
	if(1) {
	  //if(it->getType()=="gamma") {
	  //	std::cout << "Found particle " << it->getType() << std::endl;
	  particle = particleTable->FindParticle(it->getType());
	  particleGun->SetParticleDefinition(particle);
	  
	  double *p = it->get4Momentum();
	  const double *pvtx = it->getPVtx();
	  
	  particleGun->SetParticleEnergy(p[0]);
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(p[1],p[2],p[3]).unit());
	  //      G4double position = -0.5*(myDetector->GetWorldFullLength());
	  
	  G4double  position = -0.*m; //center of the target for the moment
	  //	particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
	  
	  std::cout << "Particle " << it->getType() << "  Production vertex: " 
		    << pvtx[0] << "\t" << pvtx[1] << "\t" << pvtx[2] << "\t" 
		    << std::endl;
	  particleGun->SetParticlePosition(G4ThreeVector(pvtx[0]*m+StartX,pvtx[1]*m+StartY,pvtx[2]*m ));
	  //add the randoms
	  particleGun->GeneratePrimaryVertex(anEvent);
	  //	G4cout<<"Ciao "<<pvtx[0]<<G4endl;
	}
	it++;
      }
    }
  }

  
  MyEvent *TheEvent = MyEvent::GetInstance();
  
  G4cout << "ExN02PrimaryGeneratorAction::GeneratePrimaries  Number of primaries in that event: "<< TheEvent->GetGenEvent()->getParticles()->size() << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

