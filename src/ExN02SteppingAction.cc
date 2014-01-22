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
// $Id: ExN02SteppingAction.cc,v 1.1 2014/01/22 15:35:03 veni Exp $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4TrackStatus.hh"
#include "HistoManager.hh"
#include "MyEvent.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN02SteppingAction::ExN02SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  //  //Devi capire che fare con il monitor e che cosa succede con piu' di un interazione!!!!   |-------|
  G4Track* track = step->GetTrack();
  

  MySimEvent *evt = (MyEvent::GetInstance())->GetSimEvent();


  //Storing everything that comes out of the target
  if(step->GetPostStepPoint()->GetPhysicalVolume()!=0){
    if(step->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="Target") {
      G4String lastProc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
      const G4TrackVector * trVec = step->GetSecondary ();
      //std::cout << "Number of secondaries produced:  " << trVec->size() << "  Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;

      //      std::vector<G4Track*>::iterator  it;
     
      for(int i=0;i<trVec->size();i++) {
	G4Track *trk = trVec->at(i);
	//G4cout << "Particle type:  " << trk->GetDefinition()->GetParticleName () << G4endl;
	//Get the info for the secondary:

	MyParticle part(trk->GetDefinition()->GetPDGMass ());
	part.setProcess(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().data());
	part.setType(trk->GetDefinition()->GetParticleName ().data());
	part.setEnergy(trk->GetTotalEnergy());
	part.calcMomentum(trk->GetMomentumDirection ().x(), 
			  trk->GetMomentumDirection ().y(),
			  trk->GetMomentumDirection ().z());
	part.setPVtx(trk->GetPosition ().x(),
		     trk->GetPosition ().y(),
		     trk->GetPosition ().z());

	evt->AddParticle(part);
      }
      
      
    }
  }

  //Cerca il primario  
  if(track->GetTrackID()==1){ //primary particle
    if(track->GetParticleDefinition() == G4Positron::PositronDefinition()){
      if(step->GetPostStepPoint()->GetPhysicalVolume()!=0){
	if(step->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="Target") {
	  //      G4cout<<"track->GetParticleDefinition() "<<track->GetParticleDefinition()<<G4endl;
	  G4String lastProc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	  //	  G4cout<<lastProc<<G4endl;
	  if(lastProc=="eBrem"){
	    ProcID=1.;
	  }else if(lastProc=="annihil"){ 
	    ProcID=2.;
	  }else if(lastProc=="eIoni"){
	    ProcID=3.;
	  }else{
	    ProcID=0.;
	  }
	  //	  double BeamPartE  = track->GetTotalEnergy();
	  //	  G4cout<<"Process: "<<lastProc<< " code " << ProcID << " " <<BeamPartE<<G4endl;
	  //	  G4cout<<"Volume: "<<step->GetPostStepPoint()->GetPhysicalVolume()<<G4endl;
	  //	  G4cout<<"Volume: "<<step->GetPostStepPoint()->GetPhysicalVolume()->GetName()<<G4endl;
	}
      }
    }
  }

  //Cerca il primo gamma in uscita   
  if(track->GetParticleDefinition() == G4Gamma::GammaDefinition()){
    if(step->GetPostStepPoint()->GetPhysicalVolume()!=0){
      if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Target") {
	if(track->GetCurrentStepNumber()==1 && track->GetParentID()==1){
	  //      G4cout<<"track->GetParticleDefinition() "<<track->GetParticleDefinition()<<G4endl;
	  //	  G4String lastProc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	  //	  G4int parent      = track->GetParentID();
	  //	  G4int TrID        = track->GetTrackID();
	  //	  G4bool Primo      = step->IsFirstStepInVolume();
	  GammaE    = track->GetTotalEnergy();
	  VertexPos = track->GetVertexPosition();
	  GammaDir  = track->GetVertexMomentumDirection();
	  G4ThreeVector BeamDir;
	  BeamDir[0]=0.;
	  BeamDir[1]=0.;
	  BeamDir[2]=1.;
	  ThetaGamma= SetGammaAngle(GammaDir,BeamDir);

	  //  for(int i=0; i<4; i++) P4Miss[i]=TargetEleMom[i]+BeamMom[i]-GMom[i];
	  //  G4double Mmiss2 = P4Miss[3]*P4Miss[3]-P4Miss[2]*P4Miss[2]-P4Miss[1]*P4Miss[1]-P4Miss[0]*P4Miss[0];
	  //  NChild++;
	  //  G4cout<<"Theta Gamma " <<NChild<<G4endl;
	}
      }
    }
  }  
}

G4double ExN02SteppingAction::SetGammaAngle(G4ThreeVector GammaDir,G4ThreeVector BeamDir)
{
  G4double product=0;
  G4double Theta;
  for (int i=0; i<3; i++)  product+= GammaDir[i]*BeamDir[i];
  Theta = acos (product) * 180.0 / 3.14159265;  //in degreees
  return Theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

