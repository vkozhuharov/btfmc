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
// $Id: ExN02EventAction.cc,v 1.2 2014/01/22 16:34:31 veni Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "Randomize.hh"
 
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "ExN02RunAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02SteppingAction.hh"
#include "ExN02DetectorConstruction.hh"
#include "HistoManager.hh"
#include "Constants.hh"
#include "MyEvent.hh"
#include <numeric>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN02EventAction::ExN02EventAction(ExN02RunAction* run, HistoManager* histo)
:fRunAct(run),fHistoManager(histo)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN02EventAction::~ExN02EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN02EventAction::BeginOfEventAction(const G4Event*)
{
  // Get current run and event numbers
  //  G4int nRun = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  //  G4int nEvent = evt->GetEventID();
  //  G4cout << ">>> Run " << nRun << " Event " << nEvent << " started" << G4endl;
  ETotCal  = 0;
  EtotMRod = 0;
  EtotTRod = 0;
  ECalHitT = 0;
  CalEvtT  = 0;
  ClPosX   = 0;
  ClPosY   = 0; 
  ProcID   = 0; 
  NClusters= 0;
  NTracks  = 0;

  memset(&(fHistoManager->myEvt),0,sizeof(NTEvent));
  myStepping->SetPhysProc(0);

  for(G4int i=0;i<ECalNCells;i++){Used[i]=0;}
  for(G4int i=0;i<TrackerNRings;i++){UsedRing[i]=0;}
  for(G4int i=0;i<ECalNCells;i++){Empty[i]=0;}

  //Clear completely the event:
  //  G4cout << "BeginOfEventAction:   " << G4endl; 
  
  MyEvent *TheEvent = MyEvent::GetInstance();
  TheEvent->GetSimEvent()->ClearEvent();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  // Periodic printing
  //if (event_id < 1 || event_id%NPrint == 0) 
	G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  MyEvent *TheEvent = MyEvent::GetInstance();
  MySimEvent *simEvt = TheEvent->GetSimEvent();
  MyEventGenerator *genEvt = TheEvent->GetGenEvent();
  std::vector<MyParticle>::iterator it;

  G4cout << "================================================" << G4endl  ;
  G4cout << "               Event                            " << G4endl  ;
  G4cout << "Number of primaries in that event: "<< TheEvent->GetGenEvent()->getParticles()->size() << G4endl;
  it = genEvt->getParticles()->begin();
  while (it != genEvt->getParticles()->end()) {
    G4cout << "Particle:  " << it->getType() << G4endl;
    it++;
  }

  G4cout << "Number of secondaries in that event: "<< simEvt->GetParticles()->size() << G4endl;
  it = simEvt->GetParticles()->begin();
  while (it != simEvt->GetParticles()->end()) {
    G4cout << "Particle:  " << it->getType() << G4endl;
    it++;
  }

  

  G4cout << "================================================" << G4endl  ;

  G4HCofThisEvent* LHC = evt->GetHCofThisEvent(); //list of Hit collections
  G4int nHC = LHC->GetNumberOfCollections();
  //  G4cout<<"N collections "<<nHC<<G4endl;
  for(G4int iHC=0; iHC<nHC; iHC++) {
    G4String HCname = LHC->GetHC(iHC)->GetName();  //nome della collezione
    //    G4cout << "RootIO: Found hits collection " << HCname << G4endl;
    if (HCname == "ECryCollection") {
      AddECryHits((ExN02ECalHitsCollection*) (LHC->GetHC(iHC)));
    } else if (HCname == "MRodCollection") {
      AddMRodHits((MRodHitsCollection*)(LHC->GetHC(iHC)));
    } else if (HCname == "TRodCollection") {
      AddTRodHits((TRodHitsCollection*)(LHC->GetHC(iHC)));
    } else if (HCname == "EVetoCollection") {
      AddEVetoHits((EVetoHitsCollection*)(LHC->GetHC(iHC)));
    }else if (HCname == "TrackerCollection") {
      AddTrackerHits((TrackerHitsCollection*)(LHC->GetHC(iHC)));
    }
  }
  int Ncells=0;
  for(G4int i=0;i<ECalNRow*ECalNCol;i++){
    if(ETotCry[i]>CellsE) Ncells++;
  }
  //  G4cout<<evt->GetEventID()<<G4endl;
    
  //fill histograms
  //  G4cout << "Proc ID    " <<myStepping->GetPhysProc()<<G4endl; 

  //Retrieve beam Infos!
  G4double BeamE  = myStepping->GetPositronE();
  G4double BeamX  = myStepping->GetPositronX();
  G4double BeamY  = myStepping->GetPositronY();
  G4double BeamZ  = myStepping->GetPositronZ();
  ProcID = myStepping->GetPhysProc();   
  //  G4cout<<Etrack[0]<<G4endl;
  fHistoManager->FillHisto(1, ETotCal);
  fHistoManager->FillHisto(2, BeamE);
  fHistoManager->FillHisto(3, Etrack[0]);
  fHistoManager->FillHisto(4,Ncells);
  fHistoManager->FillHisto(14,ProcID);
  fHistoManager->FillHisto(15,BeamZ);
  fHistoManager->FillHisto2(36,BeamX,BeamY,1.);

  //  G4cout<<"Proc "<<ProcID<<G4endl;
  //  if(NClusters!=0) G4cout<<"Ciao "<<NClusters <<G4endl;
  
  fHistoManager->FillHisto2(30,myStepping->GetGammaEnergy(),myStepping->GetGammaAngle(),1.);
  if(ProcID==1) fHistoManager->FillHisto2(31,myStepping->GetGammaEnergy(),myStepping->GetGammaAngle(),1.);
  if(ProcID==2) fHistoManager->FillHisto2(32,myStepping->GetGammaEnergy(),myStepping->GetGammaAngle(),1.);
  
  if(NClusters>0){
    fHistoManager->FillHisto2(33,myStepping->GetGammaEnergy(),myStepping->GetGammaAngle(),1.);
    if(ProcID==1) fHistoManager->FillHisto2(34,myStepping->GetGammaEnergy(),myStepping->GetGammaAngle(),1.);
    if(ProcID==2) fHistoManager->FillHisto2(35,myStepping->GetGammaEnergy(),myStepping->GetGammaAngle(),1.);
  }
  
  //  fill ntuple
  fHistoManager->myEvt.NTNevent   = evt->GetEventID();
  fHistoManager->myEvt.NTNCluster = NClusters;
  fHistoManager->myEvt.NTNTracks  = NTracks;
  fHistoManager->myEvt.NTEtot     = ETotCal;
  fHistoManager->myEvt.NTPBeam    = BeamE;
  fHistoManager->myEvt.NTXBeam    = BeamX;
  fHistoManager->myEvt.NTYBeam    = BeamY;
  // G4cout<<"ETOT Track"<<ETotTra<<" "<<GetETrack()<<G4endl;
  fHistoManager->myEvt.NTNtotCell = Ncells;
  fHistoManager->myEvt.NTIDProc   = ProcID;
  //  fHistoManager->myEvt.PrimE      = BeamPartE;
  //  G4cout<<"ProcID "<<ProcID<<" "<<GetETrack()<<G4endl;
  for(int i=0;i<NTracks;i++){	  
    if(i>19) break;              
    fHistoManager->myEvt.NTETracker[i] = Etrack[i]*TrackCh[i];
  }
 
  for(int i=0;i<NClusters;i++){	
    if(i>19) break;              
    fHistoManager->myEvt.NTECluster[i] = EneCl[i];
    fHistoManager->myEvt.NTXCluster[i] = XCl[i];   
    fHistoManager->myEvt.NTYCluster[i] = YCl[i];   
    fHistoManager->myEvt.NTTheta[i]    = ThCl[i];  
    fHistoManager->myEvt.NTMmiss2[i]   = MM2[i];  
  }
  for(int i=0;i< ECalNCells;i++){
    fHistoManager->myEvt.NTECell[i]=ETotCry[i];
  }
  //    G4cout<<"Fuck"<<G4endl; 
  // if(ETotCal>0.01) fHistoManager->FillNtuple(&(fHistoManager->myEvt));
  if(ETotCal>EMinSaveNT) fHistoManager->FillNtuple(&(fHistoManager->myEvt));
  //fHistoManager->FillNtuple(&(fHistoManager->myEvt));
}
void ExN02EventAction::AddECryHits(ExN02ECalHitsCollection* hcont)
{
  G4double ETotEvt=0;
  G4int nHits = hcont->entries();
  G4int SeedInd;
  for(G4int jj=0;jj<ECalNCells;jj++){ETotCry[jj]=0.0;}
  
  // Copy all hits from G4 to ROOT event
  for (G4int h=0; h<nHits; h++) {
    ExN02ECalHit* hit = (*hcont)[h]; //prende l'elemento h del vettore hit
    if ( hit != 0 ) {
      ETotCry[hit->GetCryNb()] += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristallo
      ETotEvt += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristalli
      //      CalHitT  = hit->GetTime();
    }
  }//end of loop

  while( (SeedInd=GetSeedCell()) !=-1 ){
    GetEClus(SeedInd);
    double Charge = GetCharge(EClus);
    //  G4cout<<"Charge"<<Charge<<G4endl;
    if (NClusters<10) EneCl[NClusters] = EClus;
    if (NClusters<10) XCl[NClusters]   = ClPosX;
    if (NClusters<10) YCl[NClusters]   = ClPosY;
    if (NClusters<10) ThCl[NClusters]  = Theta;
    if (NClusters<10) MM2[NClusters]   = Mmiss2;
    
    //    G4cout<<"Cl PosX "<<ClPosX<<" Cl PosY "<<ClPosY<<" Nclus "<<EClus<<G4endl;
    fHistoManager->FillHisto2(1,GetCryPosX(SeedInd),GetCryPosY(SeedInd),ETotCry[SeedInd]);  
    fHistoManager->FillHisto2(2,ClPosX,ClPosY,EClus);  
    G4double ProcID = myStepping->GetPhysProc();
 
    if(ProcID==1) fHistoManager->FillHisto2(8,EClus,Theta,1.); //Nclus==1
    if(ProcID==2) fHistoManager->FillHisto2(9,EClus,Theta,1.); //Nclus==2
    NClusters++;
  }
  fHistoManager->FillHisto(11,NClusters);
  ETotCal = ETotEvt;
}

void ExN02EventAction::AddMRodHits(MRodHitsCollection* hcont)
{
  G4double METotRod[100];
  G4double METotRodEvt=0;

  G4int nHits = hcont->entries();
  for(G4int jj=0;jj<100;jj++){METotRod[jj]=0.0;}

  // Copy all hits from G4 to ROOT event
  for (G4int h=0; h<nHits; h++) {
    MRodHit* hit = (*hcont)[h]; //prende l'elemento h del vettore hit
    if ( hit != 0 ) {
      //      ECryHit* rhit = event->AddProfilometerHit();
      //      rhit->SetTrackID(hit->GetTrackID());
      //      rhit->SetFiberNb(hit->GetFiberNb());
      METotRod[hit->GetMRodNb()] += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristallo
      METotRodEvt += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristalli
      //      rhit->SetTime(hit->GetTime());
      //      G4ThreeVector p = hit->GetPos();
      //      rhit->SetPos(TVector3(p.x(),p.y(),p.z()));
      //      G4ThreeVector lp = hit->GetLocPos();
      //      rhit->SetLocPos(TVector3(lp.x(),lp.y(),lp.z()));
    }
  }//end of loop
  EtotMRod = METotRodEvt;    
}

void ExN02EventAction::AddTRodHits(TRodHitsCollection* hcont)
{
  G4double ETotRod[100];
  G4double ETotRodEvt=0;
  G4int nHits = hcont->entries();
  for(G4int jj=0;jj<100;jj++){ETotRod[jj]=0.0;}

  // Copy all hits from G4 to ROOT event
  for (G4int h=0; h<nHits; h++) {
    TRodHit* hit = (*hcont)[h]; //prende l'elemento h del vettore hit
    if ( hit != 0 ) {
      //      ECryHit* rhit = event->AddProfilometerHit();
      //      rhit->SetTrackID(hit->GetTrackID());
      //      rhit->SetFiberNb(hit->GetFiberNb());
      ETotRod[hit->GetTRodNb()] += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristallo
      ETotRodEvt += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristalli
      EtotTRod = ETotRodEvt;
    }
  }//end of loop
  EtotTRod = ETotRodEvt; 
}

void ExN02EventAction::AddTrackerHits(TrackerHitsCollection* hcont)
{				
  const int MaxTracks=10;
  G4double ETotTrackerEvt=0;										
  G4int nHits = hcont->entries();
  //  G4int TrackCh[MaxTracks];
  G4int IDHit[MaxTracks];
  //  G4int Etrack[MaxTracks];
  G4int SeedInd=0;
  G4int LastID=0;
  for(G4int jj=0;jj<200;jj++){ETotRing[jj]=0.0;}
  for(G4int jj=0;jj<MaxTracks;jj++){Etrack[jj]=0.0;}							
  
  for (G4int h=0; h<nHits; h++) {									
    TrackerHit* hit = (*hcont)[h]; //prende l'elemento h del vettore hit					
    if ( hit != 0 ) {											
      ETotRing[hit->GetTrackerNb()]  += hit->GetEdep();
      //*hit->GetTrackCh(); //somma le energie su tutti gli hit di ogni cristal
      ETotTrackerEvt += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristalli	
      if(hit->GetTrackID()!=0 && hit->GetTrackID()!=LastID && hit->GetEdep()>5.) {	
	TrackCh[NTracks] = hit->GetTrackCh();
	IDHit[NTracks]   = hit->GetTrackID();
	Etrack[NTracks]  = hit->GetEdep();
	//	G4cout<<"EtotTracker i "<<Etrack[NTracks]<<" "<<hit->GetTrackCh()<<" "<<nHits<<" "<<hit->GetTrackID()<<G4endl;
	NTracks++;
      }
      if(NTracks>MaxTracks) continue; 
      LastID = hit->GetTrackID();
      // G4cout<<"EtotTracker i "<<hit->GetTrackerNb()<<" "<<ETotRing[hit->GetTrackerNb()]<<" "<<hit->GetTrackCh()<<" "<<nHits<<" "<<hit->GetTrackID()<<G4endl;
      // EtotTra = ETotTrackerEvt;										
    }													
  }//end of loop
  

  // G4cout<<"SI O "<<GetSeedRing()<<G4endl;
//  while( (SeedInd=GetSeedRing()) !=-1 ){
//    G4cout<<"SI "<<SeedInd<<" "<<NTracks<<G4endl;
//    SeedInd=GetSeedRing();
//    if(SeedInd != -1){
//      GetETrack(SeedInd);
//    } 
//  }
//  G4cout<<"Ntr "<<NTracks<<G4endl;
//  G4cout<<"end envent"<<G4endl;
  ETotTra = ETotTrackerEvt;                                                                              
}

void ExN02EventAction::AddEVetoHits(EVetoHitsCollection* hcont)
{
  G4double ETotVeto[100];										
  G4double ETotVetoEvt=0;										
  G4int nHits = hcont->entries();									
  for(G4int jj=0;jj<100;jj++){ETotVeto[jj]=0.0;}							
  													
  for (G4int h=0; h<nHits; h++) {									
    EVetoHit* hit = (*hcont)[h]; //prende l'elemento h del vettore hit					
    if ( hit != 0 ) {											
      ETotVeto[hit->GetEVetoNb()] += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristallo
      ETotVetoEvt += hit->GetEdep(); //somma le energie su tutti gli hit di ogni cristalli		
      //  EtotEVeto = ETotVetoEvt;										
      //      G4cout<<"EtotEVeto "<<ETotFilt[0]<<" "<<hit->GetEdep()<<" Filt Numb "<<hit->GetEVetoNb()<<G4endl;
    }													
  }//end of loop											
  //  G4cout<<"EtotEVetoEvt"<< ETotFiltEvt<< G4endl;							
  EtotEVeto = ETotVetoEvt;                                                                              
}

// Return the index of crystal with maximum energy				  
G4int ExN02EventAction::GetSeedCell()						  
{										  
  G4double Emax=0.;								  
  G4int Imax=-1;								  
  for(G4int i=0;i<ECalNRow*ECalNCol;i++){					  
    if(ETotCry[i]>SeedE && ETotCry[i]>Emax && Used[i]==0) {			  
      Emax=ETotCry[i];								  
      Imax=i;									  
      //      Used[i]=1;							  
    }										  
  }										  
  //  if(Imax !=-1) G4cout<<"I max "<<Imax<<" E max "<<Emax<<" Ncells "<<Ncells<<G4endl;
  //  fHistoManager->FillHisto2(3,Ncells,ETotCal,1.);				  
  return Imax;									  
}                                                                                 

// Return the index of crystal with maximum energy				  
G4int ExN02EventAction::GetSeedRing()						  
{										  
  G4double Emax=0.;								  
  G4int Imax=-1;
  //G4cout<<TrackerNRings<<G4endl;
  for(G4int i=0;i<TrackerNRings;i++){  				  
    if(ETotRing[i]>0 && ETotRing[i]>Emax && UsedRing[i]==0) {			  
      Emax=ETotRing[i];								  
      Imax=i;      	       	       								  
    }										  
  }										  
  return Imax;									  
}                                                                                 

G4double ExN02EventAction::GetETrack(G4int SeedTrack){
  G4double EtotTra=0;
  // G4cout<<"azz" <<SeedTrack <<G4endl;
  if(SeedTrack == -1) return -1;
  //for(G4int kk=10;kk<20;kk++){
  EtotTra+=ETotRing[SeedTrack]+ETotRing[SeedTrack+1]+ETotRing[SeedTrack-1];
  //G4cout<<"EtotTra "<<EtotTra<<" "<<SeedTrack<<G4endl;
  UsedRing[SeedTrack]=1;
  UsedRing[SeedTrack-1]=1;
  UsedRing[SeedTrack+1]=1;

    //}
  return ETotTra;
}

G4double ExN02EventAction::GetEClus(G4int SeedCell)
{
  G4int NcellsCl=0;
  G4int NNeig = GetNeigPos(SeedCell);
  EClus=0;
  ClPosX=0;
  ClPosY=0;
  
  for(G4int j=0;j<NNeig;j++){
    if(ETotCry[Neig[j]]>CellsE && Used[Neig[j]]!=1) {//some fraction of E(SeedCell)
      EClus+=ETotCry[Neig[j]];
      //       G4cout<<"EClus "<<EClus<<G4endl;
      ClPosX+= GetCryPosX(Neig[j])*ETotCry[Neig[j]];
      ClPosY+= GetCryPosY(Neig[j])*ETotCry[Neig[j]];
      //    ClTime+=CryT*ETotCry[i];
      NcellsCl++;
      Used[Neig[j]]=1;
      //    }else if(ETotCry[Neig[j]]<CellsE){  //don't use it in the iterated neig search
      //      Empty[Neig[j]]=1;
    }
  }
  ClPosX=ClPosX/EClus;
  ClPosY=ClPosY/EClus;

  G4double DxDz=ClPosX/(ECalPosiZ-ECalSizeZ/2);
  G4double DyDz=ClPosY/(ECalPosiZ-ECalSizeZ/2);
  G4double norma=sqrt(DxDz*DxDz+DyDz*DyDz+1.);
  G4double GDir[3],BDir[3];
  G4double GMom[4],BeamMom[4],TargetEleMom[4],P4Miss[4];

  GDir[0]=DxDz/norma;
  GDir[1]=DyDz/norma;
  GDir[2]=1./norma;

  GMom[0]=GDir[0]*EClus;  
  GMom[1]=GDir[1]*EClus;
  GMom[2]=GDir[2]*EClus;
  GMom[3]=EClus;        

  BDir[0]=0.;
  BDir[1]=0.;
  BDir[2]=1.;
 
  BeamMom[0]=BDir[0]*BeamEnergy;
  BeamMom[1]=BDir[1]*BeamEnergy;
  BeamMom[2]=BDir[2]*BeamEnergy;
  BeamMom[3]=sqrt(BeamEnergy*BeamEnergy+0.511*0.511);         

  //The target electron is at rest
  TargetEleMom[0]=0.;
  TargetEleMom[1]=0.;
  TargetEleMom[2]=0.;
  TargetEleMom[3]=0.511;

  for (int i=0; i<4; i++){ 
    P4Miss[i]=TargetEleMom[i]+BeamMom[i]-GMom[i];
  }
  Mmiss2 = P4Miss[3]*P4Miss[3]-P4Miss[2]*P4Miss[2]-P4Miss[1]*P4Miss[1]-P4Miss[0]*P4Miss[0];
  //  G4cout<<"***************MMiss "<<Mmiss2<<G4endl;
  //  G4cout<<"***************Gmom "<<GMom[3]<<G4endl;
  //  G4cout<<"***************MMiss "<<Mmiss2<<G4endl;
  G4double product=0;
  for (int i=0; i<3; i++)  product+= GDir[i]*BDir[i];
  ClRadius=sqrt(ClPosX*ClPosX+ClPosY*ClPosY);
  Theta = acos (product) * 180.0 / 3.14159265;  //in degreees
  //  G4cout<<" Scalat product "<<product<<" norma "<<norma<<" Theta "<<Theta<<G4endl;
  // G4cout<<"Eclus out "<<EClus<< " Ecry Seed "<<ETotCry[SeedCell]<<" ClPosX "<<ClPosX<<G4endl;
  //  if(Imax !=-1) G4cout<<"I max "<<Imax<<" E max "<<Emax<<" Ncells "<<Ncells<<G4endl;
  
  fHistoManager->FillHisto(5, NcellsCl);
  fHistoManager->FillHisto(6, EClus);
  fHistoManager->FillHisto(7, SeedCell);
  fHistoManager->FillHisto(8, Theta);
  fHistoManager->FillHisto(9, Mmiss2);
  if(EClus>50. && NcellsCl>1 && ClRadius>3. && EClus<400.) fHistoManager->FillHisto(10,Mmiss2);

  fHistoManager->FillHisto2(4,NcellsCl,EClus,1.);
  fHistoManager->FillHisto2(5,ClPosX,ClPosY,1.);
  fHistoManager->FillHisto2(6,EClus,Theta,1.);
  fHistoManager->FillHisto2(10,Mmiss2,ClRadius,1.);
  if(EClus>50.) fHistoManager->FillHisto2(7,EClus,Theta,1.);
  return EClus;
}

G4double ExN02EventAction::GetCryPosX(G4int CryNb)
{
  G4int Xind=int(CryNb/ECalNRow);
  G4double CryX  = ECalSizeX/ECalNCol;
  G4double Xcoord=-ECalSizeX*0.5+0.5*CryX+Xind*CryX;
  //  G4cout<<"Ncry "<<CryNb<<" Xind " << Xind <<" CryX "<<CryX<<" Xcoord "<<Xcoord<<G4endl;
  return Xcoord;
}

G4double ExN02EventAction::GetCryPosY(G4int CryNb)
{
  G4int Yind=int(CryNb%ECalNRow);
  G4double CryY  = ECalSizeY/ECalNRow;
  G4double Ycoord=-ECalSizeY*0.5+0.5*CryY+Yind*CryY;
  //  G4cout<<"Ncry "<<CryNb<<" Yind " << Yind <<" CryY "<<CryY<<" Ycoord "<<Ycoord<<G4endl;
  return Ycoord;
}

G4double ExN02EventAction::GetCharge(G4double Energia)
{
  G4double LY      = 1800.;
  G4double Gain    =   1E6;
  G4double Echarge = 1.6E-19;
  G4double CollEff =   0.1;
  G4double NPh = Energia*LY;
  G4double QE = 0.2;
  //  G4double rnd = G4RandPoisson::shoot(NPh);
  G4double rnd = G4RandGauss::shoot(NPh,sqrt(NPh));
  //  G4long rnd = G4Poisson(NPh);
  G4double Charge = rnd*CollEff*QE*Gain*Echarge/1E-12;
  //  G4cout<<"Energia "<<Energia<<" rnd " << rnd <<" nph "<<NPh<<" Charge "<<Charge<<" pC"<<G4endl;
  return Charge;
}

G4int ExN02EventAction::GetNeigPos(G4int CryNb)
{
  G4int NNeig=0;
  G4int SeedCell=CryNb;
  
  for(int jj=0;jj<ECalNRow*ECalNCol;jj++) Neig[jj]=0;
  G4double SX= GetCryPosX(SeedCell);
  G4double SY= GetCryPosY(SeedCell);
  
  for(int jj=0;jj<ECalNRow*ECalNCol;jj++){
    G4double NX= GetCryPosX(jj);
    G4double NY= GetCryPosY(jj);
    //    G4cout<<"NX Y "<<NX<<" "<<NY<<" SX Y "<<SX<<" "<<SY<<G4endl;
    //    G4cout<<fabs(NX-SX)<<" "<<RadiusCl<<G4endl;
    if( fabs(NX-SX)<RadiusCl && fabs(NY-SY)<RadiusCl){
      Neig[NNeig]=jj;
      //      G4cout<<" Dx "<< fabs(NX-SX)<<" DY "<< fabs(NY-SY)<<" "<<jj<<" "<<NNeig<<G4endl;
      NNeig++;
    }
  }
  //  G4cout<<" NNeig "<<NNeig<<" SeedCell "<<SeedCell<<G4endl;
  return NNeig;
}
