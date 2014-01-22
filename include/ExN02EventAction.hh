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
// $Id: ExN02EventAction.hh,v 1.3 2014/01/22 17:34:55 veni Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef ExN02EventAction_h
#define ExN02EventAction_h 1

#include "G4UserEventAction.hh"
#include "ExN02ECalHit.hh"
#include "MRodHit.hh"
#include "TrackerHit.hh"
#include "TRodHit.hh"
#include "EVetoHit.hh"
#include "GFiltHit.hh"
#include "ExN02DetectorConstruction.hh"
//#include "Constants.hh";

class G4Event;
class ExN02RunAction;
class HistoManager;
class ExN02SteppingAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN02EventAction : public G4UserEventAction
{
  public:
  ExN02EventAction(ExN02RunAction*, HistoManager*);
   ~ExN02EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
  ExN02SteppingAction * myStepping;
  private:
  void  AddECryHits(ExN02ECalHitsCollection*);
  void  AddMRodHits(MRodHitsCollection*);
  void  AddTRodHits(TRodHitsCollection*);
  void  AddEVetoHits(EVetoHitsCollection*);
  void  AddTrackerHits(TrackerHitsCollection*);
  G4double GetCryPosX(G4int CryNb);
  G4double GetCryPosY(G4int CryNb);
  G4int    GetSeedCell();
  G4int    GetSeedRing();
  G4double GetEClus(G4int SeedCell);
  G4double GetETrack(){return 0.;};
  G4double GetETrack(G4int SeedTrack);
  G4double GetCharge(G4double Energia);
  G4int GetNeig(G4int SeedCell);
  G4int GetNeigPos(G4int SeedCell);

  //  G4double Dot(G4double* v1,G4double* v1,G4int N)
  private:
   ExN02RunAction*    fRunAct;
   HistoManager* fHistoManager;
  //che devo fare ce debbo mettere il detector?
   G4double ETotCal;
   G4double EtotTRod; 
   G4double EtotMRod;
   G4double EtotEVeto;
   G4double ProcID;
   G4double ECalHitT,CalEvtT,EtotFiltEvt; 
   G4double ClPosX,ClPosY;
   G4double ClTime,EClus,Theta,ClRadius,Mmiss2,ETotTra;
   G4int NcellsCl,NClusters,NTracks;
   G4double EneCl[20];
   G4double Etrack[20];
   G4int TrackCh[20];
   G4double XCl[20];//For the Ntuple
   G4double YCl[20];//For the Ntuple
   G4double ThCl[20];//For the Ntuple
   G4double MM2[20];//For the Ntuple
   G4double ETotCry[2000]; 
   G4double ETotRing[1000];
   G4int Used[2000];
   G4int UsedRing[2000];	
   G4int Empty[2000];
   G4int Neig[2000];
   G4double CalEvtPosX,CalEvtPosY;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
