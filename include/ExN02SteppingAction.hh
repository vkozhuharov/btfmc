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
// $Id: ExN02SteppingAction.hh,v 1.2 2014/01/22 17:29:54 veni Exp $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN02SteppingAction_h
#define ExN02SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN02SteppingAction : public G4UserSteppingAction
{
  public:
    ExN02SteppingAction();
   ~ExN02SteppingAction(){};

  void UserSteppingAction(const G4Step*);
  double GetPhysProc()   {return ProcID;};
  void   SetPhysProc(float value);

  double GetPositronE()  {return BeamPartE;};
  double GetPositronX()  {return VertexPos[0];};
  double GetPositronY()  {return VertexPos[1];};
  double GetPositronZ()  {return VertexPos[2];};

  double GetGammaEnergy(){return GammaE;};
  double SetGammaAngle(G4ThreeVector GammaDir,G4ThreeVector BeamDir);
  double GetGammaAngle(){return ThetaGamma;};

private:
  double ProcID;
  double PositronE;
  double ElectronE;	
  double BeamPartE;
  double GammaE;
  double ThetaGamma;
  G4ThreeVector GammaDir;
  G4ThreeVector VertexPos;
  G4int NChild;
  //   HistoManager* fHistoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
