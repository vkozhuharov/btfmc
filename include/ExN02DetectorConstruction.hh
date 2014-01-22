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
// $Id: ExN02DetectorConstruction.hh,v 1.2 2014/01/22 17:25:08 veni Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN02DetectorConstruction_h
#define ExN02DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "ExN02MagneticField.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4PhysicalConstants.hh"



class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPVParameterisation;
class G4UserLimits;
class G4SubtractionSolid;
class ExN02DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     ExN02DetectorConstruction();
    ~ExN02DetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
  //     const
     //     G4VPhysicalVolume* GetTracker() {return physiTracker;};
     //     G4double GetTrackerFullLength() {return fTrackerLength;};
     //     G4double GetTargetFullLength()  {return fTargetLength;};
     G4double GetWorldFullLength()   {return fWorldLength;}; 
     G4double GetCryPosX(G4int CryInd){return CryX[CryInd];}; 
     G4double GetCryPosY(G4int CryInd){return CryY[CryInd];}; 
     G4double GetCryPosZ(G4int CryInd){return CryZ[CryInd];}; 
     
     void setTargetMaterial(G4String);
     void SetupDetectors();
  //     void SetMagField(G4double);
     void SetMaxStep (G4double);     
  
  private:
     void DefineMaterials();

     G4Box*             solidWorld;    
     G4LogicalVolume*   logicWorld;    
     G4VPhysicalVolume* physiWorld;    

     G4Box*             solidWall;    
     G4LogicalVolume*   logicWall;    
     G4VPhysicalVolume* physiWall;    

     G4Box*             solidTarget;   
     G4LogicalVolume*   logicTarget;   
     G4VPhysicalVolume* physiTarget;   

  // G4SubtractionSolid* solidGSubFilt;   
     G4VSolid*           solidEVeto;   
     G4LogicalVolume*    logicEVeto;   
     G4VPhysicalVolume*  physiEVeto;   
 
     G4Box*             solidTXRod;   
     G4LogicalVolume*   logicTXRod;   
     G4VPhysicalVolume* physiTXRod;   

     G4Box*             solidTYRod;   
     G4LogicalVolume*   logicTYRod;   
     G4VPhysicalVolume* physiTYRod;   

     G4Box*             solidMonitor;  
     G4LogicalVolume*   logicMonitor;  
     G4VPhysicalVolume* physiMonitor;  

     G4Box*             solidMXRod;   
     G4LogicalVolume*   logicMXRod;   
     G4VPhysicalVolume* physiMXRod;   

     G4Box*             solidMYRod;   
     G4LogicalVolume*   logicMYRod;   
     G4VPhysicalVolume* physiMYRod;   

     G4Box*             solidEcal;   // pointer to the solid Target
     G4LogicalVolume*   logicEcal;   // pointer to the logical Target
     G4VPhysicalVolume* physiEcal;   // pointer to the physical Target
 
     G4Box*             solidCry;   // pointer to the solid Target
     G4LogicalVolume*   logicCry;   // pointer to the logical Target
     G4VPhysicalVolume* physiCry;   // pointer to the physical Target
 
     G4Tubs*            solidMagIron;
     G4LogicalVolume*   logicMagIron;
     G4VPhysicalVolume* physiMagIron;

     G4Tubs*            solidMagInnJoke;
     G4LogicalVolume*   logicMagInnJoke;
     G4VPhysicalVolume* physiMagInnJoke;

     G4Tubs*            solidMagOutJoke;
     G4LogicalVolume*   logicMagOutJoke;
     G4VPhysicalVolume* physiMagOutJoke;

     G4Tubs*            solidMagBArea;
     G4LogicalVolume*   logicMagBArea;
     G4VPhysicalVolume* physiMagBArea;

     G4Box*             solidSwepMag;  // pointer to the solid Tracker
     G4LogicalVolume*   logicSwepMag;  // pointer to the logical Tracker
     G4VPhysicalVolume* physiSwepMag;  // pointer to the physical Tracker

     G4Tubs*             solidTracker;  // pointer to the solid Chamber
     G4LogicalVolume*   logicTracker;  // pointer to the logical Chamber
     G4VPhysicalVolume* physiTracker;  // pointer to the physical Chamber
     
     G4Material*         TargetMater;  // pointer to the target  material
     G4Material*         WorldMater; // pointer to the chamber material

     G4UserLimits* stepLimit;             // pointer to user step limits

     G4UniformMagField*              fMagField;   // pointer to the magnetic field 
     F03FieldSetup*              fEmFieldSetup;     
     ExN02DetectorMessenger* detectorMessenger;  // pointer to the Messenger
       
     G4double fWorldLength;            // Full length of the world volume

     private:
     G4double CryX[1000];
     G4double CryY[1000];
     G4double CryZ[1000];
     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
