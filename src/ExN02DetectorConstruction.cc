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
// $Id: ExN02DetectorConstruction.cc,v 1.1 2014/01/22 15:35:03 veni Exp $
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//#include "InputParam.input"
#include "Constants.hh" 
#include "ExN02DetectorConstruction.hh"
#include "ExN02DetectorMessenger.hh"
#include "ExN02MagneticField.hh"
#include "ExN02ECalSD.hh"
#include "TRodSD.hh"
#include "MRodSD.hh"
#include "TrackerSD.hh"
#include "GFiltSD.hh"
 
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4PSMinKinEAtGeneration.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN02DetectorConstruction::ExN02DetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 solidTarget(0), logicTarget(0), physiTarget(0), 
 TargetMater(0),
 stepLimit(0), fMagField(0), fEmFieldSetup(0), //added M. Raggi
 fWorldLength(0.)
{
  fEmFieldSetup = new F03FieldSetup();
  detectorMessenger = new ExN02DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN02DetectorConstruction::~ExN02DetectorConstruction()
{
  delete fMagField;
  if (fEmFieldSetup) delete fEmFieldSetup ;
  delete stepLimit;
  delete detectorMessenger;             
}

G4VPhysicalVolume* ExN02DetectorConstruction::Construct()
{
//--------- Material definition ---------
  G4double a, z, density;
  G4Material* Vacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,universe_mean_density,kStateGas,2.73*kelvin,3.e-18*pascal);
  //LSO crystals for the calorimeter L2S1O5
  G4Element*  O = new G4Element("Oxygen",    "O",  z=8.,  a=16.00*g/mole);
  G4Element* Lu = new G4Element("Lutetium",  "Lu",  z=71., a=174.97*g/mole);
  G4Element* Si = new G4Element("Silicon",    "Si", z=14., a=28.09*g/mole);
  G4Material* LSO = new G4Material("LSO", density=7.4*g/cm3,3);
  LSO->AddElement(Lu,2);
  LSO->AddElement(Si,1);
  LSO->AddElement(O,5);

  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  G4Material* SiO2   = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* Air    = man->FindOrBuildMaterial("G4_AIR");
  G4Material* Pb     = man->FindOrBuildMaterial("G4_Pb");
  G4Material* PbWO4  = man->FindOrBuildMaterial("G4_PbWO4");
  G4Material* elC    = man->FindOrBuildMaterial("G4_C");
  G4Material* W      = man->FindOrBuildMaterial("G4_W");
  G4Material* Concrete = man->FindOrBuildMaterial("G4_CONCRETE");
  G4Material* Iron     = man->FindOrBuildMaterial("G4_Fe");

  // Print all the materials defined.
  //  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  //------------------------------ 
  // World Volume
  //------------------------------  
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  //  G4cout<<"Computed tolerance = "<<G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm<<" mm" << G4endl;
  WorldMater=Vacuum;
  G4double HalfWorldLength = 0.5*WorldLength*m;

  solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  logicWorld= new G4LogicalVolume( solidWorld, WorldMater, "World", 0, 0, 0);
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
                                 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // copy number

  if(IsWallON==1){  
    G4ThreeVector positionWall = G4ThreeVector(WallPosiX*cm,WallPosiY*cm,WallPosiZ*cm); 
    solidWall = new G4Box("wall",WallSizeX*0.5*cm,WallSizeY*0.5*cm,WallSizeZ*0.5*cm);
    logicWall = new G4LogicalVolume(solidWall,Concrete, "Wall", 0, 0, 0);
    physiWall = new G4PVPlacement(0,               // no rotation
				  positionWall,   // at (0,0,0)
				  logicWall,      // its logical volume
				  "Wall",         // its name
				  logicWorld,     // its mother  volume
				  false,          // no boolean operations
				  0,              // copy number
				  false);         // is overlap    
  }

//  if(IsMagIronON==1){ 
//    G4double innerRad = 131.8*cm;
//    G4double outerRad = 204.6*cm;
//    G4double hz       = 60.*cm;
//    G4double startAngle = 0.*deg;
//    G4double spanningAngle = 45.*deg;
//
//    G4ThreeVector positionMagIron = G4ThreeVector(MagIronX*cm,MagIronY*cm,MagIronZ*cm); 
//
//    solidMagIron = new G4Tubs("magIron",innerRad,outerRad,hz,startAngle,spanningAngle);
//    logicMagIron = new G4LogicalVolume(solidMagIron,Iron, "MagIron", 0, 0, 0);
//    physiMagIron = new G4PVPlacement(0,              // no rotation
//				     positionMagIron,   // at (0,0,0)
//				     logicMagIron,      // its logical volume
//				     "MagIron",         // its name
//				     logicWorld,        // its mother  volume
//				     false,             // no boolean operations
//				     0);                // copy number
//  }

 if(IsTargetON==1){
  //------------------------------------------------- 
  // Target Defintion two layers of fused silica rods 
  //-------------------------------------------------
  G4ThreeVector positionTarget = G4ThreeVector(0,0,0); 

  G4double TargetX      = TargetSizeX*cm;
  G4double TargetY      = TargetSizeY*cm;
  G4double TargetLength = TargetSizeZ*cm;

  solidTarget = new G4Box("target",TargetX*0.5,TargetY*0.5,TargetLength*0.5);
  logicTarget = new G4LogicalVolume(solidTarget,elC,"Target",0,0,0);
  physiTarget = new G4PVPlacement(0,               // no rotation
				  positionTarget,  // at (x,y,z)
				  logicTarget,     // its logical volume                                  
				  "Target",        // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  0,               // copy number 
				  false);          //Check for overlaps
//  //Start Rods Description for Monitor station
//  //Start target Rods Description
//  G4int    NRodRows=5;
//  G4double TXRodLength = 2*mm;
//  G4double TXRodX      = TargetX;
//  G4double TXRodY      = 2*mm;
//
//  G4double TYRodLength = 2*mm;
//  G4double TYRodX      = 2*mm;
//  G4double TYRodY      = TargetY;
//
//  solidTXRod  = new G4Box("TXRod",TXRodX*0.5,TXRodY*0.5,TXRodLength*0.5);
//  logicTXRod  = new G4LogicalVolume(solidTXRod,SiO2,"TXRod");
//
//  solidTYRod  = new G4Box("TYRod",TYRodX*0.5,TYRodY*0.5,TYRodLength*0.5);
//  logicTYRod  = new G4LogicalVolume(solidTYRod,SiO2,"TYRod");
//
//  for (G4int i=0;i<NRodRows;i++){
//    G4ThreeVector positionTXRod = G4ThreeVector(0.   ,-TargetY*0.5+0.5*TXRodY+1*cm+i*TXRodY,+0.*cm);
//    G4ThreeVector positionTYRod = G4ThreeVector(-TargetX*0.5+0.5*TYRodX+i*TYRodX+1*cm,0.,1.*cm);
//    
//    physiTXRod  = new G4PVPlacement(0,             // no rotation
//				    positionTXRod,  // at (x,y,z)
//				    logicTXRod,      // its logical volume                                  
//				    "TXRod",        // its name
//				    logicTarget,     // its mother  volume
//				    false,         // no boolean operations
//				    i);            // copy number 
//    
//    physiTYRod  = new G4PVPlacement(0,               // no rotation
//				    positionTYRod,   // at (x,y,z)
//				    logicTYRod,      // its logical volume                                  
//				    "TYRod",         // its name
//				    logicTarget,     // its mother  volume
//				    false,           // no boolean operations
//				    i+NRodRows);     // copy number 
//  }
 }
 //------------------------------------------------- 
 // Direction monitor fused silica rods 
 //-------------------------------------------------
 if(IsMonitorON==1){
   G4ThreeVector positionMonitor = G4ThreeVector(MonitorPosiX*cm,MonitorPosiY*cm,MonitorPosiZ*cm); 
   G4double MonitorX      = MonitorSizeX*cm;
   G4double MonitorY      = MonitorSizeY*cm;
   G4double MonitorLength = MonitorSizeZ*cm;
   
   solidMonitor = new G4Box("target",MonitorX*0.5,MonitorY*0.5,MonitorLength*0.5);
   logicMonitor = new G4LogicalVolume(solidMonitor,WorldMater,"Monitor",0,0,0);
   physiMonitor = new G4PVPlacement(0,              // no rotation
				    positionMonitor,   // at (x,y,z)
				    logicMonitor,      // its logical volume                                  
				    "Monitor",         // its name
				    logicWorld,        // its mother  volume
				    false,             // no boolean operations
				    0,              // copy number 
				    false);          //Check for overlaps
   
   G4double MXRodLength = 2*mm;
   G4double MXRodX      = MonitorX-1;
   G4double MXRodY      = 2*mm;
   
   G4double MYRodLength = 2*mm;
   G4double MYRodX      = 2*mm;
   G4double MYRodY      = MonitorY-1;
   
   solidMXRod  = new G4Box("MXRod",MXRodX*0.5,MXRodY*0.5,MXRodLength*0.5);
   logicMXRod  = new G4LogicalVolume(solidMXRod,SiO2,"MXRod");
   
   solidMYRod  = new G4Box("MYRod",MYRodX*0.5,MYRodY*0.5,MYRodLength*0.5);
   logicMYRod  = new G4LogicalVolume(solidMYRod,SiO2,"MYRod");
   //
   G4int NMRodRows=10;
   for (G4int i=0;i<NMRodRows;i++){
     G4ThreeVector positionMXRod = G4ThreeVector(0.   ,-MonitorY*0.5+0.5*MXRodY+1*cm+i*MXRodY,+0.*cm);
     G4ThreeVector positionMYRod = G4ThreeVector(-MonitorX*0.5+0.5*MYRodX+i*MYRodX+1*cm,0.,1.*cm);
     
     physiMXRod  = new G4PVPlacement(0,             // no rotation
				     positionMXRod, // at (x,y,z)
				     logicMXRod,    // its logical volume                                  
				     "MXRod",       // its name
				     logicMonitor,  // its mother  volume
				     false,         // no boolean operations
				     i);            // copy number 
     
     physiMYRod  = new G4PVPlacement(0,              // no rotation
				     positionMYRod,  // at (x,y,z)
				     logicMYRod,     // its logical volume                                  
				     "MYRod",        // its name
				     logicMonitor,   // its mother  volume
				     false,          // no boolean operations
				     i+NMRodRows);   // copy number 
   }
 }
 
 if(IsEcalON==1){  
   //------------------------------ 
   // ECal Defintion
   //------------------------------  
   G4ThreeVector positionEcal = G4ThreeVector(ECalPosiX*cm,ECalPosiY*cm,ECalPosiZ*cm); 
   G4double ECalX      = ECalSizeX*cm;
   G4double ECalY      = ECalSizeY*cm;
   G4double ECalLength = ECalSizeZ*cm;

   solidEcal = new G4Box("ECalSolid",ECalX*0.5,ECalY*0.5,ECalLength*0.5);
   logicEcal = new G4LogicalVolume(solidEcal,Air,"ECalLogic",0, 0, 0);
   physiEcal  = new G4PVPlacement(0,             // no rotation
				  positionEcal,  // at (x,y,z)
				  logicEcal,      // its logical volume                                  
				  "ECal",        // its name
				  logicWorld,     // its mother  volume
				  false,         // no boolean operations
				  0,     // copy number 
				  false);          //Check for overlaps 
   
   G4double ECryLength = ECalLength;
   G4double ECryX      = ECalX/ECalNRow;
   G4double ECryY      = ECalY/ECalNCol;
   
   solidCry  = new G4Box("Ecry",ECryX*0.5,ECryY*0.5,ECryLength*0.5);
   logicCry  = new G4LogicalVolume(solidCry,LSO,"ECry",0, 0, 0);
   
   G4int ncry=0;
   for (G4int i=0;i<ECalNRow;i++){
     for (G4int j=0;j<ECalNCol;j++){
       G4double PosXCry=-ECalX*0.5+0.5*ECryX+i*ECryX;
       G4double PosYCry=-ECalY*0.5+0.5*ECryY+j*ECryY;
       G4ThreeVector positionCry = G4ThreeVector(PosXCry,PosYCry,0.);
       G4int HoleFlag=0; //should be parametric in NRow NCol
       
       //       if( fabs(PosXCry)<ECalInnHole*cm && fabs(PosYCry)<ECalInnHole*cm ){
       if( (sqrt( PosXCry*PosXCry + PosYCry*PosYCry ) <ECalInnHole*cm ) || (sqrt( PosXCry*PosXCry + PosYCry*PosYCry ) > ECalSizeX*cm/2.) ){
	 HoleFlag=1;
       }else{
	 HoleFlag=0;
       }
       //      G4cout<<" i "<<i<<" "<<j<<" "<<PosXCry)<<" "<<"InnHolt"<<ECalInnHole<<" "<<HoleFlag<<G4endl;
       //      G4cout<<" i "<<i<<" "<<j<<" "<<PosXCry<<" "<<"Y"<<PosYCry<<" "<<HoleFlag<<G4endl;
       if(HoleFlag!=1){
	 physiCry  = new G4PVPlacement(0,             // no rotation
				       positionCry,  // at (x,y,z)
				       logicCry,      // its logical volume                                  
				       "ECry",        // its name
				       logicEcal,     // its mother  volume
				       false,         // no boolean operations
				       i+j*ECalNCol,     // copy number 
				       false);          //Check for overlaps
	 ncry ++;
       }
     }
   }//end of crystal placements 
   G4cout << "Total number of LYSO crystals:  " << ncry << G4endl;
 }


 //------------------------------ 
 // low Energy gamma filter linked to Calo position
 //------------------------------  
 if(IsFilterON==1){
   G4ThreeVector positionGfilt = G4ThreeVector(GFiltPosiX*cm,GFiltPosiY*cm,GFiltPosiZ*cm); 
    
    solidHole  = new G4Box("hole", GFiltHoleSizeX*cm*0.5,GFiltHoleSizeY*cm*0.5,GFiltHoleSizeZ*cm*0.5);
    solidGfilt = new G4Box("gfilt",GFiltSizeX*cm*0.5,GFiltSizeY*cm*0.5,GFiltSizeZ*cm*0.5);
    solidGSubFilt = new G4SubtractionSolid("GSubFilt",solidGfilt,solidHole,0,G4ThreeVector(0,0,GFiltSizeZ*cm/2.));
    
    logicGfilt = new G4LogicalVolume(solidGSubFilt,Pb,"GSubFilt-H",0,0,0);
    physiGfilt = new G4PVPlacement(0,             // no rotation
				   positionGfilt,  // at (x,y,z)
				   logicGfilt,     // its logical volume    
				   "GFilt-Hole",   // its name
				   logicWorld,     // its mother volume
				   false,          // no boolean operations
				   0,
				   false);       
  }
  //  for(int i=0;i<100;i++) G4cout<<"i "<<i<<" PosX "<<GetCryPosX(i)<<" PosY "<<GetCryPosY(i)<<" PosZ 0 "<<GetCryPosZ(i)<<G4endl;
  
  //------------------------------ 
  // Magnet
  //------------------------------
  G4ThreeVector positionSwepMag = G4ThreeVector(MagnetPosiX*cm,MagnetPosiY*cm,MagnetPosiZ*cm); 
  G4double SwepMagDx=  MagnetSizeX *cm;
  G4double SwepMagDy=  MagnetSizeY *cm;
  G4double SwepMagDz=  MagnetSizeZ *cm;
  G4bool allLocal=true;
  solidSwepMag = new G4Box("swepMag",SwepMagDx*0.5,SwepMagDy*0.5,SwepMagDz*0.5);
  logicSwepMag = new G4LogicalVolume(solidSwepMag,WorldMater,"SwepMag",0,0,0);
  logicSwepMag ->SetFieldManager(fEmFieldSetup->GetLocalFieldManager(),allLocal);
  physiSwepMag = new G4PVPlacement(0,             // no rotation
				   positionSwepMag,  // at (x,y,z)
				   logicSwepMag,     // its logical volume                                  
				   "SwepMag",           // its name
				   logicWorld,       // its mother  volume
				   false,            // no boolean operations
				   0,                // copy number 
				   false);           // Overlap check    

  //------------------------------------------------- 
  // Spectrometer chambers
  //-------------------------------------------------

  if(IsTrackerON==1){
    solidTracker= new G4Tubs("TrackerBox",TrackerInnerRad*cm,TrackerOuterRad*cm,TrackerHz*cm,0.*deg,360.*deg);
    logicTracker= new G4LogicalVolume(solidTracker,Air,"LogTracker",0,0,0);
    for(int kk=0;kk<TrackerNRings;kk++){
      G4ThreeVector positionTracker = G4ThreeVector(TrackerPosiX*cm,TrackerPosiY*cm,TrackerPosiZ*cm+kk*TrackerHz*cm); 
      physiTracker= new G4PVPlacement(0,               // no rotation
				      positionTracker,  // at (x,y,z)
				      logicTracker,     // its logical volume    
				      "Tracker",        // its name
				      logicWorld,       // its mother volume
				      false,            // no boolean operations
				      kk);               // copy number 
    }
  }
 
  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 
  G4SDManager* SDman     = G4SDManager::GetSDMpointer();
  G4String ECrySDname    = "ECrySD";
  G4String TrackerSDname = "TraSD";
  G4String TRodSDname    = "TRodSD";     //target rods
  G4String MRodSDname    = "MRodSD";     //monitor rods
  G4String GFiltSDname   = "GFiltSD";   //Gamma filter
  
  if(IsEcalON==1){
    ExN02ECalSD* ECrySD = new ExN02ECalSD( ECrySDname );
    logicCry->SetSensitiveDetector( ECrySD );
    SDman->AddNewDetector( ECrySD );
  }
  if(IsTrackerON==1){
    TrackerSD* TrackSD = new TrackerSD( TrackerSDname );
    SDman->AddNewDetector( TrackSD );
    logicTracker->SetSensitiveDetector( TrackSD );
  }
  //  TRodSD* TRodSDet = new TRodSD( TRodSDname );
  //  SDman->AddNewDetector( TRodSDet );
  if(IsMonitorON==1){
    MRodSD* MRodSDet = new MRodSD( MRodSDname );
    SDman->AddNewDetector( MRodSDet );
    logicMXRod->SetSensitiveDetector( MRodSDet );
    logicMYRod->SetSensitiveDetector( MRodSDet );
  }
  if(IsFilterON==1){
    GFiltSD* GFiltSDet = new GFiltSD( GFiltSDname );
    SDman->AddNewDetector( GFiltSDet );
    logicGfilt->SetSensitiveDetector( GFiltSDet );
  }
  //  logicTXRod->SetSensitiveDetector( TRodSDet );
  //  logicTYRod->SetSensitiveDetector( TRodSDet );


//--------- Visualization attributes -------------------------------

  logicWorld  ->SetVisAttributes(G4VisAttributes::Invisible);
  //  if(IsTargetON)  logicTarget ->SetVisAttributes(G4VisAttributes::Invisible);
  //  if(IsMonitorON) logicMonitor->SetVisAttributes(G4VisAttributes::Invisible);
  if(IsEcalON)    logicEcal   ->SetVisAttributes(G4VisAttributes::Invisible);
  // logicSwepMag   ->SetVisAttributes(G4VisAttributes::Invisible);
  //  if(IsEcalON)    logicCry    ->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicCry    ->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicGFilt   ->SetVisAttributes(G4VisAttributes::Invisible);

  //  logicTXRod  ->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicTYRod  ->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  //  logicTarget ->SetVisAttributes(BoxVisAtt);

  return physiWorld;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02DetectorConstruction::SetupDetectors()
{
//  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
//  G4String filterName, particleName;
//
//  G4SDParticleFilter* gammaFilter   = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
//  G4SDParticleFilter* electronFilter= new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
//  G4SDParticleFilter* positronFilter= new G4SDParticleFilter(filterName="positronFilter",particleName="e+");
//  G4SDParticleFilter* epFilter      = new G4SDParticleFilter(filterName="epFilter");
//  epFilter->add(particleName="e-");
//  epFilter->add(particleName="e+");
//
//  G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector("SDcalo");
//  G4VPrimitiveScorer* primitive;
//  primitive = new G4PSEnergyDeposit("eDep",0);
//  det->RegisterPrimitive(primitive);
//  primitive = new G4PSNofSecondary("nGamma",0);
//  primitive->SetFilter(gammaFilter);
//  det->RegisterPrimitive(primitive);
//  G4SDManager::GetSDMpointer()->AddNewDetector(det);
//  logicEcal->SetSensitiveDetector(det);
}
 

void ExN02DetectorConstruction::DefineMaterials()
{


}


void ExN02DetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {TargetMater = pttoMaterial;
      logicTarget->SetMaterial(pttoMaterial); 
      //      G4cout << "\n----> The target is " << fTargetLength/cm << " cm of "
      //             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void ExN02DetectorConstruction::setChamberMaterial(G4String materialName)
//{
//  // search the material by its name 
//  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
//  if (pttoMaterial)
//     {ChamberMater = pttoMaterial;
//      logicChamber->SetMaterial(pttoMaterial); 
//      G4cout << "\n----> The chambers are " << ChamberWidth/cm << " cm of "
//             << materialName << G4endl;
//     }             
//}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
//void ExN02DetectorConstruction::SetMagField(G4double fieldValue)
//{
//  MagField->SetMagFieldValue(fieldValue);  
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN02DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
