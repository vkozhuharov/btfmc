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
// $Id: exampleN02.cc,v 1.2 2014/01/22 16:12:54 veni Exp $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02RunAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02SteppingAction.hh"
#include "ExN02SteppingVerbose.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
//#include "global.hh"


#include "MyEvent.hh"



#ifdef  G4MULTITHREADED
#include "G4MTHepRandom.hh"
#else
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
#include <signal.h>

void sighandler(int sig){
    G4cerr << G4endl << "********************************************************************************" << G4endl;
    G4cerr << "Killed with Signal " << sig << G4endl << "Closing ROOT files ..." << G4endl; 
    G4RunManager::GetRunManager()->AbortRun();
    //    histo->save();
    //    RootIOManager::GetInstance()->Close();
    G4cerr << "... Done" << G4endl;
    G4cerr << G4endl << "********************************************************************************" << G4endl;
    //  G4GeometryManager::GetInstance()->OpenGeometry();

#ifdef G4VIS_USE
    delete G4VVisManager::GetConcreteInstance();
#endif
    delete G4RunManager::GetRunManager();
    exit(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  //  int SeedNum = -1;
  //  if (argc == 3) SeedNum = atoi(argv[2]);
  
  // User Verbose output class
  G4VSteppingVerbose* verbosity = new ExN02SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  //choose the Random engine
  //  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  //G4MTHepRandom::setTheEngine(new CLHEP::Ranlux64Engine);  
#ifdef  G4MULTITHREADED
  G4MTHepRandom::setTheEngine(new CLHEP::RanecuEngine);
#else 
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
#endif

  //Get the event storage. Should be instantiated before the HISTO!
  MyEvent *TheEvent = MyEvent::GetInstance(); //



  // set an HistoManager for Root IO
  HistoManager*  histo = new HistoManager();
  // Run manager
  G4RunManager * runManager = new G4RunManager;
  
  // User Initialization classes (mandatory)
  //
  ExN02DetectorConstruction* detector = new ExN02DetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  //  G4VUserPhysicsList* physics = new ExN02PhysicsList;
  G4VUserPhysicsList* physics = new PhysicsList;
  runManager->SetUserInitialization(physics);
 
  // User Action classes
  //
  G4cout << "MAIN: Constructing the generator action" << G4endl;
  G4VUserPrimaryGeneratorAction* gen_action = new ExN02PrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);
  //
  //  G4UserRunAction* run_action = new ExN02RunAction(histo);
  ExN02RunAction* run_action = new ExN02RunAction(histo);
  runManager->SetUserAction(run_action); 
  //
  ExN02EventAction* myevent_action = new ExN02EventAction(run_action,histo);
  G4UserEventAction* event_action  = myevent_action;
  //  ExN02EventAction* event_action = new ExN02EventAction(run_action,histo);
  runManager->SetUserAction(event_action);
  //
  ExN02SteppingAction* mystepping_action = new ExN02SteppingAction;
  G4UserSteppingAction* stepping_action = mystepping_action ;
  runManager->SetUserAction(stepping_action);
  myevent_action->myStepping = mystepping_action; //setting the pointer to stepping action into Event action
  // Initialize G4 kernel
  //
  runManager->Initialize();

  // trap signals to close datafiles in case of abnormal termination
//  void (*prev_fn)(int);
//  signal(SIGUSR1,sighandler);
//  signal(SIGXCPU,sighandler);
//  signal(SIGINT,sighandler);
//  signal(SIGTERM,sighandler);
//  signal(127,sighandler);

//  if (prev_fn==SIGTERM) run_action::EndOfRunAction();
      
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif    
     
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  
  
  if (argc!=1)   // batch mode  
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else           // interactive mode : define UI session
    { 
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
     
#ifdef G4VIS_USE
      delete visManager;
#endif     
    }
  
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  delete runManager;
  delete verbosity;
  delete histo;

  //delete TheEvent;
  return 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
