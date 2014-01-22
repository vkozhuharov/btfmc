// GENERAL FLAGS
const int   NPrint      = 1000;
const bool  fAutoSeed   = true; //random seed with clock

// DETECTORS FLAGS
const int      IsCalibRun  = 0;
const double   EMinSaveNT  = 1.;

// DETECTORS FLAGS
const double WorldLength  = 6.; //in meters
const int   IsEVetoON     = 0;
const int   IsDipoleON    = 1;
const int   IsMagIronON   = 0;  // BTF magnet
const int   IsEcalON      = 1;
const int   IsEcalRoundON = 1;
const int   IsMonitorON   = 0;
const int   IsTargetON    = 1;
const int   IsWallON      = 0;
const int   IsTrackerON   = 1;

// Beam parameter
const int    NPrimaries     = 100000;    //Maximum 10^5 on the macbook
const int    BeamSpot       = 1;
const int    BeamESpread    = 1;
const int    BeamEmitance   = 0;    //not used still
const int    IsSincrotronON = 1;    

const double BeamEnergy    = 550.; //in MeV
const double SigmaBeamE    = 5.5; //in MeV
const double SigmaBeamX    = 0.2;
const double SigmaBeamY    = 0.2;

const double UMass         = 17.5;  //in MeV 
const double epsilon       = 1E-3; //in 10-33

const double MagneticField = 0.6; //in Tesla  0.8 or even 0.4

//  Monitor Dimension
const double MonitorSizeX =4.; //in cm
const double MonitorSizeY =4.; //in cm
const double MonitorSizeZ =4.; //in cm
				
const double MonitorPosiX =0.;	  //in cm
const double MonitorPosiY =0.;	  //in cm
const double MonitorPosiZ =-50.;  //in cm

//  Target Dimension
const double TargetSizeX =4.;	//in cm
const double TargetSizeY =4.;	//in cm
const double TargetSizeZ =0.005;//in cm

const double TargetPosiX =0.;	//in cm
const double TargetPosiY =0.;	//in cm
const double TargetPosiZ =0.;	//in cm

//  ECAL Dimension
const double ECalSizeX =30.;   //in cm
const double ECalSizeY =30.;   //in cm
const double ECalSizeZ =15.;   //in cm

const double ECalPosiX =0.;	//in cm
const double ECalPosiY =0.;	//in cm
const double ECalPosiZ =175.+ECalSizeZ*0.5;	//in cm

const int ECalNRow       = 30;
const int ECalNCol       = 30;
const int ECalNCells     = ECalNRow*ECalNCol;
const double ECalInnHole = 4.5; //2.5 cm this is radius

//Cluster reconstruction variables
const double SeedE    = 10.;     //in MeV
const double CellsE   = 0.1;     //in MeV
const double RadiusCl = 4.6;     //in cm due righe due colonne

//  Low Energy Gamma EVeto Dimension
const double EVetoPosiX = 0.;	//in cm
const double EVetoPosiY =-85.;	//in cm
const double EVetoPosiZ =-0.5*ECalSizeZ+0.5+ECalPosiZ+20;	//in cm
	     
const double EVetoSizeX =10.;   //in cm
const double EVetoSizeY =80.;   //in cm
const double EVetoSizeZ =2.;         //in cm

//const double GFiltHoleSizeX =6.;            //in cm
//const double GFiltHoleSizeY =6.;            //in cm
//const double GFiltHoleSizeZ =GFiltSizeZ+0.1;    //in cm

//const double MagIronX = 168.2;
const double MagIronX = 0;
const double MagIronY = 0;
const double MagIronZ = 20;

// Tracker geometry same position of the magnet

const double TrackerPosiX = 0.;	        //in cm
const double TrackerPosiY = 0.;	        //in cm
const double TrackerPosiZ = 20.;	//in cm

//  Magnet Dimension
const double TrackerNRings   = 100.;
const double TrackerInnerRad = 20.;
const double TrackerOuterRad = 25.;
const double TrackerHz       = 1.; //questo per due
	     
const double MagnetSizeX =TrackerInnerRad*2.;   //in cm
const double MagnetSizeY =TrackerInnerRad*2.;   //in cm
const double MagnetSizeZ =100.;   //in cm

const double MagnetPosiX =0.;	//in cm
const double MagnetPosiY =0.;	//in cm
const double MagnetPosiZ =TrackerPosiZ+MagnetSizeZ*0.5;	//in cm

//const double startAngle = 0.;
//const double spanningAngle = 360.;
const int MaxTracks = 20.;

//  Beam Dump wall Dimension
const double WallPosiX =0.;	//in cm
const double WallPosiY =0.;	//in cm
const double WallPosiZ =279;	//in cm
	     
const double WallSizeX =200.;   //in cm
const double WallSizeY =300.;   //in cm
const double WallSizeZ =40.;   //in cm
