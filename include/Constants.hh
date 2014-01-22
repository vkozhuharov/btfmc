// GENERAL FLAGS
const int   NPrint      = 10;
const bool  fAutoSeed   = true; //random seed with clock

// DETECTORS FLAGS
const int   IsCalibRun  = 0;
const double   EMinSaveNT  = 0.;

// DETECTORS FLAGS
const double WorldLength = 6.; //in meters
const int   IsFilterON  = 0;
const int   IsEcalON    = 1;
const int   IsMonitorON = 0;
const int   IsTargetON  = 0;
const int   IsMagIronON = 0;
const int   IsWallON    = 0;
const int   IsTrackerON = 0;

// Beam parameter
const int    NPrimaries    = 1; //Maximum 10^5 on the macbook
const int    BeamSpot      = 1;
const double BeamEnergy    = 550.; //in MeV 
const double UMass         = 10.; //in MeV 
const double epsilon       = 1E-3;   //in 10-33
const double SigmaBeamX    = 0.2;
const double SigmaBeamY    = 0.2;

const double MagneticField = 0.8; //in Tesla  0.8

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
const double ECalPosiZ =180.+ECalSizeZ*0.5;	//in cm

const int ECalNRow       = 30;
const int ECalNCol       = 30;
const int ECalNCells     = ECalNRow*ECalNCol;
const double ECalInnHole = 3.; //2.5 cm this is radius

//Cluster reconstruction variables
const double SeedE    = 10.;     //in MeV
const double CellsE   = 0.1;     //in MeV
const double RadiusCl = 4.6;     //in cm due righe due colonne

//  Low Energy Gamma Filter Dimension
const double GFiltPosiX =ECalPosiX;	//in cm
const double GFiltPosiY =ECalPosiY;	//in cm
const double GFiltPosiZ =-0.5*ECalSizeZ+0.5+ECalPosiZ;	//in cm

const double GFiltSizeX =ECalSizeX;   //in cm
const double GFiltSizeY =ECalSizeY;   //in cm
const double GFiltSizeZ =2.5;         //in cm

const double GFiltHoleSizeX =6.;            //in cm
const double GFiltHoleSizeY =6.;            //in cm
const double GFiltHoleSizeZ =GFiltSizeZ+0.1;    //in cm

//const double MagIronX = 168.2;
const double MagIronX = 38.2;
const double MagIronY = 0;
const double MagIronZ = 10;

//  Magnet Dimension
const double MagnetPosiX =0.;	//in cm
const double MagnetPosiY =0.;	//in cm
const double MagnetPosiZ =45.;	//in cm
	     
const double MagnetSizeX =50.;   //in cm
const double MagnetSizeY =50.;   //in cm
const double MagnetSizeZ =50.;   //in cm

// Tracker geometry same position of the magnet

const double TrackerPosiX = 0.;	        //in cm
const double TrackerPosiY = 0.;	        //in cm
const double TrackerPosiZ = 20.;	//in cm

const double TrackerNRings   = 90.;
const double TrackerInnerRad = 18.;
const double TrackerOuterRad = 25.;
const double TrackerHz       = 1.; //questo per due
//const double startAngle = 0.;
//const double spanningAngle = 360.;

//  Beam Dump wall Dimension
const double WallPosiX =0.;	//in cm
const double WallPosiY =0.;	//in cm
const double WallPosiZ =279;	//in cm
	     
const double WallSizeX =200.;   //in cm
const double WallSizeY =300.;   //in cm
const double WallSizeZ =40.;   //in cm
