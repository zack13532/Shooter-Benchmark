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
// $Id: shooter.cc,v 1.1 2007-10-11 13:01:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// shooter - perform test shots.
//

/*
  gmake CXXFLAGS="-g -pg -a -lc_p " CPPVERBOSE=1 for the library

  gmake CXXFLAGS="-g -pg -a -lc_p -static " CPPVERBOSE=1 G4TARGET=shooter
  shooter
  ((line by line profiling in detail)
  gprof -i -p -q -x -A -l `which shooter` > profile.shooter
  (normal profiling)
  gprof -p `which shooter` > profile.shooter
*/

#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <limits>

#include "G4ios.hh"

#include "BuildBoxWorld.hh"
#include "BuildCalorimeter.hh"
#include "Shoot.hh"

#include "G4Timer.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"

#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"

#include "voxeldefs.hh"
#include "Randomize.hh"

G4int numShoot ;

G4double epsilon = 0.0001; //used when comparing doubles

const G4bool optimise= true;
const G4double x0=1.12343*cm;

G4VPhysicalVolume* BuildReplicaCal(G4Material* Air)
{
  G4double xr=x0/4;
  G4int nr=20;
  G4int nl=20;

  G4double r1=20*xr,r2=21*xr,z1=10*x0,z2=11*x0;
  // Container
  G4Tubs *ecalTube=new G4Tubs("ECAL",0,r2,z2,0,360*deg);
  // End cap
  G4Tubs *leakTube=new G4Tubs("LEAK",0,r2,0.5*x0,0,360*deg);
  // Wrapper
  G4Tubs *latrTube=new G4Tubs("LATR",r1,r2,z1,0,360*deg);
  // Main calorimeter block
  G4Tubs *blocTube=new G4Tubs("BLOC",0,r1,z1,0,360*deg);
  G4Tubs *blocTubeR=new G4Tubs("BLOCR",0,r1/nr,z1,0,360*deg);
  G4Tubs *blocTubeRZ=new G4Tubs("BLOCRZ",0,r1/nr,z1/nl,0,360*deg);

  G4double zc=0.5*(z1+z2);

  G4LogicalVolume *ecalLog=new G4LogicalVolume(ecalTube,Air,
					       "ecalLog",0,0,0);
  G4LogicalVolume *leakLog=new G4LogicalVolume(leakTube,Air,
					       "leakLog",0,0,0);
  G4LogicalVolume *latrLog=new G4LogicalVolume(latrTube,Air,
					       "latrLog",0,0,0);
  G4LogicalVolume *blocLog=new G4LogicalVolume(blocTube,Air,
					       "blocLog",0,0,0);
  G4LogicalVolume *blocRLog=new G4LogicalVolume(blocTubeR,Air,
					       "blocRLog",0,0,0);
  G4LogicalVolume *blocRZLog=new G4LogicalVolume(blocTubeRZ,Air,
						 "blocRZLog",0,0,0);


  G4PVPlacement *ecalPhys=new G4PVPlacement(0,G4ThreeVector(),
					    "ecalPhys",
					    ecalLog,
					    0,false,0);
  // Position end caps, wrapper and calorimeter bloc within ecal
  // G4PVPlacement *leakPhys=
  new G4PVPlacement(0,G4ThreeVector(0,0,-zc),
		    "leakPhys",
		    leakLog,
		    ecalPhys,false,0);
  // G4PVPlacement *leakPhys2=
  new G4PVPlacement(0,G4ThreeVector(0,0,zc),
		    "leakPhys",
		    leakLog,
		    ecalPhys,false,1);
  // G4PVPlacement *latrPhys=
  new G4PVPlacement(0,G4ThreeVector(),
		    "latrPhys",
		    latrLog,
		    ecalPhys,false,0);
  G4PVPlacement *blocPhys=new G4PVPlacement(0,G4ThreeVector(),
					    "blocPhys",
					    blocLog,
					    ecalPhys,false,0);

  // Create replicas
  G4PVReplica *blocPhysR=new G4PVReplica("blocRepR",
					 blocRLog,blocPhys,
					 kRho,nr,r1/nr,0);
  // G4PVReplica *blocPhysRZ=
  new G4PVReplica("blocRepRZ",
		  blocRZLog,blocPhysR,
		  kZAxis,nl,2*z1/nl);

  return ecalPhys;
}


/*
   to change accuracy :
	fieldMgr->GetChordFinder()->SetDeltaChord( G4double newValue);
*/

int main(int argc,char *argv[])
{
  G4ThreeVector origin(0,0,0),pMX(-500,0,0);
  G4ThreeVector vx(1,0,0);
  G4ThreeVector vy(0,1,0);
  G4ThreeVector source(0,0,0);
  G4double xloc = 0, yloc = 0, zloc = 0;

  G4ThreeVector shootdir(0,0,0);

  G4VPhysicalVolume *myTopNode=0;

  G4double Field = 0.*tesla;

  kVoxelVolumeLimits* kvvl = kVoxelVolumeLimits::GetInstance();

  //set up randomization engine for particle generation
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  //set the seeds for particle generation
  long* seeds = (long*)calloc(2, sizeof(long));
  seeds[0] = 100;
  seeds[1] = 200;


  //------------------------------------------------------

  enum GeomType { BOX, CALOLOOP, CALOREP, CALOVARYLOOP} GeomType;
  G4int i;
  G4bool useField = false;

  /* Default Value */

  //starting, ending, and incrementing values for voxelization parameters
  //defaults if user does not enter them as command line arguments
  G4int numruns=0;
  G4int maxnodesinit = kvvl->GetMaxNodes(), maxnodesend = maxnodesinit, maxnodesincr = maxnodesinit;
  G4int minnodeslv1init = kvvl->GetMinVolsLevel1(), minnodeslv1end = minnodeslv1init, minnodeslv1incr = minnodeslv1init;
  G4int minnodeslv2init = kvvl->GetMinVolsLevel2(), minnodeslv2end = minnodeslv2init, minnodeslv2incr = minnodeslv2init;
  G4double ratioinit = kvvl->GetMinVolsLevel3(), ratioend = ratioinit, ratioincr = ratioinit;
  G4int orglimit = 0;

  G4bool randomdirs = false;

  GeomType = CALOLOOP ;
  numShoot = 1000000 ;
  numruns = 10;
  /* Command line parsing */
  
  for (i=1;i<argc;i++) {
    if ((i < (argc-1)) && (strcmp(argv[i],"-event") == 0)) {
      sscanf(argv[i+1],"%d",&numShoot);      
    }
    else if ((i < (argc-1)) && (strcmp(argv[i],"-geom") == 0)) {
      if (strcmp(argv[i+1],"box") == 0) {
	GeomType = BOX;
      } else    
      if (strcmp(argv[i+1],"caloloop") == 0) {
	GeomType = CALOLOOP;
      } else      
      if (strcmp(argv[i+1],"calorep") == 0) {
	GeomType = CALOREP;
      } else
      if (strcmp(argv[i+1],"calovaryloop") == 0) {
	GeomType = CALOVARYLOOP;
      } else {
	G4cerr << argv[i+1] << " is not a known geometry (box,caloloop,calorep,calovaryloop)" << G4endl ;
	exit (0);
	
      }
    }
    else if ((i < (argc-1)) && (strcmp(argv[i],"-magn") == 0)) {
      sscanf(argv[i+1],"%lf",&Field);
      G4cout << " Mag Field = " << Field << G4endl ;
      
      Field = Field * tesla ;
      useField = true;
    }
    else if ((i < (argc-3)) && (strcmp(argv[i],"-ratio") == 0)) {
      sscanf(argv[i+1],"%lf",&ratioinit);
      sscanf(argv[i+2],"%lf",&ratioend);
      sscanf(argv[i+3],"%lf",&ratioincr);
      G4cout << " Ratio range is: [" << ratioinit << ".." << ratioend << ".." << ratioincr << "]" << G4endl ;
    }
    else if ((i < (argc-3)) && (strcmp(argv[i],"-maxnodes") == 0)) {
      sscanf(argv[i+1],"%d",&maxnodesinit);
      sscanf(argv[i+2],"%d",&maxnodesend);
      sscanf(argv[i+3],"%d",&maxnodesincr);
      G4cout << " Max Nodes range is: [" << maxnodesinit << ".." << maxnodesend << ".." << maxnodesincr << "]" << G4endl ;
    }
    else if ((i < (argc-3)) && (strcmp(argv[i],"-minnodeslv1") == 0)) {
      sscanf(argv[i+1],"%d",&minnodeslv1init);
      sscanf(argv[i+2],"%d",&minnodeslv1end);
      sscanf(argv[i+3],"%d",&minnodeslv1incr);
      G4cout << " Min Nodes Lv. 1 range is: [" << minnodeslv1init << ".." << minnodeslv1end << ".." << minnodeslv1incr << "]" << G4endl ;
    }
    else if ((i < (argc-3)) && (strcmp(argv[i],"-minnodeslv2") == 0)) {
      sscanf(argv[i+1],"%d",&minnodeslv2init);
      sscanf(argv[i+2],"%d",&minnodeslv2end);
      sscanf(argv[i+3],"%d",&minnodeslv2incr);
      G4cout << " Min Nodes Lv. 2 range is: [" << minnodeslv2init << ".." << minnodeslv2end << ".." << minnodeslv2incr << "]" << G4endl ;
    }
    else if ((i < (argc-2)) && (strcmp(argv[i],"-randomize") == 0)) {
      randomdirs = true;
      sscanf(argv[i+1],"%ld",&seeds[0]);
      sscanf(argv[i+2],"%ld",&seeds[1]);
      G4cout << " Shooting direction will be random." << G4endl ;
    }
    else if ((i < (argc-1)) && (strcmp(argv[i],"-dirs") == 0)) {
      sscanf(argv[i+1],"%d",&numruns);
      G4cout << " Number of directions to shoot in = " << numruns << G4endl ;
    }
    else if ((i < (argc)) && (strcmp(argv[i],"-testorg") == 0)) {
      orglimit = 1;
      G4cout << " Testing with organization off and on" << G4endl ;
    }
  }
  
  
  //output introduction and description of run-time arguments
  G4cout << "***  Navigation Performance Tester - E.Medernach 30.10.00  ***" << G4endl
	 << ">>>  Based on original benchmark test by P.Kent" << G4endl << G4endl;

  G4cout << "Options (as arguments):" << G4endl
         << "-event <number_of_events>" << G4endl
         << "      number of events for the test. Default is 1000000" << G4endl
         << "-dirs <number_of_directions>" << G4endl
         << "      number of directions to shoot in for the test. Default is 10" << G4endl
         << "      if randomization is not on, then the actual number of runs per set is given by:" << G4endl
         << "      #dirs * #dirs - #dirs + 1. this samples angles in dphi and dcos(theta), avoiding double counting" << G4endl
         << "-randomize <seed 1> <seed 2>" << G4endl
         << "      shoot in random directions, with the engine seeded by the required inputs" << G4endl
         << "-testorg" << G4endl
         << "      if this option is specified, run all sets a second time with daughter organization on" << G4endl
         << "-geom <geometry_type>" << G4endl
         << "      where <geometry_type> can be:" << G4endl
         << "      box      - simple box (default)" << G4endl
         << "      caloloop - calorimeter made by a loop of placements" << G4endl
         << "      calorep  - calorimeter made of replicas" << G4endl
		 << "      calovaryloop – calorimeter made by a loop of placements of several different solids" << G4endl
         << "-magn <magnetic_field_value>" << G4endl
         << "      activates magnetic field (value in tesla units). Default is OFF" << G4endl
         << "-maxnodes <max_nodes_initial> <max_nodes_final> <max_nodes_increment>" << G4endl
         << "      runs a set of simulations for Max Nodes set to values" << G4endl
         << "      in the given range, incremented by the given amount" << G4endl
         << "-minnodeslv1 <min_nodes_lv1_initial> <min_nodes_lv1_final> <min_nodes_lv1_increment>" << G4endl
         << "      see -maxnodes"  << G4endl
         << "-minnodeslv2 <min_nodes_lv2_initial> <min_nodes_lv2_final> <min_nodes_lv2_increment>" << G4endl
         << "      see -maxnodes"  << G4endl
         << "-minnodeslv3 <min_nodes_lv3_initial> <min_nodes_lv3_final> <min_nodes_lv3_increment>" << G4endl
         << "      see -maxnodes"  << G4endl
         << "-ratio <ratio_intial> <ratio_final> <ratio_increment>" << G4endl
         << "      see -maxnodes"  << G4endl
         << G4endl;

  // Build the geometry
  G4cout << "Geometry type:";
  switch (GeomType) {

  case BOX:
    G4cout << " Box only." << G4endl ;
    myTopNode=BuildBoxWorld();
    break;

  case CALOLOOP:
    G4cout << " Calorimeter made of placements." << G4endl ;
    myTopNode=BuildCalorimeter();
    break;

  case CALOREP:
    G4cout << "  Calorimeter made of replicas." << G4endl ;
    myTopNode=BuildReplicaCal(0);
    break;
    
  case CALOVARYLOOP:
	G4cout << " Calorimeter made of various different placements." << G4endl;
	myTopNode=BuildVaryCalorimeter();
	break;

  }
  
  if(!myTopNode)
	  exit(1);

  G4GeometryManager* geoman = G4GeometryManager::GetInstance();

  if (!useField)
  {
    G4cout << "--> Magnetic Field is disabled !" << G4endl ;

    //Setting up variables for runs
    //----------------------------------------------------------
    
    G4double avgtime, errtime;

    G4int actualnumruns = randomdirs ? numruns : numruns * numruns - numruns + 1;
    //array of simulation times for a run with a single set of parameters
    G4double* times = (G4double*)malloc(actualnumruns*sizeof(G4double));

    //----------------------------------------------------------
    //Simulation loops - tests ranges of voxelization parameters
    //----------------------------------------------------------

    //loop over having organization on and off
    for(G4int org = 0 ; org <= orglimit ; org++){
      kvvl->SetOrganization(org);
    	
    //loop over Max Nodes values
	for(G4int maxnodes = maxnodesinit; maxnodes <= maxnodesend ; maxnodes+=maxnodesincr){
      kvvl->SetMaxNodes(maxnodes);
    
    //loop over Min Nodes Lv. 1 values
	for(G4int minnodeslv1 = minnodeslv1init; minnodeslv1 <= minnodeslv1end ; minnodeslv1+=minnodeslv1incr){
      kvvl->SetMinVolsLevel1(minnodeslv1);
    	
    //loop over Min Nodes Lv. 2 values
    for(G4int minnodeslv2 = minnodeslv2init; minnodeslv2 <= minnodeslv2end ; minnodeslv2+=minnodeslv2incr){
      kvvl->SetMinVolsLevel2(minnodeslv2);

    //loop over Ratio values
	for(G4double ratio = ratioinit; ratio < ratioend + epsilon ; ratio+=ratioincr){
      kvvl->SetDefaultSmartless(ratio);
    	
      	// Resets voxels to apply the voxel parameters
		geoman->OpenGeometry();
    	geoman->CloseGeometry(true);

		CLHEP::HepRandom::setTheSeeds(seeds);
		avgtime=0;
		errtime=0;

		G4double jlim = (randomdirs ? 1.:numruns );
		for(int j = 0 ; j < jlim ; j++){

		G4int ilim = (!randomdirs && !j) ? 1:numruns;
		for(int i = 0 ; i < ilim ; i++){ //runs one time w/ i=0 if randomdirs are on

			G4double phi, costheta, sintheta;

			//randomly generated particle source
			xloc = 1800 * CLHEP::mm * (G4UniformRand() - .5);
			yloc = 2000 * CLHEP::mm * (G4UniformRand() - .5);
			zloc = 2100 * CLHEP::mm * (G4UniformRand() - .5);
			source = G4ThreeVector(xloc, yloc, zloc);

			if(randomdirs)
			{
				costheta = 2 * G4UniformRand() - 1.; //generates z component of direction between -1 and 1
				phi = .00011 * zloc + 1.3 + G4UniformRand() * .4; //generates polar angle in a rhombus skewed in z-axis
			}
			else
			{
				costheta = 2.*i/(G4double)numruns - 1.;
				phi = 2. * j * CLHEP::pi/(G4double)numruns;
			}

			sintheta = sqrt(1. - costheta * costheta);

			shootdir = G4ThreeVector(sintheta*cos(phi),
									 sintheta*sin(phi),
									 costheta);
			G4cout << "Shooting from " << source << " along " << shootdir << G4endl;
			times[i] = Shoot(numShoot,myTopNode,source,shootdir);
			avgtime+=times[i];
		}}
	
		//calculating statistics
		avgtime/=actualnumruns;
		for(int i = 0 ; i < actualnumruns ; i++)
			errtime+=pow(times[i]-avgtime,2);
		errtime = sqrt(errtime/(actualnumruns-1));

		G4cout << "~~~Org on?:" << org << ":Ratio:" << ratio << ":Max Nodes:" << maxnodes
			   << ":Min Nodes Lv1:" << minnodeslv1 << ":Result for:" << numShoot
			   << ":Particles:" << avgtime << ":wError:" << errtime
			   << ":and average time per particle:" << avgtime/numShoot
			   << ":for # directions:" << actualnumruns << G4endl;
	}
    }
	}
	}
    }
  }
  else
  {
    //    Field = 0.1 * tesla ;
    G4double DeltaChord = 1.0e-2 * mm ;
    
    G4cout << "--> Magnetic Field test with " << Field/tesla << " Tesla" << G4endl ;

    G4cout << G4endl << "Shooting from " << origin << " along " << vx << G4endl;
    ShootVerbose(myTopNode,origin,vx);
    MagneticShoot(numShoot,myTopNode,origin,vx,Field,DeltaChord);

    G4cout << G4endl << "Shooting from " << origin << " along " << vy << G4endl;
    ShootVerbose(myTopNode,origin,vy);
    MagneticShoot(numShoot,myTopNode,origin,vy,Field,DeltaChord);

    G4cout << G4endl << "Shooting from " << pMX << " along " << vx << G4endl;
    ShootVerbose(myTopNode,pMX,vx);
    MagneticShoot(numShoot,myTopNode,pMX,vx,Field,DeltaChord);
  }

  G4GeometryManager::GetInstance()->OpenGeometry();
  return EXIT_SUCCESS;
}

