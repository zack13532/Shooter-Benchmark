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
// $Id: BuildCalorimeter.hh,v 1.1 2007-10-11 13:01:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef BUILDCALORIMETER_HH
#define BUILDCALORIMETER_HH

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Para.hh"
#include "G4Cons.hh"

//#define SHOOTERGEOMETRYDEBUG

G4VPhysicalVolume* BuildCalorimeter()
{
    G4double offset=22.5,xTlate,yTlate;
    G4int i,j,copyNo;

    G4Box *myWorldBox= new G4Box ("WBox",10000.*CLHEP::mm,10000.*CLHEP::mm,10000.*CLHEP::mm);
    G4Box *myCalBox = new G4Box ("CBox",1500.*CLHEP::mm,1500.*CLHEP::mm,1000.*CLHEP::mm);
    G4Tubs *myTargetTube = new G4Tubs ("TTube",0.,22.5*CLHEP::mm,1000.*CLHEP::mm,0.,360.*CLHEP::deg);

    G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,0,
						   "WLog",0,0,0);
    G4LogicalVolume *myCalLog=new G4LogicalVolume(myCalBox,0,
						  "CLog",0,0,0);
    G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetTube,0,
						     "TLog",0,0,0);

    G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
    G4PVPlacement *myCalPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "CalPhys",
					       myCalLog,
					       myWorldPhys,false,0);

    G4String tName1("TPhys1");	// Allow all target physicals to share
				// same name (delayed copy)
    copyNo=0;
    for (j=1;j<=25;j++)
      {
	yTlate=-1000.0-40.0+j*80.0;
	
	for (i=1;i<=50;i++)
	  {
	    copyNo++;
	    xTlate=-1000.0-20.0+i*45.0-offset;
	    // G4PVPlacement *myTargetPhys=
	    new G4PVPlacement(0,G4ThreeVector(xTlate*CLHEP::mm,yTlate*CLHEP::mm,0),
					      tName1,
					      myTargetLog,
					      myCalPhys,
					      false,
					      copyNo);
	  }
      }
    
    G4String tName2("TPhys2");	// Allow all target physicals to share
    // same name (delayed copy)
    copyNo=0;
    for (j=1;j<=26;j++)
      {
	yTlate=-1000.0-80.0+j*80.0;
	for (i=1;i<=50;i++)
	  {
	    copyNo++;
	    xTlate=-1000.0-20.0+i*45.0;
	    // G4PVPlacement *myTargetPhys=
            new G4PVPlacement(0,G4ThreeVector(xTlate*CLHEP::mm,yTlate*CLHEP::mm,0),
					      tName2,
					      myTargetLog,
					      myCalPhys,
					      false,
					      copyNo);
	  }
      }

    return myWorldPhys;
}

G4VPhysicalVolume* BuildVaryCalorimeter()
{
    G4double offset=22.5,xTlate,yTlate;
    G4int i,j,copyNo;

    G4Box *myWorldBox= new G4Box ("WBox",10000.*CLHEP::mm,10000.*CLHEP::mm,10000.*CLHEP::mm);
    G4Box *myCalBox = new G4Box ("CBox",1500.*CLHEP::mm,1500.*CLHEP::mm,1000.*CLHEP::mm);
    G4Tubs *myTargetTube = new G4Tubs ("TTube",0.,22.5*CLHEP::mm,1000.*CLHEP::mm,0.,360.*CLHEP::deg);
    G4Cons *myTargetConeSection = new G4Cons ("TCone",0.,0.,22.*CLHEP::mm,15.*CLHEP::mm,1000.*CLHEP::mm,0.,360.*CLHEP::deg);
    G4Box *myTargetBox = new G4Box ("TBox",22.*CLHEP::mm,22.*CLHEP::mm,1000.*CLHEP::mm);
    G4Para *myTargetPara = new G4Para ("TPara",15.*CLHEP::mm,15.*CLHEP::mm,800.*CLHEP::mm,15.*CLHEP::deg,15.*CLHEP::deg,15.*CLHEP::deg);

    G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,0,
						   "WLog",0,0,0);
    G4LogicalVolume *myCalLog=new G4LogicalVolume(myCalBox,0,
						  "CLog",0,0,0);

    G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetTube,0,"TTubeLog",0,0,0);
    G4LogicalVolume *myTargetHalfCone=new G4LogicalVolume(myTargetConeSection,0,"TConeLog",0,0,0);
    G4LogicalVolume *myTargetRPrism=new G4LogicalVolume(myTargetBox,0,"TBoxLog",0,0,0);
    G4LogicalVolume *myTargetPPed=new G4LogicalVolume(myTargetPara,0,"TParaLog",0,0,0);

    G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
    G4PVPlacement *myCalPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "CalPhys",
					       myCalLog,
					       myWorldPhys,false,0);

    G4PVPlacement* placedTarget = NULL;

    G4String tName1;

    copyNo=0;
    for (j=1;j<=25;j++)
      {
	yTlate=-1000.0-40.0+j*80.0;

	for (i=1;i<=50;i++)
	  {
	    copyNo++;
	    xTlate=-1000.0-20.0+i*45.0-offset;

	    G4LogicalVolume *myTargetLogical=NULL;
	    switch((j*50+i+1)%4){

			case 1 :
			case 2 : myTargetLogical = myTargetRPrism; tName1 = G4String("TBoxPhys1"); break;
			case 0 :
			default : myTargetLogical = myTargetLog; tName1 = G4String("TTubePhys1"); break;

	    }
	    placedTarget = new G4PVPlacement(0,G4ThreeVector(xTlate*CLHEP::mm,yTlate*CLHEP::mm,0),
					      tName1,
					      myTargetLogical,
					      myCalPhys,
					      false,
					      copyNo);
#ifdef SHOOTERGEOMETRYDEBUG
	    G4cout<<"Placement of a " << tName1 << " : (" << xTlate << "," << yTlate << ",0)" << G4endl;

	    if(placedTarget->CheckOverlaps()){G4cout<<"Overlapping in first loop!" <<G4endl; return NULL;}
#endif

	  }
      }

    G4String tName2;

    copyNo=0;
    for (j=1;j<=26;j++)
      {
	yTlate=-1000.0-80.0+j*80.0;
	for (i=1;i<=50;i++)
	  {
	    copyNo++;
	    xTlate=-1000.0-20.0+i*45.0;

	    G4LogicalVolume *myTargetLogical=NULL;
		switch((j*50+i)%4){

			case 1 : myTargetLogical = myTargetHalfCone; tName2 = G4String("TPhysHalfCone2"); break;
			case 2 : myTargetLogical = myTargetRPrism; tName2 = G4String("TBoxPhys2"); break;
			case 3 : myTargetLogical = myTargetPPed; tName2 = G4String("TPhysPPed2"); break;
			case 0 :
			default : myTargetLogical = myTargetLog; tName2 = G4String("TTubePhys2"); break;

		}
		placedTarget = new G4PVPlacement(0,G4ThreeVector(xTlate*CLHEP::mm,yTlate*CLHEP::mm,0),
					  tName2,
					  myTargetLog,
					  myCalPhys,
					  false,
					  copyNo);

#ifdef SHOOTERGEOMETRYDEBUG
	    G4cout<<"Placement of a " << tName1 << " : (" << xTlate << "," << yTlate << ",0)" << G4endl;

	    if(placedTarget->CheckOverlaps()){G4cout<<"Overlapping in first loop!" <<G4endl; return NULL;}
#endif

	  }
      }

#ifdef SHOOTERGEOMETRYDEBUG
    if(myWorldPhys->CheckOverlaps()){
    	G4cout<<"Overlapping world!" <<G4endl; return NULL;}

    if(myCalPhys->CheckOverlaps()){
    	G4cout<<"Overlapping cal!" <<G4endl; return NULL;}
#endif

    return myWorldPhys;
}

#endif

