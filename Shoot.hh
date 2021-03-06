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
// $Id: Shoot.hh,v 1.2 2013-08-16 zsiegel $
// $Id: Shoot.hh,v 1.1 2007-10-11 13:01:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef SHOOT_HH
#define SHOOT_HH

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Timer.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

#ifdef SHOOTERDEBUG
#include <map>
#endif

G4double Shoot(const G4int numShoot,
	   G4VPhysicalVolume *pTopNode,
	   const G4ThreeVector& pSource,
	   const G4ThreeVector& pVec)
{
  G4double physStep=kInfinity;
  G4int i;
  G4double safety;
  G4double Step;
  G4Navigator myNav;
  G4Timer timer;
  G4ThreeVector partLoc;
  G4VPhysicalVolume *located=0;

#ifdef SHOOTERDEBUG
  std::map<G4String, G4int> stepmap;
  stepmap["CalPhys"] = 0;
  stepmap["TBoxPhys1"] = 0;
  stepmap["TBoxPhys2"] = 0;
  stepmap["TPhysHalfCone2"] = 0;
  stepmap["TPhysPPed2"] = 0;
  stepmap["TTubePhys1"] = 0;
  stepmap["TTubePhys2"] = 0;
  stepmap["WPhys"] = 0;
#endif
  myNav.SetWorldVolume(pTopNode);

  timer.Start();

  for (i=numShoot;i>0;i--)
    {
      //      G4cout << "#Loop  " << i << G4endl ;
      
      partLoc=pSource;
      located=myNav.LocateGlobalPointAndSetup(partLoc,0,false,true);
#ifdef SHOOTERDEBUG
      G4cout << "Source: " << partLoc << " in " << located->GetName()
    		 << G4endl << "Path Vertices:" << G4endl;
#endif
      while (located)
	{

#ifdef SHOOTERDEBUG
	  G4cout << "Loc = " << partLoc << " Vec = " << pVec << G4endl ;
	  G4cout << "Safety = " << safety << G4endl ;

	  stepmap[located->GetName()]++;
#endif
	  Step=myNav.ComputeStep(partLoc,pVec,physStep,safety);

	  partLoc+=Step*pVec;
	  myNav.SetGeometricallyLimitedStep();
	  located=myNav.LocateGlobalPointAndSetup(partLoc);

#ifdef SHOOTERDEBUG
      G4cout << partLoc;
      if(located)
    	  G4cout << " in " << located->GetName() << " --- " << G4endl;
      else
    	  G4cout << " out of world." <<G4endl;
#endif

	};
    }
  timer.Stop();
  G4cout << "Shots = " << numShoot << " " << timer << G4endl;

#ifdef SHOOTERDEBUG
  std::map<G4String, G4int>::iterator stepmapiter;

  G4int num;
  G4cout << "Volume Illuminated\tFrequency\n-------------------------------" << G4endl;
  for(stepmapiter = stepmap.begin(), num = 0; stepmapiter != stepmap.end() ; stepmapiter++, num++)
	  G4cout << stepmapiter->first << "\t" << stepmapiter->second << G4endl;
  for(num = num ; num < 8 ; num++)
	  G4cout << G4endl;
#endif

  return timer.GetUserElapsed() + timer.GetSystemElapsed();
}

#include "G4TransportationManager.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

void MagneticShoot(const G4int numShoot,
		   G4VPhysicalVolume *pTopNode,
		   const G4ThreeVector& pSource,
		   const G4ThreeVector& pVec,
		   const G4double fieldValue, // ** already in tesla **
		   const G4double DeltaChord) // ** already in mm **
{
  /** Setting up the Magnetic field **/

  G4UniformMagField magField (0.,0.,fieldValue);
  G4Navigator *myNav = G4TransportationManager::
                       GetTransportationManager()-> GetNavigatorForTracking();
  myNav->SetWorldVolume(pTopNode);

  G4double momentum = 0.05*CLHEP::proton_mass_c2;
  G4double kineticEnergy = momentum*momentum /
       (std::sqrt(momentum*momentum+CLHEP::proton_mass_c2*CLHEP::proton_mass_c2)+CLHEP::proton_mass_c2);
  G4double velocity = momentum / (CLHEP::proton_mass_c2+kineticEnergy);
  G4double labTof= 10.*CLHEP::ns;
  G4double properTof= 0.1*CLHEP::ns;
  
  /* Field Properties */

  G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&magField);
  
  /* Choose your stepper here */
  /* G4ClassicalRK4 is the default one */
  G4MagIntegratorStepper*  pStepper = new G4ClassicalRK4( fEquation );

  /*
    pStepper = new G4ExplicitEuler( fEquation );
    pStepper = new G4ImplicitEuler( fEquation );
    pStepper = new G4SimpleRunge( fEquation );
    pStepper = new G4SimpleHeum( fEquation );
    pStepper = new G4ClassicalRK4( fEquation );
    pStepper = new G4HelixExplicitEuler( fEquation );
    pStepper = new G4HelixImplicitEuler( fEquation );
    pStepper = new G4HelixSimpleRunge( fEquation );
    pStepper = new G4RKG3_Stepper( fEquation );
  */

  G4FieldManager* pFieldMgr = G4TransportationManager::
                  GetTransportationManager()->GetFieldManager();
  pFieldMgr->SetDetectorField( &magField );

  G4ChordFinder* pChordFinder = new G4ChordFinder( &magField,DeltaChord,pStepper);
  pFieldMgr->SetChordFinder( pChordFinder );

  G4PropagatorInField *pMagFieldPropagator= G4TransportationManager::
    GetTransportationManager()-> GetPropagatorInField ();

  pChordFinder->SetChargeMomentumMass(1.,               // charge in e+ units
				      momentum,         // Momentum in Mev/c ?
				      CLHEP::proton_mass_c2 );
  G4Timer timer;
  timer.Start();

  for (G4int i=numShoot;i>0;i--)
    {
      G4VPhysicalVolume *located;
      G4ThreeVector Vec = pVec ;
      G4ThreeVector partLoc = pSource ;	    
      /*
	G4cout << "#Loop  " << i << G4endl ;
	G4cout << "Loc = " << partLoc << " Vec = " << Vec << G4endl << G4endl ;
      */
      momentum = (0.5+i*0.1) * CLHEP::proton_mass_c2;
      kineticEnergy = momentum*momentum /
       (std::sqrt(momentum*momentum+CLHEP::proton_mass_c2*CLHEP::proton_mass_c2)+CLHEP::proton_mass_c2);
      velocity = momentum / (CLHEP::proton_mass_c2+kineticEnergy);

      pFieldMgr->GetChordFinder()
	       ->SetChargeMomentumMass(1,               // charge in e+ units
				       momentum,        // Momentum in Mev/c 
				       CLHEP::proton_mass_c2); // Mass
      
      located=myNav->LocateGlobalPointAndSetup(partLoc);
      while (located)
	{
	  G4double physStep= kInfinity; // 2.5*mm ;
	  G4double safety = 1.0*CLHEP::m;
	  G4double Step = 0.0*CLHEP::m;
	  /*
	  G4cout << "Loc = " << partLoc << " Vec = " << Vec << G4endl ;
	  G4cout << "Safety = " << safety << G4endl ;
	  */
	  
	  G4FieldTrack initTrack(partLoc,pVec,0.,kineticEnergy,
	                         CLHEP::proton_mass_c2,velocity,labTof,properTof,0);
	  Step=pMagFieldPropagator->ComputeStep(initTrack,physStep,safety);

	  myNav->SetGeometricallyLimitedStep();
	  
	  partLoc = pMagFieldPropagator->EndPosition();
	  Vec = pMagFieldPropagator->EndMomentumDir(); 
	  
	  located=myNav->LocateGlobalPointAndSetup(partLoc);
	};
    }
  timer.Stop();
  //  G4cout << "Shots = " << numShoot << " " << timer << G4endl;
}


void ShootVerbose(G4VPhysicalVolume *pTopNode,
		  const G4ThreeVector& pSource,
		  const G4ThreeVector& pVec)
{
  const G4double physStep=kInfinity;
  G4double safety,Step;
  G4Navigator myNav;
  G4ThreeVector partLoc;
  G4VPhysicalVolume *located=0;

  myNav.SetWorldVolume(pTopNode);

  partLoc=pSource;
  located=myNav.LocateGlobalPointAndSetup(partLoc);
  while (located)
    {
      Step=myNav.ComputeStep(partLoc,pVec,physStep,safety);
      G4cout << "Physical Location=" << located->GetName()
	     << " #" << located->GetCopyNo() << G4endl
	     << "   Step=" << Step << "  Safety=" << safety
	     << "  ---->" << G4endl;

      partLoc+=Step*pVec;
      myNav.SetGeometricallyLimitedStep();
      located=myNav.LocateGlobalPointAndSetup(partLoc);
    };
}

#endif

