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

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "Analysis.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

//....

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
}

//....

SteppingAction::~SteppingAction()
{ 
}

//....

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 
  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) {
     fEventNumber = eventNumber;
     fScintillationCounter = 0;
     fCerenkovCounter = 0;
  }

  G4Track* track = step->GetTrack();
  G4int ID = track->GetTrackID();
  G4int ltime = track->GetLocalTime();
} 
  //const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();

//  if (secondaries->size()>0) {
  //   for(unsigned int i=0; i<secondaries->size(); ++i) {
    //      G4cout << "secondary particle: " << secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() << G4endl;
      //    G4cout << "Creator process: " << secondaries->at(i)->GetCreatorProcess()->GetProcessName() << G4endl;
     //}
  //}

//}

//....
