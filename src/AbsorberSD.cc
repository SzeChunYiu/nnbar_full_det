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

#include "AbsorberSD.hh"
#include "TubeSD.hh"

#include "NNbarHit.hh"

#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4UserEventAction.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"

#include "G4Cerenkov.hh"


//......
AbsorberSD::AbsorberSD(G4String name):
G4VSensitiveDetector(name)
{
	error_count = 0;
    G4String HCname;
    collectionName.insert(HCname="AbsorberHitCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
}

//......
AbsorberSD::~AbsorberSD()
{}

//......
void AbsorberSD::Initialize(G4HCofThisEvent*)
{
    
    HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,collectionName[0]);
}

//.....
G4bool AbsorberSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{

    if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "Abso") return false;
    
    
    // Get Direction
    G4Track * theTrack = aStep  ->  GetTrack();
    G4ThreeVector stepDelta = aStep->GetDeltaPosition();
    G4double direction = stepDelta.getZ();
    
    
    //Get particle name
    G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
    G4String particleName =  particleDef -> GetParticleName();
    
    // Get particle PDG code
    G4int pdg = particleDef ->GetPDGEncoding();
    
    // Get unique track_id (in an event)
    G4int trackID = theTrack -> GetTrackID();
    
    // Get Energy deposited
    G4double energyDeposit = aStep -> GetTotalEnergyDeposit();
   
    // Get Step Length 
    G4double DX = aStep -> GetStepLength();
    G4StepPoint* PreStep = aStep->GetPreStepPoint();


    // Get Position
    G4ThreeVector pos = PreStep->GetPosition();
    G4double z = pos.getZ();

    G4ThreeVector vertex = theTrack->GetVertexPosition();
    G4double origin = vertex.getZ();
    G4double tracklength = z - origin;

    // Read voxel indexes: i is the x index, k is the z index
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int k  = touchable->GetReplicaNumber(0);
    G4int group_ID = touchable->GetReplicaNumber(4);
    G4int module_ID = touchable->GetReplicaNumber(3);

    // Get Time 
    G4double time = theTrack->GetGlobalTime() / CLHEP::ns;

    // Get Local Time
    G4double localTime = theTrack->GetLocalTime() / CLHEP::ns;

    // Get Name
    G4String name = theTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    
    G4TouchableHandle touchPreStep = PreStep->GetTouchableHandle();
    G4VPhysicalVolume* volumePre = touchPreStep->GetVolume();
    G4String namePre = volumePre->GetName();
    
    G4int photons = 0;
    G4ParticleDefinition* particle;
    if (particleName != "opticalphoton") {
            const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();
            for (int j = 0; j < (*secondary).size(); j++) {
                    particle = (*secondary)[j]->GetDefinition();
                    if (particle->GetParticleName() == "opticalphoton" && (*secondary)[j]->GetCreatorProcess()->GetProcessName() == "Cerenkov") { photons++; } // Cerenkov exists in scintillator
            }
    }

    // Get Process
    G4int parentID = 0;
    G4String proc = "";
    // Getting Process ond parentID of primary causes seg fault
    if (trackID > 1 && theTrack->GetOriginTouchable()->GetVolume()->GetName() != "World"){
		//std::cout << name << " " <<theTrack->GetOriginTouchable()->GetVolume()->GetName() << std::endl;
		parentID = theTrack->GetParentID();
		if (parentID != 0) { proc = theTrack->GetCreatorProcess()->GetProcessName(); }
		else { proc = "primary"; }
		//std::cout << " Abso hit : " << name << " ID : " << trackID << theTrack->GetOriginTouchable()->GetVolume()->GetName() << "  " << proc << " the parent ID is : " << parentID << std::endl;
    } 

	else {
        proc = "primary";
		parentID = 0;
    }
	/***
    if (proc=="Decay") {
        G4cout << "Killing particle " << name << G4endl;
        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }
	***/

    //if (DX) { // *** tried to cancel the statement to see if things work
	    
    // Get the pre-step kinetic energy
    G4double eKinPre = aStep -> GetPreStepPoint() -> GetKineticEnergy();
    // Get the post-step kinetic energy
    G4double eKinPost = aStep -> GetPostStepPoint() -> GetKineticEnergy();
    // Get the step average kinetic energy
    G4double eKinMean = (eKinPre + eKinPost) * 0.5;

    NNbarHit* detectorHit = new NNbarHit();
    detectorHit -> SetLocalTime(localTime);
    detectorHit -> SetParentID(parentID);
    detectorHit -> SetProcess(proc);
    detectorHit -> SetTime(time);
    detectorHit -> SetName(name);
    detectorHit -> SetTrackID(trackID);
    detectorHit -> SetXID(k);
    detectorHit -> SetGroup_ID(group_ID);
    detectorHit -> SetMod_ID(module_ID);
    detectorHit -> SetPosZ(tracklength);
    detectorHit -> SetEDep(energyDeposit);
    detectorHit -> SetKinEn(eKinPost);
    detectorHit-> SetPhotons(photons);

    HitsCollection -> insert(detectorHit);

	if (parentID == 0) {
		G4double x_o = (aStep->GetPreStepPoint()->GetPosition())[0];
		G4double x_p = (aStep->GetPostStepPoint()->GetPosition())[0];
		G4double y_o = (aStep->GetPreStepPoint()->GetPosition())[1];
		G4double y_p = (aStep->GetPostStepPoint()->GetPosition())[1];
		G4double z_o = (aStep->GetPreStepPoint()->GetPosition())[2];
		G4double z_p = (aStep->GetPostStepPoint()->GetPosition())[2];

		if (x_o == x_p && y_o == y_p && z_o == z_p && DX < 10^-10 && energyDeposit < 10^-10)
                {error_count++;
                //std::cout << "- - - - Error with the propagation " << error_count << std::endl;
                //std::cout << DX << " Edep : " << energyDeposit << std::endl;
                if (error_count > 100) { std::cout << " *** - - - - Error with the propagation " << error_count << " *** - - - - -"<< std::endl;
                theTrack->SetTrackStatus(fStopAndKill);
                error_count = 0;}
		}

	}

    return true;
}

void AbsorberSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    static G4int HCID = -1;
    if(HCID < 0)
    {HCID = GetCollectionID(0);}
    HCE -> AddHitsCollection(HCID,HitsCollection);
}

