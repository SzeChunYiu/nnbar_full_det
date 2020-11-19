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

#include "PrimaryGeneratorAction.hh"
#include <iomanip>
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"
#include "Analysis.hh"
//.....

using namespace std;

//int particle_name_file_index;

PrimaryGeneratorAction::PrimaryGeneratorAction()
    :fParticleGun(nullptr)

{fParticleGun = new G4ParticleGun();}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{delete fParticleGun;}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{
    G4double x; G4double y; G4double z;
    G4double t; G4double px; G4double py; G4double pz;
    G4double KE;

    G4String particleName;

    x = 0.0 * m; y = 0.0 * m; z = 0.0 * m; t = 14.52;
    KE = 250.95 * MeV; px = 0.0; py = 0.8; pz = -0.2;

    G4String all_type[7] = { "neutron","proton","gamma","electron","muon","pion","kaon" };
    G4double name_ID = 99;

    for (int i = 0; i <= 6; i++) { if (particleName == all_type[i]) { name_ID = i; break; } }

    fParticleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
    fParticleGun->SetParticleEnergy(KE);
    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz)); //(*vect)[j]->w())
    fParticleGun->SetParticleTime(t);

    G4AnalysisManager* analysis = G4AnalysisManager::Instance();
    analysis->FillH1(21, name_ID); // particle ID
    analysis->FillH1(22, x / m);
    analysis->FillH1(23, y / m);
    analysis->FillH1(24, z / m);
    analysis->FillH1(25, t);
    analysis->FillH1(26, KE / MeV);
    analysis->FillH1(27, px);
    analysis->FillH1(28, py);
    analysis->FillH1(29, pz);
    analysis->FillH1(30, -1);

    G4int a = 0;
    analysis->FillNtupleDColumn(0, name_ID);
    analysis->FillNtupleDColumn( 1, x / m);
    analysis->FillNtupleDColumn( 2, y / m);
    analysis->FillNtupleDColumn( 3, z / m);
    analysis->FillNtupleDColumn( 4, t);
    analysis->FillNtupleDColumn( 5, KE / MeV);
    analysis->FillNtupleDColumn( 6, px);
    analysis->FillNtupleDColumn( 7, py);
    analysis->FillNtupleDColumn( 8, pz);
    analysis->FillNtupleDColumn( 9, -1);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    std::cout<<" === Particle Fired === " << std::endl;
}


//....

