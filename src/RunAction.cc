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

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include <iostream>
#include <fstream>
#include <string>
//....

int particle_name_file_index;
string particle_name_list[7] = { "", "_neutron", "_proton", "_gamma","_electron","_muon","_pion" };

RunAction::RunAction()
 : G4UserRunAction()
{ 
    fMessenger = new G4GenericMessenger(this, "/particle_generator/", "Name the particle for the file name");
    G4GenericMessenger::Command& filenameCMD = fMessenger->DeclareProperty("Particle_index", particle_name_file_index, "Index of the particle in the order");
    filenameCMD.SetParameterName("Particle Index", true);
    filenameCMD.SetDefaultValue("99");

  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager, ROOT specified in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  // Book histograms

  G4cout << "Booking histograms " << G4endl;
// Scintillator Histograms Edeposited in individual sheets for low energy gamma backgrounds keV
/***
  analysisManager->CreateH1("Scint_Bin_1", "Energy Deposit in Scint Bin", 100, 0, 1000*keV); // 0
  analysisManager->CreateH1("Scint_Bin_2", "", 100, 0, 1000*keV); // 1
  analysisManager->CreateH1("Scint_Bin_3", "", 100, 0, 1000*keV); // 2
  analysisManager->CreateH1("Scint_Bin_4", "", 100, 0, 1000*keV); // 3
  analysisManager->CreateH1("Scint_Bin_5", "", 100, 0, 1000*keV); // 4
  analysisManager->CreateH1("Scint_Bin_6", "", 100, 0, 1000*keV); // 5
  analysisManager->CreateH1("Scint_Bin_7", "", 100, 0, 1000*keV); // 6
  analysisManager->CreateH1("Scint_Bin_8", "", 100, 0, 1000*keV); // 7
  analysisManager->CreateH1("Scint_Bin_9", "", 100, 0, 1000*keV); // 8
  analysisManager->CreateH1("Scint_Bin_10", "", 100, 0, 1000*keV); //9
***/


// Scintillator Histograms Edeposited in individual sheets for Primary Particles
  analysisManager->CreateH1("Scint_Bin_1", "Energy Deposit in Scint Bin", 100, 0, 100*MeV); // 0
  analysisManager->CreateH1("Scint_Bin_2", "", 100, 0, 100*MeV); // 1
  analysisManager->CreateH1("Scint_Bin_3", "", 100, 0, 100*MeV); // 2
  analysisManager->CreateH1("Scint_Bin_4", "", 100, 0, 100*MeV); // 3
  analysisManager->CreateH1("Scint_Bin_5", "", 100, 0, 100*MeV); // 4
  analysisManager->CreateH1("Scint_Bin_6", "", 100, 0, 100*MeV); // 5
  analysisManager->CreateH1("Scint_Bin_7", "", 100, 0, 100*MeV); // 6
  analysisManager->CreateH1("Scint_Bin_8", "", 100, 0, 100*MeV); // 7
  analysisManager->CreateH1("Scint_Bin_9", "", 100, 0, 100*MeV); // 8
  analysisManager->CreateH1("Scint_Bin_10", "", 100, 0, 100*MeV); //9

  analysisManager->CreateH1("Scint_photon_1", "", 100, 0, 10000); // 10
  analysisManager->CreateH1("Scint_photon_2", "", 100, 0, 10000); // 11
  analysisManager->CreateH1("Scint_photon_3", "", 100, 0, 10000); // 12
  analysisManager->CreateH1("Scint_photon_4", "", 100, 0, 10000); // 13
  analysisManager->CreateH1("Scint_photon_5", "", 100, 0, 10000); // 14
  analysisManager->CreateH1("Scint_photon_6", "", 100, 0, 10000); // 15
  analysisManager->CreateH1("Scint_photon_7", "", 100, 0, 10000); // 16
  analysisManager->CreateH1("Scint_photon_8", "", 100, 0, 10000); // 17
  analysisManager->CreateH1("Scint_photon_9", "", 100, 0, 10000); // 18
  analysisManager->CreateH1("Scint_photon_10", "", 100, 0, 10000); //19

  analysisManager->CreateH1("NumCerenkov", "Num Cerenkov Photons", 80, 0, 30000);//20
  analysisManager->CreateH1("Particle_ID", "", 6, 0, 5); //21
  analysisManager->CreateH1("x", "", 100, -5.0, 5.0); // 22
  analysisManager->CreateH1("y", "", 100, -5.0, 5.0 ); // 23
  analysisManager->CreateH1("z", "", 100, -5.0, 5.0 ); // 24
  analysisManager->CreateH1("t", "", 100, 0, 1000.0); // 25
  analysisManager->CreateH1("KE", "", 1000, 0, 10000.0); // 26
  analysisManager->CreateH1("px", "", 50, -1.0, 1.0); // 27
  analysisManager->CreateH1("py", "", 50, -1.0, 1.0); // 28
  analysisManager->CreateH1("pz", "", 50, -1.0, 1.0); // 29
  analysisManager->CreateH1("charge", "", 3,-1, 1); // 30  
  
  // 1-D Histos

  analysisManager->CreateH1("PhotonTime", "Cerenkov Photon Production", 50, 0, 10*ns); 		// 	11
  analysisManager->CreateH1("DecayTime", "Primary Decay Time", 50, 0, 3*ns); 			//     	12 
  analysisManager->CreateH1("Range", "Primary Particle Range", 55, 0, 55); 			//	13
  analysisManager->CreateH1("EdepScint", "Energy Deposited in Scintillators", 50, 0, 250*MeV);  //      14
  analysisManager->CreateH1("EdepAbs", "Energy Deposited in Lead-glass", 50, 0, 250*MeV);       //      15
  analysisManager->CreateH1("EdepTube", "Energy Deposited in Vacuum Tube", 50, 0, 250*MeV);     //      16

  // 2-D Histos
  // name, title, nxbins, xmin, xmax, nybins, ymin, ymax
  analysisManager->CreateH2("KinE","Kinetic Energy", 550, 0, 55, 350, 0, 350); 			//	0 
  analysisManager->CreateH2("eDepvRange", "Energy Deposited", 55, 0, 55, 50, 0, 250*MeV );      //      1
  analysisManager->CreateH2("eDepvCerenkov", "Energy Deposited v Cerenkov", 80, 0, 30000, 50, 0, 250*MeV);     //      2

  G4int a = 0; 
  G4int b = 1;

  analysisManager->CreateNtuple("ntuple_0","info");
  analysisManager->CreateNtupleDColumn("Particle_ID"); //0
  analysisManager->CreateNtupleDColumn("x"); // 1
  analysisManager->CreateNtupleDColumn("y"); // 2
  analysisManager->CreateNtupleDColumn("z"); // 3
  analysisManager->CreateNtupleDColumn("t"); // 4
  analysisManager->CreateNtupleDColumn("KE"); // 5
  analysisManager->CreateNtupleDColumn("px"); // 6
  analysisManager->CreateNtupleDColumn("py"); // 7
  analysisManager->CreateNtupleDColumn("pz"); // 8
  analysisManager->CreateNtupleDColumn("charge"); // 9
  analysisManager->CreateNtupleDColumn("Scint_Edep_1"); // 10
  analysisManager->CreateNtupleDColumn("Scint_Edep_2"); // 11
  analysisManager->CreateNtupleDColumn("Scint_Edep_3"); // 12
  analysisManager->CreateNtupleDColumn("Scint_Edep_4"); // 13
  analysisManager->CreateNtupleDColumn("Scint_Edep_5"); // 14
  analysisManager->CreateNtupleDColumn("Scint_Edep_6"); // 15
  analysisManager->CreateNtupleDColumn("Scint_Edep_7"); // 16
  analysisManager->CreateNtupleDColumn("Scint_Edep_8"); // 17
  analysisManager->CreateNtupleDColumn("Scint_Edep_9"); // 18
  analysisManager->CreateNtupleDColumn("Scint_Edep_10"); //19
  analysisManager->CreateNtupleIColumn("Scint_photon_1"); // 20
  analysisManager->CreateNtupleIColumn("Scint_photon_2"); // 21
  analysisManager->CreateNtupleIColumn("Scint_photon_3"); // 22
  analysisManager->CreateNtupleIColumn("Scint_photon_4"); // 23
  analysisManager->CreateNtupleIColumn("Scint_photon_5"); // 24
  analysisManager->CreateNtupleIColumn("Scint_photon_6"); // 25
  analysisManager->CreateNtupleIColumn("Scint_photon_7"); // 26
  analysisManager->CreateNtupleIColumn("Scint_photon_8"); // 27
  analysisManager->CreateNtupleIColumn( "Scint_photon_9"); // 28
  analysisManager->CreateNtupleIColumn( "Scint_photon_10"); //29
  analysisManager->CreateNtupleIColumn( "NumCerenkov"); //30
  analysisManager->CreateNtupleIColumn( "Hit"); //31
  analysisManager->CreateNtupleDColumn( "Edep_abs"); //32
  analysisManager->FinishNtuple();
}

//....

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//.....

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  std::cout << " **** The particle file index is " << particle_name_file_index << " " << particle_name_list[particle_name_file_index] << std::endl;

  G4String fileName = "./output/calo-sim" + particle_name_list[particle_name_file_index] +".root";
  //G4String Dir = "../output/";
  //G4String fileName = Dir + "scint-calo-sim.root";
  //G4String fileName = Dir + "abs-calo-sim.root";
  //G4cout << "Printing to " << fileName << G4endl;

  analysisManager->OpenFile(fileName);
}

//....

void RunAction::EndOfRunAction(const G4Run*)
{
  
  auto analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....
