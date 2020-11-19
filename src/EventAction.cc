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


#include "EventAction.hh"
#include "Analysis.hh"
#include "ScintillatorSD.hh"
#include "AbsorberSD.hh"
#include "TubeSD.hh"
#include "TPCSD.hh"
#include "NNbarHit.hh"

#include "G4VHitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....

EventAction::EventAction(): 
    G4UserEventAction(),
    fAbsoEdepHCID(-1),
    fGapEdepHCID(-1),
    fAbsoTrackLengthHCID(-1),
    fCerenkovHCID(-1),
    fScintTrackLengthHCID(-1),
    scintHitsCollectionID(-1),
    absHitsCollectionID(-1),
    tubeHitsCollectionID(-1),
    TPCHitsCollectionID(-1)
{}

//....

EventAction::~EventAction()
{}

//....

G4THitsMap<G4double>* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<G4THitsMap<G4double>*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....

G4double EventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0.;
  for ( auto it : *hitsMap->GetMap() ) {
    // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
    sumValue += *(it.second);
  }
  return sumValue;  
}  

//.....

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    if(scintHitsCollectionID == -1) {
       scintHitsCollectionID = pSDManager->GetCollectionID("ScintillatorHitCollection");
       absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");
       tubeHitsCollectionID = pSDManager->GetCollectionID("TubeHitCollection");
       TPCHitsCollectionID = pSDManager->GetCollectionID("TPCHitCollection");
  }
}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
 
    if(scintHitsCollectionID  < 0) {return;}

	G4HCofThisEvent* HCE = event->GetHCofThisEvent();

    int CHCID = -1;
    if (CHCID<0) {CHCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitCollection");}
    
    int CHCID2 = -1;
    if (CHCID2<0) {CHCID2 = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitCollection");}

    int CHCID3 = -1;
    if (CHCID3<0) {CHCID3 = G4SDManager::GetSDMpointer()->GetCollectionID("TubeHitCollection");}

    int CHCID4 = -1;
    if (CHCID4<0) {CHCID4 = G4SDManager::GetSDMpointer()->GetCollectionID("TPCHitCollection");}

    NNbarHitsCollection* ScintHits = 0;
    NNbarHitsCollection* AbsHits   = 0;
    NNbarHitsCollection* TubeHits  = 0;
    NNbarHitsCollection* TPCHits  = 0;
	
    if (HCE) {

        G4AnalysisManager* analysis = G4AnalysisManager::Instance();
        G4int b = 1;
        G4int ltime     = 0.;
    	G4int parentID  = 0;
        G4String proc   = "";
        G4String name   = "";
        G4double time   = 0.;
        G4int trID      = 0;
        G4int i         = 0;
        G4double kinEn  = 0.;
        G4double eDep   = 0.;
        G4double trackl = 0.;	
        G4int hitCount  = 0;
        G4int org_replica = 99;
        G4int group_ID = 999;
        G4int module_ID = 999;

        ScintHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID));
        // Book vector to keep track of Edep in each Scintillator Sheet
        G4double EdepPerSheet[10] = {0., 0., 0., 0., 0.,0., 0., 0., 0., 0.};
        G4int ScintPerSheet[10] = { 0,0,0,0,0,0,0,0,0,0 };
        G4int ScintPerSheet_new[10] = { 0,0,0,0,0,0,0,0,0,0 };
        G4double totEdep   = 0.;     
        G4double eDepScint = 0.;
        G4double eDepAbs   = 0.;
        G4double eDepTube  = 0.; 
        G4double extraEdep = 0.;
        G4double eDepCompt = 0.;
        G4double eDepInelastic= 0.;
        G4double eDephIoni = 0.;
        G4double eDepHadElas = 0.;
        G4double eDepPrimary = 0.;
        G4double eDepOther = 0.;
        G4int cerenkovCounter = 0;
        G4int scint_photons = 0;
        G4int scint_photons_check = 0;

        if (ScintHits) {

           hitCount = ScintHits->entries();
           for (G4int h=0; h<hitCount; h++) {

               ltime    = ((*ScintHits)[h]) -> GetLocalTime();
               parentID = ((*ScintHits)[h]) -> GetParentID();
               proc     = ((*ScintHits)[h]) -> GetProcess();
               name     = ((*ScintHits)[h]) -> GetName();
               time     = ((*ScintHits)[h]) -> GetTime();
               trID     = ((*ScintHits)[h]) -> GetTrackID();
               i        = ((*ScintHits)[h]) -> GetXID();
               group_ID = ((*ScintHits)[h]) -> GetGroup_ID();
               module_ID = ((*ScintHits)[h]) -> GetMod_ID();
               kinEn    = ((*ScintHits)[h]) -> GetKinEn();
               eDep     = ((*ScintHits)[h]) -> GetEdep();
               trackl   = ((*ScintHits)[h]) -> GetPosZ();
               org_replica = ((*ScintHits)[h])->GetOrigin();
               G4int scint_photons_per_hit = ((*ScintHits)[h])->GetPhotons();


               std::cout<< " Scint HIT :: Module index: " << module_ID << " Group ID " << group_ID << " Layer : " << i << std::endl;

               if (proc == "Decay") {continue;}

               // ** An issue here is that the cerenkov light is also created inside the scintillator ?!

               ScintPerSheet[i] = ScintPerSheet[i] + scint_photons_per_hit;
               // Sum eDep for each scintillator sheet
               EdepPerSheet[i] = EdepPerSheet[i] + eDep;
               // Sum totEdep
               eDepScint += eDep;
               totEdep += eDep;

               if (proc != "primary" & eDep > 0) {
                   extraEdep += eDep;
                   if (proc == "compt") eDepCompt += eDep;
                   else if (proc == "pi+Inelastic") eDepInelastic += eDep;
                   else if (proc == "hIoni") eDephIoni += eDep;
                   else if (proc == "hadElastic") eDepHadElas += eDep;
                   else eDepOther += eDep;
                   //continue;
               }

               if (trID ==1) {
                   // Sum eDep for each scintillator sheet
                   //EdepPerSheet[i] += eDep;
                   eDepPrimary += eDep;

                   //G4cout << "Kinetic Energy: " << kinEn << " " << G4endl;
                   analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);

                   // When primary particle stops
                   if (kinEn == 0) {
                       analysis->FillH1(33, trackl/CLHEP::cm);
                       analysis->FillH1(32, time/CLHEP::ns);
                       // Filling only when trackID==1
                       analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV); }
                    }
                }

            // Fill scint bins with Energy Dep
            for (G4int i=0; i<10; i++) {
                analysis->FillH1(i, EdepPerSheet[i]/CLHEP::MeV);
                analysis->FillNtupleDColumn(i+10, EdepPerSheet[i] / CLHEP::MeV);}

            for (G4int i = 0; i < 10; i++) {
                analysis->FillH1(i+10, ScintPerSheet[i] );
                analysis->FillNtupleIColumn(i + 20, ScintPerSheet[i]);}
            
            analysis->FillH1(34, eDepScint/CLHEP::MeV);

        }
 
        AbsHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID2));

        if(AbsHits) {
            hitCount = AbsHits->entries();
            for (G4int h=0; h<hitCount; h++) {
                ltime           = ((*AbsHits)[h]) -> GetLocalTime();
                parentID	= ((*AbsHits)[h]) -> GetParentID();
                proc            = ((*AbsHits)[h]) -> GetProcess();
                G4String name   = ((*AbsHits)[h]) -> GetName();
                G4double time   = ((*AbsHits)[h]) -> GetTime();
                G4int trID      = ((*AbsHits)[h]) -> GetTrackID();
                G4int i         = -99;
                group_ID = ((*AbsHits)[h]) -> GetGroup_ID();
                module_ID = ((*AbsHits)[h]) -> GetMod_ID();
                G4double kinEn  = ((*AbsHits)[h]) -> GetKinEn();
                G4double eDep   = ((*AbsHits)[h]) -> GetEdep();
                G4double trackl = ((*AbsHits)[h]) -> GetPosZ();
                G4double photons_cerenkov = ((*AbsHits)[h])->GetPhotons();
                //if (name == "opticalphoton" && proc == "Cerenkov"){ std::cout << "*** Abs Hit ! " << parentID << " Name " << name << " proc " << proc << std::endl;}
                cerenkovCounter = cerenkovCounter + photons_cerenkov;
                if (proc == "Decay") {continue;}
                eDepAbs += eDep;
                totEdep += eDep;

                std::cout<< " ABS HIT :: Module index: " << module_ID << " Group ID " << group_ID << std::endl;

                if (proc != "primary" & eDep > 0) {
                    extraEdep += eDep;
                    if (proc == "compt") eDepCompt += eDep;
                    else if (proc == "pi+Inelastic") eDepInelastic += eDep;
                    else if (proc == "hIoni") eDephIoni += eDep;
                    else if (proc == "hadElastic") eDepHadElas += eDep;
                    else eDepOther += eDep;
                    //continue; // whats the reason?
                }
 	   
                if (trID == 1){
                    eDepPrimary += eDep;
                    analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
                    if (kinEn == 0) {
                        analysis->FillH1(33, trackl/CLHEP::cm);
                        analysis->FillH1(32, time/CLHEP::ns);
                        analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV);
                    }
                }

                if (name == "opticalphoton" && proc == "Cerenkov"){analysis->FillH1(31,time);}
            }
	    
            if (eDepAbs>0){analysis->FillH1(35, eDepAbs/CLHEP::MeV);}
            analysis->FillH1(20, cerenkovCounter);
            analysis->FillNtupleIColumn(30, cerenkovCounter);
        }
       
        TubeHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID3));
        if (TubeHits) {
	    hitCount = TubeHits->entries();

            for (G4int h=0; h<hitCount; h++) {
                ltime    = ((*TubeHits)[h]) -> GetLocalTime();
	        parentID = ((*TubeHits)[h]) -> GetParentID();
    		proc     = ((*TubeHits)[h]) -> GetProcess();
       	        name     = ((*TubeHits)[h]) -> GetName();
       	        time     = ((*TubeHits)[h]) -> GetTime(); 
	        trID     = ((*TubeHits)[h]) -> GetTrackID();
                group_ID = ((*TubeHits)[h]) -> GetGroup_ID();
                module_ID= ((*TubeHits)[h]) -> GetMod_ID();
		i        = ((*TubeHits)[h]) -> GetXID();
	        kinEn    = ((*TubeHits)[h]) -> GetKinEn();
	        eDep     = ((*TubeHits)[h]) -> GetEdep();
                trackl   = ((*TubeHits)[h]) -> GetPosZ();	

                std::cout<< " Tube HIT :: Module index: " << module_ID << " Group ID " << group_ID << std::endl;


                if (proc == "Decay") {continue;}

                // Sum totEdep
                eDepTube += eDep;  
                totEdep += eDep;

                if (proc != "primary" & eDep > 0) {
                    extraEdep += eDep;
                    if (proc == "compt") eDepCompt += eDep;
                    else if (proc == "pi+Inelastic") eDepInelastic += eDep;
                    else if (proc == "hIoni") eDephIoni += eDep;
                    else if (proc == "hadElastic") eDepHadElas += eDep;
                    else eDepOther += eDep;
                    //continue;
                }

                if (trID ==1) {
	            eDepPrimary += eDep;
     	            analysis->FillH2(0, trackl/CLHEP::cm, kinEn/CLHEP::MeV);
	            if (kinEn == 0) {
                        analysis->FillH1(33, trackl/CLHEP::cm);
                        analysis->FillH1(32, time/CLHEP::ns);
                        analysis->FillH2(1, trackl/CLHEP::cm, totEdep/CLHEP::MeV);
                    }
	        }
            }
            // Fill total Edep in Vacuum Tube
            analysis->FillH1(36, eDepTube/CLHEP::MeV);	
            //G4cout << "Total Edep in tube: " << eDepTube/CLHEP::MeV << G4endl;         
        }

        TPCHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID4));

        if (TPCHits) {
        hitCount = TPCHits->entries();
        for (G4int h=0; h<hitCount; h++) {

            ltime    = ((*TPCHits)[h]) -> GetLocalTime();
            parentID = ((*TPCHits)[h]) -> GetParentID();
            proc     = ((*TPCHits)[h]) -> GetProcess();
            name     = ((*TPCHits)[h]) -> GetName();
            time     = ((*TPCHits)[h]) -> GetTime();
            trID     = ((*TPCHits)[h]) -> GetTrackID();
            group_ID = ((*TPCHits)[h]) -> GetGroup_ID();
            module_ID = ((*TPCHits)[h]) -> GetMod_ID();
            kinEn    = ((*TPCHits)[h]) -> GetKinEn();
            eDep     = ((*TPCHits)[h]) -> GetEdep();
            trackl   = ((*TPCHits)[h]) -> GetPosZ();


            std::cout<< " TPC HIT :: Module index: " << module_ID << " Group ID " << group_ID << std::endl;
            /***
            eDepShield += eDep;
            totEdep += eDep;
            if (proc != "primary" & eDep > 0) {
               extraEdep += eDep;
               if (proc == "compt") eDepCompt += eDep;
               else if (proc == "pi+Inelastic") eDepInelastic += eDep;
               else if (proc == "hIoni") eDephIoni += eDep;
               else if (proc == "hadElastic") eDepHadElas += eDep;
               else eDepOther += eDep;
               //continue;
           }
            if (trID ==1) {eDepPrimary += eDep;}
            ***/
        }
        /***
        analysis->FillH1(37, eDepShield/CLHEP::MeV);
        analysis->FillNtupleDColumn(33, eDepShield/CLHEP::MeV);
        std::cout<< "Energy absorbed by shield: "<<eDepShield << " MeV *** " << eDepShield/CLHEP::MeV << std::endl;
        ***/
        }


        /// End of all hit processing, now doing the remaining things...

        std::cout << " number of scint photons " << scint_photons << " :: old method " << scint_photons_check << std::endl;
        for (int j = 0; j < 10; j++) { std::cout << "Edep " << j << ": " << EdepPerSheet[j] << " Photons: " << ScintPerSheet[j] << std::endl; }
        std::cout << "Edep Abs: " <<  eDepAbs<<" number of cerenkov " << cerenkovCounter << std::endl;

         if (totEdep > 0) {analysis->FillNtupleIColumn(31, 1);}
         else {analysis->FillNtupleIColumn(31, 0); std::cout << " === No hit for the detector since Edep = 0 " << std::endl;}

         /***
         G4cout << G4endl;
         G4cout << "---------Energy Depostied by Volume----------" << G4endl;
         G4cout << "Total Edep in tube: " << eDepTube / CLHEP::MeV << " MeV" << G4endl;
         G4cout << "Total Edep in scint: " << eDepScint / CLHEP::MeV << " MeV" << G4endl;
         G4cout << "Total Edep in abs: " << eDepAbs / CLHEP::MeV << " MeV" << G4endl;
         G4cout << "Missing Energy: " << (totEdep - eDepTube - eDepScint - eDepAbs) / CLHEP::MeV << " MeV" << G4endl << G4endl;
         G4cout << "---------Energy Deposited by Particles--------" << G4endl;
         G4cout << "Total Edep by non primary particles: " << extraEdep / CLHEP::MeV << " MeV" << G4endl;
         G4cout << "Total Edep by primary particle: " << eDepPrimary / CLHEP::MeV << " MeV" << G4endl;
         G4cout << "Total Edep: " << totEdep / CLHEP::MeV << " MeV" << G4endl;
         G4cout << "Missing Energy = " << (totEdep - extraEdep - eDepPrimary) / CLHEP::MeV << " MeV" << G4endl << G4endl;
         G4cout << "-------- Energy Deposited by Non-primary Process-----------" << G4endl;
         G4cout << "compt: " << eDepCompt / CLHEP::MeV << G4endl;
         G4cout << "pi+Inelastic: " << eDepInelastic / CLHEP::MeV << G4endl;
         G4cout << "hIoni: " << eDephIoni / CLHEP::MeV << G4endl;
         G4cout << "hadElastic: " << eDepHadElas / CLHEP::MeV << G4endl;
         G4cout << "Other: " << eDepOther / CLHEP::MeV << G4endl;
         G4cout << "Missing Energy: " << (totEdep - eDepPrimary - eDepCompt - eDepInelastic - eDephIoni - eDepHadElas - eDepOther) / CLHEP::MeV << " MeV" << G4endl << "---------------------------------" << G4endl << G4endl;
         ***/

        analysis->AddNtupleRow();
    }
    
    else {G4cout << "No HCE" << G4endl;}
    auto eventID = event->GetEventID();
    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {G4cout << "---> End of event: " << eventID << G4endl;}

}  
