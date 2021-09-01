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

int particle_name_file_index;
extern std::ofstream Particle_outFile;
extern std::vector<std::vector<G4double>> particle_gun_record;
extern G4double event_number;

PrimaryGeneratorAction::PrimaryGeneratorAction()
	:fParticleGun(nullptr)
{fParticleGun = new G4ParticleGun();}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{
	G4double x; G4double y; G4double z;
	G4double t; G4double px; G4double py; G4double pz;
	G4double KE;

	G4String particleName;

	std::vector<G4double> particle_gun_record_row;

	x = 0.0* m;
	y = 0.1*m;
	z = 0.0* m; // (*vect)[j]->z() * m
	KE = 0 *MeV; //250.0 * MeV;
	px = 0.;
	py = 1.;
	pz = 0.;
	t = 0.; 

	G4String all_type[7] = { "neutron","proton","gamma","electron","muon","pion","kaon" };
	G4double name_ID = 99;

	for (int i = 0; i <= 6; i++) { if (particleName == all_type[i]) { name_ID = i; break; } }

	fParticleGun->SetParticleDefinition(particleTable->FindParticle("pi0"));
	//fParticleGun->SetParticleEnergy(KE);
	
	int n = event_number/100; // every energy, it has 50 entries
	fParticleGun->SetParticleEnergy(0.00000001*MeV);	
	
	fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz)); //(*vect)[j]->w())
	fParticleGun->SetParticleTime(0.0);
	Particle_outFile <<  event_number << ",";
	Particle_outFile <<  0 << ","; //PID
	Particle_outFile <<  particleTable -> FindParticle("pi0") -> GetPDGMass() << ","; //Mass
	Particle_outFile <<  200 << ","; // 100+n*10 
	Particle_outFile <<  x << ",";
	Particle_outFile <<  y << ",";
	Particle_outFile <<  z << ",";
	Particle_outFile <<  t << ",";
	Particle_outFile <<  px << ",";
	Particle_outFile <<  py << ",";
	Particle_outFile <<  pz << G4endl;

	fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....

