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

#include "DetectorConstruction.hh"
#include "ScintillatorSD.hh"
#include "AbsorberSD.hh"
#include "TubeSD.hh"
#include "TPCSD.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RotationMatrix.hh"

#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4PSPopulation.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include<string>
//#include "G4GDMLParser.hh"
#include <boost/lexical_cast.hpp>

using boost::lexical_cast;

// read the position of the lead glass blocks stored in csv file
std::string filename_dataxy = "/home/billy/nnbar/simulation/nnbar-full-det-test-source/lead_glass_position/lead_glass_posx.csv";
std::string filename_dataz = "/home/billy/nnbar/simulation/nnbar-full-det-test-source/lead_glass_position/lead_glass_posz.csv";
std::string filename_data_non_ax = "/home/billy/nnbar/simulation/nnbar-full-det-test-source/lead_glass_position/non_axial.csv";
std::string filename_data_non_ax_fb = "/home/billy/nnbar/simulation/nnbar-full-det-test-source/lead_glass_position/non_axial_fb.csv";
std::vector<std::vector<double>> data_x;
std::vector<std::vector<double>> data_z;
std::vector<std::vector<double>> data_non_axial;
std::vector<std::vector<double>> data_non_axial_fb;
std::vector<G4RotationMatrix *> rot_array;
std::vector<G4RotationMatrix *> rot_array2;
std::vector<G4RotationMatrix *> rot_arrayZ;
std::vector<G4RotationMatrix *> rot_arrayZ2;


void import_lead_glass_pos(std::string file_name, std::vector<std::vector<double> >& data) {
	
	std::string row;
	std::ifstream init_file(file_name.c_str());

	// open file
	if (init_file.is_open()) {
		std::cerr << "Opening Position file : "<< file_name << " ... " << std::endl;
		// loop in each line in the file
		int count_line = 0;
		while (getline(init_file, row)) {
			count_line++;
			std::istringstream iss(row);
			// initialize a vector to store the row elements
			std::vector<double> row;
			std::string token;
			// get each element by splitting the row string by commas
			while (std::getline(iss, token, ',')) {
				// convert the string to int (or use stof for floats)
				row.push_back(boost::lexical_cast<double>(token.c_str()));
			}
			data.push_back(row);
			//std::cerr << "Reading the " << count_line << " th line in file " << file_name <<" ... "<< std::endl;
		}
		init_file.close();
		std::cerr << "Lead_glass position loaded " << std::endl;
	}
	else
		std::cerr << "ERROR: Unable to open file" << std::endl;
	return;
}


DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(), fCheckOverlaps(true)
{
  import_lead_glass_pos(filename_dataxy, data_x);
  import_lead_glass_pos(filename_dataz, data_z);
  import_lead_glass_pos(filename_data_non_ax, data_non_axial);
  import_lead_glass_pos(filename_data_non_ax_fb, data_non_axial_fb);
}
//....
DetectorConstruction::~DetectorConstruction()
{}
//....
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  // Define volumes
  return DefineVolumes();
}
//....
void DetectorConstruction::DefineMaterials()
{ 

    std::cout << " define mat" << std::endl;
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  G4double a; G4double z; G4double density;

  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elTi = nistManager->FindOrBuildElement("Ti");
  G4Element* elAs = nistManager->FindOrBuildElement("As");
  G4Element* elPb = nistManager->FindOrBuildElement("Pb");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elSi = nistManager->FindOrBuildElement("Si");
  G4Element* elN = nistManager->FindOrBuildElement("N");
  G4Element* elAr = nistManager->FindOrBuildElement("Ar");

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Tube
  G4Material* Al = new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);

  // Silicon
  G4Material* Silicon = new G4Material("Silicon", z=14., a= 28.0855*g/mole, density = 2.33*g/cm3);

  // --------Air
  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, 2);
  Air->AddElement(elN, 0.7);
  Air->AddElement(elO, 0.3);

  // ---------FR4
  //----- Epoxy
  G4Material* Epoxy = new G4Material("Epoxy" , density=1.2*g/cm3, 2);
  Epoxy->AddElement(elH, 2);
  Epoxy->AddElement(elC, 2);
  //----- SiO2 (Quarz)
  G4Material* SiO2 = new G4Material("SiO2",density= 2.200*g/cm3, 2);
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO , 2);
  //FR4 (Glass + Epoxy)
  G4Material* FR4 = new G4Material("FR4" , density=1.86*g/cm3, 2);
  FR4->AddMaterial(Epoxy, 0.472);
  FR4->AddMaterial(SiO2, 0.528);

  // ----------TPC
  // CO2
  G4Material* CO2 = new G4Material("CO2", density= 1.98*g/cm3, 2);
  CO2->AddElement(elO, 2);
  CO2->AddElement(elC, 1);
  // Ar/CO2 80/20
  G4Material* Gas = new G4Material("Gas", density=1.3954*g/cm3, 2);
  Gas->AddElement(elAr, .8);
  Gas->AddMaterial(CO2, .2);

  // BC-408 taken from datasheet
  G4Material* Scint = new G4Material("Scint", 1.023*g/cm3, 2);
  Scint->AddElement(elH, 0.524573);
  Scint->AddElement(elC, 1 - 0.524573);

  // Lead-glass (taken from PDG)
  G4Material* Abs = new G4Material("Abs", 3.86*g/cm3, 5);
  Abs->AddElement(elO, 0.156453);
  Abs->AddElement(elSi, 0.080866);
  Abs->AddElement(elTi, 0.008092);
  Abs->AddElement(elAs, .002651);
  Abs->AddElement(elPb, 0.751938);

// ----------------- Generate and Add Material Properties Table ----------------
 G4double PhotonWavelength[] =
      { 2325.4*nm, 1970.1*nm, 1529.6*nm, 1060.0*nm,
        1014.0*nm, 852.10*nm, 706.50*nm, 656.30*nm,
        643.80*nm, 632.80*nm, 589.30*nm, 587.60*nm,
        546.10*nm, 486.10*nm, 480.00*nm, 435.80*nm,
        404.70*nm, 365.00*nm};

 const G4int nEntries = sizeof(PhotonWavelength)/sizeof(G4double);
 G4double PhotonEnergy[nEntries];
 for (int i=0; i < nEntries; ++i) {PhotonEnergy[i] = (1240.*nm/PhotonWavelength[i])*eV;};

 //// Lead Glass Schott SF5
 G4double refractiveIndex[] =
        { 1.63289, 1.63785, 1.64359, 1.65104,
          1.65206, 1.65664, 1.66327, 1.66661,
          1.66756, 1.66846, 1.67252, 1.67270,
          1.67764, 1.68750, 1.68876, 1.69986,
          1.71069, 1.73056};

G4MaterialPropertiesTable* absMPT = new G4MaterialPropertiesTable();
absMPT->AddProperty("RINDEX", PhotonEnergy, refractiveIndex, nEntries)->SetSpline(true);
Abs->SetMaterialPropertiesTable(absMPT);
//G4cout << "Absorber Properties -------" << G4endl;
//G4cout << *(G4Material::GetMaterialTable()) << G4endl;
absMPT->DumpTable();

//Scintillator Optical Properties

 
  const G4int nEntries2 = 12;

  G4double ScintPhotonEnergy[nEntries2] =

  { 2.08*eV, 2.38*eV, 2.58*eV, 2.7*eV, 2.76*eV,
    2.82*eV, 2.92*eV, 2.95*eV, 3.02*eV, 3.1*eV,
    3.26*eV, 3.44*eV};

  G4double rindex_scint[nEntries2] =
    {1.58, 1.58, 1.58, 1.58, 1.58,
     1.58, 1.58, 1.58, 1.58, 1.58,
     1.58, 1.58};

  G4double atten_scint[nEntries2] =
    {210*cm, 210*cm, 210*cm, 210*cm, 210*cm,
     210*cm, 210*cm, 210*cm, 210*cm, 210*cm,
     210*cm, 210*cm};

  G4double scintilFast[nEntries2] = {.04, .07, .20, .49, .84, 1.0, .83, .55, .40, .17,.03, 0.};
  G4double scintilSlow[nEntries2] = {.04, .07, .20, .49, .84, 1.0, .83, .55, .40, .17,.03, 0.};

  G4MaterialPropertiesTable *scintMPT = new G4MaterialPropertiesTable();
  scintMPT->AddProperty("RINDEX", ScintPhotonEnergy, rindex_scint, nEntries2)->SetSpline(true);
  scintMPT->AddProperty("ABSLENGTH", ScintPhotonEnergy, atten_scint, nEntries2)->SetSpline(true);
  scintMPT->AddProperty("FASTCOMPONENT", ScintPhotonEnergy, scintilFast, nEntries2)->SetSpline(true);
  scintMPT->AddProperty("SLOWCOMPONENT", ScintPhotonEnergy, scintilSlow, nEntries2)->SetSpline(true);
  // 64% of Antracene: 17400
  scintMPT->AddConstProperty("SCINTILLATIONYIELD", 1000./ MeV); //original 11136000.
  scintMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  scintMPT->AddConstProperty("FASTTIMECONSTANT", 1.0*ns); // org: 0.9
  scintMPT->AddConstProperty("SLOWTIMECONSTANT", 1.0*ns); // org: 2.1
  scintMPT->AddConstProperty("YIELDRATIO", 1.);
  Scint->SetMaterialPropertiesTable(scintMPT);
  scintMPT->DumpTable();
}

//....
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("Abs");
  auto scintMaterial = G4Material::GetMaterial("Scint");
  auto tubeMaterial = G4Material::GetMaterial("Aluminum");
  auto FR4Material = G4Material::GetMaterial("FR4");
  auto TPCMaterial = G4Material::GetMaterial("Gas");
  auto SiliconMaterial = G4Material::GetMaterial("Silicon");

  if ( ! defaultMaterial || ! absorberMaterial || ! scintMaterial) {
    G4ExceptionDescription msg; msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()","MyCode0001", FatalException, msg);}

  // World
  auto worldSizeXY = 10.0 * m; auto worldSizeZ = 10.0 * m; //1 * calorThickness;
  auto worldS = new G4Box("WorldS",worldSizeXY/2., worldSizeXY/2., worldSizeZ/2.);
  auto worldLV = new G4LogicalVolume(worldS,defaultMaterial,"World");
  auto worldPV = new G4PVPlacement(0,G4ThreeVector(),worldLV,"WorldPV",0,false,0,fCheckOverlaps);

  // Aluminum Tube
  G4double tube_radius = 1.0*m; G4double tube_len = 3.0*m; G4double tube_thickness = 2.0*cm; G4double tube_angle = 360. * deg;
  auto tubeS = new G4Cons("TubeS", tube_radius,tube_radius+tube_thickness,tube_radius,tube_radius+tube_thickness,tube_len,0.,tube_angle);
  auto tubeLV = new G4LogicalVolume(tubeS,tubeMaterial,"TubeLV");
  auto tubePV = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),tubeLV,"TubePV",worldLV,false,0,fCheckOverlaps);

  // TPC
  G4double FR4Thickness_1 = 0.16*cm; G4double FR4Thickness_2 = 0.02*cm;
  G4double gasThickness_1 = 2.0 *cm; G4double gasThickness_2 = 50.*cm;
  G4double airThickness = 1.34 *cm;
  G4double TPCThickness = FR4Thickness_1*2.+gasThickness_2+airThickness;

  // --- TPC container
  G4double TPC_h_1 = 1.87*m; G4double TPC_h_2 = 2.04*m; G4double TPC_t = 0.85*m; G4double TPC_len = 2.0*m;
  G4double TPC_total_length = 2.0 * TPC_t + TPC_h_2;
  auto TPCS_1 = new G4Box("TPCS_1",TPC_t/2.,TPC_h_1/2.,TPC_len/2.);
  auto TPCS_2 = new G4Box("TPCS_2",TPC_h_2/2.,TPC_t/2.,TPC_len/2.);
  auto TPCLV_1 = new G4LogicalVolume(TPCS_1,defaultMaterial,"TPCLV_1");
  auto TPCLV_2 = new G4LogicalVolume(TPCS_2,defaultMaterial,"TPCLV_2");

  //------ front TPCs
  auto TPC_pos1 = G4ThreeVector(-(tube_radius+tube_thickness+TPC_t/2.),TPC_h_1/2.,-TPC_len/2.);
  auto TPC_pos2 = G4ThreeVector(0.,tube_radius+tube_thickness+TPC_t/2.,-TPC_len/2.);
  auto TPC_pos3 = G4ThreeVector((tube_radius+tube_thickness+TPC_t/2.),TPC_h_1/2.,-TPC_len/2.);
  auto TPC_pos4 = G4ThreeVector(-(tube_radius+tube_thickness+TPC_t/2.),-TPC_h_1/2.,-TPC_len/2.);
  auto TPC_pos5 = G4ThreeVector(0.,-(tube_radius+tube_thickness+TPC_t/2.),-TPC_len/2.);
  auto TPC_pos6 = G4ThreeVector((tube_radius+tube_thickness+TPC_t/2.),-TPC_h_1/2.,-TPC_len/2.);

  //------ back TPCs
  auto TPC_pos7 = G4ThreeVector(-(tube_radius+tube_thickness+TPC_t/2.),TPC_h_1/2.,TPC_len/2.);
  auto TPC_pos8 = G4ThreeVector(0.,tube_radius+tube_thickness+TPC_t/2.,TPC_len/2.);
  auto TPC_pos9 = G4ThreeVector((tube_radius+tube_thickness+TPC_t/2.),TPC_h_1/2.,TPC_len/2.);
  auto TPC_pos10 = G4ThreeVector(-(tube_radius+tube_thickness+TPC_t/2.),-TPC_h_1/2.,TPC_len/2.);
  auto TPC_pos11 = G4ThreeVector(0.,-(tube_radius+tube_thickness+TPC_t/2.),TPC_len/2.);
  auto TPC_pos12 = G4ThreeVector((tube_radius+tube_thickness+TPC_t/2.),-TPC_h_1/2.,TPC_len/2.);

  new G4PVPlacement(0,TPC_pos1,TPCLV_1,"TPCPV1",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos2,TPCLV_2,"TPCPV2",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos3,TPCLV_1,"TPCPV3",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos4,TPCLV_1,"TPCPV4",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos5,TPCLV_2,"TPCPV5",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos6,TPCLV_1,"TPCPV6",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos7,TPCLV_1,"TPCPV7",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos8,TPCLV_2,"TPCPV8",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos9,TPCLV_1,"TPCPV9",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos10,TPCLV_1,"TPCPV10",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos11,TPCLV_2,"TPCPV11",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos12,TPCLV_1,"TPCPV12",worldLV,false,0,fCheckOverlaps);

  // Scintillator
  G4double scint_layer_t = 3.*cm; G4double scint_layers = 10.; int nofLayers=10;
  G4double scint_t = scint_layers * scint_layer_t; G4double scint_length = 1.6*m;
  double n_scint_module = 3.0; G4double dx = 5.0*cm; G4double dy=5.0*cm; G4double scint_w = (TPC_total_length + dy + scint_t - 2.*dx)/n_scint_module;
  G4double dz = 0.05*m;
  
  auto scint_groupS = new G4Box("Scint_groupS",(3*scint_w+2*dx)/2.,scint_t/2.,(3*scint_length+2*dz)/2.);
  auto scint_groupLV = new G4LogicalVolume(scint_groupS,defaultMaterial,"Scint_groupLV");

  auto scintS = new G4Box("ScintS",scint_w/2., scint_t/2., scint_length/2.);
  auto scintLV = new G4LogicalVolume(scintS,defaultMaterial,"ScintLV");

  auto scint_layerS = new G4Box("ScintS",scint_w/2., scint_layer_t/2., scint_length/2.);
  auto scint_layerLV = new G4LogicalVolume(scint_layerS,scintMaterial,"Scint_layerLV");
  auto scint_layerPV = new G4PVReplica("Scint_layerPV",scint_layerLV,scintLV,kYAxis,nofLayers,scint_layer_t);

  auto scint_pos1 = G4ThreeVector(-(3*scint_w+2*dx)/2. + scint_w/2., 0. , -(3*scint_length+2*dz)/2.+scint_length/2.);
  auto scint_pos2 = scint_pos1 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos3 = scint_pos2 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos4 = scint_pos1 + G4ThreeVector(0.,0., dz + scint_length);
  auto scint_pos5 = scint_pos4 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos6 = scint_pos5 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos7 = scint_pos4 + G4ThreeVector(0.,0., dz + scint_length);
  auto scint_pos8 = scint_pos7 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos9 = scint_pos8 + G4ThreeVector(dx + scint_w,0.,0.);
  
  new G4PVPlacement(0,scint_pos1,scintLV,"ScintPV1",scint_groupLV,false,1,fCheckOverlaps);
  new G4PVPlacement(0,scint_pos2,scintLV,"ScintPV2",scint_groupLV,false,2,fCheckOverlaps);
  new G4PVPlacement(0,scint_pos3,scintLV,"ScintPV3",scint_groupLV,false,3,fCheckOverlaps);
  new G4PVPlacement(0,scint_pos4,scintLV,"ScintPV4",scint_groupLV,false,4,fCheckOverlaps);
  new G4PVPlacement(0,scint_pos5,scintLV,"ScintPV5",scint_groupLV,false,5,fCheckOverlaps);
  new G4PVPlacement(0,scint_pos6,scintLV,"ScintPV6",scint_groupLV,false,6,fCheckOverlaps);
  new G4PVPlacement(0,scint_pos7,scintLV,"ScintPV7",scint_groupLV,false,7,fCheckOverlaps);
  new G4PVPlacement(0,scint_pos8,scintLV,"ScintPV8",scint_groupLV,false,8,fCheckOverlaps);
  new G4PVPlacement(0,scint_pos9,scintLV,"ScintPV9",scint_groupLV,false,9,fCheckOverlaps);

  auto scint_group_pos1 = G4ThreeVector((3.*scint_w+2.*dx-1.0*TPC_total_length)/2.,TPC_total_length/2.+scint_t/2.+dy,0.);
  auto scint_group_pos2 = G4ThreeVector((scint_t/2.+dy+TPC_total_length/2.),-(-TPC_total_length/2.+(3.*scint_w+2.*dx)/2.),0.);
  auto scint_group_pos3 = G4ThreeVector(-(3.*scint_w+2*dx-1.0*TPC_total_length)/2.,-(TPC_total_length/2.+scint_t/2.+dy),0.);
  auto scint_group_pos4 = G4ThreeVector(-(scint_t/2.+dy+TPC_total_length/2.),(-TPC_total_length/2.+(3.*scint_w+2.*dx)/2.),0.);

  G4RotationMatrix* zRot2 = new G4RotationMatrix; zRot2 -> rotateZ(90.*deg);
  G4RotationMatrix* zRot3 = new G4RotationMatrix; zRot3 -> rotateZ(180.*deg);
  G4RotationMatrix* zRot4 = new G4RotationMatrix; zRot4 -> rotateZ(270.*deg);

  auto scintPV1 = new G4PVPlacement(0,scint_group_pos1,scint_groupLV,"ScintPV_Group1",worldLV,false,1,fCheckOverlaps);
  auto scintPV2 = new G4PVPlacement(zRot2,scint_group_pos2,scint_groupLV,"ScintPV_Group2",worldLV,false,2,fCheckOverlaps);
  auto scintPV3 = new G4PVPlacement(zRot3,scint_group_pos3,scint_groupLV,"ScintPV_Group3",worldLV,false,3,fCheckOverlaps);
  auto scintPV4 = new G4PVPlacement(zRot4,scint_group_pos4,scint_groupLV,"ScintPV_Group4",worldLV,false,4,fCheckOverlaps);
  
  // - - - front and back scintillators
  G4double scint_w_fb = TPC_total_length/2. + dy + scint_t;
  G4double scint_h_fb = scint_w_fb - tube_thickness - tube_radius;
  G4double scint_w_fb2 = 2.*scint_w_fb - 2*scint_h_fb;
  
  // virtual volume for the layers
  auto scint_fb_S = new G4Box("ScintS",scint_w_fb/2.,scint_h_fb/2.,scint_t/2.);
  auto scint_fb_LV = new G4LogicalVolume(scint_fb_S,defaultMaterial,"ScintLV");
  auto scint_fb_S2 = new G4Box("ScintS",scint_h_fb/2.,scint_w_fb2/2.,scint_t/2.);
  auto scint_fb_LV2 = new G4LogicalVolume(scint_fb_S2,defaultMaterial,"ScintLV");
  
  // defining the layers
  auto scint_layer_fb_S = new G4Box("ScintS", scint_w_fb/2.,scint_h_fb/2., scint_layer_t/2.);
  auto scint_layer_fb_LV = new G4LogicalVolume(scint_layer_fb_S,scintMaterial,"Scint_layerLV");
  auto scint_layer_fb_S2 = new G4Box("ScintS", scint_h_fb/2.,scint_w_fb2/2., scint_layer_t/2.);
  auto scint_layer_fb_LV2 = new G4LogicalVolume(scint_layer_fb_S2,scintMaterial,"Scint_layerLV");
  
  auto scint_layer_fb_PV1 = new G4PVReplica("Scint_layerPV",scint_layer_fb_LV,scint_fb_LV,kZAxis,nofLayers,scint_layer_t);
  auto scint_layer_fb_PV2 = new G4PVReplica("Scint_layerPV",scint_layer_fb_LV2,scint_fb_LV2,kZAxis,nofLayers,scint_layer_t);
  
  // rotation is introduced to make a better indexing of the layers
  G4RotationMatrix * scint_fb_rot =  new G4RotationMatrix; scint_fb_rot -> rotateX(180.0*deg);

  auto scint_fb_pos1 = G4ThreeVector(scint_w_fb/2.,tube_radius+tube_thickness+scint_h_fb/2.,-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos2 = G4ThreeVector(-scint_w_fb/2.,tube_radius+tube_thickness+scint_h_fb/2.,-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos3 = G4ThreeVector(scint_w_fb/2.,-(tube_radius+tube_thickness+scint_h_fb/2.),-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos4 = G4ThreeVector(-scint_w_fb/2.,-(tube_radius+tube_thickness+scint_h_fb/2.),-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos5 = G4ThreeVector(-(tube_radius+tube_thickness+scint_h_fb/2.),0.,-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos6 = G4ThreeVector((tube_radius+tube_thickness+scint_h_fb/2.),0.,-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos7 = G4ThreeVector(scint_w_fb/2.,tube_radius+tube_thickness+scint_h_fb/2.,(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos8 = G4ThreeVector(-scint_w_fb/2.,tube_radius+tube_thickness+scint_h_fb/2.,(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos9 = G4ThreeVector(scint_w_fb/2.,-(tube_radius+tube_thickness+scint_h_fb/2.),(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos10 = G4ThreeVector(-scint_w_fb/2.,-(tube_radius+tube_thickness+scint_h_fb/2.),(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos11 = G4ThreeVector(-(tube_radius+tube_thickness+scint_h_fb/2.),0.,(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos12 = G4ThreeVector((tube_radius+tube_thickness+scint_h_fb/2.),0.,(3*scint_length+2*dz)/2.+scint_t/2.);
  
  auto scint_fb_PV1 = new G4PVPlacement(scint_fb_rot,scint_fb_pos1,scint_fb_LV,"ScintPVfb_Group1",worldLV,false,1,fCheckOverlaps);
  auto scint_fb_PV2 = new G4PVPlacement(scint_fb_rot,scint_fb_pos2,scint_fb_LV,"ScintPVfb_Group2",worldLV,false,2,fCheckOverlaps);
  auto scint_fb_PV3 = new G4PVPlacement(scint_fb_rot,scint_fb_pos3,scint_fb_LV,"ScintPVfb_Group3",worldLV,false,3,fCheckOverlaps);
  auto scint_fb_PV4 = new G4PVPlacement(scint_fb_rot,scint_fb_pos4,scint_fb_LV,"ScintPVfb_Group4",worldLV,false,4,fCheckOverlaps);
  auto scint_fb_PV5 = new G4PVPlacement(scint_fb_rot,scint_fb_pos5,scint_fb_LV2,"ScintPVfb_Group5",worldLV,false,5,fCheckOverlaps);
  auto scint_fb_PV6 = new G4PVPlacement(scint_fb_rot,scint_fb_pos6,scint_fb_LV2,"ScintPVfb_Group6",worldLV,false,6,fCheckOverlaps);
  auto scint_fb_PV7 = new G4PVPlacement(0,scint_fb_pos7,scint_fb_LV,"ScintPVfb_Group7",worldLV,false,7,fCheckOverlaps);
  auto scint_fb_PV8 = new G4PVPlacement(0,scint_fb_pos8,scint_fb_LV,"ScintPVfb_Group8",worldLV,false,8,fCheckOverlaps);
  auto scint_fb_PV9 = new G4PVPlacement(0,scint_fb_pos9,scint_fb_LV,"ScintPVfb_Group9",worldLV,false,9,fCheckOverlaps);
  auto scint_fb_PV10 = new G4PVPlacement(0,scint_fb_pos10,scint_fb_LV,"ScintPVfb_Group10",worldLV,false,10,fCheckOverlaps);
  auto scint_fb_PV11 = new G4PVPlacement(0,scint_fb_pos11,scint_fb_LV2,"ScintPVfb_Group11",worldLV,false,11,fCheckOverlaps);
  auto scint_fb_PV12 = new G4PVPlacement(0,scint_fb_pos12,scint_fb_LV2,"ScintPVfb_Group12",worldLV,false,12,fCheckOverlaps);

  // Lead Glass
  G4double lead_glass_xy = 8.*cm; G4double lead_glass_z = 25.*cm;
  G4double height_increment = 1.5*cm;

  auto lead_glass_y_level = TPC_total_length/2. +dy + scint_t;
  int lead_index = 0;
  
  auto absorber_groupS = new G4Box("AbsoS",lead_glass_y_level, (lead_glass_z+height_increment)/2., lead_glass_xy/2.);
  auto absorber_groupLV = new G4LogicalVolume(absorber_groupS,defaultMaterial,"Abso_GroupLV");
  
  auto absorberS = new G4Box("AbsoS",lead_glass_xy/2., lead_glass_z/2., lead_glass_xy/2.);
  auto absorberLV = new G4LogicalVolume(absorberS,absorberMaterial,"AbsoLV");
  

  ///***calculate the position directly from python
  std::vector<G4RotationMatrix *> rot_array_dir11;
  std::vector<G4RotationMatrix *> rot_array_dir21;
  std::vector<G4RotationMatrix *> rot_array_dir31;
  std::vector<G4RotationMatrix *> rot_array_dir41;
  
  for (int i = 0 ; i < data_non_axial.size();i++){rot_array_dir11.push_back(new G4RotationMatrix);}
  for (int i = 0 ; i < data_non_axial.size();i++){rot_array_dir21.push_back(new G4RotationMatrix);}
  for (int i = 0 ; i < data_non_axial.size();i++){rot_array_dir31.push_back(new G4RotationMatrix);}
  for (int i = 0 ; i < data_non_axial.size();i++){rot_array_dir41.push_back(new G4RotationMatrix);}

  G4double offset_lead_glass = 0.0*mm;

  // side 1 
  std::cout << std::setprecision(15) << lead_glass_y_level << std::endl;
  for (int i = 0 ; i < data_non_axial.size();i++){ // data_non_axial.size()

      rot_array_dir11[i]-> rotateX(-1.0*data_non_axial[i][3]*deg); rot_array_dir11[i]-> rotateZ(-1.0*data_non_axial[i][5]*deg);
      new G4PVPlacement(rot_array_dir11[i],
              G4ThreeVector(-1.0*data_non_axial[i][0]*cm, data_non_axial[i][1]*cm + lead_glass_y_level + offset_lead_glass,data_non_axial[i][2]*cm)
              ,absorberLV,"AbsoPV",worldLV,false,lead_index,fCheckOverlaps);
        
      lead_index ++ ;
  }

  // side 2 
  for (int i = 0 ; i < data_non_axial.size();i++){ // data_non_axial.size()

        rot_array_dir21[i]-> rotateY(data_non_axial[i][3]*deg); rot_array_dir21[i]-> rotateZ((-1.0*data_non_axial[i][5]+90.0)*deg);
        
        new G4PVPlacement(rot_array_dir21[i],
        G4ThreeVector(data_non_axial[i][1]*cm + lead_glass_y_level + offset_lead_glass,data_non_axial[i][0]*cm,data_non_axial[i][2]*cm)
        ,absorberLV,"AbsoPV",worldLV,false,lead_index,fCheckOverlaps);
    
        lead_index ++;
  }
  
  // side 3
  for (int i = 0 ; i < data_non_axial.size();i++){ // data_non_axial.size()

        rot_array_dir31[i]-> rotateX(data_non_axial[i][3]*deg); rot_array_dir31[i]-> rotateZ(data_non_axial[i][5]*deg+180.0*deg);
        
        new G4PVPlacement(rot_array_dir31[i],
        G4ThreeVector(-1.0*data_non_axial[i][0]*cm,-1.0*(data_non_axial[i][1]*cm + lead_glass_y_level+offset_lead_glass),data_non_axial[i][2]*cm)
        ,absorberLV,"AbsoPV",worldLV,false,lead_index,fCheckOverlaps);
        lead_index ++ ;        
  }

  // side 4
  for (int i = 0 ; i < data_non_axial.size();i++){ // data_non_axial.size()

        rot_array_dir41[i]-> rotateY(-data_non_axial[i][3]*deg); rot_array_dir41[i]-> rotateZ(-data_non_axial[i][5]*deg+270.0*deg);

        new G4PVPlacement(rot_array_dir41[i],
        G4ThreeVector(-(data_non_axial[i][1]*cm + lead_glass_y_level + offset_lead_glass),-data_non_axial[i][0]*cm,data_non_axial[i][2]*cm)
        ,absorberLV,"AbsoPV",worldLV,false,lead_index,fCheckOverlaps);
        lead_index ++ ;
  }
  
  
  // for the Front and back lead glass
  G4double lead_glass_y_level_fb = (3*scint_length+2*dz)/2.+scint_t;
  
  std::vector<G4RotationMatrix *> rot_array_dir111;
  std::vector<G4RotationMatrix *> rot_array_dir211;

  for (int i = 0 ; i < data_non_axial_fb.size();i++){rot_array_dir111.push_back(new G4RotationMatrix);}
  for (int i = 0 ; i < data_non_axial_fb.size();i++){rot_array_dir211.push_back(new G4RotationMatrix);}
  
  for (int i = 0 ; i < data_non_axial_fb.size();i++){
        rot_array_dir111[i]-> rotateX(data_non_axial_fb[i][3]*deg+90.*deg); rot_array_dir111[i]-> rotateZ(-data_non_axial_fb[i][5]*deg); //rot_array_dir111[i]-> rotateY(-data_non_axial_fb[i][5]*deg);
        new G4PVPlacement(rot_array_dir111[i],
        G4ThreeVector(data_non_axial_fb[i][0]*cm,data_non_axial_fb[i][2]*cm,(data_non_axial_fb[i][1]*cm + lead_glass_y_level_fb))
        ,absorberLV,"AbsoPV",worldLV,false,lead_index,fCheckOverlaps);
        lead_index ++ ;        
  }
  
  for (int i = 0 ; i < data_non_axial_fb.size();i++){
        rot_array_dir211[i]-> rotateX(-data_non_axial_fb[i][3]*deg-90.*deg); rot_array_dir211[i]-> rotateZ(-data_non_axial_fb[i][5]*deg); //rot_array_dir111[i]-> rotateY(-data_non_axial_fb[i][5]*deg);
        new G4PVPlacement(rot_array_dir211[i],
        G4ThreeVector(data_non_axial_fb[i][0]*cm,data_non_axial_fb[i][2]*cm,-(data_non_axial_fb[i][1]*cm + lead_glass_y_level_fb))
        ,absorberLV,"AbsoPV",worldLV,false,lead_index,fCheckOverlaps);
        lead_index ++ ;        
  }
  

  /***
  G4GDMLParser fParser;
  fParser.Write("geometry.gdml",worldPV,false);
  ***/
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  
  return worldPV;
}

//....

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
}
  /***
  // declare Scintillator as SinctillatorSD

  G4String scintDetectorName = "ScintLV" ;
  ScintillatorSD* scintDetector = new ScintillatorSD(scintDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);
  SetSensitiveDetector("ScintLV", scintDetector);
 
  // declare absorber as AbsorberSD
  G4String absorberDetectorName = "AbsoLV" ;
  AbsorberSD* absorberDetector = new AbsorberSD(absorberDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(absorberDetector);
  SetSensitiveDetector("AbsoLV", absorberDetector);

  // declare vacuum as TubeSD
  G4String tubeDetectorName = "TubeLV" ;
  TubeSD* tubeDetector = new TubeSD(tubeDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(tubeDetector);
  SetSensitiveDetector("TubeLV", tubeDetector);

  // declare vacuum as TPCSD
  G4String TPCDetectorName = "TPCLV" ;
  TPCSD* TPCDetector = new TPCSD(TPCDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(TPCDetector);
  SetSensitiveDetector("gasLV", TPCDetector);
  ***/


//....

