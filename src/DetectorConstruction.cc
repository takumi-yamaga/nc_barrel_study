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
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "CDHParameterisation.hh"
#include "HodoscopeSD.hh"
#include "Constants.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), 
  cdh_logical_(nullptr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
  //delete fMessenger;

  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Construct materials
  ConstructMaterials();
  auto air = G4Material::GetMaterial("G4_AIR");
  auto scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto liquid_he3 = new G4Material("Liquid_HeThree",2.,3.016029*g/mole,81.2*mg/cm3,kStateLiquid);
  //auto vacuum = G4Material::GetMaterial("G4_Galactic");
  //auto argonGas = G4Material::GetMaterial("G4_Ar");
  //auto csI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  //auto lead = G4Material::GetMaterial("G4_Pb");
  //auto carbon = G4Material::GetMaterial("G4_C");
  //auto hydrogen = new G4Material("hydrogne", 1., 1.01*g/mole, 1.*g/cm3);

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool kCheckOverlaps = true;

  // geometries --------------------------------------------------------------

  // experimental hall (world volume) -------------------------------
  auto world_solid 
    = new G4Box("world_solid",10.*m,10.*m,10.*m);
  auto world_logical
    = new G4LogicalVolume(world_solid,air,"world_logical");
  auto world_physical
    = new G4PVPlacement(0,G4ThreeVector(),world_logical,"world_physical",0,
        false,0,kCheckOverlaps);
  // ----------------------------------------------------------------

  // target (tube) --------------------------------------------------
  auto target_radius = 35*mm;
  auto target_length   = 150*mm;
  auto target_position = G4ThreeVector(0.*mm,0.*mm,0.*mm);
  auto target_solid
    = new G4Tubs("target_solid",0.*mm,target_radius,target_length,0.*deg,360.*deg);
  auto target_logical
    = new G4LogicalVolume(target_solid,liquid_he3,"target_logical");
  //new G4PVPlacement(0,target_position,target_logical,"target_physical",
  //    world_logical,false,0,kCheckOverlaps);
  // ----------------------------------------------------------------
  
  // cdh mother volume ----------------------------------------------
  auto cdh_mother_radius = 3*m;
  auto cdh_mother_length   = 3*m;
  auto cdh_mother_position = G4ThreeVector(0.*mm,0.*mm,0.*mm);
  auto cdh_mother_solid
    = new G4Tubs("cdh_mother_solid",0.*mm,cdh_mother_radius,cdh_mother_length,0.*deg,360.*deg);
  auto cdh_mother_logical
    = new G4LogicalVolume(cdh_mother_solid,liquid_he3,"cdh_mother_logical");
  new G4PVPlacement(0,cdh_mother_position,cdh_mother_logical,"cdh_mother_physical",
      world_logical,false,0,kCheckOverlaps);
  // ----------------------------------------------------------------

  // cdh ------------------------------------------------------------
  auto cdh_number_of_segments = 36;
  auto cdh_position  = G4ThreeVector(0.*mm,0.*mm,0.*mm);
  auto cdh_radius    = 80.*cm;
  auto cdh_thickness = 3.*cm;
  auto cdh_length    = 100.*cm;

  //dummy solid volume for G4Box (modified by parameterised volume)
  auto cdh_solid 
    = new G4Box("cdh_solid",cdh_thickness/2.,1.*cm/2.,cdh_length/2.);
  cdh_logical_
    = new G4LogicalVolume(cdh_solid,scintillator,"cdh_logical");
  G4VPVParameterisation* cdh_parametrization
    = new CDHParameterisation(cdh_number_of_segments,cdh_position,
        cdh_radius,cdh_thickness,cdh_length);
  // kZAxis is dummy for CDH.
  new G4PVParameterised("cdh_physical",cdh_logical_,
      cdh_mother_logical,kZAxis,cdh_number_of_segments,cdh_parametrization,kCheckOverlaps);
  // ----------------------------------------------------------------

  // visualization attributes ------------------------------------------------
  auto visAttributes = new G4VisAttributes(G4Colour::White());
  visAttributes->SetVisibility(false);
  world_logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(MyColour::Scintillator());
  cdh_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  // return the world physical volume ----------------------------------------
  return world_physical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String sensitive_detector_name;

  // sensitive detectors -----------------------------------------------------
  auto cdh = new HodoscopeSD(sensitive_detector_name="/cdh");
  sdManager->AddNewDetector(cdh);
  cdh_logical_->SetSensitiveDetector(cdh);

}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  nistManager->FindOrBuildMaterial("G4_AIR");

  // Argon gas
  nistManager->FindOrBuildMaterial("G4_Ar");

  // Scintillator
  // (PolyVinylToluene, C_9H_10)
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // CsI
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");

  // Lead
  nistManager->FindOrBuildMaterial("G4_Pb");

  // Carbon
  nistManager->FindOrBuildMaterial("G4_C");

  // Vacuum "Galactic"
  nistManager->FindOrBuildMaterial("G4_Galactic");


  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
