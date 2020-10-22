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
#include "FNCParameterisation.hh"
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
  fnc_logical_(nullptr)
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
    = new G4Box("world_solid",1.5*m,1.5*m,1.5*m);
  auto world_logical
    = new G4LogicalVolume(world_solid,air,"world_logical");
  auto world_physical
    = new G4PVPlacement(0,G4ThreeVector(),world_logical,"world_physical",0,
        false,0,kCheckOverlaps);
  // ----------------------------------------------------------------



  // ======================================================
  // fiber nutron counter =================================
  // ======================================================
  // parameters for fnc segment ------------------------------
  auto fnc_number_of_layers_in_stack = 2;
  auto fnc_number_of_segments_in_layer = 101;
  auto fnc_fiber_width = 3.*mm;
  auto fnc_fiber_thickness = 3.*mm;
  auto fnc_fiber_length = 50.*cm;
  // ---------------------------------------------------------

  // fnc stack volume ---------------------------------------
  auto fnc_stack_width   = (G4double)fnc_number_of_segments_in_layer * fnc_fiber_width;
  auto fnc_stack_thickness   = (fnc_fiber_thickness*fnc_number_of_layers_in_stack);
  auto fnc_stack_length  = fnc_fiber_length;
  auto fnc_stack_position = G4ThreeVector(0.*mm,0.*mm,0.*mm);
  auto fnc_stack_solid
    = new G4Box("fnc_stack_solid",fnc_stack_width/2.,fnc_stack_length/2.,fnc_stack_thickness/2.);
  auto fnc_stack_logical
    = new G4LogicalVolume(fnc_stack_solid,scintillator,"fnc_stack_logical");
  new G4PVPlacement(0,fnc_stack_position,fnc_stack_logical,"fnc_stack_0",
      world_logical,false,0,kCheckOverlaps);
  // ----------------------------------------------------------------

  // fnc layer volume ---------------------------------------
  auto fnc_layer_width   = fnc_stack_width - kSpace;
  auto fnc_layer_thickness   = fnc_fiber_thickness - kSpace;
  auto fnc_layer_length  = fnc_stack_length - kSpace;
  auto fnc_layer_offset = 0.;
  auto fnc_layer_solid
    = new G4Box("fnc_layer_solid",fnc_layer_width/2.,fnc_layer_length/2.,fnc_layer_thickness/2.);
  auto fnc_layer_logical
    = new G4LogicalVolume(fnc_layer_solid,scintillator,"fnc_layer_logical");

  if(fnc_number_of_layers_in_stack%2){ // odd
    fnc_layer_offset = - (G4double)(fnc_number_of_layers_in_stack-1)/2. * fnc_fiber_thickness;
  }
  else{ // even
    fnc_layer_offset = - ((G4double)fnc_number_of_layers_in_stack/2. - 0.5) * fnc_fiber_thickness;
  }
  // along Z axis 
  for(auto i_layer = 0; i_layer<fnc_number_of_layers_in_stack; ++i_layer){
    auto z_position = fnc_layer_offset + fnc_fiber_thickness * (G4double)i_layer;
    G4cout << "z_position ::: " << z_position << G4endl;
    auto fnc_layer_position = G4ThreeVector(0.*mm,0.*mm,z_position);
    char fnc_layer_physical_name[50];
    std::sprintf(fnc_layer_physical_name,"fnc_layer_%d",i_layer);
    new G4PVPlacement(0,fnc_layer_position,fnc_layer_logical,fnc_layer_physical_name,
        fnc_stack_logical,false,0,kCheckOverlaps);
  }
  // ----------------------------------------------------------------

  //// fnc segment ----------------------------------------------------
  auto fnc_segment_width   = fnc_fiber_width - kSpace;
  auto fnc_segment_thickness   = fnc_layer_thickness - kSpace;
  auto fnc_segment_length  = fnc_layer_length - kSpace;
  auto fnc_solid 
    = new G4Box("fnc_solid",fnc_segment_width/2.,fnc_segment_length/2.,fnc_segment_thickness/2.);
  fnc_logical_
    = new G4LogicalVolume(fnc_solid,scintillator,"fnc_logical");
  G4VPVParameterisation* fnc_parametrization
    = new FNCParameterisation(fnc_number_of_segments_in_layer,fnc_segment_width);
  new G4PVParameterised("fnc_physical",fnc_logical_,
      fnc_layer_logical,kXAxis,fnc_number_of_segments_in_layer,fnc_parametrization,kCheckOverlaps);
  // kXAxis is dummy.
  //// ----------------------------------------------------------------

  // visualization attributes ------------------------------------------------
  auto visAttributes = new G4VisAttributes(G4Colour::White());
  visAttributes->SetVisibility(false);
  world_logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(MyColour::Scintillator());
  fnc_logical_->SetVisAttributes(visAttributes);
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
  auto fnc = new HodoscopeSD(sensitive_detector_name="/fnc");
  sdManager->AddNewDetector(fnc);
  fnc_logical_->SetSensitiveDetector(fnc);

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
