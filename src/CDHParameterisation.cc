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
/// \file CDHParameterisation.cc
/// \brief Implementation of the CDHParameterisation class

#include "CDHParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CDHParameterisation::CDHParameterisation(  
    G4int    number_of_segments, 
    G4ThreeVector position,
    G4double radius,
    G4double thickness,
    G4double length
    )
: G4VPVParameterisation()
{
  number_of_segments_ =  number_of_segments;
  position_ = position;
  radius_ = radius;
  thickness_ = thickness;
  length_ = length;
  d_phi_ = 360.*deg/(G4double)number_of_segments;
  width_ = (2.*radius-thickness)*std::tan(d_phi_/2.) - 2.*kCDHSpace;

  if( width_<0 ){
    G4Exception("CDHParameterisation::CDHParameterisation()",
        "InvalidSetup", FatalException,
        "too much cdh segments");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CDHParameterisation::~CDHParameterisation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDHParameterisation::ComputeTransformation
(const G4int copy_number, G4VPhysicalVolume* physical) const
{
  // copy_number will start with zero.
  G4double phi_position = d_phi_ * (G4double)copy_number;
  G4double x_position = position_.x() + radius_*std::cos(phi_position);
  G4double y_position = position_.y() + radius_*std::sin(phi_position);
  G4double z_position = position_.z();
  G4ThreeVector current_position(x_position,y_position,z_position);
  G4RotationMatrix* current_rotation = new G4RotationMatrix();
  current_rotation->rotateZ(-phi_position-90.*deg);
  physical->SetTranslation(current_position);
  physical->SetRotation(current_rotation);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CDHParameterisation::ComputeDimensions
(G4Box& box, const G4int copy_number, const G4VPhysicalVolume*) const
{
  // cdh size is adjusted to number of segments.
  // all segments have the same structure and size.
  box.SetXHalfLength(width_/2.);
  box.SetYHalfLength(thickness_/2.);
  box.SetZHalfLength(length_/2.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
