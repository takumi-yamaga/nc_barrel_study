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
/// \file FNCParameterisation.cc
/// \brief Implementation of the FNCParameterisation class

#include "FNCParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FNCParameterisation::FNCParameterisation(  
    G4int number_of_segments,
    G4double width
    )
: G4VPVParameterisation()
{
  number_of_segments_ =  number_of_segments;
  width_ = width;
  if(number_of_segments_%2){ // odd
  offset_ = -width_*((G4double)number_of_segments_-1.)/2.;
  }
  else{ // even
  offset_ = -width_*((G4double)number_of_segments_/2.-0.5);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FNCParameterisation::~FNCParameterisation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FNCParameterisation::ComputeTransformation
(const G4int copy_number, G4VPhysicalVolume* physical) const
{
  // copy_number will start with zero.
  G4double x_position = offset_ + width_ * (G4double)copy_number;
  G4double y_position = 0.;
  G4double z_position = 0.;
  G4ThreeVector current_position(x_position,y_position,z_position);
  physical->SetTranslation(current_position);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

