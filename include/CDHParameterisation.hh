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
/// \file CDHParameterisation.hh
/// \brief Definition of the CDHParameterisation class

#ifndef CDHParameterisation_h
#define CDHParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

class G4VPhysicalVolume;
class G4Box;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Tubs;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

///  A parameterisation that describes a series of boxes along Z.
///
///  The boxes have equal width, & their lengths are a linear equation.
///  They are spaced an equal distance apart, starting from given location.


// static constant
static const int kCDHSpace = 1.*mm;


class CDHParameterisation : public G4VPVParameterisation
{ 
  public:

    CDHParameterisation(
        G4int    number_of_segment,
        G4ThreeVector position,
        G4double radius,
        G4double thickness,
        G4double length
        );

    virtual ~CDHParameterisation();

    void ComputeTransformation (const G4int copy_number, G4VPhysicalVolume* physical) const;

    void ComputeDimensions (G4Box& box, const G4int copy_number,
        const G4VPhysicalVolume* physical) const;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Tubs&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trd&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
        const G4VPhysicalVolume*) const {}

  private:

    G4int number_of_segments_;   
    G4ThreeVector position_;
    G4double radius_;
    G4double thickness_;
    G4double length_;
    G4double d_phi_;
    G4double width_;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
