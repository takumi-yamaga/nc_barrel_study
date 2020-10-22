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
/// \copied from B5Constants.hh
/// \brief Definition of constants.

#ifndef Constants_h
#define Constants_h 1

#include <array>
#include "globals.hh"
#include "G4Colour.hh"

using std::array;

// hodoscopes
namespace Hodoscope{
  constexpr G4int kTotalNumber = 2;
  const array<G4String, kTotalNumber> detector_name
    = {{ "hodoscope1", "hodoscope2" }};
}

namespace MyColour{
  // G4Colour(red, green, blue, alpha)
  // alpha = 1. - transparency
  
  inline G4Colour Scintillator(){ return G4Colour(0.2,1.0,1.0,0.2); }
  inline G4Colour ScintillatorHasHit(){ return G4Colour(1.0,0.0,0.0,0.2); }

  inline G4Colour Hit(){ return G4Colour(1.0,0.0,0.0,1.0); }
}

#endif
