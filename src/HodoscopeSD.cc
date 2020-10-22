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
/// \copied from B5HodoscopeSD.cc
/// \brief Implementation of the Hodoscope class

#include "HodoscopeSD.hh"
#include "HodoscopeHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HodoscopeSD::HodoscopeSD(G4String name)
: G4VSensitiveDetector(name), 
  hits_collection_(nullptr), hits_collection_id_(-1)
{
  collectionName.insert("hodoscope_hitscollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HodoscopeSD::~HodoscopeSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HodoscopeSD::Initialize(G4HCofThisEvent* collection)
{
  hits_collection_
    = new HodoscopeHitsCollection(SensitiveDetectorName,collectionName[0]);

  if (hits_collection_id_<0) { 
    hits_collection_id_ = G4SDManager::GetSDMpointer()->GetCollectionID(hits_collection_); 
  }
  collection->AddHitsCollection(hits_collection_id_,hits_collection_);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HodoscopeSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto energy_deposit = step->GetTotalEnergyDeposit();
  if(energy_deposit==0.) return true;

  auto track = step->GetTrack();
  auto pre_steppoint = step->GetPreStepPoint();
  auto touchable = pre_steppoint->GetTouchable();
  auto physical = touchable->GetVolume(0);
  auto transform = touchable->GetHistory()->GetTopTransform();
  auto logical = physical->GetLogicalVolume();

  auto track_id = track->GetTrackID();
  auto parent_id = track->GetParentID();
  auto particle_id = track->GetParticleDefinition()->GetPDGEncoding();
  auto segment_id = physical->GetCopyNo();
  auto hit_time = pre_steppoint->GetGlobalTime();
  auto global_position = pre_steppoint->GetPosition();
  auto local_position = transform.TransformPoint(global_position);
  auto momentum = pre_steppoint->GetMomentum();
  auto polarization = track->GetPolarization();

  auto position = transform.NetTranslation();
  auto rotation = transform.NetRotation();

  G4bool segment_has_hit = false;
  for(auto i_hit=0; i_hit<hits_collection_->entries(); i_hit++){
    auto hit = (*hits_collection_)[i_hit];
    if(hit->GetSegmentID()==segment_id){
      segment_has_hit = true;
      hit->PushTotalHits();
      hit->PushTrackID(track_id);
      hit->PushParentID(parent_id);
      hit->PushParticleID(particle_id);
      hit->PushHitTime(hit_time);
      hit->PushEnergyDeposit(energy_deposit);
      hit->PushGlobalPosition(global_position);
      hit->PushLocalPosition(local_position);
      hit->PushMomentum(momentum);
      hit->PushPolarization(polarization);
    }
  }
  // if there is no hit in the segment, create new hit.
  if(!segment_has_hit){
    auto hit = new HodoscopeHit();
    hit->SetSegmentID(segment_id);
    hit->SetLogicalVolume(logical);
    hit->SetPosition(position);
    hit->SetRotation(rotation);
    hit->PushTotalHits();
    hit->PushTrackID(track_id);
    hit->PushParentID(parent_id);
    hit->PushParticleID(particle_id);
    hit->PushHitTime(hit_time);
    hit->PushEnergyDeposit(energy_deposit);
    hit->PushGlobalPosition(global_position);
    hit->PushLocalPosition(local_position);
    hit->PushMomentum(momentum);
    hit->PushPolarization(polarization);
    hits_collection_->insert(hit);
  }

  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
