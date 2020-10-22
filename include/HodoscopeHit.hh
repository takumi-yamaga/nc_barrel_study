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
/// \copied from B5HodoscopeHit.hh
/// \brief Definition of the HodoscopeHit class

#ifndef HodoscopeHit_h
#define HodoscopeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include <vector>

class G4AttDef;
class G4AttValue;

using std::vector;

/// Drift chamber hit
///
/// It records:
/// - the layer ID
/// - the particle time
/// - the particle local and global positions

class HodoscopeHit : public G4VHit
{
  public:
    HodoscopeHit();
    HodoscopeHit(const HodoscopeHit &right);
    virtual ~HodoscopeHit();

    const HodoscopeHit& operator=(const HodoscopeHit &right);
    G4bool operator==(const HodoscopeHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

    inline void SetSegmentID(const G4int id) { segment_id_ = id; }
    inline G4int GetSegmentID() const { return segment_id_; }

    inline void SetLogicalVolume(G4LogicalVolume* logical) { logical_ = logical; }
    inline G4LogicalVolume* GetLogicalVolume() const { return logical_; }

    inline void SetPosition(const G4ThreeVector position) { position_ = position; }
    inline G4ThreeVector GetPosition() const { return position_; }

    inline void SetRotation(const G4RotationMatrix rotation) { rotation_ = rotation; }
    inline G4RotationMatrix GetRotation() const { return rotation_; }

    inline void PushTotalHits() { total_hits_++; }
    inline void SetTotalHits(const G4int total_hits) { total_hits_ = total_hits; }
    inline G4int GetTotalHits() const { return total_hits_; }

    inline void PushTrackID(const G4int id) { track_id_.push_back(id); }
    inline void SetTrackID(const G4int i_hit, const G4int id);
    inline vector<G4int> GetTrackID() const { return track_id_; }
    inline G4int GetTrackID(const G4int i_hit) const;

    inline void PushParentID(const G4int id) { parent_id_.push_back(id); }
    inline void SetParentID(const G4int i_hit, const G4int id);
    inline vector<G4int> GetParentID() const { return parent_id_; }
    inline G4int GetParentID(const G4int i_hit) const;

    inline void PushParticleID(const G4int id) { particle_id_.push_back(id); }
    inline void SetParticleID(const G4int i_hit, const G4int id);
    inline vector<G4int> GetParticleID() const { return particle_id_; }
    inline G4int GetParticleID(const G4int i_hit) const;

    inline void PushHitTime(const G4double time) { hit_time_.push_back(time); }
    inline void SetHitTime(const G4int i_hit, G4double time);
    inline vector<G4double> GetHitTime() const { return hit_time_; }
    inline G4double GetHitTime(const G4int i_hit) const;

    inline void PushEnergyDeposit(const G4double de) { energy_deposit_.push_back(de); }
    inline void SetEnergyDeposit(const G4int i_hit, const G4double de);
    inline vector<G4double> GetEnergyDeposit() const { return energy_deposit_; }
    inline G4double GetEnergyDeposit(const G4int i_hit) const;

    inline void PushLocalPosition(const G4ThreeVector position) { local_position_.push_back(position); }
    inline void SetLocalPosition(const G4int i_hit, const G4ThreeVector position);
    inline vector<G4ThreeVector> GetLocalPosition() const { return local_position_; }
    inline G4ThreeVector GetLocalPosition(const G4int i_hit) const;

    inline void PushGlobalPosition(const G4ThreeVector position) { global_position_.push_back(position); }
    inline void SetGlobalPosition(const G4int i_hit, const G4ThreeVector position);
    inline vector<G4ThreeVector> GetGlobalPosition() const { return global_position_; }
    inline G4ThreeVector GetGlobalPosition(const G4int i_hit) const;

    inline void PushMomentum(const G4ThreeVector momentum) { momentum_.push_back(momentum); }
    inline void SetMomentum(const G4int i_hit, const G4ThreeVector momentum);
    inline vector<G4ThreeVector> GetMomentum() const { return momentum_; }
    inline G4ThreeVector GetMomentum(const G4int i_hit) const;

    inline void PushPolarization(const G4ThreeVector polarization) { polarization_.push_back(polarization); }
    inline void SetPolarization(const G4int i_hit, const G4ThreeVector polarization);
    inline vector<G4ThreeVector> GetPolarization() const { return polarization_; }
    inline G4ThreeVector GetPolarization(const G4int i_hit) const;

  private:
    G4int segment_id_;
    G4LogicalVolume* logical_;
    G4ThreeVector position_;
    G4RotationMatrix rotation_;
    G4int total_hits_;

    vector<G4int> track_id_;
    vector<G4int> parent_id_;
    vector<G4int> particle_id_;
    vector<G4double> hit_time_;
    vector<G4double> energy_deposit_;
    vector<G4ThreeVector> local_position_;
    vector<G4ThreeVector> global_position_;
    vector<G4ThreeVector> momentum_;
    vector<G4ThreeVector> polarization_;
};

using HodoscopeHitsCollection = G4THitsCollection<HodoscopeHit>;

extern G4ThreadLocal G4Allocator<HodoscopeHit>* HodoscopeHitAllocator;

inline void* HodoscopeHit::operator new(size_t)
{
  if (!HodoscopeHitAllocator) {
    HodoscopeHitAllocator = new G4Allocator<HodoscopeHit>;
  }
  return (void*)HodoscopeHitAllocator->MallocSingle();
}

inline void HodoscopeHit::operator delete(void* aHit)
{
  HodoscopeHitAllocator->FreeSingle((HodoscopeHit*) aHit);
}

inline void HodoscopeHit::SetTrackID(const G4int i_hit, const G4int id){
  if(i_hit<total_hits_){
    track_id_[i_hit] = id;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetTrackID()",
        "Code002", JustWarning, msg);
  }
}
inline G4int HodoscopeHit::GetTrackID(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return track_id_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetTrackID(G4int)",
        "Code002", JustWarning, msg);
    return -1;
  }
}

inline void HodoscopeHit::SetParentID(const G4int i_hit, const G4int id){
  if(i_hit<total_hits_){
    parent_id_[i_hit] = id;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetParentID()",
        "Code002", JustWarning, msg);
  }
}
inline G4int HodoscopeHit::GetParentID(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return parent_id_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetParentID(G4int)",
        "Code002", JustWarning, msg);
    return -1;
  }
}

inline void HodoscopeHit::SetParticleID(const G4int i_hit, const G4int id){
  if(i_hit<total_hits_){
    parent_id_[i_hit] = id;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetParticleID()",
        "Code002", JustWarning, msg);
  }
}
inline G4int HodoscopeHit::GetParticleID(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return particle_id_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetParticleID(G4int)",
        "Code002", JustWarning, msg);
    return -1;
  }
}

inline void HodoscopeHit::SetHitTime(const G4int i_hit, const G4double hit_time){
  if(i_hit<total_hits_){
    hit_time_[i_hit] = hit_time;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetHitTime()",
        "Code002", JustWarning, msg);
  }
}
inline G4double HodoscopeHit::GetHitTime(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return hit_time_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetHitTime(G4int)",
        "Code002", JustWarning, msg);
    return -1;
  }
}

inline void HodoscopeHit::SetEnergyDeposit(const G4int i_hit, const G4double energy_deposit){
  if(i_hit<total_hits_){
    energy_deposit_[i_hit] = energy_deposit;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetEnergyDeposit()",
        "Code002", JustWarning, msg);
  }
}
inline G4double HodoscopeHit::GetEnergyDeposit(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return energy_deposit_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetEnergyDeposit(G4int)",
        "Code002", JustWarning, msg);
    return -1;
  }
}

inline void HodoscopeHit::SetLocalPosition(const G4int i_hit, const G4ThreeVector local_position){
  if(i_hit<total_hits_){
    local_position_[i_hit] = local_position;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetLocalPosition()",
        "Code002", JustWarning, msg);
  }
}
inline G4ThreeVector HodoscopeHit::GetLocalPosition(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return local_position_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetLocalPosition(G4int)",
        "Code002", JustWarning, msg);
    return G4ThreeVector(0);
  }
}

inline void HodoscopeHit::SetGlobalPosition(const G4int i_hit, const G4ThreeVector global_position){
  if(i_hit<total_hits_){
    global_position_[i_hit] = global_position;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetGlobalPosition()",
        "Code002", JustWarning, msg);
  }
}
inline G4ThreeVector HodoscopeHit::GetGlobalPosition(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return global_position_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetGlobalPosition(G4int)",
        "Code002", JustWarning, msg);
    return G4ThreeVector(0);
  }
}

inline void HodoscopeHit::SetMomentum(const G4int i_hit, const G4ThreeVector momentum){
  if(i_hit<total_hits_){
    momentum_[i_hit] = momentum;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetMomentum()",
        "Code002", JustWarning, msg);
  }
}
inline G4ThreeVector HodoscopeHit::GetMomentum(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return momentum_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetMomentum(G4int)",
        "Code002", JustWarning, msg);
    return G4ThreeVector(0);
  }
}

inline void HodoscopeHit::SetPolarization(const G4int i_hit, const G4ThreeVector polarization){
  if(i_hit<total_hits_){
    polarization_[i_hit] = polarization;
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::SetPolarization()",
        "Code002", JustWarning, msg);
  }
}
inline G4ThreeVector HodoscopeHit::GetPolarization(const G4int i_hit) const{
  if(i_hit<total_hits_){
    return polarization_[i_hit];
  }
  else{
    G4ExceptionDescription msg;
    msg << "No hits found." << G4endl; 
    G4Exception("HodoscopeHit::GetPolarization(G4int)",
        "Code002", JustWarning, msg);
    return G4ThreeVector(0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
