#pragma once

#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "globals.hh"

#include <vector>

struct SecondaryInfo
{
    G4int    pdg            = 0;
    G4int    trackID        = -1;
    G4int    parentID       = -1;
    G4int    Z              = 0;
    G4int    A              = 0;
    G4double kineticEnergy  = 0.0; // stored numerically in MeV
    G4String name;
    G4String creatorProcess;
};

class NeutronSimHit : public G4VHit
{
public:
    NeutronSimHit()                                = default;
    ~NeutronSimHit() override                      = default;
    NeutronSimHit(const NeutronSimHit&)            = default;
    NeutronSimHit& operator=(const NeutronSimHit&) = default;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // ---- step kinematics ------------------------------------------------
    G4double      edep       = 0.0;
    G4double      stepLength = 0.0;

    // position is retained as the post-step position for backward
    // compatibility.  The explicit pre/post positions allow EventAction to
    // use the step midpoint for shower moments and the post-step point for
    // interaction vertices.
    G4ThreeVector position;
    G4ThreeVector prePosition;
    G4ThreeVector postPosition;

    // ---- track identity -------------------------------------------------
    G4int  copyNo   = -1;
    G4int  trackID  = -1;
    G4int  parentID = -1;
    G4int  pdg      = 0;
    G4int  Z        = 0;
    G4int  A        = 0;

    G4double preKineticEnergy  = 0.0;
    G4double postKineticEnergy = 0.0;
    G4double globalTime        = 0.0; // post-step global time

    G4String creatorProcess;
    G4String processName;

    // ---- volume tags ----------------------------------------------------
    // 1 = crystal, 2 = BiCoating, 0 = unknown/world
    G4int originVolume = 0;

    // Birth volume is determined directly from G4Track::GetOriginTouchable().
    // 1 = crystal, 2 = BiCoating, 0 = elsewhere/primary.
    G4int birthVolume = 0;

    // ---- step-level flags -----------------------------------------------
    G4bool firstStepInVolume = false;
    G4bool lastStepInVolume  = false;
    G4bool interactionStep   = false;

    // ---- secondaries produced at this step ------------------------------
    std::vector<SecondaryInfo> secondaries;
};

using NeutronSimHitsCollection = G4THitsCollection<NeutronSimHit>;

extern G4ThreadLocal G4Allocator<NeutronSimHit>* gNeutronSimHitAllocator;

inline void* NeutronSimHit::operator new(size_t)
{
    if (!gNeutronSimHitAllocator)
        gNeutronSimHitAllocator = new G4Allocator<NeutronSimHit>;
    return gNeutronSimHitAllocator->MallocSingle();
}

inline void NeutronSimHit::operator delete(void* hit)
{
    gNeutronSimHitAllocator->FreeSingle(
        static_cast<NeutronSimHit*>(hit));
}
