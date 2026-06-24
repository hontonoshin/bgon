#pragma once

#include "NeutronSimHit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4HCofThisEvent.hh"

// CrystalSD covers both Crystal and BiCoating logical volumes.
// Birth volume is filled via G4Track::GetOriginTouchable() -- no
// TrackingAction or static map required.
class CrystalSD : public G4VSensitiveDetector
{
public:
    CrystalSD(const G4String& name, const G4String& hcName);
    ~CrystalSD() override = default;

    void   Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* rOhist) override;
    void   EndOfEvent(G4HCofThisEvent*) override;

private:
    NeutronSimHitsCollection* fHitsCollection = nullptr;
    G4int                     fHCID           = -1;
};
