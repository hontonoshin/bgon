#pragma once
// BiCoatingSD.hh — sensitive detector for Bi-209 coating

#include "G4VSensitiveDetector.hh"
#include "NeutronSimHit.hh"

class BiCoatingSD : public G4VSensitiveDetector
{
public:
    BiCoatingSD(const G4String& name, const G4String& hcName);
    ~BiCoatingSD() override = default;

    void   Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;
    void   EndOfEvent(G4HCofThisEvent* hce) override;

private:
    NeutronSimHitsCollection* fHitsCollection = nullptr;
    G4int fHCID = -1;
};
