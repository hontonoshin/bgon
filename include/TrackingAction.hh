#pragma once
#include "G4UserTrackingAction.hh"

// Minimal TrackingAction -- kept for future extension.
// Birth-volume identification is now handled directly in CrystalSD
// via G4Track::GetOriginTouchable(), which is simpler and thread-safe.
class TrackingAction : public G4UserTrackingAction
{
public:
    TrackingAction()  = default;
    ~TrackingAction() override = default;
    void PreUserTrackingAction(const G4Track*) override {}
    void PostUserTrackingAction(const G4Track*) override {}
};
