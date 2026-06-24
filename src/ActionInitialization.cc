#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"

void ActionInitialization::BuildForMaster() const
{
    // Master RunAction books the same histograms and ntuples as workers,
    // allowing Geant4 ROOT histogram and ntuple merging.
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new RunAction());
    SetUserAction(new EventAction());

    // Kept for source compatibility and future track-level extensions.
    // Current birth-volume tagging is performed directly in CrystalSD via
    // G4Track::GetOriginTouchable().
    SetUserAction(new TrackingAction());
}
