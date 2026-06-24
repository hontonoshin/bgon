#include "BiCoatingSD.hh"
#include "NeutronSimHit.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

// ---------------------------------------------------------------------------
namespace {

// Birth-volume tag — mirrors CrystalSD exactly.
G4int BirthVolumeTag(const G4Track* track)
{
    if (!track) return 0;
    const G4VTouchable* originTouch = track->GetOriginTouchable();
    if (!originTouch) return 0;
    const G4VPhysicalVolume* pv = originTouch->GetVolume();
    if (!pv || !pv->GetLogicalVolume()) return 0;
    const auto& name = pv->GetLogicalVolume()->GetName();
    if (name == "Crystal")   return 1;
    if (name == "BiCoating") return 2;
    return 0;
}

bool IsTransportation(const G4VProcess* process)
{
    return !process ||
           process->GetProcessName() == "Transportation";
}

} // namespace

// ---------------------------------------------------------------------------
BiCoatingSD::BiCoatingSD(const G4String& name, const G4String& hcName)
    : G4VSensitiveDetector(name)
{
    collectionName.insert(hcName);
}

// ---------------------------------------------------------------------------
void BiCoatingSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new NeutronSimHitsCollection(
                          SensitiveDetectorName, collectionName[0]);
    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()
                    ->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(fHCID, fHitsCollection);
}

// ---------------------------------------------------------------------------
G4bool BiCoatingSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    if (!step) return false;

    const auto* track = step->GetTrack();
    const auto* pre   = step->GetPreStepPoint();
    const auto* post  = step->GetPostStepPoint();
    if (!track || !pre || !post || !track->GetDefinition()) return false;

    const G4double edep            = step->GetTotalEnergyDeposit();
    const auto*    process         = post->GetProcessDefinedStep();
    const G4bool   interactionStep = !IsTransportation(process);
    const auto*    created         = step->GetSecondaryInCurrentStep();
    const G4bool   hasSecondaries  = created && !created->empty();
    const G4bool   firstStep       = step->IsFirstStepInVolume();
    const G4bool   lastStep        = step->IsLastStepInVolume();

    // Keep step if it deposits energy, is an interaction, produced secondaries,
    // or is the first/last step in the volume (for track entry/exit accounting).
    // CRITICAL: do NOT drop zero-edep interaction steps — hadronic inelastic
    // reactions often deposit zero energy on the step that creates fragments.
    if (edep <= 0.0 && !interactionStep && !hasSecondaries
        && !firstStep && !lastStep)
        return false;

    auto* hit = new NeutronSimHit();

    // ---- Energy / geometry -----------------------------------------------
    hit->edep       = edep;
    hit->stepLength = step->GetStepLength();
    hit->position   = post->GetPosition();
    hit->copyNo     = pre->GetTouchableHandle()
                          ? pre->GetTouchableHandle()->GetCopyNumber(0) : -1;

    // ---- Track identity --------------------------------------------------
    hit->trackID  = track->GetTrackID();
    hit->parentID = track->GetParentID();

    const auto* pd = track->GetDefinition();
    hit->pdg = pd->GetPDGEncoding();
    hit->Z   = pd->GetAtomicNumber();
    hit->A   = pd->GetAtomicMass();

    // ---- Kinematics ------------------------------------------------------
    hit->preKineticEnergy  = pre->GetKineticEnergy();
    hit->postKineticEnergy = post->GetKineticEnergy();
    hit->globalTime        = track->GetGlobalTime();

    // ---- Process names ---------------------------------------------------
    const auto* creator  = track->GetCreatorProcess();
    hit->creatorProcess  = creator ? creator->GetProcessName() : "primary";
    hit->processName     = process ? process->GetProcessName()
                                   : "Transportation";

    // ---- Volume tags -----------------------------------------------------
    // originVolume == 2 means BiCoating throughout EventAction.
    hit->originVolume      = 2;
    hit->firstStepInVolume = firstStep;
    hit->lastStepInVolume  = lastStep;
    hit->interactionStep   = interactionStep;

    // birthVolume: which volume was this track created in?
    hit->birthVolume = BirthVolumeTag(track);

    // ---- Secondaries (needed for fragment identification) ----------------
    if (created) {
        hit->secondaries.reserve(created->size());
        for (const auto* sec : *created) {
            if (!sec || !sec->GetDefinition()) continue;
            SecondaryInfo info;
            info.pdg            = sec->GetDefinition()->GetPDGEncoding();
            info.trackID        = sec->GetTrackID();
            info.parentID       = sec->GetParentID();
            info.Z              = sec->GetDefinition()->GetAtomicNumber();
            info.A              = sec->GetDefinition()->GetAtomicMass();
            info.kineticEnergy  = sec->GetKineticEnergy() / MeV;
            info.name           = sec->GetDefinition()->GetParticleName();
            const auto* cp      = sec->GetCreatorProcess();
            info.creatorProcess = cp ? cp->GetProcessName()
                                     : hit->processName;
            hit->secondaries.push_back(info);
        }
    }

    fHitsCollection->insert(hit);
    return true;
}

// ---------------------------------------------------------------------------
void BiCoatingSD::EndOfEvent(G4HCofThisEvent*) {}
