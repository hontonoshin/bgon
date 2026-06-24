#include "CrystalSD.hh"

#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4StepPoint.hh"

namespace {

G4int RegionTag(const G4StepPoint* point)
{
    if (!point) return 0;
    const auto* pv = point->GetPhysicalVolume();
    if (!pv || !pv->GetLogicalVolume()) return 0;

    const auto& name = pv->GetLogicalVolume()->GetName();
    if (name == "Crystal")   return 1;
    if (name == "BiCoating") return 2;
    return 0;
}

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

G4bool IsTransportation(const G4VProcess* process)
{
    return !process || process->GetProcessName() == "Transportation";
}

} // namespace

CrystalSD::CrystalSD(const G4String& name, const G4String& hcName)
    : G4VSensitiveDetector(name)
{
    collectionName.insert(hcName);
}

void CrystalSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new NeutronSimHitsCollection(
        SensitiveDetectorName, collectionName[0]);

    if (fHCID < 0) {
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(
            SensitiveDetectorName + "/" + collectionName[0]);
    }

    hce->AddHitsCollection(fHCID, fHitsCollection);
}

G4bool CrystalSD::ProcessHits(G4Step* step, G4TouchableHistory*)
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

    // Store zero-edep steps only when they carry interaction, boundary, or
    // secondary-production information required by the event analysis.
    if (edep <= 0.0 && !interactionStep && !hasSecondaries
        && !firstStep && !lastStep) {
        return false;
    }

    auto* hit = new NeutronSimHit();

    hit->edep       = edep;
    hit->stepLength = step->GetStepLength();
    hit->prePosition  = pre->GetPosition();
    hit->postPosition = post->GetPosition();
    hit->position     = hit->postPosition;

    hit->copyNo = pre->GetTouchableHandle()
        ? pre->GetTouchableHandle()->GetCopyNumber(0) : -1;

    hit->trackID  = track->GetTrackID();
    hit->parentID = track->GetParentID();
    hit->pdg      = track->GetDefinition()->GetPDGEncoding();
    hit->Z        = track->GetDefinition()->GetAtomicNumber();
    hit->A        = track->GetDefinition()->GetAtomicMass();

    hit->preKineticEnergy  = pre->GetKineticEnergy();
    hit->postKineticEnergy = post->GetKineticEnergy();
    hit->globalTime        = post->GetGlobalTime();

    const auto* creator = track->GetCreatorProcess();
    hit->creatorProcess = creator ? creator->GetProcessName() : "primary";
    hit->processName = process ? process->GetProcessName() : "Transportation";

    hit->originVolume      = RegionTag(pre);
    hit->birthVolume       = BirthVolumeTag(track);
    hit->firstStepInVolume = firstStep;
    hit->lastStepInVolume  = lastStep;
    hit->interactionStep   = interactionStep;

    if (created) {
        hit->secondaries.reserve(created->size());

        for (const auto* sec : *created) {
            if (!sec || !sec->GetDefinition()) continue;

            SecondaryInfo info;
            info.pdg           = sec->GetDefinition()->GetPDGEncoding();
            info.trackID       = sec->GetTrackID();
            info.parentID      = sec->GetParentID();
            info.Z             = sec->GetDefinition()->GetAtomicNumber();
            info.A             = sec->GetDefinition()->GetAtomicMass();
            info.kineticEnergy = sec->GetKineticEnergy() / MeV;
            info.name          = sec->GetDefinition()->GetParticleName();

            const auto* cp = sec->GetCreatorProcess();
            info.creatorProcess = cp ? cp->GetProcessName() : hit->processName;

            hit->secondaries.push_back(info);
        }
    }

    fHitsCollection->insert(hit);
    return true;
}

void CrystalSD::EndOfEvent(G4HCofThisEvent*) {}
