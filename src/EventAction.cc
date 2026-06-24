#include "EventAction.hh"
#include "NeutronSimHit.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace {

constexpr G4double kPromptWindow = 10.0 * ns;
constexpr G4double kFirstInteractionSecondaryMinKEMeV = 0.1;
constexpr G4int kNumberOfCrystals =
    DetectorConstruction::NX * DetectorConstruction::NY * DetectorConstruction::NZ;

std::string Lower(std::string s)
{
    for (char& c : s) {
        c = static_cast<char>(
            std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

bool ContainsNoCase(const std::string& s, const std::string& token)
{
    return Lower(s).find(Lower(token)) != std::string::npos;
}

bool IsInelastic(const std::string& process)
{
    return ContainsNoCase(process, "inelastic");
}

bool IsHadronicSelected(const std::string& process)
{
    return IsInelastic(process)
        || ContainsNoCase(process, "capture")
        || ContainsNoCase(process, "fission");
}

bool IsTransportation(const std::string& process)
{
    return process.empty() || ContainsNoCase(process, "transportation");
}

// 0=other, 1=neutron, 2=proton, 3=light ion A<=4,
// 4=heavy ion A>4, 5=e/gamma, 6=pion/kaon.
int ClassifyPDG(G4int pdg, G4int Z, G4int A)
{
    if (pdg == 2112) return 1;
    if (pdg == 2212) return 2;

    if (pdg == 11 || pdg == -11 || pdg == 22) return 5;

    if (pdg == 211 || pdg == -211 || pdg == 111
        || pdg == 321 || pdg == -321
        || pdg == 130 || pdg == 310) {
        return 6;
    }

    if (std::abs(pdg) > 1000000000 && Z > 0 && A > 0) {
        return (A <= 4) ? 3 : 4;
    }

    return 0;
}

// Integer process code for compact ROOT storage.
// 0 none, 1 elastic, 2 inelastic, 3 capture, 4 fission,
// 5 decay, 6 electromagnetic/ionisation, 99 other.
int ProcessCode(const std::string& process)
{
    if (process.empty()) return 0;
    if (ContainsNoCase(process, "inelastic")) return 2;
    if (ContainsNoCase(process, "capture"))   return 3;
    if (ContainsNoCase(process, "fission"))   return 4;
    if (ContainsNoCase(process, "elastic"))   return 1;
    if (ContainsNoCase(process, "decay"))     return 5;

    if (ContainsNoCase(process, "ioni")
        || ContainsNoCase(process, "compt")
        || ContainsNoCase(process, "phot")
        || ContainsNoCase(process, "conv")
        || ContainsNoCase(process, "brem")
        || ContainsNoCase(process, "msc")) {
        return 6;
    }

    return 99;
}

struct SecondaryCounters
{
    G4int neutrons  = 0;
    G4int protons   = 0;
    G4int gammas    = 0;
    G4int electrons = 0;
    G4int positrons = 0;
    G4int pions     = 0;
    G4int kaons     = 0;
    G4int lightIons = 0;
    G4int heavyIons = 0;
    G4int other     = 0;
};

void CountSecondary(const SecondaryInfo& secondary, SecondaryCounters& counts)
{
    const G4int pdg = secondary.pdg;

    if (pdg == 2112) { ++counts.neutrons; return; }
    if (pdg == 2212) { ++counts.protons;  return; }
    if (pdg == 22)   { ++counts.gammas;   return; }
    if (pdg == 11)   { ++counts.electrons; return; }
    if (pdg == -11)  { ++counts.positrons; return; }

    if (pdg == 211 || pdg == -211 || pdg == 111) {
        ++counts.pions;
        return;
    }

    if (pdg == 321 || pdg == -321 || pdg == 130 || pdg == 310) {
        ++counts.kaons;
        return;
    }

    if (std::abs(pdg) > 1000000000 && secondary.Z > 0 && secondary.A > 0) {
        if (secondary.A <= 4) ++counts.lightIons;
        else                  ++counts.heavyIons;
        return;
    }

    ++counts.other;
}

G4double CrystalFrontZ()
{
    constexpr G4double crystalHalfZ = 1.5 * cm;
    constexpr G4double gapZ = 0.1 * mm;
    constexpr G4double pitchZ = 2.0 * crystalHalfZ + gapZ;

    return 0.5 * (DetectorConstruction::NZ - 1) * pitchZ + crystalHalfZ;
}

G4double CalorimeterActiveSpanZ()
{
    constexpr G4double crystalFullZ = 3.0 * cm;
    constexpr G4double gapZ = 0.1 * mm;

    return DetectorConstruction::NZ * crystalFullZ
         + (DetectorConstruction::NZ - 1) * gapZ;
}

G4double DepthFromFront(const G4ThreeVector& position)
{
    return CrystalFrontZ() - position.z();
}

} // namespace

EventAction::EventAction() : G4UserEventAction() {}

void EventAction::BeginOfEventAction(const G4Event*) {}

void EventAction::EndOfEventAction(const G4Event* event)
{
    if (!event) return;

    G4double primaryEnergyGeV = 0.0;
    G4double sourceX = 0.0;
    G4double sourceY = 0.0;

    if (const auto* vertex = event->GetPrimaryVertex()) {
        sourceX = vertex->GetX0();
        sourceY = vertex->GetY0();

        if (const auto* primary = vertex->GetPrimary()) {
            primaryEnergyGeV = primary->GetKineticEnergy() / GeV;
        }
    }

    const G4double eventWeight = PrimaryGeneratorAction::GetCurrentWeight();
    const G4int eventID = event->GetEventID();

    NeutronSimHitsCollection* crystalHC = nullptr;
    if (auto* hce = event->GetHCofThisEvent()) {
        if (fCrystalHCID < 0) {
            fCrystalHCID = G4SDManager::GetSDMpointer()->GetCollectionID(
                "CrystalSD/CrystalHitsCollection");
        }

        if (fCrystalHCID >= 0) {
            crystalHC = static_cast<NeutronSimHitsCollection*>(
                hce->GetHC(fCrystalHCID));
        }
    }

    auto* analysis = G4AnalysisManager::Instance();

    G4double edepTotal    = 0.0;
    G4double edepNP       = 0.0;
    G4double edepLightIon = 0.0;
    G4double edepHeavyIon = 0.0;
    G4double edepEM       = 0.0;
    G4double edepPionKaon = 0.0;
    G4double edepOther    = 0.0;

    G4double totalTrackLength = 0.0;
    G4double maxStepEdep      = 0.0;

    G4double weightedDepthSum   = 0.0;
    G4double weightedDepth2Sum  = 0.0;
    G4double weightedRadiusSum  = 0.0;
    G4double weightedRadius2Sum = 0.0;

    G4double firstDepositTime = std::numeric_limits<G4double>::max();
    G4double lastDepositTime  = -1.0;

    G4double edepFront  = 0.0;
    G4double edepMiddle = 0.0;
    G4double edepRear   = 0.0;

    G4double edepFirstLayer  = 0.0;
    G4double edepSecondLayer = 0.0;

    bool anyHadronicInteraction   = false;
    bool primaryNeutronHadronic   = false;
    bool primaryNeutronInelastic  = false;
    G4int nPrimaryNeutronInelastic = 0;

    const NeutronSimHit* firstPrimaryHadronicHit = nullptr;
    G4double firstPrimaryHadronicTime = std::numeric_limits<G4double>::max();

    SecondaryCounters secondaryCounts;
    G4int nSecondaries = 0;
    G4double sumSecondaryKE = 0.0;

    std::unordered_set<G4int> uniqueTrackIDs;
    std::unordered_map<G4int, const NeutronSimHit*> lastHitByTrack;
    std::array<G4double, kNumberOfCrystals> crystalEdep{};

    const G4double activeSpanZ = CalorimeterActiveSpanZ();
    const G4double oneThirdZ = activeSpanZ / 3.0;

    const G4int nHits = crystalHC ? crystalHC->entries() : 0;

    for (G4int i = 0; i < nHits; ++i) {
        const auto* hit = (*crystalHC)[i];
        if (!hit) continue;

        uniqueTrackIDs.insert(hit->trackID);
        totalTrackLength += hit->stepLength;
        maxStepEdep = std::max(maxStepEdep, hit->edep);

        auto lastIt = lastHitByTrack.find(hit->trackID);
        if (lastIt == lastHitByTrack.end()
            || !lastIt->second
            || hit->globalTime >= lastIt->second->globalTime) {
            lastHitByTrack[hit->trackID] = hit;
        }

        const bool selectedHadronic = IsHadronicSelected(hit->processName);
        if (selectedHadronic) anyHadronicInteraction = true;

        if (hit->trackID == 1 && hit->pdg == 2112) {
            if (selectedHadronic) {
                primaryNeutronHadronic = true;
                if (hit->globalTime < firstPrimaryHadronicTime) {
                    firstPrimaryHadronicTime = hit->globalTime;
                    firstPrimaryHadronicHit = hit;
                }
            }

            if (IsInelastic(hit->processName)) {
                primaryNeutronInelastic = true;
                ++nPrimaryNeutronInelastic;
            }
        }

        for (const auto& secondary : hit->secondaries) {
            ++nSecondaries;
            sumSecondaryKE += secondary.kineticEnergy * MeV;
            CountSecondary(secondary, secondaryCounts);
        }

        edepTotal += hit->edep;

        switch (ClassifyPDG(hit->pdg, hit->Z, hit->A)) {
            case 1:
            case 2: edepNP       += hit->edep; break;
            case 3: edepLightIon += hit->edep; break;
            case 4: edepHeavyIon += hit->edep; break;
            case 5: edepEM       += hit->edep; break;
            case 6: edepPionKaon += hit->edep; break;
            default: edepOther   += hit->edep; break;
        }

        if (hit->copyNo >= 0 && hit->copyNo < kNumberOfCrystals) {
            crystalEdep[static_cast<std::size_t>(hit->copyNo)] += hit->edep;

            const G4int layer = hit->copyNo /
                (DetectorConstruction::NX * DetectorConstruction::NY);
            if (layer == 0) edepFirstLayer += hit->edep;
            else if (layer == 1) edepSecondLayer += hit->edep;
        }

        if (hit->edep <= 0.0) continue;

        firstDepositTime = std::min(firstDepositTime, hit->globalTime);
        lastDepositTime  = std::max(lastDepositTime,  hit->globalTime);

        const G4ThreeVector midpoint =
            0.5 * (hit->prePosition + hit->postPosition);
        const G4double depth = DepthFromFront(midpoint);
        const G4double dx = midpoint.x() - sourceX;
        const G4double dy = midpoint.y() - sourceY;
        const G4double radius = std::sqrt(dx*dx + dy*dy);

        weightedDepthSum   += hit->edep * depth;
        weightedDepth2Sum  += hit->edep * depth * depth;
        weightedRadiusSum  += hit->edep * radius;
        weightedRadius2Sum += hit->edep * radius * radius;

        if (depth < oneThirdZ) edepFront += hit->edep;
        else if (depth < 2.0 * oneThirdZ) edepMiddle += hit->edep;
        else edepRear += hit->edep;

        // Energy-weighted shower profiles. Histogram content is deposited
        // energy in MeV, already multiplied by the event weight.
        if (depth >= 0.0 && depth <= activeSpanZ) {
            analysis->FillH2(4, primaryEnergyGeV, depth / mm,
                             (hit->edep / MeV) * eventWeight);
        }
        analysis->FillH2(5, primaryEnergyGeV, radius / mm,
                         (hit->edep / MeV) * eventWeight);
    }

    G4double promptEdep  = 0.0;
    G4double delayedEdep = 0.0;

    if (crystalHC && firstDepositTime < std::numeric_limits<G4double>::max()) {
        for (G4int i = 0; i < nHits; ++i) {
            const auto* hit = (*crystalHC)[i];
            if (!hit || hit->edep <= 0.0) continue;

            if (hit->globalTime - firstDepositTime <= kPromptWindow)
                promptEdep += hit->edep;
            else
                delayedEdep += hit->edep;
        }
    }

    G4double depthCentroid  = -1.0;
    G4double depthRMS       = -1.0;
    G4double radialCentroid = -1.0;
    G4double radialRMS      = -1.0;

    if (edepTotal > 0.0) {
        depthCentroid  = weightedDepthSum / edepTotal;
        radialCentroid = weightedRadiusSum / edepTotal;

        const G4double depthVariance = std::max(
            0.0, weightedDepth2Sum / edepTotal - depthCentroid * depthCentroid);
        const G4double radialVariance = std::max(
            0.0, weightedRadius2Sum / edepTotal - radialCentroid * radialCentroid);

        depthRMS  = std::sqrt(depthVariance);
        radialRMS = std::sqrt(radialVariance);
    }

    G4double leakageTotal    = 0.0;
    G4double leakageNP       = 0.0;
    G4double leakageEM       = 0.0;
    G4double leakageLightIon = 0.0;
    G4double leakageHeavyIon = 0.0;
    G4double leakagePionKaon = 0.0;

    bool primaryNeutronEscaped = false;
    G4double primaryNeutronExitKE = 0.0;

    for (const auto& item : lastHitByTrack) {
        const auto* hit = item.second;
        if (!hit || !hit->lastStepInVolume
            || !IsTransportation(hit->processName)
            || hit->postKineticEnergy <= 0.0) {
            continue;
        }

        const G4double escapingKE = hit->postKineticEnergy;
        leakageTotal += escapingKE;

        switch (ClassifyPDG(hit->pdg, hit->Z, hit->A)) {
            case 1:
            case 2: leakageNP       += escapingKE; break;
            case 3: leakageLightIon += escapingKE; break;
            case 4: leakageHeavyIon += escapingKE; break;
            case 5: leakageEM       += escapingKE; break;
            case 6: leakagePionKaon += escapingKE; break;
            default: break;
        }

        if (hit->trackID == 1 && hit->pdg == 2112) {
            primaryNeutronEscaped = true;
            primaryNeutronExitKE = escapingKE;
        }
    }

    G4int nHitCrystals = 0;
    G4int maxCrystalID = -1;
    G4double maxCrystalEdep = 0.0;

    for (G4int crystalID = 0; crystalID < kNumberOfCrystals; ++crystalID) {
        const G4double e = crystalEdep[static_cast<std::size_t>(crystalID)];
        if (e <= 0.0) continue;

        ++nHitCrystals;
        if (e > maxCrystalEdep) {
            maxCrystalEdep = e;
            maxCrystalID = crystalID;
        }

        analysis->FillH2(6, primaryEnergyGeV,
                         static_cast<G4double>(crystalID),
                         (e / MeV) * eventWeight);
    }

    const G4double responseFraction =
        (primaryEnergyGeV > 0.0)
        ? (edepTotal / (primaryEnergyGeV * GeV)) : 0.0;

    G4int firstInteractionProcessCode = 0;
    G4double firstInteractionDepth = -1.0;
    G4double firstInteractionRadius = -1.0;
    G4double firstInteractionTime = -1.0;
    G4double firstInteractionPreKE = -1.0;
    G4int firstInteractionMultiplicity = 0;
    G4double firstInteractionSecondaryKE = 0.0;

    if (firstPrimaryHadronicHit) {
        firstInteractionProcessCode = ProcessCode(
            firstPrimaryHadronicHit->processName);
        firstInteractionDepth = DepthFromFront(
            firstPrimaryHadronicHit->postPosition);

        const G4double dx = firstPrimaryHadronicHit->postPosition.x() - sourceX;
        const G4double dy = firstPrimaryHadronicHit->postPosition.y() - sourceY;
        firstInteractionRadius = std::sqrt(dx*dx + dy*dy);
        firstInteractionTime = firstPrimaryHadronicHit->globalTime;
        firstInteractionPreKE = firstPrimaryHadronicHit->preKineticEnergy;
        firstInteractionMultiplicity = static_cast<G4int>(
            firstPrimaryHadronicHit->secondaries.size());

        for (const auto& secondary : firstPrimaryHadronicHit->secondaries) {
            firstInteractionSecondaryKE += secondary.kineticEnergy * MeV;
        }
    }

    const G4double componentSum = edepNP + edepLightIon + edepHeavyIon
                                + edepEM + edepPionKaon + edepOther;
    const G4double closureCheck = edepTotal - componentSum;

    // Existing histograms: IDs 0-8 remain unchanged.
    analysis->FillH1(0, edepTotal / MeV, eventWeight);
    analysis->FillH1(1, primaryEnergyGeV, eventWeight);
    analysis->FillH1(2, edepNP / MeV, eventWeight);
    analysis->FillH1(3, edepHeavyIon / MeV, eventWeight);
    analysis->FillH1(4, edepEM / MeV, eventWeight);
    analysis->FillH1(5, edepPionKaon / MeV, eventWeight);
    analysis->FillH1(6, edepLightIon / MeV, eventWeight);
    analysis->FillH1(8, primaryEnergyGeV, eventWeight);
    if (anyHadronicInteraction)
        analysis->FillH1(7, primaryEnergyGeV, eventWeight);

    // New histograms.
    if (primaryNeutronHadronic)
        analysis->FillH1(9, primaryEnergyGeV, eventWeight);
    if (primaryNeutronInelastic)
        analysis->FillH1(10, primaryEnergyGeV, eventWeight);
    if (responseFraction > 0.0)
        analysis->FillH1(11, responseFraction, eventWeight);
    if (firstInteractionDepth >= 0.0)
        analysis->FillH1(12, firstInteractionDepth / mm, eventWeight);

    analysis->FillH1(13, static_cast<G4double>(nSecondaries), eventWeight);
    analysis->FillH1(14, static_cast<G4double>(uniqueTrackIDs.size()), eventWeight);
    analysis->FillH1(15, totalTrackLength / mm, eventWeight);

    if (depthCentroid >= 0.0)
        analysis->FillH1(16, depthCentroid / mm, eventWeight);
    if (radialRMS >= 0.0)
        analysis->FillH1(17, radialRMS / mm, eventWeight);
    if (promptEdep > 0.0)
        analysis->FillH1(18, promptEdep / MeV, eventWeight);
    if (delayedEdep > 0.0)
        analysis->FillH1(19, delayedEdep / MeV, eventWeight);
    if (leakageTotal > 0.0)
        analysis->FillH1(20, leakageTotal / MeV, eventWeight);
    if (primaryNeutronExitKE > 0.0)
        analysis->FillH1(24, primaryNeutronExitKE / MeV, eventWeight);

    analysis->FillH2(0, primaryEnergyGeV, edepTotal / MeV, eventWeight);
    if (responseFraction > 0.0)
        analysis->FillH2(1, primaryEnergyGeV, responseFraction, eventWeight);
    if (firstInteractionDepth >= 0.0)
        analysis->FillH2(2, primaryEnergyGeV,
                         firstInteractionDepth / mm, eventWeight);
    analysis->FillH2(3, primaryEnergyGeV,
                     static_cast<G4double>(nSecondaries), eventWeight);
    if (leakageTotal > 0.0)
        analysis->FillH2(8, primaryEnergyGeV,
                         leakageTotal / MeV, eventWeight);
    if (depthCentroid >= 0.0)
        analysis->FillH2(9, primaryEnergyGeV,
                         depthCentroid / mm, eventWeight);
    if (radialRMS >= 0.0)
        analysis->FillH2(10, primaryEnergyGeV,
                         radialRMS / mm, eventWeight);

    // Ntuple 0: preserve columns 0-23 exactly, then append new observables.
    analysis->FillNtupleIColumn(0,  0, eventID);
    analysis->FillNtupleDColumn(0,  1, primaryEnergyGeV);
    analysis->FillNtupleDColumn(0,  2, edepTotal / MeV);
    analysis->FillNtupleDColumn(0,  3, 0.0);
    analysis->FillNtupleDColumn(0,  4, 0.0);
    analysis->FillNtupleDColumn(0,  5, 0.0);
    analysis->FillNtupleDColumn(0,  6, edepNP / MeV);
    analysis->FillNtupleDColumn(0,  7, edepHeavyIon / MeV);
    analysis->FillNtupleDColumn(0,  8, edepEM / MeV);
    analysis->FillNtupleDColumn(0,  9, edepPionKaon / MeV);
    analysis->FillNtupleDColumn(0, 10, edepLightIon / MeV);
    analysis->FillNtupleDColumn(0, 11, closureCheck / MeV);
    analysis->FillNtupleIColumn(0, 12, anyHadronicInteraction ? 1 : 0);
    analysis->FillNtupleIColumn(0, 13, 0);
    analysis->FillNtupleIColumn(0, 14, 0);
    analysis->FillNtupleIColumn(0, 15, 0);
    analysis->FillNtupleIColumn(0, 16, 0);
    analysis->FillNtupleIColumn(0, 17, 0);
    analysis->FillNtupleIColumn(0, 18, 0);
    analysis->FillNtupleIColumn(0, 19, 0);
    analysis->FillNtupleIColumn(0, 20, 0);
    analysis->FillNtupleIColumn(0, 21, 0);
    analysis->FillNtupleIColumn(0, 22, 0);
    analysis->FillNtupleDColumn(0, 23, eventWeight);

    analysis->FillNtupleDColumn(0, 24, edepOther / MeV);
    analysis->FillNtupleDColumn(0, 25, responseFraction);
    analysis->FillNtupleDColumn(0, 26, sourceX / mm);
    analysis->FillNtupleDColumn(0, 27, sourceY / mm);
    analysis->FillNtupleIColumn(0, 28, nHits);
    analysis->FillNtupleIColumn(0, 29, static_cast<G4int>(uniqueTrackIDs.size()));
    analysis->FillNtupleIColumn(0, 30, nSecondaries);
    analysis->FillNtupleIColumn(0, 31, secondaryCounts.neutrons);
    analysis->FillNtupleIColumn(0, 32, secondaryCounts.protons);
    analysis->FillNtupleIColumn(0, 33, secondaryCounts.gammas);
    analysis->FillNtupleIColumn(0, 34, secondaryCounts.electrons);
    analysis->FillNtupleIColumn(0, 35, secondaryCounts.positrons);
    analysis->FillNtupleIColumn(0, 36, secondaryCounts.pions);
    analysis->FillNtupleIColumn(0, 37, secondaryCounts.kaons);
    analysis->FillNtupleIColumn(0, 38, secondaryCounts.lightIons);
    analysis->FillNtupleIColumn(0, 39, secondaryCounts.heavyIons);
    analysis->FillNtupleIColumn(0, 40, secondaryCounts.other);
    analysis->FillNtupleDColumn(0, 41, sumSecondaryKE / MeV);
    analysis->FillNtupleIColumn(0, 42, primaryNeutronHadronic ? 1 : 0);
    analysis->FillNtupleIColumn(0, 43, primaryNeutronInelastic ? 1 : 0);
    analysis->FillNtupleIColumn(0, 44, nPrimaryNeutronInelastic);
    analysis->FillNtupleIColumn(0, 45, firstInteractionProcessCode);
    analysis->FillNtupleDColumn(0, 46, firstInteractionDepth >= 0.0
        ? firstInteractionDepth / mm : -1.0);
    analysis->FillNtupleDColumn(0, 47, firstInteractionRadius >= 0.0
        ? firstInteractionRadius / mm : -1.0);
    analysis->FillNtupleDColumn(0, 48, firstInteractionTime >= 0.0
        ? firstInteractionTime / ns : -1.0);
    analysis->FillNtupleDColumn(0, 49, firstInteractionPreKE >= 0.0
        ? firstInteractionPreKE / MeV : -1.0);
    analysis->FillNtupleIColumn(0, 50, firstInteractionMultiplicity);
    analysis->FillNtupleDColumn(0, 51, firstInteractionSecondaryKE / MeV);
    analysis->FillNtupleDColumn(0, 52, totalTrackLength / mm);
    analysis->FillNtupleDColumn(0, 53, maxStepEdep / MeV);
    analysis->FillNtupleDColumn(0, 54, depthCentroid >= 0.0
        ? depthCentroid / mm : -1.0);
    analysis->FillNtupleDColumn(0, 55, depthRMS >= 0.0
        ? depthRMS / mm : -1.0);
    analysis->FillNtupleDColumn(0, 56, radialCentroid >= 0.0
        ? radialCentroid / mm : -1.0);
    analysis->FillNtupleDColumn(0, 57, radialRMS >= 0.0
        ? radialRMS / mm : -1.0);
    analysis->FillNtupleDColumn(0, 58, edepFront / MeV);
    analysis->FillNtupleDColumn(0, 59, edepMiddle / MeV);
    analysis->FillNtupleDColumn(0, 60, edepRear / MeV);

    const bool hasDepositTime =
        firstDepositTime < std::numeric_limits<G4double>::max();
    analysis->FillNtupleDColumn(0, 61,
        hasDepositTime ? firstDepositTime / ns : -1.0);
    analysis->FillNtupleDColumn(0, 62,
        hasDepositTime ? lastDepositTime / ns : -1.0);
    analysis->FillNtupleDColumn(0, 63,
        hasDepositTime ? (lastDepositTime - firstDepositTime) / ns : -1.0);
    analysis->FillNtupleDColumn(0, 64, promptEdep / MeV);
    analysis->FillNtupleDColumn(0, 65, delayedEdep / MeV);
    analysis->FillNtupleDColumn(0, 66, leakageTotal / MeV);
    analysis->FillNtupleDColumn(0, 67, leakageNP / MeV);
    analysis->FillNtupleDColumn(0, 68, leakageEM / MeV);
    analysis->FillNtupleDColumn(0, 69, leakageLightIon / MeV);
    analysis->FillNtupleDColumn(0, 70, leakageHeavyIon / MeV);
    analysis->FillNtupleDColumn(0, 71, leakagePionKaon / MeV);
    analysis->FillNtupleIColumn(0, 72, primaryNeutronEscaped ? 1 : 0);
    analysis->FillNtupleDColumn(0, 73, primaryNeutronExitKE / MeV);
    analysis->FillNtupleIColumn(0, 74, maxCrystalID);
    analysis->FillNtupleDColumn(0, 75, maxCrystalEdep / MeV);
    analysis->FillNtupleIColumn(0, 76, nHitCrystals);
    analysis->FillNtupleDColumn(0, 77, edepFirstLayer / MeV);
    analysis->FillNtupleDColumn(0, 78, edepSecondLayer / MeV);
    analysis->AddNtupleRow(0);

    // Ntuple 1: secondaries from the first selected primary-neutron hadronic
    // interaction. Event-level multiplicities remain unfiltered in ntuple 0;
    // only this row-wise diagnostic ntuple uses a 0.1 MeV size-control cut.
    if (firstPrimaryHadronicHit) {
        G4int secondaryIndex = 0;

        for (const auto& secondary : firstPrimaryHadronicHit->secondaries) {
            if (secondary.kineticEnergy < kFirstInteractionSecondaryMinKEMeV) {
                ++secondaryIndex;
                continue;
            }

            const G4int particleClass = ClassifyPDG(
                secondary.pdg, secondary.Z, secondary.A);
            const G4int isFragment =
                (std::abs(secondary.pdg) > 1000000000
                 && secondary.Z > 0 && secondary.A > 1) ? 1 : 0;

            analysis->FillNtupleIColumn(1,  0, eventID);
            analysis->FillNtupleDColumn(1,  1, primaryEnergyGeV);
            analysis->FillNtupleDColumn(1,  2, eventWeight);
            analysis->FillNtupleIColumn(1,  3, firstInteractionProcessCode);
            analysis->FillNtupleIColumn(1,  4, secondaryIndex);
            analysis->FillNtupleIColumn(1,  5, secondary.pdg);
            analysis->FillNtupleIColumn(1,  6, secondary.Z);
            analysis->FillNtupleIColumn(1,  7, secondary.A);
            analysis->FillNtupleDColumn(1,  8, secondary.kineticEnergy);
            analysis->FillNtupleDColumn(1,  9,
                firstPrimaryHadronicHit->postPosition.x() / mm);
            analysis->FillNtupleDColumn(1, 10,
                firstPrimaryHadronicHit->postPosition.y() / mm);
            analysis->FillNtupleDColumn(1, 11,
                firstPrimaryHadronicHit->postPosition.z() / mm);
            analysis->FillNtupleDColumn(1, 12,
                firstInteractionDepth / mm);
            analysis->FillNtupleDColumn(1, 13,
                firstPrimaryHadronicHit->globalTime / ns);
            analysis->FillNtupleIColumn(1, 14,
                firstPrimaryHadronicHit->trackID);
            analysis->FillNtupleIColumn(1, 15, secondary.trackID);
            analysis->FillNtupleIColumn(1, 16, particleClass);
            analysis->FillNtupleIColumn(1, 17, isFragment);
            analysis->AddNtupleRow(1);

            if (secondary.kineticEnergy > 0.0) {
                analysis->FillH1(21, secondary.kineticEnergy, eventWeight);
            }

            if (isFragment) {
                analysis->FillH1(22, static_cast<G4double>(secondary.A), eventWeight);
                analysis->FillH1(23, static_cast<G4double>(secondary.Z), eventWeight);
                analysis->FillH2(7,
                    static_cast<G4double>(secondary.Z),
                    static_cast<G4double>(secondary.A),
                    eventWeight);
            }

            ++secondaryIndex;
        }
    }
}
