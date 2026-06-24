#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

G4ThreadLocal G4double PrimaryGeneratorAction::fgCurrentEnergy = 0.0;
G4ThreadLocal G4double PrimaryGeneratorAction::fgCurrentWeight = 1.0;

namespace SourceModel {

constexpr G4double kMinimumEnergyMeV = 100.0;
constexpr G4double kMaximumEnergyMeV = 10000.0;


enum class RigidityModel { R60_MV, R125_MV, R200_MV };
constexpr RigidityModel kRigidityModel = RigidityModel::R125_MV;

struct SpectrumPoint {
    G4double energyMeV;
    G4double differentialFlux; // relative units; proportional to n cm^-2 MeV^-1
};

// Approximate digitization of the R0 = 60 MV curve in Fig. 1.6.
constexpr std::array<SpectrumPoint, 13> kSpectrumR60 = {{
    {15.0, 6.0e-8}, {20.0, 1.2e-7}, {30.0, 2.8e-7},
    {50.0, 6.5e-7}, {80.0, 1.05e-6}, {120.0, 1.35e-6},
    {180.0, 1.55e-6}, {250.0, 1.35e-6}, {350.0, 9.0e-7},
    {500.0, 4.5e-7}, {700.0, 2.0e-7}, {900.0, 9.0e-8},
    {1050.0, 5.0e-8}
}};

// Approximate digitization of the R0 = 125 MV curve in Fig. 1.6.
constexpr std::array<SpectrumPoint, 19> kSpectrumR125 = {{
    {15.0, 7.0e-8}, {20.0, 1.5e-7}, {30.0, 4.0e-7},
    {50.0, 9.0e-7}, {80.0, 1.65e-6}, {120.0, 2.6e-6},
    {180.0, 3.8e-6}, {250.0, 4.8e-6}, {350.0, 5.3e-6},
    {500.0, 4.4e-6}, {700.0, 3.0e-6}, {1000.0, 1.55e-6},
    {1500.0, 6.5e-7}, {2200.0, 2.8e-7}, {3000.0, 1.2e-7},
    {4500.0,  4.5e-8}, {6000.0,  1.8e-8}, {8000.0,  6.0e-9},
    {10000.0, 2.0e-9}
}};

// Approximate digitization of the R0 = 200 MV curve in Fig. 1.6.
constexpr std::array<SpectrumPoint, 17> kSpectrumR200 = {{
    {15.0, 8.0e-8}, {20.0, 1.8e-7}, {30.0, 5.0e-7},
    {50.0, 1.25e-6}, {80.0, 2.4e-6}, {120.0, 4.0e-6},
    {180.0, 6.3e-6}, {250.0, 8.3e-6}, {350.0, 9.8e-6},
    {500.0, 9.2e-6}, {700.0, 7.2e-6}, {1000.0, 4.8e-6},
    {1500.0, 2.7e-6}, {2200.0, 1.45e-6}, {3200.0, 7.5e-7},
    {4500.0, 3.8e-7}, {6000.0, 1.8e-7}
}};

template <std::size_t N>
G4double LogLogInterpolate(const std::array<SpectrumPoint, N>& points,
                          G4double energyMeV)
{
    if (energyMeV < points.front().energyMeV ||
        energyMeV > points.back().energyMeV) {
        return 0.0;
    }

    auto upper = std::upper_bound(
        points.begin(), points.end(), energyMeV,
        [](G4double value, const SpectrumPoint& point) {
            return value < point.energyMeV;
        });

    if (upper == points.begin()) {
        return points.front().differentialFlux;
    }
    if (upper == points.end()) {
        return points.back().differentialFlux;
    }

    const auto& p1 = *(upper - 1);
    const auto& p2 = *upper;

    const G4double x  = std::log(energyMeV);
    const G4double x1 = std::log(p1.energyMeV);
    const G4double x2 = std::log(p2.energyMeV);
    const G4double y1 = std::log(p1.differentialFlux);
    const G4double y2 = std::log(p2.differentialFlux);

    const G4double fraction = (x - x1) / (x2 - x1);
    return std::exp(y1 + fraction * (y2 - y1));
}

G4double DifferentialSpectrum(G4double energyMeV)
{
    switch (kRigidityModel) {
        case RigidityModel::R60_MV:
            return LogLogInterpolate(kSpectrumR60, energyMeV);
        case RigidityModel::R125_MV:
            return LogLogInterpolate(kSpectrumR125, energyMeV);
        case RigidityModel::R200_MV:
            return LogLogInterpolate(kSpectrumR200, energyMeV);
    }
    return 0.0;
}

G4double MinimumEnergyMeV()
{
    return kMinimumEnergyMeV;
}

G4double MaximumEnergyMeV()
{
    return kMaximumEnergyMeV;
}

struct TabulatedCdf {
    std::vector<G4double> energyMeV;
    std::vector<G4double> cdf;
};

const TabulatedCdf& GetCdf()
{
    // A dense logarithmic grid preserves both the low-energy rise and high-energy tail.
    static const TabulatedCdf table = [] {
        constexpr std::size_t numberOfGridPoints = 4097;

        TabulatedCdf result;
        result.energyMeV.resize(numberOfGridPoints);
        result.cdf.resize(numberOfGridPoints, 0.0);

        const G4double emin = MinimumEnergyMeV();
        const G4double emax = MaximumEnergyMeV();
        const G4double logRatio = std::log(emax / emin);

        for (std::size_t i = 0; i < numberOfGridPoints; ++i) {
            const G4double fraction =
                static_cast<G4double>(i) /
                static_cast<G4double>(numberOfGridPoints - 1);
            result.energyMeV[i] = emin * std::exp(fraction * logRatio);
        }

        for (std::size_t i = 1; i < numberOfGridPoints; ++i) {
            const G4double e1 = result.energyMeV[i - 1];
            const G4double e2 = result.energyMeV[i];
            const G4double f1 = DifferentialSpectrum(e1);
            const G4double f2 = DifferentialSpectrum(e2);

            result.cdf[i] = result.cdf[i - 1]
                + 0.5 * (f1 + f2) * (e2 - e1);
        }

        const G4double normalization = result.cdf.back();
        if (!(normalization > 0.0)) {
            throw std::runtime_error("Solar-neutron source spectrum has zero integral.");
        }

        for (auto& value : result.cdf) {
            value /= normalization;
        }
        result.cdf.back() = 1.0;

        return result;
    }();

    return table;
}

} // namespace SourceModel

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(new G4ParticleGun(1))
{
    auto* neutron = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    if (neutron == nullptr) {
        throw std::runtime_error("Geant4 neutron particle definition was not found.");
    }

    fParticleGun->SetParticleDefinition(neutron);
    fParticleGun->SetParticleEnergy(SourceModel::MinimumEnergyMeV() * MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    G4double eventWeight = 1.0;
    const G4double energy = SampleEnergy(eventWeight);

    fgCurrentEnergy = energy / GeV;
    fgCurrentWeight = eventWeight;

    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0));

    // Preserve the detector footprint and the existing array pitch exactly.
    constexpr G4double crystalHX = 2.5 * cm / 2.0;
    constexpr G4double formerCoatingAllowance = 5.0 * um;
    constexpr G4double gap = 0.1 * mm;
    constexpr G4double pitchXY =
        2.0 * (crystalHX + formerCoatingAllowance) + gap;
    constexpr G4double caloHX = DetectorConstruction::NX * pitchXY / 2.0;
    constexpr G4double caloHY = DetectorConstruction::NY * pitchXY / 2.0;

    // Solar neutrons at the detector are treated as a parallel beam from the Sun.
    // Their impact points are uniform over the complete rectangular front face.
    const G4double x0 = (2.0 * G4UniformRand() - 1.0) * caloHX;
    const G4double y0 = (2.0 * G4UniformRand() - 1.0) * caloHY;

    fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, fSourceZ));
    fParticleGun->GeneratePrimaryVertex(event);

    // The tabulated spectrum is sampled directly, so no importance correction
    // is required. Results are detector response per incident solar neutron.
    if (auto* vertex = event->GetPrimaryVertex()) {
        vertex->SetWeight(eventWeight);
    }
}

G4double PrimaryGeneratorAction::SampleEnergy(G4double& weight) const
{
    const auto& table = SourceModel::GetCdf();
    const G4double randomCdf = G4UniformRand();

    const auto upper = std::lower_bound(
        table.cdf.begin(), table.cdf.end(), randomCdf);

    if (upper == table.cdf.begin()) {
        weight = 1.0;
        return table.energyMeV.front() * MeV;
    }
    if (upper == table.cdf.end()) {
        weight = 1.0;
        return table.energyMeV.back() * MeV;
    }

    const std::size_t i2 = static_cast<std::size_t>(upper - table.cdf.begin());
    const std::size_t i1 = i2 - 1;

    const G4double c1 = table.cdf[i1];
    const G4double c2 = table.cdf[i2];
    const G4double e1 = table.energyMeV[i1];
    const G4double e2 = table.energyMeV[i2];

    const G4double localFraction =
        (c2 > c1) ? (randomCdf - c1) / (c2 - c1) : 0.0;

    // Interpolate in log(E), consistent with the logarithmic source grid.
    const G4double sampledEnergyMeV =
        e1 * std::exp(localFraction * std::log(e2 / e1));

    weight = 1.0;
    return sampledEnergyMeV * MeV;
}

G4ThreeVector PrimaryGeneratorAction::SampleHemisphere() const
{
    // Retained for header/source compatibility; not used for the parallel beam.
    return G4ThreeVector(0.0, 0.0, -1.0);
}
