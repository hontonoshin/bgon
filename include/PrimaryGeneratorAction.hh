#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4Types.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

    static G4double GetCurrentEnergy() { return fgCurrentEnergy; }
    static G4double GetCurrentWeight() { return fgCurrentWeight; }

private:
    G4ParticleGun* fParticleGun = nullptr;

    static G4ThreadLocal G4double fgCurrentEnergy;
    static G4ThreadLocal G4double fgCurrentWeight;

    G4double fSourceZ    = 26.0*CLHEP::cm;  

    G4double      SampleEnergy(G4double& weight) const;
    G4ThreeVector SampleHemisphere() const;
};
