#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Array layout — identical to the GAGG simulation so that
    // PrimaryGeneratorAction::GeneratePrimaries() needs no changes.
    static constexpr G4int NX = 4;
    static constexpr G4int NY = 4;
    static constexpr G4int NZ = 2;

    // Z position of the particle gun upstream face (set in Construct()).
    G4double GetSourceZ() const { return fSourceZ; }

private:
    void DefineMaterials();

    G4LogicalVolume* fWorldLogical   = nullptr;
    G4LogicalVolume* fCrystalLogical = nullptr;  // BGO crystal logical volume

    // Source plane: just upstream of the calorimeter front face.
    G4double fSourceZ = 0.0;
};
