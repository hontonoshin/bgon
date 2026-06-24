#include "DetectorConstruction.hh"
#include "CrystalSD.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

// No G4Region / G4RegionStore needed — BGO has no thin coating requiring
// a special production cut region. The global default (0.1 mm) is fine.

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction() {}

// ---------------------------------------------------------------------------
void DetectorConstruction::DefineMaterials()
{
    auto* nist = G4NistManager::Instance();

    // BGO — Bi4Ge3O12, density 7.13 g/cm3.
    // Built from NIST elements so the composition is exact.
    auto* Bi = nist->FindOrBuildElement("Bi");
    auto* Ge = nist->FindOrBuildElement("Ge");
    auto* O  = nist->FindOrBuildElement("O");

    auto* BGO = new G4Material("BGO", 7.13*g/cm3, 3, kStateSolid);
    BGO->AddElement(Bi,  4);
    BGO->AddElement(Ge,  3);
    BGO->AddElement(O,  12);

    nist->FindOrBuildMaterial("G4_AIR");
}

// ---------------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();

    auto* air = G4Material::GetMaterial("G4_AIR");
    auto* BGO = G4Material::GetMaterial("BGO");

    // Crystal dimensions — identical to the GAGG simulation.
    // No coating, so pitch = crystal size + inter-crystal gap only.
    const G4double crystalX = 2.5*cm / 2.0;
    const G4double crystalY = 2.5*cm / 2.0;
    const G4double crystalZ = 3.0*cm / 2.0;
    const G4double gap      = 0.1*mm;

    // Pitch: crystal full width + gap (no coating term).
    const G4double pitchXY = 2.0*crystalX + gap;
    const G4double pitchZ  = 2.0*crystalZ + gap;

    const G4double caloHX = NX * pitchXY / 2.0;
    const G4double caloHY = NY * pitchXY / 2.0;
    const G4double caloHZ = NZ * pitchZ  / 2.0;

    // World
    auto* worldBox = new G4Box("World",
                               3.0*caloHX + 10.0*cm,
                               3.0*caloHY + 10.0*cm,
                               3.0*caloHZ + 30.0*cm);
    fWorldLogical = new G4LogicalVolume(worldBox, air, "World");
    auto* worldPhys = new G4PVPlacement(nullptr, {}, fWorldLogical,
                                        "World", nullptr, false, 0, true);

    // Calorimeter envelope (air)
    auto* caloBox = new G4Box("Calorimeter", caloHX, caloHY, caloHZ);
    auto* caloLog = new G4LogicalVolume(caloBox, air, "Calorimeter");
    new G4PVPlacement(nullptr, {}, caloLog, "Calorimeter",
                      fWorldLogical, false, 0, true);

    // BGO crystal logical volume (shared across all placements)
    auto* crystalBox = new G4Box("Crystal", crystalX, crystalY, crystalZ);
    fCrystalLogical  = new G4LogicalVolume(crystalBox, BGO, "Crystal");

    // Place 4×4×2 array of BGO crystals — layout identical to GAGG simulation.
    G4int copyNum = 0;
    for (G4int iz = 0; iz < NZ; ++iz) {
        for (G4int iy = 0; iy < NY; ++iy) {
            for (G4int ix = 0; ix < NX; ++ix) {
                const G4double xPos = (ix - (NX - 1)/2.0) * pitchXY;
                const G4double yPos = (iy - (NY - 1)/2.0) * pitchXY;
                const G4double zPos = (iz - (NZ - 1)/2.0) * pitchZ;
                const G4ThreeVector pos(xPos, yPos, zPos);

                new G4PVPlacement(nullptr, pos, fCrystalLogical, "Crystal",
                                  caloLog, false, copyNum, true);
                ++copyNum;
            }
        }
    }

    // Source plane: just upstream of the calorimeter front face.
    // The calorimeter is centred at z=0, so its front face is at z = -caloHZ.
    fSourceZ = -(caloHZ + 1.0*cm);

    // Visualisation
    fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
    caloLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto* crystalVis = new G4VisAttributes(G4Colour(0.20, 0.45, 0.85, 0.75));
    crystalVis->SetForceSolid(true);
    fCrystalLogical->SetVisAttributes(crystalVis);

    return worldPhys;
}

// ---------------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField()
{
    auto* sdMgr = G4SDManager::GetSDMpointer();

    // Single sensitive detector on BGO crystals.
    // No BiCoatingSD needed — BGO has no coating.
    auto* crystalSD = new CrystalSD("CrystalSD", "CrystalHitsCollection");
    sdMgr->AddNewDetector(crystalSD);
    SetSensitiveDetector("Crystal", crystalSD);
}
