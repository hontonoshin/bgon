#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4INCLXXInterfaceStore.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
    SetVerboseLevel(1);

    RegisterPhysics(new G4EmStandardPhysics(1));
    RegisterPhysics(new G4DecayPhysics(1));
    RegisterPhysics(new G4RadioactiveDecayPhysics(1));
    RegisterPhysics(new G4HadronElasticPhysicsHP(1));
    RegisterPhysics(new G4HadronPhysicsINCLXX("inclxx", true, true, false));
    RegisterPhysics(new G4IonINCLXXPhysics(1));

    auto* neutronTrackingCut = new G4NeutronTrackingCut(1);
    neutronTrackingCut->SetTimeLimit(1.0 * ms);
    RegisterPhysics(neutronTrackingCut);
}

void PhysicsList::ConstructParticle()
{
    G4BaryonConstructor().ConstructParticle();
    G4BosonConstructor().ConstructParticle();
    G4LeptonConstructor().ConstructParticle();
    G4MesonConstructor().ConstructParticle();
    G4IonConstructor().ConstructParticle();
    G4ShortLivedConstructor().ConstructParticle();
    G4VModularPhysicsList::ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

    auto* deexParams = G4NuclearLevelData::GetInstance()->GetParameters();
    deexParams->SetDeexChannelsType(G4DeexChannelType(4));

    // Preserve explicit light-cluster production from INCL++ while allowing
    // heavier residual nuclei from Bi and Ge reactions to be de-excited and
    // transported normally.
    G4INCLXXInterfaceStore::GetInstance()->SetMaxClusterMass(8);

    G4cout << "\n[PhysicsList] BGO calorimeter configuration:"
           << "\n  Energy range: 100 MeV to 10 GeV primary neutrons"
           << "\n  EM:           G4EmStandardPhysics"
           << "\n  Elastic:      G4HadronElasticPhysicsHP"
           << "\n  Inelastic:    G4HadronPhysicsINCLXX + de-excitation"
           << "\n  Ions:         G4IonINCLXXPhysics"
           << "\n  Neutron tracking time limit: 1 ms"
           << "\n  Global production cut: 0.1 mm\n"
           << G4endl;
}

void PhysicsList::SetCuts()
{
    // Homogeneous BGO geometry: no thin coating or special low-cut region is
    // present. A 0.1 mm global range cut avoids the excessive low-energy EM
    // track production caused by the former 1 um coating-specific setup.
    constexpr G4double defaultCut = 0.1 * mm;

    SetDefaultCutValue(defaultCut);
    SetCutsWithDefault();
    SetCutValue(defaultCut, "gamma");
    SetCutValue(defaultCut, "e-");
    SetCutValue(defaultCut, "e+");
    SetCutValue(defaultCut, "proton");
}
