#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

int main(int argc, char** argv)
{
    // ── Random engine ──────────────────────────────────────────────────────
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
    CLHEP::HepRandom::setTheSeed(12345, 0);

    // ── Run manager ────────────────────────────────────────────────────────
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new PhysicsList());

    runManager->SetUserInitialization(new ActionInitialization());

    //runManager->Initialize();

    // ── UI ─────────────────────────────────────────────────────────────────
    auto* ui = G4UImanager::GetUIpointer();

    if (argc > 1) {
        // Batch mode: execute supplied macro
        G4String cmd = "/control/execute ";
        ui->ApplyCommand(cmd + argv[1]);
    } else {
        // Interactive mode
        G4UIExecutive* uiExec = new G4UIExecutive(argc, argv);
        G4VisManager* visManager = new G4VisExecutive();
        visManager->Initialize();
        ui->ApplyCommand("/control/execute vis.mac");
        uiExec->SessionStart();
        delete visManager;
        delete uiExec;
    }

    delete runManager;
    return 0;
}
