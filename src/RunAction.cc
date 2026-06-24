#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4ios.hh"

RunAction::RunAction() : G4UserRunAction()
{
    auto* analysis = G4AnalysisManager::Instance();
    analysis->SetDefaultFileType("root");
    analysis->SetVerboseLevel(1);
    analysis->SetNtupleMerging(true);

    // ------------------------------------------------------------------
    // Original histograms: IDs 0-8 are unchanged so the existing analysis
    // macro can continue to read the same objects.
    // ------------------------------------------------------------------
    analysis->CreateH1(
        "EdepBGO_total",
        "Total BGO energy deposition per event;"
        "E_{dep} (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H0

    analysis->CreateH1(
        "PrimaryNeutronEnergy",
        "Primary neutron energy;E_{n} (GeV);Weighted primaries / log bin",
        300, 0.1, 10.0, "none", "none", "log");                    // H1

    analysis->CreateH1(
        "EdepBGO_np",
        "BGO Edep from n+p;E_{dep} (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H2

    analysis->CreateH1(
        "EdepBGO_heavyIon",
        "BGO Edep from heavy ions (A>4);"
        "E_{dep} (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H3

    analysis->CreateH1(
        "EdepBGO_em",
        "BGO Edep from e^{#pm}/#gamma;"
        "E_{dep} (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H4

    analysis->CreateH1(
        "EdepBGO_piK",
        "BGO Edep from pions and kaons;"
        "E_{dep} (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H5

    analysis->CreateH1(
        "EdepBGO_lightIon",
        "BGO Edep from light ions (A#leq4, excluding proton);"
        "E_{dep} (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H6

    analysis->CreateH1(
        "HadronicInBGO_num",
        "Events with any selected hadronic process in BGO;"
        "E_{n} (GeV);Weighted events",
        50, 0.1, 10.0, "none", "none", "log");                      // H7

    analysis->CreateH1(
        "HadronicInBGO_den",
        "All incident-neutron events;E_{n} (GeV);Weighted events",
        50, 0.1, 10.0, "none", "none", "log");                      // H8

    // ------------------------------------------------------------------
    // New one-dimensional histograms.
    // ------------------------------------------------------------------
    analysis->CreateH1(
        "PrimaryNeutronHadronic_num",
        "Primary neutron with inelastic/capture/fission process in BGO;"
        "E_{n} (GeV);Weighted events",
        50, 0.1, 10.0, "none", "none", "log");                      // H9

    analysis->CreateH1(
        "PrimaryNeutronInelastic_num",
        "Primary neutron with inelastic process in BGO;"
        "E_{n} (GeV);Weighted events",
        50, 0.1, 10.0, "none", "none", "log");                      // H10

    analysis->CreateH1(
        "BGO_ResponseFraction",
        "BGO deposited-energy fraction;E_{dep}/E_{n};Weighted events",
        400, 1.0e-6, 2.0, "none", "none", "log");                   // H11

    analysis->CreateH1(
        "FirstInteractionDepth",
        "First selected primary-neutron hadronic interaction depth;"
        "Depth from front BGO face (mm);Weighted events",
        350, 0.0, 70.0, "none", "none", "linear");                  // H12

    analysis->CreateH1(
        "SecondaryMultiplicity",
        "All secondaries produced in BGO per event;"
        "N_{secondary};Weighted events",
        500, -0.5, 499.5, "none", "none", "linear");                // H13

    analysis->CreateH1(
        "UniqueTrackMultiplicity",
        "Unique tracks producing stored BGO hits per event;"
        "N_{tracks};Weighted events",
        1000, -0.5, 999.5, "none", "none", "linear");               // H14

    analysis->CreateH1(
        "TotalTrackLengthBGO",
        "Total stored track length in BGO per event;"
        "Track length (mm);Weighted events",
        500, 0.0, 1.0e5, "none", "none", "linear");                 // H15

    analysis->CreateH1(
        "ShowerDepthCentroid",
        "Energy-weighted shower depth centroid;"
        "Depth from front BGO face (mm);Weighted events",
        350, 0.0, 70.0, "none", "none", "linear");                  // H16

    analysis->CreateH1(
        "ShowerRadialRMS",
        "Energy-weighted transverse shower RMS;"
        "Radial RMS about incident axis (mm);Weighted events",
        400, 0.0, 100.0, "none", "none", "linear");                 // H17

    analysis->CreateH1(
        "PromptEdep10ns",
        "BGO energy deposited within 10 ns of first deposit;"
        "Prompt E_{dep} (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H18

    analysis->CreateH1(
        "DelayedEdepAfter10ns",
        "BGO energy deposited later than 10 ns after first deposit;"
        "Delayed E_{dep} (MeV);Weighted events / log bin",
        500, 1.0e-6, 1.0e4, "none", "none", "log");                 // H19

    analysis->CreateH1(
        "EscapingKineticEnergy",
        "Approximate kinetic-energy leakage from final BGO exits;"
        "Escaping kinetic energy (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H20

    analysis->CreateH1(
        "FirstInteractionSecondaryKE",
        "Secondary kinetic energy at first primary-neutron hadronic interaction;"
        "Secondary kinetic energy (MeV);Weighted secondaries / log bin",
        500, 1.0e-4, 1.0e4, "none", "none", "log");                 // H21

    analysis->CreateH1(
        "FirstInteractionFragmentA",
        "Fragment mass number from first primary-neutron hadronic interaction;"
        "A;Weighted fragments",
        220, 0.5, 220.5, "none", "none", "linear");                 // H22

    analysis->CreateH1(
        "FirstInteractionFragmentZ",
        "Fragment atomic number from first primary-neutron hadronic interaction;"
        "Z;Weighted fragments",
        90, 0.5, 90.5, "none", "none", "linear");                   // H23

    analysis->CreateH1(
        "PrimaryNeutronExitKE",
        "Kinetic energy of primary neutrons leaving their final BGO visit;"
        "Exit kinetic energy (MeV);Weighted events / log bin",
        500, 1.0e-3, 1.0e4, "none", "none", "log");                 // H24

    // ------------------------------------------------------------------
    // Two-dimensional histograms.
    // ------------------------------------------------------------------
    analysis->CreateH2(
        "PrimaryE_vs_EdepBGO",
        "Primary neutron energy vs total BGO Edep;"
        "Primary neutron energy (GeV);Total BGO E_{dep} (MeV)",
        100, 0.1, 10.0,
        100, 1.0e-3, 1.0e4,
        "none", "none", "none", "none", "log", "log");            // H2-0

    analysis->CreateH2(
        "PrimaryE_vs_ResponseFraction",
        "Primary neutron energy vs BGO response fraction;"
        "Primary neutron energy (GeV);E_{dep}/E_{n}",
        100, 0.1, 10.0,
        120, 1.0e-6, 2.0,
        "none", "none", "none", "none", "log", "log");            // H2-1

    analysis->CreateH2(
        "PrimaryE_vs_FirstInteractionDepth",
        "Primary energy vs first primary-neutron hadronic depth;"
        "Primary neutron energy (GeV);Depth from front face (mm)",
        100, 0.1, 10.0,
        140, 0.0, 70.0,
        "none", "none", "none", "none", "log", "linear");         // H2-2

    analysis->CreateH2(
        "PrimaryE_vs_SecondaryMultiplicity",
        "Primary energy vs total secondary multiplicity;"
        "Primary neutron energy (GeV);N_{secondary}",
        100, 0.1, 10.0,
        250, -0.5, 499.5,
        "none", "none", "none", "none", "log", "linear");         // H2-3

    analysis->CreateH2(
        "PrimaryE_vs_DepthProfile",
        "Energy-weighted longitudinal BGO shower profile;"
        "Primary neutron energy (GeV);Depth from front face (mm)",
        100, 0.1, 10.0,
        140, 0.0, 70.0,
        "none", "none", "none", "none", "log", "linear");         // H2-4

    analysis->CreateH2(
        "PrimaryE_vs_RadialProfile",
        "Energy-weighted transverse BGO shower profile;"
        "Primary neutron energy (GeV);Radius from incident axis (mm)",
        100, 0.1, 10.0,
        200, 0.0, 100.0,
        "none", "none", "none", "none", "log", "linear");         // H2-5

    analysis->CreateH2(
        "PrimaryE_vs_CrystalID",
        "Energy deposited by crystal and primary energy;"
        "Primary neutron energy (GeV);Crystal copy number",
        100, 0.1, 10.0,
        32, -0.5, 31.5,
        "none", "none", "none", "none", "log", "linear");         // H2-6

    analysis->CreateH2(
        "FirstInteractionFragmentZA",
        "Fragment map at first primary-neutron hadronic interaction;Z;A",
        90, 0.5, 90.5,
        220, 0.5, 220.5,
        "none", "none", "none", "none", "linear", "linear");      // H2-7

    analysis->CreateH2(
        "PrimaryE_vs_LeakageKE",
        "Primary neutron energy vs approximate kinetic-energy leakage;"
        "Primary neutron energy (GeV);Escaping kinetic energy (MeV)",
        100, 0.1, 10.0,
        120, 1.0e-3, 1.0e4,
        "none", "none", "none", "none", "log", "log");            // H2-8

    analysis->CreateH2(
        "PrimaryE_vs_DepthCentroid",
        "Primary neutron energy vs shower depth centroid;"
        "Primary neutron energy (GeV);Depth centroid (mm)",
        100, 0.1, 10.0,
        140, 0.0, 70.0,
        "none", "none", "none", "none", "log", "linear");         // H2-9

    analysis->CreateH2(
        "PrimaryE_vs_RadialRMS",
        "Primary neutron energy vs transverse shower RMS;"
        "Primary neutron energy (GeV);Radial RMS (mm)",
        100, 0.1, 10.0,
        200, 0.0, 100.0,
        "none", "none", "none", "none", "log", "linear");         // H2-10

    // ------------------------------------------------------------------
    // Ntuple 0: one row per incident neutron.
    // Columns 0-23 preserve the previous layout exactly.
    // ------------------------------------------------------------------
    analysis->CreateNtuple(
        "NeutronStudy",
        "Per-event BGO calorimeter response to solar neutrons");

    analysis->CreateNtupleIColumn("eventID");                    //  0
    analysis->CreateNtupleDColumn("primaryE_GeV");               //  1
    analysis->CreateNtupleDColumn("edepGAGG_MeV");               //  2 = BGO
    analysis->CreateNtupleDColumn("edepBi_MeV");                 //  3 = 0
    analysis->CreateNtupleDColumn("edepBiDesc_GAGG_MeV");        //  4 = 0
    analysis->CreateNtupleDColumn("edepBiFrag_GAGG_MeV");        //  5 = 0
    analysis->CreateNtupleDColumn("edepNP_GAGG_MeV");            //  6
    analysis->CreateNtupleDColumn("edepHeavyIon_GAGG_MeV");      //  7
    analysis->CreateNtupleDColumn("edepEM_GAGG_MeV");            //  8
    analysis->CreateNtupleDColumn("edepPiK_GAGG_MeV");           //  9
    analysis->CreateNtupleDColumn("edepLightIon_GAGG_MeV");      // 10
    analysis->CreateNtupleDColumn("closureCheck_MeV");           // 11
    analysis->CreateNtupleIColumn("hadronicInGAGG");             // 12 = any BGO
    analysis->CreateNtupleIColumn("hadronicInBi");               // 13 = 0
    analysis->CreateNtupleIColumn("anyNeutronInelasticInBi");    // 14 = 0
    analysis->CreateNtupleIColumn("nNeutronInelasticInBi");      // 15 = 0
    analysis->CreateNtupleIColumn("biFragmentationEvent");       // 16 = 0
    analysis->CreateNtupleIColumn("nBiFragsEnteringGAGG");       // 17 = 0
    analysis->CreateNtupleIColumn("nBiFragsDepositingGAGG");     // 18 = 0
    analysis->CreateNtupleIColumn("biFragEnteredGAGG");          // 19 = 0
    analysis->CreateNtupleIColumn("biFragDepositedGAGG");        // 20 = 0
    analysis->CreateNtupleIColumn("nBiFragsTotal");              // 21 = 0
    analysis->CreateNtupleIColumn("nBiDescendants");             // 22 = 0
    analysis->CreateNtupleDColumn("weight");                     // 23

    analysis->CreateNtupleDColumn("edepOther_BGO_MeV");          // 24
    analysis->CreateNtupleDColumn("edepFraction");               // 25
    analysis->CreateNtupleDColumn("sourceX_mm");                 // 26
    analysis->CreateNtupleDColumn("sourceY_mm");                 // 27
    analysis->CreateNtupleIColumn("nHits");                      // 28
    analysis->CreateNtupleIColumn("nUniqueTracks");              // 29
    analysis->CreateNtupleIColumn("nSecondaries");               // 30
    analysis->CreateNtupleIColumn("nSecondaryNeutrons");         // 31
    analysis->CreateNtupleIColumn("nSecondaryProtons");          // 32
    analysis->CreateNtupleIColumn("nSecondaryGammas");           // 33
    analysis->CreateNtupleIColumn("nSecondaryElectrons");        // 34
    analysis->CreateNtupleIColumn("nSecondaryPositrons");        // 35
    analysis->CreateNtupleIColumn("nSecondaryPions");            // 36
    analysis->CreateNtupleIColumn("nSecondaryKaons");            // 37
    analysis->CreateNtupleIColumn("nSecondaryLightIons");        // 38
    analysis->CreateNtupleIColumn("nSecondaryHeavyIons");        // 39
    analysis->CreateNtupleIColumn("nSecondaryOther");            // 40
    analysis->CreateNtupleDColumn("sumSecondaryKE_MeV");         // 41
    analysis->CreateNtupleIColumn("primaryNeutronHadronic");     // 42
    analysis->CreateNtupleIColumn("primaryNeutronInelastic");    // 43
    analysis->CreateNtupleIColumn("nPrimaryNeutronInelastic");   // 44
    analysis->CreateNtupleIColumn("firstInteractionProcessCode");// 45
    analysis->CreateNtupleDColumn("firstInteractionDepth_mm");   // 46
    analysis->CreateNtupleDColumn("firstInteractionR_mm");       // 47
    analysis->CreateNtupleDColumn("firstInteractionTime_ns");    // 48
    analysis->CreateNtupleDColumn("firstInteractionPreKE_MeV");  // 49
    analysis->CreateNtupleIColumn("firstInteractionMultiplicity");//50
    analysis->CreateNtupleDColumn("firstInteractionSecondaryKE_MeV");//51
    analysis->CreateNtupleDColumn("totalTrackLength_mm");        // 52
    analysis->CreateNtupleDColumn("maxStepEdep_MeV");            // 53
    analysis->CreateNtupleDColumn("depthCentroid_mm");           // 54
    analysis->CreateNtupleDColumn("depthRMS_mm");                // 55
    analysis->CreateNtupleDColumn("radialCentroid_mm");          // 56
    analysis->CreateNtupleDColumn("radialRMS_mm");               // 57
    analysis->CreateNtupleDColumn("edepFrontThird_MeV");         // 58
    analysis->CreateNtupleDColumn("edepMiddleThird_MeV");        // 59
    analysis->CreateNtupleDColumn("edepRearThird_MeV");          // 60
    analysis->CreateNtupleDColumn("firstHitTime_ns");            // 61
    analysis->CreateNtupleDColumn("lastHitTime_ns");             // 62
    analysis->CreateNtupleDColumn("timeSpan_ns");                // 63
    analysis->CreateNtupleDColumn("promptEdep10ns_MeV");         // 64
    analysis->CreateNtupleDColumn("delayedEdepAfter10ns_MeV");   // 65
    analysis->CreateNtupleDColumn("leakKE_total_MeV");           // 66
    analysis->CreateNtupleDColumn("leakKE_np_MeV");              // 67
    analysis->CreateNtupleDColumn("leakKE_em_MeV");              // 68
    analysis->CreateNtupleDColumn("leakKE_lightIon_MeV");        // 69
    analysis->CreateNtupleDColumn("leakKE_heavyIon_MeV");        // 70
    analysis->CreateNtupleDColumn("leakKE_piK_MeV");             // 71
    analysis->CreateNtupleIColumn("primaryNeutronEscaped");      // 72
    analysis->CreateNtupleDColumn("primaryNeutronExitKE_MeV");   // 73
    analysis->CreateNtupleIColumn("maxCrystalID");               // 74
    analysis->CreateNtupleDColumn("maxCrystalEdep_MeV");         // 75
    analysis->CreateNtupleIColumn("nHitCrystals");               // 76
    analysis->CreateNtupleDColumn("edepLayer0_MeV");             // 77
    analysis->CreateNtupleDColumn("edepLayer1_MeV");             // 78
    analysis->FinishNtuple(0);

    // ------------------------------------------------------------------
    // Ntuple 1: one row per secondary above 0.1 MeV from the first selected
    // primary-neutron hadronic interaction. Complete multiplicities remain
    // in ntuple 0 and are not subject to this size-control threshold.
    // ------------------------------------------------------------------
    analysis->CreateNtuple(
        "FirstInteractionSecondaries",
        "Secondaries from first primary-neutron hadronic interaction");
    analysis->CreateNtupleIColumn("eventID");                //  0
    analysis->CreateNtupleDColumn("primaryE_GeV");           //  1
    analysis->CreateNtupleDColumn("weight");                 //  2
    analysis->CreateNtupleIColumn("processCode");            //  3
    analysis->CreateNtupleIColumn("secondaryIndex");         //  4
    analysis->CreateNtupleIColumn("pdg");                    //  5
    analysis->CreateNtupleIColumn("Z");                      //  6
    analysis->CreateNtupleIColumn("A");                      //  7
    analysis->CreateNtupleDColumn("kineticEnergy_MeV");      //  8
    analysis->CreateNtupleDColumn("x_mm");                   //  9
    analysis->CreateNtupleDColumn("y_mm");                   // 10
    analysis->CreateNtupleDColumn("z_mm");                   // 11
    analysis->CreateNtupleDColumn("depth_mm");               // 12
    analysis->CreateNtupleDColumn("time_ns");                // 13
    analysis->CreateNtupleIColumn("parentTrackID");          // 14
    analysis->CreateNtupleIColumn("secondaryTrackID");       // 15
    analysis->CreateNtupleIColumn("particleClass");          // 16
    analysis->CreateNtupleIColumn("isFragment");             // 17
    analysis->FinishNtuple(1);
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    auto* analysis = G4AnalysisManager::Instance();
    const G4String fileName =
        "bgo_run" + std::to_string(run->GetRunID()) + ".root";

    analysis->OpenFile(fileName);
    G4cout << "[RunAction] Output file: " << fileName << G4endl;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    auto* analysis = G4AnalysisManager::Instance();
    analysis->Write();
    analysis->CloseFile();

    G4cout << "[RunAction] Run " << run->GetRunID()
           << " (" << run->GetNumberOfEvent()
           << " events) complete. File closed." << G4endl;
}
