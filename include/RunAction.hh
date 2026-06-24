#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction  (const G4Run* run) override;
};
