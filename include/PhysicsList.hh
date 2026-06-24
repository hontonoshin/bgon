#pragma once
#include "G4VModularPhysicsList.hh"

class PhysicsList : public G4VModularPhysicsList
{
public:
    PhysicsList();
    ~PhysicsList() override = default;

    void ConstructParticle() override;
    void ConstructProcess() override;
    void SetCuts() override;
};
