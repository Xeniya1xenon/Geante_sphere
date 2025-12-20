#include "G4VModularPhysicsList.hh"
#include "G4OpticalPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4SystemOfUnits.hh"


class PhysicsList : public G4VModularPhysicsList
{
public:
    PhysicsList();
    ~PhysicsList() override = default;

    void ConstructProcess() override;
};

PhysicsList::PhysicsList()
{
    // Стандартная электромагнитная физика
    RegisterPhysics(new G4EmStandardPhysics());

    // Распады
    RegisterPhysics(new G4DecayPhysics());

    // Оптические процессы (ОБЯЗАТЕЛЬНО для сцинтилляции и ФЭУ)
    RegisterPhysics(new G4OpticalPhysics());

    SetVerboseLevel(1);
}

void PhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();
}

//#include "G4EmStandardPhysics.hh"
//#include "G4OpticalPhysics.hh"
//#include "G4DecayPhysics.hh"
//#include "G4RadioactiveDecayPhysics.hh"
//#include "G4VModularPhysicsList.hh"
//#include "G4SystemOfUnits.hh"



//void PhysicsList::ConstructProcess()
//{
  //  G4VModularPhysicsList::ConstructProcess();

    //auto opticalPhysics = new G4OpticalPhysics();
    //RegisterPhysics(opticalPhysics);
//}