/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"

namespace B1
{

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}


G4double efficiency(G4double energy) {

  G4double photonEnergy[] =

  {4.6, 4.32, 4.2724, 4.19989,
   4.062, 3.94856, 3.62987, 3.394,
   3.217, 2.905, 2.744, 2.54,
   2.37546, 2.304, 2.20, 2.051,
   2.0131, 1.957, 1.908, 1.86867,
   1.884, 1.81388, 1.778, 1.75596,
    1.73442, 1.722, 1.70455};


        G4double efficiencyPMT[] =
  {0.00215, 0.04005, 0.068275, 0.1023,
   0.18765, 0.2472, 0.33765, 0.35,
   0.3434557, 0.31288, 0.269869, 0.204579,
   0.157976, 0.113397, 0.06647, 0.033,
   0.0255, 0.0158, 0.00893, 0.0049557,
   0.0032455, 0.00183, 0.0008015, 0.00044485, 
   0.000237976, 0.000181, 0.01119};

        G4double dE = abs(photonEnergy[0]-energy);
  int j=0;
  for (int i=1; i<23; i++)
  {
    if (abs(photonEnergy[i]-energy)<dE)
       {
    dE = photonEnergy[i]-energy;
          j = i;
             }
  }
  return efficiencyPMT[j];
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  // get volume of the current step
  G4LogicalVolume* volume =
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
}


}  // namespace B1
