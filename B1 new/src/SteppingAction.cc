/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4OpticalPhoton.hh"

namespace B1
{

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}


G4double efficiency(G4double energy) {


  G4double photonEnergy[] =
{1.70455, 1.722, 1.73442, 1.75596, 1.778, 1.81388, 1.86867,
  1.884, 1.908, 1.957, 2.0131, 2.051, 2.20, 2.304, 2.37546,
  2.54, 2.744, 2.905, 3.217, 3.394, 3.62987, 3.94856, 4.062,
  4.19989, 4.2724, 4.32, 4.6
};


   G4double efficiencyPMT[] =
{
  0.01119, 0.000181, 0.000237976, 0.00044485, 0.0008015, 0.00183, 0.0032455,
  0.0049557, 0.00893, 0.0158, 0.0255, 0.033, 0.06647, 0.113397,
  0.157976, 0.204579, 0.269869, 0.31288, 0.3434557, 0.35,
  0.33765, 0.2472, 0.18765, 0.1023, 0.068275, 0.04005, 0.00215
};

        G4double dE = abs(photonEnergy[0]-energy);
  int j=0;
  for (int i=1; i<27; i++)
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
   G4Track *aTrack=step->GetTrack();

    G4LogicalVolume* volume 
      = step->GetPreStepPoint()->GetTouchableHandle()
        ->GetVolume()->GetLogicalVolume();

  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  
  // if ((volume == fScoringVolume) && (aTrack->GetDynamicParticle()->GetCharge()!=0)){
  //   // collect energy deposited in this step
  //     G4double edepStep = step->GetTotalEnergyDeposit();
  //     fEventAction->AddEdep(edepStep);
  // }

  if (aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
      return;
    
    // get volume of the current step
    G4StepPoint* postPoint = step->GetPostStepPoint();
    G4VPhysicalVolume* postVolume = postPoint->GetPhysicalVolume();
    G4LogicalVolume* pvolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary 
        && pvolume->GetName()=="logic_detector")
    {
        G4double Ephoton = aTrack->GetKineticEnergy();
        G4double prob_photon = efficiency(Ephoton*1e6);
        if (G4UniformRand() < prob_photon) {
          fEventAction->AddPhotoElectron();
        }    
    } 

}


}  // namespace B1
