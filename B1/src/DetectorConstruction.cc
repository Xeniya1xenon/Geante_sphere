//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();



  // Envelope parameters
  //
  G4double env_sizeXY = 400 * cm, env_sizeZ = 400 * cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;



  //
  // World
  //
  G4double world_sizeXY = 1.0 * env_sizeXY;
  G4double world_sizeZ = 1.0 * env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld =
    new G4Box("World",  // its name
              0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        world_mat,  // its material
                                        "World");  // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                     G4ThreeVector(),  // at (0,0,0)
                                     logicWorld,  // its logical volume
                                     "World",  // its name
                                     nullptr,  // its mother  volume
                                     false,  // no boolean operation
                                     0,  // copy number
                                     checkOverlaps);  // overlaps checking

  //
  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",  // its name
                            0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
                                      env_mat,  // its material
                                      "Envelope");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicEnv,  // its logical volume
                    "Envelope",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking



//LAB scintallor               
 G4Element* elH = nist->FindOrBuildElement("H");
 G4Element* elC = nist->FindOrBuildElement("C");

 G4double density = 0.856 * g/cm3;
 G4Material* LAB = new G4Material("LAB", density, 2);
 LAB->AddElement(elH, 0.12);
 LAB->AddElement(elC, 0.88);

const G4int NUMENTRIES = 9; // scint properties
 G4double LAB_PP[NUMENTRIES] = 
{6.6*eV,6.7*eV,6.8*eV,6.9*eV,7.0*eV, 7.1*eV,7.2*eV,7.3*eV,7.4*eV}; // энергия выделившегося при сцинтилляции фотона
 G4double LAB_RIND[NUMENTRIES] = 
{ 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57}; // индекс преломления
 G4double LAB_ABSL[NUMENTRIES] = 
{ 35.*cm, 35.*cm, 35.*cm, 35.*cm, 35.*cm, 35.*cm, 35.*cm, 35.*cm, 35*cm }; // длина пробега волны

 G4MaterialPropertiesTable* LAB_MPT = new G4MaterialPropertiesTable();

LAB_MPT->AddProperty("RINDEX", LAB_PP, LAB_RIND, NUMENTRIES);
LAB_MPT->AddProperty("ABSLENGTH", LAB_PP, LAB_ABSL, NUMENTRIES);
LAB_MPT->AddConstProperty("SCINTILLATIONYIELD", 8000./MeV);


 LAB -> SetMaterialPropertiesTable(LAB_MPT);


//PMMA
 G4Element* elO = nist->FindOrBuildElement("O");

 G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
    PMMA -> AddElement(elC, 5);
    PMMA -> AddElement(elO, 2);
    PMMA -> AddElement(elH, 8);


//Сцинтиллятор LAB + PPO + Bis-MSB + Gd


G4Element* elGd = nist->FindOrBuildElement("Gd");

//Создаем LAB с добавками

G4Material* LAB_sc = new G4Material("LAB_Scintillator", 0.856 * g/cm3, 3);
LAB_sc->AddElement(elC, 0.8788);
LAB_sc->AddElement(elH, 0.12);
LAB_sc->AddElement(elGd, 0.0012); // добавка гадолиния (1 г/л)

//Оптические свойства
const G4int NUM = 9;
G4double photonEnergy[NUM] = {6.6*eV, 6.7*eV, 6.8*eV, 6.9*eV, 7.0*eV, 7.1*eV, 7.2*eV, 7.3*eV, 7.4*eV};

//Индекс преломления
G4double refractiveIndex[NUM] = {1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57};

//Длина поглощения света 
G4double absorptionLength[NUM] = {600*cm, 600*cm, 600*cm, 600*cm, 600*cm, 600*cm, 600*cm, 600*cm, 600*cm};

//Эмиссионный спектр 
G4double scintEmission[NUM] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

//Таблица свойств
G4MaterialPropertiesTable* MPT_LAB = new G4MaterialPropertiesTable();
MPT_LAB->AddProperty("RINDEX", photonEnergy, refractiveIndex, NUM);
MPT_LAB->AddProperty("ABSLENGTH", photonEnergy, absorptionLength, NUM);
MPT_LAB->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintEmission, NUM);
MPT_LAB->AddConstProperty("SCINTILLATIONYIELD", 8000./MeV);
MPT_LAB->AddConstProperty("RESOLUTIONSCALE", 1.0);
MPT_LAB->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 3.0*ns);
MPT_LAB->AddConstProperty("SCINTILLATIONYIELD1", 1.0);

LAB_sc->SetMaterialPropertiesTable(MPT_LAB);



  //
  // Shape 1
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4ThreeVector pos1 = G4ThreeVector(0, 0 * cm, 0 * cm);

  // Conical section shape
  G4double shape1_rmina = 0. * cm, shape1_rmaxa = 929. * mm;
  G4double shape1_rminb = 0. * cm, shape1_rmaxb = 929. * mm;
  G4double shape1_hz = 929. * mm;
  G4double shape1_phimin = 0. * deg, shape1_phimax = 360. * deg;
  auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
                                shape1_hz, shape1_phimin, shape1_phimax);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
                                         shape1_mat,  // its material
                                         "Shape1");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    pos1,  // at position
                    logicShape1,  // its logical volume
                    "Shape1",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking

  // 
  // Shape 2
  //
  G4Material* shape2_mat = LAB;
  G4ThreeVector pos2 = G4ThreeVector(0, 0 * cm, 0 * cm);

  // Conical section shape
  G4double shape2_rmina = 0. * cm, shape2_rmaxa = 927. * mm;
  G4double shape2_rminb = 0. * cm, shape2_rmaxb = 927. * mm;
  G4double shape2_hz = 927. * mm;
  G4double shape2_phimin = 0. * deg, shape2_phimax = 360. * deg;
  auto solidShape2 = new G4Cons("Shape2", shape2_rmina, shape2_rmaxa, shape2_rminb, shape2_rmaxb,
                                shape2_hz, shape2_phimin, shape2_phimax);

  auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
                                         shape2_mat,  // its material
                                         "Shape2");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    pos2,  // at position
                    logicShape2,  // its logical volume
                    "Shape2",  // its name
                    logicShape1,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking


  
  // 
  // Shape 3
  //
  G4Material* shape3_mat = PMMA; // это сосуд для сферического лаба, нцжно замениить материал и прописать его выше
  G4ThreeVector pos3 = G4ThreeVector(0, 0 * cm, 0 * cm);

  G4double shape3_Rmax = 630.35 * mm;
  auto solidShape3 = new G4Orb("Shape3", shape3_Rmax);

  auto logicShape3 = new G4LogicalVolume(solidShape3,  // its solid
                                         shape3_mat,  // its material
                                         "Shape3");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    pos3,  // at position
                    logicShape3,  // its logical volume
                    "Shape3",  // its name
                    logicShape2,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking

  // 
  // Shape 4
  //
   G4Material* shape4_mat = LAB_sc;
  G4ThreeVector pos4 = G4ThreeVector(0, 0 * cm, 0 * cm);

  G4double shape4_Rmax = 620.35 * mm; //cltkfnm gjvtymit
  auto solidShape4 = new G4Orb("Shape4", shape4_Rmax);

  auto logicShape4 = new G4LogicalVolume(solidShape4,  // its solid
                                         shape4_mat,  // its material
                                         "Shape4");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    pos4,  // at position
                    logicShape4,  // its logical volume
                    "Shape4",  // its name
                    logicShape3,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking


  //
  // Shape 2
  // //
  // G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  // G4ThreeVector pos2 = G4ThreeVector(0, -1 * cm, 7 * cm);

  // // Trapezoid shape
  // G4double shape2_dxa = 12 * cm, shape2_dxb = 12 * cm;
  // G4double shape2_dya = 10 * cm, shape2_dyb = 16 * cm;
  // G4double shape2_dz = 6 * cm;
  // auto solidShape2 =
  //   new G4Trd("Shape2",  // its name
  //             0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
  //             0.5 * shape2_dz);  // its size

  // auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
  //                                        shape2_mat,  // its material
  //                                        "Shape2");  // its name

  // new G4PVPlacement(nullptr,  // no rotation
  //                   pos2,  // at position
  //                   logicShape2,  // its logical volume
  //                   "Shape2",  // its name
  //                   logicEnv,  // its mother  volume
  //                   false,  // no boolean operation
  //                   0,  // copy number
  //                   checkOverlaps);  // overlaps checking

  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicShape4;
  
  //
  // always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
