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

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include <vector>
#include <cmath>


class PMTSD : public G4VSensitiveDetector
{
public:
    PMTSD(const G4String& name)
        : G4VSensitiveDetector(name) {}

    ~PMTSD() override = default;

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override
    {
        auto track = step->GetTrack();

        // Считаем только оптические фотоны
        if (track->GetDefinition()->GetParticleName() == "opticalphoton") {
            track->SetTrackStatus(fStopAndKill);
        }
        return true;
    }
};

namespace B1
{
//новое
  class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction() = default;
    virtual ~DetectorConstruction() = default;

    virtual G4VPhysicalVolume* Construct();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

private:
    G4LogicalVolume* fScoringVolume = nullptr;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

const G4int NENT = 9;
static G4double photonEnergy[NENT] = 
{1.70455*eV, 1.73442*eV, 1.81388*eV, 1.908*eV, 2.304*eV, 2.744*eV, 3.394*eV, 4.062*eV, 4.6*eV};


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

auto logicEnv = new G4LogicalVolume(
    solidEnv,
    env_mat,
    "Envelope"
);

auto physEnv = new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    logicEnv,
    "Envelope",
    logicWorld,
    false,
    0,
    checkOverlaps
);



// Воздух (мир)

static G4double RINDEX_AIR[NENT] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
G4MaterialPropertiesTable* MPT_AIR = new G4MaterialPropertiesTable();
MPT_AIR->AddProperty("RINDEX", photonEnergy, RINDEX_AIR, NENT);
env_mat->SetMaterialPropertiesTable(MPT_AIR);
world_mat->SetMaterialPropertiesTable(MPT_AIR);



//LAB scintallor               
 G4Element* elH = nist->FindOrBuildElement("H");
 G4Element* elC = nist->FindOrBuildElement("C");

 G4double density = 0.856 * g/cm3;
 G4Material* LAB = new G4Material("LAB", density, 2);
 LAB->AddElement(elH, 0.12);
 LAB->AddElement(elC, 0.88);

const G4int NUMENTRIES = 9; // scint properties
static G4double LAB_PP[NUMENTRIES] = 
{1.70455*eV, 1.73442*eV, 1.81388*eV, 1.908*eV, 2.304*eV, 2.744*eV, 3.394*eV, 4.062*eV, 4.6*eV}; // энергия выделившегося при сцинтилляции фотона
static G4double LAB_RIND[NUMENTRIES] = 
{ 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57}; // индекс преломления
static G4double LAB_ABSL[NUMENTRIES] = 
{ 1200.*cm, 1200.*cm, 1200.*cm, 1200.*cm, 1200.*cm, 1200.*cm, 1200.*cm, 1200.*cm, 1200.*cm }; // длина пробега волны

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

  static G4double RINDEX_PMMA[NENT] = {1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49};
  static G4double ABSL_PMMA[NENT]   = {200*cm,200*cm,200*cm,200*cm,200*cm,200*cm,200*cm,200*cm,200*cm};

G4MaterialPropertiesTable* MPT_PMMA = new G4MaterialPropertiesTable();
MPT_PMMA->AddProperty("RINDEX", photonEnergy, RINDEX_PMMA, NENT);
MPT_PMMA->AddProperty("ABSLENGTH", photonEnergy, ABSL_PMMA, NENT);
PMMA->SetMaterialPropertiesTable(MPT_PMMA);


//Сцинтиллятор LAB + PPO + Bis-MSB + Gd


G4Element* elGd = nist->FindOrBuildElement("Gd");

//Создаем LAB с добавками

G4Material* LAB_sc = new G4Material("LAB_Scintillator", 0.856 * g/cm3, 3);
LAB_sc->AddElement(elC, 0.8788);
LAB_sc->AddElement(elH, 0.12);
LAB_sc->AddElement(elGd, 0.0012); // добавка гадолиния (1 г/л)


//Индекс преломления
static G4double refractiveIndex[NUMENTRIES] = {1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57};

//Длина поглощения света 
static G4double absorptionLength[NUMENTRIES] = {600*cm, 600*cm, 600*cm, 600*cm, 600*cm, 600*cm, 600*cm, 600*cm, 600*cm};

//Эмиссионный спектр 
static G4double scintEmission[NUMENTRIES] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

//Таблица свойств
G4MaterialPropertiesTable* MPT_LAB = new G4MaterialPropertiesTable();
MPT_LAB->AddProperty("RINDEX", photonEnergy, refractiveIndex, NUMENTRIES);
MPT_LAB->AddProperty("ABSLENGTH", photonEnergy, absorptionLength, NUMENTRIES);
MPT_LAB->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintEmission, NUMENTRIES);
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

//  Нержавеющая сталь (отражающая поверхность)

G4MaterialPropertiesTable* MPT_STEEL = new G4MaterialPropertiesTable();
static G4double RINDEX_STEEL[NENT] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5};
MPT_STEEL->AddProperty("RINDEX", photonEnergy, RINDEX_STEEL, NENT);
shape1_mat->SetMaterialPropertiesTable(MPT_STEEL);

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
    auto physShape2 = new G4PVPlacement(
    nullptr,
    pos2,
    logicShape2,
    "Shape2",
    logicShape1,
    false,
    0,
    checkOverlaps
);



// физические размеры (Hamamatsu R5912)
G4double PMT_outer_diam = 204.0 * mm;    // диаметр колбы
G4double PMT_outer_rad  = PMT_outer_diam/2.0; // ~102 mm
G4double Photocathode_diam = 190.0 * mm; // активная область
G4double Photocathode_rad  = Photocathode_diam/2.0; // ~95 mm
G4double Photocathode_thickness = 0.2 * mm; // тонкая пленка

// 
const G4int PMT_N = 9;
static G4double pmtPhotonEnergy[PMT_N] = {
    1.70455*eV, 1.73442*eV, 1.81388*eV, 1.908*eV, 2.304*eV, 2.744*eV, 3.394*eV, 4.062*eV, 4.6*eV
}; 

// QE
static G4double pmtQE[PMT_N] = { 0.01, 0.03, 0.10, 0.20, 0.26, 0.28, 0.26, 0.20, 0.08
};

// REFLECTIVITY
static G4double pmtRefl[PMT_N];
for (int i=0;i<PMT_N;i++) pmtRefl[i] = 0.0; // можно выставить небольшую отраж. при желании

// 
G4double glassRindex[PMT_N];
G4double glassABSL[PMT_N];
for (int i=0;i<PMT_N;i++){
    glassRindex[i] = 1.47;    
    glassABSL[i]    = 50.*cm; 
}
G4Material* glassMat = nist->FindOrBuildMaterial("G4_BOROSILICATE_GLASS");
G4MaterialPropertiesTable* glassMPT = new G4MaterialPropertiesTable();
glassMPT->AddProperty("RINDEX", pmtPhotonEnergy, glassRindex, PMT_N);
glassMPT->AddProperty("ABSLENGTH", pmtPhotonEnergy, glassABSL, PMT_N);
glassMat->SetMaterialPropertiesTable(glassMPT);

// Photocathode (тонкая пленка) 
//
G4Element* elK = nist->FindOrBuildElement("K");
G4Material* photocathodeMat = new G4Material("PhotocathodeMat", 5.0*g/cm3, 1);
photocathodeMat->AddElement(elK, 1); // не обязательно хим. корректно — имя условное
static G4double pcRindex[PMT_N];
for (int i=0;i<PMT_N;i++) pcRindex[i] = 1.0;

G4MaterialPropertiesTable* pcMPT = new G4MaterialPropertiesTable();
pcMPT->AddProperty("RINDEX", pmtPhotonEnergy, pcRindex, PMT_N);
photocathodeMat->SetMaterialPropertiesTable(pcMPT);
// 

// Optical surface на фотокатоде:
G4OpticalSurface* photocathodeSurface = new G4OpticalSurface("PhotocathodeSurface");
photocathodeSurface->SetType(dielectric_metal);    // фотокатод — поглощающий/эмиссионный слой
photocathodeSurface->SetModel(unified);
photocathodeSurface->SetFinish(polished);          

G4MaterialPropertiesTable* photocMPT = new G4MaterialPropertiesTable();
photocMPT->AddProperty("EFFICIENCY", pmtPhotonEnergy, pmtQE, PMT_N);   // QE(λ)
photocMPT->AddProperty("REFLECTIVITY", pmtPhotonEnergy, pmtRefl, PMT_N);
photocathodeSurface->SetMaterialPropertiesTable(photocMPT);

// Геометрия PMT 
//
G4double rmin = 0.0*mm;
G4double rmax = PMT_outer_rad;
G4double startPhi = 0.*deg, deltaPhi = 360.*deg;
G4double startTheta = 0.*deg, deltaTheta = 180.*deg; 
auto solidGlass = new G4Sphere("PMTGlass_solid", rmin, rmax, startPhi, deltaPhi, 0.*deg, 180.*deg);

// Логический объем стекла
auto logicGlass = new G4LogicalVolume(solidGlass, glassMat, "PMTGlass_log");

// Фотокатод — тонкая оболочка: внутренний радиус = Photocathode_rad, толщина = Photocathode_thickness
auto pci_rmin = Photocathode_rad;
auto pci_rmax = Photocathode_rad + Photocathode_thickness;
auto solidPhotocat = new G4Sphere("Photocat_solid", pci_rmin, pci_rmax, startPhi, deltaPhi, 0.*deg, 180.*deg);
auto logicPhotocat = new G4LogicalVolume(solidPhotocat, photocathodeMat, "Photocat_log");

// Связываем skin surface с логическим объёмом фотокатода
new G4LogicalSkinSurface("PhotocatSkin", logicPhotocat, photocathodeSurface);


//  Размещение и ориентация PMT (12 вершин икосаэдра)
G4double R_PMT_place = 780.0*mm; // радиус размещения (вблизи PMMA)
G4double phi = (1.0 + std::sqrt(5.0)) / 2.0;
std::vector<G4ThreeVector> vertices = {
    { -1,  phi,  0}, {  1,  phi,  0}, { -1, -phi,  0}, {  1, -phi,  0},
    {  0, -1,  phi}, {  0,  1,  phi}, {  0, -1, -phi}, {  0,  1, -phi},
    {  phi,  0, -1}, {  phi,  0,  1}, { -phi,  0, -1}, { -phi,  0,  1}
};
for (auto& v : vertices) v = v.unit() * R_PMT_place;

  

// поворот фэу
for (size_t i=0;i<vertices.size(); ++i) {

auto rot = new G4RotationMatrix();

// Разворачиваем PMT "внутрь"
rot->rotateY(180.0 * deg);

// Поворачиваем вокруг оси Z по азимутальному углу
rot->rotateZ(vertices[i].phi());
  
   G4ThreeVector pos = vertices[i];

    /* G4ThreeVector w = -pos.unit();
    G4ThreeVector u = G4ThreeVector(0,0,1).cross(w);
    if (u.mag() == 0) u = G4ThreeVector(0,1,0);
    u = u.unit();
    G4ThreeVector v = w.cross(u);

    auto rot = new G4RotationMatrix();
    rot->set(u, v, w); */

    auto physGlass = new G4PVPlacement(
        rot,
        pos,
        logicGlass,
        "PMT_glass",
        logicShape1,
        false,
        (G4int)i,
        checkOverlaps
    );

    auto physPhotocat = new G4PVPlacement(
        nullptr,
        G4ThreeVector(),
        logicPhotocat,
        "PMT_photocat",
        logicGlass,
        false,
        (G4int)i,
        checkOverlaps
    );

   // new G4LogicalBorderSurface(
     //   "GlassToPhotocat",
       // physGlass,
        //physPhotocat,
        //photocathodeSurface
    //);

   /* G4ThreeVector pos = vertices[i];

    // ориентация: ось Z локального PMT направлена в сторону центра
    G4ThreeVector w = -pos.unit();   // направление на центр
    G4ThreeVector u = G4ThreeVector(0,0,1).cross(w);
    if (u.mag() == 0) u = G4ThreeVector(0,1,0);
    u = u.unit();
    G4ThreeVector v = w.cross(u);
    auto rot = new G4RotationMatrix();
    rot->set(u, v, w);

    // glass placement (внутри logicShape2, т.е. в LAB без добавок)
    new G4PVPlacement(rot,
                      pos, 
                      logicGlass, 
                      "PMT_glass", 
                      logicShape2, 
                      false, 
                      (G4int)i, 
                      checkOverlaps);

    // фотокатод размещаем как дочка glass в тех же локальных координатах (центр совпадает)
    auto physPhotocat = new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    logicPhotocat,
    "PMT_photocat",
    logicGlass,
    false,
    i,
    checkOverlaps
); */
}

// чувствительный детектор только к фотокатоду
G4SDManager* SDman = G4SDManager::GetSDMpointer();
PMTSD* pmtSD = new PMTSD("R5912_PMT_SD");
SDman->AddNewDetector(pmtSD);
logicPhotocat->SetSensitiveDetector(pmtSD);


  
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
auto physShape3 = new G4PVPlacement(
    nullptr,
    pos3,
    logicShape3,
    "Shape3",
    logicShape2,
    false,
    0,
    checkOverlaps
);

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

  auto physShape4 = new G4PVPlacement(
    nullptr,
    pos4,
    logicShape4,
    "Shape4",
    logicShape3,
    false,
    0,
    checkOverlaps
);

// LAB_sc → PMMA
auto scint_to_pmma = new G4OpticalSurface("ScintToPMMA");
scint_to_pmma->SetType(dielectric_dielectric);
scint_to_pmma->SetModel(unified);
scint_to_pmma->SetFinish(polished);

// PMMA → LAB
auto pmma_to_lab = new G4OpticalSurface("PMMAToLAB");
pmma_to_lab->SetType(dielectric_dielectric);
pmma_to_lab->SetModel(unified);
pmma_to_lab->SetFinish(ground);

// LAB → air (steel tank)
auto lab_to_air = new G4OpticalSurface("LABToAir");
lab_to_air->SetType(dielectric_metal);
lab_to_air->SetModel(unified);
lab_to_air->SetFinish(groundfrontpainted);

// ---------- Border surfaces ----------

new G4LogicalBorderSurface("ScintToPMMA_Surface",
                           physShape4,   // LAB_sc
                           physShape3,   // PMMA
                           scint_to_pmma);

new G4LogicalBorderSurface("PMMAToLAB_Surface",
                           physShape3,   // PMMA
                           physShape2,   // LAB
                           pmma_to_lab);

new G4LogicalBorderSurface("LABToAir_Surface",
                           physShape2,   // LAB
                           physEnv,      // воздух
                           lab_to_air);


/*
//  Граница между сцинтиллятором и PMMA
G4OpticalSurface* scint_to_pmma = new G4OpticalSurface("ScintToPMMA");
scint_to_pmma->SetType(dielectric_dielectric);
scint_to_pmma->SetModel(unified);
scint_to_pmma->SetFinish(polished);

new G4LogicalBorderSurface(
    "ScintToPMMA_Surface",
    physShape4,   // внутренняя сфера (LAB_sc)
    physShape3,   // PMMA
    scint_to_pmma
);

// 2. Граница между PMMA и внешним LAB
G4OpticalSurface* pmma_to_lab = new G4OpticalSurface("PMMAToLAB");
pmma_to_lab->SetType(dielectric_dielectric);
pmma_to_lab->SetModel(unified);
pmma_to_lab->SetFinish(ground);

new G4LogicalBorderSurface(
    "PMMAToLAB_Surface",
    physShape3,   // PMMA
    physShape2,   // LAB
    pmma_to_lab
);

// 3. Граница LAB — воздух (стенки бака)
G4OpticalSurface* lab_to_air = new G4OpticalSurface("LABToAir");
lab_to_air->SetType(dielectric_metal);
lab_to_air->SetFinish(groundfrontpainted);
lab_to_air->SetModel(unified);

new G4LogicalBorderSurface(
    "LABToAir_Surface",
    physShape2,
    physEnv,
    lab_to_air   
);
*/

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
