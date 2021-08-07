 
#include "DetectorConstruction.hh"

#include "MCEvent.hh"
#include "DataManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <vector>
#include <algorithm>

// void myfunction (int i) {
//   G4cout << i<<G4endl;
// }

//-----------------------------------------------------------------------------
 
DetectorConstruction::DetectorConstruction(G4double worldSize)
 
{
  solidWorld =0;
  logicWorld =0; 
  physiWorld =0;
  
  solidWrapping = 0;
  logicWrapping = 0;
  physiWrapping = 0;   

  solidScint = 0;
  logicScint = 0;
  physiScint = 0;   

  solidSensor = 0;
  logicSensor = 0;
  physiSensor = 0;   

  stepLimit = 0;

  NISTManager         = G4NistManager::Instance();

  World_VisAtt   = 0;
  VisAttWrapping = 0;
  VisAttScint    = 0;
  VisAttSensor   = 0;  
 
  fWorldLength        = worldSize;
  //field size devided by (fiber_diameter+fiber_gap) has to be even number to obtain symetric arangement of fibers
  field_size          = 40.455*cm;//10.44*cm;  //40.455*cm for 931 fibers per row; 
                                   //4.35*cm for 101 fibers per row; 
                                   //2.61*cm for 61 fibers per row; 
                                   //2.61*mm for 7 fibers per row; 
                                   //100.*cm;
  scint_x             = 2.0*cm;
  scint_y             = 2.0*cm;
  scint_z             = 4.0*cm;
  wrapping_thinkness  = 3.0*mm;
  sensor_thickness    = 1.0*cm;
  //-----------------------------------------------------------------;
  //DataManager* dataManager = DataManager::GetInstance();
}

//-----------------------------------------------------------------------------
 
DetectorConstruction::~DetectorConstruction()
{
  // if(World_VisAtt!=0)  {delete World_VisAtt;   World_VisAtt   = 0;}
  // if(VisAttWrapping!=0){delete VisAttWrapping; VisAttWrapping = 0;}
  // if(VisAttScint!=0)   {delete VisAttScint;    VisAttScint    = 0;}
  // if(VisAttSensor!=0)  {delete VisAttSensor;   VisAttSensor   = 0;}  

  // if(solidWrapping!=0){delete solidWrapping; solidWrapping = 0;}
  // if(logicWrapping!=0){delete logicWrapping; logicWrapping = 0;}
  // if(physiWrapping!=0){delete physiWrapping; physiWrapping = 0;}

  // if(solidScint!=0){delete solidScint; solidScint = 0;}
  // if(logicScint!=0){delete logicScint; logicScint = 0;}
  // if(physiScint!=0){delete physiScint; physiScint = 0;}

  // if(solidSensor!=0){delete solidSensor; solidSensor = 0;}
  // if(logicSensor!=0){delete logicSensor; logicSensor = 0;}
  // if(physiSensor!=0){delete physiSensor; physiSensor = 0;}
  
  if(solidWorld!=0){delete solidWorld; solidWorld = 0;}
  if(logicWorld!=0){delete logicWorld; logicWorld = 0;}
  if(physiWorld!=0){delete physiWorld; physiWorld = 0;}

  if(stepLimit!=0){delete stepLimit; stepLimit =0;}
}

//-----------------------------------------------------------------------------
 
G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  //------------------------------
  // Construct Materials
  //------------------------------
  ConstructMaterials();

  //------------------------------ 
  // World
  //------------------------------ 
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  G4double HalfWorldLength = 0.5*fWorldLength;

  auto solidWorld = new G4Box("sWorld",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  auto logicWorld = new G4LogicalVolume(solidWorld, NISTManager->FindOrBuildMaterial("G4_AIR"), "lWorld", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  auto physiWorld = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*mm,0.*mm), logicWorld, "pvWorld", 0, false, 0 ,true );   
  //G4cout<<"world center position:  "<<physiWorld->GetObjectTranslation()<<G4endl; 

  //Start building the detector
  //------------------------------------------------------------------
  auto solidWrapping = new G4Box("sWrapping",scint_x/2.+wrapping_thinkness,scint_y/2.+wrapping_thinkness,scint_z/2.+wrapping_thinkness/2.);
  auto logicWrapping = new G4LogicalVolume(solidWrapping, NISTManager->FindOrBuildMaterial("TEFLON_OPTICAL"), "lWrapping", 0, 0, 0);
  auto physiWrapping = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*mm,0.*mm), logicWrapping, "pvWrapping", logicWorld, false, 0);   

  //------------------------------------------------------------------
  solidScint = new G4Box("sScint",scint_x/2.,scint_y/2.,scint_z/2.);
  logicScint = new G4LogicalVolume(solidScint, NISTManager->FindOrBuildMaterial("plasticScint"), "lScint", 0, 0, 0);
  // G4double scintZ = (physiWrapping->GetTranslation()).z()+wrapping_thinkness/2.
  physiScint = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*mm,wrapping_thinkness/2.), logicScint, "pvScint", logicWrapping, false, 0);   

  //--------------------- Sensor ----------------------------------
  solidSensor = new G4Box("sSensor",scint_x/2.,scint_y/2.,sensor_thickness/2.);
  logicSensor = new G4LogicalVolume(solidSensor, NISTManager->FindOrBuildMaterial("G4_AIR"),"lSensor", 0, 0, 0);
  physiSensor = new G4PVPlacement(0, G4ThreeVector(0.*mm,0.*mm,scint_z/2.+sensor_thickness/2.+wrapping_thinkness/2.), 
				  logicSensor, "pvSensor", logicWorld, false, 0);   

  
  //--------- Visualization attributes -------------------------------

  VisualizationAttributes();

  //--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in PhysicsList how to setup the processes
  // G4StepLimiter or G4UserSpecialCuts).
    
  // Sets a max Step length in the tracker region, with G4StepLimiter
  //
  // G4double maxStep = 0.5*cm;
  // stepLimit = new G4UserLimits(maxStep);
  // logicWorld->SetUserLimits(stepLimit);
  
  // Set additional contraints on the track, with G4UserSpecialCuts
  //
  // G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
  return physiWorld;
}
 
//======================================================================
void DetectorConstruction::ConstructMaterials()
{
  NISTManager=G4NistManager::Instance();

  //-----------------------------
  G4Material* scintMat = new G4Material ("plasticScint", 1.06*g/cm3, 1, kStateSolid);
  scintMat->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE"), 1.0);

  //------------------- Material Properties Table for Scint -----------------------------

  //========================================================================//
  //======               Material properties for the Scintillation      ====//
  //====== Refractive index values are also used for cherenkov generation ==//
  //========================================================================//

  const G4int ScintENTRIES = 13;
  G4double ppscint[ScintENTRIES] = {2.48*eV, 2.43*eV, 2.41*eV, 2.40*eV, 2.38*eV, 2.34*eV, 2.26*eV, 
				    2.21*eV, 2.18*eV, 2.10*eV, 2.07*eV, 1.98*eV, 1.91*eV};
    
  G4double rindex_scint[ScintENTRIES]={1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6};
  G4double absorption[ScintENTRIES]={3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m, 3.5*m};
  G4double ScintilFast[ScintENTRIES]={0.019, 0.037, 0.074, 0.111, 0.148, 0.185, 0.148, 
				      0.110, 0.074, 0.037, 0.033, 0.018, 0.006};

  G4MaterialPropertiesTable* MPT_scint = new G4MaterialPropertiesTable();
  MPT_scint->AddProperty("RINDEX", ppscint, rindex_scint, ScintENTRIES);
  MPT_scint->AddProperty("ABSLENGTH",ppscint,absorption, ScintENTRIES);
  MPT_scint->AddProperty("FASTCOMPONENT",ppscint, ScintilFast, ScintENTRIES);
  MPT_scint->AddConstProperty("SCINTILLATIONYIELD",7100./MeV);
  MPT_scint->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPT_scint->AddConstProperty("FASTTIMECONSTANT", 7.0*ns);
  scintMat->SetMaterialPropertiesTable(MPT_scint);
  scintMat->GetIonisation()->SetBirksConstant(0.07943*mm/MeV);
  G4cout<<"Polystyrene radLength = "<<scintMat->GetRadlen()/cm<<G4endl;


  //----====================
  // build TEFLON
  //----====================
  G4Material * Teflon = NISTManager->FindOrBuildMaterial("G4_TEFLON");
  G4Material* Teflon_optical = new G4Material("TEFLON_OPTICAL", Teflon->GetDensity(), Teflon, kStateSolid);
  const G4int numEntriesTeflon = 16;

  G4double TeflonEnergies[numEntriesTeflon] = { 1.90769*eV, 1.9375*eV, 1.96825*eV, 2*eV, 2.03279*eV, 2.06667*eV, 
					    2.10169*eV, 2.13793*eV, 2.17544*eV, 2.21429*eV, 2.25455*eV, 
					    2.2963*eV, 2.33962*eV, 2.38462*eV, 2.43137*eV, 2.48*eV};

  G4double teflonrindices[numEntriesTeflon] = { 1.29999, 1.30009, 1.3002, 1.30031, 1.30042, 1.30053, 1.30064, 1.30096, 1.30128, 
					      1.30161, 1.30193, 1.30226, 1.30251, 1.30274, 1.30297, 1.3032}; 
  //DOI:10.1117/1.2965541   //J. Micro/Nanolith. MEMS MOEMS 7 3, 033010 ͑Jul–Sep 2008

  G4MaterialPropertiesTable* teflonprop = new G4MaterialPropertiesTable();
  teflonprop->AddProperty("RINDEX", TeflonEnergies, teflonrindices, numEntriesTeflon);
  Teflon_optical->SetMaterialPropertiesTable(teflonprop);

}
//======================================================================

void DetectorConstruction::VisualizationAttributes()
{

  World_VisAtt = new G4VisAttributes();
  World_VisAtt->SetForceWireframe(true);  
  // logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  logicWorld->SetVisAttributes(World_VisAtt);

  VisAttWrapping = new G4VisAttributes(true, G4Colour(0.5, 0.5, 0.5, 0.4));
  VisAttWrapping->SetForceWireframe(false);
  logicWrapping->SetVisAttributes(VisAttWrapping);

  VisAttScint = new G4VisAttributes(true, G4Colour::Blue());
  VisAttScint->SetForceWireframe(false);
  logicScint->SetVisAttributes(VisAttScint);

  // VisAttSensor = new G4VisAttributes(true, G4Colour::White());  
  // VisAttSensor->SetForceWireframe(true);  
  // logicSensor->SetVisAttributes(VisAttSensor);

}
