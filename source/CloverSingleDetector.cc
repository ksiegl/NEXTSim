/************************************************************************
 * \file CloverSingleDetector.hh
 *
 * \class CloverSingleDetector   
 * \brief Design of a HPGe detector crystal (Clover 40x50x70)
 * 
 * \author Christophe Sotty : csotty@gmail.com
 *                            Christophe.Sotty@fys.kuleuven.be
 * 
 *************************************************************************/
//R. Lica - removed det_env


#include "CloverSingleDetector.hh"
//#include "CloverSingleHit.hh"
//#include "CloverSingleSD.hh"
#include "CADMesh.hh"
#include <TText.h>
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

using namespace std;



CloverSingleDetector::CloverSingleDetector(G4VPhysicalVolume* p_mother, G4Transform3D p_transformation, G4int p_cl_nb, G4int p_cr_nb, G4VisAttributes* pvis_att)
  :mother(p_mother), cl_nb(p_cl_nb), cr_nb(p_cr_nb), transformation(p_transformation), det_vis_att(pvis_att)
{
	
}
CloverSingleDetector::~CloverSingleDetector() 
{; }



G4VPhysicalVolume* CloverSingleDetector::Construct()
{
    // Material
    //
    G4NistManager* NISTmgr = G4NistManager::Instance();

    G4Material* Vacuum          = NISTmgr->FindOrBuildMaterial("G4_Galactic");
    G4Material* HPGe   		= NISTmgr->FindOrBuildMaterial("G4_Ge"); 

    // Detector Construction
    //
    // Meshing
    // 18022015 mesh_CloverSingle =  new CADMesh(const_cast<char*>(Form("../STL_export/Clover_KULeuven/Clover_Assembly_Crystal2_%i.stl",cr_nb+4)), mm,  G4ThreeVector(-0.2524*cm + 75*mm), false);

//    mesh_CloverSingle =  new CADMesh(const_cast<char*>(Form("/ARCHIVE/Ddata/geant4_stl/vandle/isolde/Clover_KULeuven_RefFace/Clover_Assembly_RefModif_Crystal2_%i.stl",cr_nb+4)), mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false); // back to 08122015
  
    mesh_CloverSingle =  new CADMesh(const_cast<char*>(Form("/home/ksiegl/NEXTSim/stl/Clover_Assembly_RefModif_Crystal2_%i.stl",cr_nb+4)), mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false); // back to 08122015

  
    //G4cout << Form("/home/jheidema/opt/NEXTSim/STL_export/Clover_KULeuven_RefFace/Clover_Assembly_RefModif_Crystal2_%i.stl",cr_nb+4) << G4endl; 
    CloverSingle_sol = mesh_CloverSingle->TessellatedMesh();
    CloverSingle_log = new G4LogicalVolume(CloverSingle_sol, HPGe, Form("/Clover%i_Crystal%i_log",cl_nb ,cr_nb));
    CloverSingle_log -> SetVisAttributes(det_vis_att);

   
    CloverSingle_phys = new G4PVPlacement(transformation, Form("/Clover%iCrystal%i_phys", cl_nb, cr_nb), CloverSingle_log, mother, false, 0);

    //G4SDManager* SDman = G4SDManager::GetSDMpointer();
    //cloverSingleSD = new CloverSingleSD(Form("/Clover%i_Crystal%i_SD", cl_nb, cr_nb));
    //SDman->AddNewDetector(cloverSingleSD);
    //CloverSingle_log->SetSensitiveDetector(cloverSingleSD);

//    G4cout << "CloverSingle Sensitive Detector Name = " << cloverSingleSD->GetName() << flush << G4endl;
    
    return(CloverSingle_phys);
}







































































