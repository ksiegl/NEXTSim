/************************************************************************
 * \file CloverQuadDetector.hh
 *
 * \class CloverQuadDetector   
 * \brief Design of a HPGe detector crystal (Clover 40x50x70) (KU Leuven)
 *
 * \author Christophe Sotty : csotty@gmail.com
 *                            Christophe.Sotty@fys.kuleuven.be
 *
 *************************************************************************/
//R. Lica on 01.03.2016 - Fixed the G4Transform3D bug, removed det_env

#include "CloverSingleDetector.hh"
//#include "CloverSingleHit.hh"
#include "CloverQuadDetector.hh"
//#include "CloverSingleSD.hh"
#include "CADMesh.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include <TText.h>

using namespace std;



CloverQuadDetector::CloverQuadDetector(G4VPhysicalVolume* p_mother, 
                                       G4double p_rho, 
                                       G4double p_theta, 
                                       G4double p_phi, 
                                       G4double p_spin, 
                                       G4int p_cl_nb)
  :mother(p_mother), rho(p_rho), theta(p_theta), phi(p_phi), spin(p_spin), cl_nb(p_cl_nb)
{
	//G4cout	<<"\nClover Quad (KU Leuven) DETECTOR#######################\n"	<<flush
	//	<<"\nRho =\t\t"		<<rho/mm	<< " mm"	<<flush
	//	<<"\nTheta =\t\t"	<<theta/deg	<< " deg" 	<<flush
	//	<<"\nPhi =\t\t"		<<phi/deg	<< " deg"       <<flush
	//	<<"\nSpin =\t\t"	<<spin/deg	<< " deg"       <<flush
	//	<<"\nCloverNb =\t"	<<cl_nb 			<<G4endl;
}

CloverQuadDetector::~CloverQuadDetector()
{;}


G4VPhysicalVolume* CloverQuadDetector::Construct()
{
    // Naming
    //
    G4String name="/Clover"+cl_nb;
    
    // Material
    //
    G4NistManager* NISTmgr = G4NistManager::Instance();

    G4Material* Vacuum          = NISTmgr->FindOrBuildMaterial("G4_Galactic");
    G4Material* Air             = NISTmgr->FindOrBuildMaterial("G4_AIR");
    G4Material* Al_mat          = NISTmgr->FindOrBuildMaterial("G4_Al");
    G4Element* H_el  = new G4Element("Hydrogen_el"   ,  "H_el", 1 ,     1.01*g/mole);
    G4Element* C_el  = new G4Element("Carbone_el"    ,  "C_el", 6 ,    12.01*g/mole);
    G4Element* O_el  = new G4Element("Oxygene_el"    ,  "O_el", 8 ,    16.00*g/mole);
    G4Material* CEpoxy = new G4Material("CEpoxy",1.6*g/cm3,3);
    CEpoxy->AddElement(H_el, 11);
    CEpoxy->AddElement(O_el, 3);
    CEpoxy->AddElement(C_el, 16);

    // Visualisation
    // 
    vector<G4VisAttributes*> CrysVisAttM;
    CrysVisAttM.push_back(new G4VisAttributes(G4Colour(0.0,0.0,1.0)));
    CrysVisAttM.push_back(new G4VisAttributes(G4Colour(0.0,1.0,0.0)));
    CrysVisAttM.push_back(new G4VisAttributes(G4Colour(1.0,0.0,0.0)));
    CrysVisAttM.push_back(new G4VisAttributes(G4Colour(0.0,1.0,1.0)));

    // Visualisation
    //
    Crystal0_vis_att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));   
    Crystal1_vis_att = new G4VisAttributes(G4Colour(0.0,1.0,0.0));   
    Crystal2_vis_att = new G4VisAttributes(G4Colour(1.0,0.0,0.0));   
    Crystal3_vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));   
    G4VisAttributes* Al_vis_att 	= new G4VisAttributes(G4Colour(0.65,0.65,0.65));
    G4VisAttributes* CEpoxy_vis_att 	= new G4VisAttributes(G4Colour(0.345,0.07059,0.07059));
    Capsule_vis_att       = Al_vis_att;    
    CapsuleCap_vis_att    = Al_vis_att;
    CarboneWindow_vis_att = CEpoxy_vis_att;
	
	

    // Detector Construction
    //
    // Meshing
 
    mesh_Capsule	= new CADMesh(const_cast<char*>("/ARCHIVE/Ddata/geant4_stl/vandle/isolde/Clover_KULeuven_RefFace/Clover_Assembly_RefModif_capsule2_1.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false); // back to 08122015 
    mesh_CapsuleCap	= new CADMesh(const_cast<char*>("/ARCHIVE/Ddata/geant4_stl/vandle/isolde/Clover_KULeuven_RefFace/Clover_Assembly_RefModif_Capsule_cap_2.stl"),    mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false); // back to 08122015 
    mesh_CarboneWindow	= new CADMesh(const_cast<char*>("/ARCHIVE/Ddata/geant4_stl/vandle/isolde/Clover_KULeuven_RefFace/Clover_Assembly_RefModif_Carbone_window_3.stl"), mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false); // back to 08122015 

/*	
    mesh_Capsule	= new CADMesh(const_cast<char*>("/home/ksiegl/NEXTSim/stl/Clover_Assembly_RefModif_capsule2_1.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false); // back to 08122015 
    mesh_CapsuleCap	= new CADMesh(const_cast<char*>("/home/ksiegl/NEXTSim/stl/Clover_Assembly_RefModif_Capsule_cap_2.stl"),    mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false); // back to 08122015 
    mesh_CarboneWindow	= new CADMesh(const_cast<char*>("/home/ksiegl/NEXTSim/stl/Clover_Assembly_RefModif_Carbone_window_3.stl"), mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false); // back to 08122015 
*/
	  cout<<"./stl/Clover_Assembly_RefModif_capsule2_1.stl"<<endl;
    Capsule_sol		= mesh_Capsule->TessellatedMesh();
    cout<<"./stl/Clover_Assembly_RefModif_Capsule_cap_2.stl"<<endl;
    CapsuleCap_sol	= mesh_CapsuleCap->TessellatedMesh();
    cout<<"./stl/Clover_Assembly_RefModif_Carbone_window_3.stl"<<endl;
    CarboneWindow_sol	= mesh_CarboneWindow->TessellatedMesh();

    // Logical Volume
    Capsule_log		= new G4LogicalVolume(Capsule_sol       , Al_mat, name+"/Capsule_log"           );
    CapsuleCap_log	= new G4LogicalVolume(CapsuleCap_sol    , Al_mat, name+"/CapsuleCap_log"        );
    CarboneWindow_log	= new G4LogicalVolume(CarboneWindow_sol , CEpoxy, name+"/CarboneWindow_log"     );
    Capsule_log		->SetVisAttributes(Capsule_vis_att);
    CapsuleCap_log	->SetVisAttributes(CapsuleCap_vis_att);
    CarboneWindow_log	->SetVisAttributes(CarboneWindow_vis_att);
    
    // Position and Orientation
    rho = -rho ;

    G4RotationMatrix rotation = G4RotationMatrix();
 
      if(0 != spin) rotation.rotateZ(spin);
      rotation.rotateY(phi);
      rotation.rotateX(theta);   

    
      G4ThreeVector translation(0., 0., rho);
      if(0 != spin) translation.rotateZ(spin);
      translation.rotateY(phi);
      translation.rotateX(theta);   
      
      G4Transform3D transformation(rotation, translation);
	
	
	 // Capsule, CapsuleCap, CarboneWindow
    Capsule_phys       = new G4PVPlacement(transformation, name+"/Capsule"      , Capsule_log      , mother, false, 0);
    CapsuleCap_phys    = new G4PVPlacement(transformation, name+"/CapsuleCap"   , CapsuleCap_log   , mother, false, 0);
    CarboneWindow_phys = new G4PVPlacement(transformation, name+"/CarboneWindow", CarboneWindow_log, mother, false, 0);
	
	 // Crystals
    Crystal0	= new CloverSingleDetector(mother, transformation, cl_nb, 0, Crystal0_vis_att);
    Crystal1	= new CloverSingleDetector(mother, transformation, cl_nb, 1, Crystal1_vis_att);
    Crystal2	= new CloverSingleDetector(mother, transformation, cl_nb, 2, Crystal2_vis_att);
    Crystal3	= new CloverSingleDetector(mother, transformation, cl_nb, 3, Crystal3_vis_att);
    
    
    G4VPhysicalVolume* cl0 = Crystal0->Construct();
    G4VPhysicalVolume* cl1 = Crystal1->Construct();
    G4VPhysicalVolume* cl2 = Crystal2->Construct();
    G4VPhysicalVolume* cl3 = Crystal3->Construct();



    return(mother);

}






















