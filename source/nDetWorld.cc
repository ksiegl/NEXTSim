
#include "nDetWorld.hh"
#include "nDetWorldObject.hh"
#include "nDetWorldMessenger.hh"
#include "nDetMaterials.hh"

#include "gdmlSolid.hh"
#include "optionHandler.hh" // split_str
#include "termColors.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "CADMesh.hh"
//CERN
#include "Tape.hh"
#include "IS530_Chamber.hh"
#include "IS530_Plastic.hh"
#include "Polyhedron.hh"
#include "CloverSingleDetector.hh"
#include "CloverQuadDetector.hh"
#include "CloverSingleBuchDetector.hh"
#include "CloverQuadBuchDetector.hh"
#include "global.hh"

#include "CERNFrame.hh"
#include "CERNFloor.hh"
#include "CERNTapeBox.hh"
#include "CERNSupport.hh"

#define DEFAULT_FLOOR_MATERIAL "G4_CONCRETE"

nDetWorld::nDetWorld() : solidV(NULL), logV(NULL), physV(NULL), fillMaterial("air"), floorMaterial(), floorThickness(0), floorSurfaceY(0) {
	// Set the default size of the experimental hall
	hallSize = G4ThreeVector(10*m, 10*m, 10*m);
	
	messenger = new nDetWorldMessenger(this);
}
	
bool nDetWorld::setWorldFloor(const G4String &input){
	// Expects a space-delimited string of the form:
	//  "centerY(cm) thickness(cm) [material=G4_CONCRETE]"
	std::vector<std::string> args;
	unsigned int Nargs = split_str(input, args);
	if(Nargs < 2){
		std::cout << " nDetConstruction: Invalid number of arguments given to ::SetWorldFloor(). Expected 2, received " << Nargs << ".\n";
		std::cout << " nDetConstruction:  SYNTAX: <centerY> <thickness> [material=G4_CONCRETE]\n";
		return false;
	}
	floorSurfaceY = strtod(args.at(0).c_str(), NULL)*cm;
	floorThickness = strtod(args.at(1).c_str(), NULL)*cm;
	if(Nargs < 3) // Defaults to concrete
		floorMaterial = DEFAULT_FLOOR_MATERIAL;
	else
		floorMaterial = args.at(2);
	return true;
}

void nDetWorld::buildExpHall(nDetMaterials *materials){
	solidV = new G4Box("solidV", hallSize.getX()/2, hallSize.getY()/2, hallSize.getZ()/2);

	G4Material *expHallFill = materials->getMaterial(fillMaterial);
	if(!expHallFill){ // Use the default material, if
		std::cout << Display::WarningStr("nDetWorld") << "Failed to find user-specified world material (" << fillMaterial << ")!" << Display::ResetStr() << std::endl;
		std::cout << Display::WarningStr("nDetWorld") << " Defaulting to filling world volume with air" << Display::ResetStr() << std::endl;
		expHallFill = materials->fAir;
	}

	logV = new G4LogicalVolume(solidV, expHallFill, "expHallLogV", 0, 0, 0);
	logV->SetVisAttributes(G4VisAttributes::Invisible);

	// Add a floor to the experimental hall (disabled by default)
	if(!floorMaterial.empty() && floorThickness > 0){
		G4Material *mat = materials->getMaterial(floorMaterial);
		if(mat){
			G4Box *floorBox = new G4Box("floor", hallSize.getX()/2, floorThickness/2, hallSize.getZ()/2);
			G4LogicalVolume *floor_logV;
			if(floorPitSize.getX() > 0 && floorPitSize.getX() > 0 && floorPitSize.getX() > 0){ // Dig a pit
				G4double pitCenterOffsetY = 0.5*(floorThickness - floorPitSize.getY());
				G4Box *pitBox = new G4Box("pitBox", floorPitSize.getX()/2, floorPitSize.getY()/2, floorPitSize.getZ()/2);
				G4SubtractionSolid *floorWithPit = new G4SubtractionSolid("floorWithPit", floorBox, pitBox, NULL, G4ThreeVector(0, pitCenterOffsetY, 0));
				floor_logV = new G4LogicalVolume(floorWithPit, mat, "floor_logV");
			}
			else{
				floor_logV = new G4LogicalVolume(floorBox, mat, "floor_logV");
			}
			floor_logV->SetVisAttributes(materials->visShadow);
			logV->SetVisAttributes(materials->visAssembly);
			new G4PVPlacement(NULL, G4ThreeVector(0, -(floorSurfaceY+floorThickness/2), 0), floor_logV, "floorBox_physV", logV, 0, 0, false);
		}
		else{
			std::cout << Display::WarningStr("nDetWorld") << "Failed to find user-specified floor material (" << floorMaterial << ")!" << Display::ResetStr() << std::endl;
			std::cout << Display::WarningStr("nDetWorld") << " Disabling the use of a floor" << Display::ResetStr() << std::endl;
			floorMaterial = "";
			floorThickness = 0;
		}
	}

	// Add additional objects
	for(auto obj : objects){
		obj->placeObject(logV, materials);
	}
	
	if(!expName.empty() && expName!="isolde" && expName!="RIKEN" && expName!="ORNL2016" && expName!="Argonne" && expName!="FDSi" && expName!="e14060") 
		cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n" <<
				"<<<<<<<<<<<<<<<<<<<<<<<<< Unrecognizable expriment name. Please check for appropriate naming schemes. >>>>>>>>>>>>>>>>>>>>>>>>>\n" <<
				"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
			
	if(expName=="isolde") BuildCERNStructures();

	// Place the experimental hall into the world
	physV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logV, "expHallPhysV", 0, false, 0);
	
	if(expName=="isolde") BuildCERNElements();

	if(expName=="FDSi") BuildFDSi(materials);

	if(expName=="e14060") BuildE14060(materials);

	return;
}

void nDetWorld::BuildCERNStructures(){

   CERNFloor* cernFloor = new CERNFloor();
   G4RotationMatrix* floorRot = new G4RotationMatrix();
   G4double floorXPos = -126.5*cm;
   G4ThreeVector floorPosition = G4ThreeVector(floorXPos,0., 0.);
   cernFloor->Place(floorRot, floorPosition, "cernFloor", logV); 
   
   CERNFrame* cernFrame = new CERNFrame();
   G4RotationMatrix* rotFrame = new G4RotationMatrix();
   rotFrame->rotateX(0*degree);
   G4double frameXPos = -22*cm;
   G4double frameZPos = 25*cm;
   G4ThreeVector framePosition = G4ThreeVector(frameXPos, 0, frameZPos);
   cernFrame->Place(rotFrame, framePosition, "cernFrame", logV); 

   CERNTapeBox* tapeBox = new CERNTapeBox();
   G4RotationMatrix* rotTapeBox = new G4RotationMatrix();
   G4double tapeBoxXPos = -86.55*cm;
   G4ThreeVector tapeBoxPosition = G4ThreeVector(tapeBoxXPos, 0, 0);
   tapeBox->Place(rotTapeBox, tapeBoxPosition, "cernTapeBox", logV);                 
    
   CERNSupport* cernSupport = new CERNSupport();
   G4RotationMatrix* rotSupport = new G4RotationMatrix(); 
   G4ThreeVector supportPos(0.0,0.0, 5.7*cm);  
   cernSupport->Place(rotSupport, supportPos, "cernSupport", logV);        

   
	return;
}

void nDetWorld::BuildCERNElements(){
  vector<CloverQuadDetector*> 		clquad_array;
  vector<CloverQuadBuchDetector*> 	clquadbuch_array;
  vector<Tape*>						tape_array;
  vector<Polyhedron*>				poly_array;
  //vector<Cubic_Chamber*>			cubic_chamber_array;
  vector<IS530_Chamber*>			IS530_chamber_array;
  vector<IS530_Plastic*>			IS530_plastic_array;

	///setup copied from Ola's config file
    /* 1	81		 45		-35.26		15		KU-Blue-103
       2	75		-45		-35.26		-15		Buch-Red-99
       1	82		-45		 35.26		15		Ku-Yellow-102
       2	72		 45		 35.26	 	-15		Buch-Green-101
       #2	60		180		0		0		Buch-Extra
       8	0		0		0		0	        Tape
       #10	0		0		0		0		OSIRIS
       #9	0		0		0		0		TPiece-Chamber
       12	0		0		0		0		IS530-Chamber
       13	0		0		0		0		IS530-Plastic
       #14	49		180		0		0		VetoPlastic */

	G4int    gType[8]     = {1,2,1,2,8,12,13,10};
	G4double gDistance[8] = {81,75,82,72,0,0,0,0};
	G4double gTheta[8]    = {45,-45,-45,45,0,0,0,0};
	G4double gPhi[8]      = {-35.26,-35.26,35.26,35.26,0,0,0,0};
	G4double gSpin[8]     = {15,-15,15,-15,0,0,0,0};
	G4int gLines = 8;

	for(int l=0;l<gLines;l++){
		
		if(1==gType[l]){   									// Clover KU
			clquad_array.push_back(new CloverQuadDetector(	physV,(G4double) gDistance[l]*mm,(G4double) (gTheta[l]*deg),(G4double) (gPhi[l]*deg),(G4double) (gSpin[l]*deg),clquad_array.size()));
			
		}
		if(2==gType[l]){  								// Clover Buch
			clquadbuch_array.push_back(	new CloverQuadBuchDetector(	physV, (G4double) gDistance[l]*mm, (G4double) (gTheta[l]*deg), (G4double) (gPhi[l]*deg), (G4double) (gSpin[l]*deg), clquadbuch_array.size()));
		}

		if(8==gType[l]){  								// Tape
			tape_array.push_back(new Tape(physV));
			if(1<tape_array.size()){
				cout<<"\nYou can only contruct the tape one time!!!\n";
				exit(0);
			}
		}
		if(10==gType[l]){  								// Polyhedron 
			poly_array.push_back(new Polyhedron(physV, (G4double) gDistance[l]*mm, (G4double) (gTheta[l]*deg), (G4double) (gPhi[l]*deg), (G4double) (gSpin[l]*deg)));
			if(1<poly_array.size()){
				cout<<"\nYou can only contruct the Polyhedron frame one time!!!\n";
				exit(0);
			}
		}		
		if(12==gType[l]){  				// IS530 Chamber 
			IS530_chamber_array.push_back(new IS530_Chamber(physV, (G4double) gDistance[l]*mm, (G4double) (gTheta[l]*deg), (G4double) (gPhi[l]*deg), (G4double) (gSpin[l]*deg)));
			if(1<IS530_chamber_array.size()){
				cout<<"\nYou can only contruct the IS530 chamber frame one time!!!\n"; 
				exit(0);
			}
			
		}
		if(13==gType[l]){  							// IS530 Plastic - non-sensitive detector 
			IS530_plastic_array.push_back(new IS530_Plastic(physV, (G4double) gDistance[l]*mm, (G4double) (gTheta[l]*deg), (G4double) (gPhi[l]*deg), (G4double) (gSpin[l]*deg)));
			if(1<IS530_plastic_array.size()){
				cout<<"\nYou can only contruct the IS530 Plastic one time!!!\n"; 
				exit(0);
			}
		}
			
	}

  //Construction

  // 1. Clover KU Leuven
  for (int clq=0; clq<clquad_array.size(); clq++){
    clquad_array.at(clq)->Construct();
  }
  // 2. Clover Bucharest
  for (int clq=0; clq<clquadbuch_array.size(); clq++){
    clquadbuch_array.at(clq)->Construct();
  }

  // 8. Tape
  for (int t=0; t<tape_array.size(); t++){
    tape_array.at(t)->Construct();
  }

  // 10. OSIRIS chamber
  for (int pl=0; pl<poly_array.size(); pl++){
    poly_array.at(pl)->Construct();
  }
  /*
  // 11. Cubic chamber
  for (int cc=0; cc<cubic_chamber_array.size(); cc++){
    cubic_chamber_array.at(cc)->Construct();
  }*/
  // 12. IS530 chamber
  for (int is=0; is<IS530_chamber_array.size(); is++){
    IS530_chamber_array.at(is)->Construct();
  }
  // 13. IS530 plastic - non sensitive detector
  for (int is=0; is<IS530_plastic_array.size(); is++){
    IS530_plastic_array.at(is)->Construct();
  }
 
  cout << "CERN setup - DONE" <<endl;
}

void nDetWorld::BuildFDSi(nDetMaterials *materials){ 
	G4RotationMatrix rotation = G4RotationMatrix();
	rotation.rotateX(90*degree);
	//CADMesh* Ball_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/FDS/stl/FDSI_Clarionet_VANDLE_Ball_Shell.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* VANDLEFrame_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/FDS/stl/FDSI_Clarionet_VANDLE_moreLow.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	//G4VSolid* Ball_sol = Ball_mesh->TessellatedMesh();
	G4VSolid* VANDLEFrame_sol = VANDLEFrame_mesh->TessellatedMesh();
	G4Material* Frame_mat = materials->fAluminum;
	//G4LogicalVolume* Ball_log= new G4LogicalVolume(Ball_sol, Frame_mat,"Ball log");
	G4LogicalVolume* VANDLEFrame_log= new G4LogicalVolume(VANDLEFrame_sol, Frame_mat,"Ball log");
    	G4VisAttributes* Al_vis_att 	= new G4VisAttributes(G4Colour(0.65,0.65,0.65));
	//Ball_log->SetVisAttributes(Al_vis_att);
	VANDLEFrame_log->SetVisAttributes(Al_vis_att);
	G4Transform3D transformation(rotation, G4ThreeVector(0, 0, 0));
	//G4VPhysicalVolume* Ball_phys = new G4PVPlacement(transformation, Ball_log, "Ball_PhysV", logV, false, 0);
	G4VPhysicalVolume* VANDLEFrame_phys = new G4PVPlacement(transformation, VANDLEFrame_log, "VANDLEFrame_PhysV", logV, false, 0);
}

void nDetWorld::BuildE14060(nDetMaterials *materials){ 
	G4RotationMatrix frame_rotation = G4RotationMatrix();
	frame_rotation.rotateX(90*degree);
	G4RotationMatrix stack_rotation = G4RotationMatrix();
	stack_rotation.rotateZ(-90*degree);

	CADMesh* VANDLEFrame_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_full_frame_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* Foam_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_double_styrofoam_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* HAGFrame_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_double_HAGRID_mount_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* HAGShell_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_HAGRiD_Crystal_shells_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* HAGCrystal_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_HAGRiD_Crystal_labr_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* HAGPMT_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_HAGRiD_PMTs_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	//CADMesh* NaIShell_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_NaI_Shell_Stack_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* NaICrystal_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_NaI_Crystal_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	
	CADMesh* ImplantFrame_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/Implant_Frame_Al_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* ImplantYAP_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/E14060_YAP_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* ImplantLG_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/Quartz_LightGuide_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* ImplantVeto_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/VETO_Plastic_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* ImplantAngr_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/Angr_board_round_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	CADMesh* ImplantPLA_mesh = new CADMesh(const_cast<char*>("/SCRATCH/DScratch3/e14060/stl/Implant_Frame_Pla_mm_high.stl"),       mm,  G4ThreeVector(0*cm, 0*cm, 0*cm), false);
	
	
	G4VSolid* VANDLEFrame_sol = VANDLEFrame_mesh->TessellatedMesh();
	G4VSolid* Foam_sol = Foam_mesh->TessellatedMesh();
	G4VSolid* HAGFrame_sol = HAGFrame_mesh->TessellatedMesh();
	G4VSolid* HAGShell_sol = HAGShell_mesh->TessellatedMesh();
	G4VSolid* HAGCrystal_sol = HAGCrystal_mesh->TessellatedMesh();
	G4VSolid* HAGPMT_sol = HAGPMT_mesh->TessellatedMesh();
	//G4VSolid* NaIShell_sol = NaIShell_mesh->TessellatedMesh();
	G4VSolid* NaICrystal_sol = NaICrystal_mesh->TessellatedMesh();

	G4VSolid* ImplantFrame_sol = ImplantFrame_mesh->TessellatedMesh();
	G4VSolid* ImplantYAP_sol = ImplantYAP_mesh->TessellatedMesh();
	G4VSolid* ImplantVeto_sol = ImplantVeto_mesh->TessellatedMesh();
	G4VSolid* ImplantAngr_sol = ImplantAngr_mesh->TessellatedMesh();
	G4VSolid* ImplantPLA_sol = ImplantPLA_mesh->TessellatedMesh();
	G4VSolid* ImplantLG_sol = ImplantLG_mesh->TessellatedMesh();
	


	G4Material* Frame_mat = materials->fAluminum;
	G4Material* Foam_mat = materials->fPLA;
	G4Material* HAG_mat = materials->fLaBr3;
	G4Material* PMT_mat = materials->fMuMetal;
	G4Material* NaIShell_mat = materials->fSteel;
	G4Material* NaI_mat = materials->fNaI;
	G4Material* YAP_mat = materials->fYAP;
	G4Material* LG_mat = materials->fSiO2;
	G4Material* Veto_mat = materials->fEJ200;
	
	G4LogicalVolume* VANDLEFrame_log= new G4LogicalVolume(VANDLEFrame_sol, Frame_mat,"Frame log");
	G4LogicalVolume* Foam_log= new G4LogicalVolume(Foam_sol, Foam_mat,"Foam log");
	G4LogicalVolume* HAGFrame_log= new G4LogicalVolume(HAGFrame_sol, Frame_mat,"HAGFrame log");
	G4LogicalVolume* HAGShell_log= new G4LogicalVolume(HAGShell_sol, Frame_mat,"HAGShell log");
	G4LogicalVolume* HAGCrystal_log= new G4LogicalVolume(HAGCrystal_sol, HAG_mat,"HAGCrystal log");
	G4LogicalVolume* HAGPMT_log= new G4LogicalVolume(HAGPMT_sol, PMT_mat,"HAGPMT log");
	//G4LogicalVolume* NaIShell_log= new G4LogicalVolume(VANDLEFrame_sol, NaIShell_mat,"NaIShell log");
	G4LogicalVolume* NaICrystal_log= new G4LogicalVolume(NaICrystal_sol, NaI_mat,"NaI Crystal log");

	G4LogicalVolume* ImplantFrame_log= new G4LogicalVolume(ImplantFrame_sol, Frame_mat,"ImplantFrame log");
	G4LogicalVolume* ImplantYAP_log= new G4LogicalVolume(ImplantYAP_sol, YAP_mat,"ImplantYAP log");
	G4LogicalVolume* ImplantVeto_log= new G4LogicalVolume(ImplantVeto_sol, Veto_mat,"ImplantVeto log");
	G4LogicalVolume* ImplantAngr_log= new G4LogicalVolume(ImplantFrame_sol, NaIShell_mat,"ImplantAngr log");
	G4LogicalVolume* ImplantPLA_log= new G4LogicalVolume(ImplantPLA_sol, Foam_mat,"ImplantPLA log");
	G4LogicalVolume* ImplantLG_log= new G4LogicalVolume(ImplantFrame_sol, LG_mat,"ImplantLG log");


    G4VisAttributes* Al_vis_att 	= new G4VisAttributes(G4Colour(0.65,0.65,0.65));
    G4VisAttributes* Steel_vis_att 	= new G4VisAttributes(G4Colour(0.45,0.45,0.45));
	Steel_vis_att->SetForceSolid(true);
    G4VisAttributes* Det_vis_att 	= new G4VisAttributes(G4Colour(0.0,0.95,0.0));
    G4VisAttributes* YAP_vis_att 	= new G4VisAttributes(G4Colour(0.0,0,0,0.95));
	YAP_vis_att->SetForceSolid();
    G4VisAttributes* PMT_vis_att 	= new G4VisAttributes(G4Colour(0.15,0.15,0.15));
    G4VisAttributes* Foam_vis_att 	= new G4VisAttributes(G4Colour(0.85,0.85,0.85));
	
	VANDLEFrame_log->SetVisAttributes(Al_vis_att);
	Foam_log->SetVisAttributes(Foam_vis_att);
	HAGFrame_log->SetVisAttributes(Al_vis_att);
	HAGShell_log->SetVisAttributes(Al_vis_att);
	HAGCrystal_log->SetVisAttributes(Det_vis_att);
	//NaIShell_log->SetVisAttributes(Steel_vis_att);
	NaICrystal_log->SetVisAttributes(Det_vis_att);
	HAGPMT_log->SetVisAttributes(PMT_vis_att);

	ImplantFrame_log->SetVisAttributes(Al_vis_att);
	ImplantYAP_log->SetVisAttributes(YAP_vis_att);
	ImplantVeto_log->SetVisAttributes(Det_vis_att);
	ImplantAngr_log->SetVisAttributes(Steel_vis_att);
	ImplantPLA_log->SetVisAttributes(Foam_vis_att);
	ImplantLG_log->SetVisAttributes(Foam_vis_att);
	
	G4Transform3D frame_transformation(frame_rotation, G4ThreeVector(0,0,19*cm));
	
	G4VPhysicalVolume* VANDLEFrame_phys = new G4PVPlacement(frame_transformation, VANDLEFrame_log, "VANDLEFrame_PhysV", logV, false, 0);
	G4VPhysicalVolume* Foam_phys = new G4PVPlacement(frame_transformation, Foam_log, "Foam_PhysV", logV, false, 0);
	G4VPhysicalVolume* HAGFrame_phys = new G4PVPlacement(frame_transformation, HAGFrame_log, "HAGFrame_PhysV", logV, false, 0);
	G4VPhysicalVolume* HAGShell_phys = new G4PVPlacement(frame_transformation, HAGShell_log, "HAGShell_PhysV", logV, false, 0);
	G4VPhysicalVolume* HAGCrystal_phys = new G4PVPlacement(frame_transformation, HAGCrystal_log, "HAGCrystal_PhysV", logV, false, 0);
	G4VPhysicalVolume* HAGPMT_phys = new G4PVPlacement(frame_transformation, HAGPMT_log, "HAGPMT_PhysV", logV, false, 0);
	//G4VPhysicalVolume* NaIShell_phys = new G4PVPlacement(frame_transformation, NaIShell_log, "NaIShell_PhysV", logV, false, 0);
	G4VPhysicalVolume* NaICrystal_phys = new G4PVPlacement(frame_transformation, NaICrystal_log, "NaICrystal_PhysV", logV, false, 0);

	G4Transform3D stack_transformation(stack_rotation,G4ThreeVector(-2.3*cm,0,0));

	G4VPhysicalVolume* ImplantFrame_phys = new G4PVPlacement(stack_transformation, ImplantFrame_log, "ImplantFrame_PhysV", logV, false, 0);
	G4VPhysicalVolume* ImplantYAP_phys = new G4PVPlacement(stack_transformation, ImplantYAP_log, "ImplantYAP_PhysV", logV, false, 0);
	G4VPhysicalVolume* ImplantVeto_phys = new G4PVPlacement(stack_transformation, ImplantVeto_log, "ImplantVeto_PhysV", logV, false, 0);
	G4VPhysicalVolume* ImplantAngr_phys = new G4PVPlacement(stack_transformation, ImplantAngr_log, "ImplantAngr_PhysV", logV, false, 0);
	G4VPhysicalVolume* ImplantPLA_phys = new G4PVPlacement(stack_transformation, ImplantPLA_log, "ImplantPLA_PhysV", logV, false, 0);
	G4VPhysicalVolume* ImplantLG_phys = new G4PVPlacement(stack_transformation, ImplantLG_log, "ImplantLG_PhysV", logV, false, 0);

	ImplantYAP_phys->CheckOverlaps();
}

nDetWorldPrimitive *nDetWorld::addNewPrimitive(const G4String &str){
	nDetWorldPrimitive *obj = new nDetWorldPrimitive(str);
	obj->decodeString();
	if(!obj->decodeString()){
		std::cout << " nDetWorld: Invalid number of arguments given to ::addNewPrimitive(). Expected at least " << obj->getNumRequiredArgs() << " but received " << obj->getNumSuppliedArgs() << ".\n";
		std::cout << " nDetWorld:  SYNTAX: " << obj->syntaxStr() << std::endl;
		delete obj;
		return NULL;
	}
	objects.push_back(obj);
	return obj;
}

void nDetWorld::listAllPrimitives(){
	nDetWorldPrimitive dummy("");
	dummy.listAllPrimitives();
}

void nDetWorld::printDefinedObjects(){
	if(!objects.empty()){
		for(auto obj : objects){
			std::cout << "***********************************************************\n";
			obj->print();
		}
		std::cout << "***********************************************************\n";
	}
	else
		std::cout << " nDetWorldPrimitive: No 3d primitives are currently defined\n";
}

void nDetWorld::reset(){
	for(auto obj : objects)
		delete obj;
	objects.clear();
}

gdmlObject *nDetWorld::loadGDML(const G4String &input){
	gdmlObject *obj = new gdmlObject(input);
	obj->decodeString();
	if(!obj->decodeString()){
		std::cout << " nDetWorld: Invalid number of arguments given to ::loadGDML(). Expected " << obj->getNumRequiredArgs() << " but received " << obj->getNumSuppliedArgs() << ".\n";
		std::cout << " nDetWorld:  SYNTAX: " << obj->syntaxStr() << std::endl;
		delete obj;
		return NULL;
	}
	objects.push_back(obj);
	return obj;
}
