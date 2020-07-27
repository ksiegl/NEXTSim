#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include "termColors.hh"
#include "nDetMaterials.hh"
#include "nDetMaterialsMessenger.hh"
#include "nDetDynamicMaterial.hh"

nDetMaterials::nDetMaterials() : messenger(), isInitialized(false), scintsAreDefined(false), lightYieldScale(1) { 
	messenger = new nDetMaterialsMessenger(this);
}

nDetMaterials::~nDetMaterials(){
	if(isInitialized){
		// Optical surfaces
		delete fTeflonOpSurf;
		delete fSiliconOpSurf;
		delete fPerfectOpSurf;
		delete fGreaseOpSurf;	
		delete fMylarOpSurf;
		delete fEsrOpSurf;
		
		// Materials
		delete fGrease;
		
		// Material properties tables
		delete fPerfectMPT;
		delete fGreaseMPT;
		delete fEsrMPT;

		// Vis attributes
		delete visAssembly;
		delete visSensitive;
		delete visWindow;
		delete visGrease;
		delete visWrapping;
		delete visScint;
		delete visShadow;
	}
	if(scintsAreDefined){
		delete fEJ200;
		delete fEJ276;
		delete fEJ200MPT;
		delete fEJ276MPT;
	}
	delete messenger;
}

void nDetMaterials::initialize(){
	if(isInitialized) return;
	
	defineMaterials();
	defineScintillators();

	elementList["H"] = fH;
	elementList["C"] = fC;
	elementList["O"] = fO;
	elementList["F"] = fF;
	elementList["Si"] = fSi;
	elementList["Y"] = fAl;
	elementList["La"] = fLa;
	elementList["Br"] = fBr;
	elementList["Na"] = fNa;
	elementList["I"] = fI;
	elementList["Cr"] = fCr;
	elementList["Fe"] = fFe;
	elementList["Ni"] = fNi;
	
	materialList["air"] = fAir;
	materialList["vacuum"] = fVacuum;
	materialList["teflon"] = fTeflon;
    materialList["ej200"] = fEJ200;
    materialList["ej276"] = fEJ276; 
	materialList["grease"] = fGrease;
	materialList["quartz"] = fSiO2;
	materialList["silicon"] = fSilicon;
	materialList["mylar"] = fMylar;
	materialList["acrylic"] = fAcrylic;
	materialList["aluminum"] = fAluminum;
	materialList["labr3"] = fLaBr3;
	materialList["nai"] = fNaI;
	materialList["yap"] = fYAP;
	materialList["yso"] = fYSO;
	materialList["steel"] = fSteel;
	materialList["mumetal"] = fMuMetal;
	materialList["hdpe"] = fHDPE;

	opticalSurfaceList["teflon"] = fTeflonOpSurf;
	opticalSurfaceList["silicon"] = fSiliconOpSurf;
	opticalSurfaceList["mylar"] = fMylarOpSurf;
	opticalSurfaceList["esr"] = fEsrOpSurf;
	opticalSurfaceList["perfect"] = fPerfectOpSurf;
	opticalSurfaceList["grease"] = fGreaseOpSurf;

	//visAttributesList["air"] = &G4VisAttributes::Invisible;
	//visAttributesList["vacuum"] = &G4VisAttributes::Invisible;
	visAttributesList["teflon"] = visWrapping;
	visAttributesList["ej200"] = visScint;
	visAttributesList["ej276"] = visScint;
	visAttributesList["grease"] = visGrease;
	visAttributesList["quartz"] = visWindow;
	visAttributesList["silicon"] = visSensitive;
	visAttributesList["mylar"] = visWrapping;
	visAttributesList["acrylic"] = visWindow;
	visAttributesList["aluminum"] = visShadow;
}

G4Material* nDetMaterials::getUserDetectorMaterial(const G4String &name){
	G4Material *mat = getMaterial(name);
	
	if(!mat) // Material was not found
		std::cout << Display::ErrorStr("nDetDynamicMaterial") << "Detector material named \"" << name << "\" was not found in list!\n" << Display::ResetStr();
		
	return mat;
}

G4Material* nDetMaterials::getUserSurfaceMaterial(const G4String &name){
	G4Material *mat = NULL;
	if(name == "mylar" || name == "esr" || name == "perfect")
		mat = fMylar;
	else
		mat = getMaterial(name);

	if(!mat) // Material was not found
		std::cout << Display::ErrorStr("nDetDynamicMaterial") << "Wrapping material named \"" << name << "\" was not found in list!\n" << Display::ResetStr();

	return mat;
}

G4OpticalSurface* nDetMaterials::getUserOpticalSurface(const G4String &name){
	G4OpticalSurface *opt = getOpticalSurface(name);

	if(!opt) // Surface was not found
		std::cout << Display::ErrorStr("nDetDynamicMaterial") << "Optical surface named \"" << name << "\" was not found in list!\n" << Display::ResetStr();

    return opt; // default
}

G4VisAttributes* nDetMaterials::getUserVisAttributes(const G4String &name){
	G4VisAttributes *visatt = NULL;
	if(name == "mylar" || name == "teflon" || name == "esr" || name == "perfect")
		visatt = visWrapping;
	else
		visatt = getVisualAttributes(name);

	if(!visatt){ // Visual attributes was not found (use a default)
		std::cout << Display::WarningStr("nDetDynamicMaterial") << "Visual attributes named \"" << name << "\" was not found in list!\n" << Display::ResetStr();
		visatt = visAssembly;
	}
		
	return visatt;
}

void nDetMaterials::setLightYield(const G4double &yield){ 
	lightYieldScale = yield; 
	
	if(scintsAreDefined) // Re-define the scintillators as the light yield has changed
		defineScintillators();
}

bool nDetMaterials::searchForElement(const G4String &name){
	size_t dummy;
	return nist.searchElementList(name, dummy);
}

G4Element* nDetMaterials::getElement(const G4String &name){
	std::map<G4String, G4Element*>::iterator iter = elementList.find(name);
	if(iter != elementList.end()) // Found element in the dictionary
		return iter->second;
		
	// No element found in the dictionary, search for it in the NIST database
	G4Element *ptr = nist.searchForElement(name);
	if(ptr){ // Found element in the database, add it to the dictionary
		elementList[name] = ptr;
		return ptr;
	}
	
	return NULL;
}

bool nDetMaterials::searchForMaterial(const G4String &name){
	size_t dummy;
	return nist.searchMaterialList(name, dummy);
}

G4Material* nDetMaterials::getMaterial(const G4String &name){
	std::map<G4String, G4Material*>::iterator iter = materialList.find(name);
	if(iter != materialList.end()) // Found material in the dictionary
		return iter->second;
	// No material found in the dictionary, search for it in the NIST database
	G4Material *ptr = nist.searchForMaterial(name);
	if(ptr){ // Found material in the database, add it to the dictionary
		materialList[name] = ptr;
		return ptr;
	}
	
	return NULL;
}

G4OpticalSurface* nDetMaterials::getOpticalSurface(const G4String &name){
	std::map<G4String, G4OpticalSurface*>::iterator iter = opticalSurfaceList.find(name);
	if(iter != opticalSurfaceList.end()) // Found optical surface in the dictionary
		return iter->second;
	return NULL;
}
	
G4VisAttributes* nDetMaterials::getVisualAttributes(const G4String &name){
	std::map<G4String, G4VisAttributes*>::iterator iter = visAttributesList.find(name);
	if(iter != visAttributesList.end()) // Found optical surface in the dictionary
		return iter->second;
	return NULL;
}

void nDetMaterials::listMaterials() const {
	std::cout << "/////////////////////////////\n";
	std::cout << "// Defined Materials\n";
	std::cout << "/////////////////////////////\n\n";
	for(auto material : materialList){
		std::cout << "  " << material.first << std::endl;
	}
	std::cout << std::endl;
}

void nDetMaterials::listVisAttributes() const {
	std::cout << "/////////////////////////////\n";
	std::cout << "// Defined Visual Attributes\n";
	std::cout << "/////////////////////////////\n\n";
	for(auto visatt : visAttributesList){
		std::cout << "  " << visatt.first << std::endl;
	}
	std::cout << std::endl;
}

void nDetMaterials::listOptSurfaces() const {
	std::cout << "/////////////////////////////\n";
	std::cout << "// Defined Optical Surfaces\n";
	std::cout << "/////////////////////////////\n\n";
	for(auto optsurf : opticalSurfaceList){
		std::cout << "  " << optsurf.first << std::endl;
	}
	std::cout << std::endl;
}

bool nDetMaterials::buildNewMaterial(const G4String &filename){
	nDetDynamicMaterial* dmat = new nDetDynamicMaterial();
	dmat->setScalingFactor(lightYieldScale);
	if(!dmat->read(filename, this)) // Read the material file
		return false;
	std::string name = dmat->getFilePrefix();
	int nameCounter = 2;
	while(materialList.find(name) != materialList.end()){ // Check for material existing in material list
		std::stringstream stream;
		stream << name << "-" << nameCounter++;
		
		name = stream.str();
	}
	if(nameCounter != 2){ // Print a warning about existing material name
		std::cout << Display::WarningStr("nDetDynamicMaterial") << "Material named \"" << dmat->getFilePrefix() << "\" is already in material list!\n" << Display::ResetStr();
		std::cout << Display::WarningStr("nDetDynamicMaterial") << " Renaming new material to \"" << name << "\"\n" << Display::ResetStr();
	}
	materialList[name] = dmat->getMaterial();
	visAttributesList[name] = dmat->getVisAttributes();
	opticalSurfaceList[name] = dmat->getOpticalSurface();
	return true;
}

void nDetMaterials::printMaterial(const G4String &name){
	G4Material *mat = getMaterial(name);
	if(mat){
		std::cout << "/////////////////////////////\n";
		std::cout << "// Material (" << name << ")\n";
		std::cout << "/////////////////////////////\n\n";
		std::cout << (*mat);
		// Re-writing G4MaterialPropertiesTable::DumpTable() because it prints to G4cout only
		// This is messy, but I have to do a lot of workarounds to deal with the strange setup of the G4MaterialPropertiesTable class
		std::cout << "\n/////////////////////////////\n";
		std::cout << "// Properties\n";
		std::cout << "/////////////////////////////\n\n";
		G4MaterialPropertiesTable* MPT = mat->GetMaterialPropertiesTable();
		if(MPT){
#ifndef GEANT_OLDER_VERSION
			std::vector<G4String> propNames = MPT->GetMaterialPropertyNames(); // Ugh...
			std::vector<G4String> cpropNames = MPT->GetMaterialConstPropertyNames(); // UGH... WHY, GEANT???
			for(auto prop : (*MPT->GetPropertyMap())){ // Iterate over variable properties
				G4PhysicsVector *vec = prop.second;
				G4String propertyName = propNames[prop.first];
				std::cout << std::string(propertyName.length(), '-') << std::endl;
				std::cout << propertyName << std::endl;
				std::cout << std::string(propertyName.length(), '-') << std::endl;
				if(vec != NULL){
					// Now I re-write G4PhysicsVector::DumpValues() for the same reason
					for (size_t i = 0; i < vec->GetVectorLength(); i++) // Iterate over all values in the vectors
						std::cout << vec->Energy(i) << "\t" << (*vec)[i] << std::endl;
				}
			}
			std::cout << "\n/////////////////////////////\n";
			std::cout << "// Constant Properties\n";
			std::cout << "/////////////////////////////\n\n";
			for(auto cprop : (*MPT->GetConstPropertyMap())) // Iterate over constant properties
				std::cout << cpropNames[cprop.first] << " = " << cprop.second << std::endl;
#else
			const std::map<G4String, G4MaterialPropertyVector*, std::less< G4String > >* pMap = MPT->GetPropertiesMap();
			const std::map<G4String, G4double, std::less< G4String > >* cpMap = MPT->GetPropertiesCMap();
			for(auto prop : (*pMap)){ // Iterate over variable properties
				G4PhysicsVector *vec = prop.second;
				std::cout << std::string(prop.first.length(), '-') << std::endl;
				std::cout << prop.first << std::endl;
				std::cout << std::string(prop.first.length(), '-') << std::endl;
				if(vec != NULL){
					// Now I re-write G4PhysicsVector::DumpValues() for the same reason
					for (size_t i = 0; i < vec->GetVectorLength(); i++) // Iterate over all values in the vectors
						std::cout << vec->Energy(i) << "\t" << (*vec)[i] << std::endl;
				}
			}
			std::cout << "\n/////////////////////////////\n";
			std::cout << "// Constant Properties\n";
			std::cout << "/////////////////////////////\n\n";
			for(auto cprop : (*cpMap)) // Iterate over constant properties
				std::cout << cprop.first << " = " << cprop.second << std::endl;
#endif
		}
		else{
			std::cout << " No properties table defined for this material\n";
		}
		std::cout << std::endl;
	}
}

void nDetMaterials::defineMaterials(){
	if(isInitialized) return;

	/////////////////////////////////////////////////////////////////
	// Visual Attributes
	/////////////////////////////////////////////////////////////////

	visAssembly = new G4VisAttributes();

	visSensitive = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75)); // grey
	visSensitive->SetForceSolid(true);

	visWindow = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0)); // cyan	

	visGrease = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red

	visWrapping = new G4VisAttributes();
	visWrapping->SetColor(1, 0, 1, 0.5); // Alpha=50%
	visWrapping->SetForceSolid(true);

	visScint = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue

	visShadow = new G4VisAttributes();
	visShadow->SetColor(0, 1, 0, 0.5); // green, Alpha=50%
	visShadow->SetForceSolid(true);	

	/////////////////////////////////////////////////////////////////
	// Elements & Molecular Components
	/////////////////////////////////////////////////////////////////

	// Elements
	fH = nist.searchForElement("H");
	fC = nist.searchForElement("C");
	fO = nist.searchForElement("O");
	fF = nist.searchForElement("F");
	fSi = nist.searchForElement("Si");
	fAl = nist.searchForElement("Al");
	fY = nist.searchForElement("Y");
	fLa = nist.searchForElement("La");
	fBr = nist.searchForElement("Br");
	fNa = nist.searchForElement("Na");
	fI = nist.searchForElement("I");
	fCr = nist.searchForElement("Cr");
	fFe = nist.searchForElement("Fe");
	fNi = nist.searchForElement("Ni");
	fCu = nist.searchForElement("Cu");

	// Air
    fAir = nist.searchForMaterial("G4_AIR");

	// Lab vacuum
	fVacuum = nist.searchForMaterial("G4_Galactic");

	//Sodium Iodide
	fNaI = nist.searchForMaterial("G4_SODIUM_IODIDE");
	
	//High-Density Polyethylene
	fHDPE = nist.searchForMaterial("G4_POLYETHYLENE");

	//Stainless Steel
	fSteel = nist.searchForMaterial("G4_STAINLESS-STEEL");

	G4int ncomp = 3;
    	fYSO =new G4Material("YSO", 4.5*g/cm3, ncomp); 
    	G4int nAtoms;
    	fYSO->AddElement(fY,  nAtoms=2);
    	fYSO->AddElement(fSi,  nAtoms=1);
   	 fYSO->AddElement(fO,  nAtoms=5);

	//YAP
	fYAP =new G4Material("YAP", 5.37*g/cm3, ncomp); 
    	fYAP->AddElement(fY,  nAtoms=1);
    	fYAP->AddElement(fAl,  nAtoms=1);
   	fYAP->AddElement(fO,  nAtoms=3);
    
    	// Lanthanum Bromide
    	G4int ncomp_labr3 = 2;
    	fLaBr3 = new G4Material("LaBr3", 5.06*g/cm3, ncomp_labr3);
    	fLaBr3->AddElement(fLa,  nAtoms=1);
    	fLaBr3->AddElement(fBr,  nAtoms=3);

    	// PLA
    	G4int ncomp_PLA = 3;
    	fPLA = new G4Material("PLA", 1.24*g/cm3, ncomp_PLA);
    	fPLA->AddElement(fH,  nAtoms = 1);
    	fPLA->AddElement(fO,  nAtoms = 3);
    	fPLA->AddElement(fC,  nAtoms = 3);
    
    	materialList["PLA"] = fPLA;

	// Mu-Metal
    	G4int ncomp_MM = 3;
    	fMuMetal = new G4Material("MuMetal", 1.24*g/cm3, ncomp_PLA);
    	fMuMetal->AddElement(fNi,  nAtoms = 77);
    	fMuMetal->AddElement(fFe,  nAtoms = 16);
    	fMuMetal->AddElement(fCu,  nAtoms = 5);

	/////////////////////////////////////////////////////////////////
	// Teflon (C2F4)n
	/////////////////////////////////////////////////////////////////

    fTeflon = nist.searchForMaterial("G4_TEFLON");

    G4double photonEnergy_teflon[9] = {1.607*eV, 1.743*eV, 1.908*eV, 2.108*eV, 2.354*eV, 2.664*eV, 3.070*eV, 3.621*eV, 4.413*eV};
    G4double reflectivity_teflon[9] = {0.514, 0.583, 0.656, 0.727, 0.789, 0.836, 0.868, 0.887, 0.892}; // https://www.osti.gov/servlets/purl/1184400 (1 layer)
    G4double efficiency_teflon[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    G4double absorption_teflon[9] =  {0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm};
    G4double refIndex_teflon[9] = {1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315};
    
    fTeflonMPT = new G4MaterialPropertiesTable();
    fTeflonMPT->AddProperty("REFLECTIVITY", photonEnergy_teflon, reflectivity_teflon, 9);
    fTeflonMPT->AddProperty("EFFICIENCY", photonEnergy_teflon, efficiency_teflon, 9);
    fTeflonMPT->AddProperty("RINDEX", photonEnergy_teflon, refIndex_teflon, 9);
    fTeflonMPT->AddProperty("ABSLENGTH", photonEnergy_teflon, absorption_teflon, 9);
    fTeflon->SetMaterialPropertiesTable(fTeflonMPT);

	/////////////////////////////////////////////////////////////////
	// Silicone Optical Grease (C2H6OSi)n
	/////////////////////////////////////////////////////////////////

    fGrease = new G4Material("Grease", 1.06*g/cm3, 4);

    fGrease->AddElement(fC, 2);
    fGrease->AddElement(fH, 6);
    fGrease->AddElement(fO, 1);
    fGrease->AddElement(fSi, 1);

    /*const G4int nEntries_Grease = 5;
    G4double Photon_Energy[nEntries_Grease] = { 2.757*eV, 3.102*eV, 3.312*eV, 3.545*eV, 4.136*eV };
    G4double RefractiveIndex_grease[nEntries_Grease] = {1.468, 1.473, 1.477, 1.482, 1.496};
    G4double Absorption_grease[nEntries_Grease] = { 195*mm,  195*mm, 195*mm, 195*mm, 195*mm };
    fGreaseMPT=new G4MaterialPropertiesTable();

    fGreaseMPT->AddProperty("RINDEX", Photon_Energy, RefractiveIndex_grease, nEntries_Grease);
    fGreaseMPT->AddProperty("ABSLENGTH", Photon_Energy, Absorption_grease, nEntries_Grease);*/

    G4double greasePhotonEnergy[3] = { 2*eV, 3*eV, 4*eV};
    G4double refIndexGrease[3] = {1.465, 1.465, 1.465};
    G4double absorptionGrease[3] = {195*mm, 195*mm, 195*mm};

    fGreaseMPT = new G4MaterialPropertiesTable();
    fGreaseMPT->AddProperty("RINDEX", greasePhotonEnergy, refIndexGrease, 3);
    fGreaseMPT->AddProperty("ABSLENGTH", greasePhotonEnergy, absorptionGrease, 3);
    fGrease->SetMaterialPropertiesTable(fGreaseMPT);

	/////////////////////////////////////////////////////////////////
	// Quartz (SiO2)
	/////////////////////////////////////////////////////////////////

    fSiO2 = nist.searchForMaterial("G4_SILICON_DIOXIDE");

    //optical properties of SiO2 - fused silica or fused quartz
    G4double PhotonEnergy[5] = { 2.484*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.105*eV };
    /*G4double RefractiveIndex_SiO2[5] = { 1.54, 1.54, 1.54, 1.54, 1.54 };
    G4double Absorption_SiO2[5] = {125.*cm, 123.5*cm, 122.*cm, 121.*cm, 120.*cm};

    fSiO2MPT = new G4MaterialPropertiesTable();
    fSiO2MPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_SiO2, 5);
    fSiO2MPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_SiO2, 5);*/
    
    G4double refIndexSiO2[3] = {1.458, 1.458, 1.458};
    G4double absorptionSiO2[3] = {125*cm, 125*cm, 125*cm};

    fSiO2MPT = new G4MaterialPropertiesTable();
    fSiO2MPT->AddProperty("RINDEX", greasePhotonEnergy, refIndexSiO2, 3);
    fSiO2MPT->AddProperty("ABSLENGTH", greasePhotonEnergy, absorptionSiO2, 3);   
    fSiO2->SetMaterialPropertiesTable(fSiO2MPT);

	/////////////////////////////////////////////////////////////////
	// Silicon (Si)
	/////////////////////////////////////////////////////////////////

    fSilicon = nist.searchForMaterial("G4_Si");

    // optical properties,    
    fSiliconMPT = new G4MaterialPropertiesTable();
    
    /*G4double RefractiveReal_Si[nEntries_Sil] = { 4.293, 4.453, 4.676, 5.008, 5.587 };
    G4double RefractiveImg_Si[nEntries_Sil] = { 0.045, 0.060, 0.091, 0.150, 0.303 };
    fSiliconMPT->AddProperty("REALRINDEX", PhotonEnergy, RefractiveReal_Si, nEntries_Sil);
    fSiliconMPT->AddProperty("IMAGINARYRINDEX", PhotonEnergy, RefractiveImg_Si, nEntries_Sil);*/
    
    // Ideal detector
    G4double EfficiencyIndex_Si[5] = { 1., 1., 1., 1., 1. };
    G4double Reflective_Si[5] = { 0., 0., 0., 0., 0.};
    
    // Non-ideal
    //G4double EfficiencyIndex_Si[nEntries_Sil] = { 0.37, 0.42, 0.39, 0.36, 0.32 };        
    //G4double Reflective_Si[nEntries_Sil] = { 0.49, 0.45, 0.42, 0.40, 0.39};

    fSiliconMPT->AddProperty("EFFICIENCY", PhotonEnergy, EfficiencyIndex_Si, 5);
    fSiliconMPT->AddProperty("REFLECTIVITY", PhotonEnergy, Reflective_Si, 5);

    fSilicon->SetMaterialPropertiesTable(fSiliconMPT);

	/////////////////////////////////////////////////////////////////
	// ACRYLIC (C5O2H8)n -- CRT
	/////////////////////////////////////////////////////////////////

    fAcrylic = nist.searchForMaterial("G4_PLEXIGLASS");

	// Photon energy (eV)
	G4double ENERGY_ACRYLIC[11] = {6.1992*eV, 4.95936*eV, 4.1328*eV, 3.5424*eV, 3.0996*eV, 2.7552*eV, 2.47968*eV, 2.25426*eV, 2.0664*eV, 1.90745*eV, 1.7712*eV};
	
	// Refractive index
	G4double RINDEX_ACRYLIC[11] = {1.57237, 1.54724, 1.5286, 1.51533, 1.50629, 1.50033, 1.49633, 1.49313, 1.4896, 1.48461, 1.47702};
	
	G4MaterialPropertiesTable *MPT_Acrylic = new G4MaterialPropertiesTable();
	MPT_Acrylic->AddProperty("RINDEX", ENERGY_ACRYLIC, RINDEX_ACRYLIC, 11);

	// Photon energy (eV)
	G4double energy2[25] = {12.3984*eV, 5.0292*eV, 4.75755*eV, 4.69897*eV, 4.66072*eV, 4.61377*eV, 4.56776*eV, 4.5316*eV, 4.48721*eV, 
	                        4.43507*eV, 4.13105*eV, 3.87258*eV, 3.64454*eV, 3.43671*eV, 3.25129*eV, 3.10158*eV, 2.94219*eV, 2.81212*eV, 
	                        2.69307*eV, 2.58078*eV, 2.47479*eV, 2.38212*eV, 2.29384*eV, 2.21614*eV, 1.7712*eV};
	                       
	// Acrylic absorption length (mm)
	G4double abslength[25] = {0, 0, 1.02102, 1.28345, 1.8604, 2.44271, 3.20801, 4.15751, 5.55214, 6.99127, 13.0334, 19.3961, 27.6547, 
	                          32.87, 34.1447, 35.5173, 34.1447, 32.87, 32.87, 34.1447, 34.1447, 35.5173, 35.5173, 34.1447, 33.7719};
	
	MPT_Acrylic->AddProperty("ABSLENGTH", energy2, abslength, 25);

	fAcrylic->SetMaterialPropertiesTable(MPT_Acrylic);

	/////////////////////////////////////////////////////////////////
	// Natural Aluminum
	/////////////////////////////////////////////////////////////////

	fAluminum = nist.searchForMaterial("G4_Al");

	double AlEnergies[3] = {2.0*eV, 3.0*eV, 4.0*eV};	                       
	double AlRefIndex[3] = {0.86, 0.50, 0.28};
	double AlAbsLength[3] = {1*mm, 1*mm, 1*mm};

	// H. EHRENREICH et al. Phys. Rev. 132, 5 (1963)
	double aluminumEnergy[8] = {0.000*eV, 0.697*eV, 1.172*eV, 1.350*eV, 1.504*eV, 2.305*eV, 3.174*eV, 4.000*eV};
	double aluminumReflectivity[8] = {1.000, 0.977, 0.950, 0.911, 0.869, 0.914, 0.921, 0.922};

	fAluminumMPT = new G4MaterialPropertiesTable();
	fAluminumMPT->AddProperty("RINDEX", AlEnergies, AlRefIndex, 3);
	fAluminumMPT->AddProperty("ABSLENGTH", AlEnergies, AlAbsLength, 3);
	fAluminumMPT->AddProperty("REFLECTIVITY", aluminumEnergy, aluminumReflectivity, 8);
	fAluminum->SetMaterialPropertiesTable(fAluminumMPT);

	/////////////////////////////////////////////////////////////////
	// Optical Surfaces
	/////////////////////////////////////////////////////////////////

	// Teflon	
    fTeflonOpSurf = new G4OpticalSurface("TeflonSurface", glisur, ground, dielectric_metal, 0.1); // polish level

    fTeflonOpSurf->SetFinish(ground);
    //fTeflonOpSurf->SetType(dielectric_metal);
    fTeflonOpSurf->SetType(dielectric_dielectric);
    fTeflonOpSurf->SetMaterialPropertiesTable(fTeflonMPT);

    // Silicon
    fSiliconOpSurf = new G4OpticalSurface("SiPMSurface");
    fSiliconOpSurf->SetType(dielectric_metal);
    fSiliconOpSurf->SetFinish(polished);
    fSiliconOpSurf->SetModel(glisur);
    fSiliconOpSurf->SetMaterialPropertiesTable(fSiliconMPT);

	// Aluminized mylar
    G4Material *Mylar = nist.searchForMaterial("G4_MYLAR");

    fMylar = new G4Material("AluminizedMylar", 1.39*g/cm3, 2);
    fMylar->AddMaterial(Mylar, 0.8);
    fMylar->AddMaterial(fAluminum, 0.2);

    //G4double RefractiveReal_Mylar[5] = {0.81257,0.72122,0.63324,0.55571,0.48787};
    //G4double RefractiveImg_Mylar[5] = {6.0481,5.7556,5.4544,5.1464,4.8355};

	//fMylarMPT = new G4MaterialPropertiesTable();
	//fMylarMPT->AddProperty("REALRINDEX", PhotonEnergy, RefractiveReal_Mylar, 5);
	//fMylarMPT->AddProperty("IMAGINARYRINDEX", PhotonEnergy, RefractiveImg_Mylar, 5);

    fMylarOpSurf = new G4OpticalSurface("MylarSurface");
	fMylarOpSurf->SetType(dielectric_metal);
	fMylarOpSurf->SetFinish(polished); // dielectric_metal only allows polished or ground. Polished dielectric_metal uses only reflectivity or absorption.
	fMylarOpSurf->SetModel(glisur);
	fMylarOpSurf->SetMaterialPropertiesTable(fAluminumMPT);

	// 100% reflectivity
	double perfectEfficiency[3] = {0.0, 0.0, 0.0};
	double perfectReflectivity[3] = {1.0, 1.0, 1.0};
	double perfectSpecularSpike[3] = {1.0, 1.0, 1.0};
	
	// ESR (98% reflectivity)
	double esrReflectivity[3] = {0.98, 0.98, 0.98};
    
    fEsrMPT = new G4MaterialPropertiesTable();
	fEsrMPT->AddProperty("EFFICIENCY", AlEnergies, perfectEfficiency, 3);
	fEsrMPT->AddProperty("REFLECTIVITY", AlEnergies, esrReflectivity, 3);    
    
    // ESR film (built in look-up-table)
    fEsrOpSurf = new G4OpticalSurface("EsrSurface");
    fEsrOpSurf->SetType(dielectric_LUT);
    fEsrOpSurf->SetModel(LUT);    
    //fEsrOpSurf->SetFinish(polishedvm2000air);
    fEsrOpSurf->SetFinish(polishedvm2000glue);
    fEsrOpSurf->SetMaterialPropertiesTable(fEsrMPT);

	fPerfectMPT = new G4MaterialPropertiesTable();
	fPerfectMPT->AddProperty("EFFICIENCY", AlEnergies, perfectEfficiency, 3);
	fPerfectMPT->AddProperty("REFLECTIVITY", AlEnergies, perfectReflectivity, 3);
	fPerfectMPT->AddProperty("SPECULARSPIKECONSTANT", AlEnergies, perfectSpecularSpike, 3);

	fPerfectOpSurf = new G4OpticalSurface("PerfectReflector");
	fPerfectOpSurf->SetType(dielectric_metal);
	fPerfectOpSurf->SetFinish(polished);
	fPerfectOpSurf->SetModel(glisur);
	fPerfectOpSurf->SetMaterialPropertiesTable(fPerfectMPT);

	fGreaseOpSurf = new G4OpticalSurface("GreaseSurface");
	fGreaseOpSurf->SetType(dielectric_dielectric);
	fGreaseOpSurf->SetFinish(ground);
	fGreaseOpSurf->SetModel(unified); // Defaults to Lambertian reflection (i.e. rough surface) --CRT
	fGreaseOpSurf->SetMaterialPropertiesTable(fGreaseMPT);	

	isInitialized = true;
}

void nDetMaterials::defineScintillators(){
	if(scintsAreDefined){
		delete fEJ200;
		delete fEJ276;
		delete fEJ200MPT;
		delete fEJ276MPT;
	}

	/////////////////////////////////////////////////////////////////
	// EJ200 N(H)=52.4%, N(C)=47.6%
	/////////////////////////////////////////////////////////////////

    fEJ200 = new G4Material("EJ200", 1.023*g/cm3, 2);
    fEJ200->AddElement(fH, 0.08457);
    fEJ200->AddElement(fC, 0.91543);

	G4double photonEnergy_Ej200[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV, 
		                               2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV, 
		                               2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV, 
		                               2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV, 
		                               3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};

	G4double ScintilFast_EJ200[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060, 
		                              0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594, 
		                              0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798, 
		                              0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170, 
		                              0.105, 0.028, 0.004, 0.000};
                                  
	G4double photonEnergy_Ej200_2[2] = {2.004*eV, 3.191*eV};
	G4double RefIndex_EJ200[2] = {1.580, 1.580};
	G4double Absorption_EJ200[2] = {400*cm, 400*cm};

    fEJ200MPT = new G4MaterialPropertiesTable();
    fEJ200MPT->AddProperty("RINDEX", photonEnergy_Ej200_2, RefIndex_EJ200, 2);
    fEJ200MPT->AddProperty("ABSLENGTH", photonEnergy_Ej200_2, Absorption_EJ200, 2);
    fEJ200MPT->AddProperty("FASTCOMPONENT", photonEnergy_Ej200, ScintilFast_EJ200, 44);

    //fEJ200MPT->AddConstProperty("SCINTILLATIONYIELD", 0.64*17400/MeV); // 64% of Anthracene
    fEJ200MPT->AddConstProperty("SCINTILLATIONYIELD", 10000/MeV); // Scintillation efficiency as per Eljen specs
    fEJ200MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

    //fEJ200MPT->AddConstProperty("RISETIMECONSTANT", 0.9*ns); Geant4 10.1 TODO
    fEJ200MPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.9*ns);
    fEJ200MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
    fEJ200MPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	std::cout << "nDetConstruction: Photon yield is set to " << lightYieldScale << "x scale\n";
    G4double pEF = lightYieldScale; 
    G4double pSF = pEF * 1.35;

    //light yield - data taken form V.V. Verbinski et al, Nucl. Instrum. & Meth. 65 (1968) 8-25
	G4double particleEnergy[36] = {1E-3, 0.10*MeV, 0.13*MeV, 0.17*MeV, 0.20*MeV, 0.24*MeV, 0.30*MeV, 0.34*MeV, 0.40*MeV, 
		                           0.48*MeV, 0.60*MeV, 0.72*MeV, 0.84*MeV, 1.00*MeV, 1.30*MeV, 1.70*MeV, 2.00*MeV, 2.40*MeV, 
		                           3.00*MeV, 3.40*MeV, 4.00*MeV, 4.80*MeV, 6.00*MeV, 7.20*MeV, 8.40*MeV, 10.00*MeV, 11.00*MeV, 
		                           12.00*MeV, 13.00*MeV, 14.00*MeV, 15.00*MeV, 16.00*MeV, 17.00*MeV, 18.00*MeV, 19.00*MeV, 20.00*MeV};
		                           
	G4double electronYield[36] = {0.0E+00*pEF, 1.0E+03*pEF, 1.3E+03*pEF, 1.7E+03*pEF, 2.0E+03*pEF, 2.4E+03*pEF, 3.0E+03*pEF, 3.4E+03*pEF, 4.0E+03*pEF, 
		                          4.8E+03*pEF, 6.0E+03*pEF, 7.2E+03*pEF, 8.4E+03*pEF, 1.0E+04*pEF, 1.3E+04*pEF, 1.7E+04*pEF, 2.0E+04*pEF, 2.4E+04*pEF, 
		                          3.0E+04*pEF, 3.4E+04*pEF, 4.0E+04*pEF, 4.8E+04*pEF, 6.0E+04*pEF, 7.2E+04*pEF, 8.4E+04*pEF, 1.0E+05*pEF, 1.1E+05*pEF, 
		                          1.2E+05*pEF, 1.3E+05*pEF, 1.4E+05*pEF, 1.5E+05*pEF, 1.6E+05*pEF, 1.7E+05*pEF, 1.8E+05*pEF, 1.9E+05*pEF, 2.0E+05*pEF};

    fEJ200MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield, 36)->SetSpline(true);
	G4double protonYield[36] = {0.6*pSF, 67.1*pSF, 88.6*pSF, 120.7*pSF, 146.5*pSF, 183.8*pSF, 246.0*pSF, 290.0*pSF, 365.0*pSF, 
		                        483.0*pSF, 678.0*pSF, 910.0*pSF, 1175.0*pSF, 1562.0*pSF, 2385.0*pSF, 3660.0*pSF, 4725.0*pSF, 6250.0*pSF, 
		                        8660.0*pSF, 10420.0*pSF, 13270.0*pSF, 17180.0*pSF, 23100.0*pSF, 29500.0*pSF, 36200.0*pSF, 45500.0*pSF, 51826.7*pSF, 
		                        58313.7*pSF, 65047.2*pSF, 72027.4*pSF, 79254.2*pSF, 86727.6*pSF, 94447.6*pSF, 102414.2*pSF, 110627.4*pSF, 119087.2*pSF};

    fEJ200MPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield, 36)->SetSpline(true);

	G4double ionYield[36] = {0.2*pEF, 10.4*pEF, 12.7*pEF, 15.7*pEF, 17.9*pEF, 20.8*pEF, 25.1*pEF, 27.9*pEF, 31.9*pEF, 
		                     36.8*pEF, 43.6*pEF, 50.2*pEF, 56.9*pEF, 65.7*pEF, 81.3*pEF, 101.6*pEF, 116.5*pEF, 136.3*pEF, 
		                     166.2*pEF, 187.1*pEF, 218.6*pEF, 260.5*pEF, 323.5*pEF, 387.5*pEF, 451.5*pEF, 539.9*pEF, 595.5*pEF, 
		                     651.8*pEF, 708.7*pEF, 766.2*pEF, 824.2*pEF, 882.9*pEF, 942.2*pEF, 1002.1*pEF, 1062.6*pEF, 1123.7*pEF}; 
                         
    fEJ200MPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield, 36)->SetSpline(true);

    fEJ200->SetMaterialPropertiesTable(fEJ200MPT);

	/////////////////////////////////////////////////////////////////
	// EJ276 N(H)=48.1%, N(C)=51.9%
	/////////////////////////////////////////////////////////////////

    fEJ276 = new G4Material("EJ276", 1.096*g/cm3, 2);
    fEJ276->AddElement(fH, 0.07216);
    fEJ276->AddElement(fC, 0.92784);

	G4double photonEnergy_Ej276[36] = {3.131*eV, 3.087*eV, 3.060*eV, 3.044*eV, 3.029*eV, 3.017*eV, 3.010*eV, 3.001*eV, 2.993*eV, 2.984*eV, 
		                               2.976*eV, 2.967*eV, 2.959*eV, 2.950*eV, 2.941*eV, 2.910*eV, 2.857*eV, 2.838*eV, 2.821*eV, 2.802*eV, 
		                               2.784*eV, 2.764*eV, 2.739*eV, 2.705*eV, 2.671*eV, 2.646*eV, 2.625*eV, 2.599*eV, 2.567*eV, 2.533*eV, 
		                               2.500*eV, 2.468*eV, 2.437*eV, 2.406*eV, 2.377*eV, 2.350*eV};

	G4double ScintilFast_EJ276[36] = {0.000, 0.010, 0.088, 0.157, 0.225, 0.293, 0.354, 0.415, 0.492, 0.570, 
	                                  0.649, 0.730, 0.807, 0.882, 0.934, 1.000, 0.890, 0.826, 0.761, 0.692, 
	                                  0.629, 0.569, 0.509, 0.445, 0.388, 0.326, 0.263, 0.200, 0.144, 0.100, 
	                                  0.068, 0.038, 0.024, 0.012, 0.002, 0.000};

    fEJ276MPT = new G4MaterialPropertiesTable();
    fEJ276MPT->AddProperty("RINDEX", photonEnergy_Ej200_2, RefIndex_EJ200, 2);
    fEJ276MPT->AddProperty("ABSLENGTH", photonEnergy_Ej200_2, Absorption_EJ200, 2);
    fEJ276MPT->AddProperty("FASTCOMPONENT", photonEnergy_Ej276, ScintilFast_EJ276, 36);

    fEJ276MPT->AddConstProperty("SCINTILLATIONYIELD", 8600/MeV); // Scintillation efficiency as per Eljen specs
    fEJ276MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

    fEJ276MPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.9*ns);
    fEJ276MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
    fEJ276MPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	G4double electronYield_EJ276[36];
	G4double protonYield_EJ276[36];
	G4double ionYield_EJ276[36];

	// Produce the scaled light-yield for EJ276 (scaled down from EJ200 by 14%).
	for(size_t i = 0; i < 36; i++){
		electronYield_EJ276[i] = 0.86 * electronYield[i];
		protonYield_EJ276[i] = 0.86 * protonYield[i];
		ionYield_EJ276[i] = 0.86 * ionYield[i];
	}

    fEJ276MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield_EJ276, 36)->SetSpline(true);
    fEJ276MPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield_EJ276, 36)->SetSpline(true);
    fEJ276MPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield_EJ276, 36)->SetSpline(true);

	fEJ276->SetMaterialPropertiesTable(fEJ276MPT);
	
	// Update the material dictionary
	materialList["ej200"] = fEJ200;
	materialList["ej276"] = fEJ276;
	
	scintsAreDefined = true;
}
