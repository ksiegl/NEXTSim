#ifndef NDET_MATERIALS_HH
#define NDET_MATERIALS_HH

#include "globals.hh"

#include "nistDatabase.hh"

class G4Element;
class G4Material;
class G4MaterialPropertiesTable;
class G4OpticalSurface;
class G4VisAttributes;

class nDetMaterialsMessenger;

/** @class nDetMaterials
  * @brief Geant materials and elements used for detector construction
  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu)
  * @date July 3, 2019
  *
  * This class is designed as a handler for common Geant materials which are used by
  * other NEXTSim classes. The materials and elements defined in this class are expandable
  * by using the nistDatabase class.
  */

class nDetMaterials{
  public:
    // Materials and elements
    G4Element* fH; ///< Hydrogen
    G4Element* fC; ///< Carbon
    G4Element* fO; ///< Oxygen
    G4Element* fF; ///< Flourine
    G4Element* fSi; ///< Silicon
    G4Element* fAl; ///< Aluminium
    G4Element* fY; ///<Yttrium
    G4Element* fLa; ///<Lanthanum
    G4Element* fBr; ///<Bromine
    G4Element* fNa; ///<Sodium
    G4Element* fI; ///<Iodine
    G4Element* fCr; ///<Chromium
    G4Element* fFe; ///<Iron
    G4Element* fNi; ///<Nickel
    G4Element* fCu; ///<Copper

    G4Material* fAir; ///< Material corresponding to air
    G4Material* fVacuum; ///< Material corresponding to natural vacuum
    G4Material* fTeflon; ///< Material corresponding to teflon
    G4Material* fEJ200; ///< Material corresponding to EJ-200 scintillator
    G4Material* fEJ276; ///< Material corresponding to EJ-276 scintillator
    G4Material* fGrease; ///< Material corresponding to optical grease
    G4Material* fSiO2; ///< Material corresponding to quartz
    G4Material* fSilicon; ///< Material corresponding to silicon
    G4Material* fMylar; ///< Material corresponding to aluminized mylar
    G4Material* fAcrylic; ///< Material corresponding to acrylic
    G4Material* fAluminum; ///< Material corresponding to aluminum
    G4Material* fLaBr3; ///< Lanthanum Bromide detector material
    G4Material* fNaI; ///< Sodium Iodide detector material
    G4Material* fYAP; ///< Yttrium Aluminum Perovskite detecor material
    G4Material* fYSO; ///< Yttrium Orthosilicate detecor material
    G4Material* fSteel; ///< Stainless Steel material
    G4Material* fMuMetal; ///< Mu Metal material
    G4Material* fHDPE; ///< High Density Polyethylene material
    G4Material* fPLA; ///<Material correspondint to PLA


    // Material table properties
    G4MaterialPropertiesTable* fAirMPT; ///< Material properties table for air
    G4MaterialPropertiesTable* fTeflonMPT; ///< Material properties table for teflon
    G4MaterialPropertiesTable* fEJ200MPT; ///< Material properties table for EJ-200 scintillator
    G4MaterialPropertiesTable* fEJ276MPT; ///< Material properties table for EJ-276 scintillator
    G4MaterialPropertiesTable* fGreaseMPT; ///< Material properties table for optical grease
    G4MaterialPropertiesTable* fSiO2MPT; ///< Material properties table for quartz
    G4MaterialPropertiesTable* fSiliconMPT; ///< Material properties table for silicon
    G4MaterialPropertiesTable* fMylarMPT; ///< Material properties table for aluminized mylar
    G4MaterialPropertiesTable* fPerfectMPT; ///< Material properties table for a perfect reflector
    G4MaterialPropertiesTable* fAluminumMPT; ///< Material properties table for aluminum
    G4MaterialPropertiesTable* fEsrMPT; ///< Material properties table for 3M Enhanced Specular Reflector

    // Optical Surfaces
    G4OpticalSurface* fTeflonOpSurf; ///< Optical surface for teflon
    G4OpticalSurface* fSiliconOpSurf; ///< Optical surface for silicon
    G4OpticalSurface* fMylarOpSurf; ///< Optical surface for aluminized mylar
    G4OpticalSurface* fEsrOpSurf; ///< Optical surface for 3M Enhanced Specular Reflector
	G4OpticalSurface* fPerfectOpSurf; ///< Optical surface for a perfect reflector
	G4OpticalSurface* fGreaseOpSurf; ///< Optical surface for optical grease

	// Visual attributes
	G4VisAttributes *visAssembly; ///< Visual attributes for the mother assembly
	G4VisAttributes *visSensitive; ///< Visual attributes for the photo-sensitive surface
	G4VisAttributes *visWindow; ///< Visual attributes for the optical window
	G4VisAttributes *visGrease; ///< Visual attributes for the optical grease
	G4VisAttributes *visWrapping; ///< Visual attributes for the inner/outer detector wrappint
	G4VisAttributes *visScint; ///< Visual attributes for the scintillator material
	G4VisAttributes *visShadow; ///< Visual attributes for the shadow object
	
	/** Default constructor
	  */
	nDetMaterials();
	
	/** Destructor
	  */
	~nDetMaterials();

	/** Define all materials, scintillators, and optical surfaces
	  */
	void initialize();

	/** Return true if all materials and surfaces are defined and return false otherwise
	  */
	bool materialsAreDefined() const { return isInitialized; }

	/** Get a pointer to the scintillator material corresponding to @a name
	  * 
	  * Valid material names are shown in the table below
	  * | Name  | Description |
	  * |-------|-------------|
	  * | ej200 | <a href="https://eljentechnology.com/products/plastic-scintillators/ej-200-ej-204-ej-208-ej-212">Eljen EJ-200</a> plastic scintillator
	  * | ej276 | <a href="https://eljentechnology.com/products/plastic-scintillators/ej-276">Eljen EJ-276</a> plastic scintillator
	  *
	  * @return A pointer to the scintillator material if the name is valid or a pointer to EJ200 if it is invalid
	  */
	G4Material* getUserDetectorMaterial(const G4String &name);

	/** Get a pointer to the surface material corresponding to @a name
	  *
	  * Valid material names are shown in the table below
	  * | Name    | Description |
	  * |---------|-------------|
	  * | mylar   | Aluminized mylar (80% mylar, 20% aluminum by mass)
	  * | teflon  | Polytetrafluoroethylene (C2F4)n (PTFE)
	  * | esr     | <a href="https://www.3m.com/3M/en_US/company-us/all-3m-products/~/3M-Enhanced-Specular-Reflector-3M-ESR-/?N=5002385+3293061534&rt=rud">3M Enhanced Specular Reflector</a>
	  * | silicon | Natural silicon (from the NIST database)
	  * | perfect | Perfect optical reflector (reflectivity = 1, specular spike = 1)
	  * 
	  * @return A pointer to the surface material if the name is valid or a pointer to Mylar if it is invalid
	  */
	G4Material* getUserSurfaceMaterial(const G4String &name);

	/** Get a pointer to the optical surface corresponding to @a name
	  * @note See getUserSurfaceMaterial() for a description of valid optical surface names
	  * @return A pointer to the optical surface if the name is valid or a pointer to Mylar optical surface if it is invalid	
	  */
	G4OpticalSurface* getUserOpticalSurface(const G4String &name);

	/** Get a pointer to the visual attributes corresponding to @a name
	  */
	G4VisAttributes* getUserVisAttributes(const G4String &name);

	/** Get the current light yield multiplier used for scintillator photon production
	  */
	G4double getLightYield() const { return lightYieldScale; }

	/** Get a pointer to the database for pre-defined NIST element and materials
	  */
	nistDatabase* getNistDatabase(){ return &nist; }

	/** Get a pointer to the messenger for this class
	  */
	nDetMaterialsMessenger* getMessenger(){ return messenger; }

	/** Set the light yield multiplier for scintillator photon production
	  * @param yield The fraction of the light yield to use for optical photon production in scintillators (default is 1)
	  */	
	void setLightYield(const G4double &yield);

	/** Search for an element in the NIST database 
	  * @return True if the element is found in the NIST database and return false otherwise
	  */
	bool searchForElement(const G4String &name);

	/** Search for an element name in the pre-defined dictionary
	  * @note If the specified name corresponds to a pre-defined element, a pointer to that element is returned.
	  *       Otherwise the name is searched for within the Geant NIST element database.
	  * @return A pointer to the Geant element if it is found in the pre-defined dictionary OR the NIST database and return NULL otherwise
	  */
	G4Element* getElement(const G4String &name);

	/** Search for a material in the NIST database 
	  * @return True if the material is found in the NIST database and return false otherwise
	  */
	bool searchForMaterial(const G4String &name);
	
	/** Search for a material name in the pre-defined dictionary
	  * @note If the specified name corresponds to a pre-defined material, a pointer to that material is returned.
	  *       Otherwise the name is searched for within the Geant NIST material database.
	  * @return A pointer to the Geant material if it is found in the pre-defined dictionary OR the NIST database and return NULL otherwise
	  */
	G4Material* getMaterial(const G4String &name);

	/** Search for an optical surface name in the pre-defined dictionary
	  * @return A pointer to the Geant optical surface if it is found in the pre-defined dictionary and return NULL otherwise
	  */
	G4OpticalSurface* getOpticalSurface(const G4String &name);

	/** Search for an object visual attributes name in the pre-defined dictionary
	  * @return A pointer to the Geant visual attributes if it is found in the pre-defined dictionary and return NULL otherwise
	  */	
	G4VisAttributes* getVisualAttributes(const G4String &name);

	/** List all materials in the materials dictionary
	  */
	void listMaterials() const ;

	/** List all visual attributes in the dictionary
	  */
	void listVisAttributes() const ;

	/** List all optical surfaces in the dictionary
	  */
	void listOptSurfaces() const ;

	/** Read an input material file and build a new material to add to the list
	  * @param filename The path to the material file to read
	  * @return True if the file is read successfully and return false otherwise
	  */
	bool buildNewMaterial(const G4String &filename);

	/**
	  */
	void printMaterial(const G4String &name);

  private:
	nDetMaterialsMessenger* messenger; ///< Pointer to the messenger for this class

	bool isInitialized; ///< Flag indicating that all materials have been defined
	bool scintsAreDefined; ///< Flag indicating that the scintillator materials have been defined
	
    G4double lightYieldScale; ///< Scaling parameter used for scintillation light yield (default=1)	

	nistDatabase nist; ///< Database for pre-defined NIST element and material lookups

	std::map<G4String, G4Element*> elementList; ///< List of elements which are available for use
	
	std::map<G4String, G4Material*> materialList; ///< List of materials which are available for use
	
	std::map<G4String, G4OpticalSurface*> opticalSurfaceList; ///< List of optical surfaces which are available for use
	
	std::map<G4String, G4VisAttributes*> visAttributesList; ///< List of visual attributes which are available for use

	/** Define materials and optical surfaces
	  */
	void defineMaterials();

	/** Define scintillator materials EJ200 and EJ276
	  */
	void defineScintillators();
};

#endif
