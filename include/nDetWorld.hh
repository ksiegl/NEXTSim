#ifndef NDET_WORLD_HH
#define NDET_WORLD_HH

#include <vector>

#include "globals.hh"

#include "G4ThreeVector.hh"

class nDetWorldObject;
class nDetWorldPrimitive;
class nDetWorldMessenger;
class nDetMaterials;
class gdmlObject;

class G4Box;
class G4LogicalVolume;
class G4PVPlacement;

/** @class nDetWorld
  * @brief Build experimental hall setup area
  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu)
  * @date August 14, 2019
  */

class nDetWorld{
  public:
	/** Default constructor
	  */
	nDetWorld();
	
	/** Destructor
	  */
	~nDetWorld(){ }

	/** Get a pointer to the geometry of the experimental setup area
	  */
	G4Box *getGeometry(){ return solidV; }
	
	/** Get a pointer to the logical volume of the experimental setup area
	  */
	G4LogicalVolume *getLogicalVolume(){ return logV; }

	/** Get a pointer to the physical volume of the experimental setup area
	  */	
	G4PVPlacement *getPhysicalVolume(){ return physV; }

	/** Get a pointer to the messenger used for this class
	  */
	nDetWorldMessenger *getMessenger(){ return messenger; }

	/** Set the bounding box size of the experimental hall along the X, Y, and Z axes (all in mm)
	  * @note The new world size will take effect the next time nDetConstruction::ConstructDetector() is called
	  */
	void setWorldSize(const G4ThreeVector &size){ hallSize = size; }

	/** Set the name of the material to use to fill the experimental hall
	  * @note Defaults to "air" which corresponds to nDetMaterials::fAir
	  */
	void setWorldMaterial(const G4String &material){ fillMaterial = material; }

	/** Set the size of the pit in the floor of the experimental hall
	  */
	void setWorldFloorPitSize(const G4ThreeVector &size){ floorPitSize = size; }

	/** Setup a floor for the bottom of the experimental hall using parameters from a space-delimited input string
	  * @note String syntax: <surfaceY> <thickness> [material=G4_CONCRETE]
	  * | Parameter | Description |
	  * |-----------|-------------|
	  * | surfaceY  | The vertical distance to the surface of the floor with respect to the world origin (in cm)
	  * | thickness | The thickness of the floor slab (in cm)
	  * | material  | The name of the material to use (default=G4_CONCRETE)
	  * @return True if the correct number of arguments are found in the input string and return false otherwise
	  */
	bool setWorldFloor(const G4String &input);

	/** Build the experiment hall physical and logical volumes
	  */
	void buildExpHall(nDetMaterials *materials);

	/** Set Name for experimental setup to be constructed
	 */
	void SetExp(std::string const expName_){expName = expName_;}

	/** Build the CERN hall structures
	 */
	void BuildCERNStructures();

	/** Build the CERN IDS elements
	 */
	void BuildCERNElements();

	/** Build FDSi Setup
	 */
	void BuildFDSi(nDetMaterials *materials);

	/** Build NSCL e14060 setup
	 */
	void BuildE14060(nDetMaterials *materials);

	/** Add a new Geant primitive object to the world
	  * @param str Space-delimited input string with format specified in nDetWorldPrimitive::decodeArgs()
	  * @return A pointer to the new object if construction was successful and return NULL if construction failed
	  */
	nDetWorldPrimitive *addNewPrimitive(const G4String &str);

	/** List all available Geant shapes which may be added to the experimental setup area
	  */
	void listAllPrimitives();

	/** Print the list of all user-defined 3d objects in the setup area
	  */
	void printDefinedObjects();

	/** Reset the experimental setup area and delete all defined primitive objects
	  */
	void reset();

	/** Load a GDML model from a file using parameters from a space-delimited input string
	  * @note String syntax: <filename> <posX> <posY> <posZ> <rotX> <rotY> <rotZ> <matString>
	  * | Parameter | Description |
	  * |-----------|-------------|
	  * | filename  | Filename of the input GDML file
	  * | posX(Y,Z) | X, Y, and Z position of the center of the model (in mm)
	  * | rotX(Y,Z) | Rotation about the X, Y, and Z axes (in degrees)
	  * | matString | The NIST database name of the material to use for the model
	  * @return A pointer to the gdmlObject containing the model
	  */
	gdmlObject *loadGDML(const G4String &input);
	
  private:
	G4Box *solidV; ///< Geometry of the experimental setup area
	G4LogicalVolume* logV; ///< Logical volume of the experimental setup area
	G4PVPlacement* physV; ///< Physical volume of the experimental setup area

	G4String fillMaterial; ///< The name of the material which fills the experimental hall
	G4String floorMaterial; ///< The name of the material which comprises the floor of the experimental hall

	G4double floorThickness; ///< Thickness of the floor of the experimental hall (in mm)
	G4double floorSurfaceY; ///< Vertical distance to the surface of the floor of the experimental hall (with respect to the origin, in mm)

	G4ThreeVector floorPitSize; ///< The size of the floor pit along the X, Y, and Z axes (all in mm)
	G4ThreeVector hallSize; ///< Size of the experimental hall along the X, Y, and Z axes (all in mm)

	std::vector<nDetWorldObject*> objects; ///< Vector of objects to add to the experimental setup area
	
	std::string expName;

	nDetWorldMessenger *messenger; ///< Pointer to the messenger used for this class
};

#endif
