#include <iostream>
#include <fstream>

#include "G4Step.hh"

#include "nDetDetector.hh"
#include "centerOfMass.hh"
#include "vertilon.hh"

const double coeff = 1.23984193E-3; // hc = Mev * nm

///////////////////////////////////////////////////////////////////////////////
// class centerOfMass
///////////////////////////////////////////////////////////////////////////////

centerOfMass::~centerOfMass(){
}

centerOfMass centerOfMass::clone() const {
	centerOfMass retval;
	retval.response = response.clone();
	for(size_t i = 0; i < 4; i++)
		retval.anodeResponse[i] = anodeResponse[i].clone();
	retval.gainMatrix = gainMatrix;
	retval.countMatrix = countMatrix;
	return retval;
}

G4ThreeVector centerOfMass::getCenter() const {
	return (totalMass > 0 ? (1/totalMass)*center : center);
}

double centerOfMass::getCenterX() const {
	return (totalMass > 0 ? (1/totalMass)*center.getX() : 0);
}

double centerOfMass::getCenterY() const {
	return (totalMass > 0 ? (1/totalMass)*center.getY() : 0);
}

double centerOfMass::getCenterZ() const {
	return (totalMass > 0 ? (1/totalMass)*center.getZ() : 0);
}

bool centerOfMass::getCenterSegment(G4ThreeVector &pos, short &col, short &row) const {
	double xpos=0;
	double ypos=0;
	if(pixelWidth>0 && pixelHeight>0){
	xpos = (pos.getX()+activeWidth/2)/pixelWidth;
	ypos = (pos.getY()+activeHeight/2)/pixelHeight;
	}
	else{
	//std::cout << "Error zero pixel size" << std::endl;
	}
	col = (short)floor(xpos);
	row = (short)floor(ypos);
	
	//std::cout << " activeWidth=" << activeWidth << ", activeHeight=" << activeHeight << std::endl;
	//std::cout << " pixelWidth=" << pixelWidth << ", pixelHeight=" << pixelHeight << std::endl;
	//std::cout << " x=" << xpos << ", y=" << ypos << ", col=" << col << ", row=" << row << std::endl;
	
	return ((col >= 0 && col < Ncol) && (row >= 0 && row < Nrow));
}

bool centerOfMass::getCenterSegment(short &col, short &row) const {
	G4ThreeVector pos = this->getCenter();
	return this->getCenterSegment(pos, col, row);
}	

void centerOfMass::getAnodeCurrents(double *array) const {
	for(size_t i = 0; i < 4; i++){
		array[i] = anodeCurrent[i];
	}
}

double centerOfMass::getReconstructedX() const {
	return ((anodeCurrent[0]+anodeCurrent[1])-(anodeCurrent[2]+anodeCurrent[3]))/(anodeCurrent[0]+anodeCurrent[1]+anodeCurrent[2]+anodeCurrent[3]);
}

double centerOfMass::getReconstructedY() const {
	return ((anodeCurrent[1]+anodeCurrent[2])-(anodeCurrent[3]+anodeCurrent[0]))/(anodeCurrent[0]+anodeCurrent[1]+anodeCurrent[2]+anodeCurrent[3]);
}

short centerOfMass::setNumColumns(const short &col_){ 
	Ncol = col_;
	pixelWidth = activeWidth / (Ncol > 0 ? Ncol : 1);
	return Ncol; 
}

short centerOfMass::setNumRows(const short &row_){
	Nrow = row_;
	pixelHeight = activeHeight / (Nrow > 0 ? Nrow : 1);
	return Nrow; 
}

double centerOfMass::setActiveAreaWidth(const double &width_){
	activeWidth = width_;
	pixelWidth = activeWidth / (Ncol > 0 ? Ncol : 1);
	return activeWidth;
}

double centerOfMass::setActiveAreaHeight(const double &height_){
	activeHeight = height_;
	pixelHeight = activeHeight / (Nrow > 0 ? Nrow : 1);
	return activeHeight;
}

void centerOfMass::setSegmentedPmt(const nDetDetectorParams *params){
	if(params->GetNumPmtColumns() <= 0 || params->GetNumPmtRows() <= 0)
		return;

	Ncol = params->GetNumPmtColumns();
	Nrow = params->GetNumPmtRows();
	activeWidth = params->GetPmtWidth();
	activeHeight = params->GetPmtHeight();
	pixelWidth = activeWidth / Ncol;
	pixelHeight = activeHeight / Nrow;
	
	// Setup the anode gain matrix.
	gainMatrix.clear();
	countMatrix.clear();
	for(short i = 0; i < Ncol; i++){
		gainMatrix.push_back(std::vector<double>(Nrow, 100));
		countMatrix.push_back(std::vector<int>(Nrow, 0));
	}
}

bool centerOfMass::loadSpectralResponse(const char *fname){
	return response.loadSpectralResponse(fname);
}

bool centerOfMass::loadGainMatrix(const char *fname){
	if(gainMatrix.empty() || Ncol*Nrow == 0) return false;
	std::ifstream gainFile(fname);
	if(!gainFile.good()) return false;
	
	double readval;
	for(short col = 0; col < Ncol; col++){
		for(short row = 0; row < Nrow; row++){
			gainFile >> readval;
			if(gainFile.eof()){
				gainFile.close();
				return false;
			}
			gainMatrix[col][row] = readval;
		}
	}
	
	gainFile.close();
	
	return true;
}

void centerOfMass::copySpectralResponse(centerOfMass *other){
	response.copySpectralResponse(other->getPmtResponse()->getSpectralResponse());
}

void centerOfMass::copySpectralResponse(const centerOfMass *other){
	response.copySpectralResponse(other->getConstPmtResponse()->getConstSpectralResponse());
}

void centerOfMass::copyGainMatrix(centerOfMass *other){
	other->getGainMatrix(gainMatrix);
}

void centerOfMass::copyGainMatrix(const centerOfMass *other){
	other->getGainMatrix(gainMatrix);
}

void centerOfMass::clear(){
	Npts = 0;
	NnotDetected = 0;
	tSum = 0;
	lambdaSum = 0;
	totalMass = 0;
	center = G4ThreeVector();
	t0 = std::numeric_limits<double>::max();	
	response.clear();
	for(size_t i = 0; i < 4; i++){
		anodeCurrent[i] = 0;
		anodeResponse[i].clear();
	}
	for(short i = 0; i < Ncol; i++){
		for(short j = 0; j < Nrow; j++){
			countMatrix[i][j] = 0;
		}
	}
}

bool centerOfMass::addPoint(const double &energy, const double &time, const G4ThreeVector &position, const double &mass/*=1*/){
	double wavelength = coeff/energy; // in nm
	if(Ncol < 0 && Nrow < 0){ // Default behavior
		center += mass*position;	
		
		// Add the PMT response to the "digitized" trace
		response.addPhoton(time, wavelength);
		
		// Add the "mass" to the others
		totalMass += mass;		
	}
	else{ // Segmented PMT behavior
		short xpos, ypos;
		G4ThreeVector pos = position;
		//std::cout << "pre: " << pos.getX() << "\t" << pos.getY() << std::endl;
		if(this->getCenterSegment(pos, xpos, ypos)){
			// Get the gain of this anode.
			double gain = getGain(xpos, ypos);
			increment(xpos, ypos);

			// Add the anger logic currents to the anode outputs.
			double totalCurrent = 0;
			double *current = getCurrent(xpos, ypos);
			if(current){
				for(size_t i = 0; i < 4; i++){
					anodeCurrent[i] += gain*mass*current[i];
					totalCurrent += current[i];
				}
				for(size_t i = 0; i < 4; i++){
					anodeResponse[i].addPhoton(time, 0, gain*mass*(current[i]/totalCurrent));
				}
			}
			
			// Compute resistor network leakage current. This is unnecessary for symmetric leakage... CRT
			/*const double leakage[3][3] = {{1E-3, 1E-2, 1E-3},
			                              {1E-2, 1.00, 1E-2},
			                              {1E-3, 1E-2, 1E-3}};
			for(short anodeX = -1; anodeX <= 1; anodeX++){
				for(short anodeY = -1; anodeY <= 1; anodeY++){
					double *current = getCurrent(xpos+anodeX, ypos+anodeY);
					if(current){ // Add the anger logic currents to the anode outputs.
						for(size_t i = 0; i < 4; i++){
							anodeCurrent[i] += gain*leakage[anodeX+1][anodeY+1]*current[3-i];
						}
					}
				}
			}*/
			
			// Add the PMT response to the "digitized" trace
			response.addPhoton(time, wavelength, gain);

			// Add the "mass" to the others weighted by the individual anode gain
			center += pos;
			totalMass += mass;
		}
	}
	
	tSum += time;
	lambdaSum += wavelength;
	if(time < t0) t0 = time;

	Npts++;
	
	return true;
}

void centerOfMass::printCounts() const {
	for(short i = Nrow-1; i >= 0; i--){
		for(short j = 0; j < Ncol; j++){
			std::cout << countMatrix[j][i] << "\t";
		}
		std::cout << std::endl;
	}		
}

void centerOfMass::print() const {
	if(!empty()){
		std::cout << "M=" << totalMass << ", c=(" << getCenterX() << ", " << getCenterY() << ", " << getCenterZ() << ")\n";
		std::cout << " t0=" << t0 << ", tAvg=" << tSum/Npts << std::endl;
	}
}

void centerOfMass::increment(const int &x, const int &y){
	if((x < 0 || x >= Ncol) || (y < 0 || y >= Nrow) || countMatrix.empty()) return;
	countMatrix[x][y]++;
}

double centerOfMass::getGain(const int &x, const int &y){
	if((x < 0 || x >= Ncol) || (y < 0 || y >= Nrow) || gainMatrix.empty()) return 0;
	return gainMatrix[x][y]/100;
}

double *centerOfMass::getCurrent(const int &x, const int &y){
	if((x < 0 || x >= 8) || (y < 0 || y >= 8)) return NULL;
	return vertilon::currents[x][y];
}
