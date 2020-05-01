
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TGraph.h"

#include "Randomize.hh"

#include "pmtResponse.hh"

const double sqrt2pi = 2.5066282746;

void copyTGraph(TGraph *g, std::vector<double> &xvec, std::vector<double> &yvec){
	double x, y;
	xvec.clear();
	yvec.clear();
	for(int i = 0; i < g->GetN(); i++){
		g->GetPoint(i, x, y);
		xvec.push_back(x);
		yvec.push_back(y);
	}
}

///////////////////////////////////////////////////////////////////////////////
// class spectralResponse
///////////////////////////////////////////////////////////////////////////////

/// Destructor.
spectralResponse::~spectralResponse(){
	clear();
}

spectralResponse spectralResponse::clone() const {
	spectralResponse retval;
	retval.xmin = xmin;
	retval.xmax = xmax;
	retval.size = size;
	retval.wavelength = wavelength;
	retval.percentage = percentage;
	return retval;
}

void spectralResponse::getRange(double &min_, double &max_) const {
	min_ = xmin;
	max_ = xmax;
}

void spectralResponse::copy(const spectralResponse *other){
	clear();
	wavelength = other->getWavelength();
	percentage = other->getPercentage();
	xmin = wavelength.front();
	xmax = wavelength.back();
	size = wavelength.size();
}

/// Load response function from a file.
bool spectralResponse::load(const std::string &fname){
	clear();

	if(fname.find(".root") != std::string::npos){ // Root file. Expects a TGraph named "spec"
		TFile *f = new TFile(fname.c_str(), "READ");
		if(!f->IsOpen()) return false;
		TGraph *g1 = (TGraph*)f->Get("spec");
		if(!g1){
			f->Close();
			return false;
		}
		copyTGraph(g1, wavelength, percentage);
	
		xmin = wavelength.front();
		xmax = wavelength.back();
		size = wavelength.size();
	
		f->Close();
		delete f;	
	}
	else{ // Ascii file
		std::ifstream f(fname.c_str());
		if(!f.good()) return false;
		double xval, yval;
		while(true){
			f >> xval >> yval;
			if(f.eof())
				break;
			wavelength.push_back(xval);
			percentage.push_back(yval);
		}
		
		xmin = wavelength.front();
		xmax = wavelength.back();
		size = wavelength.size();
		
		f.close();
	}
	
	return true;
}

/// Return the PMT quantum efficiency for a given wavelength.
double spectralResponse::eval(const double &lambda_) const {
	if(size == 0) return 0;
	if(lambda_ >= xmin && lambda_ <= xmax){ // Interpolate the distribution.
		for(size_t i = 0; i < size-1; i++){
			if(lambda_ == wavelength.at(i))
				return percentage.at(i);
			else if(lambda_ > wavelength.at(i) && lambda_ <= wavelength.at(i+1))
				return (percentage.at(i)+(percentage.at(i+1)-percentage.at(i))*(lambda_-wavelength.at(i))/(wavelength.at(i+1)-wavelength.at(i)));
		}
	}
	return 0;
}

void spectralResponse::clear(){
	wavelength.clear();
	percentage.clear();
	xmin = 0;
	xmax = 0;
	size = 0;
}

///////////////////////////////////////////////////////////////////////////////
// class pmtResponse
///////////////////////////////////////////////////////////////////////////////

pmtResponse::pmtResponse() : risetime(4.0), falltime(20.0), timeSpread(0), traceDelay(50), gain(1E4), maximum(-9999), baseline(-9999),
                             baselineFraction(0), baselineJitterFraction(0), polyCfdFraction(0.5), adcClockTick(4), tLatch(0), pulseIntegralLow(5), pulseIntegralHigh(10),
                             maxIndex(0), adcBins(4096), pulseLength(100), isDigitized(false), isTriggered(false), useSpectralResponse(false), pulseIsSaturated(false),
                             printTrace(false), pulseArray(), spec(), minimumArrivalTime(0), functionType(EXPO) {
	this->setPulseLength(pulseLength);
}

pmtResponse::pmtResponse(const double &risetime_, const double &falltime_) : risetime(risetime_), falltime(falltime_), timeSpread(0), traceDelay(50), gain(1E4), maximum(-9999), baseline(-9999),
                                                                             baselineFraction(0), baselineJitterFraction(0), polyCfdFraction(0.5), adcClockTick(4), tLatch(0), pulseIntegralLow(5), pulseIntegralHigh(10),
                                                                             maxIndex(0), adcBins(4096), pulseLength(100), isDigitized(false), isTriggered(false), useSpectralResponse(false), pulseIsSaturated(false),
                                                                             printTrace(false), pulseArray(), spec(), minimumArrivalTime(0), functionType(EXPO) {
	this->setPulseLength(pulseLength);
}

pmtResponse::~pmtResponse(){
}

pmtResponse pmtResponse::clone() const {
	pmtResponse retval;
	retval.risetime = risetime;
	retval.falltime = falltime;
	retval.timeSpread = timeSpread;
	retval.traceDelay = traceDelay;
	retval.gain = gain;
	retval.baselineFraction = baselineFraction;
	retval.baselineJitterFraction = baselineJitterFraction;
	retval.polyCfdFraction = polyCfdFraction;
	retval.adcClockTick = adcClockTick;
	retval.pulseIntegralLow = pulseIntegralLow;
	retval.pulseIntegralHigh = pulseIntegralHigh;
	retval.adcBins = adcBins;
	retval.pulseLength = pulseLength;
	retval.useSpectralResponse = useSpectralResponse;
	retval.printTrace = printTrace;
	retval.pulseArray = pulseArray;
	retval.spec = spec.clone();
	retval.functionType = functionType;
	return retval;
}

/// Return the gain-weighted arrival time of the photon pulse.
double pmtResponse::getWeightedPhotonArrivalTime() const {
	double weightedAverage = 0;
	double totalWeight = 0;
	std::vector<double>::const_iterator iter1, iter2;
	for(auto arrival : arrivalTimes){
		weightedAverage += arrival.time * arrival.gain;
		totalWeight += arrival.gain;
	}
	return weightedAverage/totalWeight;
}

void pmtResponse::setRisetime(const double &risetime_){ 
	risetime = risetime_; 
}

void pmtResponse::setFalltime(const double &falltime_){ 
	falltime = falltime_; 
}

/// Set the length of the pulse in ADC bins.
void pmtResponse::setPulseLength(const size_t &len){
	pulseLength = len;
	pulseArray = std::vector<unsigned short>(pulseLength, 0);
	this->clear();
}

/// Set the length of the pulse in nanoseconds.
void pmtResponse::setPulseLengthInNanoSeconds(const double &length){
	size_t len = (size_t)ceil(length/adcClockTick);
	this->setPulseLength(len);
}

/// Set the dynamic bit range of the ADC.
void pmtResponse::setBitRange(const size_t &len){
	adcBins = (unsigned int)std::pow(2, len);
}

/// Load PMT spectral response from root file.
bool pmtResponse::loadSpectralResponse(const char *fname){
	return (useSpectralResponse = spec.load(fname));
}

void pmtResponse::addPhoton(const double &arrival, const double &wavelength/*=0*/, const double &gain_/*=1*/){
	double efficiency = 1;
	if(useSpectralResponse){ // Compute the quantum efficiency of the PMT for this wavelength.
		efficiency = spec.eval(wavelength)/100;
	}

	// Compute the offset of the response function due to the trace delay and the time spread of the PMT
	double dt = arrival + traceDelay; // Arrival time is the leading edge of the pulse.
	if(timeSpread > 0){ // Smear the time offset based on the photo-electron transit time spread.
		dt += (G4UniformRand()-0.5)*timeSpread;
	}

	// Add the photon to the list of arrival times
	arrivalTimes.push_back(photonArrivalTime(arrival, gain_*efficiency, (dt >= 0 ? dt : 0)));

	// Check if this is the first photon
	if(arrival < minimumArrivalTime)
		minimumArrivalTime = arrival;
}

double pmtResponse::sample(const double &time) const {
	double retval = 0;
	for(auto arrival : arrivalTimes)
		retval += arrival.gain * func(time, arrival.dt);
	return retval;
}

void pmtResponse::getRawPulse(std::vector<double> &rawPulse) const {
	rawPulse = std::vector<double>(pulseLength, 0);
	double time = tLatch + adcClockTick/2;
	double prevAmp = 0;
	for(size_t index = 0; index < pulseLength; index++){
		rawPulse[index] = sample(time); // Sample the total light response spectrum
		if(rawPulse[index] < prevAmp && rawPulse[index] < 1){ // Amplitude is less than one ADC bin
			rawPulse[index] = 0;
			break;
		}
		prevAmp = rawPulse[index];
		time += adcClockTick;
	}
}

void pmtResponse::digitize(const double &baseline_, const double &jitter_){
	if(isDigitized) 
		return;
	pulseIsSaturated = false;
	
	// Digitize the light pulse
	unsigned int value, bin;
	double time = tLatch + adcClockTick/2;
	for(size_t i = 0; i < pulseLength; i++){
		bin = (unsigned int)floor(sample(time)); // Sample the total light response spectrum
		if(bin >= adcBins) bin = adcBins-1;
		value = bin;
		value += baseline_*adcBins;
		if(jitter_ != 0) 
			value += (-jitter_ + 2*G4UniformRand()*jitter_)*adcBins;
		if(value <= adcBins-1) // Pulse is not saturated.
			pulseArray[i] = (unsigned short)value;
		else{ // Pulse is saturated.
			pulseIsSaturated = true;
			pulseArray[i] = (unsigned short)(adcBins-1);
		}
		time += adcClockTick;
	}
	isDigitized = true;
}

void pmtResponse::digitize(){
	this->digitize(baselineFraction, baselineJitterFraction);
}

/// Integrate the baseline corrected trace for QDC in the range [start_, stop_] and return the result.
double pmtResponse::integratePulse(const size_t &start_, const size_t &stop_){
	if(pulseLength == 0 || pulseArray.empty()) return -9999;
	size_t stop = (stop_ > pulseLength ? pulseLength : stop_);

	// Check for start index greater than stop index.
	if(start_+1 >= stop) return -9999;

	double qdc = 0.0;
	for(size_t i = start_+1; i < stop; i++){ // Integrate using trapezoidal rule.
		qdc += 0.5*(pulseArray[i-1] + pulseArray[i]) - baseline;
	}
	if(qdc>50) isTriggered = true;
	return qdc;
}

/// Integrate the baseline corrected trace for QDC in the range [0, pulseLength) and return the result.
double pmtResponse::integratePulse(){
	return this->integratePulse(0, pulseLength);
}

/// Integrate the baseline corrected trace for QDC in the range [maxIndex-start_, maxIndex+stop_] and return the result.
double pmtResponse::integratePulseFromMaximum(const short &start_, const short &stop_){
	size_t low = (maxIndex > start_ ? maxIndex-start_ : 0);
	return this->integratePulse(low, maxIndex+stop_);
}

/// Integrate the baseline corrected trace for QDC in the range [maxIndex-pulseIntegralLow, maxIndex+pulseIntegralHigh] and return the result.
double pmtResponse::integratePulseFromMaximum(){
	return this->integratePulseFromMaximum(pulseIntegralLow, pulseIntegralHigh);
}

/// Perform traditional CFD analysis on the waveform.
double pmtResponse::analyzeCFD(const double &F_/*=0.5*/, const size_t &D_/*=1*/, const size_t &L_/*=1*/){
	if(pulseLength == 0 || pulseArray.empty()) return -9999;
	double cfdMinimum = 9999;
	size_t cfdMinIndex = 0;
	
	double phase = -9999;
	double *cfdvals = new double[pulseLength];

	// Compute the cfd waveform.
	for(size_t cfdIndex = 0; cfdIndex < pulseLength; ++cfdIndex){
		cfdvals[cfdIndex] = 0.0;
		if(cfdIndex >= L_ + D_ - 1){
			for(size_t i = 0; i < L_; i++)
				cfdvals[cfdIndex] += F_ * (pulseArray[cfdIndex - i]-baseline) - (pulseArray[cfdIndex - i - D_]-baseline);
		}
		if(cfdvals[cfdIndex] < cfdMinimum){
			cfdMinimum = cfdvals[cfdIndex];
			cfdMinIndex = cfdIndex;
		}
	}

	// Find the zero-crossing.
	if(cfdMinIndex > 0){
		// Find the zero-crossing.
		for(size_t cfdIndex = cfdMinIndex-1; cfdIndex > 0; cfdIndex--){
			if(cfdvals[cfdIndex] >= 0.0 && cfdvals[cfdIndex+1] < 0.0){
				phase = cfdIndex - cfdvals[cfdIndex]/(cfdvals[cfdIndex+1]-cfdvals[cfdIndex]);
				break;
			}
		}
	}
	
	delete[] cfdvals;

	return (phase*adcClockTick-traceDelay+tLatch);
}


/// Perform polynomial CFD analysis on the waveform.
double pmtResponse::analyzePolyCFD(const double &F_){
	if(pulseLength == 0 || pulseArray.empty()) return -9999;
	if(maximum <= 0 && findMaximum() <= 0) return -9999;

	double threshold = F_*maximum + baseline;

	double phase = -9999;
	for(unsigned short cfdIndex = maxIndex; cfdIndex > 0; cfdIndex--){
		if(pulseArray[cfdIndex-1] < threshold && pulseArray[cfdIndex] >= threshold){
			if(false){
			calculateP2(cfdIndex-1, &pulseArray[cfdIndex-1], &cfdPar[4]);
			// Calculate the phase of the trace.
			if(cfdPar[6] != 0)
				phase = (-cfdPar[5]+std::sqrt(cfdPar[5]*cfdPar[5] - 4*cfdPar[6]*(cfdPar[4] - threshold)))/(2*cfdPar[6]);
			else
				phase = (threshold-cfdPar[4])/cfdPar[5];
			}
			else if (true){
				double x[2] = {(double)(cfdIndex-1), (double)cfdIndex};
				double p1 = (pulseArray[x[1]]-pulseArray[x[0]]) / (x[1]-x[0]);
				double p0 = pulseArray[x[1]]-p1*x[1];
			
				phase = ( threshold - p0 ) / p1;
			}
			break;
		}
	}

	return (phase*adcClockTick-traceDelay+tLatch);
}

/// Perform polynomial CFD analysis on the waveform using F=polyCfdFraction
double pmtResponse::analyzePolyCFD(){
	return this->analyzePolyCFD(polyCfdFraction);
}

void pmtResponse::copyTrace(unsigned short *arr, const size_t &len){
	if(!isDigitized) this->digitize();
	size_t stop = (len <= pulseLength ? len : pulseLength);
	for(size_t i = 0; i < stop; i++){
		arr[i] = pulseArray[i];
	}
}

/// Copy the digitized trace into a vector.
void pmtResponse::copyTrace(std::vector<unsigned short> &vec){
	vec.clear();
	vec.reserve(pulseLength);
	for(size_t i = 0; i < pulseLength; i++){
		vec.push_back(pulseArray[i]);
	}
}

void pmtResponse::clear(){
	for(size_t i = 0; i < 7; i++){
		cfdPar[i] = 0;
	}

	maximumTime = -9999;
	maximum = -9999;
	baseline = -9999;

	minimumArrivalTime = 1E6;

	arrivalTimes.clear();

	maxIndex = 0;
	
	isDigitized = false;
	isTriggered = false;
	for(size_t i = 0; i < pulseLength; i++)
		pulseArray[i] = 0;
}

void pmtResponse::print(){
	std::cout << "***********************************************************\n";
	if(functionType == EXPO){
		std::cout << "* response : exponential\n";
		std::cout << "* risetime : " << risetime << " ns" << std::endl;
		std::cout << "* falltime : " << falltime << " ns" <<  std::endl;
	}
	else if(functionType == VANDLE){
		std::cout << "* response : vandle\n";
		std::cout << "* beta     : " << risetime << " ns" << std::endl;
		std::cout << "* gamma    : " << falltime << " ns" <<  std::endl;		
	}
	else if(functionType == GAUSS){
		std::cout << "* response : gaussian\n";
		std::cout << "* sigma    : " << risetime << " ns" << std::endl;
	}
	std::cout << "* spread   : " << timeSpread << " ns" << std::endl;
	std::cout << "* delay    : " << traceDelay << " ns" << std::endl;
	std::cout << "* gain     : " << gain << "x" << std::endl;
	std::cout << "* baseline : " << baselineFraction*100 << "% (" << (int)(baselineFraction*adcBins) << " channels)" << std::endl;
	std::cout << "* jitter   : " << baselineJitterFraction*100 << "% (+-" << (int)(baselineJitterFraction*adcBins) << " channels)" << std::endl;
	std::cout << "* CfdF     : " << polyCfdFraction << std::endl;
	std::cout << "* sampling : " << getAdcClockFrequency() << " MSPS (" << adcClockTick << " ns)\n";
	std::cout << "* Ilow     : " << pulseIntegralLow << " clock ticks (" << pulseIntegralLow*adcClockTick << " ns)" << std::endl;
	std::cout << "* Ihigh    : " << pulseIntegralHigh << " clock ticks (" << pulseIntegralHigh*adcClockTick << " ns)" << std::endl;
	std::cout << "* range    : " << adcBins << " channels" << std::endl;
	std::cout << "* length   : " << pulseLength << " clock ticks (" << pulseLength*adcClockTick << " ns)" << std::endl;
	std::cout << "***********************************************************\n";	
}

void pmtResponse::printRaw(){
	/*for(size_t i = 0; i < pulseLength; i++){
		std::cout << i << "\t" << rawPulse[i] << std::endl;
	}*/
}

double pmtResponse::calculateP2(const short &x0, unsigned short *y, double *p){
	double x1[3], x2[3];
	for(size_t i = 0; i < 3; i++){
		x1[i] = (x0+i);
		x2[i] = std::pow(x0+i, 2);
	}

	double denom = (x1[1]*x2[2]-x2[1]*x1[2]) - x1[0]*(x2[2]-x2[1]*1) + x2[0]*(x1[2]-x1[1]*1);
	if (denom==0){
		std::cout << "calculateP2 denom zero" << std::endl;
		return -999;
		}
	p[0] = (y[0]*(x1[1]*x2[2]-x2[1]*x1[2]) - x1[0]*(y[1]*x2[2]-x2[1]*y[2]) + x2[0]*(y[1]*x1[2]-x1[1]*y[2]))/denom;
	p[1] = ((y[1]*x2[2]-x2[1]*y[2]) - y[0]*(x2[2]-x2[1]*1) + x2[0]*(y[2]-y[1]*1))/denom;
	p[2] = ((x1[1]*y[2]-y[1]*x1[2]) - x1[0]*(y[2]-y[1]*1) + y[0]*(x1[2]-x1[1]*1))/denom;
	if (p[2]==0){
		std::cout << "calculateP2 p2 zero" << std::endl;
		return -999;
		}
	// Calculate the maximum of the polynomial.
	return (p[0] - p[1]*p[1]/(4*p[2]));
}

double pmtResponse::calculateP3(const short &x, unsigned short *y, double *p, double &xmax){
	double x1[4], x2[4], x3[4];
	for(size_t i = 0; i < 4; i++){
		x1[i] = (x+i);
		x2[i] = std::pow(x+i, 2);
		x3[i] = std::pow(x+i, 3);
	}

	double denom = (x1[1]*(x2[2]*x3[3]-x2[3]*x3[2]) - x1[2]*(x2[1]*x3[3]-x2[3]*x3[1]) + x1[3]*(x2[1]*x3[2]-x2[2]*x3[1])) - (x1[0]*(x2[2]*x3[3]-x2[3]*x3[2]) - x1[2]*(x2[0]*x3[3]-x2[3]*x3[0]) + x1[3]*(x2[0]*x3[2]-x2[2]*x3[0])) + (x1[0]*(x2[1]*x3[3]-x2[3]*x3[1]) - x1[1]*(x2[0]*x3[3]-x2[3]*x3[0]) + x1[3]*(x2[0]*x3[1]-x2[1]*x3[0])) - (x1[0]*(x2[1]*x3[2]-x2[2]*x3[1]) - x1[1]*(x2[0]*x3[2]-x2[2]*x3[0]) + x1[2]*(x2[0]*x3[1]-x2[1]*x3[0]));
	if (denom==0){
		std::cout << "calculateP3 denom zero" << std::endl;
		return -999;
		}
	p[0] = (y[0]*(x1[1]*(x2[2]*x3[3]-x2[3]*x3[2]) - x1[2]*(x2[1]*x3[3]-x2[3]*x3[1]) + x1[3]*(x2[1]*x3[2]-x2[2]*x3[1])) - y[1]*(x1[0]*(x2[2]*x3[3]-x2[3]*x3[2]) - x1[2]*(x2[0]*x3[3]-x2[3]*x3[0]) + x1[3]*(x2[0]*x3[2]-x2[2]*x3[0])) + y[2]*(x1[0]*(x2[1]*x3[3]-x2[3]*x3[1]) - x1[1]*(x2[0]*x3[3]-x2[3]*x3[0]) + x1[3]*(x2[0]*x3[1]-x2[1]*x3[0])) - y[3]*(x1[0]*(x2[1]*x3[2]-x2[2]*x3[1]) - x1[1]*(x2[0]*x3[2]-x2[2]*x3[0]) + x1[2]*(x2[0]*x3[1]-x2[1]*x3[0]))) / denom;
	p[1] = ((y[1]*(x2[2]*x3[3]-x2[3]*x3[2]) - y[2]*(x2[1]*x3[3]-x2[3]*x3[1]) + y[3]*(x2[1]*x3[2]-x2[2]*x3[1])) - (y[0]*(x2[2]*x3[3]-x2[3]*x3[2]) - y[2]*(x2[0]*x3[3]-x2[3]*x3[0]) + y[3]*(x2[0]*x3[2]-x2[2]*x3[0])) + (y[0]*(x2[1]*x3[3]-x2[3]*x3[1]) - y[1]*(x2[0]*x3[3]-x2[3]*x3[0]) + y[3]*(x2[0]*x3[1]-x2[1]*x3[0])) - (y[0]*(x2[1]*x3[2]-x2[2]*x3[1]) - y[1]*(x2[0]*x3[2]-x2[2]*x3[0]) + y[2]*(x2[0]*x3[1]-x2[1]*x3[0]))) / denom;
	p[2] = ((x1[1]*(y[2]*x3[3]-y[3]*x3[2]) - x1[2]*(y[1]*x3[3]-y[3]*x3[1]) + x1[3]*(y[1]*x3[2]-y[2]*x3[1])) - (x1[0]*(y[2]*x3[3]-y[3]*x3[2]) - x1[2]*(y[0]*x3[3]-y[3]*x3[0]) + x1[3]*(y[0]*x3[2]-y[2]*x3[0])) + (x1[0]*(y[1]*x3[3]-y[3]*x3[1]) - x1[1]*(y[0]*x3[3]-y[3]*x3[0]) + x1[3]*(y[0]*x3[1]-y[1]*x3[0])) - (x1[0]*(y[1]*x3[2]-y[2]*x3[1]) - x1[1]*(y[0]*x3[2]-y[2]*x3[0]) + x1[2]*(y[0]*x3[1]-y[1]*x3[0]))) / denom;
	p[3] = ((x1[1]*(x2[2]*y[3]-x2[3]*y[2]) - x1[2]*(x2[1]*y[3]-x2[3]*y[1]) + x1[3]*(x2[1]*y[2]-x2[2]*y[1])) - (x1[0]*(x2[2]*y[3]-x2[3]*y[2]) - x1[2]*(x2[0]*y[3]-x2[3]*y[0]) + x1[3]*(x2[0]*y[2]-x2[2]*y[0])) + (x1[0]*(x2[1]*y[3]-x2[3]*y[1]) - x1[1]*(x2[0]*y[3]-x2[3]*y[0]) + x1[3]*(x2[0]*y[1]-x2[1]*y[0])) - (x1[0]*(x2[1]*y[2]-x2[2]*y[1]) - x1[1]*(x2[0]*y[2]-x2[2]*y[0]) + x1[2]*(x2[0]*y[1]-x2[1]*y[0]))) / denom;

	if(p[3] == 0){
		if (p[2]==0){
			std::cout << "calculateP3 p2 zero" << std::endl;
			return -999;
		}
		// Handle the case of p[3] == 0.
		return (p[0] - p[1]*p[1]/(4*p[2]));
	}
	// Calculate the maximum of the polynomial.
	if (4*p[2]*p[2]-12*p[3]*p[1]<0){
			std::cout << "calculateP3 max calc error" << std::endl;
			return -999;
		}
	double xmax1 = (-2*p[2]+std::sqrt(4*p[2]*p[2]-12*p[3]*p[1]))/(6*p[3]);
	double xmax2 = (-2*p[2]-std::sqrt(4*p[2]*p[2]-12*p[3]*p[1]))/(6*p[3]);

	if((2*p[2]+6*p[3]*xmax1) < 0) // The second derivative is negative (i.e. this is a maximum).
		xmax = xmax1;
	else
		xmax = xmax2;
	
	return (p[0] + p[1]*xmax + p[2]*xmax*xmax + p[3]*xmax*xmax*xmax);
}

double pmtResponse::func(const double &t, const double &dt/*=0*/) const {
	if(functionType == EXPO){ // Standard single photon response function.
		if(t-dt <= 0) return 0;
		return gain*(1/(falltime-risetime))*(std::exp(-(t-dt)/falltime)-std::exp(-(t-dt)/risetime));
	}
	else if(functionType == VANDLE){ // "Vandle" function. The risetime is BETA and the decay time is GAMMA.
		if(t-dt <= 0) return 0;
		return (100*gain*std::exp((dt-t)*risetime)*(1 - std::exp(-std::pow((t-dt)*falltime,4))));
	}
	else if(functionType == GAUSS){ // Normalized gaussian function. The risetime is sigma and the decay time is not used.
		return gain*((1/(risetime*sqrt2pi))*std::exp(-0.5*std::pow((t-dt)/risetime, 2.0)));
	}
	return 0;
}

double pmtResponse::findMaximum(){
	if(pulseLength == 0 || pulseArray.empty()) return -9999;

	// Find the baseline.
	double tempbaseline = 0.0;
	size_t sample_size = floor(traceDelay/adcClockTick);
	sample_size = (sample_size <= pulseLength ? 15:pulseLength);
	for(size_t i = 0; i < sample_size; i++){
		tempbaseline += pulseArray[i];
	}
	tempbaseline = tempbaseline/sample_size;
	baseline = double(tempbaseline);

	// Find the maximum ADC value and the maximum bin.
	unsigned short maxADC = 0;
	for(size_t i = 0; i < pulseLength; i++){
		if(pulseArray[i] - baseline > maxADC){ 
			maxADC = pulseArray[i] - baseline;
			maxIndex = i;
		}
	}

	// Find the pulse maximum by fitting with a third order polynomial.
	if(pulseArray[maxIndex-1] >= pulseArray[maxIndex+1]) // Favor the left side of the pulse.
		maximum = calculateP3(maxIndex-2, &pulseArray[maxIndex-2], cfdPar, maximumTime) - baseline;
	else // Favor the right side of the pulse.
		maximum = calculateP3(maxIndex-1, &pulseArray[maxIndex-1], cfdPar, maximumTime) - baseline;

	return maximum;
}
