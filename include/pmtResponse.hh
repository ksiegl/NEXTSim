#ifndef PMT_RESPONSE_HH
#define PMT_RESPONSE_HH

#include <vector>

/** @class spectralResponse
  * @brief Class used to interpolate PMT anode quantum efficiency from an input spectrum.
  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu)
  * @date May 23, 2019
  */

class spectralResponse{
  public:
	/** Default constructor
	  */
	spectralResponse() : xmin(0), xmax(0), size(0) { }

	/** Destructor
	  */
	~spectralResponse();

	/** Clone this object
	  */
	spectralResponse clone() const ;

	/** Get the range of the quantum efficiency spectrum
	  * @param min_ The minimum wavelength defined in the spectrum (in nm)
	  * @param max_ The maximum wavelength defined in the spectrum (in nm)
	  */
	void getRange(double &min_, double &max_) const ;

	/** Get a vector containing all wavelengths contained in the spectrum (in nm)
	  */
	std::vector<double> getWavelength() const { return wavelength; }
	
	/** Get a vector containing all quantum efficiencies contained in the spectrum (in percent)
	  */
	std::vector<double> getPercentage() const { return percentage; }

	/** Get the number of elements in the spectrum
	  */
	size_t getSize() const { return size; }

	/** Copy the quantum efficiency spectrum from another spectralResponse object
	  * @param other Pointer to another spectralResponse from which the spectrum will be copied
	  */
	void copy(const spectralResponse *other);

	/** Load the quantum efficiency spectrum from a root file
	  * @param fname Path to a root file containing a TGraph named 'spec' or to an ascii file with two columns of values
	  * @return True if all objects are loaded successfully and return false otherwise
	  */
	bool load(const std::string &fname);

	/** Interpolate the quantum efficiency for a given wavelength
	  * @param wavelength Optical photon wavelength (in nm)
	  * @return Return the quantum efficiency for a given wavelength
	  */
	double eval(const double &wavelength) const ;

  private:
	double xmin; ///< Minimum optical photon wavelenth defined in the quantum efficiency spectrum (in nm)
	double xmax; ///< Maximum optical photon wavelenth defined in the quantum efficiency spectrum (in nm)
	
	size_t size; ///< The number of elements in the spectrum

	std::vector<double> wavelength; ///< The vector of all wavelength values (in nm)
	std::vector<double> percentage; ///< The vector of all quantum efficiencies (in percent)

	/** Delete the spectrum and clear all values
	  */	
	void clear();
};

/** @class pmtResponse
  * @brief Class used to simulate the light response due to detection of multiple optical photons.
  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu)
  * @date May 23, 2019
  */

class pmtResponse{
public:
	/** @class photonArrivalTime
	  * @brief Simple structure to contain an optical photon's arrival time, its gain, and its time offset
	  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu)
	  * @date September 19, 2019
	  */
	class photonArrivalTime{
	public:
		double time; ///< Optical photon arrival time (in ns)
		double gain; ///< The gain of the single-photon response
		double dt; ///< The offset of the response function due to PMT time spread (in ns)

		/** Default constructor
		  */
		photonArrivalTime() : time(0), gain(1), dt(0) { }
	
		/** Arrival time constructor
		  */
		photonArrivalTime(const double &time_, const double &gain_, const double &dt_) : time(time_), gain(gain_), dt(dt_) { }
	};

	/** Types of single-photon light response functions
	  */
	enum photonResponseType {EXPO, VANDLE, GAUSS};

	/** Default constructor
	  */
	pmtResponse();
	
	/** Pulse constructor
	  * @param risetime_ The rise time of the single photon response function (in ns)
	  * @param falltime_ The decay time of the single photon response function (in ns)
	  */
	pmtResponse(const double &risetime_, const double &falltime_);

	/** Destructor
	  */
	~pmtResponse();

	/** Clone this object
	  */
	pmtResponse clone() const ;

	/** Get the length of the light pulse (in ADC clock ticks)
	  */
	size_t getPulseLength() const { return pulseLength; }

	/** Get the maximum value of the pulse (in ADC channels)
	  */
	double getMaximum() const { return maximum; }

	/** Get the baseline of the pulse (in ADC channels)
	  */
	double getBaseline() const { return baseline; }

	/** Get the pulse phase at the maximum (in ns)
	  */
	double getMaximumTime() const { return maximumTime; }

	/** Get the fast trigger from the pmt signal
	  */
	bool getTrigger() const {return isTriggered; }
	
	/** Get a pointer to the array of Poly CFD parameters
	  * @return Pointer to an array of 7 doubles 
	  */
	double *getCfdParams() { return cfdPar; }

	/** Get the maximum ADC index of the pulse (in ADC clock ticks)
	  */
	unsigned short getMaximumIndex() const { return maxIndex; }

	/** Get a pointer to the digitized pulse array
	  */
	unsigned short *getDigitizedPulse(){ return pulseArray.data(); }

	/** Return true if the digitized pulse is saturated and return false otherwise
	  */
	bool getPulseIsSaturated() const { return pulseIsSaturated; }

	/** Get the gain-weighted arrival time of the photon pulse (in ns)
	  */
	double getWeightedPhotonArrivalTime() const ;

	/** Get the total light response spectrum by sampling the single-photon response
	  * spectra of each detected optical photon
	  */
	void getRawPulse(std::vector<double> &rawPulse) const ;

	/** Get the minimum photon arrival time of the pulse (in ns)
	  */
	double getMinimumPhotonArrivalTime() const { return minimumArrivalTime; }

	/** Return true if the digitized pulse is to be printed and return false otherwise
	  */
	bool getPrintTrace() const { return printTrace; }

	/** Get the period of the ADC clock (in ns)
	  */
	double getAdcClockInNanoseconds() const { return adcClockTick; }

	/** Get the frequency of the ADC clock (in MHz)
	  */
	double getAdcClockFrequency() const { return ((1E9/adcClockTick)/1E6); }

	/** Return true if the anode quantum efficiency will be used, and return false otherwise
	  */
	bool getSpectralResponseEnabled() const { return useSpectralResponse; }

	/** Get a pointer to the PMT quantum efficiency spectrum
	  */
	spectralResponse* getSpectralResponse(){ return &spec; }

	/** Get a const pointer to the PMT quantum efficiency spectrum
	  */
	const spectralResponse* getConstSpectralResponse() const { return (const spectralResponse*)(&spec); }	

	/** Set the rise time of the single photon pulse
	  * @param risetime_ The rise time of the single photon response function (in ns)
	  */
	void setRisetime(const double &risetime_);
	
	/** Set the decay time of the single photon pulse
	  * @param falltime_ The decaytime of the single photon response function (in ns)
	  */
	void setFalltime(const double &falltime_);

	/** Set the FWHM spread in the photo-electron transit time
	  * @param spread_ The uncertainty spread in the arrival time of the single photon response (in ns)
	  */
	void setTransitTimeSpread(const double &spread_){ timeSpread = spread_; }

	/** Set the time delay of the pulse
	  * @param traceDelay_ The offset in the arrival time of the single photon response (in ns)
	  */
	void setTraceDelay(const double &traceDelay_){ traceDelay = traceDelay_; }

	/** Set the gain
	  * @param gain_ The multiplicitive gain of the single photon response
	  */
	void setGain(const double &gain_){ gain = gain_; }

	/** Set the length of the pulse
	  * @param len_ The length of the light response pulse (in ADC clock ticks)
	  */
	void setPulseLength(const size_t &len);

	/** Set the length of the pulse 
	  * @param length The length of the light response pulse (in ns)
	  */
	void setPulseLengthInNanoSeconds(const double &length);

	/** Set the dynamic bit range of the ADC
	  * @param len_ The number of bits of dynamic range the ADC has
	  */
	void setBitRange(const size_t &len);

	/** Set the single photon response function type. The available function types with their usage of the fall time and rise time 
	  * parameters (set by setFalltime() and setRisetime() respectively) are described in the table below
	  *
	  * | Type           | Description | Fall Time | Rise Time |
	  * |----------------|-------------|-----------|-----------|
	  * | EXPO (default) | [Double exponential function](https://iopscience.iop.org/article/10.1088/0031-9155/54/21/004) | Decay time of pulse (in ns) | Rise time of pulse (in ns) |
	  * | VANDLE         | [Vandle function](https://www.sciencedirect.com/science/article/pii/S0168900213015672) | Gamma (in ns) | Beta (in ns) |
	  * | GAUSS          | Normalized gaussian | Not used | Sigma (in ns) |
	  *
	  * @param type The single photon response function to use
	  */
	void setFunctionType(const photonResponseType &type){ functionType = type; }

	/** Set the height of the digitized pulse baseline
	  * @param percentage The amplitude of the pulse baseline (percentage of the total ADC dynamic range)
	  */	
	void setBaselinePercentage(const double &percentage){ baselineFraction = percentage/100; }

	/** Set the noise in the digitized pulse baseline
	  * @param percentage The maximum amplitude of noise above or below the true pulse baseline (percentage of the total ADC dynamic range)
	  */	
	void setBaselineJitterPercentage(const double &percentage){ baselineJitterFraction = percentage/100; }

	/** Set the PolyCFD fraction parameter (F)
	  * @param frac The fraction of the baseline corrected pulse height where the PolyCFD crossing point will be extracted
	  */
	void setPolyCfdFraction(const double &frac){ polyCfdFraction = frac; }
	
	/** Set the lower limit for pulse integral TQDC 
	  * @param low The number of ADC clock ticks to the left of the peak maximum to start TQDC integral (inclusive)
	  */
	void setPulseIntegralLow(const short &low){ pulseIntegralLow = low; }

	/** Set the upper limit for pulse integral TQDC
	  * @param low The number of ADC clock ticks to the right of the peak maximum to stop TQDC integral (inclusive)
	  */
	void setPulseIntegralHigh(const short &high){ pulseIntegralHigh = high; }

	/** Set the flag which indicates whether or not the digitized pulse should be printed
	  */
	void setPrintTrace(const bool &enabled){ printTrace = enabled; }

	/** Set the period of the ADC clock (in ns)
	  */
	void setAdcClockInNanoseconds(const double &period){ adcClockTick = period; }

	/** Set the frequency of the ADC clock (in MHz)
	  * @return The period of the ADC clock (in ns)
	  */
	double setAdcClockFrequency(const double &frequency){ return (adcClockTick = (1/(frequency*1E6))*1E9); }

	/** Set the ADC latching time (in clock ticks)
	  * @return The latching time (in ns)
	  */
	double setAdcLatchTicks(const double &latch){ return (tLatch = adcClockTick*latch); }

	/** Disable use of PMT quantum efficiency spectrum
	  */
	void disableSpectralResponse(){ useSpectralResponse = false; }

	/** Load PMT spectral response from root file
	  * @param fname Path to a root file containing the anode quantum efficiency spectrum
	  * @return True if the spectrum is loaded successfully and return false otherwise
	  */
	bool loadSpectralResponse(const char *fname);

	/** Copy quantum efficiency spectrum from another pmtResponse object
	  * @param spec_ Pointer to a spectralResponse from which the spectrum will be copied
	  */
	void copySpectralResponse(const spectralResponse *spec_){ spec.copy(spec_); }

	/** Add a photon signal to the raw pulse
	  * @param arrival Arrival time of the optical photon (in ns)
	  * @param wavelength Wavelength of the optical photon (in nm)
	  * @param gain_ Additional multiplicitive gain to use for single photon response
	  */
	void addPhoton(const double &arrival, const double &wavelength=0, const double &gain_=1);

	/** Sample the total photon response by iterating over the list of single-photon arrival times
	  * @param time The time at which the full photon light response will be sampled (ns)
	  * @return The full photon light response at the specified time
	  */
	double sample(const double &time) const ;

	/** Build the raw light response pulse by stepping through the trace and sampling the sum of the 
	  * response of all detected optical photons and then digitize the spectrum
	  * @param baseline The light pulse baseline as a fraction of the full ADC dynamic range
	  * @param jitter The random jitter in the light pulse baseline as a fraction of the full ADC dynamic range
	  */
	void digitize(const double &baseline, const double &jitter);

	/** Digitize the raw light pulse (uses baseline=baselineFraction and jitter=baselineJitter)
	  */
	void digitize();

	/** Integrate the light pulse in a range and return the result
	  * @param start_ Inclusive start index for integration (in ADC clock ticks)
	  * @param stop_ Inclusive stop index for integration (in ADC clock ticks)
	  * @return The integral of the baseline corrected light pulse in the range [start_, stop_]
	  */
	double integratePulse(const size_t &start_, const size_t &stop_);

	/** Integrate the light pulse in a range and return the result
	  * @return The integral of the baseline corrected light pulse in the range [0, pulseLength)
	  */
	double integratePulse();

	/** Integrate the light pulse in a range about the maximum and return the result
	  * @param start_ The number of ADC clock ticks to the left of the light pulse max index to begin integration
	  * @param stop_ The number of ADC clock ticks to the right of the light pulse max index to stop integration
	  * @return The integral of the baseline corrected light pulse in the range [maxIndex-start_, maxIndex+stop_]
	  */
	double integratePulseFromMaximum(const short &start_, const short &stop_);

	/** Integrate the light pulse in a range about the maximum and return the result
	  * @return The integral of the baseline corrected light pulse in the range [maxIndex-pulseIntegralLow, maxIndex+pulseIntegralHigh]
	  */
	double integratePulseFromMaximum();

	/** Perform traditional CFD analysis on the waveform using [XIA CFD Algorithm](https://www.xia.com/Papers/TimeRes_DigConstFracDiscrimination.pdf)
	  * @param F_ The CFD crossing point as a function of the pulse height
	  * @param D_ CFD delay length (in ADC clock ticks)
	  * @param L_ Running average length (in ADC clock ticks)
	  * @return The phase of the pulse, i.e. the CFD crossing point (in ns)
	  */
	double analyzeCFD(const double &F_=0.5, const size_t &D_=1, const size_t &L_=1);
	
	/** Perform polynomial CFD analysis (PolyCFD) on the waveform
	  * @param F_ The PolyCFD crossing point as a function of the pulse height
	  * @return The phase of the pulse, i.e. the PolyCFD crossing point (in ns)
	  */
	double analyzePolyCFD(const double &F_);

	/** Perform polynomial CFD analysis (PolyCFD) on the waveform using F=polyCfdFraction
	  * @return The phase of the pulse, i.e. the PolyCFD crossing point (in ns)
	  */
	double analyzePolyCFD();

	/** Copy the digitized trace into an array
	  * @param arr The array of unsigned shorts into which the digitized trace will be copied
	  * @param len The length of the array pointed at by arr
	  */
	void copyTrace(unsigned short *arr, const size_t &len);

	/** Copy the digitized trace into a vector
	  * @param vec The vector of unsigned shorts into which the digitized trace will be copied
	  */
	void copyTrace(std::vector<unsigned short> &vec);

	/** Clear the ADC trace and all analysis variables
	  */
	void clear();

	/** Print pulse parameters to stdout
	  */
	void print();

	/** Print the raw pulse to stdout
	  */
	void printRaw();

	/** Calculate the parameters for a second order polynomial which passes through 3 points
	  * @param x0 Initial x value. Sequential x values are assumed to be x0, x0+1, and x0+2
	  * @param y Pointer to the beginning of the array of unsigned shorts containing the three y values
	  * @param p Pointer to the array of doubles for storing the three polynomial parameters
	  * @return The maximum/minimum of the polynomial
	  */
	static double calculateP2(const short &x0, unsigned short *y, double *p);

	/** Calculate the parameters for a third order polynomial which passes through 4 points
	  * @param x0 Initial x value. Sequential x values are assumed to be x0, x0+1, x0+2, and x0+3
	  * @param y Pointer to the beginning of the array of unsigned shorts containing the four y values
	  * @param p Pointer to the array of doubles for storing the three polynomial parameters
	  * @param xmax  X value at the peak of the pulse
	  * @return The local maximum/minimum of the polynomial
	  */
	static double calculateP3(const short &x0, unsigned short *y, double *p, double &xmax);

private:
	double cfdPar[7]; ///< Cfd fitting parameters.	

	double risetime; ///< Risetime of the single photon response function
	double falltime; ///< Decay time of the single photon response function
	double timeSpread; ///< Spread in the single photon response time offset based on the photo-electron transit time spread (in ns)
	double traceDelay; ///< Time offset of single photon response function (in ns)
	double gain; ///< Multiplicitive gain of the single photon response function

	double maximumTime; ///< The time at the trace maximum (in ns)
	double maximum; ///< The maximum value of the trace
	double baseline; ///< Baseline offset of the trace

	double baselineFraction; ///< Digitized light pulse baseline as a fraction of the total ADC dynamic range
	double baselineJitterFraction; ///< Digitized light pulse baseline jitter as a fraction of the total ADC dynamic range
	double polyCfdFraction; ///< PolyCFD F parameter as a fraction of max pulse height minus the baseline

	double adcClockTick; ///< The period of the ADC clock (in ns)
	double tLatch; ///< ADC clock latching time (in ns)
	
	short pulseIntegralLow; ///< Number of ADC clock ticks for start of digitized pulse TQDC integral wrt pulse maximum
	short pulseIntegralHigh; ///< Number of ADC clock ticks for stop of digitized pulse TQDC integral wrt pulse maximum

	unsigned short maxIndex; ///< The index of the maximum trace bin (in ADC clock ticks)

	unsigned int adcBins; ///< Dynamic range of ADC in discrete channels

	size_t pulseLength; ///< Length of the raw and digitized ADC light pulses (in ADC clock ticks)
	
	bool isDigitized; ///< Flag indicating that the light response pulse has been digitized
	bool isTriggered; ///< Flag to indicate if signal would trigger a filter
	bool useSpectralResponse; ///< Flag indicating that PMT anode quantum efficiency spectrum should be used to compute anode gain
	bool pulseIsSaturated; ///< Flag indicating that the digitized pulse has exceeded the maximum ADC dynamic range
	bool printTrace; ///< Flag indicating that the left and right digitized PMT traces will be printed to the screen
	
	std::vector<unsigned short> pulseArray; ///< Array to store digitized light response pulse

	spectralResponse spec; ///< Anode quantum efficiency

	double minimumArrivalTime; ///< Minimum photon arrival time

	photonResponseType functionType; ///< Integer indicating the single photon response function to use to build the light response pulse

	std::vector<photonArrivalTime> arrivalTimes; ///< Vector of all optical photon arrival times and their individual single-photon response gains
	
	/** Evaluate the single-photon-response function for a given time and time offset
	  * @param t The time to along the pulse at which to evaluate the single-photon-response function (in ns)
	  * @param dt The offset of the pulse along the time axis (in ns)
	  * @return The value of the single-photon-response function at the user specified time
	  */
	double func(const double &t, const double &dt=0) const ;
	
	/** Compute the baseline and maximum of the light pulse
	  * @return The maximum of the light pulse if the array is properly initialized and return -9999 otherwise
	  */
	double findMaximum();
};

#endif
