#include "AnitaEventSummary.h"

#include "RawAnitaHeader.h" 
#include "UsefulAdu5Pat.h"
#include "TruthAnitaEvent.h"
#include "TBuffer.h"
#include "TClass.h"
//are these even necessary anymore? Who knows! 
//ClassImp(AnitaEventSummary)
//ClassImp(AnitaEventSummary::PointingHypothesis)
//ClassImp(AnitaEventSummary::SourceHypothesis)
  

const double C_IN_M_NS = 0.299792;


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaEventSummary::AnitaEventSummary()
    : anitaLocation()
{
  zeroInternals();
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * Takes care of copying the header info into the event summary
 */
AnitaEventSummary::AnitaEventSummary(const RawAnitaHeader* header)
    : anitaLocation()
{
  zeroInternals();

  setTriggerInfomation(header);
  eventNumber = header->eventNumber;
  run = header->run;
  realTime = header->realTime;

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * Takes care of copying the header and GPS info into the event summary
 */
AnitaEventSummary::AnitaEventSummary(const RawAnitaHeader* header, UsefulAdu5Pat* pat, const TruthAnitaEvent * truth)
    : anitaLocation(dynamic_cast<Adu5Pat*>(pat))
{

  zeroInternals();
  setTriggerInfomation(header);
  setSourceInformation(pat,truth);
  eventNumber = header->eventNumber;
  run = header->run;
  realTime = header->realTime;
}








/** 
 * Set everything to zero, should be called in constructor
 */
void AnitaEventSummary::zeroInternals(){

  fLastEventNumber = -1; 
  run = 0;
  eventNumber = 0;
  resetNonPersistent();


  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    nPeaks[polInd] = 0; 
    for(Int_t dir=0; dir < maxDirectionsPerPol; dir++)
    {
      memset(&peak[polInd][dir],0,sizeof(PointingHypothesis)); 
      memset(&coherent[polInd][dir],0,sizeof(WaveformInfo)); 
      memset(&deconvolved[polInd][dir],0,sizeof(WaveformInfo)); 
      memset(&coherent_filtered[polInd][dir],0,sizeof(WaveformInfo)); 
      memset(&deconvolved_filtered[polInd][dir],0,sizeof(WaveformInfo));

    }
    for(Int_t ant=0; ant < 48; ant++){
      memset(&channels[polInd][ant],0,sizeof(ChannelInfo));
      channels[polInd][ant].pol = pol;
      channels[polInd][ant].ant = ant;
    }
  }
  memset(&flags,0, sizeof(EventFlags)); 
  flags.pulser = EventFlags::NONE; 

  sun.reset();
  wais.reset(); 
  ldb.reset(); 
  mc.reset(); 

}



/** 
 * Utility function to get the polarisation of the highest map peak
 * Useful for TTree::Draw() 
 * 
 * @return the polarisation of the largest interferometric map peak
 */
AnitaPol::AnitaPol_t AnitaEventSummary::highestPol() const{
  findHighestPeak();
  return fHighestPol;
}


/** 
 * Utility function to get the polarisation of the highest map peak
 * Useful for TTree::Draw()
 * 
 * @return the polarisation of the largest interferometric map peak
 */
int AnitaEventSummary::highestPolAsInt() const{
  findHighestPeak();
  return int(fHighestPol);
}


/** 
 * Utility function to get the index of the highest map peak
 * Useful for TTree::Draw()
 * 
 * @return the index of the largest interferometric map peak
 */
int AnitaEventSummary::highestPeakInd() const{
  findHighestPeak();
  return fHighestPeakIndex;
}



/** 
 * Utility function to return a const reference to the higher map peak. 
 * Useful for TTree::Draw() 
 * 
 * @return the peak with largest value
 */
const AnitaEventSummary::PointingHypothesis& AnitaEventSummary::highestPeak() const{
  findHighestPeak();
  return peak[fHighestPol][fHighestPeakIndex];
}


/** 
 * Utility function to return a const reference to the the unfiltered coherently summed waveform info of the polarisation of the highest map peak
 * Useful for TTree::Draw() 
 * 
 * @return the unfiltered coherently summed waveform info corresponding to the highest map peak
 */
const AnitaEventSummary::WaveformInfo& AnitaEventSummary::highestCoherent() const{
  findHighestPeak();
  return coherent[fHighestPol][fHighestPeakIndex];
}


/** 
 * Utility function to return a const reference to the the unfiltered, deconvolved coherently summed waveform info of the polarisation of the highest map peak
 * Useful for TTree::Draw() 
 * 
 * @return the unfiltered, deconvolved coherently summed waveform info corresponding to the highest map peak
 */
const AnitaEventSummary::WaveformInfo& AnitaEventSummary::highestDeconvolved() const{
  findHighestPeak();
  return deconvolved[fHighestPol][fHighestPeakIndex];
}




/** 
 * Utility function to return a const reference to the the filtered coherently summed waveform info of the polarisation of the highest map peak
 * Useful for TTree::Draw() 
 * 
 * @return the filtered coherently summed waveform info corresponding to the highest map peak
 */
const AnitaEventSummary::WaveformInfo& AnitaEventSummary::highestCoherentFiltered() const{
  findHighestPeak();
  return coherent_filtered[fHighestPol][fHighestPeakIndex];
}


/** 
 * Utility function to return a const reference to the the filtered, deconvolved coherently summed waveform info of the polarisation of the highest map peak
 * Useful for TTree::Draw() 
 * 
 * @return the filtered, deconvolved coherently summed waveform info corresponding to the highest map peak
 */
const AnitaEventSummary::WaveformInfo& AnitaEventSummary::highestDeconvolvedFiltered() const{
  findHighestPeak();
  return deconvolved_filtered[fHighestPol][fHighestPeakIndex];
}






/** 
 * Set trigger information in EventFlags directly from the header
 * 
 * @param header is a pointer to the event header
 */
void AnitaEventSummary::setTriggerInfomation(const RawAnitaHeader* header){  

  // setting to zero, for testing this function.
  flags.isVPolTrigger = 0;
  flags.isHPolTrigger = 0;
  
  // need two adjacent L3 phi-sectors to trigger ANITA-3
  for(Int_t phi=0; phi<NUM_PHI; phi++){
    Int_t phi2 = (phi+1)%NUM_PHI; // adjacent phi-sector

    Int_t trigBit = (header->l3TrigPattern >> phi) & 1;
    Int_t trigBit2 = (header->l3TrigPattern >> phi2) & 1;

    // std::cout << header->l3TrigPattern << "\t"
    // 	      << phi << "\t" << phi2 << "\t"
    // 	      << trigBit << "\t" << trigBit2 << std::endl;

    if(trigBit > 0 && trigBit2 > 0){
      flags.isVPolTrigger = 1;      
    }

    Int_t trigBitH = (header->l3TrigPatternH >> phi) & 1;
    Int_t trigBitH2 = (header->l3TrigPatternH >> phi2) & 1;

    if(trigBitH > 0 && trigBitH2 > 0){
      flags.isHPolTrigger = 1;      
    }
  }

  flags.isRF = header->getTriggerBitRF();
  flags.isSoftwareTrigger = header->getTriggerBitADU5() | header->getTriggerBitG12() | header->getTriggerBitSoftExt();
  flags.isAdu5Trigger = header->getTriggerBitADU5();
  flags.isG12Trigger = header->getTriggerBitG12();
  flags.isSoftwareTrigger = header->getTriggerBitSoftExt();
  flags.isMinBiasTrigger = flags.isAdu5Trigger || flags.isG12Trigger || header->getTriggerBitG12() || flags.isSoftwareTrigger;

  
}




/** 
 * Set the source information using the GPS info (and MC truth if non-NULL)
 * 
 * @param pat is a pointer to the event GPS information
 * @param truth is a pointer to the MC Truth, default value is NULL
 */
void AnitaEventSummary::setSourceInformation(UsefulAdu5Pat* pat, const TruthAnitaEvent * truth){


  pat->getSunPosition(sun.phi, sun.theta);
  sun.distance = 150e9;  // I guess in theory we could compute this! 

  pat->getThetaAndPhiWaveWaisDivide(wais.theta, wais.phi);
  wais.theta *= 180/ TMath::Pi(); 
  wais.phi *= 180/ TMath::Pi(); 
  wais.distance = pat->getWaisDivideTriggerTimeNs() * C_IN_M_NS; 

  pat->getThetaAndPhiWaveLDB(ldb.theta, ldb.phi);
  ldb.distance = pat->getLDBTriggerTimeNs() * C_IN_M_NS; 
  ldb.theta *= 180/ TMath::Pi(); 
  ldb.phi *= 180/ TMath::Pi(); 



  if (truth) 
  {
    pat->getThetaAndPhiWave(truth->sourceLon, truth->sourceLat, truth->sourceAlt, mc.theta,mc.phi);
    mc.theta*=TMath::RadToDeg();
    mc.phi*=TMath::RadToDeg();
    mc.weight = truth->weight; 
    mc.distance = pat->getTriggerTimeNsFromSource(truth->sourceLat, truth->sourceLon, truth->sourceAlt);
    mc.energy = truth->nuMom; // I guess this won't be true if icemc ever simulates non-relativistic neutrinos :P
  }
}



/** 
 * Set default values for MC truth
 */
void AnitaEventSummary::MCTruth::reset()
{
  SourceHypothesis::reset();
  memset(&wf[0],0,sizeof(WaveformInfo));
  memset(&wf[1],0,sizeof(WaveformInfo));
  weight = 0;
  energy = 0;
}





/** 
 * Utility function to return the linear polarization fraction so I can stop typing it
 * Useful for TTree::Draw() 
 * 
 * 
 * @return the linear polarization fraction
 */

double AnitaEventSummary::WaveformInfo::linearPolFrac() const {

  double value = TMath::Sqrt( pow(Q,2) + pow(U,2) ) / I;

  return value;
}


/** 
 * Utility function to return the linear polarization angle so I can stop typing it
 * Useful for TTree::Draw() 
 * 
 * 
 * @return the linear polarization angle is degrees
 */

double AnitaEventSummary::WaveformInfo::linearPolAngle() const {
  
  double value = (TMath::ATan(U/Q)/2)*TMath::RadToDeg();

  return value;

}



/** 
 * Utility function to return the circular polarization fraction so I can stop typing it
 * Useful for TTree::Draw() 
 * 
 * 
 * @return the circular polarization fraction
 */

double AnitaEventSummary::WaveformInfo::circPolFrac() const {
  
  double value = TMath::Abs(V)/I;

  return value;

}

/** 
 * Utility function to return the total polarization fraction so I can stop typing it
 * Useful for TTree::Draw() 
 * 
 * 
 * @return the circular polarization fraction
 */

double AnitaEventSummary::WaveformInfo::totalPolFrac() const {
  
  double value = TMath::Sqrt(pow(Q,2) + pow(U,2) + pow(V,2))/I;

  return value;

}




/** 
 * Return the standardized values of the peakMoments.
 * See https://en.wikipedia.org/wiki/Standardized_moment
 * 
 * @param k is the degree (as shown on wikipedia)
 * k = 1 would be zero if the origin of the moments was the mean
 * k = 2 is 1 is 1, definitionally
 * k = 3 is a measure of skewness about the origin
 * k = 4 is a measure of kurtosis about the origin
 * @return the normalized peak moment
 */
double AnitaEventSummary::WaveformInfo::standardizedPeakMoment(int k) const {

  // the first moment k=1, is stored in peakMoments[0]...  
  if(k <= 0){
    return -1;
  }
  else if (k <= 5){
    return peakMoments[k-1]/pow(peakMoments[1], 0.5*double(k));
  }
  else{
    return -1;
  }
}



/** 
 * Get the phi in degree of the current channel
 */
double AnitaEventSummary::ChannelInfo::getPhi() const {
  return AnitaGeomTool::Instance()->getPhiFromAnt(ant);
}





/** 
 * Make the source hypothesis go back to nonsense thats easy to recognize
 */
void AnitaEventSummary::SourceHypothesis::reset() {

  theta = -999;
  phi = -999;
  distance = -999;

  memset(mapHistoryVal,0,NUM_POLS*sizeof(Double_t));
  memset(mapValue,0,NUM_POLS*sizeof(Double_t));
}
  



/** 
 * Let the constructor do the hard work.
 * 
 * @param pat is ANITA's gps data
 * 
 */
AnitaEventSummary::PayloadLocation::PayloadLocation(const Adu5Pat* pat){
  update(pat);
}


/** 
 * Copy all the values of interest from the Adu5Pat
 * 
 * @param pat is ANITA's gps data
 */
void AnitaEventSummary::PayloadLocation::update(const Adu5Pat* pat){
  if(pat){
    latitude = pat->latitude;
    longitude = pat->longitude;
    altitude = pat->altitude;
    
    prevHeading = heading;
    heading = pat->heading;
    
  }
  else{
    reset();
  }
}


/** 
 * Get the angle between this peak and the source
 * 
 * @param source is a SourceHypothesis for the event (e.g. WAIS)
 * 
 * @return the phi angle between the selected peak and source
 */
double AnitaEventSummary::PointingHypothesis::dPhiSource(const SourceHypothesis& source) const{
  return dPhi(source.phi);
}


/** 
 * Get the angle between this peak and an arbitrary phi
 * 
 * @return the phi angle between the selected peak and source
 */
double AnitaEventSummary::PointingHypothesis::dPhi(double phi2) const{
  double dPhi = phi - phi2;
  if(dPhi < -180){
    dPhi += 360;
  }
  else if(dPhi >= 180){
    dPhi -= 360;
  }

  if(dPhi < -180 || dPhi >= 180){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " dPhi = " << dPhi
              << " is outside expected range, peak phi = " << phi << ", source phi = "
              << phi2 << std::endl;
  }

  return dPhi;
}




/** 
 * Use the summary information to get the difference in angle between two peak thetas.
 * 
 * Note:
 * The sign convention for theta varies between ANITA libraries for silly historical reasons.
 * Acclaim/UCorrelator have +ve theta is up.
 * AnitaEventCorrelator (i.e. all usefulAdu5Pat functions) have +ve theta means down.
 * 
 * @param theta2 an angle to find the difference between
 * @param different_sign_conventions is a boolian to invert theta2, set to true if interfacing between theta derived in anitaEventCorrelator and Acclaim/UCorrelator (default is false)
 * 
 * @return the theta angle between the peak and the source
 */
double AnitaEventSummary::PointingHypothesis::dTheta(double theta2, bool different_sign_conventions) const{
  int factor = different_sign_conventions ? -1 : 1;
  double dTheta = theta - factor*theta2; // + instead of - due to sign convention difference
  return dTheta;
}


/** 
 * Use the summary information to get the angle difference in theta peak.theta - source.theta.
 * 
 * Note:
 * Accounts for the silly sign convention difference, UsefulAdu5Pat has +ve theta is down, the UCorrelator/anitaAnalysisTools have +ve theta is up
 * 
 * @param source is a SourceHypothesis for the event (e.g. WAIS)
 * 
 * @return the theta angle between the peak and the source
 */
double AnitaEventSummary::PointingHypothesis::dThetaSource(const SourceHypothesis& source) const{
  return dTheta(source.theta, true);
}



/** 
 * Get the bearing of the peak (i.e. phi angle from north increasing clockwise)
 * 
 * @return the peak bearing
 */
double AnitaEventSummary::PointingHypothesis::bearing() const{
  // heading increases clockwise, payload phi increases anti-clockwise so we subtract it from heading.
  const AnitaEventSummary* sum = getContainer(__PRETTY_FUNCTION__);
  double bearing = -9999;
  if(sum){
    bearing = double(sum->anitaLocation.heading) - phi;

    bearing = bearing < 0 ? bearing + 360 : bearing;
    bearing = bearing >= 360 ? bearing - 360 : bearing;
  
    if(bearing < 0 || bearing >= 360){
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", peak bearing = "
                << bearing << ", phi = " << phi << ", heading = " << sum->anitaLocation.heading << std::endl;
    }
  }
  return bearing;
}


/** 
 * Get the polarisation of the peak best corresponding to the MC peak
 * In the case of data, just returns the highestPeak
 * 
 * @return the polarisation of the peak closest to MC truth (or highest peak if data)
 */
AnitaPol::AnitaPol_t AnitaEventSummary::mcPol() const {
  findMC();
  return fMCPol;
}

int AnitaEventSummary::mcPolAsInt() const {
  findMC();
  return int(fMCPol);
}

int AnitaEventSummary::mcPeakInd() const {
  findMC();
  return fMCPeakIndex;
}

const AnitaEventSummary::PointingHypothesis& AnitaEventSummary::mcPeak() const {
  findMC();
  return peak[fMCPol][fMCPeakIndex];
}

const AnitaEventSummary::WaveformInfo& AnitaEventSummary::mcCoherent() const {
  findMC();
  return coherent[fMCPol][fMCPeakIndex];
}

const AnitaEventSummary::WaveformInfo& AnitaEventSummary::mcDeconvolved() const {
  findMC();
  return deconvolved[fMCPol][fMCPeakIndex];
}

const AnitaEventSummary::WaveformInfo& AnitaEventSummary::mcCoherentFiltered() const {
  findMC();
  return coherent_filtered[fMCPol][fMCPeakIndex];
}

const AnitaEventSummary::WaveformInfo& AnitaEventSummary::mcDeconvolvedFiltered() const {
  findMC();
  return deconvolved_filtered[fMCPol][fMCPeakIndex];
}


/** 
 * Print warning if fContainer is NULL as the majority the utility functions that rely on it will print nonsense
 * 
 * @param funcName should be the __PRETTY_FUNCTION__ macro for nice debugging info
 * 
 * @return the fContainer pointer
 */
const AnitaEventSummary* AnitaEventSummary::PointingHypothesis::getContainer(const char* funcName) const{
  if(!fContainer){
    std::cerr << "Error in " << funcName
              << " don't have access to AnitaEventSummary that contains me!"
              << " Was the AnitaEventSummary constructor called?"
              << std::endl;
  }
  return fContainer;
}


/** 
 * Return the smaller of the two hardware angles, hwAngle / hwAngleXPol
 * To know which one was smaller, same pol as peak or xpol, see 
 * @return 
 */
double AnitaEventSummary::PointingHypothesis::minAbsHwAngle() const {
  if(TMath::Abs(hwAngle) < TMath::Abs(hwAngleXPol)){
    return hwAngle;
  }
  else{
    return hwAngleXPol;
  }
}

/** 
 * Return the smaller of the two hardware angles, hwAngle / hwAngleXPol
 * To know which one was smaller, same pol as peak or xpol, see 
 * @return 
 */
Bool_t AnitaEventSummary::PointingHypothesis::absHwAngleLessThanAbsHwAngleXPol() const {
  return (TMath::Abs(hwAngle) < TMath::Abs(hwAngleXPol));
}




double AnitaEventSummary::PointingHypothesis::dPhiWais() const {
  return getContainer(__PRETTY_FUNCTION__) ? dPhiSource(fContainer->wais) : dPhi(-9999); // should trigger warning message
}
double AnitaEventSummary::PointingHypothesis::dThetaWais() const {
  return getContainer(__PRETTY_FUNCTION__) ? dThetaSource(fContainer->wais) : dPhi(-9999); // should trigger warning message
}
double AnitaEventSummary::PointingHypothesis::dPhiSun() const {
  return getContainer(__PRETTY_FUNCTION__) ? dPhiSource(fContainer->sun) : dPhi(-9999); // should trigger warning message
}
double AnitaEventSummary::PointingHypothesis::dThetaSun() const {
  return getContainer(__PRETTY_FUNCTION__) ? dThetaSource(fContainer->sun) : dPhi(-9999); // should trigger warning message
}
double AnitaEventSummary::PointingHypothesis::dPhiLDB() const {
  return getContainer(__PRETTY_FUNCTION__) ? dPhiSource(fContainer->ldb) : dPhi(-9999); // should trigger warning message
}
double AnitaEventSummary::PointingHypothesis::dThetaLDB() const {
  return getContainer(__PRETTY_FUNCTION__) ? dThetaSource(fContainer->ldb) : dPhi(-9999); // should trigger warning message
}
double AnitaEventSummary::PointingHypothesis::dPhiMC() const {
  return getContainer(__PRETTY_FUNCTION__) && fContainer->mc.weight > 0 ? dPhiSource(fContainer->mc) : -9999;
}
double AnitaEventSummary::PointingHypothesis::dThetaMC() const {
  return getContainer(__PRETTY_FUNCTION__) && fContainer->mc.weight > 0 ? dThetaSource(fContainer->mc) : -9999;
}


void AnitaEventSummary::findHighestPeak() const {
  resetNonPersistent();
  if(fHighestPeakIndex < 0){ // then we've not done this before
    double highestVal = -1e99;
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(int peakInd=0; peakInd < nPeaks[polInd]; peakInd++){
        if(peak[polInd][peakInd].value > highestVal){
          highestVal = peak[polInd][peakInd].value;
          fHighestPeakIndex = peakInd;
          fHighestPol = pol;
        }
      }
    }
  }
}


/** 
 * Workhorse function to find the peak closest to the MC
 * Caches the result in the mutable, non-ROOT-persistent members fMCPol and fMCPeakIndex
 * In the case of non-MC data, sets the indices to fHighestPol and fHighestPeakIndex
 */
void AnitaEventSummary::findMC() const {
  resetNonPersistent();
  if(mc.weight <= 0){
    findHighestPeak();
    fMCPeakIndex = fHighestPeakIndex;
    fMCPol = fHighestPol;
  }
  else if(fMCPeakIndex < 0){ // then we've not done this before
    double minDeltaAngleSq = 1e99;
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(int peakInd=0; peakInd < nPeaks[polInd]; peakInd++){
        double dPhi = peak[polInd][peakInd].dPhiMC();
        double dTheta = peak[polInd][peakInd].dThetaMC();
        double deltaAngleSq = dPhi*dPhi + dTheta*dTheta;
        if(deltaAngleSq < minDeltaAngleSq){
          minDeltaAngleSq = deltaAngleSq;
          fMCPeakIndex = peakInd;
          fMCPol = pol;
        }
      }
    }
  }
}





/** 
 * Reset the mutable "interesting" indices to defaults
 * The default values trigger the loop through peaks in findHighestPeak() and findMC().
 */
void AnitaEventSummary::resetNonPersistent() const{
  if(fLastEventNumber!=eventNumber){
    fHighestPeakIndex = -1;
    fHighestPol = AnitaPol::kNotAPol;
    fMCPeakIndex = -1;
    fMCPol = AnitaPol::kNotAPol;
    for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      for(Int_t dir=0; dir < maxDirectionsPerPol; dir++){
        peak[polInd][dir].fContainer = const_cast<AnitaEventSummary*>(this); // Set non-persistent pointer to container in hacky fashion.
      }
    }
    fLastEventNumber=eventNumber;
  }  
}


