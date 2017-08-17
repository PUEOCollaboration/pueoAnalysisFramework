#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h" 
#include "AnitaConventions.h" 
class Adu5Pat;
class UsefulAdu5Pat;
class RawAnitaHeader; 
class TruthAnitaEvent; 

/** Common analysis output format 
 *  
 *  Needless to say, there's no guarantee that everything will be filled, so be wary if something is 0 (it may not have been filled).   
 *
 * */ 


class AnitaEventSummary 
{
public: 

  /** The maximum number of hypotheses storable per polarization */ 
  static const Int_t maxDirectionsPerPol = 5; 
  static const Int_t peaksPerSpectrum = 3; 

  /** A pointing hypothesis stores the results of interferometric pointing */ 
  class PointingHypothesis 
  {
  public: 
    PointingHypothesis() { ; }
    Double_t phi;  /// peak phi, degrees
    Double_t theta; /// peak theta, degrees
    Double_t value; /// peak value
    Double_t snr; /// snr of peak
    Double_t mapRMS; /// rms of interferometric map
    Double_t mapHistoryVal; /// value of average of the peak location over the past 60 min-bias events
    Double_t hwAngle; /// angle with respect to triggering phi sector

    Double_t latitude;/// on continent, or -9999 if doesn't intersect
    Double_t longitude;/// on continent, or -9999 if doesn't intersect
    Double_t altitude;/// on continent, or -9999 if doesn't intersect
    Double_t distanceToSource; /// on continent, or -9999 if doesn't intersect
    
    Double_t sigma_theta;  ///error on theta
    Double_t sigma_phi;  /// error on phi
    Double_t rho;  ///correlation coefficient between theta and phi
    Double_t chisq; /// chisq/ndof of peak finding process, if available (otherwise zero)


    Double_t theta_adjustment_needed; /// If an event barely missed the ground, it is useful to see the coordinates at which it would hit if theta adjustment by a small amount. This is the calculated small amount that leads to it hitting the ground. 
    Double_t phi_separation; //angular separation from higher value peak in same event. 1000 if highest value event (i.e. first hypothesis) 

    Double_t dphi_rough;  /// phi - phi rough
    Double_t dtheta_rough; /// theta - theta rough 


    Bool_t triggered; /// was this in a triggered phi sector? 
    Bool_t triggered_xpol; /// was this in a triggered xpol phi sector?  
    Bool_t masked; /// was this in a masked phi sector? 
    Bool_t masked_xpol; /// was this in a masked phi xpol sector? 

    ClassDefNV(PointingHypothesis,13); 
  }; 



  /** Stores information about a waveform (coherent or deconvolve) */ 
  class WaveformInfo
  {

  public: 
    WaveformInfo() {; } 
    Double_t snr; ///Signal to Noise of waveform 
    Double_t peakHilbert; /// peak of hilbert envelope
    Double_t peakVal;  /// peak value
    Double_t xPolPeakVal;  // Peak of xpol trace
    Double_t xPolPeakHilbert;  // Peak of xpol hilbert Envelope

    Double_t I,Q,U,V;  // Stokes Parameters

    //some utilities for polarization info
    double linearPolFrac() const;
    double linearPolAngle() const;
    double circPolFrac() const;
    double totalPolFrac() const;


    Double_t totalPower;  ///Total power in waveform
    Double_t totalPowerXpol;  ///Total power in xPol waveform

    //spectrum info 
    Double_t bandwidth[peaksPerSpectrum];  /// bandwidth of each peak (implementation defined, may not be comparable between analyses) 
    Double_t peakFrequency[peaksPerSpectrum]; //peak frequency of power spectrum 
    Double_t peakPower[peaksPerSpectrum]; //power within +/- bandwidth of each peak 
    Double_t spectrumSlope;  ///  Slope of line fit to spectrum (in log-space, so this is spectral-index) 
    Double_t spectrumIntercept; /// Intercept of line fit to spectrum (in log-space) 

    //Shape parameters, computed using hilbert envelope 
    // This should probably taken out into its own class 
    Double_t riseTime_10_90;  /// Rise time of hilbert env from 10% to 90% of peak
    Double_t riseTime_10_50;  /// Rise time of hilbert env from 10% to 50% of peak
    Double_t fallTime_90_10;  /// Fall time of hilbert env from 90% to 10% of peak
    Double_t fallTime_50_10;  /// Fall time of hilbert env from 50% to 10% of peak
    Double_t width_50_50;   /// Width from first envelope crossing of 50 percent of peak to last 
    Double_t width_10_10;  /// Width from first envelope crossing of 10 percent of peak to last 
    Double_t power_10_10;  /// Power enclosed within 10_10 width
    Double_t power_50_50;  /// Power enclosed within 50_50 width
    Double_t peakTime;  // Time that peak hilbert env occurs
    Double_t peakMoments[5];  // moments about Peak  (1st - 5th moments) 


    //See a number that has something to do with how impulsive it is 
    Double_t impulsivityMeasure; 

    Int_t numAntennasInCoherent; // number of antennas used to make this 

    Double_t localMaxToMin; /// Largest value of local max to neighbouring local min (see Acclaim::RootTools::getLocalMaxToMin)
    Double_t localMaxToMinTime; /// Time between local maxima and minima +ve means max is before min, -ve means min is before max
    Double_t globalMaxToMin; /// Difference between maximum and minimum voltage
    Double_t globalMaxToMinTime; /// Time between maximum and minimum volts, +ve means max is before min, -ve means min is before max 

    ClassDefNV(WaveformInfo, 9);
  }; 


  /** Stores various event flags */
  class EventFlags
  {
  public: 
    EventFlags() {; }
    /** Is this event from a cal pulser? */ 
    enum CalPulser 
      {
        NONE, 
        WAIS, 
        LDB, 
        SIPLE,
        TD
      }; 

    Int_t isGood;
    Int_t isRF;
    Int_t isAdu5Trigger;
    Int_t isG12Trigger;
    Int_t isSoftwareTrigger;
    Int_t isMinBiasTrigger;
    Int_t isPayloadBlast;
    Int_t nadirFlag;
    Int_t strongCWFlag;
    Int_t isHPolTrigger;
    Int_t isVPolTrigger;

    CalPulser pulser;
    Bool_t isVarner;
    Bool_t isVarner2;

    /** These are used to cut out payload blasts and stuf like that. 
     *  The first element is the total, and then the next are by ring 
     *  So to get the top ring, do 1 + AnitaRing::kTopRing, etc. 
     */
    Double_t meanPower[1+AnitaRing::kNotARing]; 
    Double_t medianPower[1+AnitaRing::kNotARing]; 
    Double_t meanPowerFiltered[1+AnitaRing::kNotARing]; 
    Double_t medianPowerFiltered[1+AnitaRing::kNotARing]; 

    Double_t maxBottomToTopRatio[AnitaPol::kNotAPol]; 
    int maxBottomToTopRatioSector[AnitaPol::kNotAPol]; 
    Double_t minBottomToTopRatio[AnitaPol::kNotAPol]; 
    int minBottomToTopRatioSector[AnitaPol::kNotAPol]; 

    ClassDefNV(EventFlags,7); 
  };

  /** A Source Hypothesis tells us about different potential sources of signals (e.g. calibration pulser) */ 
  class SourceHypothesis
  {
    public:
      SourceHypothesis() { reset(); }
      Double_t theta;
      Double_t phi;
      Double_t distance;

      Double_t mapValue[NUM_POLS];  ///what the instantaneous map value is at this source hypothesis
      Double_t mapHistoryVal[NUM_POLS]; /// a history of the interferometric map value for the source location

      void reset(); /// sets all the values to nonsense.  Sorry, mapHistoryVal means this is in source now 
      

      ClassDefNV(SourceHypothesis,3);
  };

  class MCTruth : public SourceHypothesis
  {
    public: 
      MCTruth() { reset(); } 
    WaveformInfo wf[AnitaPol::kNotAPol]; 
    double weight; 
    void reset(); 

    ClassDefNV(MCTruth,4);
  }; 

  
  /** Adu5Pat has 16 ints stored with it, but really only 4 are ever important. They are SUPER IMPORTANT, so they should be here */
  class PayloadLocation
  {
  public:
    PayloadLocation() { reset(); }
    PayloadLocation(const Adu5Pat* pat); //!< Slightly more useful constructor

    Float_t latitude;
    Float_t longitude;
    Float_t altitude;
    Float_t heading;

    Float_t prevHeading; //useful for determining rotation rate

    void reset() { latitude = -999; longitude = -999; altitude = -999; heading = -999; prevHeading = -999;};
    void update(const Adu5Pat* pat); //!< Copy the data from the pat into the object

    ClassDefNV(PayloadLocation,2);
  };  
  
  PayloadLocation anitaLocation;

 
  Int_t run;
  UInt_t eventNumber;
  UInt_t realTime;
  
  
  Int_t nPeaks[AnitaPol::kNotAPol]; ///Number of peaks actually found; this might be less than maxDirectionsPerPol 

  PointingHypothesis peak[AnitaPol::kNotAPol][maxDirectionsPerPol]; 

  WaveformInfo coherent[AnitaPol::kNotAPol][maxDirectionsPerPol]; 
  WaveformInfo deconvolved[AnitaPol::kNotAPol][maxDirectionsPerPol];
  WaveformInfo coherent_filtered[AnitaPol::kNotAPol][maxDirectionsPerPol]; 
  WaveformInfo deconvolved_filtered[AnitaPol::kNotAPol][maxDirectionsPerPol];


  // WaveformInfo maxWaveform;  // do we want this for all ? 

  EventFlags flags;

  SourceHypothesis sun;
  SourceHypothesis wais;
  SourceHypothesis ldb;
  MCTruth mc; 

  AnitaEventSummary();
  AnitaEventSummary(const RawAnitaHeader* header);
  AnitaEventSummary(const RawAnitaHeader* header, UsefulAdu5Pat* pat, const TruthAnitaEvent * truth = 0 );  
  void setTriggerInfomation(const RawAnitaHeader* header);
  void setSourceInformation(UsefulAdu5Pat* pat, const TruthAnitaEvent * truth = 0);  
  void zeroInternals();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Utility functions to save precious key strokes when doing a TTree::Draw()
  // Basically, these figure out the higher map peak and return reconstruction variables of that
  // polarisation to you.
  // If you don't want to decide on the primary polarisation of an event based on the higher map peak then
  // these might not be for you...
  //////////////////////////////////////////////////////////////// /////////////////////////////////////////
  const PointingHypothesis& higherPeak(int peakInd=0) const;
  AnitaPol::AnitaPol_t higherPeakPol(int peakInd=0) const; // peak[vpol][peakInd] >= peak[hpol][peakInd] ? vpol : hpol
  const WaveformInfo& higherCoherent(int peakInd=0) const;
  const WaveformInfo& higherDeconvolved(int peakInd=0) const;
  const WaveformInfo& higherCoherentFiltered(int peakInd=0) const;
  const WaveformInfo& higherDeconvolvedFiltered(int peakInd=0) const;

  // Resolution utility functions, for more keystoke saving
  double dPhiSource(const SourceHypothesis& source, int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const; 
  double dThetaSource(const SourceHypothesis& source, int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const;

  double dPhiWais(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const{
    return dPhiSource(wais, peakInd, pol);
  }
  double dThetaWais(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const{
    return dThetaSource(wais, peakInd, pol);
  }
  double dPhiLDB(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const{
    return dPhiSource(ldb, peakInd, pol);
  }
  double dThetaLDB(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const{
    return dThetaSource(ldb, peakInd, pol);
  }
  double dPhiSun(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const{
    return dPhiSource(sun, peakInd, pol);
  }
  double dThetaSun(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const{
    return dThetaSource(sun, peakInd, pol);
  }
  double dPhiMC(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const{
    return dPhiSource(mc, peakInd, pol);
  }
  double dThetaMC(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const{
    return dThetaSource(mc, peakInd, pol);
  }

  double peakBearing(int peakInd=0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol) const;
  double weight(){return mc.weight > 0 ? mc.weight : 1;}
  
  private: 

  ClassDefNV(AnitaEventSummary, 22);
}; 





#endif 
