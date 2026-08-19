#include "pueo/Nav.h"
#include "pueo/Conventions.h"
#include "TMath.h"
#include "Rtypes.h"

/** 
 * @class PueoEventSummary
 * @brief Reduced Data Format for PUEO analysis. Heavily lifted from AnitaEventSummary
 *
 * Lets see how well this does
 */

namespace pueo
{
class EventSummary
{
public:


  //------------------------------------------------------------------------------------
  /*************************************************************************************
     * Static members (for array sizes inside class)
     *************************************************************************************/
  static const Int_t maxDirectionsPerPol = 3; /// The maximum number of hypotheses storable per polarization
  static const Int_t peaksPerSpectrum = 3; /// The maximum number of frequency peaks per waveform spectrum
  static const Int_t numFracPowerWindows = 5; /// The number of considered windows to integrate around for power estimates

  //------------------------------------------------------------------------------------
  /*************************************************************************************
     * Sub-classes
     *************************************************************************************/


  /**
     * @class PointingHypothesis
     * @brief Stores some summary information about the peak of an interferometric map.
     */

  class PointingHypothesis 
  {

  public:
    PointingHypothesis() {; }
    Double_t az;  /// peak az, degrees
    Double_t el; /// peak el, degrees
    Double_t value; /// peak value
    Double_t snr; /// snr of peak
    Double_t mapRMS; /// rms of interferometric map
    Double_t mapHistoryVal; /// value of average of the peak location over the past 60 min-bias events
    Double_t latitude;/// on continent, or -999 if doesn't intersect
    Double_t longitude;/// on continent, or -999 if doesn't intersect
    Double_t altitude;/// on continent, or -999 if doesn't intersect
    Double_t distanceToSource; /// on continent, or -999 if doesn't intersect
    Double_t el_adjustment_needed; /// If an event barely missed the ground, it is useful to see the coordinates at which it would hit if el adjustment by a small amount. This is the calculated small amount that leads to it hitting the ground.
    Double_t az_separation; /// angular separation from higher value peak in same event. -999 if highest value event (i.e. first hypothesis)
    Double_t daz_rough;  /// az - az rough
    Double_t del_rough; /// el - el rough

    Double_t width_el;  /// width in el
    Double_t width_az;  /// width in az
    Double_t chisq; /// chisq/ndof of peak finding process, if available (otherwise zero)

    Double_t hwAngle; /// angle in degrees with respect to closest triggering phi sector in same polarization, or -999 if no trigger in same polarization
    Double_t hwAngleXPol; /// angle in degrees with respect to closest triggering phi sector in opposite polarization, or -999 if not trigger in opposite polarization
    Bool_t triggered; /// was this in a triggered phi sector?
    Bool_t triggered_xpol; /// was this in a triggered xpol phi sector?
    Bool_t masked; /// was this in a masked phi sector?
    Bool_t masked_xpol; /// was this in a masked phi xpol sector?

  private:

    ClassDefNV(PointingHypothesis,1);
  }; //end PointingHypothesis

  /**
     * @class WaveformInfo
     * @brief Stores information about a coherently summed waveform (filtered/unfiltered/deconvolved)
     * The coherent summing of the waveform corresponds to a direction stored in a PointingHypothesis
     * Units are assumed to be in mV, ns,and GHz. 
     */
  class WaveformInfo
  {

  public:
    WaveformInfo() {; }

    Double_t snr; // Signal to Noise of waveform

    Double_t pk2pk; // peak 2 peak value
    Double_t rms; // rms of coherentsum
    Double_t peakHilbert;// peak of hilbert envelope
    Double_t rmsHilbert; // rms of hilbert envelope

    Double_t xPolPk2pk; // Peak of xpol trace
    Double_t xPolRms; // RMS of xpol trace
    Double_t xPolPeakHilbert; // Peak of xpol hilbert Envelope
    Double_t xPolRmsHilbert; // RMS of xpol hilbert Envelope

    Double_t totalPower;  /// Total power in waveform
    Double_t totalPowerXpol;  /// Total power in xPol waveform

    Double_t I,Q,U,V;  /// Integral Stokes Parameters (over the entire waveform) 
    Double_t max_dI,max_dQ,max_dU,max_dV,polErr; /// instantanteous stokes parameters (computed near max_dI).
    Int_t NPointsMaxStokes; /// The number of points used in the above estimates 

    //spectrum info
    Double_t bandwidth[peaksPerSpectrum]; /// bandwidth of each peak 
    Double_t peakFrequency[peaksPerSpectrum];  /// peak frequency of power spectrum
    Double_t peakPower[peaksPerSpectrum]; //power within +/- bandwidth of each peak
    Double_t spectrumSlope; ///  Slope of line fit to spectrum (in log-space, so this is spectral-index)
    Double_t spectrumIntercept;  /// Intercept of line fit to spectrum (in log-space)

    //Shape parameters, computed using hilbert envelope
    Double_t scRiseTime_10_65; // Fraction of window length it takes to integrate from 10% total power to 65% total power
    Double_t scRiseTime_20_65; // Fraction of window length it takes to integrate from 20% total power to 65% total power
    Double_t scRiseTime_30_65; // Fraction of window length it takes to integrate from 30% total power to 65% total power

    Double_t numPeaks; // Number of peaks in the waveform

    Double_t riseTime_20_80; // Rise time from 20% to 80% of peak
    Double_t riseTime_50_80; // Rise time from 50% to 80% of peak
    Double_t riseTime_20_50; // Rise time from 20% to 50% of peak

    Double_t fallTime_80_20; // Fall time from 80% to 20% of peak
    Double_t fallTime_80_50; // Fall time from 50% to 80% of peak
    Double_t fallTime_50_20; // Fall time from 50% to 20% of peak

    Double_t width_50_50; // Width from first envelope crossing of 50 percent of peak to last
    Double_t width_20_20;  // Width from first envelope crossing of 20 percent of peak to last

    Double_t power_20_20;  /// Power enclosed within 20_20 width
    Double_t power_50_50;  /// Power enclosed within 50_50 width

    Double_t peakTime;  // Time that max peak occurs

    Double_t fracPowerWindowBegins[numFracPowerWindows]; // Narrowest width containing {10%, 20%, 30%, 40%, 50%} of the total power. Will be 0 if percentage is never reached (possible for low peak HE)
    Double_t fracPowerWindowEnds[numFracPowerWindows]; // Narrowest width containing {10%, 20%, 30%, 40%, 50%} of the total power. Will be length of envelope if percentage is never reached (possible for low peak HE)

    Int_t numAntennasInCoherent; /// Number of antennas used to make this

    //some utilities for polarization info
    Double_t linearPolFrac() const { return TMath::Sqrt( pow(Q,2) + pow(U,2) ) / I; };
    Double_t linearPolAngle() const { return (TMath::ATan(U/Q)/2)*TMath::RadToDeg(); };
    Double_t circPolFrac() const { return TMath::Abs(V)/I; };
    Double_t totalPolFrac() const { return TMath::Sqrt(pow(Q,2) + pow(U,2) + pow(V,2))/I; };

    //utilities for instantaneous polarization info
    Double_t instantaneousLinearPolFrac() const { return TMath::Sqrt( pow(max_dQ,2) + pow(max_dU,2) ) / max_dI; };
    Double_t instantaneousLinearPolAngle() const { return (TMath::ATan(max_dU/max_dQ)/2)*TMath::RadToDeg(); };
    Double_t instantaneousCircPolFrac() const { return TMath::Abs(max_dV)/max_dI; };
    Double_t instantaneousTotalPolFrac() const { return TMath::Sqrt(pow(max_dQ,2) + pow(max_dU,2) + pow(max_dV,2))/max_dI; };

  private:

    ClassDefNV(WaveformInfo, 1);

  }; //end WaveformInfo

  /**
     * @class ChannelInfo
     * @brief Stores brief information of a channel's waveform
     */
  class ChannelInfo
  {

  public:

    ChannelInfo() : pol(pueo::pol::kNotAPol), ant(-1) {; }

    Double_t rms; 
    Double_t avgPower;
    Double_t snr;
    Double_t peakHilbert;
    std::vector<Double_t> sinsubtract_freq; //in GHz
    std::vector<Double_t> sinsubtract_phase; //from -pi to +pi
    std::vector<Double_t> sinsubtract_amp; //in nominal mV

    inline Int_t getAnt() const {return ant;} 
    inline pueo::pol::pol_t getPol() const {return pol;} 

  private:

    pueo::pol::pol_t pol; //! DOES NOT PERSIST IN ROOT! the polarization
    Int_t ant; //! DOES NOT PERSIST IN ROOT! the antenna

    ClassDefNV(ChannelInfo, 1);
  }; //end ChannelInfo

  /**
     * @class EventFlags
     * @brief Stores simple integer numbers based on quality cuts, calibration pulser timing, trigger type, etc.
     */
  class EventFlags
  {
  public:
    EventFlags() {; }
    /** Is this event from a cal pulser? */
    enum CalPulser
    {
      NONE,
      LDB,
      TD,
      S200,
      HICALA,
      HICALB
    };

    Int_t isRF;
    Int_t isSoftwareTrigger;
    Int_t isPPSTrigger;
    Int_t isHPolTrigger;
    Int_t isVPolTrigger;

    //Quality Cut Flags
    Int_t isPayloadBlast;
    Double_t blastImpulsivity; //impulsivity measure for blast checking
    Double_t blastCW; //CW measure for blast checking
    Int_t isSaturated;
    Int_t isShortRun; //for short runs
    Int_t isMissingFT; //for runs with no force triggers
    Int_t isMissingPPS; //for runs with no PPS
    Int_t unsualFT; //force trigger / pps rates too high for this run
    Int_t noRF; //no RF triggers in run
    Int_t noDaqHsk; //no DAQ HSK information for htis run
    Int_t lowThermal; //Median of all channels have low RMS values
    Int_t lowChannelThermal;  //at least one channel has low RMS values
    Int_t strongCWFlag;
    Int_t hasGlitch; //Data is flatlined in some way
    Int_t isInDeadtime; //Event occurs in first 2 seconds of a run, before AGC and masking have settled
    Int_t allTrigsMasked; //All triggered channels are in masked off phi sectors
    Int_t isMUOSContaminated; //If MUOS is present in trigger

    CalPulser pulser;

    /** The fraction of nearby events that are payload blasts */
    Double_t blastFraction;

    ClassDefNV(EventFlags,2);
  }; //end EventFlags

  /**
     * @class SourceHypothesis
     * For known sources, such as a calibration pulser
     */
  class SourceHypothesis
  {
  public:
    SourceHypothesis() { reset(); }
    Double_t el;
    Double_t az;
    Double_t distance;

    Double_t mapValue[pueo::k::NUM_POLS]; ///what the instantaneous map value is at this source hypothesis

    void reset() {el = -999;
      az = -999;
      distance = -999;

      memset(mapValue,0,pueo::k::NUM_POLS*sizeof(Double_t));
    };

    ClassDef(SourceHypothesis,1);
  }; //end SourceHypothesis

  /**
     * @class PayloadLocation
     * Contains the most important 4 numbers in the GPS data.
     */
  class PayloadLocation
  {
  public:
    PayloadLocation() { reset(); }
    PayloadLocation(const pueo::nav::Attitude* att); /// Slightly more useful constructor

    Float_t latitude;
    Float_t longitude;
    Float_t altitude;
    Float_t heading;

    Float_t prevHeading; //useful for determining rotation rate

    void reset() { latitude = -999; longitude = -999; altitude = -999; heading = -999; prevHeading = -999;};
    void update(const pueo::nav::Attitude* att) {
      if(att){
        latitude = att->latitude;
        longitude = att->longitude;
        altitude = att->altitude;

        prevHeading = heading;
        heading = att->heading;

      }
      else{
        reset();
      }
    }; /// Copy the data from the pat into the object

    /**
       * Convert to an Attitude object,
       * @return an Attitude object only with location information
       */
    pueo::nav::Attitude att () const {
      pueo::nav::Attitude att;
      att.longitude = longitude;
      att.latitude = latitude;
      att.altitude = altitude;
      att.heading = heading;
      return att;
    }

    ClassDefNV(PayloadLocation,1);
  };//end PayloadLocation

  //------------------------------------------------------------------------------------
  /*************************************************************************************
     * Public member variables
     *************************************************************************************/
  Int_t run; /// Run
  UInt_t eventNumber; /// Event number
  Double_t realTime; /// Time of the event
  std::array<int, pueo::pol::kNotAPol> nPeaks; /// Number of peaks stored in this EventSummary (might be less than maxDirectionsPerPol)

  PointingHypothesis peak[pueo::pol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the event peak directions (indices of all WaveformInfo member arrays match peak index)
  WaveformInfo coherent[pueo::pol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the (unfiltered) coherently summed waveforms, array index correponds to entry in peak[][]
  WaveformInfo dedispersed[pueo::pol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the (unfiltered) de-dispersed coherently summed waveforms, array index correponds to entry in peak[][]
  WaveformInfo coherent_filtered[pueo::pol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the filtered, coherently summed waveforms, array index correponds to entry in peak[][]
  WaveformInfo dedispersed_filtered[pueo::pol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the filtered, de-dispersed, coherently summed waveforms, array index correponds to entry in peak[][]

  ChannelInfo channels[pueo::pol::kNotAPol][pueo::k::NUM_ANTS]; /// Summaries of each channel's waveform.

  EventFlags flags; /// Flags corresponding the event quality, trigger type, calibration pulser timing, etc.

  SourceHypothesis sun; /// Contains location of sun in map coordinates at time of event
  SourceHypothesis muos1; /// Contains location of MUOS 1 in map coordinates at time of event
  SourceHypothesis muos2; /// Contains location of MUOS 2 in map coordinates at time of event
  SourceHypothesis muos3; /// Contains location of MUOS 3 in map coordinates at time of event
  SourceHypothesis muos4; /// Contains location of MUOS 4 in map coordinates at time of event
  SourceHypothesis muos5; /// Contains location of MUOS 5 in map coordinates at time of event

  SourceHypothesis td; /// Contains location of TD cal pulser in map coordinates at time of event
  SourceHypothesis ldb; /// Contains location of LDB cal pulser in map coordinates at time of event
  SourceHypothesis s200; /// Contains location of S+200 cal pulser in map coordinates at time of event

  SourceHypothesis hicalA; /// Contains location of HiCalA in map coordinates at time of event
  SourceHypothesis hicalB; /// Contains location of HiCalB in map coordinates at time of event

  //MCTruth mc; /// Contains summary information about MC truth, if real data then this filled with constant, unphysical values.
  PayloadLocation location; /// Reduced GPS data

private:

  ClassDefNV(EventSummary, 1);
}; //end class definition
} //end name space
