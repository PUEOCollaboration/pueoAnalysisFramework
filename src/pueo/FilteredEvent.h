#ifndef _PUEO_FILTERED_EVENT_H_
#define _PUEO_FILTERED_EVENT_H_
#include "TObject.h"
#include "pueo/Conventions.h"
#include "pueo/UsefulAttitude.h"
#include "pueo/AnalysisWaveform.h"

class TGraph;
class TCanvas;



namespace pueo 
{


 namespace nav{ class Attitude;}
class UsefulEvent;
class RawHeader;
class FilterStrategy;


/**
  \brief
  This class is intended to store all the necessary data about an ANITA event for filtering and analysis.
  It stores the raw and filtered waveforms as well as auxilliary information.
*/

class FilteredEvent
{

  friend class FilterOperation;
  friend class FilterStrategy;

  public:


   /** Create a FilteredEvent from a useful event, a filter strategy, GPS info and the event header. If keep_all_stages is true, then  all stages of filtering are kept. Note that filtering occurs during construction, not later!*/
   FilteredEvent(const UsefulEvent * event, FilterStrategy * strategy, const nav::Attitude * att, const RawHeader * header, bool keep_all_stages = false);
   /**
      Created a FilteredEvent from another FilteredEvent.
      !THIS IS NOT A COPY CONSTRUCTOR!
      The output of fEv will be used as an input to create this filtered event.
      Useful for doing things like deconvolving the output of another filter.
      If keep_all_stages is true, then  all stages of filtering are kept. Note that filtering occurs during construction, not later!*/
   FilteredEvent(const FilteredEvent* fEv, FilterStrategy * strategy, bool keep_all_stages = false);
  

   /** Empty FilteredEvent. Not particularly useful */
   FilteredEvent();

   /** Destructor */
   virtual ~FilteredEvent();

   /** Accessor for raw waveform based on index. Mostly useful if you want to iterate over everything */
   const AnalysisWaveform * getRawGraph(UInt_t i) const { return rawGraphs[i]; }

   /** Accessor for raw waveform based on antenna number and polarization. */
   const AnalysisWaveform * getRawGraph(UInt_t ant, pol::pol_t pol) const { return rawGraphsByAntPol[pol][ant]; }

   /** Accessor for raw waveform based on phi, ring, and pol */
   const AnalysisWaveform * getRawGraph(UInt_t phi, ring::ring_t ring, pol::pol_t pol) const;

   /** Accessor for filtered waveform based on index. Mostly useful if you want to iterate over everything */
   const AnalysisWaveform * getFilteredGraph(UInt_t i) const { return filteredGraphs[i]; }

   /** Accessor for filtered waveform based on antenna number and polarization. */
   const AnalysisWaveform * getFilteredGraph(UInt_t ant, pol::pol_t pol) const { return filteredGraphsByAntPol[pol][ant]; }

   /** Accessor for raw waveform based on phi, ring, and pol */
   const AnalysisWaveform * getFilteredGraph(UInt_t phi, ring::ring_t ring, pol::pol_t pol) const;


   /* If the FilteredEvent was told to keep_all_stages, you can get the waveform at each stage. Stage 0 is after the first filter operation. */
   const AnalysisWaveform * getFilteredGraphAtStage(UInt_t ant, pol::pol_t pol, UInt_t stage) const;

   /** Return the UsefulEvent this event was constructed with. */
   const UsefulEvent* getUsefulEvent() const { return useful; }
   /** Return the Useful GPS info */
   const UsefulAttitude * getGPS() const { return &att; }

   /** Return the header */
   const RawHeader * getHeader() const { return header; }

   /** Return the strategy */
   const FilterStrategy* getStrategy() const {return strategy;}

   void plotSummary(TCanvas * chpol = 0, TCanvas * cvpol = 0) const;

   int checkSaturation(ULong64_t *save_hsat  =0, ULong64_t* save_vsat = 0, double threshold=1500) const; 
   
	 int checkStepFunction(Int_t lab = 1, ring::ring_t ring = ring::kTopRing, Int_t phiSector = 8, pol::pol_t pol = pol::kVertical) const; 
	 int checkSurfForGlitch(Int_t surf = 0, Int_t lab = -1, double glitchThreshold=800) const; 

   /** Various calculations. Don't necessarily have to be in this class. */
   void getAverageSpectrum (TGraph * target, pol::pol_t pol = pol::kNotAPol ) const;
   void getMedianSpectrum  (TGraph * target, pol::pol_t pol = pol::kNotAPol , double pctile = 0.5) const;
   double getAveragePower(pol::pol_t pol = pol::kNotAPol, ring::ring_t ring = ring::kNotARing, bool filtered = false) const;
   double getMedianPower(pol::pol_t pol = pol::kNotAPol,  ring::ring_t ring = ring::kNotARing, bool filtered = false) const;
   void getMinMaxRatio(pol::pol_t pol, double * max_ratio, double * min_ratio, int* max_sector, int* min_sector, ring::ring_t ring1 = ring::kBottomRing, ring::ring_t ring2 = ring::kTopRing, int nth = 0, int * n_greater = 0) const;

  int getPueoVersion() const { return pueoVersion; }

  
  protected:
   AnalysisWaveform *rawGraphs[k::NUM_ANTS*pol::kNotAPol];
   AnalysisWaveform *rawGraphsByAntPol[pol::kNotAPol][k::NUM_ANTS];
   AnalysisWaveform *filteredGraphs[k::NUM_ANTS*pol::kNotAPol];
   AnalysisWaveform *filteredGraphsByAntPol[pol::kNotAPol][k::NUM_ANTS];


   const UsefulEvent * useful;
   const FilterStrategy * strategy;
   UsefulAttitude att;
   const RawHeader * header;

   int pueoVersion;

   bool keep_all_stages;
   std::vector<AnalysisWaveform *> all_stages[pol::kNotAPol][k::NUM_ANTS] ;
   void saveStage(int nreserve);
};



}


#endif
