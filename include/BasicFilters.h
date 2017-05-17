#ifndef ANITA_ANALYSIS_BASIC_FILTERS
#define ANITA_ANALYSIS_BASIC_FILTERS

#include "FilterOperation.h"
#include "TString.h"
#include "SystemResponse.h"

/** A set of basic filter operations that serve as examples and also should be quite useful */ 


/** SimplePassBandFilter just cuts everything outside the pass band in
 * fourier space like a brick wall, with all the ensuing acausal goodness.
 *
 */ 


namespace FFTtools
{
  class DigitalFilter; 
}

// 
namespace AnitaResponse
{
  class DeconvolutionMethod;
  class ResponseManager;
}

class SimplePassBandFilter : public UniformFilterOperation 
{ 
  public: 
    const char * tag() const { return "SimplePassBandFilter"; } 

    const char * description() const { return descStr.Data(); }

    /** Filter everything outside the pass band. Numbers are given in GHz. */ 
    SimplePassBandFilter(double low = 0.2, double high = 1.2)  
      : low(low), high(high)
    {
        descStr = TString::Format("SimplePassbandFilter(%g,%g)", low,high); 
    }

    virtual void processOne(AnalysisWaveform *) ;

  private: 
    TString descStr; 
    double low; 
    double high;

}; 

/** Brick wall notch filter.  */ 
class SimpleNotchFilter : public UniformFilterOperation
{
  public: 

    /** Ghz*/ 
    SimpleNotchFilter(double minfreq, double maxfreq) 
      : min(minfreq), max(maxfreq) 
    {
      desc.Form("SimpleNotchFilter(%g,%g)",min,max); 
    }

    const char * tag() const { return "SimpleNotchFilter"; } 
    const char * description() const { return desc.Data(); } 
    virtual void processOne(AnalysisWaveform *); 
  private: 
    TString desc;  
    double min,max; 
}; 



/** Brickwall ALFA filter */ 
class ALFASincFilter : public FilterOperation
{
  public: 
    ALFASincFilter(double cutoff = 0.7)
      : pb(0,cutoff) {descStr.Form("ALFA Filter with cutoff at %f GHz",cutoff); } 


    virtual void process(FilteredAnitaEvent * event) 
    {
      if (AnitaVersion::get()!=3) return; 

      pb.processOne( getWf(event, 4, AnitaPol::kHorizontal) ); 
      //cross talk is strong in this one 
      pb.processOne( getWf(event, 12, AnitaPol::kHorizontal) ); 
    }


    const char * tag() const { return "ALFAFilter"; } 
    const char * description() const { return descStr.Data(); } 
  private:
    SimplePassBandFilter pb; 
    double power_before; 
    double power_after;  
    TString descStr;

}; 

//backwards compatibility? 
typedef ALFASincFilter ALFAFilter; 


class HybridFilter : public FilterOperation
{

  public:

    const char * tag () const { return "HybridFilter"; } 
    const char * description () const { return "Hybrid Filter"; } 
    virtual void process(FilteredAnitaEvent * event); 
};

class SumDifferenceFilter : public FilterOperation
{

  public:

    const char * tag () const { return "SumDifferenceFilter"; } 
    const char * description () const { return "SumDifference Filter"; } 
    virtual void process(FilteredAnitaEvent * event); 
};


class DigitalFilterOperation : public UniformFilterOperation 
{
  public: 
    /**Digital filter based on digi. 
     * If correct_delay is true, the waveform will be adjusted by average group delay over the band. 
     * The band is defined by delay_min_freq and delay_max_freq (given in terms of fnyq). If supersampled, you may need to adjust
     * these from the values here. 
     *
     */
    DigitalFilterOperation(const FFTtools::DigitalFilter *digi, bool correct_delay = true, 
                          double delay_min_freq = 0.18/1.3, double delay_max_freq = 1); 
    const char * tag () const { return "DigitalFilter"; } 
    const char * description () const { return "DigitalFilter"; } 
    virtual void processOne(AnalysisWaveform* wf); 

  private: 
    const FFTtools::DigitalFilter * digi; 
    double delay; 


}; 

class ALFALanczosFilter : public FilterOperation
{
  public:
    ALFALanczosFilter(double cutoff = 0.6, int a= 3); 
    virtual ~ALFALanczosFilter(); 
    const char * tag() const { return "ALFAFilter"; } 
    const char * description() const { return descStr.Data(); } 
    virtual void process(FilteredAnitaEvent * event) ; 

  private: 
    DigitalFilterOperation *pb; 
    FFTtools::DigitalFilter * filt;
    TString descStr; 

}; 

class ALFAButterworthFilter : public FilterOperation
{
  public: 
    ALFAButterworthFilter(double cutoff = 0.55); /** The cutoff is scaled by 1.3, so if you oversampled, this won't be right anymore */ 
    virtual ~ALFAButterworthFilter(); 

    virtual void process(FilteredAnitaEvent * event) ; 

    const char * tag() const { return "ALFAFilter"; } 
    const char * description() const { return descStr.Data(); } 
  private:
    DigitalFilterOperation *pb; 
    FFTtools::DigitalFilter * filt;
    double delay; 
    double power_before; 
    double power_after;  
    TString descStr;

};

namespace AnitaResponse{
class DeconvolveFilter : public FilterOperation
{
 public: 

  DeconvolveFilter(const AnitaResponse::ResponseManager *rm, const AnitaResponse::DeconvolutionMethod * dm) 
      : rm(rm), dm(dm)  
  {;} 
      
  

  virtual const char * tag() const { return "DeconvolveFilter"; } 
  virtual const char * description() const { return "DeconvolveFilter"; } 

  virtual void process(FilteredAnitaEvent * ev); 
      

 private: 
  const AnitaResponse::ResponseManager *rm; 
  const AnitaResponse::DeconvolutionMethod *dm; 

}; 
}


#endif 
