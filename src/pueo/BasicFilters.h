#ifndef PUEO_ANALYSIS_BASIC_FILTERS
#define PUEO_ANALYSIS_BASIC_FILTERS

#include "pueo/FilterOperation.h"
#include "TString.h"
#include "pueo/SystemResponse.h"
#include "pueo/RawHeader.h"

/** A set of basic filter operations that serve as examples and also should be quite useful */ 


/** SimplePassBandFilter just cuts everything outside the pass band in
 * fourier space like a brick wall, with all the ensuing acausal goodness.
 *
 */ 


namespace FFTtools
{
  class DigitalFilter; 
}

 


namespace pueo 
{
  class DeconvolutionMethod;
  class ResponseManager;

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

    virtual void processOne(AnalysisWaveform *, const RawHeader * = 0, int = 0, int pol = 0) ;

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
    virtual void processOne(AnalysisWaveform *, const RawHeader * = 0, int = 0, int pol = 0); 
  private: 
    TString desc;  
    double min,max; 
}; 




class HybridFilter : public FilterOperation {

  public:

    const char * tag () const { return "HybridFilter"; } 
    const char * description () const { return "Hybrid Filter"; } 
    virtual void process(FilteredEvent * event); 
    virtual void processOne(AnalysisWaveform * awf, const RawHeader * header = 0, int ant=0, int pol = 0);
};


class SumDifferenceFilter : public FilterOperation {

  public:

    const char * tag () const { return "SumDifferenceFilter"; } 
    const char * description () const { return "SumDifference Filter"; } 
    virtual void process(FilteredEvent * event);
    virtual void processOne(AnalysisWaveform * awf, const RawHeader * header = 0, int ant = 0, int pol = 0);
};


class FlipHVFilter : public FilterOperation {

  public:

    const char * tag () const { return "FlipHVFilter"; } 
    const char * description () const { return "FlipHVFilter"; } 
    virtual void process(FilteredEvent * event);
    virtual void processOne(AnalysisWaveform * awf, const RawHeader * header = 0, int ant = 0, int pol = 0);
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
    virtual void processOne(AnalysisWaveform* wf, const RawHeader * = 0, int ant = 0, int pol = 0); 

  private: 
    const FFTtools::DigitalFilter * digi; 
    double delay; 


}; 

; 


/** Removes glitches */

class DeglitchFilter : public UniformFilterOperation
{
  public:

    enum RemoveAction { DELETE, AVERAGE} ; 
    /** 
     * Construct the filter. 
     * @param threshold_above_neighbors the minimum abs(v) threshold above neigbhors to be considered a glitch
     * @param n_neighbors number of neighbors (on both sides, so up to 2 * n_neighbors + 1 are considered 
     * @param action What to do when removing a point. If DELETE is used, will function on uneven waveform. AVERAGE will function on even waveform. 
     * @param max_remove The maximum number of points to remove. If more than this are attempted, None are removed (since probably something silly is going on)
     *
     * */
    DeglitchFilter(double threshold_above_neighbors = 100, int n_neighbors = 9, RemoveAction action = DELETE, int max_remove = 10); 
    virtual ~DeglitchFilter() {;}
    virtual void processOne(AnalysisWaveform *wf, const RawHeader * = 0, int ant = 0, int pol = 0); 
    const char * tag() const { return "DeglitchFilter"; } 
    const char * description() const { return descStr.Data(); } 
       
  private: 
  TString descStr; 
  RemoveAction action; 
  double thresh; 
  int neighbors; 
  int nremoved; 
  int max_remove;



}; 


class DeconvolveFilter : public FilterOperation
{
 public: 

  DeconvolveFilter(const ResponseManager *rm, const DeconvolutionMethod * dm) 
      : rm(rm), dm(dm)  
  {;} 
      
  

  virtual const char * tag() const { return "DeconvolveFilter"; } 
  virtual const char * description() const { return "DeconvolveFilter"; } 

  virtual void process(FilteredEvent * ev); 
	virtual void processOne(AnalysisWaveform * awf, const RawHeader * header=0, int ant=0, int pol = 0);
      

 private: 
  const ResponseManager *rm; 
  const DeconvolutionMethod *dm; 

}; 


/**
 *  Filter corresponding to differintegral operator in Fourier domain.
 *  Modifying scipy.fftpack.diff from SciPy
 *  <https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.fftpack.diff.html>
 *  for usage with FFTtools. Assumes the waveform is uniformly sampled, zero-meaned when
 *  order is negative, isn't constant after differentiation, and the result is expected to be
 *  real over the input domain, so don't use it otherwise.
 **/ 
class IFFTDiffFilter : public UniformFilterOperation
{
  public:

    IFFTDiffFilter(double order = 1, int branchOrder = 0) 
      : order(order), branchOrder(branchOrder) 
    {
      desc.Form("IFFTDiffFilter(%g,%d)", order, branchOrder); 
    }
    /** "order" refers to order of operation, with positive values corresponding to orders of differation,
     *  negative to integration. Default is simple differentiation.
     *  "branchOrder" provides a unique result when "order" is noninteger, then "branchOrder" corresponds to
     *  a branch cut in complex analysis.
     **/

    const char * tag() const { return "IFFTDiffFilter"; } 
    const char * description() const { return desc.Data(); } 
    virtual void processOne(AnalysisWaveform *, const RawHeader * = 0, int = 0, int pol = 0);

  private: 
    TString desc;  
    double order;
    int branchOrder; 
}; 


}

#endif 
