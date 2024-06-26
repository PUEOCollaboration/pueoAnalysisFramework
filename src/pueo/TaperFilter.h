#ifndef PUEO_ANALYSIS_TAPER_FILTER_H
#define PUEO_ANALYSIS_TAPER_FILTER_H

#include "pueo/FilterOperation.h" 
#include "TString.h" 

/* These are filters that taper (window?) the waveform. */

namespace pueo
{
class GaussianTaper : public UniformFilterOperation
{
  public: 

    const char * tag () const { return "GaussianTaper"; } 

    const char * description() const { return descStr.Data(); }

    GaussianTaper(double distance_ns = 5, double nsigma = 2, bool filter_beginning = true, bool filter_end = true) 
    : filter_beginning(filter_beginning), filter_end(filter_end), mean(distance_ns), sigma(distance_ns/nsigma)
    {
        descStr = TString::Format("GaussianTaper%s%s(%g,%g)",
                                  filter_beginning ? "B" : "",
                                  filter_end ? "E" : "",
                                  distance_ns,nsigma); 
    }

    virtual void processOne(AnalysisWaveform *, const RawHeader *, int, int) ;

  private: 
    TString descStr; 
    bool filter_beginning; 
    bool filter_end; 
    double mean; 
    double sigma;



}; 
}

#endif

