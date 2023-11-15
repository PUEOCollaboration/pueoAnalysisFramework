#ifndef PUEO_ANALYSIS_DIODE_FILTER_H
#define PUEO_ANALYSIS_DIODE_FILTER_H

#include "pueo/FilterOperation.h"
#include "TString.h" 
#include "pueo/AnalysisWaveform.h" 

/* Implements the time domain diode model */ 


namespace pueo 
{

class DiodeFilter : public UniformFilterOperation
{

  public: 

    DiodeFilter(); 

    const char * tag () const { return "DiodeFilter"; } 
    const char * description () const { return "DiodeFilter"; } 

    virtual void processOne(AnalysisWaveform *, const RawHeader *, int, int) ;

  private: 
    AnalysisWaveform response; 





}; 
}



#endif


