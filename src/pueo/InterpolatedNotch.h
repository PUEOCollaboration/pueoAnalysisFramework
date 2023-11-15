#ifndef PUEO_ANALYSIS_INTERPOLATED_NOTCH_H
#define PUEO_ANALYSIS_INTERPOLATED_NOTCH_H

#include "TString.h"
#include "AnalysisWaveform.h"

#include "pueo/FilterOperation.h"

namespace pueo 
{
class InterpolatedNotchFilter : public UniformFilterOperation
{

 public:

  InterpolatedNotchFilter();

  const char * tag () const { return "InterpolatedNotch"; }
  const char * description () const { return "InterpolatedNotchFilter"; }

  virtual void processOne(AnalysisWaveform *, const RawHeader *,int unused , int unsused2) ;
  virtual void process(FilteredEvent* ev);
  virtual void cutParams(bool notchLower,Double_t lowBinCenter,Double_t lowBinSigma,bool notchHigher,
			 Double_t highBinCenter,Double_t highBinSigma);
  TGraph * interpolatedFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);
  
 private:
  Double_t lowBin,highBin,lowSigma,highSigma;
  bool notchLower,notchHigher;

  void init();
};
}

#endif 
