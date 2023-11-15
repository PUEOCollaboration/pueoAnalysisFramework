#ifndef PUEO_NOISE_SUMMARY
#define PUEO_NOISE_SUMMARY

#include "math.h"
#include <sstream>

#include "TObject.h"
#include "TH2D.h"
#include "TProfile2D.h"

#include "pueo/Conventions.h"
#include "pueo/FilteredEvent.h"

namespace pueo 
{
/*===================
  A class to store information about the thermal environment */
class NoiseSummary
{
 public:
  NoiseSummary();

  virtual ~NoiseSummary();

  //functions
  void zeroInternals();
  void deleteHists();

  //size
  int fifoLength; //comes from AnitaNoiseMachine
  static const int nPhi = 180; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
  static const int nTheta = 100; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic


  //flags
  bool isMinBias; //is it a min bias event?  Can use this to determine where the updates are
  bool mapFifoFillFlag; //has the fifo filled up yet?
  bool rmsFifoFillFlag; //has the fifo filled up yet?

  //data
  double avgRMSNoise[k::NUM_PHI][ring::kNotARing][k::NUM_POLS]; //averaged rms of min bias waveforms
  TH2D *avgMapProf[k::NUM_POLS]; //Averaged histograms over fifoLength minBias samples


 private:

  ClassDefNV(NoiseSummary,1);
  
};
}

/*--------------*/

#endif
