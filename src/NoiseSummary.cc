#include "pueo/NoiseSummary.h"



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
pueo::NoiseSummary::NoiseSummary() {
  for (int poli=0; poli<k::NUM_POLS; poli++) {
    avgMapProf[poli] = NULL;
  }
  zeroInternals();
}

pueo::NoiseSummary::~NoiseSummary() {

  deleteHists();
}

void pueo::NoiseSummary::zeroInternals() {

  fifoLength=0;

  memset(avgRMSNoise,0,k::NUM_ANTS*k::NUM_POLS*sizeof(double));

  for (int poli=0; poli<k::NUM_POLS; poli++) {
    if (avgMapProf[poli] != NULL) {
      delete avgMapProf[poli];
      avgMapProf[poli] = NULL;
    }
  }

  return;
}


void pueo::NoiseSummary::deleteHists() {
  
  for (int poli=0; poli<k::NUM_POLS; poli++) {
    if (avgMapProf[poli] != NULL) {
      delete avgMapProf[poli];
      avgMapProf[poli] = NULL;
    }
  }  

  return;
}

