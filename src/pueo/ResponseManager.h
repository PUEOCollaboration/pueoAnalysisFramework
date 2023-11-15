#ifndef PUEO_RESPONSE_MANAGER_H 
#define PUEO_RESPONSE_MANAGER_H

/* This class does the dirty work of loading appopriate responses
 * for each antenna. See README.response for more details.
 *
 *
 *   Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 * */ 

//#include "AnalysisConfig.h" 
#include <vector>
#include "pueo/Conventions.h"

namespace pueo
{

class AbstractResponse; 
class DeconvolutionMethod; 
class ResponseManager
{

 public: 

  // ResponseManager(const UCorrelator::AnalysisConfig * cfg); 
  ResponseManager(const char * responseDir, int npad, const DeconvolutionMethod* methodPtr=NULL, unsigned int evTime = 0); 

  const AbstractResponse * response(int pol, int iant) const { return responses[iant][pol]; } 
  const DeconvolutionMethod * getDeconvolutionMethod() const  { return method; }
	
  // check time is used to change which response you are using for time dependent responses (such as the TUFF ones)
  void checkTime(unsigned int evTime);

  static std::string getConfigurationFromIndex(const char * indexFile, unsigned int evTime); 

  virtual ~ResponseManager(); 

      
 private: 
  int loadResponsesFromDir(const char * dir, int npad, unsigned int evTime = 0); 

  const AbstractResponse* responses[k::NUM_ANTS][2]; 
  std::vector<AbstractResponse*> response_store; 
  const DeconvolutionMethod * method; 
  int lastTime;
  const char* whichDir;
  int savePad;
  bool hasIndex;


};


}


#endif

