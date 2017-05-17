#ifndef UCORRELATOR_RESPONSE_MANAGER_H 
#define UCORRELATOR_RESPONSE_MANAGER_H

/* This class does the dirty work of loading appopriate responses
 * for each antenna. See README.response for more details.
 *
 *
 *   Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 * */ 

//#include "AnalysisConfig.h" 
#include <vector>
#include "AnitaConventions.h"

// namespace UCorrelator
// {
namespace AnitaResponse
{

class AbstractResponse; 
class DeconvolutionMethod; 
class ResponseManager
{

 public: 

  // ResponseManager(const UCorrelator::AnalysisConfig * cfg); 
  ResponseManager(const char * responseDir, int npad, const AnitaResponse::DeconvolutionMethod* methodPtr=NULL); 

  const AbstractResponse * response(int pol, int iant) const { return responses[iant][pol]; } 
  const DeconvolutionMethod * getDeconvolutionMethod() const  { return method; }

  virtual ~ResponseManager(); 

      
 private: 
  int loadResponsesFromDir(const char * dir, int npad); 
  const AbstractResponse* responses[NUM_SEAVEYS][2]; 
  std::vector<AbstractResponse*> response_store; 
  const DeconvolutionMethod * method; 


};


}


#endif

