#ifndef _PUEO_FILTER_OPERATION_H
#define _PUEO_FILTER_OPERATION_H

#include "pueo/Conventions.h" 
#include "pueo/RawHeader.h"

namespace pueo 
{
class FilteredEvent; 
class AnalysisWaveform; 

/** A FilteredOperation does things to the waveforms inside a FilteredEvent 
 *
 *
 **/ 

class FilterOperation
{
  public: 
    /** short name for operation, will be used for output tree name, if there is one  */ 
    virtual const char * tag () const = 0; 

    /** human readable description, should provide sufficient information to understand what was done  */ 
    virtual const char * description () const = 0; 
 
    /**  operate on the FilteredEvent */ 
    virtual void process(FilteredEvent * event)= 0; 
    
		/**  operate on one waveform (ABL added moved this from UniformFitlerOperation so that FilterStrategy could call its this on single waveforms, and the weird argument structure is just so it matches ad sinsub) */ 
    virtual void processOne(AnalysisWaveform* awf, const RawHeader * header = 0, int ant = 0, int pol = 0)= 0; 


    // If you want to output stuff to trees, define the output here 
    
    /** The number of output variables (doubles or double arrays) defined by this operation */
    virtual unsigned nOutputs() const  { return 0; } 

    /** The name of the ith output variable */ 
    virtual const char *  outputName(unsigned i) const  { (void) i; return ""; } 

    /** The length of the ith output variable  (it's a double array of this size)*/ 
    virtual unsigned outputLength(unsigned i) const { (void) i; return 0; } 

    /** Fill the ith output */ 
    virtual void fillOutput(unsigned i, double * v) const{ (void) v; (void) i;  return; } 

    /** Destructor */ 
    virtual ~FilterOperation(); 

  protected: 
    /** Accessor for waveform */ 
    AnalysisWaveform * getWf(FilteredEvent *ev, int i); 

    /** Accessor for waveform */ 
    AnalysisWaveform * getWf(FilteredEvent *ev, int ant, pol::pol_t pol); 
}; 


/** For filter operations that do the same thing to each waveform */ 
class UniformFilterOperation : public FilterOperation
{

  public: 
    /** Processes an event, calling processOne on each waveform */ 
    virtual void process(FilteredEvent * event) ; 
    virtual void processOne(AnalysisWaveform* awf, const RawHeader * header = 0, int ant = 0, int pol = 0)=0; 
    virtual ~UniformFilterOperation() {; } 

}; 



/** A ConditionalFilterOperation only applies the passed FilterOperation 
 *  if the condition is true. The tag and description are combinations of the 
 *  passed operation and the condition tag and description provided. 
 */ 
class ConditionalFilterOperation : public FilterOperation
{

  public: 
    /** Convert a UniformFilterOperation to a conditional operation. The UniformFilterOperation will only be applied
     * to the event if the passed condition function returns true */ 

    ConditionalFilterOperation(UniformFilterOperation * operation, 
                               bool (*condition)(FilteredEvent * ev, int ant, pol::pol_t pol), 
                               const char * condition_tag, const char * condition_description, bool own = false) ; 
    
    

    virtual const char * tag() const { return condition_tag; } 
    virtual const char * description () const { return condition_desc; } 

    virtual ~ConditionalFilterOperation(); 
    virtual void process(FilteredEvent * event); 
    virtual void processOne(AnalysisWaveform* awf, const RawHeader * header = 0, int ant = 0, int pol = 0); 

  protected:
    bool (*fn)(FilteredEvent *, int, pol::pol_t);
    char * condition_tag; 
    char * condition_desc; 
    UniformFilterOperation * fo; 
    bool own; 
};


}




#endif 
