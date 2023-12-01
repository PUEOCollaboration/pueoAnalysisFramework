#include "pueo/FilterOperation.h" 
#include "TGraph.h" 
#include "pueo/FilteredEvent.h" 
#include <string.h> 
#include <assert.h>


pueo::FilterOperation::~FilterOperation()
{
}

pueo::ConditionalFilterOperation::~ConditionalFilterOperation()
{
 free (condition_tag); 
 free (condition_desc); 
 if (own) delete fo; 
}

pueo::ConditionalFilterOperation::ConditionalFilterOperation(UniformFilterOperation * operation, 
                                                       bool (*condition)(FilteredEvent * ev, int ant, pol::pol_t), 
                                                       const char * condition_tag_suffix, const char * condition_description_suffix, bool should_own_operation) 
                                                       : fn(condition), fo(operation), own(should_own_operation)
{
  int ret = asprintf(&condition_tag, "%s_%s", fo->tag(), condition_tag_suffix) ; 
  assert(ret > 0); 

  ret = asprintf(&condition_desc, "%s (if %s) ", fo->description(), condition_description_suffix); 
  assert( ret > 0); 
}




pueo::AnalysisWaveform* pueo::FilterOperation::getWf(FilteredEvent *ev, int i) 
{ 
  return ev->filteredGraphs[i]; 
}

pueo::AnalysisWaveform* pueo::FilterOperation::getWf(FilteredEvent *ev, int ant, pol::pol_t pol) 
{ 
  return ev->filteredGraphsByAntPol[pol][ant]; 
}


void pueo::ConditionalFilterOperation::process(FilteredEvent * ev) 
{
  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int ant = 0; ant <k::NUM_ANTS; ant++) 
    {
      if (fn(ev,ant, (pol::pol_t)pol))
      {
        fo->processOne(getWf(ev,ant, (pol::pol_t)pol)); 
      }
    }
  }
}

void pueo::ConditionalFilterOperation::processOne(AnalysisWaveform* awf, const RawHeader* header, int ant, int pol) 
{
	printf("processOne not implemented yet, sorry!\n");
}

void pueo::UniformFilterOperation::process(FilteredEvent * ev) 
{
#ifdef USE_OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < k::NUM_ANTS * 2; i++) 
  {
    processOne(getWf(ev,i)); 
  }
}
