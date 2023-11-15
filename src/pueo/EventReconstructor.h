#ifndef _PUEO_EVENT_RECONSTRUCTOR_H
#define _PUEO_EVENT_RECONSTRUCTOR_H

namespace  pueo 
{

class FilteredEvent; 
class EventSummary;
class UsefulAttitude;

class EventReconstructor 
{
  public:
    virtual void process(const FilteredEvent * ev, UsefulAttitude* usefulPat, EventSummary * summary) const = 0; 
    virtual ~EventReconstructor() { ; }
}; 
}

#endif 










