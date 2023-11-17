#include "pueo/FilteredEvent.h"
#include "pueo/UsefulEvent.h"
#include "pueo/FilterStrategy.h"
#include "pueo/RawHeader.h"
#include "pueo/AnalysisWaveform.h"
#include "pueo/Version.h"
#include "pueo/GeomTool.h"
#include <algorithm>
#include "TCanvas.h"
#include "TStyle.h"

// maybe this should be configurable? :) 
static const int pad_to_size = pueo::k::MAX_NUMBER_SAMPLES; 

//ClassImp(pueo::FilteredEvent);

pueo::FilteredEvent::FilteredEvent()
{
  for (unsigned i = 0; i < 2 * k::NUM_ANTS; i++ )
    {
    filteredGraphs[i] = NULL;
    rawGraphs[i] = NULL;
  }
}


pueo::FilteredEvent::FilteredEvent(const pueo::FilteredEvent* fEv, FilterStrategy * strategy, bool save_stages )
  : useful(fEv->getUsefulEvent()),
    strategy(strategy),
    att(fEv->getGPS()),
    header(fEv->getHeader()),
    keep_all_stages(save_stages)
{

#ifdef MULTIVERSION_PUEO_ENABLED
  pueoVersion = version::getVersionFromUnixTime(header->realTime);
#else
  pueoVersion = 0;
#endif

  if (pueoVersion <= 0) pueoVersion = version::get();

  // Initialize the filtered graphs with the raw graphs from Raw Anita Event


  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {
    for (unsigned ant = 0; ant < k::NUM_ANTS; ant++)
    {
      int k = pol * k::NUM_ANTS  + ant;
      // int i = geom->getChanIndexFromAntPol(ant, (pol::pol_t) pol);
      // filteredGraphs[k] = new AnalysisWaveform (useful->fNumPoints[i], useful->fTimes[i], useful->fVolts[i]);
      // rawGraphs[k] = new AnalysisWaveform (useful->fNumPoints[i], useful->fTimes[i], useful->fVolts[i]);

      const AnalysisWaveform* rawIn = fEv->getRawGraph(k);
      const AnalysisWaveform* filteredIn = fEv->getFilteredGraph(k);
      
      rawGraphs[k] = new AnalysisWaveform (*rawIn);      
      filteredGraphs[k] = new AnalysisWaveform(*filteredIn);
      
      filteredGraphs[k]->forceEvenSize(pad_to_size); // do this for correlations
//      rawGraphs[k]->forceEvenSize(max_size); // do this for correlations
      filteredGraphsByAntPol[pol][ant] = filteredGraphs[k];
      rawGraphsByAntPol[pol][ant] = rawGraphs[k];
    }
  }


  //tell the strategy to process this
  strategy->process(this);

  
  int max_size = pad_to_size;
  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {
    for (unsigned ant = 0; ant < k::NUM_ANTS; ant++)
    {
      int k = pol * k::NUM_ANTS  + ant;
      if (filteredGraphs[k]->Neven() > pad_to_size) 
      {
        max_size = filteredGraphs[k]->Neven(); 
      }
    }
  }

  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {
    for (unsigned ant = 0; ant < k::NUM_ANTS; ant++)
    {
      int k = pol * k::NUM_ANTS  + ant;
      filteredGraphs[k]->forceEvenSize(max_size); // do this for correlations
    }
  }
}


pueo::FilteredEvent::FilteredEvent(const UsefulEvent * event, FilterStrategy * strategy, const nav::Attitude * rawatt, const RawHeader * header, bool save_stages )
  : useful(event),
    strategy(strategy),
    att(rawatt),
    header(header),
    keep_all_stages(save_stages)
{

#ifdef MULTIVERSION_ANITA_ENABLED
  pueoVersion = version::getVersionFromUnixTime(header->realTime);
#else
  pueoVersion = 0;
#endif

  if (pueoVersion <= 0) pueoVersion = version::get();

  auto geom = &pueo::GeomTool::Instance();
  // Initialize the filtered graphs with the raw graphs from Raw Anita Event

  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {
    for (unsigned ant = 0; ant < k::NUM_ANTS; ant++)
    {
      int k = pol * k::NUM_ANTS  + ant;
      int i = geom->getChanIndexFromAntPol(ant, (pol::pol_t) pol);
      filteredGraphs[k] = new AnalysisWaveform (useful->volts[i].size(), &useful->volts[i][0], useful->dt[i], useful->t0[i]);
      rawGraphs[k] = new AnalysisWaveform (useful->volts[i].size(), &useful->volts[i][0], useful->dt[i], useful->t0[i]); 
      filteredGraphs[k]->forceEvenSize(pad_to_size); // do this for correlations
//      rawGraphs[k]->forceEvenSize(260); // do this for correlations
      filteredGraphsByAntPol[pol][ant] = filteredGraphs[k];
      rawGraphsByAntPol[pol][ant] = rawGraphs[k];
    }
  }


  //tell the strategy to process this
  strategy->process(this);

  
  int max_size = pad_to_size; 

  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {
    for (unsigned ant = 0; ant < k::NUM_ANTS; ant++)
    {
      int k = pol * k::NUM_ANTS  + ant;
      if (filteredGraphs[k]->Neven() > max_size) 
      {
        max_size = filteredGraphs[k]->Neven(); 
      }
    }
  }

  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {
    for (unsigned ant = 0; ant < k::NUM_ANTS; ant++)
    {
      int k = pol * k::NUM_ANTS  + ant;
      filteredGraphs[k]->forceEvenSize(max_size); // do this for correlations
    }
  }
}



const pueo::AnalysisWaveform * pueo::FilteredEvent::getFilteredGraphAtStage(UInt_t ant, pol::pol_t pol, UInt_t stage) const
{
  if (!keep_all_stages)
  {
    fprintf(stderr,"You didn't ask to save all the stages!\n");
    return 0;
  }

  if (stage >= all_stages[pol][ant].size()) return filteredGraphsByAntPol[pol][ant];
  else return all_stages[pol][ant][stage];
}


void pueo::FilteredEvent::saveStage(int nreserve)
{

  for (int pol = 0; pol < 2; pol++)
  {
    for (int ant = 0; ant < k::NUM_ANTS; ant++)
    {
      all_stages[pol][ant].reserve(nreserve);
      all_stages[pol][ant].push_back(new AnalysisWaveform(*filteredGraphsByAntPol[pol][ant]));
    }
  }
}

pueo::FilteredEvent::~FilteredEvent()
{
  for (unsigned pol = 0; pol < 2; pol++)
  {
    for (unsigned ant = 0; ant <  k::NUM_ANTS; ant++ )
    {
      delete filteredGraphsByAntPol[pol][ant];
      delete rawGraphsByAntPol[pol][ant];

      for (unsigned j = 0; j < all_stages[pol][ant].size() ; j++)
      {
        delete all_stages[pol][ant][j];
      }
    }
  }
}


double pueo::FilteredEvent::getAveragePower(pol::pol_t pol, ring::ring_t ring, bool filtered) const
{

  pol::pol_t pol_start =pol == pol::kVertical ? pol::kVertical : pol::kHorizontal;
  pol::pol_t pol_end = pol == pol::kHorizontal ? pol::kHorizontal : pol::kVertical;

  double sum = 0;

  int n = k::NUM_ANTS * (int(pol_end) - int(pol_start) +1);
  if (ring != ring::kNotARing) n/=ring::kNotARing;

  for (int pol = (int) pol_start; pol <=(int) pol_end; pol++ )
  {
    for (int ant = 0; ant < k::NUM_ANTS; ant++)
    {
      if (ring != ring::kNotARing && pueo::GeomTool::getRingFromAnt(ant) != ring) continue;
      sum += ( filtered ?  filteredGraphsByAntPol[pol][ant]->even() : rawGraphsByAntPol[pol][ant]->uneven())->getSumV2();
    }
  }

  return sum/n;

}


double pueo::FilteredEvent::getMedianPower(pol::pol_t pol, ring::ring_t ring, bool filtered) const
{

  pol::pol_t pol_start =pol == pol::kVertical ? pol::kVertical : pol::kHorizontal;
  pol::pol_t pol_end = pol == pol::kHorizontal ? pol::kHorizontal : pol::kVertical;

  int n = k::NUM_ANTS * (int(pol_end) - int(pol_start) +1);

  if (ring != ring::kNotARing) n/=ring::kNotARing;

  double vals[n];


  int i = 0;
  for (int pol = (int) pol_start; pol <=(int) pol_end; pol++ )
  {
    for (int ant = 0; ant < k::NUM_ANTS; ant++)
    {
      if (ring != ring::kNotARing && pueo::GeomTool::getRingFromAnt(ant) != ring) continue;
      vals[i++] = ( filtered ?  filteredGraphsByAntPol[pol][ant]->even() : rawGraphsByAntPol[pol][ant]->uneven())->getSumV2();
    }
  }

  std::nth_element(vals, vals + n/2, vals +n);
  return vals[n/2];
}


void pueo::FilteredEvent::getMedianSpectrum(TGraph * target, pol::pol_t pol, double frac) const
{
  target->Set(131);

  pol::pol_t pol_start =pol == pol::kVertical ? pol::kVertical : pol::kHorizontal;
  pol::pol_t pol_end = pol == pol::kHorizontal ? pol::kHorizontal : pol::kVertical;

//  memset(target->GetY(),0, target->GetN() * sizeof(double));
  memcpy(target->GetX(),filteredGraphsByAntPol[pol_start][0]->power()->GetX(), target->GetN() * sizeof(double));

  int n = k::NUM_ANTS * (int(pol_end) - int(pol_start) +1);


  double vals[n];
  int i = 0;

  //Can paralellize this loop if it's helpful
  for (int j = 0; j < 131; j++)
  {
    i = 0;
    for (int pol = (int) pol_start; pol <=(int) pol_end; pol++ )
    {
      for (int ant = 0; ant < k::NUM_ANTS; ant++)
      {
        //TODO: Is it faster just to put things into a set?
        vals[i++] = filteredGraphsByAntPol[pol][ant]->powerdB()->GetY()[j];
      }
    }
    std::nth_element(vals, vals + int(n*frac), vals +n);
    target->GetY()[j] = vals[int(n*frac)];
  }
}


void pueo::FilteredEvent::getAverageSpectrum(TGraph * target, pol::pol_t pol) const
{
  target->Set(131);

  pol::pol_t pol_start =pol == pol::kVertical ? pol::kVertical : pol::kHorizontal;
  pol::pol_t pol_end = pol == pol::kHorizontal ? pol::kHorizontal : pol::kVertical;

  memset(target->GetY(),0, target->GetN() * sizeof(double));
  memcpy(target->GetX(),filteredGraphsByAntPol[pol_start][0]->power()->GetX(), target->GetN() * sizeof(double));

  int n = k::NUM_ANTS * (int(pol_end) - int(pol_start) +1);

  for (int pol = (int) pol_start; pol <=(int) pol_end; pol++ )
  {
    for (int ant = 0; ant < k::NUM_ANTS; ant++)
    {
      for (int j = 0; j < 131; j++)
      {
        target->GetY()[j] += filteredGraphsByAntPol[pol][ant]->power()->GetY()[j]/n;
      }
    }
  }

  for (int j = 0; j < 131; j++)
  {
    target->GetY()[j] = 10 * TMath::Log10(target->GetY()[j]);
  }


}

void pueo::FilteredEvent::getMinMaxRatio(pol::pol_t pol, double * max_ratio, double * min_ratio, int* max_sector, int* min_sector, ring::ring_t ring1 , ring::ring_t ring2, int nth, int * n_greater )  const
{

  double max = 0;
  double min = 999;
  int imax = -1;
  int imin = -1;

  int greater = 0; 


  for (int i = 0; i < k::NUM_PHI; i++)
  {
    int ant1 = pueo::GeomTool::getAntFromPhiRing(i, ring1);


    int ant2 = pueo::GeomTool::getAntFromPhiRing(i, ring2);


//    printf("%d %d %d\n", pol, ant1,ant2);
//    printf("%p\n", rawGraphsByAntPol[pol][ant1]);
//    printf("%p\n", rawGraphsByAntPol[pol][ant2]);
    double peak1 = rawGraphsByAntPol[pol][ant1]->uneven()->pk2pk(nth,nth);
    double peak2 = rawGraphsByAntPol[pol][ant2]->uneven()->pk2pk(nth,nth);

    double ratio = peak1/peak2;
    if (ratio > 1) greater++; 

    if ( imax < 0 || ratio > max )
    {
      imax = i;
      max = ratio;
    }
    if ( imax < 0 || ratio < min )
    {
      imin = i;
      min = ratio;
    }
  }

  if (max_ratio) *max_ratio = max;
  if (min_ratio) *min_ratio = min;
  if (max_sector) *max_sector= imax;
  if (min_sector) *min_sector= imin;
  if (n_greater) *n_greater = greater; 
}


void pueo::FilteredEvent::plotSummary(TCanvas * ch, TCanvas * cv) const
{

  if (!ch) ch = new TCanvas("pueo::FilteredEvent_hpol","HPol", 1000,1000);
  if (!cv) cv = new TCanvas("pueo::FilteredEvent_vpol","VPol", 1000,1000);

  ch->Clear();
  cv->Clear();

  ch->Divide(8,k::NUM_ANTS/8);
  cv->Divide(8,k::NUM_ANTS/8);


  for (int i =0; i < k::NUM_ANTS; i++)
  {
    ch->cd(i+1);
    getRawGraph(i, pol::kHorizontal)->drawUneven("",1);
    getFilteredGraph(i, pol::kHorizontal)->drawEven("lsame",2);

    cv->cd(i+1);
    getRawGraph(i, pol::kVertical)->drawUneven("",1);
    getFilteredGraph(i, pol::kVertical)->drawEven("lsame",2);
  }

}
int pueo::FilteredEvent::checkSaturation(std::bitset<k::NUM_HORNS> * save_hsat, std::bitset<k::NUM_HORNS> * save_vsat, double thresh) const
{
  
  std::bitset<k::NUM_HORNS> hsat = 0; 
  std::bitset<k::NUM_HORNS> vsat = 0; 

  int totalsat = 0; 
  for (int i = 0; i < k::NUM_ANTS; i++) 
  {
      int hindex = pueo::GeomTool::getChanIndexFromAntPol(i,pol::kHorizontal); 
      int vindex = pueo::GeomTool::getChanIndexFromAntPol(i,pol::kVertical); 
      const double *yh = &useful->volts[hindex][0]; 
      const double *yv = &useful->volts[vindex][0]; 

      for (size_t j = 0; j < useful->volts[hindex].size(); j++)
      {
        if (fabs(yh[j]) > thresh)
        {
          hsat.set(i); 
          totalsat ++; 
          break; 
        }
      }

      for (size_t j = 0; j < useful->volts[vindex].size(); j++)
      {
        if (fabs(yv[j]) > thresh)
        {
          vsat.set(i); 
          totalsat ++; 
          break; 
        }
      }
  }


  if (save_hsat) *save_hsat = hsat; 
  if (save_vsat) *save_vsat = vsat; 

  return totalsat; 
}

int pueo::FilteredEvent::checkStepFunction(Int_t lab, ring::ring_t ring, Int_t phiSector, pol::pol_t pol)  const {
  int ant0 = pueo::GeomTool::getAntFromPhiRing(phiSector, ring);
  if(rawGraphsByAntPol[pol][ant0]->uneven()->pk2pk(0,0) < 1000) return 0;
	for (int i = 0; i < k::NUM_PHI; i++)
  {
    int ant1 = pueo::GeomTool::getAntFromPhiRing(i, ring::kBottomRing);
    int ant2 = pueo::GeomTool::getAntFromPhiRing(i, ring::kLowerMiddleRing);
    int ant3 = pueo::GeomTool::getAntFromPhiRing(i, ring::kUpperMiddleRing);
    int ant4 = pueo::GeomTool::getAntFromPhiRing(i, ring::kTopRing);

    double peak1 = rawGraphsByAntPol[pol][ant1]->uneven()->pk2pk(0,0);
    double peak2 = rawGraphsByAntPol[pol][ant2]->uneven()->pk2pk(0,0);
    double peak3 = rawGraphsByAntPol[pol][ant3]->uneven()->pk2pk(0,0);
    double peak4 = rawGraphsByAntPol[pol][ant4]->uneven()->pk2pk(0,0);
		if((peak1 > 700) && (ant1 != ant0)) return 0;
		if((peak2 > 700) && (ant2 != ant0)) return 0;
		if((peak3 > 700) && (ant3 != ant0)) return 0;
		if((peak4 > 700) && (ant4 != ant0)) return 0;
  }
	return 1;
}

int pueo::FilteredEvent::checkSurfForGlitch(Int_t surf, Int_t lab, double glitchThreshold)  const {
	for (int i = 0; i < k::NUM_DIGITZED_CHANNELS; i++)
  {
		int ant;
		pol::pol_t pol;
		pueo::GeomTool::getAntPolFromSurfChan(surf, i, ant, pol);
		double maxVal = abs(rawGraphsByAntPol[pol][ant]->uneven()->peakVal(0,0,-1,true));
		double p2p = rawGraphsByAntPol[pol][ant]->uneven()->pk2pk();
		double minVal = p2p - maxVal;
		double asymmetry = abs(maxVal - abs(minVal));
		if(asymmetry > glitchThreshold) return 1;
  }
	return 0;
}


const pueo::AnalysisWaveform *pueo::FilteredEvent::getRawGraph(UInt_t phi,
							ring::ring_t ring,
							pol::pol_t pol) const {

    // this works for A4 and should work for A3.
    const Int_t ant = 16*ring + phi;
    return rawGraphsByAntPol[pol][ant];
}




const pueo::AnalysisWaveform * pueo::FilteredEvent::getFilteredGraph(UInt_t phi, 
							      ring::ring_t ring, 
							      pol::pol_t pol) const {
    // this works for A4 and should work for A3.
    const Int_t ant = 16*ring + phi;
    return filteredGraphsByAntPol[pol][ant];
}
