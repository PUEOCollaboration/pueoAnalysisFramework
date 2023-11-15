#include "pueo/BasicFilters.h"
#include "pueo/AnalysisWaveform.h" 
#include "DigitalFilter.h" 
#include "FFTWComplex.h"
#include "pueo/ResponseManager.h"

void pueo::SimplePassBandFilter::processOne(AnalysisWaveform* g, const RawHeader * header, int ant, int pol) 
{
  int nfreq = g->Nfreq(); 
  double df = g->deltaF(); 

//  printf("SimplePassBandFilter::processOne!\n"); 
  FFTWComplex * vals = g->updateFreq(); 
  for (int i = 0; i < nfreq; i++) 
  {
    if (i * df < low ||  i * df > high) 
    {
      vals[i].re = 0; 
      vals[i].im = 0; 
    }
  }
}


void pueo::SimpleNotchFilter::processOne(AnalysisWaveform* g, const RawHeader * header, int ant, int pol) 
{
//  printf("SimpleNotchFilter::processOne!\n"); 
  int nfreq = g->Nfreq(); 
  double df = g->deltaF(); 

  FFTWComplex * vals = g->updateFreq(); 
  for (int i = 0; i < nfreq; i++) 
  {
    if (i * df >= min &&  i * df <= max) 
    {
      vals[i].re = 0; 
      vals[i].im = 0; 
    }
  }
}


void pueo::HybridFilter::process(FilteredEvent * event) {

#ifdef USE_OMP
#pragma omp  parallel for 
#endif

  for (int i = 0; i < k::NUM_ANTS; i++) AnalysisWaveform::basisChange(getWf(event,i, pueo::pol::kHorizontal), getWf(event,i,pueo::pol::kVertical));
}


void pueo::HybridFilter::processOne(AnalysisWaveform* g, const RawHeader* header, int ant, int pol) {}


void pueo::SumDifferenceFilter::process(FilteredEvent * event) {

#ifdef USE_OMP
#pragma omp  parallel for 
#endif

  for (int i = 0; i < k::NUM_ANTS; i++) AnalysisWaveform::sumDifference(getWf(event,i, pueo::pol::kHorizontal), getWf(event,i,pueo::pol::kVertical));
}


void pueo::SumDifferenceFilter::processOne(AnalysisWaveform* g, const RawHeader* header, int ant, int pol) {}


void pueo::FlipHVFilter::process(FilteredEvent * event) {

#ifdef USE_OMP
#pragma omp  parallel for 
#endif

  for (int i = 0; i < k::NUM_ANTS; i++) AnalysisWaveform::flipHV(getWf(event,i, pueo::pol::kHorizontal), getWf(event,i,pueo::pol::kVertical));
}


void pueo::FlipHVFilter::processOne(AnalysisWaveform* g, const RawHeader* header, int ant, int pol) {}


pueo::DigitalFilterOperation::DigitalFilterOperation(const FFTtools::DigitalFilter * digi, bool correct, double fmin, double fmax)
  : digi(digi), delay(0) 
{
  if (correct) delay = digi->avgDelay(fmin,fmax,201); //TODO don't hardcode number of points? 

}

void pueo::DigitalFilterOperation::processOne(AnalysisWaveform * wf, const RawHeader * header, int ant, int pol) 
{
  TGraphAligned * g = wf->updateEven(); 

  digi->filterGraph(g); 

  if (delay)
  {
    double dt = (g->GetX()[1]-g->GetX()[0]) * delay; 
    for (int i = 0; i < g->GetN(); i++)
    {
      g->GetX()[i] -= dt; 
    }
  }
}




void pueo::DeconvolveFilter::process(FilteredEvent * ev) 
{

#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for 
#endif
  for (int i = 0; i < 2*k::NUM_ANTS; i++) 
  {
    pueo::pol::pol_t pol = pueo::pol::pol_t( i %2); 
    int ant = i /2; 
    rm->response(pol,ant)->deconvolveInPlace(getWf(ev,ant,pol), dm); 
  }
}

void pueo::DeconvolveFilter::processOne(AnalysisWaveform* awf, const RawHeader* header, int ant, int pol)
{
	printf("processOne not implemented for this yet, sorry!\n");
}



pueo::DeglitchFilter::DeglitchFilter(double th, int n, RemoveAction ac,int mr) 
   :action(ac), thresh(th), neighbors(n), nremoved(0), max_remove(mr) 
{
  descStr.Form("DeglitchFilter with thresh=%g, n_neighbors=%d, action = %s, max_remove = %d", thresh, neighbors, action == DELETE ? "DELETE" : "AVERAGE", max_remove); 
}

static double max_abs(int n, const double * y) 
{
  double max = 0; 

  for (int i =0; i < n; i++)
  {
    if (fabs(y[i]) > max) max = fabs(y[i]); 
  }
  return max; 
}

void pueo::DeglitchFilter::processOne(AnalysisWaveform * wf, const RawHeader * header, int ant, int pol) 
{

  TGraphAligned * g = action == DELETE ? wf->updateUneven() : wf->updateEven(); 

  //circular buffer of values, initialize with first
  int ndelete = 0; 
  int delete_list[max_remove]; // only delete up to 10 values 

  for (int i = 0; i < g->GetN(); i++)
  {
    double max = 0; 

    if (i > 0)
    {
      int start = TMath::Max(0, i - neighbors); 
      int end = i-1; 
      max = max_abs(end-start+1, g->GetY() + start); 
    }

    if (i < g->GetN()-1)
    {
      int start = i+1; 
      int end = TMath::Min(g->GetN()-1, i+neighbors); 
      max = TMath::Max(max,max_abs(end-start+1, g->GetY() + start)); 
    }

    if (g->GetY()[i] > max + thresh)
    {
      if (ndelete > max_remove) 
      {
        fprintf(stderr,"DeglitchFilter: ALREADY REMOVED MORE THAN %d points in his waveform. Giving up.\n", max_remove); 
        return; 
      }

      delete_list[ndelete++] = i; 
    }
  }

  nremoved+=ndelete; 

  for (int i = 0; i < ndelete; i++)
  {
    if (action == DELETE) 
    {
      g->RemovePoint(delete_list[i]-ndelete); 
    }
    else
    {
      g->GetY()[delete_list[i]] =  delete_list[i] == 0 ? g->GetY()[1] : 
                                  delete_list[i] == g->GetN()-1 ? g->GetY()[g->GetN()-2] : 
                                  0.5 * g->GetY()[delete_list[i]-1] + 0.5 * g->GetY()[delete_list[i] +1]; 
    }
  }
}


void pueo::IFFTDiffFilter::processOne(AnalysisWaveform* g, const RawHeader * header, int ant, int pol) {

  int nfreq = g -> Nfreq();

  double dw = 2 * M_PI * g -> deltaF();
  double cosOrder = cos(order * (0.5 + 2 * branchOrder) * M_PI);  //  i^x = cos(x * (0.5 + 2 * k) * pi) + i * sin(x * (0.5 + 2 * k) * pi), k = ..., -1, 0, 1, ...
  double sinOrder = sin(order * (0.5 + 2 * branchOrder) * M_PI);

  FFTWComplex * vals = g -> updateFreq();
  vals[0] = order ? 0 : vals[0];  //  Just in case the FFT itself (order == 0) is desired. But should be zero, anyway.
  for (int i = 1; i < nfreq; i++) {

    double wOrder = pow(i * dw, order);
    vals[i].re = wOrder * (cosOrder * vals[i].re - sinOrder * vals[i].im);
    vals[i].im = wOrder * (sinOrder * vals[i].re + cosOrder * vals[i].im);
  }
  if (fmod(order, 2) && !(nfreq % 2)) vals[nfreq - 1] = 0;  //  "For odd order and even [nfreq], the Nyquist mode is taken zero."

}
