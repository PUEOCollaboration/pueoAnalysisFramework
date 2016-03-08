#include "AnalysisWaveform.h" 
#include "FFTtools.h" 
#include "RFInterpolate.h" 
#include <assert.h>



AnalysisWaveform::AnalysisWaveform(int N, const double *x, const double * y, double nominal_dt, InterpolationType interp_type, InterpolationOptions *opt)
  : g_uneven(N,x,y),  dt(nominal_dt), fft(0), interpolation_type(interp_type), must_update_uneven(false), must_update_freq(true), must_update_even(true), uneven_equals_even(false)  
{

  if (opt) interpolation_options = *opt; 
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true;
}

AnalysisWaveform::AnalysisWaveform(int N, const double * y, double dt, double t0)
  : g_even(N), dt(dt), fft(0),  interpolation_type(AKIMA), must_update_uneven(false), must_update_freq(true), must_update_even(false), uneven_equals_even(true) 
{
  for (int i  = 0; i < N; i++) 
  {
    g_even.GetX()[i] = t0  + dt * i;  
    g_even.GetY()[i] = y[i]; 
  }
  fft_len = N/2 +1; 
  df = 1./(N * dt); 
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true;
}

AnalysisWaveform::AnalysisWaveform(int N, const FFTWComplex * f, double df, double t0)
  :  g_even(N), df(df), interpolation_type(AKIMA), must_update_uneven(false), must_update_freq(false), must_update_even(true), uneven_equals_even(true) 
{
  fft_len = N/2 +1; 
  fft = new FFTWComplex[fft_len]; 
  memcpy(fft, f, fft_len * sizeof(FFTWComplex)); 
  dt = 1./ (N*df); 

  for (int i  = 0; i < N; i++) 
  {
    g_even.GetX()[i] = t0  + dt * i;  
    g_even.GetY()[i] = 0; 
  }

  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
}



const TGraph * AnalysisWaveform::uneven() const
{

  if (uneven_equals_even) return even(); 

  if (must_update_uneven) 
  {
    calculateUnevenFromEven(); 
  }

  return &g_uneven; 
}

const TGraph * AnalysisWaveform::even()  const
{

  //kaboom 
  assert (!(must_update_even && must_update_freq &&  must_update_uneven)); 

  if (must_update_even && must_update_freq) 
  {
    calculateEvenFromUneven(); 
  }
  else if (must_update_even && must_update_even)
  {
    calculateEvenFromFreq(); 
  }
  return &g_even; 
}

const FFTWComplex * AnalysisWaveform::freq()  const
{


  if (must_update_freq) 
  {
    calculateFreqFromEven(); 
  }

  return fft; 
}




void AnalysisWaveform::updateEven(const TGraph * replace)
{
  g_even.Set(replace->GetN()); 
  memcpy(g_even.GetX(), replace->GetX(), replace->GetN() * sizeof(double)); 
  memcpy(g_even.GetY(), replace->GetY(), replace->GetN() * sizeof(double)); 
  must_update_uneven = !uneven_equals_even; 
  must_update_freq = true; 
}


TGraph * AnalysisWaveform::updateEven()
{
  must_update_uneven = !uneven_equals_even; 
  must_update_freq = true; 
  return (TGraph *) even(); 
}

TGraph * AnalysisWaveform::updateUneven()
{
  if (uneven_equals_even) return updateEven(); 

  must_update_even = true; 
  must_update_freq = true; 
  return (TGraph *) uneven(); 
}

FFTWComplex * AnalysisWaveform::updateFreq() 
{
  must_update_even = true; 
  must_update_uneven = !uneven_equals_even; 

  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 

  return (FFTWComplex *) freq(); 
}




void AnalysisWaveform::updateUneven(const TGraph * replace)
{
  if (uneven_equals_even) updateEven(replace); 

  g_uneven.Set(replace->GetN()); 
  memcpy(g_uneven.GetX(), replace->GetX(), replace->GetN() * sizeof(double)); 
  memcpy(g_uneven.GetY(), replace->GetY(), replace->GetN() * sizeof(double)); 

  must_update_even = true; 
  must_update_freq = true; 
}

void AnalysisWaveform::calculateUnevenFromEven() const
{
  const TGraph * g = even(); 

  if (interpolation_type == AKIMA) 
  {
    ROOT::Math::Interpolator irp(g->GetN(), ROOT::Math::Interpolation::kAKIMA); 
    irp.SetData(g->GetN(), g->GetX(), g->GetY()); 
     
    for (int i = 0; i < g_uneven.GetN(); i++) 
    {
      if (g_uneven.GetX()[i] < g->GetX()[g->GetN()-1]) 
        g_uneven.GetY()[i] = irp.Eval(g_uneven.GetX()[i]); 
    }
  }
  else 
  {

    FFTtools::getInterpolatedGraphSparseInvert(g, &g_uneven, interpolation_options.max_distance, 0, 0, 
                                                  interpolation_type == SPARSE_YEN ? 0 : interpolation_options.mu, 
                                                  interpolation_options.regularization_order); 


  }

  must_update_uneven = false; 
}

void AnalysisWaveform::calculateEvenFromUneven()  const
{


  const TGraph * g = uneven(); 

  if (interpolation_type == AKIMA) 
  {
    ROOT::Math::Interpolator irp(g->GetN(), ROOT::Math::Interpolation::kAKIMA); 
    irp.SetData(g->GetN(), g->GetX(), g->GetY()); 
     
    int npoints = (g->GetX()[g->GetN()-1] - g->GetX()[0]) / dt; 
    g_even.Set(npoints); 

    double t0 = g->GetX()[0]; 

    for (int i = 0; i < g_even.GetN(); i++) 
    {
      double x = t0 + i * dt; 
      g_even.GetY()[i] = irp.Eval(x); 
      g_even.GetX()[i] = x;   
    }
  }
  else 
  {
    FFTtools::getInterpolatedGraphSparseInvert(g, &g_uneven, interpolation_options.max_distance, 0, 0, 
                                                  interpolation_type == SPARSE_YEN ? 0 : interpolation_options.mu, 
                                                  interpolation_options.regularization_order); 

  }

  must_update_even = false; 
}

void AnalysisWaveform::calculateEvenFromFreq() const
{
  double *y =  FFTtools::doInvFFT(g_even.GetN(), fft); 

  memcpy(g_even.GetY(), y, g_even.GetN() * sizeof(double)); 
  double t0 = g_even.GetX()[0]; 
  if (g_even.GetX()[1] - t0 != dt) 
  {
    for (int i = 1; i < g_even.GetN(); i++) 
    {
      g_even.GetX()[i] = i * dt + t0; 
    }
  }
  
  delete [] y; 

  must_update_even = false; 
}

void AnalysisWaveform::calculateFreqFromEven() const
{

  const TGraph * g = even();  // in case it must be converted from uneven

  dt = g->GetX()[1] - g->GetX()[0]; 
  fft_len = g->GetN()/2+1;  
  df = 1./ (g->GetN() * dt); 


  if (fft) 
    delete [] fft; 

  fft = FFTtools::doFFT(g->GetN(), g->GetY()); 

  must_update_freq = false; 
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
}



void AnalysisWaveform::updateFreq(int new_N, const FFTWComplex * new_fft, double new_df )
{
  if (new_N && new_N != g_even.GetN())
  {
    delete [] fft; 
    fft = new FFTWComplex[new_N/2+1]; 
    fft_len = new_N/2 + 1; 
    g_even.Set(new_N); 
  }
  memcpy(fft, new_fft, fft_len * sizeof(FFTWComplex)); 
  if (new_df)
  {
    df = new_df; 
  }
  dt = 1./ (g_even.GetN() * df); 

  must_update_even = true; 
  must_update_uneven = !uneven_equals_even; 
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
}

const TGraph * AnalysisWaveform::phase() const
{
  const FFTWComplex * the_fft = freq(); //will update if necessary
  if (phase_dirty)
  {

    g_phase.Set(fft_len); 
    for (int i = 0; i < fft_len; i++) 
    {
        g_phase.SetPoint(i, i * df, the_fft[i].getPhase()); 

    }

    phase_dirty = false;
  }

  return &g_phase; 
}


const TGraph * AnalysisWaveform::powerdB() const
{
  const TGraph * the_power = power(); //will update if necessary
  if (power_db_dirty)
  {

    g_power_db.Set(fft_len); 
    for (int i = 0; i < fft_len; i++) 
    {
        g_power_db.SetPoint(i, the_power->GetX()[i], 10 * TMath::Log10(the_power->GetY()[i])); 
    }

    power_db_dirty = false;
  }

  return &g_power_db; 
}



const TGraph * AnalysisWaveform::power() const
{
  const FFTWComplex * the_fft = freq(); //will update if necessary
  if (power_dirty)
  {

    g_power.Set(fft_len); 
    for (int i = 0; i < fft_len; i++) 
    {
      if (i == 0 || i == fft_len -1)
        g_power.SetPoint(i, i * df, the_fft[i].getAbsSq()/fft_len); 

      else
        g_power.SetPoint(i, i * df, the_fft[i].getAbsSq()*2/fft_len); 
    }

    power_dirty = false;
  }

  return &g_power; 
}

AnalysisWaveform::~AnalysisWaveform() 
{
  if (fft) 
    delete [] fft; 



}



AnalysisWaveform::AnalysisWaveform(const AnalysisWaveform & other) 
{
  

  //first copy over scalars 
  
  dt = other.dt; 
  df = other.df; 
  fft_len = other.fft_len; 

  interpolation_type = other.interpolation_type; 
  interpolation_options = other.interpolation_options; 

  uneven_equals_even = other.uneven_equals_even; 
  must_update_even = other.must_update_even; 
  must_update_freq = other.must_update_freq; 
  must_update_uneven = other.must_update_uneven; 


  //don't bother copying these, they can be regenerated if needed
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 


  // we must copy g_uneven if it's not qual to even 
  if (!uneven_equals_even)
  {
    g_uneven = other.g_uneven; 
  }
 

  if (!must_update_even)
  {
    g_even = other.g_even; 
  }
  else // still need the size! 
  {
    g_even.Set(other.g_even.GetN()); 
  }

  if (!must_update_freq)
  {
    memcpy(fft, other.fft, fft_len * sizeof(FFTWComplex)); 
  }
  else
  {
    fft = 0; 
  }

  



}

AnalysisWaveform * AnalysisWaveform::correlation(const AnalysisWaveform *A, const AnalysisWaveform *B) 
{
  if (A->Nfreq() != B->Nfreq()) 
  {
    fprintf(stderr,"correlation does not handle the case where A and B are of different lengths!\n"); 
    return 0; 
  }

  // in principle, we could check if even is zeropadded to warn people but it's kinda expensive


  AnalysisWaveform * answer = new AnalysisWaveform(A->Neven(), A->freq(), 
                                                   A->deltaF(), A->even()->GetX()[0]); 

  FFTWComplex * update = answer->updateFreq(); 

  for (int i = 0; i < B->Nfreq(); i++) 
  {
    FFTWComplex vA = update[i]; 
    FFTWComplex vB = B->freq()[i]; 
    update[i].re =  vA.re * vB.re + vA.im * vB.im ; 
    update[i].im =  vA.im * vB.re - vA.re * vB.im; 
  }

  return answer; 
}

void AnalysisWaveform::padFreq(int npad)
{
  //new even size
  int new_N = even()->GetN() * npad; 

  //allocate new memory
  FFTWComplex * new_fft = new FFTWComplex[new_N/2+1]; 

  //copy old
  memcpy(new_fft, freq(), Nfreq() * sizeof(FFTWComplex));

  //scale
  FFTWComplex scale(npad,0); 
  for (int i =0; i < Nfreq(); i++) { new_fft[i] *= scale; }


  // zero rest 
  memset(new_fft + fft_len, 0, (new_N/2 + 1 - fft_len) * sizeof(FFTWComplex));

  //set new lengths and dt
  fft_len = new_N/2 + 1; 
  g_even.Set(new_N); 
  dt = 1./ (g_even.GetN() * df); 


  //swap in new fft
  delete [] fft; 
  fft = new_fft; 

  //others need to update now! 
  must_update_even = true; 
  must_update_uneven = !uneven_equals_even; 
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
}

void AnalysisWaveform::padEven(int npad)
{
  TGraph * g = updateEven(); 
  int old_n = g->GetN(); 
  g->Set(g->GetN() * npad); 
  for (int i = old_n; i < g->GetN(); i++) 
  {
    g->GetX()[i] = g->GetX()[i-1] + dt; 
    g->GetY()[i] = 0; 
  }
}

void AnalysisWaveform::drawEven(const char * opt) const{ ((TGraph*)even())->Draw(opt); }
void AnalysisWaveform::drawUneven(const char * opt) const{ ((TGraph*)uneven())->Draw(opt); }
void AnalysisWaveform::drawPower(const char * opt)const { ((TGraph*)power())->Draw(opt); }
void AnalysisWaveform::drawPowerdB(const char * opt)const { ((TGraph*)powerdB())->Draw(opt); }
void AnalysisWaveform::drawPhase(const char * opt) const{ ((TGraph*)phase())->Draw(opt); }

