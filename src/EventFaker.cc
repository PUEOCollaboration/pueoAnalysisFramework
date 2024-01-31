#include "pueo/EventFaker.h" 
#include "pueo/UsefulEvent.h" 
#include "pueo/SystemResponse.h" 
#include "pueo/GeomTool.h" 
#include "FFTtools.h" 
#include "TRandom3.h" 
#include <fftw3.h> 

const TF2 & pueo::EventFaker::getDefaultHPolAntennaGain() 
{
  static TF2 f ("default_hpol_antenna_gain", "exp(-x*x/(2*(20*20))) * exp(-y*y/(2*(20*20)))", -90,90,-90,90); 
  return f;
}


const TF2 & pueo::EventFaker::getDefaultVPolAntennaGain() 
{
  static TF2 f ("default_vpol_antenna_gain", "exp(-x*x/(2*(20*20))) * exp(-y*y/(2*(20*20)))", -90,90,-90,90); 
  return f; 
}

const TF2 & pueo::EventFaker::getDefaultOffAxisDelay() 
{
  const double c2 = 1.45676e-8; 
  const double c1 = 5.01452e-6; 
  static TF2 f("default_off_axis_delay", "[0] * (x*x+y*y) + [1] * (x*x*x*x + y*y*y*y)", -90,90,-90,90); 
  f.SetParameter(0,c1); 
  f.SetParameter(1,c2); 
  return f; 
}

const TF1 & pueo::EventFaker::getDefaultMagResponse() 
{
  static TF1 f("default_mag_response", "1",0,2); 
  return f; 
}

pueo::EventFaker::EventFaker(const char * responseDir, const TF1 & mag_response , double delay , double signal_dt) 
      : manager( responseDir, ceil((1./2.6)/ signal_dt)), 
        offAxisGain_hpol(getDefaultHPolAntennaGain()), 
        offAxisGain_vpol(getDefaultVPolAntennaGain()), 
        offAxisDelay(getDefaultOffAxisDelay())
    {
      setSignalFromMagnitudeResponse(mag_response, delay, signal_dt); 
    }

static void normalizeSignal(pueo::AnalysisWaveform * sig) 
{
  double scale = sig->even()->getSumV2(); 
  pueo::TGraphAligned * g = sig->updateEven(); 
  for (int i = 0; i < g->GetN(); i++)
  {
    g->GetY()[i] /= scale; 
  }

}

void pueo::EventFaker::setSignal(const AnalysisWaveform & sig) 
{
  prototype = sig; 
  normalizeSignal(&prototype); 
  prototype.padEven(0); 

  for (int ant = 0; ant < k::NUM_ANTS; ant++)
  {
    for (int ipol = 0; ipol < 2; ipol++) 
    {
      signal[ant][ipol] = prototype; 
      manager.response(ipol,ant)->convolveInPlace(&signal[ant][ipol]); 
    }
  }
}


template <class Evalable>  
static pueo::AnalysisWaveform * makeSignalFromMagnitudeResponse( const Evalable  & mag_response, double delay, double signal_dt, int N ) 
{
  //find the minimum df in the magnitude response
  double sig_df = 1 / (N*signal_dt); 
  std::vector <double> G(N/2+1); 

  for (int i = 0; i < N/2+1; i++) 
  {
    G[i]=mag_response.Eval(sig_df * i); 
  }

  FFTWComplex * F = FFTtools::makeMinimumPhase(N/2+1, &G[0]); 

  pueo::AnalysisWaveform * answer =  new pueo::AnalysisWaveform (N, F, sig_df, delay); 

  delete [] F; 

  return answer; 
}


void pueo::EventFaker::setSignalFromMagnitudeResponse(const TF1 & mag_response, double delay , double signal_dt, int npoints ) 
{
  AnalysisWaveform * tmp = makeSignalFromMagnitudeResponse( mag_response, delay, signal_dt, npoints); 
  setSignal(*tmp); 
  delete tmp; 

}

void pueo::EventFaker::setSignalFromMagnitudeResponse(const TGraph& mag_response, double delay , double signal_dt, int npoints ) 
{
  AnalysisWaveform * tmp = makeSignalFromMagnitudeResponse( mag_response, delay, signal_dt, npoints); 
  setSignal(*tmp); 
  delete tmp; 
}


void pueo::EventFaker::addSignal(pueo::UsefulEvent * event, double theta, double phi, double A, double extra_delay, 
                                std::complex<double> jones_H, std::complex<double> jones_V) const
{
  

  double th_rad = theta * TMath::DegToRad(); 
  double phi_rad = phi * TMath::DegToRad(); 

  auto geom = GeomTool::Instance(); 

  for (int ich = 0; ich < k::NUM_DIGITZED_CHANNELS; ich++)
  {
    int chan, surf; 
    geom.getSurfChanFromChanIndex(ich, surf, chan); 
    if (chan == 8) continue; //this is the clock, ignore it; 
    int ant; 
    pueo::pol::pol_t pol ; 
    geom.getAntPolFromSurfChan(surf,chan, ant, pol); 
    double sig_t0 = signal[ant][pol].even()->GetX()[0]; 
    double sig_t1 = signal[ant][pol].even()->GetX()[signal[ant][pol].Neven()-1]; 

    //we need to figure out the gain and delay for this channel 
    
    double R = geom.getAntR(ant, pol); 
    double z = geom.getAntZ(ant, pol); 
    double phi0_rad =  geom.getAntPhiPositionRelToAftFore(ant,pol); 
    double phi0=  phi0_rad  * TMath::RadToDeg(); 
    double Greal = A * (pol == pueo::pol::kHorizontal ?  offAxisGain_hpol :  offAxisGain_vpol).Eval(phi-phi0, theta-10); //TODO: don't hardcode 10 
    std::complex<double> G = Greal * (pol == pueo::pol::kHorizontal ? jones_H : jones_V); 

    double ts= (z * tan(th_rad) - R * cos(phi_rad-phi0_rad)) *  1e9 * cos(th_rad) / C_LIGHT; 
    ts+= offAxisDelay.Eval( phi-phi0, theta-10); 
    ts += extra_delay; 

    for (size_t i = 0; i < event->volts[ich].size(); i++) 
    {

      if (event->t(ich,i) - ts< sig_t0) continue; 
      if (event->t(ich,i) - ts> sig_t1) continue; 

      std::complex<double> val = G * std::complex<double> ( signal[ant][pol].evalEven(event->t(ich,i) - ts),
                                                            signal[ant][pol].hilbertTransform()->evalEven(event->t(ich,i) - ts)); 
      event->volts[ich][i] += std::real(val); 
    }
  }
}





