#ifndef ANITA_ANALYSIS_BASIC_FILTERS
#define ANITA_ANALYSIS_BASIC_FILTERS

#include "FilterOperation.h"
#include "TString.h" 
#include "FFTWComplex.h"
#include "AnitaGeomTool.h"
#include "FilteredAnitaEvent.h"
#include <vector>
using namespace std;

/** A set of basic filter operations that serve as examples and also should be quite useful */ 


/** SimplePassBandFilter just cuts everything outside the pass band in
 * fourier space like a brick wall, with all the ensuing acausal goodness.
 *
 */ 


/** Geometric filter */
class GeometricFilter : public FilterOperation
{
  public: 
    typedef struct {} notchPlaceholder;  //  whatever is your notch definition data structure
    GeometricFilter() {}                 // dummy ctor for development; will be deprecated
    GeometricFilter(std::vector< vector<AnalysisWaveform*> > noise, AnitaGeomTool* geom)  
      : noiseSamples(noise), geomTool(geom) // initialize variables here as needed;   
      {
        descStr.Form("Geometric Filter"); 
      }

    // process an event
    virtual void process(FilteredAnitaEvent * event);
    virtual void processOne(AnalysisWaveform* wf); 
    void setDbCut(double dB) {dbCut = dB;}
    const char * tag() const { return "GeomFilter"; } 
    const char * description() const { return descStr.Data(); } 
    int printFlag = 0;

  protected:
    void getNotchandBandwidth(int nfreq, double dBCut, int nAntennasToUse, int nadirFlag, float meanFreqVert, float meanFreqHoriz);
    void getGroupsofAntennas(int nAntennasToUse, int nadirFlag);
    void getClosestNAntennas(int nantennasToUse, double peakPhi, vector<int>& whichAntennasToUse, int nadirFlag);
    void adaptiveFilterPartialPayload(int pol, double dBCut, int nfreq,double *frequencies,double *bandwidth,double *magPeak,
            int nantennasToUse,vector<int>& whichAntennasToUse, float &mean_freq);
    void getFrequenciestoCut(int antenna,vector< vector<double> > &antennaFreq,vector< vector<double> > &bandwidth, 
            vector< vector<double> > &PeakMag, vector<double> &uniquefreqs,vector<double> &uniquebandwidth, int nfreq, 
            vector<double> &uniquePhase, vector<double> &uniquePhase_bandwidth);
    void applyAdaptiveFilter_singleAnt(double centerFrequency, double bandWidth, int polFlag,int ant);               
    TGraph* interpolatedFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq); //interpolated filter. Do interpolation across notch    
    void GeomMethod(int ant,int pol,vector<double> Freq,vector<double> bandWidth,vector<double> cutFreqs); //geom method for phase filtering  
    double solveGamma_plus(double theta, double psi, double delta);
    double solveGamma_minus(double theta, double psi, double delta);      
  private:
    std::vector<vector<AnalysisWaveform*> > noiseSamples;         // populated in constructor
    TString descStr;
    int nFreq = 3;
    double dbCut = 2.0;
    int nAntsToUse = 9;
    int groupFlag = 0;
    std::vector<std::vector<double> > notchFreqs;        
    std::vector< std::vector<double> > notchBands;
    std::vector< std::vector<double> > notchFreqs_Horiz;
    std::vector< std::vector<double> > notchBands_Horiz;
    std::vector< std::vector<double> > notchPhase;
    std::vector< std::vector<double> > notchPhaseBands;
    std::vector< std::vector<double> > notchPhase_Horiz;
    std::vector< std::vector<double> > notchPhaseBands_Horiz;
    vector< vector<int> > antenna_group_holder;
    int antenna_groups_size;
    vector<int> unique_phis;
    AnitaGeomTool* geomTool = 0;
    FilteredAnitaEvent* currentEvent = 0;
    double peakPhi=0.;
    int saturatedChannels[2*NUM_SEAVEYS] = {0};
}; 


#endif 