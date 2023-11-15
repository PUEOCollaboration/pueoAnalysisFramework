#ifndef PUEO_TEMPLATES
#define PUEO_TEMPLATES

#include "TObject.h"
#include "pueo/EventSummary.h"
#include "FFTtools.h"
#include "TFile.h"
#include "pueo/AnalysisWaveform.h"
#include "pueo/SystemResponse.h"
#include "pueo/Version.h"


/* * * * * * * * *
   
   
   Be careful using templates, we might not know what we're looking for


   Anyway this is where all the template stuff will live


 * * * * * * * * * */


/**

   _emplateMachine_

    Class to import and hold onto the different templates
    
    Saveable to a ROOT file I think if you want to do that sort of thing!
*/


namespace pueo 
{
class TemplateSummary;

class TemplateMachine : public TObject {

 public:
  /** Default Constructor **/
  TemplateMachine(const int inLength=2048);

  /** Default Destructor **/
  ~TemplateMachine();
  
  /** in case you want to reset for some reason **/
  void zeroInternals();

  /** template waveform defaults **/
  const int length;
  const double dT;// = 0.1; //in ns
  
  /** Template Storage **/
  //impulse response
  FFTWComplex *theImpTemplateFFT;
  TGraph *theImpTemplate;
  //ANITA3 averaged wais pulses
  FFTWComplex *theWaisTemplateFFT;
  TGraph *theWaisTemplate;
  //ZHAires simulated CRs at various angles
  static constexpr int numCRTemplates = 10; //this is the place it is defined for everything
  FFTWComplex *theCRTemplateFFTs[numCRTemplates];
  TGraph *theCRTemplates[numCRTemplates];

  //deconvolved versions for all those
  //impulse response
  FFTWComplex *theImpTemplateFFT_deconv;
  TGraph *theImpTemplate_deconv;
  //ANITA3 averaged wais pulses
  FFTWComplex *theWaisTemplateFFT_deconv;
  TGraph *theWaisTemplate_deconv;
  //ZHAires simulated CRs at various angles
  FFTWComplex *theCRTemplateFFTs_deconv[numCRTemplates];
  TGraph *theCRTemplates_deconv[numCRTemplates];
  

  /** Filling the stored templates **/
  void loadTemplates( unsigned int evTime = 0, int version = pueo::version::get());
  void deconvolveTemplates(DeconvolutionMethod *deconv);
  

  /** Flags */
  //Check to see if the templates are loaded in bad english
  bool isTmpltsLoaded() { return kTmpltsLoaded; };
  //Check to see if you filled the deconvolved templates too
  bool isTmpltsDeconv() { return kTmpltsDeconv; };

  /** use the templates, and a supplied AnalysisWaveform, to fill up the summary */
  void doTemplateAnalysis(const AnalysisWaveform *waveform, int poli, int dir, TemplateSummary *summary, bool do_impulse = true, bool do_wais = true, bool do_cr = true);
  void doDeconvolvedTemplateAnalysis(const AnalysisWaveform *waveform, const DeconvolutionMethod *deconv, 
				     int poli, int dir, TemplateSummary *summary,bool do_impulse = true, bool do_wais = true, bool do_cr = true );


  /** Write templates to file **/
  void writeTemplatesToFile(TFile *outFile);

  void setUseAverageCRTemplate(bool opt) { fUseAverageCRTemplate = opt; }
  void setDoWindow(bool opt) { fDoWindow = opt; }
  std::string getNotchStr() { return fNotchStr ;}
  
 private:
  /* flags to see what you still might need to do */
  bool kTmpltsLoaded;
  bool kTmpltsDeconv;
  bool fUseAverageCRTemplate;
  bool fDoWindow;

  std::string fNotchStr;
  std::vector<int> payloadTimes;
  std::vector<std::string> notchConfigs;

  void fillNotchConfigs();
  
  /* hidden functions for filling the templates one at a time */
  void getImpulseResponseTemplate(int version);
  void getWaisTemplate(int version);
  void getCRTemplates(int version);
  
  
  ClassDefNV(TemplateMachine, 1);
  
};


/**

   _AnitaTemplateSummary_

    Class to store results of any template correlations

*/

class TemplateSummary
{
 public:

  /** Default Constructor **/
  TemplateSummary();

  /** Default Destructor **/
  virtual ~TemplateSummary();

  /** The maximum number of hypotheses storable per polarization */ 
  static const Int_t maxDirectionsPerPol = EventSummary::maxDirectionsPerPol; 

  /* Number of points on the cone */
  static const int numCRTemplates = TemplateMachine::numCRTemplates;

  /*The template correlatin values for a single coherent waveform comparison*/
  class SingleTemplateResult 
  {
  public:
    SingleTemplateResult() {; }
    //impulse response
    Double_t impulse; //peak
    Double_t impulse_loc; //location
    bool impulse_pol; //max(1) or min(0)

    //one for the WAIS template too
    Double_t wais;
    Double_t wais_loc;
    bool wais_pol;
    
    //and for the bigger multi-coherence-angle one
    Double_t cRay[numCRTemplates];
    Double_t cRay_loc[numCRTemplates];
    bool cRay_pol[numCRTemplates];
    

    ClassDefNV(SingleTemplateResult,3);
  };

  
  SingleTemplateResult coherent[k::NUM_POLS][maxDirectionsPerPol];
  
  SingleTemplateResult deconvolved[k::NUM_POLS][maxDirectionsPerPol];


  void zeroInternals();
    
 private:
  ClassDefNV(TemplateSummary, 1); 
};





//-----------------
//random utilities
// vvvvvvvvvvvv


/**

   Stuff I use that should maybe be in FFTtools?
   
*/

namespace FFTtoolsAnnex{
  TGraph *normalizeWaveform(TGraph* inGraph);
  double *getCorrelationFromFFT(int length,const FFTWComplex *theFFT1, const FFTWComplex *theFFT2);
}







/**

   windowing functions!

*/

namespace WindowingTools {

  TGraph *windowWave(TGraph*, int&, const int, const int, const int, const int);
  
  
  TGraph *windowCut(TGraph *inGraph,int length);
  
  TGraph *windowDispersed(TGraph *inGraph, int &peakHilbertLoc);
  TGraph *windowDispersed(TGraph *inGraph);
  
  TGraph *windowEField(TGraph *inGraph, int &peakHilbertLoc);
  TGraph *windowEField(TGraph *inGraph);
  

  TGraph *windowWave(TGraph *inGraph, int &peakHilbertLoc,
		     const int upRampLoc = 50,   const int downRampLoc = 600,
		     const int upRampLen = 100, const int downRampLen = 400);

}


}




#endif
