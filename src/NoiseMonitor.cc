#include "pueo/NoiseMonitor.h"
#include "pueo/FilteredEvent.h"
#include "pueo/RawHeader.h"
#include "pueo/FilterStrategy.h"
#include "pueo/FilterOperation.h"
#include "pueo/Dataset.h"
#include "pueo/Conventions.h" 
#include <numeric>

#include "TTree.h"
#include "TFile.h"
#include "TProfile2D.h"
#include <stdlib.h>

void pueo::NoiseMonitor::getRmsDirEnv(){
  fRmsDir = getenv("PUEO_RMS_DIR");
  if(!fRmsDir){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__
              << ", can't see env PUEO_RMS_DIR. "
              << "Will look for/make rms profiles in pwd"
              << std::endl;
    fRmsDir = ".";
  }
}

UInt_t pueo::NoiseMonitor::makeStratHashFromDesc(const FilterStrategy* fs){
  TString hasher;
  for(unsigned i=0; i < fs->nOperations(); i++){
    hasher += fs->getOperation(i)->description();
  }
  UInt_t hash = hasher.Hash();
  return hash;
}


TString pueo::NoiseMonitor::getFileName(int run){
  return TString::Format("%s/rms_pueo%d_run%d_hash%u.root", fRmsDir, version::get(), run, fHash);
}

void pueo::NoiseMonitor::getProfilesFromFileRun(int run){

  TString theRootPwd = gDirectory->GetPath();
  TString fileName = getFileName(run);
  if(fFile){
    fFile->Close();
  }
  fFile = TFile::Open(fileName);

  ProfPair p;
  if(fFile){
    TProfile2D* pH = (TProfile2D*) fFile->Get(getHistName(pol::kHorizontal, run));
    TProfile2D* pV = (TProfile2D*) fFile->Get(getHistName(pol::kVertical, run));
    p.set(pH, pV);
  }
  else{
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", couldn't find file " << fileName << ", will need to generate it. This may take some time." << std::endl;
  }

  if(!p.get(pol::kVertical) || !p.get(pol::kHorizontal)){
    std::cerr << "Error! unable to get the profiles I just tried to read for run " << run << ", "
              << p.get(pol::kVertical) << ", " << p.get(pol::kHorizontal) << std::endl;
    makeProfiles(run);

    fFile = TFile::Open(fileName);
    TProfile2D* pH = (TProfile2D*) fFile->Get(getHistName(pol::kHorizontal, run));
    TProfile2D* pV = (TProfile2D*) fFile->Get(getHistName(pol::kVertical, run));
    p.set(pH, pV);
  }

  if(!p.get(pol::kVertical) || !p.get(pol::kHorizontal)){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__
              << ", unable to find RMS profiles for run " << run
              << ". Something bad is about to happen!"
              << std::endl;
  }
  // fFiles[run] = fFile;

  fCurrent.set(p);

  gDirectory->cd(theRootPwd);
}


TString pueo::NoiseMonitor::getHistName(pol::pol_t pol, int run){
  TString name = TString::Format("rms?_%d_%u", run, fHash);
  const Ssiz_t namePolCharPos = name.First('?');
  name[namePolCharPos] = pol::asChar(pol);
  return name;
}


void pueo::NoiseMonitor::makeProfiles(int run){

  TString fileName = getFileName(run);

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ", generating minimim bias RMS profiles for run " << run << " in file " << fileName << std::endl;
  TFile* f = TFile::Open(fileName, "recreate");

  Dataset d(run);
  const int n = d.N();
  RawHeader* header = NULL;

  d.first();
  header = d.header();
  double startTime = header->realTime;

  d.last();
  header = d.header();
  double endTime = header->realTime + 1; // make sure last event is in bin

  const int nBin = (endTime - startTime)/double(defaultTimeScaleSeconds);

  // temp turn off output, turned back on at the end of the function
  std::vector<bool> original_enable_outputs = fFilterStrat->enable_outputs;
  for(unsigned i=0; i < fFilterStrat->enable_outputs.size(); i++){
    fFilterStrat->enable_outputs.at(i) = false;
  }

  TString name = TString::Format("rms?_%d_%u", run, fHash);
  TString title = TString::Format("RMS for each ?Pol for run %d (FilterStrategy operations descriptions hash = %u);Antenna;realTime", run, fHash);

  const Ssiz_t titlePolCharPos = title.First('?');

  std::vector<TProfile2D*> profs(pol::kNotAPol, NULL);
  for(int polInd=0; polInd < pol::kNotAPol; polInd++){
    pol::pol_t pol = (pol::pol_t) polInd;

    TString name = getHistName(pol, run);

    title[titlePolCharPos] = pol::asChar(pol);
    profs.at(pol) = new TProfile2D(name, title, k::NUM_ANTS, -0.5, k::NUM_ANTS-0.5, nBin, startTime, endTime);
  }

  const double deltaPrint = double(n)/1000;
  double nextPrint = 0;

  for(int entry=0; entry < n; entry++){
    d.getEntry(entry);
    header = d.header();
    if(!trigger::isRFTrigger(header->trigType)){

      FilteredEvent fEv(d.useful(), fFilterStrat, d.gps(), header);
      for(int polInd=0; polInd < pol::kNotAPol; polInd++){
        pol::pol_t pol = (pol::pol_t) polInd;

        for(int ant=0; ant < k::NUM_ANTS; ant++){
          const AnalysisWaveform* wf = fEv.getFilteredGraph(ant, pol);
          const TGraphAligned* gr = wf->even();
          double rms = TMath::RMS(gr->GetN(), gr->GetY());
          profs.at(pol)->Fill(ant, header->realTime, rms);
        }
      }
    }
    if(entry >= nextPrint){
      const int nm = 50;
      int m = nm*nextPrint/n;
      fprintf(stderr, "\r%4.2f %% complete", 100*nextPrint/n);
      std::cerr << "[";
      for(int i=0; i < m; i++) std::cerr << "=";
      for(int i=m; i < nm; i++) std::cerr << " ";
      std::cerr << "]";
      nextPrint += deltaPrint;
      if(nextPrint >= n){std::cerr << std::endl;}      
    }
    
    // if(entry == nextPrint){
    //   std::cerr << "\r" << 100*(entry+1)/n << "% ";
    //   nextPrint += printEvery;
    //   nextPrint = nextPrint >= n ? n - 1 : nextPrint;
    //   if(entry==n-1){
    //     std::cerr << std::endl;
    //   }
    // }
  }

  // turn the outputs back on, if they were on
  fFilterStrat->enable_outputs = original_enable_outputs;

  f->Write();
  f->Close();

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ", complete!" << std::endl;
}



/** 
 * Default constructor, requires a filter strategy
 * If you pass it a NULL, things probably won't end well
 * 
 * @param fs is the filter strategy we want to use to get the waveform RMSs
 */
pueo::NoiseMonitor::NoiseMonitor(FilterStrategy* fs)
    : fFilterStrat(fs), fFile(NULL)
{
  getRmsDirEnv();
  fHash = makeStratHashFromDesc(fs);
}



/** 
 * Default constructor, requires a filter strategy
 * If you pass it a hash which doesn't match something in PUEO_RMS_DIR, it won't end well
 * 
 * @param fs is the filter strategy we want to use to get the waveform RMSs
 */
pueo::NoiseMonitor::NoiseMonitor(UInt_t hash)
    : fFilterStrat(NULL), fFile(NULL)
{
  getRmsDirEnv();
  fHash = hash;
}




/** 
 * Destructor.
 */
pueo::NoiseMonitor::~NoiseMonitor(){

  // std::map<int, TFile*>::iterator it;
  // for(it = fFiles.begin(); it != fFiles.end(); it++){
  //   it->second->Close();
  // }
  if(fFile){
    fFile->Close();
    fFile = NULL;
  }
}



// /** 
//  * Set fCurrent if found,
//  * 
//  * @param realTime is used to find the run, to find the profile in memory
//  */
// void NoiseMonitor::findProfilesInMemoryFromTime(UInt_t realTime){

//   int run = AnitaDataset::getRunAtTime(double(realTime));
//   findProfilesInMemoryFromRun(run);
// }



// /** 
//  * Set fCurrent if found,
//  * 
//  * @param realTime is used to find the run, to find the profile in memory
//  */
// void NoiseMonitor::findProfilesInMemoryFromRun(Int_t run){

//   // first look in the map
//   std::map<int, ProfPair>::iterator it = fRunProfiles.find(run);
//   if(it!=fRunProfiles.end()){
//     // then this the map
//     fCurrent.set(it->second);
//   }
//   else{
//     getProfilesFromFileRun(run);
//   }
// }


void pueo::NoiseMonitor::getProfilesFromFileTime(UInt_t realTime){
  int run = Dataset::getRunAtTime(double(realTime));
  getProfilesFromFileRun(run); // update fCurrent
}


double pueo::NoiseMonitor::getRMS(pol::pol_t pol, Int_t ant, UInt_t realTime){

  version::setVersionFromUnixTime(realTime);

  const TProfile2D* p = fCurrent.get(pol);
  
  if(!(p && realTime >= fCurrent.startTime() && realTime < fCurrent.endTime())){
    getProfilesFromFileTime(realTime); // update fCurrent
    p = fCurrent.get(pol); // should get p again
  }

  // surely find bin should be const?
  double rms = 0;
  if(p){
    Int_t bin = ((TProfile2D*)p)->FindBin(ant, realTime);
    rms = p->GetBinContent(bin);
  }
  return rms;
}


void pueo::NoiseMonitor::ProfPair::set(const TProfile2D* h, const TProfile2D* v){

  H = const_cast<TProfile2D*>(h);
  V = const_cast<TProfile2D*>(v);
  // use both H and V so if one is NULL it fucks up
  // presumably the limits are the same.. but I'm not going to check right now
  fStartTime = H->GetYaxis()->GetBinLowEdge(1);
  fEndTime = V->GetYaxis()->GetBinLowEdge(V->GetNbinsY());

  // prettiness
  if(H && V){
    H->GetYaxis()->SetTimeDisplay(1);
    V->GetYaxis()->SetTimeDisplay(1);

    const char* timeFormat = "#splitline{%H:%M}{%d/%m/%Y}";
    H->GetYaxis()->SetTimeFormat(timeFormat);
    V->GetYaxis()->SetTimeFormat(timeFormat);

    H->GetYaxis()->SetTimeOffset(0);
    V->GetYaxis()->SetTimeOffset(0);

    H->GetYaxis()->SetLabelSize(0.02);
    V->GetYaxis()->SetLabelSize(0.02);

    fStartTime = H->GetYaxis()->GetBinLowEdge(1);
    fEndTime = V->GetYaxis()->GetBinLowEdge(V->GetNbinsY());
  }
  else{
    // something impossible
    fStartTime = 0;
    fEndTime = 0;
  }

  // H->GetXaxis()->LabelsOption("h");
  // V->GetXaxis()->LabelsOption("h");
}

void pueo::NoiseMonitor::ProfPair::set(const ProfPair& other){

  set(other.get(pol::kHorizontal), other.get(pol::kVertical));
  
}
