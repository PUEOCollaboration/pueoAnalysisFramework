#ifndef PUEO_NOISE_MONITOR_H
#define PUEO_NOISE_MONITOR_H

#include "TString.h"
#include "TObject.h"
#include <map>

#include "pueo/Conventions.h"


class TFile;
class TTree;
class TProfile2D;
namespace pueo
{

class FilterStrategy;

/** 
 * @class NoiseMonitor is designed to track the channel RMS of min bias events.
 * 
 * Define PUEO_RMS_DIR to generate profile histograms
 */
class NoiseMonitor {

 public:

  // The default timescale in seconds (bin width)
  enum {
    defaultTimeScaleSeconds = 10
  };

  NoiseMonitor(FilterStrategy* fs);
  NoiseMonitor(UInt_t hash);

  virtual ~NoiseMonitor();
  double getRMS(pol::pol_t pol, Int_t ant, UInt_t realTime);

  const TProfile2D* getProfile(pol::pol_t pol, int run){
    getProfilesFromFileRun(run);
    return fCurrent.get(pol);
  }


  /** 
   * @class Overengineered pair holder with set/get methods for pol_t
   * 
   * It shouldn't be public, and only is so ROOT doesn't throw a hissy fit.
   */
  class ProfPair{
   public:
    ProfPair() : H(NULL), V(NULL), fStartTime(0), fEndTime(0) {}
    const TProfile2D* get(pol::pol_t pol) const{
      return pol == pol::kHorizontal ? H : V;
    }
    void set(const TProfile2D* h, const TProfile2D* v);
    void set(const ProfPair& other);
    double startTime(){return fStartTime;}
    double endTime(){return fEndTime;}
   private:
    TProfile2D* H;
    TProfile2D* V;
    double fStartTime;
    double fEndTime;
    
  };

 protected:

  FilterStrategy* fFilterStrat;
  const char* fRmsDir;
  TFile* fFile; // Pointer to currently opened file
  // std::map<int, TFile*> fFiles; /// map of run to opened TFiles
  // std::map<int, ProfPair> fRunProfiles; /// map of run to pair of TProfiles (H/V)

  ProfPair fCurrent; //
  UInt_t fHash;

  // void findProfilesInMemoryFromTime(UInt_t realTime); // from map in memory
  // void findProfilesInMemoryFromRun(Int_t run); // from map in memory
  void getProfilesFromFileRun(int run);  // search for file
  void getProfilesFromFileTime(UInt_t realTime);  // search for file
  void makeProfiles(int run); // last resort, make them from scratch (generates a file)
  void getRmsDirEnv();
  TString getHistName(pol::pol_t pol, int run);
  
  
  
  static UInt_t makeStratHashFromDesc(const FilterStrategy* fs);
  TString getFileName(int run);
  
};

}

#endif //NOISE_MONITOR_H
