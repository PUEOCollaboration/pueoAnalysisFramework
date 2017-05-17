#include "ResponseManager.h" 
#include "SystemResponse.h" 
#include <stdio.h>
#include "AnalysisWaveform.h" 
#include <sys/types.h>
#include <dirent.h>
#include <assert.h>
#include "FFTtools.h"
#include "TGraph.h"
#include "AnitaGeomTool.h"


/** See README.response */ 




int AnitaResponse::ResponseManager::loadResponsesFromDir(const char * raw_dir, int npad)
{
  DIR *dp; 
  struct dirent *ep; 

  //where do we look for the dir? 
  // Then try data/responses/ 
  // Then try ${ANITA_UTIL_INSTALL_DIR}/share/UCorrelator/responses
  
  TString dir ; 

  dir.Form(raw_dir); 
  dp = opendir(dir.Data()); 
  if (!dp)
  {
    dir.Form("./data/responses/%s",raw_dir); 
    dp = opendir(dir.Data()); 
    if (!dp)
    {
      dir.Form("%s/share/AnitaAnalysisFramework/responses/%s", getenv("ANITA_UTIL_INSTALL_DIR"), raw_dir); 
      dp = opendir(dir.Data()); 

      if (!dp)
      {
          fprintf(stderr,"Could not open response dir %s\n",raw_dir);
          return 1; 
      }
    }
  }


   

  std::map<const char *, AbstractResponse *> prefix_map; 

  while ( (ep = readdir(dp)) )
  {

    char * entry = ep->d_name; 
    char * prefix = strdup(entry); 

    if (entry[0]=='.') continue; 

//    printf("ResponseManager found: %s\n",entry); 
    // find dot
   
    char * dot = strstr(prefix,"."); 
    if (!dot) 
    {
      fprintf(stderr, "Entry %s contains no .\n", entry); 
      continue; 
    }

    //pointer to the suffix
    char * suffix = dot+1; 

    //truncate the prefix
    *dot = 0;

    double angle = 0; 

    char * angstr=strstr(prefix,"_"); 

    //get angle (search for _) 
    if (angstr) 
    {
      angle = atof(angstr+1); //parse the angle
      *prefix = 0;  //chop it off from the prefix
    }
    
    AnitaResponse::AbstractResponse * r; 

    // do we have something with the same prefix but different angle? If so, we should add to it 
    if (prefix_map.count(prefix))
    {
      r = prefix_map[prefix]; 
    }
    else
    {
      r = 0; 
    }


    //now create a response from this 
    if (!strcasecmp(suffix,"imp"))
    {
      TString fname;
      fname.Form("%s/%s", dir.Data(), entry); 
      TGraph imp(fname); 
      if (!r) 
      {
        r = new Response(&imp, npad); 
        response_store.push_back(r); 
      }
      else
      {
        assert(imp.GetN()); 
        AnalysisWaveform aw (imp.GetN(), imp.GetY(), imp.GetX()[1] - imp.GetX()[0], imp.GetX()[0]); 
        aw.padEven(npad); 
        ((AnitaResponse::Response*) r)->addResponseAtAngle(angle, aw.freq()); 
      }
    }

    else if (!strcasecmp(suffix,"freq"))
    {
      // TODO writer parser

      fprintf(stderr,"Reading of .freq format not implemented yet!\n"); 

    }

    else if (!strcasecmp(suffix,"iir"))
    {
      //TODO write parser 
      fprintf(stderr,"Reading of .iir format not implemented yet!\n"); 
    }



    //now figure out what to apply it to


    int start_ant, stop_ant, start_pol, stop_pol;

    if (!strcasecmp(prefix,"all"))
    {
      start_pol = 0; start_ant = 0; 
      stop_pol = 1; stop_ant = NUM_SEAVEYS-1; 
    }

    else if (!strcasecmp(prefix,"allHpol"))
    {
      start_pol = 0; start_ant = 0; 
      stop_pol = 0; stop_ant = NUM_SEAVEYS-1; 
    }

    else if (!strcasecmp(prefix,"allVpol"))
    {
      start_pol = 1; start_ant = 0; 
      stop_pol = 1; stop_ant = NUM_SEAVEYS-1; 
    }
    else
    {
      int phi; 
      char ring; 
      char pol; 

      sscanf(prefix,"%02d%c%c", &phi,&ring,&pol); 

      if (pol == 'h' || pol == 'H')
      {
        start_pol = 0; 
        stop_pol = 0; 
      }
      else if  (pol == 'v' || pol == 'V')
      {
        start_pol = 1; 
        stop_pol = 1; 
      }
      else
      {
        fprintf(stderr,"Something wrong with %s\n",prefix); 
        continue; 
      }

      AnitaRing::AnitaRing_t the_ring; 

      if (ring == 'T' || ring == 't')
      {
        the_ring = AnitaRing::kTopRing; 
      }
      else  if (ring == 'M' || ring == 'm')
      {
        the_ring = AnitaRing::kMiddleRing; 
      }
      else if (ring == 'B' || ring == 'b')
      {
        the_ring = AnitaRing::kBottomRing; 
      }
     
      else
      {
          fprintf(stderr,"Something wrong with %s\n",prefix); 
          continue; 
      }

      int ant = AnitaGeomTool::Instance()->getAntFromPhiRing(phi-1, the_ring); 
      start_ant = ant; 
      stop_ant = ant; 
    }

//    printf("\t%d %d %d %d\n", start_ant, stop_ant, start_pol, stop_pol); 
    for (int pol = start_pol; pol <= stop_pol; pol++)
    {
      for (int ant = start_ant; ant <= stop_ant; ant++) 
      {
        responses[ant][pol] = r; 
      }
    }


    free(prefix); 
  }

  for (int ant = 0; ant < NUM_SEAVEYS; ant++)
  {
    for (int pol = 0; pol <2; pol++) 
    {
      if (! responses[ant][pol])
      {
        printf("%d %d\n",ant,pol); 
      }
      assert(responses[ant][pol]); 
    }
  }


  closedir(dp); 


  return 0; 
}

// AnitaResponse::ResponseManager::ResponseManager(const AnitaResponse::AnalysisConfig *cfg ) 
// {
//   memset(responses,0,sizeof(responses)); 

//   if (cfg->response_option)
//   {
//     loadResponsesFromDir(AnalysisConfig::getResponseString(cfg->response_option),cfg->response_npad); 
//   }

//   method = cfg->deconvolution_method; 

// }

AnitaResponse::ResponseManager::ResponseManager(const char * dir, int npad, const AnitaResponse::DeconvolutionMethod* methodPtr)
{
  memset(responses,0,sizeof(responses)); 
  loadResponsesFromDir(dir,npad);
  method = methodPtr == NULL ? &AnitaResponse::kDefaultDeconvolution : methodPtr;
  // method = &kDefaultDeconvolution;
}

AnitaResponse::ResponseManager::~ResponseManager() 
{

  for (size_t i = 0; i < response_store.size(); i++) 
  {
    delete response_store[i]; 
  }
}
