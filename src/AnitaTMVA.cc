#include "AnitaTMVA.h" 
#include "TTree.h" 
#include "TFile.h" 
#include <string>
#include "TString.h" 
#include <stdlib.h> 
#include <sstream>



struct idiocy
{
  int i;
  float f;
}; 


AnitaTMVA::MVAVarSet::MVAVarSet(int n, const char * names[], const char * expr[], const char * type, const bool * spectator) 
{
  for (int i = 0; i < n; i++) 
  {
    add(MVAVar(expr[i], names[i], type ? type[i] :'F', spectator ? spectator[i] : false)); 
  }
}

AnitaTMVA::MVAVarSet::MVAVarSet(const char * ifile) 
{

  FILE * f = fopen(ifile,"r"); 

  if (!f) 
  {
    fprintf(stderr,"Could not load %s\n",ifile); 
    return;
  }


  char buf[2048]; 
  int lineno = 0;
  while(!feof(f))
  {
    fgets(buf,sizeof(buf), f); 
    buf[strlen(buf)-1]=0; //truncate \n 
    lineno++;

    char * comment = strchr(buf,'#'); 
    if (comment) *comment = 0; 



    if(!strchr(buf,'@'))
      continue; 

    std::vector<TString> toks; 
    char * tok = strtok(buf,"@"); 
    while (tok)
    {
      toks.push_back(TString(tok).Strip(TString::kBoth,' ')); 
      tok = strtok(0,"@"); 
    }

    if (toks.size() < 2) 
    {
      fprintf(stderr,"Not enough tokens in line %d: %s\n", lineno, buf); 
    }
    
    printf("Adding: %s/%s\n", toks[1].Data(), toks[0].Data()); 
    add(MVAVar(strdup(toks[1].Data()), strdup(toks[0].Data()), toks.size() > 2 ? toks[2].Data()[0] : 'F' , toks.size() > 3 ? atoi(toks[3].Data()) : false)); 
  }


  fclose(f); 

}



TTree* AnitaTMVA::makeTMVATree(TTree * in, TFile * outf, const char * tree_name, const AnitaTMVA::MVAVarSet & vars , const char * cut) 
{
  return makeTMVATree(1,&in,outf,tree_name,vars,cut); 
}

TTree* AnitaTMVA::makeTMVATree(int ntrees, TTree ** in, TFile * outf, const char * tree_name, const AnitaTMVA::MVAVarSet & vars , const char * cut) 
{

  std::stringstream drawstr; 

  outf->cd(); 
  TTree * out = new TTree(tree_name,tree_name); 


  std::vector<idiocy> mem(vars.N()+2); 

  for (int i = 0; i < vars.N(); i++) 
  {
    drawstr << vars.at(i).expression << ":"; 

    char type = vars.at(i).type; 
    switch(type)
    {
      case 'I': 
        out->Branch(vars.at(i).name, &mem[i].i); break; 
      case 'F': 
      default:
        out->Branch(vars.at(i).name, &mem[i].f);
    }
  }

  out->Branch("entry", &mem[vars.N()].i); 
  out->Branch("iteration", &mem[vars.N()+1].i); 
  drawstr << "Entry$:Iteration$"; 

  std::vector<int> Nout(ntrees); 
//  printf("%s\n",drawstr.str().c_str()); 

  //now the real work happens... which we can parallelize! 
#ifdef USE_OMP
#pragma omp parallel for 
#endif
  for (int t = 0; t < ntrees; t++) 
  {
    in[t]->SetEstimate(in[t]->GetEntries() *10); 
    Nout[t] = in[t]->Draw(drawstr.str().c_str(),cut,"goff"); 
  }
 

  outf->cd(); 
  //This is madness. I should just use a custom selector but I'm too lazy for that right now. 
  for (int t = 0; t < ntrees; t++)
  {
    printf("Nout[%d]=%d\n", t, Nout[t]); 
    for (int j =0; j < Nout[t]; j++) 
    {

      for (int i = 0; i < vars.N()+2; i++)
      {
        char type = i < vars.N() ? vars.at(i).type : 'I'; 
        switch(type) 
        {
          case 'I':
            mem[i].i = in[t]->GetVal(i)[j]; break;
          case 'F': 
          default: 
            mem[i].f = in[t]->GetVal(i)[j]; 
        }
      }

      outf->cd(); 
      out->Fill(); 
    }
  }

  return out; 
}


int AnitaTMVA::evaluateTMVA(TTree * tree, const AnitaTMVA::MVAVarSet & vars, const char * branch_name, const char * weights_file, double aux) 
{
  std::vector<idiocy> mem(vars.N()); 
  TMVA::Reader reader; 


  for (int i = 0; i < vars.N(); i++) 
  {
    if(! tree->FindBranch(vars.at(i).name))
    {
      fprintf(stderr,"%s not found in tree... did you create this tree with the same variable set? Aborting evaluateTMVA.", vars.at(i).name); 
      return 1; 
    }

    if (vars.at(i).spectator) continue; 

    switch(vars.at(i).type)
    {
      case 'I':
         tree->SetBranchAddress(vars.at(i).name, &mem[i].i); 
         reader.AddVariable(vars.at(i).name, &mem[i].f); //ugh ??!? 
      case 'F':
      default:
         tree->SetBranchAddress(vars.at(i).name, &mem[i].f); 
         reader.AddVariable(vars.at(i).name, &mem[i].f); 
         break;
    }
  }

  float value = 0; 
  reader.BookMVA(branch_name,weights_file); 
  TBranch * b = tree->Branch( branch_name ,&value); 

  for (int j = 0; j < tree->GetEntries(); j++) 
  {
    tree->GetEntry(j); 
    if (vars.at(j).type == 'I') 
    {
      mem[j].f = mem[j].i; 
    }
    value = reader.EvaluateMVA(branch_name,aux); 
    b->Fill(); 
  }



  return 0; 

}



int AnitaTMVA::MVAVarSet::toFile(const char * ofile) 
{

  FILE * f = fopen(ofile,"w"); 
  if (!f) return 1; 

  for (int i = 0; i < N(); i++)
  {
    fprintf(f,"%s @ %s @ %c @ %d\n", at(i).name, at(i).expression, at(i).type, at(i).spectator); 
  }
  fclose(f); 
  return 0; 
}





