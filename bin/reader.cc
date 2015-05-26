#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "generator/McEventCollection_p5.h"
#include "generator/GenParticle_p5.h"

#include <vector>
#include <iostream>

int main() {

  TFile* f = TFile::Open("data/a20a20.root");
  TTree* CollectionTree = (TTree*) f->Get("CollectionTree");
  TBranch* branch = CollectionTree->GetBranch("McEventCollection_p5_GEN_EVENT");
  McEventCollection_p5* event = 0;


  CollectionTree->SetBranchAddress("McEventCollection_p5_GEN_EVENT", &event);
  std::cout << branch->GetEntry(100) << std::endl;
  std::cout << event->m_genParticles.size() << std::endl;

  
}
  
