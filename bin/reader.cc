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
  
  branch->SetAddress(&event);

  std::cout << branch->GetEntry(102) << std::endl;
  std::cout << branch->GetEntry(103) << std::endl;
  std::cout << event->m_genParticles.size() << std::endl;
  std::cout << event->m_genParticles[100].m_m << std::endl;
  //  delete event;

  
}
  
