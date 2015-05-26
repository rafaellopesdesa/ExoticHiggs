#include "fastjet/ClusterSequence.hh"
#include "generator/McEventCollection_p5.h"
#include "generator/GenParticle_p5.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <cmath>

using namespace std;

int main(int argc, char** argv){
  
  // Gets the input and stores in fastjet
  std::vector<fastjet::PseudoJet> input_particles;
  TFile* dataFile = TFile::Open("data/a20a20.root");
  TTree* CollectionTree = (TTree*) dataFile->Get("CollectionTree");
  TBranch* branch = CollectionTree->GetBranch("McEventCollection_p5_GEN_EVENT");
  McEventCollection_p5* event = 0;

  // I don't care about the other information
  CollectionTree->SetBranchStatus("*", 0);
  CollectionTree->SetBranchStatus("McEventCollection_p5_GEN_EVENT", 1);
  CollectionTree->SetBranchAddress("McEventCollection_p5_GEN_EVENT", &event);

  // Which jet we want
  double R = 0.4;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  Long64_t nentries = CollectionTree->GetEntries();
  cout << "Reading " << nentries << " events" << endl;
  for (Long64_t i=0; i<nentries; i++) {

    branch->GetEntry(0);
    TLorentzVector partvec;
    input_particles.clear();
    int index = -1;
    for (auto part : event->m_genParticles) {
      index++;
      partvec.SetXYZM(part.m_px, part.m_py, part.m_pz, part.m_m);
      
      // Stable particles
      if (part.m_status != 1) continue;

      // Not neutrinos
      if (abs(part.m_pdgId) == 12) continue;
      if (abs(part.m_pdgId) == 14) continue;
      if (abs(part.m_pdgId) == 16) continue;

      bool isHardScatterLepton = false;
      if (abs(part.m_pdgId) == 11 || abs(part.m_pdgId) == 13) {
	TLorentzVector partvec2;
	for (auto part2 : event->m_genParticles) {
	  if (part2.m_status != 23 || part.m_pdgId != part2.m_pdgId) continue;
	  partvec2.SetXYZM(part2.m_px, part2.m_py, part2.m_pz, part2.m_m);
	  if (partvec.DeltaR(partvec2) < 0.1) isHardScatterLepton = true;
	}
      }
      if (isHardScatterLepton) continue;

      input_particles.push_back(fastjet::PseudoJet(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E()));
      input_particles.back().set_user_index(index);
    }

    // cluster the jets
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);


    double ptmin = 10000.0;
    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

    cout << "Ran " << jet_def.description() << endl;

    printf("%5s %15s %15s %15s\n","jet #", "eta", "phi", "pt");
 
    for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
      printf("%5u %15.8f %15.8f %15.8f\n",
	     i, inclusive_jets[i].eta(), inclusive_jets[i].phi(),
	     inclusive_jets[i].pt()/1000.);
    }
    
  }
  return 0;

}
