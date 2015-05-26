#include "fastjet/ClusterSequence.hh"
#include "generator/McEventCollection_p5.h"
#include "generator/GenParticle_p5.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH2D.h"

#include <algorithm>
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <cmath>

using namespace std;

vector<TLorentzVector> findBS(vector<GenParticle_p5>& partList) {

  vector<TLorentzVector> retval;
  TLorentzVector partvec;
  for (auto part : partList) {
    if (part.m_status != 23 || part.m_pdgId != 5) continue;
    partvec.SetXYZM(part.m_px, part.m_py, part.m_pz, part.m_m);
    retval.push_back(partvec);
    for (auto part2 : partList) {
      if (part2.m_status != 23 || part2.m_pdgId != -5) continue;
      if (part.m_prodVtx != part2.m_prodVtx) continue;
      partvec.SetXYZM(part2.m_px, part2.m_py, part2.m_pz, part2.m_m);
      retval.push_back(partvec);
    }
  }
  return retval;

}

int main(int argc, char** argv){

  if (argc != 3) {
    cout << "analysis.exe [eventFile] [outputFile]" << endl;
    return 0;
  }
  
  TH2D* nsubjets = new TH2D("nsubjets", "", 10, -0.5, 9.5, 10, -0.5, 9.5);
  
  // Gets the input and stores in fastjet
  std::vector<fastjet::PseudoJet> input_particles;
  TFile* dataFile = TFile::Open(argv[1]);
  TTree* CollectionTree = (TTree*) dataFile->Get("CollectionTree");
  TBranch* branch = CollectionTree->GetBranch("McEventCollection_p5_GEN_EVENT");
  McEventCollection_p5* event = 0;

  // I don't care about the other information
  CollectionTree->SetBranchStatus("*", 0);
  CollectionTree->SetBranchStatus("McEventCollection_p5_GEN_EVENT", 1);
  CollectionTree->SetBranchAddress("McEventCollection_p5_GEN_EVENT", &event);

  // Which jet we want
  double R = 1.0;
  fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, R);

  Long64_t nentries = CollectionTree->GetEntries();
  cout << "Reading " << nentries << " events" << endl;
  for (Long64_t ievt=0; ievt<nentries; ievt++) {
    if (ievt % 500 == 0) cout << "Event  " << ievt << endl;
    
    branch->GetEntry(ievt);
    TLorentzVector partvec;
    input_particles.clear();
    int index = -1;

    vector<TLorentzVector> theBS = findBS(event->m_genParticles);
    
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

    // and look for subjets
    double Rsub = 0.5;
    double dcut = pow(Rsub/R,2);
    double ptmin = 20000.0;
    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

    vector<TLorentzVector> inclusive_jets_vec;
    for (auto jet : inclusive_jets) {
      TLorentzVector jetvec(jet.px(), jet.py(), jet.pz(), jet.e());
      inclusive_jets_vec.push_back(jetvec);
      vector<fastjet::PseudoJet> subjets = sorted_by_pt(jet.exclusive_subjets(dcut));
      nsubjets->Fill(inclusive_jets.size(),subjets.size());
    }
    vector<TLorentzVector> matched_jets;
    if (inclusive_jets_vec.size() > 0) {
      for (auto bs : theBS) {
	auto min = min_element(inclusive_jets_vec.begin(), inclusive_jets_vec.end(),
			       [bs](TLorentzVector v1, TLorentzVector v2)->bool {return bs.DeltaR(v1) < bs.DeltaR(v2);});
	matched_jets.push_back(*min);
      }
    }

      
  }
  dataFile->Close();
  TFile* resultFile = TFile::Open(argv[2], "RECREATE");
  nsubjets->Write();
  resultFile->Close();
  return 0;

}


    
