#include "fastjet/ClusterSequence.hh"
#include "generator/McEventCollection_p5.h"
#include "generator/GenParticle_p5.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLorentzVector.h"

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
  for (Long64_t ievt=0; ievt<nentries; ievt++) {

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


    double ptmin = 10000.0;
    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
    vector<TLorentzVector> inclusive_jets_vec;
    for (auto jet : inclusive_jets) {
      TLorentzVector jetvec(jet.px(), jet.py(), jet.pz(), jet.e());
      if (fabs(jetvec.Eta()) < 4.7 && jetvec.Pt() > 20)
	inclusive_jets_vec.push_back(jetvec);
    }

    vector<TLorentzVector> matched_jets;
    for (auto bs : theBS) {
      if (inclusive_jets_vec.size() > 0) {
	auto min = min_element(inclusive_jets_vec.begin(), inclusive_jets_vec.end(),
			       [bs](TLorentzVector v1, TLorentzVector v2)->bool {return bs.DeltaR(v1) < bs.DeltaR(v2);});
	if (bs.DeltaR(*min) < 0.4) {
	  matched_jets.push_back(*min);
	}
      } else break;      
    }

    
    // Write a report
    cout << "Using algorithm " << jet_def.description() << endl;
    cout << "Fastjet found " << inclusive_jets.size() << " jets with pT > 5 GeV" << endl;
    cout << "I found " << inclusive_jets_vec.size() << " jets with pT > 20 GeV and |eta| < 4.7" << endl;
    cout << "Resolved case!" << endl;
    cout << "First a->bb: " << endl;
    cout << "    b parton 1 : " << theBS[0].Eta() << " " << theBS[0].Phi() << " " << theBS[0].Pt()/1000. << endl;
    cout << "    b jet    1 : " << matched_jets[0].Eta() << " " << matched_jets[0].Phi() << " " << matched_jets[0].Pt()/1000. << endl;
    cout << "    b parton 2 : " << theBS[1].Eta() << " " << theBS[1].Phi() << " " << theBS[1].Pt()/1000. << endl;
    cout << "    b jet    2 : " << matched_jets[1].Eta() << " " << matched_jets[1].Phi() << " " << matched_jets[1].Pt()/1000. << endl;
    cout << "    a mass     : " << (theBS[0] + theBS[1]).M()/1000. << endl;
    cout << "    dijet mass : " << (matched_jets[0] + matched_jets[1]).M()/1000. << endl;
    cout << "Second a->bb: " << endl;
    cout << "    b parton 1 : " << theBS[2].Eta() << " " << theBS[2].Phi() << " " << theBS[2].Pt()/1000. << endl;
    cout << "    b jet    1 : " << matched_jets[2].Eta() << " " << matched_jets[2].Phi() << " " << matched_jets[2].Pt()/1000. << endl;
    cout << "    b parton 2 : " << theBS[3].Eta() << " " << theBS[3].Phi() << " " << theBS[3].Pt()/1000. << endl;
    cout << "    b jet    2 : " << matched_jets[3].Eta() << " " << matched_jets[3].Phi() << " " << matched_jets[3].Pt()/1000. << endl;
    cout << "    a mass     : " << (theBS[2] + theBS[3]).M()/1000. << endl;
    cout << "    dijet mass : " << (matched_jets[2] + matched_jets[3]).M()/1000. << endl;

  }
  return 0;

}


    
