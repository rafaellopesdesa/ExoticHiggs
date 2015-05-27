#include "fastjet/ClusterSequence.hh"
#include "generator/McEventCollection_p5.h"
#include "generator/GenParticle_p5.h"

#include "ExoticHiggs/LeptonUtils.h"
#include "ExoticHiggs/JetUtils.h"

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

  double ghost_factor = 10.e-18;
  
  Long64_t nentries = CollectionTree->GetEntries();
  cout << "Reading " << nentries << " events" << endl;
  for (Long64_t ievt=0; ievt<nentries; ievt++) {
    if (ievt % 500 == 0) cout << "Event  " << ievt << endl;
    
    branch->GetEntry(ievt);
    TLorentzVector partvec;
    input_particles.clear();
    int index = -1;

    // Again kind of an overkill, but ok.
    vector<TLorentzVector> theBS = findBS(event->m_genParticles);
    vector<vector<int>> theBHdecays = findBHdecays(event->m_genParticles);

    int theHL = findHardScatterLepton(event->m_genParticles);
    fastjet::PseudoJet Vlepton_jet;
    TLorentzVector Vlepton_vec;
    GenParticle_p5 Vlepton_part;

    for (auto part : event->m_genParticles) {
      index++;
      partvec.SetXYZM(part.m_px, part.m_py, part.m_pz, part.m_m);

      if (index == theBHdecays[0][0]) {
	cout << "add 1" << endl;
	partvec *= ghost_factor;
	fastjet::PseudoJet bghost(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	bghost.set_user_index(-50001);
	input_particles.push_back(bghost);
	continue;
      } else if (index == theBHdecays[0][1]) {
	cout << "add 2" << endl;
	partvec *= ghost_factor;
	fastjet::PseudoJet bghost(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	bghost.set_user_index(-50002);
	input_particles.push_back(bghost);
	continue;
      } else if (index == theBHdecays[1][0]) {
	cout << "add 3" << endl;
	partvec *= ghost_factor;
	fastjet::PseudoJet bghost(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	bghost.set_user_index(-50003);
	input_particles.push_back(bghost);
	continue;
      } else if (index == theBHdecays[1][1]) {
	cout << "add 4" << endl;
	partvec *= ghost_factor;
	fastjet::PseudoJet bghost(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	bghost.set_user_index(-50004);
	input_particles.push_back(bghost);
	continue;
      }

      if (isBhadron(part.m_pdgId)) {
	partvec *= ghost_factor;
	fastjet::PseudoJet bghost(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	bghost.set_user_index(-5);
	input_particles.push_back(bghost);
	continue;
      }

      // Stable particles
      if (part.m_status != 1) continue;

      // Not neutrinos
      if (abs(part.m_pdgId) == 12) continue;
      if (abs(part.m_pdgId) == 14) continue;
      if (abs(part.m_pdgId) == 16) continue;

      if (index == theHL) {
	// A bit of an overkill, but simpler.
	Vlepton_jet = fastjet::PseudoJet(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	Vlepton_vec = partvec;
	Vlepton_part = part;
	continue;
      }

      
      input_particles.push_back(fastjet::PseudoJet(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E()));
      input_particles.back().set_user_index(index);
    }

    // cluster the jets
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);

    // Jet selection
    double jet_ptmin = 20000.0;
    double jet_etamax = 2.6;

    // Lepton selection
    double lep_ptmin = 26000.0;
    double lep_etamax = 2.4;
    double lep_isomax = 0.15;
    
    // and look for subjets
    double Rsub = 0.5;
    double dcut = pow(Rsub/R,2);
    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(jet_ptmin));

    vector<particleJet> selected_jets;
    for (auto jet : inclusive_jets) {
      particleJet jetvec(TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.e()), isBjet(jet));
      if (isBpartonJet(jet,1)) {
	cout << "found 1" << endl;
	jetvec.parton.push_back(&(event->m_genParticles[theBHdecays[0][0]]));
      } else if (isBpartonJet(jet,2)) {
	cout << "found 2" << endl;
	jetvec.parton.push_back(&(event->m_genParticles[theBHdecays[0][1]]));
      }	else if (isBpartonJet(jet,3)) {
	cout << "found 3" << endl;
	jetvec.parton.push_back(&(event->m_genParticles[theBHdecays[1][0]]));
      }	else if (isBpartonJet(jet,4)) {
	cout << "found 4" << endl;
	jetvec.parton.push_back(&(event->m_genParticles[theBHdecays[1][1]]));
      }	
      if (fabs(jetvec.jet.Eta()) < jet_etamax)
	selected_jets.push_back(jetvec);
    }

    vector<TLorentzVector> selected_lepton;
    if (Vlepton_vec.Pt() > lep_ptmin &&
	fabs(Vlepton_vec.Eta()) < lep_etamax &&
	!(abs(Vlepton_part.m_pdgId) == 11 && LeptonIsolation(Vlepton_jet, inclusive_jets, 0.2, 0.4) > lep_isomax))
      selected_lepton.push_back(Vlepton_vec);
    cout << "--- Event report ---" << endl;
    cout << "b quarks" << endl;
    cout << "   " << event->m_genParticles[theBHdecays[0][0]] << endl;
    cout << "   " << event->m_genParticles[theBHdecays[0][1]] << endl;
    cout << "   " << event->m_genParticles[theBHdecays[1][0]] << endl;
    cout << "   " << event->m_genParticles[theBHdecays[1][1]] << endl;
    cout << "I found " << selected_jets.size() << "/" << selected_lepton.size() << " good jets/leptons" << endl;
    cout << "Jets: " << endl;
    for (auto jet : selected_jets)
      cout << jet;
    
    // vector<TLorentzVector> matched_jets;
    // if (selected_jets.size() > 0) {
    //   for (auto bs : theBS) {
    // 	auto min = min_element(selected_jets.begin(), selected_jets.end(),
    // 			       [bs](TLorentzVector v1, TLorentzVector v2)->bool {return bs.DeltaR(v1) < bs.DeltaR(v2);});
    // 	matched_jets.push_back(*min);
    //   }
    // }

      
  }
  dataFile->Close();
  TFile* resultFile = TFile::Open(argv[2], "RECREATE");
  nsubjets->Write();
  resultFile->Close();
  return 0;

}


    
