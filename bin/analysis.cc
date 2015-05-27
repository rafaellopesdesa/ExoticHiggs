// Analysis strategy

/*


  1) Reconstruct C/A with 1.0
  2) Do some sort of substructure tagging
  3) Reconstruct the subjets
  4) b-tag them
  5) Reconstruct the mass of the two subjets  

 */ 


#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "generator/McEventCollection_p5.h"
#include "generator/GenParticle_p5.h"

#include "ExoticHiggs/LeptonUtils.h"
#include "ExoticHiggs/JetUtils.h"

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <algorithm>
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <cmath>

using namespace std;


int main(int argc, char** argv){

  // Parameters

  // Fat jet definition
  double fat_R = 1.0;
  fastjet::JetDefinition fat_jet_def(fastjet::cambridge_algorithm, fat_R);

  // Fat jet selection
  double jet_ptmin = 40000.0;
  double jet_etamax = 2.6;

  // Substructure tagger
  double mu_thr = 0.667;
  double y_thr = 0.09;    

  // Prunning
  bool doPrunning = true;
  double zcut = 0.1;
  double rcut_factor = 0.5;
  
  // Lepton selection
  double lep_ptmin = 26000.0;
  double lep_etamax = 2.4;
  double lep_isomax = 0.5; // mini iso limit

  // Ghost factor for truth-level B-jet tagging
  double ghost_factor = 1.e-21;

  // Sub jet definition
  double sub_R = 0.4;
  fastjet::JetDefinition sub_jet_def(fastjet::antikt_algorithm, sub_R);
  double sub_jet_ptmin = 5000.0;
  bool use_only_charged = false;

  // Histograms to do
  TH1D* bparton_pt = new TH1D("bparton_pt", "bparton_pt", 100, 0., 150.);
  TH1D* bbar_dR = new TH1D("bbar_dR", "bbar_dR", 100, 0., 5.);
  TH2D* bbar_dR_vs_H_pT = new TH2D("bbar_dR_vs_H_pT", "bbar_dR_vs_H_pT", 100, 0., 5., 100, 0., 150.);
  TH1D* njets = new TH1D("njets", "njets", 10, -0.5, 9.5);
  TH1D* nbjets = new TH1D("nbjets", "nbjets", 10, -0.5, 9.5);
  TH1D* lepton_pt = new TH1D("lepton_pt", "lepton_pt", 100, 0., 150.);
  TH1D* lepton_miniiso = new TH1D("lepton_miniiso", "lepton_miniiso", 100, 0., 0.5);
  TH1D* fatjet_mass = new TH1D("fatjet_mass", "fatjet_mass", 100, 0., 150.);
  TH1D* fatjet_pt = new TH1D("fatjet_pt", "fatjet_pt", 100, 0., 150.);  
  TH1D* fatjet_pruned_mass = new TH1D("fatjet_pruned_mass", "fatjet_pruned_mass", 100, 0., 150.);
  TH1D* fatjet_pruned_pt = new TH1D("fatjet_pruned_pt", "fatjet_pruned_pt", 100, 0., 150.);  
  TH1D* nsubjet = new TH1D("nsubjet", "nsubjet", 10, -0.5, 9.5);
  TH1D* nbsubjet = new TH1D("nbsubjet", "nbsubjet", 10, -0.5, 9.5);  
  TH1D* subjet_mass_mdpieces = new TH1D("subjet_mass_mdpieces", "subjet_mass_mdpieces", 100, 0., 150.);
  TH1D* subjet_mass = new TH1D("subjet_mass", "subjet_mass", 100, 0., 150.);
  TH1D* bsubjet_mass_mdpieces = new TH1D("bsubjet_mass_mdpieces", "bsubjet_mass_mdpieces", 100, 0., 150.);
  TH1D* bsubjet_mass = new TH1D("bsubjet_mass", "bsubjet_mass", 100, 0., 150.);
  
  
  if (argc != 3) {
    cout << "analysis.exe [eventFile] [outputFile]" << endl;
    return 0;
  }

  TDatabasePDG *db= TDatabasePDG::Instance();

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
  
  Long64_t nentries = CollectionTree->GetEntries();
  cout << "Reading " << nentries << " events" << endl;
  for (Long64_t ievt=0; ievt<nentries; ievt++) {
    if (ievt % 500 == 0) cout << "Event  " << ievt << endl;
    
    branch->GetEntry(ievt);
    TLorentzVector partvec;
    input_particles.clear();
    int index = -1;

    // Again kind of an overkill, but ok.
    vector<vector<int>> theBHdecays = findBHdecays(event->m_genParticles);

    int theHL = findHardScatterLepton(event->m_genParticles);
    fastjet::PseudoJet Vlepton_jet;
    TLorentzVector Vlepton_vec;
    GenParticle_p5 Vlepton_part;

    for (auto part : event->m_genParticles) {
      index++;
      partvec.SetXYZM(part.m_px, part.m_py, part.m_pz, part.m_m);

      if (index == theBHdecays[0][0]) {
	partvec *= ghost_factor;
	fastjet::PseudoJet bghost(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	bghost.set_user_index(-50001);
	input_particles.push_back(bghost);
	continue;
      } else if (index == theBHdecays[0][1]) {
	partvec *= ghost_factor;
	fastjet::PseudoJet bghost(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	bghost.set_user_index(-50002);
	input_particles.push_back(bghost);
	continue;
      } else if (index == theBHdecays[1][0]) {
	partvec *= ghost_factor;
	fastjet::PseudoJet bghost(partvec.Px(),partvec.Py(),partvec.Pz(),partvec.E());
	bghost.set_user_index(-50003);
	input_particles.push_back(bghost);
	continue;
      } else if (index == theBHdecays[1][1]) {
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
    fastjet::ClusterSequence clust_seq(input_particles, fat_jet_def);    
    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(jet_ptmin));
    
    fastjet::MassDropTagger md_tagger(mu_thr, y_thr);        
    fastjet::Pruner pruner(fastjet::cambridge_algorithm, zcut, rcut_factor);

    vector<particleJet> selected_jets;
    for (auto jet : inclusive_jets) {

      TLorentzVector jetvec(jet.px(), jet.py(), jet.pz(), jet.e());
      if (fabs(jetvec.Eta()) > jet_etamax) continue;

      particleJet partJet;
      partJet.jet = jetvec;
      partJet.isBjet = isBjet(jet);
      partJet.pseudoJet = jet;

      if (isBpartonJet(jet,1)) {
      	partJet.parton.push_back(&(event->m_genParticles[theBHdecays[0][0]]));
      } else if (isBpartonJet(jet,2)) {
      	partJet.parton.push_back(&(event->m_genParticles[theBHdecays[0][1]]));
      } else if (isBpartonJet(jet,3)) {
      	partJet.parton.push_back(&(event->m_genParticles[theBHdecays[1][0]]));
      } else if (isBpartonJet(jet,4)) {
      	partJet.parton.push_back(&(event->m_genParticles[theBHdecays[1][1]]));
      }

      fastjet::PseudoJet tagged = md_tagger(jet);
      partJet.hasSubstructure = false;
      if (!(tagged == 0)) {
	partJet.hasSubstructure = true;

	input_particles.clear();
	fastjet::PseudoJet pruned_jet;
	if (doPrunning)
	  pruned_jet= pruner(jet);
	else
	  pruned_jet = jet;
	partJet.pruned_pseudoJet = pruned_jet;

	input_particles = pruned_jet.constituents();
	if (use_only_charged) {
	  for (auto ipart = input_particles.begin(); ipart != input_particles.end(); ipart++) {
	    if (ipart->user_index() > 0) {
	      if ((db->GetParticle(event->m_genParticles[ipart->user_index()].m_pdgId))->Charge() == 0)
		input_particles.erase(ipart);
	    }
	  }
	}

	fastjet::ClusterSequence sub_clust_seq(input_particles, sub_jet_def);    
	vector<fastjet::PseudoJet> sub_inclusive_jets = sorted_by_pt(sub_clust_seq.inclusive_jets(sub_jet_ptmin));

	vector<int> sub_jets_btag;

	for (auto subjet : sub_inclusive_jets)
	  sub_jets_btag.push_back((int) isBjet(subjet));

	partJet.subjets = sub_inclusive_jets;
	partJet.subjets_btag = sub_jets_btag;
	partJet.pieces = tagged.pieces();
      }

      selected_jets.push_back(partJet);

    }
      
    vector<TLorentzVector> selected_lepton;
    if (Vlepton_vec.Pt() > lep_ptmin &&
	fabs(Vlepton_vec.Eta()) < lep_etamax &&
	!(abs(Vlepton_part.m_pdgId) == 11 && LeptonMiniIsolation(Vlepton_part, event->m_genParticles) > lep_isomax))
      selected_lepton.push_back(Vlepton_vec);

    // now do the analysis

    // Truth level, before selection

    TLorentzVector b1, b2, b3, b4;
    b1.SetXYZM(event->m_genParticles[theBHdecays[0][0]].m_px, event->m_genParticles[theBHdecays[0][0]].m_py, event->m_genParticles[theBHdecays[0][0]].m_pz, event->m_genParticles[theBHdecays[0][0]].m_m);
    b2.SetXYZM(event->m_genParticles[theBHdecays[0][1]].m_px, event->m_genParticles[theBHdecays[0][1]].m_py, event->m_genParticles[theBHdecays[0][1]].m_pz, event->m_genParticles[theBHdecays[0][1]].m_m);
    b3.SetXYZM(event->m_genParticles[theBHdecays[1][0]].m_px, event->m_genParticles[theBHdecays[1][0]].m_py, event->m_genParticles[theBHdecays[1][0]].m_pz, event->m_genParticles[theBHdecays[1][0]].m_m);
    b4.SetXYZM(event->m_genParticles[theBHdecays[1][1]].m_px, event->m_genParticles[theBHdecays[1][1]].m_py, event->m_genParticles[theBHdecays[1][1]].m_pz, event->m_genParticles[theBHdecays[1][1]].m_m);

    bparton_pt->Fill(b1.Pt()/1000.);
    bparton_pt->Fill(b2.Pt()/1000.);
    bparton_pt->Fill(b3.Pt()/1000.);
    bparton_pt->Fill(b4.Pt()/1000.);

    bbar_dR->Fill(b1.DeltaR(b2));
    bbar_dR->Fill(b3.DeltaR(b4));

    bbar_dR_vs_H_pT->Fill(b1.DeltaR(b2), (b1+b2+b3+b4).Pt()/1000.);
    bbar_dR_vs_H_pT->Fill(b3.DeltaR(b4), (b1+b2+b3+b4).Pt()/1000.);

    lepton_pt->Fill(Vlepton_vec.Pt()/1000.);
    lepton_miniiso->Fill(LeptonMiniIsolation(Vlepton_part, event->m_genParticles));

    // Now do lepton selection
    if (selected_lepton.size() > 0) {
      njets->Fill(selected_jets.size());

      int nb = 0;
      for (auto jet : selected_jets) {       
	if (jet.isBjet) nb++;
	fatjet_mass->Fill(jet.pseudoJet.m()/1000.);
	fatjet_pruned_mass->Fill(jet.pruned_pseudoJet.m()/1000.);
	fatjet_pt->Fill(jet.pseudoJet.pt()/1000.);
	fatjet_pruned_pt->Fill(jet.pruned_pseudoJet.pt()/1000.);

	if (jet.hasSubstructure) {

	  int nbsub = 0;
	  int bsub1 = -1;
	  int bsub2 = -1;
	  nsubjet->Fill(jet.subjets.size());
	  for (int isub = 0; isub < jet.subjets_btag.size(); isub++) {
	    if (jet.subjets_btag[isub]) {
	      nbsub++;
	      bsub1 < 0 ? bsub1 = isub : bsub2 = isub;
	    }
	  }
	  nbsubjet->Fill(nbsub);

	  subjet_mass_mdpieces->Fill((jet.pieces[0]+jet.pieces[1]).m()/1000.);
	  subjet_mass->Fill((jet.subjets[0]+jet.subjets[1]).m()/1000.);

	  if (isBjet(jet.pieces[0]) && isBjet(jet.pieces[1]))
	    bsubjet_mass_mdpieces->Fill((jet.pieces[0]+jet.pieces[1]).m()/1000.);
	  if (nbsub >= 2)
	    bsubjet_mass->Fill((jet.subjets[bsub1]+jet.subjets[bsub2]).m()/1000.);
	}
      }
      nbjets->Fill(nb);
    }
      
    
      
  }
  dataFile->Close();
  TFile* resultFile = TFile::Open(argv[2], "RECREATE");

  bparton_pt->Write();
  bbar_dR->Write();
  bbar_dR_vs_H_pT->Write();
  njets->Write();
  nbjets->Write();
  lepton_pt->Write();
  lepton_miniiso->Write();
  fatjet_mass->Write();
  fatjet_pt->Write();
  fatjet_pruned_mass->Write();
  fatjet_pruned_pt->Write();
  nsubjet->Write();
  nbsubjet->Write();
  subjet_mass_mdpieces->Write();
  subjet_mass->Write();
  bsubjet_mass_mdpieces->Write();
  bsubjet_mass->Write();

  resultFile->Close();
  return 0;

}


    
