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

#include "inih/INIReader.h"

#include "ExoticHiggs/LeptonUtils.h"
#include "ExoticHiggs/JetUtils.h"

#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TString.h"
#include "TFileCollection.h"

#include <algorithm>
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <cmath>
#include <fstream>

using namespace std;


int main(int argc, char** argv){

  // Parameters
  if (argc != 2) {
    cout << "analysis.exe [configFile]" << endl;
    return 0;
  }
  INIReader reader(argv[1]);
  if (reader.ParseError() < 0) {
    std::cout << "Can't load " << argv[1] << endl;
    return 1;
  }
  
  // Fat jet definition
  double fat_R = reader.GetInteger("fatjet", "R", 1.0);
  fastjet::JetDefinition fat_jet_def(fastjet::cambridge_algorithm, fat_R);

  // Fat jet selection
  double jet_ptmin = reader.GetReal("fatjet", "ptmin", 40000.0);
  double jet_etamax = reader.GetReal("fatjet", "etamax", 2.6);

  // Substructure tagger
  double mu_thr = reader.GetReal("mdtagger", "mu_thr", 0.667);
  double y_thr = reader.GetReal("mdtagger", "y_thr", 0.09);    

  // Prunning
  bool doPrunning = reader.GetBoolean("prunning", "doPrunning", true);
  double zcut = reader.GetReal("prunning", "zcut", 0.1);
  double rcut_factor = reader.GetReal("prunning", "rcut_factor", 0.5);
  
  // Lepton selection
  double lep_ptmin = reader.GetReal("lepton", "lep_ptmin", 26000.0);
  double lep_etamax = reader.GetReal("lepton", "lep_etamax", 2.4);
  double lep_isomax = reader.GetReal("lepton", "lep_isomax", 0.5); // mini iso limit

  // Ghost factor for truth-level B-jet tagging
  double ghost_factor = reader.GetReal("btagging", "ghost_factor", 1.e-21);

  // Sub jet definition
  double sub_R = reader.GetReal("subjet", "R", 0.4);
  fastjet::JetDefinition sub_jet_def(fastjet::antikt_algorithm, sub_R);
  double sub_jet_ptmin = reader.GetReal("subjet", "ptmin", 5000.0);
  bool use_only_charged = reader.GetBoolean("subjet", "use_only_charged", false);

  // files
  const char* dataFileName = reader.Get("io", "data_file_name", "UNKNOWN").c_str();
  const char* outputFileName = reader.Get("io", "result_file_name", "UNKNOWN").c_str();
  double hypMass = reader.GetReal("io", "hyp_mass", 20.);
  bool isSignal = reader.GetBoolean("io", "isSignal", false);
  bool isPythia6 = reader.GetBoolean("io", "isPythia6", false);

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
  TH1D* cutflow_h = new TH1D("cutflow", "cutflow", 9, 0., 9.);
  
  TDatabasePDG *db= TDatabasePDG::Instance();

  // Gets the input and stores in fastjet
  std::vector<fastjet::PseudoJet> input_particles;
  TChain* CollectionTree = new TChain("CollectionTree");
  TFileCollection fc("dum","",dataFileName);
  CollectionTree->AddFileInfoList((TCollection*) fc.GetList());
  //  TBranch* branch = CollectionTree->GetBranch("McEventCollection_p5_GEN_EVENT");
  McEventCollection_p5* event = 0;

  // I don't care about the other information
  CollectionTree->SetBranchStatus("*", 0);
  CollectionTree->SetBranchStatus("McEventCollection_p5_GEN_EVENT", 1);
  CollectionTree->SetBranchAddress("McEventCollection_p5_GEN_EVENT", &event);
  
  Long64_t nentries = CollectionTree->GetEntries();
  cout << "Reading " << nentries << " events" << endl;


  double cutflow[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
  cutflow_h->GetXaxis()->SetBinLabel(1, "all");
  cutflow_h->GetXaxis()->SetBinLabel(2, "1 lep");  
  cutflow_h->GetXaxis()->SetBinLabel(3, "lep pt");
  cutflow_h->GetXaxis()->SetBinLabel(4, "lep eta");
  cutflow_h->GetXaxis()->SetBinLabel(5, "lep iso");
  cutflow_h->GetXaxis()->SetBinLabel(6, "2 jets");
  cutflow_h->GetXaxis()->SetBinLabel(7, "substructure");
  cutflow_h->GetXaxis()->SetBinLabel(8, "2 subjets");
  cutflow_h->GetXaxis()->SetBinLabel(9, "2 b subjets");
  
  // all, lepton pt, lepton eta, lepton iso, 2 jets, 2 jets with substructure, 2 jet with 2 jets, 2 jets with 2 bjets
  for (Long64_t ievt=0; ievt<nentries; ievt++) {

    // all of it
    cutflow[0]++;

    if (ievt % 500 == 0) cout << "Event  " << ievt << endl;
    
    CollectionTree->GetEntry(ievt);
    TLorentzVector partvec;
    input_particles.clear();
    int index = -1;

    // Again kind of an overkill, but ok.
    vector<vector<int>> theBHdecays = findBHdecays(event->m_genParticles);

    int theHL = findHardScatterLepton(event->m_genParticles, isPythia6);
    fastjet::PseudoJet Vlepton_jet;
    TLorentzVector Vlepton_vec;
    GenParticle_p5 Vlepton_part;

    for (auto part : event->m_genParticles) {
      index++;
      partvec.SetXYZM(part.m_px, part.m_py, part.m_pz, part.m_m);

      if (isSignal) {
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

      if (isSignal) {
	if (isBpartonJet(jet,1)) {
	  partJet.parton.push_back(&(event->m_genParticles[theBHdecays[0][0]]));
	} else if (isBpartonJet(jet,2)) {
	  partJet.parton.push_back(&(event->m_genParticles[theBHdecays[0][1]]));
	} else if (isBpartonJet(jet,3)) {
	  partJet.parton.push_back(&(event->m_genParticles[theBHdecays[1][0]]));
	} else if (isBpartonJet(jet,4)) {
	  partJet.parton.push_back(&(event->m_genParticles[theBHdecays[1][1]]));
	}
      }
      
      fastjet::PseudoJet tagged = md_tagger(jet);
      partJet.hasSubstructure = false;
      if (!(tagged == 0)) {
	partJet.hasSubstructure = true;

	std::vector<fastjet::PseudoJet> charged_input_particles;
	fastjet::PseudoJet pruned_jet;
	if (doPrunning)
	  pruned_jet= pruner(jet);
	else
	  pruned_jet = jet;
	partJet.pruned_pseudoJet = pruned_jet;

	if (use_only_charged) {
	  for (auto& ipart : pruned_jet.constituents()) {
	    if (ipart.user_index() > 0) { 
	      if ((db->GetParticle(event->m_genParticles[ipart.user_index()].m_pdgId))->Charge() != 0)
		charged_input_particles.push_back(ipart);
	    } else {
	      charged_input_particles.push_back(ipart); // keep the ghosts
	    }
	  }
	} else {
	  charged_input_particles = pruned_jet.constituents();
	}
	//	cout << "RCLSA " << charged_input_particles.size() << endl;
	fastjet::ClusterSequence sub_clust_seq(charged_input_particles, sub_jet_def);    
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
    if (theHL > -1) {
      cutflow[1]++;
      if (Vlepton_vec.Pt() > lep_ptmin) {
	cutflow[2]++;
	if (fabs(Vlepton_vec.Eta()) < lep_etamax) {
	  cutflow[3]++;
	  if (LeptonMiniIsolation(Vlepton_part, event->m_genParticles) > lep_isomax) {
	    cutflow[4]++;
	    selected_lepton.push_back(Vlepton_vec);
	  }
	}
      }
    }
    
    // Truth level, before selection
    if (isSignal) {
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
    }

    if (theHL > -1) {
      lepton_pt->Fill(Vlepton_vec.Pt()/1000.);
      lepton_miniiso->Fill(LeptonMiniIsolation(Vlepton_part, event->m_genParticles));
    }
    
    // Now do lepton selection
    if (selected_lepton.size() > 0) {
      njets->Fill(selected_jets.size());

      if (selected_jets.size() >= 2) cutflow[5]++;

      vector<int> cf_jetHasSubstructure;
      vector<int> cf_numberSubjets;
      vector<int> cf_numberBsubjets;
      for (int ijet=0; ijet < selected_jets.size(); ijet++) {
	cf_jetHasSubstructure.push_back(0);
	cf_numberSubjets.push_back(0);
	cf_numberBsubjets.push_back(0);
      }
      int ijet = -1;
      
      int nb = 0;
      for (auto jet : selected_jets) {
	ijet++;
	if (jet.isBjet) nb++;
	fatjet_mass->Fill(jet.pseudoJet.m()/1000.);
	fatjet_pruned_mass->Fill(jet.pruned_pseudoJet.m()/1000.);
	fatjet_pt->Fill(jet.pseudoJet.pt()/1000.);
	fatjet_pruned_pt->Fill(jet.pruned_pseudoJet.pt()/1000.);

	if (jet.hasSubstructure) {

	  cf_jetHasSubstructure[ijet] = 1;
	  cf_numberSubjets[ijet] = jet.subjets.size();

	  int nbsub = 0;
	  int bsub1 = -1;
	  int bsub2 = -1;
	  nsubjet->Fill(jet.subjets.size());
	  for (int isub = 0; isub < jet.subjets_btag.size(); isub++) {
	    if (jet.subjets_btag[isub]) {
	      cf_numberBsubjets[ijet]++;
	      nbsub++;
	      if (bsub1 < 0) {
		bsub1 = isub;
	      } else if (bsub2 < 0) {
		bsub2 = isub;
	      }
	    }
	  }
	  nbsubjet->Fill(nbsub);
	  
	  subjet_mass_mdpieces->Fill((jet.pieces[0]+jet.pieces[1]).m()/1000.);

	  if (jet.subjets.size() >= 2)
	    subjet_mass->Fill((jet.subjets[0]+jet.subjets[1]).m()/1000.);

	  if (isBjet(jet.pieces[0]) && isBjet(jet.pieces[1]))
	    bsubjet_mass_mdpieces->Fill((jet.pieces[0]+jet.pieces[1]).m()/1000.);
	  if (nbsub >= 2) {
	    bsubjet_mass->Fill((jet.subjets[bsub1]+jet.subjets[bsub2]).m()/1000.);
	  }
	}
      }
      nbjets->Fill(nb);

      if (count_if(cf_jetHasSubstructure.begin(), cf_jetHasSubstructure.end(), [&](int x){return (x==1);})>=2)
	cutflow[6]++;
      if (count_if(cf_numberSubjets.begin(), cf_numberSubjets.end(), [&](int x){return (x>=2);})>=2)
	cutflow[7]++;
      if (count_if(cf_numberBsubjets.begin(), cf_numberBsubjets.end(), [&](int x){return (x>=2);})>=2)
	cutflow[8]++;
      
    }         
      
  }

  for (int ibin=1; ibin<=9; ibin++) cutflow_h->SetBinContent(ibin,cutflow[ibin-1]);
  if (isSignal) {
    fatjet_mass->Fit("gaus", "", "", hypMass-15., hypMass+5.);
    fatjet_pruned_mass->Fit("gaus", "", "", hypMass-15., hypMass+5.);
    subjet_mass->Fit("gaus", "", "", hypMass-15., hypMass+5.);
    subjet_mass_mdpieces->Fit("gaus", "", "", hypMass-15., hypMass+5.);
    bsubjet_mass->Fit("gaus", "", "", hypMass-15., hypMass+5.);
    bsubjet_mass_mdpieces->Fit("gaus", "", "", hypMass-15., hypMass+5.);
  }
  
  cout << "---- REPORT ----" << endl;
  cout << "Acceptance: " << endl << "   " << 100*cutflow[8]/cutflow[0] << "%" << endl;

  if (isSignal) {
    cout << "Mass resolution: " << endl;
    cout << "   fat jet mass: " << 100.*((TF1*) fatjet_mass->GetListOfFunctions()->At(0))->GetParameter(2)/hypMass << "%" << endl;
    cout << "   pruned fat jet mass: " << 100.*((TF1*) fatjet_pruned_mass->GetListOfFunctions()->At(0))->GetParameter(2)/hypMass << "%" << endl;
    cout << "   AK04 subjet mass: " << 100.*((TF1*) subjet_mass->GetListOfFunctions()->At(0))->GetParameter(2)/hypMass << "%" << endl;
    cout << "   MD subjet mass: " << 100.*((TF1*) subjet_mass_mdpieces->GetListOfFunctions()->At(0))->GetParameter(2)/hypMass << "%" << endl;
    cout << "   btagged AK04 subjet mass: " << 100.*((TF1*) bsubjet_mass->GetListOfFunctions()->At(0))->GetParameter(2)/hypMass << "%" << endl;
    cout << "   btagged MD subjet mass: " << 100.*((TF1*) bsubjet_mass_mdpieces->GetListOfFunctions()->At(0))->GetParameter(2)/hypMass << "%" << endl;
  }
  
  TFile* resultFile = TFile::Open(outputFileName, "RECREATE");

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
  cutflow_h->Write();
  
  resultFile->Close();
  return 0;

}


    
