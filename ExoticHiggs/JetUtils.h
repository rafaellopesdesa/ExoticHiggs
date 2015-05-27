#include "fastjet/ClusterSequence.hh"
#include "generator/McEventCollection_p5.h"
#include "generator/GenParticle_p5.h"

#include <vector>
#include <utility>
#include <iostream>

#include "TLorentzVector.h"

std::vector<int> findPartonHadrons(int parton, std::vector<GenParticle_p5>& partList);
std::vector< std::vector<int> > findBHdecays(std::vector<GenParticle_p5>& partList);
bool isBhadron(int pdgId);
bool isBjet(const fastjet::PseudoJet& jet);
bool isBpartonJet(const fastjet::PseudoJet& jet, int parton);

class particleJet {

 public:

  particleJet();

  TLorentzVector jet;
  fastjet::PseudoJet pseudoJet;  
  bool isBjet;
  
  // Subjet information
  bool hasSubstructure;
  fastjet::PseudoJet pruned_pseudoJet;
  std::vector<fastjet::PseudoJet> pieces;
  std::vector<fastjet::PseudoJet> subjets;
  std::vector<int> subjets_btag;
  
  std::vector<GenParticle_p5*> parton;
  
};

// bool isBJet(std::vector<GenParticle_p5>& partList, fastjet::PseudoJet& jet);
// bool isBHdecayJet(int BHdecayidx, std::vector<GenParticle_p5>& partList, fastjet::PseudoJet& jet);

std::ostream& operator<<(std::ostream& os, const particleJet jet);
std::ostream& operator<<(std::ostream& os, const GenParticle_p5 parton);
std::ostream& operator<<(std::ostream& ostr, const fastjet::PseudoJet& jet);


