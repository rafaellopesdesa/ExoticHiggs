#include "fastjet/ClusterSequence.hh"
#include "generator/McEventCollection_p5.h"
#include "generator/GenParticle_p5.h"

#include <vector>

#include "TLorentzVector.h"

int findHardScatterLepton(std::vector<GenParticle_p5>& partList, bool isPythia6 = false);
double LeptonIsolation(fastjet::PseudoJet& lepton, std::vector<fastjet::PseudoJet>& jetList, double Rmin, double Rmax);
double LeptonMiniIsolation(GenParticle_p5& lepton, std::vector<GenParticle_p5>& partList);
