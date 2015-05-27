#include "ExoticHiggs/LeptonUtils.h"

#include <cmath>

int findHardScatterLepton(std::vector<GenParticle_p5>& partList)
{

  int index = -1;
  int prodVtx = 1;
  for (auto part : partList) {
    index++;
    if ((abs(part.m_pdgId) == 24 || abs(part.m_pdgId) == 23) && part.m_status == 22) {
      prodVtx = part.m_endVtx;
    }
    if ((abs(part.m_pdgId) == 24 || abs(part.m_pdgId) == 23) && part.m_prodVtx == prodVtx) {
      prodVtx = part.m_endVtx;
    }
    if ((abs(part.m_pdgId) == 11 || abs(part.m_pdgId) == 13) && part.m_prodVtx == prodVtx) {
      if (part.m_status == 1) {
	break;
      } else {
	prodVtx = part.m_endVtx;
      }
    }
  }
  return index;
}


double LeptonIsolation(fastjet::PseudoJet& lepton, std::vector<fastjet::PseudoJet>& jetList, double Rmin, double Rmax)
{
  double iso = 0.;
  for (auto jet : jetList) {
    if (jet.delta_R(lepton) < Rmax && jet.delta_R(lepton) > Rmin)
      iso += jet.e();
  }
  iso /= lepton.e();
  return iso;
}
    
    

