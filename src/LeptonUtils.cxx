#include "ExoticHiggs/LeptonUtils.h"

#include <cmath>
#include <iostream>


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

double LeptonMiniIsolation(GenParticle_p5& lepton, std::vector<GenParticle_p5>& partList)
{

  TLorentzVector lep_vec;
  lep_vec.SetXYZM(lepton.m_px,lepton.m_py,lepton.m_pz,lepton.m_m);
  double dr = 0.2;
  if (lep_vec.Pt() > 50)  dr = 10./lep_vec.Pt();
  if (lep_vec.Pt() > 200) dr = 0.05;

  double iso=0.;
  for (auto part : partList) {
    if (abs(part.m_pdgId) == 12) continue;
    if (abs(part.m_pdgId) == 14) continue;
    if (abs(part.m_pdgId) == 16) continue;

    TLorentzVector part_vec;
    part_vec.SetXYZM(part.m_px,part.m_py,part.m_pz,part.m_m);
    if (part_vec == lep_vec) continue;
    if (part_vec.Pt() < 500.) continue;
    if (lep_vec.DeltaR(part_vec) < dr)
      iso += part_vec.Pt();
  }
  iso /= lep_vec.Pt();
  return iso;
  
}
    

