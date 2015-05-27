#include "ExoticHiggs/JetUtils.h"

#include <algorithm>
#include <iostream>

using namespace std;

particleJet::particleJet()
{}

particleJet::particleJet(TLorentzVector _jet, bool _isBjet) :
  jet(_jet),
  isBjet(_isBjet)
{}
  

std::ostream& operator<<(std::ostream& os, const particleJet jet)
{

  os << jet.jet.Px() << " " << jet.jet.Py() << " " << jet.jet.Pz() << " " << jet.jet.E() << " " << jet.isBjet << " " << jet.jet.Eta() << " " << jet.jet.Phi() << " " << jet.jet.Pt()/1000. << endl;
  for (auto parton : jet.parton) {
    TLorentzVector _temp;
    _temp.SetXYZM(parton->m_px,parton->m_py,parton->m_pz,parton->m_m);
    os << "   matched parton: " << parton->m_px << " " << parton->m_py << " " << parton->m_pz << " " << parton->m_pdgId << " " << _temp.Eta() << " " << _temp.Phi() << " " << _temp.Pt()/1000. << endl;
  }
  return os;      
}

std::ostream& operator<<(std::ostream& os, const GenParticle_p5 parton)
{

  TLorentzVector _temp;
  _temp.SetXYZM(parton.m_px,parton.m_py,parton.m_pz,parton.m_m);
  os << parton.m_px << " " << parton.m_py << " " << parton.m_pz << " " << parton.m_pdgId << " " << _temp.Eta() << " " << _temp.Phi() << " " << _temp.Pt()/1000.;
  return os;      
}

bool isBhadron(int pdgId)
{

  if (pdgId >= 500 && pdgId <600) return true;
  else if (pdgId >= 5000 && pdgId < 6000) return true;
  else return false;
  
}

bool isBjet(fastjet::PseudoJet& jet)
{

  vector<fastjet::PseudoJet> constituents = jet.constituents();
  for (auto constituent : constituents) 
    if (constituent.user_index() == -5) return true;
  return false;
  
}

bool isBpartonJet(fastjet::PseudoJet& jet, int parton)
{

  vector<fastjet::PseudoJet> constituents = jet.constituents();
  for (auto constituent : constituents) 
    if (constituent.user_index() == (-50000-parton)) return true;
  return false;
  
}

std::vector<int> findPartonHadrons(int parton, std::vector<GenParticle_p5>& partList)
{

  vector<int> retval;
  
  vector<int> prodVtx;
  prodVtx.push_back(partList[parton].m_endVtx);

  for (int index = parton; index < partList.size(); index++) {
    //    for (auto ii : prodVtx) cout << ii << " "; cout << endl;
    if (find(prodVtx.begin(), prodVtx.end(), partList[index].m_prodVtx) != prodVtx.end()) {
      if (partList[index].m_status == 1)
	retval.push_back(index);
      else
	prodVtx.push_back(partList[index].m_endVtx);
    }
  }
  return retval;
}


std::vector< std::vector<int> > findBHdecays(std::vector<GenParticle_p5>& partList)
{

  int index = -1;
  int prodVtx1 = 1;
  int prodVtx2 = 1;
  std::vector<int> vec1;
  std::vector<int> vec2;
  std::vector<std::vector<int>> retval;
  for (auto part : partList) {
    index++;
    if ((abs(part.m_pdgId) == 36) && part.m_status == 22) {
      if (prodVtx1 > 0)
	prodVtx1 = part.m_endVtx;
      else
	prodVtx2 = part.m_endVtx;
    }
    if ((abs(part.m_pdgId) == 36) && part.m_prodVtx == prodVtx1) {
      prodVtx1 = part.m_endVtx;
    }
    if ((abs(part.m_pdgId) == 36) && part.m_prodVtx == prodVtx2) {
      prodVtx2 = part.m_endVtx;
    }
    if (abs(part.m_pdgId) == 5 && part.m_prodVtx == prodVtx1) {
      vec1.push_back(index);
    }
    if (abs(part.m_pdgId) == 5 && part.m_prodVtx == prodVtx2) {
      vec2.push_back(index);
    }
  }
  retval.push_back(vec1);
  retval.push_back(vec2);
  return retval;
}    
  


  
