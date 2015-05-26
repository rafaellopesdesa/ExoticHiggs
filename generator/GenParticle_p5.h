#ifndef GENERATOROBJECTS_GENPARTICLE_P5_H 
#define GENERATOROBJECTS_GENPARTICLE_P5_H 

#include <vector>
#include <utility> //> for std::pair

#include "TObject.h"

class GenParticle_p5 : public TObject
{ 

 public: 

  GenParticle_p5();
  ~GenParticle_p5();

  float m_px;
  float m_py;
  float m_pz;
  float m_m;
  int m_pdgId;
  int m_status;
  std::vector< std::pair<int, int> > m_flow;
  float m_thetaPolarization;
  float m_phiPolarization;
  int m_prodVtx;
  int m_endVtx;
  int m_barcode;
  float m_generated_mass;
  short m_recoMethod;

  ClassDef(GenParticle_p5, 1);
  
}; 



#endif
