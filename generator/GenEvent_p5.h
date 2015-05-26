#ifndef GENERATOROBJECTS_GENEVENT_P5_H 
#define GENERATOROBJECTS_GENEVENT_P5_H 

// STL includes
#include <vector>
#include <utility> //> for std::pair

#include "TObject.h"

class GenEvent_p5 : public TObject
{ 
 public: 

  GenEvent_p5();
  ~GenEvent_p5();

  int m_signalProcessId;
  int m_eventNbr;
  int m_mpi;
  double m_eventScale;
  double m_alphaQCD;
  double m_alphaQED;
  int m_signalProcessVtx;
  int m_beamParticle1;
  int m_beamParticle2;
  std::vector<double> m_weights;
  std::vector<long int> m_randomStates;
  std::vector<double> m_crossSection;
  std::vector<float> m_heavyIon;
  std::vector<double> m_pdfinfo;
  int m_momentumUnit;
  int m_lengthUnit;
  unsigned int m_verticesBegin;
  unsigned int m_verticesEnd;
  unsigned int m_particlesBegin;
  unsigned int m_particlesEnd;

  ClassDef(GenEvent_p5, 1);

}; 


#endif
