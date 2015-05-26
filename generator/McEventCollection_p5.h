#ifndef GENERATOROBJECTS_MCEVENTCOLLECTION_P5_H 
#define GENERATOROBJECTS_MCEVENTCOLLECTION_P5_H 

#include <vector>

#include "generator/GenEvent_p5.h"
#include "generator/GenParticle_p5.h"
#include "generator/GenVertex_p5.h"

#include "TObject.h"

class McEventCollection_p5 : public TObject {
  
 public: 

  McEventCollection_p5();
  ~McEventCollection_p5();

  std::vector<GenEvent_p5> m_genEvents;
  std::vector<GenVertex_p5> m_genVertices;
  std::vector<GenParticle_p5> m_genParticles;

  ClassDef(McEventCollection_p5, 1);
  
};


#endif
