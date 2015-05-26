#ifndef GENERATOROBJECTS_GENVERTEX_P5_H 
#define GENERATOROBJECTS_GENVERTEX_P5_H 

#include <vector>

#include "TObject.h"

class GenVertex_p5 : public TObject
{ 

 public: 

  GenVertex_p5();
  ~GenVertex_p5();

  float m_x;
  float m_y;
  float m_z;
  float m_t;
  std::vector<int> m_particlesIn;
  std::vector<int> m_particlesOut;
  int m_id;
  std::vector<float> m_weights;
  int m_barcode;

  ClassDef(GenVertex_p5, 1);

}; 


#endif
