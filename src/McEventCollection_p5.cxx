///////////////////////// -*- C++ -*- /////////////////////////////
// McEventCollectionCnv_p5.cxx 
// Implementation file for class McEventCollectionCnv_p5
// Author: S.Binet<binet@cern.ch>
/////////////////////////////////////////////////////////////////// 

// STL includes
#include <iostream>
#include <utility>
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// GeneratorObjectsTPCnv includes
#include "generator/GenParticle_p5.h"
#include "generator/GenVertex_p5.h"
#include "generator/GenEvent_p5.h"
#include "generator/McEventCollection_p5.h"

ClassImp(McEventCollection_p5)

McEventCollection_p5::McEventCollection_p5() :
  m_genEvents( ),
  m_genVertices( ),
  m_genParticles( )
{
  
}

McEventCollection_p5::~McEventCollection_p5()
{

  std::cout << "Deleting events" << std::endl;
  m_genEvents.clear();
  std::cout << "Deleting vertices" << std::endl;
  m_genVertices.clear();
  std::cout << "Deleting particles" << std::endl;
  m_genParticles.clear();

}

