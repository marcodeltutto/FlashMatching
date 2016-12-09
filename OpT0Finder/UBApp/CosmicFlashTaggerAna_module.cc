////////////////////////////////////////////////////////////////////////
// Class:       CosmicFlashTaggerAna
// Plugin Type: analyzer (art v2_05_00)
// File:        CosmicFlashTaggerAna_module.cc
//
// Generated at Fri Dec  9 09:44:39 2016 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

class CosmicFlashTaggerAna;


class CosmicFlashTaggerAna : public art::EDAnalyzer {
public:
  explicit CosmicFlashTaggerAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicFlashTaggerAna(CosmicFlashTaggerAna const &) = delete;
  CosmicFlashTaggerAna(CosmicFlashTaggerAna &&) = delete;
  CosmicFlashTaggerAna & operator = (CosmicFlashTaggerAna const &) = delete;
  CosmicFlashTaggerAna & operator = (CosmicFlashTaggerAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.

};


CosmicFlashTaggerAna::CosmicFlashTaggerAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void CosmicFlashTaggerAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(CosmicFlashTaggerAna)
