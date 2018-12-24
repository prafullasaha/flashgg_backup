#ifndef THQLeptonicFILTER_h
#define THQLeptonicFILTER_h


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;
using namespace edm;
class tHqLeptonicFilter : public edm::EDFilter {
 public:
  explicit tHqLeptonicFilter(const edm::ParameterSet&);
  ~tHqLeptonicFilter();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob(); 
 private:

// ----------member data ---------------------------

EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
  int n_thqLeptonicEvents=0;
  int n_thqmuonEvents=0;
  int n_thqeleEvents=0;
};
#endif
