#include "flashgg/MicroAOD/interface/tHqLeptonicFilter.h"

tHqLeptonicFilter::tHqLeptonicFilter(const edm::ParameterSet& iConfig) :
genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag>( "genParticleTag" ) ) )
{
}

tHqLeptonicFilter::~tHqLeptonicFilter()
{
}

// member functions

// ------------ method called on each new Event  ------------

bool
tHqLeptonicFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool accepted = false;
  bool accepted1 = false;
  Handle<View<reco::GenParticle> > genParticles;
  iEvent.getByToken( genParticleToken_, genParticles );

// -- look for Higgs decaying to two photons

for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
    Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
    if (gen->pdgId() == 25){

//cout << "Higgs found " <<endl;
         if (gen->numberOfDaughters() !=2 ) continue;
            const reco::Candidate * d1 = gen->daughter( 0 );
            const reco::Candidate * d2 = gen->daughter( 1 );

      //cout << "Higgs with status = " << gen->status() << " has two daughters with pdgId: " << endl;
      //      //cout << "d1 pdgId = " << d1->pdgId() << "   d2 pdgId = "<< d2->pdgId() <<endl;
     if (d1->pdgId()!=22 || d2->pdgId()!=22) continue;
      accepted = true;
      break;
    }
  }
 
  for( unsigned int j = 0 ; j < genParticles->size(); j++ ) {
    Ptr<reco::GenParticle> gen1 = genParticles->ptrAt(j);
    if (std::abs(gen1->pdgId()) != 24)    continue;
    if (gen1->numberOfDaughters() !=2 ) continue;
      const reco::Candidate * t1 = gen1->daughter( 0 );
      const reco::Candidate * t2 = gen1->daughter( 1 );
      if ((std::abs(t1->pdgId())==11 && std::abs( t2->pdgId())==12) || (std::abs(t1->pdgId())==12 && std::abs(t2->pdgId())==11)){
        accepted1 = true;
	n_thqeleEvents++; 

//        cout<<"t1 & t2="<<t1->pdgId()<<"   "<<t2->pdgId()<<endl;
        break;
       }
      else if ((std::abs(t1->pdgId())==13 && std::abs(t2->pdgId())==14) || (std::abs(t1->pdgId())==14 && std::abs(t2->pdgId())==13)){
        accepted1 = true;
        n_thqmuonEvents++;
//        cout<<"t1 & t2="<<t1->pdgId()<<"   "<<t2->pdgId()<<endl;

        break;        
       }
  }

if (accepted && accepted1) {
    n_thqLeptonicEvents++;
    return true; 
  } else { 
    return false;
  }
  
}
void tHqLeptonicFilter::endJob(){
  std::cout<<"Total Number of thqLeptonic Events: "<<n_thqLeptonicEvents<<std::endl;
  std::cout<<"Total Number of thqmuon Events: "<<n_thqmuonEvents<<std::endl;
  std::cout<<"Total Number of thqele Events: "<<n_thqeleEvents<<std::endl;
  }
DEFINE_FWK_MODULE(tHqLeptonicFilter);





