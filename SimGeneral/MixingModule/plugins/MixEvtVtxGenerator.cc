#ifndef HI_MixEvtVtxGenerator_H
#define HI_MixEvtVtxGenerator_H
/*
*/
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Provenance/interface/Provenance.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TMatrixD.h"
#include <iostream>

using namespace edm;
using namespace std;


namespace HepMC {
   class FourVector ;
}


class MixEvtVtxGenerator : public edm::EDProducer
{
   public:
  
  // ctor & dtor
  explicit MixEvtVtxGenerator( const edm::ParameterSet& );
  virtual ~MixEvtVtxGenerator();
  
  virtual void produce( edm::Event&, const edm::EventSetup& ) override;
  
  virtual HepMC::FourVector* getVertex(edm::Event&);
  virtual HepMC::FourVector* getRecVertex(edm::Event&);
  
  protected:
  
  HepMC::FourVector*       fVertex ;
  TMatrixD *boost_;
  
  private :
  
  edm::EDGetTokenT<reco::VertexCollection>   hiLabel;
  edm::EDGetTokenT<HepMCProduct>  signalLabel;
  edm::EDGetTokenT<CrossingFrame<HepMCProduct> >   cfLabel;
  
   bool                     useRecVertex;
   std::vector<double>      vtxOffset;
  bool useCF_;

};

MixEvtVtxGenerator::MixEvtVtxGenerator( const ParameterSet& pset ) 
	: fVertex(0), boost_(0),
	  useRecVertex(pset.exists("useRecVertex")?pset.getParameter<bool>("useRecVertex"):false)
	  
{   
   produces<bool>("matchedVertex"); 
   vtxOffset.resize(3);
   if(pset.exists("vtxOffset")) vtxOffset=pset.getParameter< std::vector<double> >("vtxOffset");

   if(useRecVertex) useCF_ = 0;
   else{
     useCF_ = pset.getUntrackedParameter<bool>("useCF",false);
     cfLabel = consumes<CrossingFrame<HepMCProduct> >(pset.getParameter<edm::InputTag>("mixLabel"));
   }

}

MixEvtVtxGenerator::~MixEvtVtxGenerator() 
{
   delete fVertex ;
   if (boost_ != 0 ) delete boost_;
   // no need since now it's done in HepMCProduct
   // delete fEvt ;
}

HepMC::FourVector* MixEvtVtxGenerator::getVertex( Event& evt){

  HepMC::GenVertex* genvtx = 0;
  const HepMC::GenEvent* inev = 0;

  //cout<<" use CF "<<useCF_<<endl;
  
  if(useCF_){
    Handle<CrossingFrame<HepMCProduct> > cf;
    evt.getByToken(cfLabel,cf);
    MixCollection<HepMCProduct> mix(cf.product());
    if(mix.size() < 2){
      cout<<"Less than 2 sub-events, mixing seems to have failed!"<<endl;
    }
    const HepMCProduct& bkg = mix.getObject(1);
    if(!(bkg.isVtxGenApplied())){
      cout<<"Input background does not have smeared vertex!"<<endl;
    }else{
      inev = bkg.GetEvent();
   } 
  }else{
    //cout<<" hiLabel "<<hiLabel<<endl;
    Handle<HepMCProduct> input;
    evt.getByToken(signalLabel,input);
    inev = input->GetEvent();
  }

  genvtx = inev->signal_process_vertex();
  if(!genvtx){
    //cout<<"No Signal Process Vertex!"<<endl;
    HepMC::GenEvent::particle_const_iterator pt=inev->particles_begin();
    HepMC::GenEvent::particle_const_iterator ptend=inev->particles_end();
    while(!genvtx || ( genvtx->particles_in_size() == 1 && pt != ptend ) ){
      //if(!genvtx) cout<<"No Gen Vertex!"<<endl;
      if(pt == ptend) cout<<"End reached, No Gen Vertex!"<<endl;
      genvtx = (*pt)->production_vertex();
      ++pt;
    }
  }
  double aX,aY,aZ,aT;

  aX = genvtx->position().x();
  aY = genvtx->position().y();
  aZ = genvtx->position().z();
  aT = genvtx->position().t();
  
  if(!fVertex){
    fVertex = new HepMC::FourVector();
    //cout<<" creating new vertex "<<endl;
  }
  cout<<" setting vertex "<<" aX "<<aX<<" aY "<<aY<<" aZ "<<aZ<<" aT "<<aT<<endl;
  fVertex->set(aX,aY,aZ,aT);
  

  return fVertex;

}


HepMC::FourVector* MixEvtVtxGenerator::getRecVertex( Event& evt){

  Handle<reco::VertexCollection> input;
  evt.getByToken(hiLabel,input);

  double aX,aY,aZ;

  aX = input->begin()->position().x() + vtxOffset[0];
  aY = input->begin()->position().y() + vtxOffset[1];
  aZ = input->begin()->position().z() + vtxOffset[2];

  /*
  std::cout << "reco::Vertex = " << input->begin()->position().x()
	    << ", " << input->begin()->position().y()
	    << ", " << input->begin()->position().z()
	    << std::endl;

  std::cout << "offset = " << vtxOffset[0]
	    << ", " << vtxOffset[1]
	    << ", " << vtxOffset[2]
	    << std::endl;

  std::cout << "embedded GEN vertex = " << aX
	    << ", " << aY << ", " << aZ << std::endl;

  */
  if(!fVertex) fVertex = new HepMC::FourVector();
  fVertex->set(10.0*aX,10.0*aY,10.0*aZ,0.0); // HepMC positions in mm (RECO in cm)
  
  return fVertex;

}

void MixEvtVtxGenerator::produce( Event& evt, const EventSetup& )
{
   
   
   Handle<HepMCProduct> HepMCEvt ;
   
   evt.getByToken( signalLabel, HepMCEvt ) ;
   
   // generate new vertex & apply the shift 
   //
   if(HepMCEvt->isVtxGenApplied()) throw cms::Exception("MatchVtx")
				      <<"Signal HepMCProduct is not compatible for embedding - it's vertex is already smeared."
				      <<std::endl;
   HepMCEvt->applyVtxGen( useRecVertex ? getRecVertex(evt) : getVertex(evt) ) ;

   //   HepMCEvt->boostToLab( GetInvLorentzBoost(), "vertex" );
   //   HepMCEvt->boostToLab( GetInvLorentzBoost(), "momentum" );
   
   // OK, create a (pseudo)product and put in into edm::Event
   //
   auto_ptr<bool> NewProduct(new bool(true)) ;      
   evt.put( NewProduct ,"matchedVertex") ;
      
   return ;

}

DEFINE_FWK_MODULE(MixEvtVtxGenerator);

#endif
