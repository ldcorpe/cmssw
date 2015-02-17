// -*- C++ -*-
//
//

// user include files
#include "FWCore/Framework/interface/MakerMacros.h"
//#include "DataFormats/EgammaReco/interface/SuperCluster.h"
//#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "Calibration/EcalCalibAlgos/interface/ElectronRecalibSuperClusterAssociator.h"

//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include <iostream>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#define DEBUG 

using namespace reco;
using namespace edm;

//#define P4_FROM_SUPER_CLUSTER 0
 
ElectronRecalibSuperClusterAssociator::ElectronRecalibSuperClusterAssociator(const edm::ParameterSet& iConfig) 
{
#ifdef DEBUG
  std::cout<< "ElectronRecalibSuperClusterAssociator::ElectronRecalibSuperClusterAssociator" << std::endl;
#endif

  //register your products
  produces<GsfElectronCollection>();
  produces<GsfElectronCoreCollection>() ;
  //  produces<SuperClusterCollection>();
  
  superClusterCollectionEB_ = iConfig.getParameter<edm::InputTag > ("superClusterCollectionEB");
  superClusterCollectionEE_ = iConfig.getParameter<edm::InputTag > ("superClusterCollectionEE");

  outputLabel_ = iConfig.getParameter<std::string>("outputLabel");
  electronSrc_ = iConfig.getParameter<edm::InputTag > ("electronSrc");

  electronToken_ = consumes<reco::GsfElectronCollection>(electronSrc_);
  ebScToken_     = consumes<reco::SuperClusterCollection>(superClusterCollectionEB_);
  eeScToken_     = consumes<reco::SuperClusterCollection>(superClusterCollectionEE_);
  
#ifdef DEBUG
  std::cout<< "ElectronRecalibSuperClusterAssociator::ElectronRecalibSuperClusterAssociator::end" << std::endl;
#endif
}

ElectronRecalibSuperClusterAssociator::~ElectronRecalibSuperClusterAssociator()
{
}

// ------------ method called to produce the data  ------------
void ElectronRecalibSuperClusterAssociator::produce(edm::Event& e, const edm::EventSetup& iSetup) 
{
#ifdef DEBUG
  std::cout<< "GEDElectronRecalibSuperClusterAssociator::produce" << std::endl;
#endif

  // Create the output collections   
  std::auto_ptr<GsfElectronCollection> pOutEle(new GsfElectronCollection);
  std::auto_ptr<GsfElectronCoreCollection> pOutEleCore(new GsfElectronCoreCollection);

  // Get SuperClusters in EB
  Handle<reco::SuperClusterCollection> superClusterEBHandle;
  e.getByToken(ebScToken_, superClusterEBHandle); 
  //const reco::SuperClusterCollection* scCollection = superClusterEBHandle.product();
  
#ifdef DEBUG
  std::cout<<"EB scCollection->size()" << superClusterEBHandle->size()<<std::endl;
#endif

  // Get SuperClusters in EE
  Handle<reco::SuperClusterCollection> superClusterEEHandle;
  e.getByToken(eeScToken_, superClusterEEHandle);
  //  const reco::SuperClusterCollection* eeScCollection = superClusterEEHandle.product();

#ifdef DEBUG
  std::cout<<"EE scCollection->size()" << superClusterEEHandle->size() << std::endl;
#endif

  // Get Electrons
  edm::Handle<reco::GsfElectronCollection> eleHandle;
  e.getByToken(electronToken_, eleHandle);
  //  const reco::GsfElectronCollection* electronCollection = eleHandle.product();

  GsfElectronCoreRefProd rEleCore=e.getRefBeforePut<GsfElectronCoreCollection>();
  edm::Ref<GsfElectronCoreCollection>::key_type idxEleCore = 0;
  
  for(reco::GsfElectronCollection::const_iterator eleIt = eleHandle->begin(); eleIt != eleHandle->end(); eleIt++)
    {
      float DeltaRMineleSCbarrel(0.15); //initial minDeltaR
      float DeltaRMineleSCendcap(0.15); 
      const reco::SuperCluster* nearestSCbarrel=0;
      const reco::SuperCluster* nearestSCendcap=0;
      int iscRef=-1, iscRefendcap=-1;
      int iSC=0;
      
      if(eleIt->trackerDrivenSeed()){
	edm::LogError("trackerDriven") << "skipping trackerDriven electrons";
	continue;
      }
      // first loop is on EB superClusters
      iSC=0;
      for(reco::SuperClusterCollection::const_iterator scIt = superClusterEBHandle->begin();
	  scIt != superClusterEBHandle->end(); scIt++, iSC++){

	double DeltaReleSC = sqrt(reco::deltaR2(eleIt->eta(), eleIt->phi(),
						scIt->eta(), scIt->phi()));
	
	if(DeltaReleSC<DeltaRMineleSCbarrel) //save the nearest SC
	  {
	    DeltaRMineleSCbarrel = DeltaReleSC;
	    nearestSCbarrel = &*scIt;
	    iscRef = iSC;
	  }
#ifdef DEBUG	
	std::cout << "EB: " << scIt - superClusterEBHandle->begin() << " " << iSC << " " << iscRef 
		  << "\t" << std::setprecision(4) << scIt->energy() 
		  << " " << scIt->eta() << " " << scIt->phi() 
		  << "\t--\t" << eleIt->energy() << " " << eleIt->eta() << " " << eleIt->phi() 
		  << "\t" << DeltaRMineleSCbarrel
		  << std::endl;
#endif
      }
      
      // second loop is on EE superClusters
      iSC=0;
      for(reco::SuperClusterCollection::const_iterator scIt = superClusterEEHandle->begin();
	  scIt != superClusterEEHandle->end(); scIt++, iSC++){
#ifdef DEBUG	
	std::cout << "EE: " << scIt - superClusterEEHandle->begin() << " " << iSC << " " << iscRef 
		  << "\t" << std::setprecision(4) << scIt->energy() 
		  << " " << scIt->eta() << " " << scIt->phi() 
		  << "\t--\t " << eleIt->energy() << " " << eleIt->eta() << " " << eleIt->phi() 
		  << "\t" << DeltaRMineleSCendcap
		  << std::endl;
#endif
	
	double DeltaReleSC = sqrt(reco::deltaR2(eleIt->eta(), eleIt->phi(),
						scIt->eta(), scIt->phi()));
	
	if(DeltaReleSC<DeltaRMineleSCendcap)
	  {
	    DeltaRMineleSCendcap = DeltaReleSC;
	    nearestSCendcap = &*scIt;
	    iscRefendcap = iSC;
	  }
      }
      ////////////////////////      
      //      if(eleIt->isEB()) assert(DeltaRMineleSCbarrel < DeltaRMineleSCendcap);
      //else assert(DeltaRMineleSCbarrel > DeltaRMineleSCendcap);
      if(eleIt->isEB() && DeltaRMineleSCbarrel > DeltaRMineleSCendcap){
	edm::LogError("ElectronRecalibAssociator") << "EB electron, but nearest SC is in EE";;
	continue;
      }

      if(eleIt->isEB() && nearestSCbarrel){
#ifdef DEBUG
	std::cout << "Starting Association is with EB superCluster "<< std::endl;
#endif  
	
	pOutEleCore->push_back(*eleIt->core());
	std::cout << "pushed" << std::endl;
	reco::GsfElectronCore & newEleCore = pOutEleCore->back();
	std::cout << "backed" << std::endl;
	newEleCore.setGsfTrack(eleIt->gsfTrack());           // gsf track
	std::cout << "set gsf track ref" << std::endl;
	newEleCore.setSuperCluster(eleIt->superCluster());   // refined supercluster
	
	reco::SuperClusterRef scRef(reco::SuperClusterRef(superClusterEBHandle, iscRef));  
#ifndef CMSSW_5_3_X
        newEleCore.setParentSuperCluster(scRef);             // mustache 
#endif

	//	continue; 
	std::cout << "take reference to core" << std::endl;
	reco::GsfElectronCoreRef newEleCoreRef(rEleCore, idxEleCore++);
	std::cout << "ref made" << std::endl;
	reco::GsfElectron newEle(*eleIt,newEleCoreRef);
	std::cout << "created ele" << std::endl;
	//	std::cout << newEleCore.gsfTrack()->pMode() << std::endl;
	//std::cout << newEleCoreRef->gsfTrack()->pMode() << std::endl;
	
	continue;
	//-- first possibility: set the new p4SC using refined SC
	
	//newEle.setP4(reco::GsfElectron::P4_FROM_SUPER_CLUSTER, eleIt->p4(reco::GsfElectron::P4_FROM_SUPER_CLUSTER)); //*newEle.superCluster()->energy()/eleIt->superCluster()->energy());
	
	//-- second possibility: set the new p4SC using mustache SC
	
	newEle.setP4(reco::GsfElectron::P4_FROM_SUPER_CLUSTER, eleIt->p4(reco::GsfElectron::P4_FROM_SUPER_CLUSTER)*newEle.parentSuperCluster()->energy()/eleIt->parentSuperCluster()->energy(), eleIt->p4Error(reco::GsfElectron::P4_FROM_SUPER_CLUSTER), false); 

	//-- update the correctedEcalEnergy
	//newEle.setCorrectedEcalEnergy(eleIt->ecalEnergy()*(newEle.superCluster()->energy()/eleIt->superCluster()->energy()));
	newEle.setCorrectedEcalEnergy(eleIt->ecalEnergy()*(newEle.parentSuperCluster()->energy()/
							   eleIt->parentSuperCluster()->energy()));
  
	newEle.setCorrectedEcalEnergyError(eleIt->ecalEnergyError()*(nearestSCbarrel->energy()/eleIt->ecalEnergy()));
	pOutEle->push_back(newEle);
      }  else if(!(eleIt->isEB()) && nearestSCendcap)
	{
	  continue;
#ifdef DEBUG
	  std::cout << "Starting Association is with EE superCluster "<< std::endl;
#endif  
	  reco::GsfElectronCore newEleCore(*(eleIt->core()));
	  std::cout << "new ele core" << std::endl;
	  newEleCore.setGsfTrack(eleIt->gsfTrack());
	  std::cout << "new gsfTrack to ele core" << std::endl;
	  reco::SuperClusterRef scRef(superClusterEEHandle, iscRefendcap);
	  std::cout << "new sc ref" << scRef->energy() << std::endl;
	  newEleCore.setSuperCluster(scRef);
	  newEleCore.setParentSuperCluster(scRef);
	  std::cout << "new sc ref to ele core" << std::endl;
	  reco::GsfElectronCoreRef newEleCoreRef = edm::Ref<GsfElectronCoreCollection>(&(*pOutEleCore), idxEleCore++);
	  //reco::GsfElectronCoreRefProd newEleCoreRef(rEleCore, idxEleCore++);
	  std::cout << "new ele core ref" << std::endl;
	  pOutEleCore->push_back(newEleCore);
	  std::cout << "pushed back" << std::endl;
	  //continue;
	  reco::GsfElectron newEle(*eleIt,newEleCoreRef);
#ifdef DEBUG
	  //	  assert(newEle.superCluster()->seed()!=NULL);
#endif
	  std::cout << "newEle:" << newEle.superCluster()->energy() << "\t" << newEle.superCluster()->seed()->energy() <<  std::endl;
	  
	//-- second possibility: set the new p4SC using mustache SC
	  //newEle.setP4(reco::GsfElectron::P4_FROM_SUPER_CLUSTER, eleIt->p4(reco::GsfElectron::P4_FROM_SUPER_CLUSTER)*newEle.parentSuperCluster()->energy()/eleIt->parentSuperCluster()->energy(),eleIt->p4Error(reco::GsfElectron::P4_FROM_SUPER_CLUSTER), false); 

	//-- update the correctedEcalEnergy
	//newEle.setCorrectedEcalEnergy(eleIt->ecalEnergy()*(newEle.superCluster()->energy()/eleIt->superCluster()->energy()));
	//newEle.setCorrectedEcalEnergy(eleIt->ecalEnergy()*(newEle.parentSuperCluster()->energy()/
	//						   eleIt->parentSuperCluster()->energy()));
	//newEle.setCorrectedEcalEnergyError(eleIt->ecalEnergyError()*(nearestSCbarrel->energy()/eleIt->ecalEnergy()));
	  //	pOutEle->push_back(newEle);
	}else{
	edm::LogError("Failed SC association") << "No SC to be associated to the electron";
      }
    }
  
  
  
#ifdef DEBUG
  std::cout << "Filled new electrons  " << pOutEle->size() << std::endl;
  std::cout << "Filled new electronsCore  " << pOutEleCore->size() << std::endl;
  //  std::cout << "Filled new endcapSC  " << pOutNewEndcapSC->size() << std::endl;
#endif  
  
  // put result into the Event

  e.put(pOutEle);
  e.put(pOutEleCore);
  
  //  e.put(pOutNewEndcapSC);
  
}

DEFINE_FWK_MODULE(ElectronRecalibSuperClusterAssociator);
