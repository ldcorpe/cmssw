#include "RecoEcal/EgammaClusterAlgos/interface/HGCALShowerBasedEmIdentificationLC2.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "TPrincipal.h"

namespace {
template<typename DETID>
std::pair<int,int> getlayer(const unsigned rawid) {
DETID id(rawid);
return std::make_pair(id.zside(),id.layer());
}
}

HGCALShowerBasedEmIdentificationLC2::HGCALShowerBasedEmIdentificationLC2 (bool withPileup,const HGCalGeometry *geom_ ): 
 withPileup_(withPileup)
{

  geom = geom_;	
  // initialize showerPos and showerDir
  showerPos_ = math::XYZPoint(0.,0.,0.);
  showerDir_ = math::XYZVector(0.,0.,0.);
  showerPosIsSet_ = false;
  showerDirIsSet_ = false;
  
  // parameters
  mip_ = 0.0000551;
  minenergy_ = 4.;
  rmax_ = 100.; // no transverse limitation for no PU case
  if (withPileup_) rmax_ = 1.5*2.27;
  hovereConesize_ = 0.05;
    
  // HGCAL average medium
  criticalEnergy_ = 0.00536; // in GeV
  radiationLength_ = 0.968; // in cm
    
  // longitudinal parameters
  // mean values
  // shower max <T> = t0 + t1*lny
  // <alpha> = alpha0 + alpha1*lny
  // shower average = alpha/beta  
  meant0_ = -1.396;
  meant1_ = 1.007;
  meanalpha0_ = -0.0433;
  meanalpha1_ = 0.540;
  // sigmas
  // sigma(lnT) = 1 /sigmalnt0 + sigmalnt1*lny; 
  // sigma(lnalpha) = 1 /sigmalnt0 + sigmalnt1*lny; 
  sigmalnt0_ = -2.506;
  sigmalnt1_ = 1.245;
  sigmalnalpha0_ = -0.08442;
  sigmalnalpha1_ = 0.7904;
  // corr(lnalpha,lnt) = corrlnalpha0_+corrlnalphalnt1_*y
  corrlnalphalnt0_ = 0.7858;
  corrlnalphalnt1_ = -0.0232;
    
  // cut values, to be moved as configurable parameters
  cutStartPosition_ = 322.5;
  cutSigmaetaeta_ = 0.0055;
  if (withPileup_) cutSigmaetaeta_ = 0.00480;
  cutHoverem_ = 0.003;
  if (withPileup_) cutHoverem_ = 0.065;
  cutLengthCompatibility_ = 4.0;
  
  /*
  std::cout << "*** HGCAL ShowerBased EmIdentification ***" << std::endl;
  std::cout << "- max transverse radius: " << rmax_ << std::endl;
  std::cout << "- hovere cone size: " << hovereConesize_ << std::endl;
  std::cout << "- cut sigmaetaeta: " << cutSigmaetaeta_ << std::endl;
  std::cout << "- cut hoverem: " << cutHoverem_ << std::endl;
  std::cout << "- cut start position: " << cutStartPosition_ << std::endl;
  std::cout << "- cut length compatibility (in sigmas): " << cutLengthCompatibility_ << std::endl;
  */
  
}

HGCALShowerBasedEmIdentificationLC2::~HGCALShowerBasedEmIdentificationLC2 ()
{
//  delete pcaShowerAnalysis_;
}

void HGCALShowerBasedEmIdentificationLC2::setShowerPosition(const math::XYZPoint &pos)
{
  showerPos_ = pos;
  showerPosIsSet_ = true;
}

void HGCALShowerBasedEmIdentificationLC2::setShowerDirection(const math::XYZVector &dir)
{
  showerDir_ = dir;
  showerDirIsSet_ = true;
}

bool HGCALShowerBasedEmIdentificationLC2::isEm(const reco::PFCluster& clu,const edm::PtrVector<HGCRecHit>& rcs )
{
  //return (cutStartPosition(clu) &&  cutSigmaetaeta(clu) && cutHadOverEm(clu));
  return (cutStartPosition(clu,rcs) &&  cutSigmaetaeta(clu,rcs) && cutLengthCompatibility(clu,rcs));
} 

double  HGCALShowerBasedEmIdentificationLC2::firstLayerWithv2(float frac,const reco::PFCluster& clu,const  edm::SortedCollection<HGCRecHit>& rcs, const  edm::SortedCollection<HGCRecHit>& rcsHEF  )
{

	math::XYZPoint firstPos;
	float cluster_energy =0;
//	float cluster_energy_HEF =0;
//	std::cout << " HGCEF rechits " << rcsHEF.size() << std::endl;
	//------------------------------------------------------------------------------
	//FIRST need to ascertain FULL shower energy fro,HEF and EE
	//Do this by going through the rechits in each subdet & summing up the energy...
	//------------------------------------------------------------------------------

	//loop over rechits
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;// get detId

		// Part oif the problem is that we need to work with HGCEE and HGCEF rechits rather than the default, so do a matching...
		int matchIndex =-1; //matchIndex in case of HGCEE match 
	//	int matchIndexHEF =-1; // matchIndex in case of HGCHEF match


		//-------------------------
		//first consider HGCEE case
		//-------------------------
		if(id_.subdetId() == HGCEE){
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_)	{ // they will have dame detID
					matchIndex = irh;
					break;
				}
			}
			if( matchIndex >-1){		
				cluster_energy = cluster_energy + rcs[matchIndex].energy();
			}
		}

		//-------------------------
		//next consider HGCHEF case
		//-------------------------
/*		if(id_.subdetId() == HGCHEF){
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_)	{ // they will have dame detID
					matchIndexHEF = irh;
	//				std::cout << "	FOUND" << std::endl;
					break;
				}
			}
			if( matchIndexHEF >-1){		
				cluster_energy = cluster_energy + rcsHEF[matchIndexHEF].energy();
				cluster_energy_HEF = cluster_energy_HEF + rcsHEF[matchIndexHEF].energy();
			//	std::cout << "adding HGCHEF energy " << rcsHEF[matchIndexHEF].energy() << std::endl;
			}
		}*/
		//end of loop over rechits.
	}
	// now cluster_energy has the total enegry of cluster in HGCEE and HGCEF
//	std::cout <<" cluster e: EE " << cluster_energy - cluster_energy_HEF << ", HEF " << cluster_energy_HEF <<",  ie " << (100*cluster_energy_HEF )/cluster_energy<< " pc " << std::endl;

	//------------------------------------------------------------------------------
	//Now that we have teh total energy of the cluster, need to go through layer by
	//layer starting from the outer layers (HGCHEF then HGCEE) and se how far many 
	//layers we need to integrate before reaching frac % of the cluster energy...
	//Start with HEF, then move onto HEE.
	//------------------------------------------------------------------------------

	float energy_fraction =0;
	int nRecHits =0;
	int match_layer =0;
	auto  match_detector = rcs[0].id().subdetId(); // will probably be in HGCEE,
	// so pick subdetID of first HGCEE rechit... (so the "auto" gets decoded ok) 


	//------------------------------------------------------------------------------
	//Start by looping over HGCEHEF layers 12->0.
	//------------------------------------------------------------------------------

/*	for ( int iLayer = 12; iLayer > -1; iLayer--){
		// loop over all avialable rechits, but we will only select those in HGCHEF
		for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

			const DetId & id_ =  clu.hitsAndFractions()[ih].first;
			//skip rechits which are not in HGCHEF
			if (id_.subdetId() != HGCHEF){ continue;}
			//get layer
			int layer = getlayer<HGCHEDetId>(id_).second;
			//	skip recHit if not in the considered layer;
			if (layer != iLayer) {continue;}
			// now we do the matching to the HGCHEF recHit collection as above...
			int matchIndexHEF =-1;
			// loop over HGCHEF rechits to find the match
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_){
					matchIndexHEF = irh;
					break;
				}
			}
			if (matchIndexHEF ==-1 ) { 
				continue;
			}
			auto refhit  = (rcsHEF[matchIndexHEF]) ;
			nRecHits++;
			energy_fraction =energy_fraction + (refhit.energy())/cluster_energy;
	//		std::cout << "adding fraction" << (refhit.energy())/cluster_energy << std::endl;
		}
		if (energy_fraction > frac) {
			match_layer =iLayer;
			match_detector = HGCHEF;
			break;
		}
 	if (iLayer ==0) std::cout << " HGCHEF layer " << iLayer << ", energy is " << energy_fraction*100 << "pc of total "<< std::endl;
	}*/
	//------------------------------------------------------------------------------
	//Continue with HGCEE layers 30->0
	//------------------------------------------------------------------------------
	float layer_e =0;
	for ( int iLayer = 30; iLayer > -1; iLayer--){
	layer_e =0;
		// loop over all avialable rechits, but we will only select those in HGCEE
		for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

			const DetId & id_ =  clu.hitsAndFractions()[ih].first;
			//skip rechits which are not in HGCEE
			if (id_.subdetId() != HGCEE){ continue;}
			//get layer
			int layer = getlayer<HGCEEDetId>(id_).second;
			//	skip recHit if not in the considered layer;
			if (layer != iLayer) {continue;}
			// now we do the matching to the HGCHEF recHit collection as above...
			int matchIndex =-1;
			// loop over HGCEE rechits to find the match
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = (rcs[matchIndex]) ;
			nRecHits++;
			energy_fraction =energy_fraction + (refhit.energy())/cluster_energy;
			layer_e += refhit.energy();
		}
		if (energy_fraction > frac) {
			match_layer =iLayer;
			match_detector = HGCEE; 
			break;
		}
//  std::cout << " HGCEE layer " << iLayer << ", energy is " << layer_e << ", fraction " << energy_fraction*100 << "pc of total "<< std::endl;
	//		std::cout << " HGCEE layer " << iLayer << ", energy is " << energy_fraction*100 << "pc of total "<< std::endl;
	}

//std::cout <<" match layer " << match_layer << std::endl;
	//------------------------------------------------------------------------------
	//  Next  find the most energetic recHit in the layer and use that position..
	//------------------------------------------------------------------------------

	double maxfirstenergy=0.; 
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;

		if(match_detector == HGCEE){
			int	layer = getlayer<HGCEEDetId>(id_).second;
			if ( layer != match_layer) {continue;}
			int matchIndex =-1;
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = rcs[matchIndex];
			const GlobalPoint pos =  geom->getPosition(refhit.id());
			if (id_.det()==DetId::Forward) {
				if (refhit.energy() > maxfirstenergy) {
					firstPos = pos;
					maxfirstenergy = refhit.energy();
				}
			}
		} /*else if(match_detector == HGCHEF){
			int	layer = getlayer<HGCHEDetId>(id_).second;
			if ( layer != match_layer) {continue;}
			int matchIndex =-1;
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = rcsHEF[matchIndex];
			const GlobalPoint pos =  geom->getPosition(refhit.id());
			if (id_.det()==DetId::Forward) {
				if (refhit.energy() > maxfirstenergy) {
					firstPos = pos;
					maxfirstenergy = refhit.energy();
				}
			}
		}*/
	}

	//------------------------------------------------------------------------------
	// finally refine firstPos x and y using the meaured direction 
	//------------------------------------------------------------------------------
	if (!showerPosIsSet_ || !showerPosIsSet_) return firstPos.z();

	double lambda = (firstPos-showerPos_).z()/showerDir_.z();
	math::XYZPoint extraPos = showerPos_ + lambda*showerDir_;	 
	firstPos = extraPos;

	return firstPos.z();

}
double  HGCALShowerBasedEmIdentificationLC2::firstLayerWithPedro(float frac,const reco::PFCluster& clu,const  edm::SortedCollection<HGCRecHit>& rcs, const  edm::SortedCollection<HGCRecHit>& rcsHEF  )
{
	double x0[42]={0.08, //layer 1
	0.620, 0.620, 0.620, 0.620, 0.620, 0.620, 0.620, 0.620, 0.620, 0.620, // layer 2-11
	0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809, // layer 12,-21
	1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239, //layer 22-30
	3.580, // layer 31
	3.103,3.103,3.103,3.103,3.103,3.103,3.103,3.103,3.103,3.103,3.103} ;//layer 32-42
	double mipEE = 0.000055;
	double mipHEF = 0.000085;

	math::XYZPoint firstPos;
	float cluster_energy =0;
	float cluster_energy_HEF =0;
//	std::cout << " HGCEF rechits " << rcsHEF.size() << std::endl;
	//------------------------------------------------------------------------------
	//FIRST need to ascertain FULL shower energy fro,HEF and EE
	//Do this by going through the rechits in each subdet & summing up the energy...
	//------------------------------------------------------------------------------

	//loop over rechits
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;// get detId

		// Part oif the problem is that we need to work with HGCEE and HGCEF rechits rather than the default, so do a matching...
		int matchIndex =-1; //matchIndex in case of HGCEE match 
		int matchIndexHEF =-1; // matchIndex in case of HGCHEF match

			int iLayer = getlayer<HGCHEDetId>(id_).second;

		//-------------------------
		//first consider HGCEE case
		//-------------------------
		if(id_.subdetId() == HGCEE){
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_)	{ // they will have dame detID
					matchIndex = irh;
					break;
				}
			}
			if( matchIndex >-1){		
				cluster_energy = cluster_energy + x0[iLayer-1]*(rcs[matchIndex].energy()/mipEE);
			}
		}

		//-------------------------
		//next consider HGCHEF case
		//-------------------------
		if(id_.subdetId() == HGCHEF){
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_)	{ // they will have dame detID
					matchIndexHEF = irh;
	//				std::cout << "	FOUND" << std::endl;
					break;
				}
			}
			if( matchIndexHEF >-1){		
				cluster_energy = cluster_energy + x0[iLayer+29]*(rcs[matchIndexHEF].energy()/mipHEF);
		//		cluster_energy = cluster_energy + rcsHEF[matchIndexHEF].energy();
				cluster_energy_HEF = cluster_energy_HEF +  x0[iLayer+29]*(rcs[matchIndexHEF].energy()/mipHEF);
			//	std::cout << "adding HGCHEF energy " << rcsHEF[matchIndexHEF].energy() << std::endl;
			}
		}
		//end of loop over rechits.
	}
	// now cluster_energy has the total enegry of cluster in HGCEE and HGCEF
	 if (cluster_energy_HEF);
//	 std::cout <<" cluster e: EE " << cluster_energy - cluster_energy_HEF << ", HEF " << cluster_energy_HEF <<",  ie " << (100*cluster_energy_HEF )/cluster_energy<< " pc " << std::endl;

	//------------------------------------------------------------------------------
	//Now that we have teh total energy of the cluster, need to go through layer by
	//layer starting from the outer layers (HGCHEF then HGCEE) and se how far many 
	//layers we need to integrate before reaching frac % of the cluster energy...
	//Start with HEF, then move onto HEE.
	//------------------------------------------------------------------------------

	float energy_fraction =0;
	int nRecHits =0;
	int match_layer =0;
	float layer_e =0;
	auto  match_detector = rcs[0].id().subdetId(); // will probably be in HGCEE,
	// so pick subdetID of first HGCEE rechit... (so the "auto" gets decoded ok) 


	//------------------------------------------------------------------------------
	//Start by looping over HGCEHEF layers 12->0.
	//------------------------------------------------------------------------------

	for ( int iLayer = 12; iLayer > 0; iLayer--){
	layer_e =0;
		// loop over all avialable rechits, but we will only select those in HGCHEF
		for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

			const DetId & id_ =  clu.hitsAndFractions()[ih].first;
			//skip rechits which are not in HGCHEF
			if (id_.subdetId() != HGCHEF){ continue;}
			//get layer
			int layer = getlayer<HGCHEDetId>(id_).second;
			//	skip recHit if not in the considered layer;
			if (layer != iLayer) {continue;}
			// now we do the matching to the HGCHEF recHit collection as above...
			int matchIndexHEF =-1;
			// loop over HGCHEF rechits to find the match
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_){
					matchIndexHEF = irh;
					break;
				}
			}
			if (matchIndexHEF ==-1 ) { 
				continue;
			}
			auto refhit  = (rcsHEF[matchIndexHEF]) ;
			nRecHits++;
			energy_fraction =energy_fraction + (x0[iLayer+29]*(refhit.energy()/mipHEF))/cluster_energy;
			layer_e += (x0[iLayer+29]*(refhit.energy()/mipHEF));
	//		std::cout << "adding fraction" << (refhit.energy())/cluster_energy << std::endl;
		}
		if (energy_fraction > frac) {
			match_layer =iLayer;
			match_detector = HGCHEF;
			break;
		}
//  std::cout << " HGCHEF layer " << iLayer << ", energy is " << layer_e << ", fraction " << energy_fraction*100 << "pc of total "<< std::endl;
	}
	//------------------------------------------------------------------------------
	//Continue with HGCEE layers 30->0
	//------------------------------------------------------------------------------

	for ( int iLayer = 30; iLayer > 0; iLayer--){
	layer_e =0;
		// loop over all avialable rechits, but we will only select those in HGCEE
		for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

			const DetId & id_ =  clu.hitsAndFractions()[ih].first;
			//skip rechits which are not in HGCEE
			if (id_.subdetId() != HGCEE){ continue;}
			//get layer
			int layer = getlayer<HGCEEDetId>(id_).second;
			//	skip recHit if not in the considered layer;
			if (layer != iLayer) {continue;}
			// now we do the matching to the HGCHEF recHit collection as above...
			int matchIndex =-1;
			// loop over HGCEE rechits to find the match
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = (rcs[matchIndex]) ;
			nRecHits++;
			energy_fraction =energy_fraction + (x0[iLayer-1]*(refhit.energy()/mipEE))/cluster_energy;
			//energy_fraction =energy_fraction + (refhit.energy())/cluster_energy;
			layer_e += (x0[iLayer-1]*(refhit.energy()/mipEE));
		}
		if (energy_fraction > frac) {
			match_layer =iLayer;
			match_detector = HGCEE; 
			break;
		}
 // std::cout << " HGCEE layer " << iLayer << ", energy is " << layer_e << ", fraction " << energy_fraction*100 << "pc of total "<< std::endl;
	//		std::cout << " HGCEE layer " << iLayer << ", energy is " << energy_fraction*100 << "pc of total "<< std::endl;
	}
//std::cout <<" match layer " << match_layer << std::endl;
	//------------------------------------------------------------------------------
	//  Next  find the most energetic recHit in the layer and use that position..
	//------------------------------------------------------------------------------

	double maxfirstenergy=0.; 
//	int finalLayer =-1;
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;

		if(match_detector == HGCEE){
			int	layer = getlayer<HGCEEDetId>(id_).second;
			if ( layer != match_layer) {continue;}
			int matchIndex =-1;
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = rcs[matchIndex];
			const GlobalPoint pos =  geom->getPosition(refhit.id());
			if (id_.det()==DetId::Forward) {
				if (refhit.energy() > maxfirstenergy) {
					firstPos = pos;
					maxfirstenergy = refhit.energy();
				}
			}
		} /*else if(match_detector == HGCHEF){
			int	layer = getlayer<HGCHEDetId>(id_).second;
			if ( layer != match_layer) {continue;}
			int matchIndex =-1;
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = rcsHEF[matchIndex];
			const GlobalPoint pos =  geom->getPosition(refhit.id());
			if (id_.det()==DetId::Forward) {
				if (refhit.energy() > maxfirstenergy) {
					firstPos = pos;
					maxfirstenergy = refhit.energy();
				}
			}
		}*/
	}

	//------------------------------------------------------------------------------
	// finally refine firstPos x and y using the meaured direction 
	//------------------------------------------------------------------------------
	if (!showerPosIsSet_ || !showerPosIsSet_) return firstPos.z();

	double lambda = (firstPos-showerPos_).z()/showerDir_.z();
	math::XYZPoint extraPos = showerPos_ + lambda*showerDir_;	 
	firstPos = extraPos;

	return firstPos.z();

}
double  HGCALShowerBasedEmIdentificationLC2::firstLayerWith(float frac,const reco::PFCluster& clu,const  edm::SortedCollection<HGCRecHit>& rcs, const  edm::SortedCollection<HGCRecHit>& rcsHEF  )
{

	math::XYZPoint firstPos;
	float cluster_energy =0;
	float cluster_energy_HEF =0;
//	std::cout << " HGCEF rechits " << rcsHEF.size() << std::endl;
	//------------------------------------------------------------------------------
	//FIRST need to ascertain FULL shower energy fro,HEF and EE
	//Do this by going through the rechits in each subdet & summing up the energy...
	//------------------------------------------------------------------------------

	//loop over rechits
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;// get detId

		// Part oif the problem is that we need to work with HGCEE and HGCEF rechits rather than the default, so do a matching...
		int matchIndex =-1; //matchIndex in case of HGCEE match 
		int matchIndexHEF =-1; // matchIndex in case of HGCHEF match


		//-------------------------
		//first consider HGCEE case
		//-------------------------
		if(id_.subdetId() == HGCEE){
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_)	{ // they will have dame detID
					matchIndex = irh;
					break;
				}
			}
			if( matchIndex >-1){		
				cluster_energy = cluster_energy + rcs[matchIndex].energy();
			}
		}

		//-------------------------
		//next consider HGCHEF case
		//-------------------------
		if(id_.subdetId() == HGCHEF){
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_)	{ // they will have dame detID
					matchIndexHEF = irh;
	//				std::cout << "	FOUND" << std::endl;
					break;
				}
			}
			if( matchIndexHEF >-1){		
				cluster_energy = cluster_energy + rcsHEF[matchIndexHEF].energy();
				cluster_energy_HEF = cluster_energy_HEF + rcsHEF[matchIndexHEF].energy();
			//	std::cout << "adding HGCHEF energy " << rcsHEF[matchIndexHEF].energy() << std::endl;
			}
		}
		//end of loop over rechits.
	}
	// now cluster_energy has the total enegry of cluster in HGCEE and HGCEF
	// std::cout <<" cluster e: EE " << cluster_energy - cluster_energy_HEF << ", HEF " << cluster_energy_HEF <<",  ie " << (100*cluster_energy_HEF )/cluster_energy<< " pc " << std::endl;
if(cluster_energy_HEF);
	//------------------------------------------------------------------------------
	//Now that we have teh total energy of the cluster, need to go through layer by
	//layer starting from the outer layers (HGCHEF then HGCEE) and se how far many 
	//layers we need to integrate before reaching frac % of the cluster energy...
	//Start with HEF, then move onto HEE.
	//------------------------------------------------------------------------------

	float energy_fraction =0;
	int nRecHits =0;
	int match_layer =0;
	float layer_e =0;
	auto  match_detector = rcs[0].id().subdetId(); // will probably be in HGCEE,
	// so pick subdetID of first HGCEE rechit... (so the "auto" gets decoded ok) 


	//------------------------------------------------------------------------------
	//Start by looping over HGCEHEF layers 12->0.
	//------------------------------------------------------------------------------

	for ( int iLayer = 12; iLayer > -1; iLayer--){
	layer_e =0;
		// loop over all avialable rechits, but we will only select those in HGCHEF
		for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

			const DetId & id_ =  clu.hitsAndFractions()[ih].first;
			//skip rechits which are not in HGCHEF
			if (id_.subdetId() != HGCHEF){ continue;}
			//get layer
			int layer = getlayer<HGCHEDetId>(id_).second;
			//	skip recHit if not in the considered layer;
			if (layer != iLayer) {continue;}
			// now we do the matching to the HGCHEF recHit collection as above...
			int matchIndexHEF =-1;
			// loop over HGCHEF rechits to find the match
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_){
					matchIndexHEF = irh;
					break;
				}
			}
			if (matchIndexHEF ==-1 ) { 
				continue;
			}
			auto refhit  = (rcsHEF[matchIndexHEF]) ;
			nRecHits++;
			energy_fraction =energy_fraction + (refhit.energy())/cluster_energy;
			layer_e += refhit.energy();
	//		std::cout << "adding fraction" << (refhit.energy())/cluster_energy << std::endl;
		}
		if (energy_fraction > frac) {
			match_layer =iLayer;
			match_detector = HGCHEF;
			break;
		}
  //std::cout << " HGCHEF layer " << iLayer << ", energy is " << layer_e << ", fraction " << energy_fraction*100 << "pc of total "<< std::endl;
	}
	//------------------------------------------------------------------------------
	//Continue with HGCEE layers 30->0
	//------------------------------------------------------------------------------

	for ( int iLayer = 30; iLayer > -1; iLayer--){
	layer_e =0;
		// loop over all avialable rechits, but we will only select those in HGCEE
		for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

			const DetId & id_ =  clu.hitsAndFractions()[ih].first;
			//skip rechits which are not in HGCEE
			if (id_.subdetId() != HGCEE){ continue;}
			//get layer
			int layer = getlayer<HGCEEDetId>(id_).second;
			//	skip recHit if not in the considered layer;
			if (layer != iLayer) {continue;}
			// now we do the matching to the HGCHEF recHit collection as above...
			int matchIndex =-1;
			// loop over HGCEE rechits to find the match
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = (rcs[matchIndex]) ;
			nRecHits++;
			energy_fraction =energy_fraction + (refhit.energy())/cluster_energy;
			layer_e += refhit.energy();
		}
		if (energy_fraction > frac) {
			match_layer =iLayer;
			match_detector = HGCEE; 
			break;
		}
//  std::cout << " HGCEE layer " << iLayer << ", energy is " << layer_e << ", fraction " << energy_fraction*100 << "pc of total "<< std::endl;
	//		std::cout << " HGCEE layer " << iLayer << ", energy is " << energy_fraction*100 << "pc of total "<< std::endl;
	}
//std::cout <<" match layer " << match_layer << std::endl;
	//------------------------------------------------------------------------------
	//  Next  find the most energetic recHit in the layer and use that position..
	//------------------------------------------------------------------------------

	double maxfirstenergy=0.; 
//	int finalLayer =-1;
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;

		if(match_detector == HGCEE){
			int	layer = getlayer<HGCEEDetId>(id_).second;
			if ( layer != match_layer) {continue;}
			int matchIndex =-1;
			for (unsigned int irh=0; irh < rcs.size() ; irh++){
				if (rcs[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = rcs[matchIndex];
			const GlobalPoint pos =  geom->getPosition(refhit.id());
			if (id_.det()==DetId::Forward) {
				if (refhit.energy() > maxfirstenergy) {
					firstPos = pos;
					maxfirstenergy = refhit.energy();
				}
			}
		} /*else if(match_detector == HGCHEF){
			int	layer = getlayer<HGCHEDetId>(id_).second;
			if ( layer != match_layer) {continue;}
			int matchIndex =-1;
			for (unsigned int irh=0; irh < rcsHEF.size() ; irh++){
				if (rcsHEF[irh].id() == id_){
					matchIndex = irh;
					break;
				}
			}
			if (matchIndex ==-1 ) { 
				continue;
			}
			auto refhit  = rcsHEF[matchIndex];
			const GlobalPoint pos =  geom->getPosition(refhit.id());
			if (id_.det()==DetId::Forward) {
				if (refhit.energy() > maxfirstenergy) {
					firstPos = pos;
					maxfirstenergy = refhit.energy();
				}
			}
		}*/
	}

	//------------------------------------------------------------------------------
	// finally refine firstPos x and y using the meaured direction 
	//------------------------------------------------------------------------------
	if (!showerPosIsSet_ || !showerPosIsSet_) return firstPos.z();

	double lambda = (firstPos-showerPos_).z()/showerDir_.z();
	math::XYZPoint extraPos = showerPos_ + lambda*showerDir_;	 
	firstPos = extraPos;

	return firstPos.z();

}
math::XYZPoint HGCALShowerBasedEmIdentificationLC2::startPosition0(const reco::PFCluster& clu)
{
	math::XYZPoint firstPos;
	double zmin = 10000.0;
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {
		const auto& refhit = clu.recHitFractions()[ih].recHitRef();
		const auto& pos = refhit->position();
		const DetId & id_(refhit->detId()) ;
		if (id_.det()==DetId::Forward) {
			if (std::abs(pos.z())<zmin) {
				firstPos = pos;
				zmin = std::abs(pos.z());
			}
		}
	}
	// refine the first position estimation, taking the max energy in the first layer
	double maxfirstenergy=0.;
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {
		const auto& refhit = clu.recHitFractions()[ih].recHitRef();
		const auto& pos = refhit->position();
		const DetId & id_(refhit->detId());
		if (id_.det()==DetId::Forward) {
			if (std::abs(pos.z()) != zmin) continue;
			if (refhit->energy() > maxfirstenergy) {
				firstPos = pos;
				maxfirstenergy = refhit->energy();
			}
		}
	}
	// finally refine firstPos x and y using the meaured direction
	if (!showerPosIsSet_ || !showerPosIsSet_) return firstPos;
	double lambda = (firstPos-showerPos_).z()/showerDir_.z();
	math::XYZPoint extraPos = showerPos_ + lambda*showerDir_;
	firstPos = extraPos;
	return firstPos;
}
math::XYZPoint HGCALShowerBasedEmIdentificationLC2::startPosition3(const reco::PFCluster& clu,const  edm::SortedCollection<HGCRecHit>& rcs )
{
bool applyMinRecHits =1;
int minRecHits =3;

	math::XYZPoint firstPos;
	double zmin = 10000.0;  
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;
		int matchIndex =-1;

		for (unsigned int irh=0; irh < rcs.size() ; irh++){

			if (rcs[irh].id() == id_)
			{
				matchIndex = irh;
				break;
			}
		}

		if (matchIndex ==-1 ) { 
			continue;
		}
		auto refhit  = (rcs[matchIndex]) ;
		const GlobalPoint pos =  geom->getPosition(refhit.id());

		if (id_.det()==DetId::Forward) {      
			if (std::abs(pos.z())<zmin) {

				if (applyMinRecHits){
					int nRecHitInLayer =0;
					int	layer = getlayer<HGCEEDetId>(id_).second;
			//		std::cout << "--> found new potential match in layer " << layer <<", need " << minRecHits << " rechits to be valid." << std::endl;
					for (unsigned int ih2=0;ih2<clu.recHitFractions().size();++ih2) {

						const DetId & id2_ =  clu.hitsAndFractions()[ih2].first;
						int	layer2 = getlayer<HGCEEDetId>(id2_).second;
						if (layer == layer2) { nRecHitInLayer++;}
					}
					if (nRecHitInLayer >= minRecHits){
				//	std::cout << "----> OK!         recHits in layer..."<< nRecHitInLayer << std:: endl;
					firstPos = pos;
					zmin = std::abs(pos.z());
					} else {
			//		std::cout << "----> insufficient recHits in layer..."<< nRecHitInLayer << std:: endl;;
					}
				} else {

					firstPos = pos;
					zmin = std::abs(pos.z());

				}
			}
		}
	}

	double maxfirstenergy=0.; 
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;
		int matchIndex =-1;

		for (unsigned int irh=0; irh < rcs.size() ; irh++){

			if (rcs[irh].id() == id_)
			{

				matchIndex = irh;
				break;

			}
		}

		if (matchIndex ==-1 ) { 
			continue;
		}

		auto refhit  = rcs[matchIndex];
		const GlobalPoint pos =  geom->getPosition(refhit.id());
		if (id_.det()==DetId::Forward) {
			if (std::abs(pos.z()) != zmin) continue;
			if (refhit.energy() > maxfirstenergy) {
				firstPos = pos;
				maxfirstenergy = refhit.energy();
			}
		}
	}

	// finally refine firstPos x and y using the meaured direction 
	if (!showerPosIsSet_ || !showerPosIsSet_) return firstPos;

	double lambda = (firstPos-showerPos_).z()/showerDir_.z();
	math::XYZPoint extraPos = showerPos_ + lambda*showerDir_;	 
	firstPos = extraPos;

	return firstPos;

}
math::XYZPoint HGCALShowerBasedEmIdentificationLC2::startPosition(const reco::PFCluster& clu,const  edm::PtrVector<HGCRecHit>& rcs )
{
bool applyMinRecHits =0;
int minRecHits =0;

	math::XYZPoint firstPos;
	double zmin = 10000.0;  
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;
		int matchIndex =-1;

		for (unsigned int irh=0; irh < rcs.size() ; irh++){

			if (rcs[irh]->id() == id_)
			{
				matchIndex = irh;
				break;
			}
		}

		if (matchIndex ==-1 ) { 
			continue;
		}
		auto refhit  = *(rcs[matchIndex]) ;
		const GlobalPoint pos =  geom->getPosition(refhit.id());

		if (id_.det()==DetId::Forward) {      
			if (std::abs(pos.z())<zmin) {

				if (applyMinRecHits){
					int nRecHitInLayer =0;
					int	layer = getlayer<HGCEEDetId>(id_).second;
			//		std::cout << "--> found new potential match in layer " << layer <<", need " << minRecHits << " rechits to be valid." << std::endl;
					for (unsigned int ih2=0;ih2<clu.recHitFractions().size();++ih2) {

						const DetId & id2_ =  clu.hitsAndFractions()[ih2].first;
						int	layer2 = getlayer<HGCEEDetId>(id2_).second;
						if (layer == layer2) { nRecHitInLayer++;}
					}
					if (nRecHitInLayer >= minRecHits){
				//	std::cout << "----> OK!         recHits in layer..."<< nRecHitInLayer << std:: endl;
					firstPos = pos;
					zmin = std::abs(pos.z());
					} else {
			//		std::cout << "----> insufficient recHits in layer..."<< nRecHitInLayer << std:: endl;;
					}
				} else {

					firstPos = pos;
					zmin = std::abs(pos.z());

				}
			}
		}
	}

	double maxfirstenergy=0.; 
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {

		const DetId & id_ =  clu.hitsAndFractions()[ih].first;
		int matchIndex =-1;

		for (unsigned int irh=0; irh < rcs.size() ; irh++){

			if (rcs[irh]->id() == id_)
			{

				matchIndex = irh;
				break;

			}
		}

		if (matchIndex ==-1 ) { 
			continue;
		}

		auto refhit  = *rcs[matchIndex];
		const GlobalPoint pos =  geom->getPosition(refhit.id());
		if (id_.det()==DetId::Forward) {
			if (std::abs(pos.z()) != zmin) continue;
			if (refhit.energy() > maxfirstenergy) {
				firstPos = pos;
				maxfirstenergy = refhit.energy();
			}
		}
	}

	// finally refine firstPos x and y using the meaured direction 
	if (!showerPosIsSet_ || !showerPosIsSet_) return firstPos;

	double lambda = (firstPos-showerPos_).z()/showerDir_.z();
	math::XYZPoint extraPos = showerPos_ + lambda*showerDir_;	 
	firstPos = extraPos;

	return firstPos;

}

double HGCALShowerBasedEmIdentificationLC2::sigmaRR(const std::vector<HGCEEDetId>& detidstack, const edm::Handle<HGCRecHitCollection> recHits_ )
{
	double sigmaRR_global=0;	
	double Etot=0;
	double sigmaRR[30];
	double E[30];

	for (unsigned int iL =0 ; iL<30 ; iL++){ //loop over every layer.
	
	double sigmaRR_layer=0., sumnrj_layer=0.;
	double r_weighted=0; 
	int n=0;
	//loop over rechits in layer to find weighted average radial position 
	for (unsigned int ih=0;ih<detidstack.size();++ih) {
		n++;
		const DetId & id_ =  detidstack[ih];
		int layer = getlayer<HGCEEDetId>(id_).second;
		if (layer != (int) iL) continue; // enforce that the detid is in the correct layer.
		
		if (id_.det()==DetId::Forward) { //ensure we are in the forward part of detetcor.
			const GlobalPoint cellPos0 =  geom->getPosition(id_);
			math::XYZPoint cellPos (cellPos0.x(), cellPos0.y(), cellPos0.z());	     
			double cellR = cellPos.z()/sinh(cellPos.eta());
    	HGCRecHitCollection::const_iterator refhit = recHits_->find(detidstack[ih]);
			// apply energy cut cut
			//if (!withPileup_ || refhit->energy()>minenergy_*mip_) 
			{

					sumnrj_layer += refhit->energy();
					r_weighted += refhit->energy()*cellR;

			}
		}
	}
	
	r_weighted=r_weighted/sumnrj_layer;
	E[iL]=sumnrj_layer;
	
	//now we have avergae R in the cluster, can get sigma RR. 
	for (unsigned int ih=0;ih<detidstack.size(); ih++) {
		const DetId & id_ =  detidstack[ih];
		int layer = getlayer<HGCEEDetId>(id_).second;
		if (layer != (int) iL) continue;
		//	std::cout << "detId pos " << id_.position() << std::endl;
		if (id_.det()==DetId::Forward) {
			const GlobalPoint cellPos0 =  geom->getPosition(detidstack[ih]);
			math::XYZPoint cellPos (cellPos0.x(), cellPos0.y(), cellPos0.z());	     
			double cellR = cellPos.z()/sinh(cellPos.eta());
    	HGCRecHitCollection::const_iterator refhit = recHits_->find(detidstack[ih]);
		//	if (!withPileup_ || refhit->energy()>minenergy_*mip_) 
			{
				// simple transversal cut, later can refine as function of depth
					const double dR = fabs(r_weighted - cellR );
					sigmaRR_layer += dR*dR*(refhit->energy());
			}
		}
	}

	//std::cout << "debug 2" << std::endl;
	//std::cout << "debug 0 " << sigmaRR_layer << ", E[iL]] " << E[iL] <<std::endl;
	if (E[iL] !=0) sigmaRR_layer /= E[iL];
//	std::cout << "debug 1 " << sigmaRR_layer << std::endl;
	sigmaRR_layer = sqrt(sigmaRR_layer);
	//std::cout << "debug 2 " << sigmaRR_layer << std::endl;
	sigmaRR[iL] = sigmaRR_layer;
//	std::cout << "debug 3 " << sigmaRR[iL] << std::endl;
	//std::cout << "debug s_ee 2 " << sigmaetaeta << std::endl;
	
//	std::cout << "debug 4 " << sigmaRR_global << std::endl;
	sigmaRR_global += sigmaRR[iL]*sigmaRR[iL]*E[iL];
//	std::cout << "debug 5 " << sigmaRR_global << std::endl;
//	std::cout << "debug 6 " << Etot<< ", " << E[iL] << std::endl;
	Etot += E[iL];

	
	}
	
	sigmaRR_global /=Etot;
	sigmaRR_global =sqrt(sigmaRR_global);


	return sigmaRR_global;

}

double HGCALShowerBasedEmIdentificationLC2::sigmaetaeta(const reco::PFCluster& clu,const edm::PtrVector<HGCRecHit>& rcs )
{

	//	std::cout << "debug 0" << std::endl;
	double sigmaetaeta=0., sumnrj=0.;
	math::XYZPoint firstPos = startPosition(clu, rcs);

	//	std::cout << "debug s_ee 0 " << sigmaetaeta << std::endl;
	//	std::cout << "debug 1" << std::endl;
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {
		//		const auto& refhit = clu.recHitFractions()[ih].recHitRef();
		//		const DetId & id_ = refhit->detId() ;       
		const DetId & id_ =  clu.hitsAndFractions()[ih].first;
		int matchIndex =-1;
		//	std::cout << "rcs size " << rcs.size() << std::endl;
		//	std::cout << "detId pos " << id_.position() << std::endl;
		for (unsigned int irh=0; irh < rcs.size() ; irh++){

			//		if (irh%1000 ==0) { std::cout << "rechit " << irh << " pos " <<rcs[irh]->position() << std::endl ;}
			if (rcs[irh]->id() == id_)
			{

				//				refhit  = rcs[irh];
				matchIndex = irh;
				//	std::cout << " add " << i <<"th rechit, subtotal "<< subtotal<< std::endl;
				break;

			}
		}

		if (matchIndex ==-1 ) { 
			//std::cout << "no match" << std::endl;
			continue;}
		auto refhit  = *rcs[matchIndex];
		if (id_.det()==DetId::Forward) {

			const GlobalPoint cellPos0 =  geom->getPosition(refhit.id());
			math::XYZPoint cellPos (cellPos0.x(), cellPos0.y(), cellPos0.z());	     
			math::XYZVector radius, longitudinal, transverse;
			radius = cellPos - firstPos;
			// distances in local coordinates
			longitudinal =  (radius.Dot(showerDir_))*showerDir_.unit()/showerDir_.R();
			transverse = radius - longitudinal;
			// apply energy cut cut
			if (!withPileup_ || refhit.energy()>minenergy_*mip_) {
				// simple transversal cut, later can refine as function of depth
				if (!withPileup_ || transverse.R() < rmax_) {
					const double deta = (cellPos.eta()-showerPos_.eta());
					sigmaetaeta += deta*deta*refhit.energy();
					//	std::cout << "debug s_ee 1 " << sigmaetaeta << std::endl;
					sumnrj += refhit.energy();
				}
			}
		}
	}

	//std::cout << "debug 2" << std::endl;
	sigmaetaeta /= sumnrj;
	sigmaetaeta = sqrt(sigmaetaeta);
	//std::cout << "debug s_ee 2 " << sigmaetaeta << std::endl;

	// now correct the eta dependency
	double feta;
	constexpr double feta_0 = 0.00964148 - 0.01078431*1.5 + 0.00495703*1.5*1.5;
	const double clu_eta = std::abs(clu.eta()); 
	feta = 0.00964148 - clu_eta*(0.0107843 - 0.00495703*clu_eta);
	sigmaetaeta *= feta_0 / feta ;
	//	std::cout << "debug s_ee 3 " << sigmaetaeta << std::endl;

	//	std::cout << "debug 3" << std::endl;
	return sigmaetaeta;

}

double HGCALShowerBasedEmIdentificationLC2::lengthCompatibility(const reco::PFCluster& clu,const  edm::PtrVector<HGCRecHit>& rcs  )
{
	bool  photon =1;
	if (photon) radiationLength_=radiationLength_*(9/7); //photon radiation length is 9/7 higher than for electrons.
	// check that showerPos and showerDir have been set
	if (!showerPosIsSet_) {
		std::cout << "[HGCALShowerBasedEmIdentificationLC2::lengthCompatibility] error, showwer position not set " << std::endl;
		std::cout << "[HGCALShowerBasedEmIdentificationLC2::lengthCompatibility] error, please invoke setShowerPos and setShowerDir before invoking this function " << std::endl;
		return 0.;
	}

	double lengthCompatibility=0., predictedLength=0., predictedSigma=0.;

	// shower length	 
	const double length =  (showerPos_ - startPosition(clu, rcs)).R();
	const double lny = clu.emEnergy()/criticalEnergy_>1. ? std::log(clu.emEnergy()/criticalEnergy_) : 0.;

	// inject here parametrization results
	const double meantmax = meant0_ + meant1_*lny;
	const double meanalpha = meanalpha0_ + meanalpha1_*lny;
	const double sigmalntmax = 1.0 / (sigmalnt0_+sigmalnt1_*lny);
	const double sigmalnalpha = 1.0 / (sigmalnalpha0_+sigmalnalpha1_*lny);
	const double corrlnalphalntmax = corrlnalphalnt0_+corrlnalphalnt1_*lny;

	const double invbeta = meantmax/(meanalpha-1.);
	predictedLength = meanalpha*invbeta;
	predictedLength *= radiationLength_;

	double sigmaalpha = meanalpha*sigmalnalpha;
	if (sigmaalpha<0.) sigmaalpha = 1.;
	double sigmatmax = meantmax*sigmalntmax;
	if (sigmatmax<0.) sigmatmax = 1.;

	predictedSigma = sigmalnalpha*sigmalnalpha/((meanalpha-1.)*(meanalpha-1.));
	predictedSigma += sigmalntmax*sigmalntmax;
	predictedSigma -= 2*sigmalnalpha*sigmalntmax*corrlnalphalntmax/(meanalpha-1.);
	predictedSigma = predictedLength*sqrt(predictedSigma);
	if (photon) {
	predictedSigma = sqrt( predictedSigma*predictedSigma + radiationLength_*radiationLength_);
	radiationLength_= (7/9)*radiationLength_ ; //return radition length to the original value.
	}

	lengthCompatibility = (predictedLength-length)/predictedSigma;
	
	
	return lengthCompatibility;

}

bool HGCALShowerBasedEmIdentificationLC2::cutSigmaetaeta(const reco::PFCluster& clu,const edm::PtrVector<HGCRecHit>& rcs )
{
	return (sigmaetaeta(clu, rcs)<cutSigmaetaeta_);
}

bool HGCALShowerBasedEmIdentificationLC2::cutStartPosition(const reco::PFCluster& clu,const edm::PtrVector<HGCRecHit>& rcs )
{
	return ( std::abs( startPosition(clu,rcs).z() ) < cutStartPosition_ );
}

bool HGCALShowerBasedEmIdentificationLC2::cutLengthCompatibility(const reco::PFCluster& clu,const  edm::PtrVector<HGCRecHit>& rcs )
{
	return (std::abs(lengthCompatibility(clu,rcs))<cutLengthCompatibility_); 
}

