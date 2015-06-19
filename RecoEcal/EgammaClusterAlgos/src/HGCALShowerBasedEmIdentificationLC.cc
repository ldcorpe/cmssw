#include "RecoEcal/EgammaClusterAlgos/interface/HGCALShowerBasedEmIdentificationLC.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "TPrincipal.h"


HGCALShowerBasedEmIdentificationLC::HGCALShowerBasedEmIdentificationLC (bool withPileup): 
 withPileup_(withPileup)
{
  
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

HGCALShowerBasedEmIdentificationLC::~HGCALShowerBasedEmIdentificationLC ()
{
//  delete pcaShowerAnalysis_;
}

void HGCALShowerBasedEmIdentificationLC::setShowerPosition(const math::XYZPoint &pos)
{
  showerPos_ = pos;
  showerPosIsSet_ = true;
}

void HGCALShowerBasedEmIdentificationLC::setShowerDirection(const math::XYZVector &dir)
{
  showerDir_ = dir;
  showerDirIsSet_ = true;
}

bool HGCALShowerBasedEmIdentificationLC::isEm(const reco::PFCluster& clu,const  edm::PtrVector<reco::PFRecHit>& rcs )
{
  //return (cutStartPosition(clu) &&  cutSigmaetaeta(clu) && cutHadOverEm(clu));
  return (cutStartPosition(clu,rcs) &&  cutSigmaetaeta(clu,rcs) && cutLengthCompatibility(clu,rcs));
} 

math::XYZPoint HGCALShowerBasedEmIdentificationLC::startPosition(const reco::PFCluster& clu,const  edm::PtrVector<reco::PFRecHit>& rcs )
{

//	std::cout << "debug 0.0" << std::endl;
	math::XYZPoint firstPos;
	double zmin = 10000.0;  
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {
//		std::cout << "debug 0.0.0" << std::endl;

		// auto& refhit;// = clu.recHitFractions()[ih].recHitRef();
		//const DetId & id_ =  clu.recHitFractions()[ih].first;
		const DetId & id_ =  clu.hitsAndFractions()[ih].first;
		int matchIndex =-1;
	//	std::cout << "rcs size " << rcs.size() << std::endl;
	//	std::cout << "detId pos " << id_.position() << std::endl;
		for (unsigned int irh=0; irh < rcs.size() ; irh++){

		//		if (irh%1000 ==0) { std::cout << "rechit " << irh << " pos " <<rcs[irh]->position() << std::endl ;}
			if (rcs[irh]->detId() == id_)
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
	  auto refhit  = rcs[matchIndex];

//		std::cout << "debug 0.0.1" << std::endl;
		const auto& pos = refhit->position();
//		std::cout << "debug 0.0.2" << std::endl;
//		std::cout << "debug 0.0.3" << std::endl;
		if (id_.det()==DetId::Forward) {      
	//		std::cout << "debug 0.0.4" << std::endl;
			if (std::abs(pos.z())<zmin) {
				firstPos = pos;
	//			std::cout << "debug 0.0.5" << std::endl;
				zmin = std::abs(pos.z());
		//		std::cout << "debug 0.0.6" << std::endl;
			}
		}
	}

//	std::cout << "debug 0.1" << std::endl;
	// refine the first position estimation, taking the max energy in the first layer 
	double maxfirstenergy=0.; 
	for (unsigned int ih=0;ih<clu.recHitFractions().size();++ih) {
	//	const auto& refhit = clu.recHitFractions()[ih].recHitRef();
	//	const DetId & id_(refhit->detId());
		const DetId & id_ =  clu.hitsAndFractions()[ih].first;
		int matchIndex =-1;
	//	std::cout << "rcs size " << rcs.size() << std::endl;
	//	std::cout << "detId pos " << id_.position() << std::endl;
		for (unsigned int irh=0; irh < rcs.size() ; irh++){

		//		if (irh%1000 ==0) { std::cout << "rechit " << irh << " pos " <<rcs[irh]->position() << std::endl ;}
			if (rcs[irh]->detId() == id_)
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
	  auto refhit  = rcs[matchIndex];
		const auto& pos = refhit->position();
		if (id_.det()==DetId::Forward) {
			if (std::abs(pos.z()) != zmin) continue;
			if (refhit->energy() > maxfirstenergy) {
				firstPos = pos;
				maxfirstenergy = refhit->energy();
			}
		}
	}

//	std::cout << "debug 0.2" << std::endl;
	// finally refine firstPos x and y using the meaured direction 
	if (!showerPosIsSet_ || !showerPosIsSet_) return firstPos;

	double lambda = (firstPos-showerPos_).z()/showerDir_.z();
	math::XYZPoint extraPos = showerPos_ + lambda*showerDir_;	 
	firstPos = extraPos;

	return firstPos;

}

double HGCALShowerBasedEmIdentificationLC::sigmaetaeta(const reco::PFCluster& clu,const  edm::PtrVector<reco::PFRecHit>& rcs )
{

//	std::cout << "debug 0" << std::endl;
	double sigmaetaeta=0., sumnrj=0.;
	math::XYZPoint firstPos = startPosition(clu, rcs);

	std::cout << "debug s_ee 0 " << sigmaetaeta << std::endl;
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
			if (rcs[irh]->detId() == id_)
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
	  auto refhit  = rcs[matchIndex];
		if (id_.det()==DetId::Forward) {
			math::XYZPoint cellPos = refhit->position();	     
			math::XYZVector radius, longitudinal, transverse;
			radius = cellPos - firstPos;
			// distances in local coordinates
			longitudinal =  (radius.Dot(showerDir_))*showerDir_.unit()/showerDir_.R();
			transverse = radius - longitudinal;
			// apply energy cut cut
			if (!withPileup_ || refhit->energy()>minenergy_*mip_) {
				// simple transversal cut, later can refine as function of depth
				if (!withPileup_ || transverse.R() < rmax_) {
					const double deta = (cellPos.eta()-showerPos_.eta());
					sigmaetaeta += deta*deta*refhit->energy();
	std::cout << "debug s_ee 1 " << sigmaetaeta << std::endl;
					sumnrj += refhit->energy();
				}
			}
		}
	}

	//std::cout << "debug 2" << std::endl;
	sigmaetaeta /= sumnrj;
	sigmaetaeta = sqrt(sigmaetaeta);
	std::cout << "debug s_ee 2 " << sigmaetaeta << std::endl;

	// now correct the eta dependency
	double feta;
	constexpr double feta_0 = 0.00964148 - 0.01078431*1.5 + 0.00495703*1.5*1.5;
	const double clu_eta = std::abs(clu.eta()); 
	feta = 0.00964148 - clu_eta*(0.0107843 - 0.00495703*clu_eta);
	sigmaetaeta *= feta_0 / feta ;
	std::cout << "debug s_ee 3 " << sigmaetaeta << std::endl;

//	std::cout << "debug 3" << std::endl;
	return sigmaetaeta;

}

double HGCALShowerBasedEmIdentificationLC::lengthCompatibility(const reco::PFCluster& clu,const  edm::PtrVector<reco::PFRecHit>& rcs  )
{

	// check that showerPos and showerDir have been set
	if (!showerPosIsSet_) {
		std::cout << "[HGCALShowerBasedEmIdentificationLC::lengthCompatibility] error, showwer position not set " << std::endl;
		std::cout << "[HGCALShowerBasedEmIdentificationLC::lengthCompatibility] error, please invoke setShowerPos and setShowerDir before invoking this function " << std::endl;
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

	lengthCompatibility = (predictedLength-length)/predictedSigma;

	return lengthCompatibility;

}

bool HGCALShowerBasedEmIdentificationLC::cutSigmaetaeta(const reco::PFCluster& clu,const  edm::PtrVector<reco::PFRecHit>& rcs )
{
	return (sigmaetaeta(clu, rcs)<cutSigmaetaeta_);
}

bool HGCALShowerBasedEmIdentificationLC::cutStartPosition(const reco::PFCluster& clu,const  edm::PtrVector<reco::PFRecHit>& rcs )
{
	return ( std::abs( startPosition(clu,rcs).z() ) < cutStartPosition_ );
}

bool HGCALShowerBasedEmIdentificationLC::cutLengthCompatibility(const reco::PFCluster& clu,const  edm::PtrVector<reco::PFRecHit>& rcs )
{
	return (std::abs(lengthCompatibility(clu,rcs))<cutLengthCompatibility_); 
}

