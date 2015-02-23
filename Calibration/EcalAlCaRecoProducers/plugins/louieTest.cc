#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
using namespace std;
using namespace edm;
// **********************************************************************
class louieTest : public edm::EDAnalyzer {
	public:
		explicit louieTest(const edm::ParameterSet&);
		~louieTest();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	private:
		edm::Service<TFileService> fs_;
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		//void initEventStructure();
		TTree *flashggTree;

		float var1;
		float var2;
		float var3;
};
// ******************************************************************************************
// ******************************************************************************************
//
// constants, enums and typedefs
//
//
// static data member definitions
//
//
// constructors and destructor
//
louieTest::louieTest(const edm::ParameterSet& iConfig)
{
}
louieTest::~louieTest()
{
}
void
louieTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	// ********************************************************************************
	// access edm objects
var1 = 1;
var2 =2;
var3 =3;

	flashggTree->Fill();
}
	void
louieTest::beginJob()
{
	flashggTree = fs_->make<TTree>("flashggTree","flashgg tree");
	flashggTree->Branch("var1", &var1, "var1/F");
	flashggTree->Branch("var2", &var2, "var2/F");
	flashggTree->Branch("var3", &var3, "var3/F");
}
	void
louieTest::endJob()
{
}
void
louieTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
typedef louieTest louieTest;
DEFINE_FWK_MODULE(louieTest);

