#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/History.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonLevel.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/RHadrons.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/SimpleWeakShowerMEs.h>
#include <Pythia8/SpaceShower.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/TimeShower.h>
#include <Pythia8/Weights.h>
#include <functional>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <functional>
#include <string>
#include <Pythia8/UserHooks.h>
#include <Pythia8/SplittingsOnia.h>
#include <Pythia8/HeavyIons.h>
#include <Pythia8/BeamShape.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*);
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>);
#endif

void bind_Pythia8_SimpleWeakShowerMEs(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::SimpleWeakShowerMEs file:Pythia8/SimpleWeakShowerMEs.h line:26
		pybind11::class_<Pythia8::SimpleWeakShowerMEs, std::shared_ptr<Pythia8::SimpleWeakShowerMEs>> cl(M("Pythia8"), "SimpleWeakShowerMEs", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SimpleWeakShowerMEs(); } ) );
		cl.def( pybind11::init( [](Pythia8::SimpleWeakShowerMEs const &o){ return new Pythia8::SimpleWeakShowerMEs(o); } ) );
		cl.def("getMEqg2qg", (double (Pythia8::SimpleWeakShowerMEs::*)(double, double, double)) &Pythia8::SimpleWeakShowerMEs::getMEqg2qg, "C++: Pythia8::SimpleWeakShowerMEs::getMEqg2qg(double, double, double) --> double", pybind11::arg("sH"), pybind11::arg("tH"), pybind11::arg("uH"));
		cl.def("getMEqq2qq", (double (Pythia8::SimpleWeakShowerMEs::*)(double, double, double, bool)) &Pythia8::SimpleWeakShowerMEs::getMEqq2qq, "C++: Pythia8::SimpleWeakShowerMEs::getMEqq2qq(double, double, double, bool) --> double", pybind11::arg("sH"), pybind11::arg("tH"), pybind11::arg("uH"), pybind11::arg("sameID"));
		cl.def("getMEgg2gg", (double (Pythia8::SimpleWeakShowerMEs::*)(double, double, double)) &Pythia8::SimpleWeakShowerMEs::getMEgg2gg, "C++: Pythia8::SimpleWeakShowerMEs::getMEgg2gg(double, double, double) --> double", pybind11::arg("sH"), pybind11::arg("tH"), pybind11::arg("uH"));
		cl.def("getMEgg2qqbar", (double (Pythia8::SimpleWeakShowerMEs::*)(double, double, double)) &Pythia8::SimpleWeakShowerMEs::getMEgg2qqbar, "C++: Pythia8::SimpleWeakShowerMEs::getMEgg2qqbar(double, double, double) --> double", pybind11::arg("sH"), pybind11::arg("tH"), pybind11::arg("uH"));
		cl.def("getMEqqbar2gg", (double (Pythia8::SimpleWeakShowerMEs::*)(double, double, double)) &Pythia8::SimpleWeakShowerMEs::getMEqqbar2gg, "C++: Pythia8::SimpleWeakShowerMEs::getMEqqbar2gg(double, double, double) --> double", pybind11::arg("sH"), pybind11::arg("tH"), pybind11::arg("uH"));
		cl.def("getMEqqbar2qqbar", (double (Pythia8::SimpleWeakShowerMEs::*)(double, double, double, bool)) &Pythia8::SimpleWeakShowerMEs::getMEqqbar2qqbar, "C++: Pythia8::SimpleWeakShowerMEs::getMEqqbar2qqbar(double, double, double, bool) --> double", pybind11::arg("sH"), pybind11::arg("tH"), pybind11::arg("uH"), pybind11::arg("sameID"));
		cl.def("getMEqg2qgZ", (double (Pythia8::SimpleWeakShowerMEs::*)(class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4)) &Pythia8::SimpleWeakShowerMEs::getMEqg2qgZ, "C++: Pythia8::SimpleWeakShowerMEs::getMEqg2qgZ(class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4) --> double", pybind11::arg("p1"), pybind11::arg("p2"), pybind11::arg("p3"), pybind11::arg("p4"), pybind11::arg("p5"));
		cl.def("getMEqq2qqZ", (double (Pythia8::SimpleWeakShowerMEs::*)(class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4)) &Pythia8::SimpleWeakShowerMEs::getMEqq2qqZ, "C++: Pythia8::SimpleWeakShowerMEs::getMEqq2qqZ(class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4) --> double", pybind11::arg("p1"), pybind11::arg("p2"), pybind11::arg("p3"), pybind11::arg("p4"), pybind11::arg("p5"));
		cl.def("assign", (class Pythia8::SimpleWeakShowerMEs & (Pythia8::SimpleWeakShowerMEs::*)(const class Pythia8::SimpleWeakShowerMEs &)) &Pythia8::SimpleWeakShowerMEs::operator=, "C++: Pythia8::SimpleWeakShowerMEs::operator=(const class Pythia8::SimpleWeakShowerMEs &) --> class Pythia8::SimpleWeakShowerMEs &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::Clustering file:Pythia8/History.h line:32
		pybind11::class_<Pythia8::Clustering, std::shared_ptr<Pythia8::Clustering>> cl(M("Pythia8"), "Clustering", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Clustering(); } ) );
		cl.def( pybind11::init( [](Pythia8::Clustering const &o){ return new Pythia8::Clustering(o); } ) );
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4){ return new Pythia8::Clustering(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5){ return new Pythia8::Clustering(a0, a1, a2, a3, a4, a5); } ), "doc" , pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"), pybind11::arg("flavRadBefIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6){ return new Pythia8::Clustering(a0, a1, a2, a3, a4, a5, a6); } ), "doc" , pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"), pybind11::arg("flavRadBefIn"), pybind11::arg("spinRadIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7){ return new Pythia8::Clustering(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc" , pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"), pybind11::arg("flavRadBefIn"), pybind11::arg("spinRadIn"), pybind11::arg("spinEmtIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7, int const & a8){ return new Pythia8::Clustering(a0, a1, a2, a3, a4, a5, a6, a7, a8); } ), "doc" , pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"), pybind11::arg("flavRadBefIn"), pybind11::arg("spinRadIn"), pybind11::arg("spinEmtIn"), pybind11::arg("spinRecIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9){ return new Pythia8::Clustering(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); } ), "doc" , pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"), pybind11::arg("flavRadBefIn"), pybind11::arg("spinRadIn"), pybind11::arg("spinEmtIn"), pybind11::arg("spinRecIn"), pybind11::arg("spinRadBefIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10){ return new Pythia8::Clustering(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); } ), "doc" , pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"), pybind11::arg("flavRadBefIn"), pybind11::arg("spinRadIn"), pybind11::arg("spinEmtIn"), pybind11::arg("spinRecIn"), pybind11::arg("spinRadBefIn"), pybind11::arg("radBefIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11){ return new Pythia8::Clustering(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); } ), "doc" , pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"), pybind11::arg("flavRadBefIn"), pybind11::arg("spinRadIn"), pybind11::arg("spinEmtIn"), pybind11::arg("spinRecIn"), pybind11::arg("spinRadBefIn"), pybind11::arg("radBefIn"), pybind11::arg("recBefIn"));
		cl.def( pybind11::init<int, int, int, int, double, int, int, int, int, int, int, int, class std::map<int, int, struct std::less<int>, class std::allocator<struct std::pair<const int, int> > >>(), pybind11::arg("emtIn"), pybind11::arg("radIn"), pybind11::arg("recIn"), pybind11::arg("partnerIn"), pybind11::arg("pTscaleIn"), pybind11::arg("flavRadBefIn"), pybind11::arg("spinRadIn"), pybind11::arg("spinEmtIn"), pybind11::arg("spinRecIn"), pybind11::arg("spinRadBefIn"), pybind11::arg("radBefIn"), pybind11::arg("recBefIn"), pybind11::arg("posIn") );

		cl.def_readwrite("emitted", &Pythia8::Clustering::emitted);
		cl.def_readwrite("emittor", &Pythia8::Clustering::emittor);
		cl.def_readwrite("recoiler", &Pythia8::Clustering::recoiler);
		cl.def_readwrite("partner", &Pythia8::Clustering::partner);
		cl.def_readwrite("pTscale", &Pythia8::Clustering::pTscale);
		cl.def_readwrite("flavRadBef", &Pythia8::Clustering::flavRadBef);
		cl.def_readwrite("spinRad", &Pythia8::Clustering::spinRad);
		cl.def_readwrite("spinEmt", &Pythia8::Clustering::spinEmt);
		cl.def_readwrite("spinRec", &Pythia8::Clustering::spinRec);
		cl.def_readwrite("spinRadBef", &Pythia8::Clustering::spinRadBef);
		cl.def_readwrite("radBef", &Pythia8::Clustering::radBef);
		cl.def_readwrite("recBef", &Pythia8::Clustering::recBef);
		cl.def_readwrite("iPosInMother", &Pythia8::Clustering::iPosInMother);
		cl.def("pT", (double (Pythia8::Clustering::*)() const) &Pythia8::Clustering::pT, "C++: Pythia8::Clustering::pT() const --> double");
		cl.def("list", (void (Pythia8::Clustering::*)() const) &Pythia8::Clustering::list, "C++: Pythia8::Clustering::list() const --> void");
	}
	{ // Pythia8::History file:Pythia8/History.h line:116
		pybind11::class_<Pythia8::History, std::shared_ptr<Pythia8::History>> cl(M("Pythia8"), "History", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int, double, class Pythia8::Event, class Pythia8::Clustering, class std::shared_ptr<class Pythia8::MergingHooks>, class Pythia8::BeamParticle, class Pythia8::BeamParticle, class Pythia8::ParticleData *, class Pythia8::Info *, class Pythia8::PartonLevel *, class Pythia8::CoupSM *, bool, bool, bool, bool, double, class Pythia8::History *>(), pybind11::arg("depthIn"), pybind11::arg("scalein"), pybind11::arg("statein"), pybind11::arg("c"), pybind11::arg("mergingHooksPtrIn"), pybind11::arg("beamAIn"), pybind11::arg("beamBIn"), pybind11::arg("particleDataPtrIn"), pybind11::arg("infoPtrIn"), pybind11::arg("showersIn"), pybind11::arg("coupSMPtrIn"), pybind11::arg("isOrdered"), pybind11::arg("isStronglyOrdered"), pybind11::arg("isAllowed"), pybind11::arg("isNextInInput"), pybind11::arg("probin"), pybind11::arg("mothin") );

		cl.def("projectOntoDesiredHistories", (bool (Pythia8::History::*)()) &Pythia8::History::projectOntoDesiredHistories, "C++: Pythia8::History::projectOntoDesiredHistories() --> bool");
		cl.def("weightCKKWL", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double)) &Pythia8::History::weightCKKWL, "C++: Pythia8::History::weightCKKWL(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"));
		cl.def("weightNL3Loop", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, double)) &Pythia8::History::weightNL3Loop, "C++: Pythia8::History::weightNL3Loop(class Pythia8::PartonLevel *, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("RN"));
		cl.def("weightNL3First", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, class Pythia8::Rndm *)) &Pythia8::History::weightNL3First, "C++: Pythia8::History::weightNL3First(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, class Pythia8::Rndm *) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"), pybind11::arg("rndmPtr"));
		cl.def("weightNL3Tree", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double)) &Pythia8::History::weightNL3Tree, "C++: Pythia8::History::weightNL3Tree(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"));
		cl.def("weightUMEPSTree", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double)) &Pythia8::History::weightUMEPSTree, "C++: Pythia8::History::weightUMEPSTree(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"));
		cl.def("weightUMEPSSubt", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double)) &Pythia8::History::weightUMEPSSubt, "C++: Pythia8::History::weightUMEPSSubt(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"));
		cl.def("weightUNLOPSTree", [](Pythia8::History &o, class Pythia8::PartonLevel * a0, class Pythia8::AlphaStrong * a1, class Pythia8::AlphaStrong * a2, class Pythia8::AlphaEM * a3, class Pythia8::AlphaEM * a4, double const & a5) -> std::vector<double, class std::allocator<double> > { return o.weightUNLOPSTree(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"));
		cl.def("weightUNLOPSTree", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, int)) &Pythia8::History::weightUNLOPSTree, "C++: Pythia8::History::weightUNLOPSTree(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, int) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"), pybind11::arg("depthIn"));
		cl.def("weightUNLOPSSubt", [](Pythia8::History &o, class Pythia8::PartonLevel * a0, class Pythia8::AlphaStrong * a1, class Pythia8::AlphaStrong * a2, class Pythia8::AlphaEM * a3, class Pythia8::AlphaEM * a4, double const & a5) -> std::vector<double, class std::allocator<double> > { return o.weightUNLOPSSubt(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"));
		cl.def("weightUNLOPSSubt", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, int)) &Pythia8::History::weightUNLOPSSubt, "C++: Pythia8::History::weightUNLOPSSubt(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, int) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"), pybind11::arg("depthIn"));
		cl.def("weightUNLOPSLoop", [](Pythia8::History &o, class Pythia8::PartonLevel * a0, class Pythia8::AlphaStrong * a1, class Pythia8::AlphaStrong * a2, class Pythia8::AlphaEM * a3, class Pythia8::AlphaEM * a4, double const & a5) -> std::vector<double, class std::allocator<double> > { return o.weightUNLOPSLoop(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"));
		cl.def("weightUNLOPSLoop", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, int)) &Pythia8::History::weightUNLOPSLoop, "C++: Pythia8::History::weightUNLOPSLoop(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, int) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"), pybind11::arg("depthIn"));
		cl.def("weightUNLOPSSubtNLO", [](Pythia8::History &o, class Pythia8::PartonLevel * a0, class Pythia8::AlphaStrong * a1, class Pythia8::AlphaStrong * a2, class Pythia8::AlphaEM * a3, class Pythia8::AlphaEM * a4, double const & a5) -> std::vector<double, class std::allocator<double> > { return o.weightUNLOPSSubtNLO(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"));
		cl.def("weightUNLOPSSubtNLO", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, int)) &Pythia8::History::weightUNLOPSSubtNLO, "C++: Pythia8::History::weightUNLOPSSubtNLO(class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, int) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"), pybind11::arg("depthIn"));
		cl.def("weightUNLOPSFirst", (class std::vector<double, class std::allocator<double> > (Pythia8::History::*)(int, class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, class Pythia8::Rndm *)) &Pythia8::History::weightUNLOPSFirst, "C++: Pythia8::History::weightUNLOPSFirst(int, class Pythia8::PartonLevel *, class Pythia8::AlphaStrong *, class Pythia8::AlphaStrong *, class Pythia8::AlphaEM *, class Pythia8::AlphaEM *, double, class Pythia8::Rndm *) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("order"), pybind11::arg("trial"), pybind11::arg("asFSR"), pybind11::arg("asISR"), pybind11::arg("aemFSR"), pybind11::arg("aemISR"), pybind11::arg("RN"), pybind11::arg("rndmPtr"));
		cl.def("foundAllowedHistories", (bool (Pythia8::History::*)()) &Pythia8::History::foundAllowedHistories, "C++: Pythia8::History::foundAllowedHistories() --> bool");
		cl.def("foundOrderedHistories", (bool (Pythia8::History::*)()) &Pythia8::History::foundOrderedHistories, "C++: Pythia8::History::foundOrderedHistories() --> bool");
		cl.def("foundCompleteHistories", (bool (Pythia8::History::*)()) &Pythia8::History::foundCompleteHistories, "C++: Pythia8::History::foundCompleteHistories() --> bool");
		cl.def("getStartingConditions", (void (Pythia8::History::*)(const double, class Pythia8::Event &)) &Pythia8::History::getStartingConditions, "C++: Pythia8::History::getStartingConditions(const double, class Pythia8::Event &) --> void", pybind11::arg("RN"), pybind11::arg("outState"));
		cl.def("getClusteredEvent", (bool (Pythia8::History::*)(const double, int, class Pythia8::Event &)) &Pythia8::History::getClusteredEvent, "C++: Pythia8::History::getClusteredEvent(const double, int, class Pythia8::Event &) --> bool", pybind11::arg("RN"), pybind11::arg("nSteps"), pybind11::arg("outState"));
		cl.def("getFirstClusteredEventAboveTMS", [](Pythia8::History &o, const double & a0, int const & a1, class Pythia8::Event & a2, int & a3) -> bool { return o.getFirstClusteredEventAboveTMS(a0, a1, a2, a3); }, "", pybind11::arg("RN"), pybind11::arg("nDesired"), pybind11::arg("process"), pybind11::arg("nPerformed"));
		cl.def("getFirstClusteredEventAboveTMS", (bool (Pythia8::History::*)(const double, int, class Pythia8::Event &, int &, bool)) &Pythia8::History::getFirstClusteredEventAboveTMS, "C++: Pythia8::History::getFirstClusteredEventAboveTMS(const double, int, class Pythia8::Event &, int &, bool) --> bool", pybind11::arg("RN"), pybind11::arg("nDesired"), pybind11::arg("process"), pybind11::arg("nPerformed"), pybind11::arg("updateProcess"));
		cl.def("nClusterings", (int (Pythia8::History::*)()) &Pythia8::History::nClusterings, "C++: Pythia8::History::nClusterings() --> int");
		cl.def("lowestMultProc", (class Pythia8::Event (Pythia8::History::*)(const double)) &Pythia8::History::lowestMultProc, "C++: Pythia8::History::lowestMultProc(const double) --> class Pythia8::Event", pybind11::arg("RN"));
		cl.def("getPDFratio", (double (Pythia8::History::*)(int, bool, bool, int, double, double, int, double, double)) &Pythia8::History::getPDFratio, "C++: Pythia8::History::getPDFratio(int, bool, bool, int, double, double, int, double, double) --> double", pybind11::arg("side"), pybind11::arg("forSudakov"), pybind11::arg("useHardPDF"), pybind11::arg("flavNum"), pybind11::arg("xNum"), pybind11::arg("muNum"), pybind11::arg("flavDen"), pybind11::arg("xDen"), pybind11::arg("muDen"));
		cl.def("getWeakProb", (double (Pythia8::History::*)()) &Pythia8::History::getWeakProb, "C++: Pythia8::History::getWeakProb() --> double");
		cl.def("getWeakProb", (double (Pythia8::History::*)(class std::vector<int, class std::allocator<int> > &, class std::vector<class Pythia8::Vec4, class std::allocator<class Pythia8::Vec4> > &, class std::vector<int, class std::allocator<int> >)) &Pythia8::History::getWeakProb, "C++: Pythia8::History::getWeakProb(class std::vector<int, class std::allocator<int> > &, class std::vector<class Pythia8::Vec4, class std::allocator<class Pythia8::Vec4> > &, class std::vector<int, class std::allocator<int> >) --> double", pybind11::arg("mode"), pybind11::arg("mom"), pybind11::arg("fermionLines"));
		cl.def("getSingleWeakProb", (double (Pythia8::History::*)(class std::vector<int, class std::allocator<int> > &, class std::vector<class Pythia8::Vec4, class std::allocator<class Pythia8::Vec4> > &, class std::vector<int, class std::allocator<int> >)) &Pythia8::History::getSingleWeakProb, "C++: Pythia8::History::getSingleWeakProb(class std::vector<int, class std::allocator<int> > &, class std::vector<class Pythia8::Vec4, class std::allocator<class Pythia8::Vec4> > &, class std::vector<int, class std::allocator<int> >) --> double", pybind11::arg("mode"), pybind11::arg("mom"), pybind11::arg("fermionLines"));
		cl.def("findStateTransfer", (void (Pythia8::History::*)(class std::map<int, int, struct std::less<int>, class std::allocator<struct std::pair<const int, int> > > &)) &Pythia8::History::findStateTransfer, "C++: Pythia8::History::findStateTransfer(class std::map<int, int, struct std::less<int>, class std::allocator<struct std::pair<const int, int> > > &) --> void", pybind11::arg("transfer"));
		cl.def("printHistory", (void (Pythia8::History::*)(const double)) &Pythia8::History::printHistory, "C++: Pythia8::History::printHistory(const double) --> void", pybind11::arg("RN"));
		cl.def("printStates", (void (Pythia8::History::*)()) &Pythia8::History::printStates, "C++: Pythia8::History::printStates() --> void");
	}
}
