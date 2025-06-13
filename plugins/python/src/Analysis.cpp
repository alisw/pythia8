#include <Pythia8/Analysis.h>
#include <Pythia8/Basics.h>
#include <Pythia8/Event.h>
#include <Pythia8/ParticleData.h>
#include <iterator>
#include <memory>
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

void bind_Pythia8_Analysis(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::Sphericity file:Pythia8/Analysis.h line:27
		pybind11::class_<Pythia8::Sphericity, std::shared_ptr<Pythia8::Sphericity>> cl(M("Pythia8"), "Sphericity", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Sphericity(); } ), "doc" );
		cl.def( pybind11::init( [](double const & a0){ return new Pythia8::Sphericity(a0); } ), "doc" , pybind11::arg("powerIn"));
		cl.def( pybind11::init<double, int>(), pybind11::arg("powerIn"), pybind11::arg("selectIn") );

		cl.def("analyze", (bool (Pythia8::Sphericity::*)(const class Pythia8::Event &)) &Pythia8::Sphericity::analyze, "C++: Pythia8::Sphericity::analyze(const class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("sphericity", (double (Pythia8::Sphericity::*)() const) &Pythia8::Sphericity::sphericity, "C++: Pythia8::Sphericity::sphericity() const --> double");
		cl.def("aplanarity", (double (Pythia8::Sphericity::*)() const) &Pythia8::Sphericity::aplanarity, "C++: Pythia8::Sphericity::aplanarity() const --> double");
		cl.def("eigenValue", (double (Pythia8::Sphericity::*)(int) const) &Pythia8::Sphericity::eigenValue, "C++: Pythia8::Sphericity::eigenValue(int) const --> double", pybind11::arg("i"));
		cl.def("eventAxis", (class Pythia8::Vec4 (Pythia8::Sphericity::*)(int) const) &Pythia8::Sphericity::eventAxis, "C++: Pythia8::Sphericity::eventAxis(int) const --> class Pythia8::Vec4", pybind11::arg("i"));
		cl.def("list", (void (Pythia8::Sphericity::*)() const) &Pythia8::Sphericity::list, "C++: Pythia8::Sphericity::list() const --> void");
		cl.def("nError", (int (Pythia8::Sphericity::*)() const) &Pythia8::Sphericity::nError, "C++: Pythia8::Sphericity::nError() const --> int");
	}
	{ // Pythia8::Thrust file:Pythia8/Analysis.h line:80
		pybind11::class_<Pythia8::Thrust, std::shared_ptr<Pythia8::Thrust>> cl(M("Pythia8"), "Thrust", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Thrust(); } ), "doc" );
		cl.def( pybind11::init<int>(), pybind11::arg("selectIn") );

		cl.def("analyze", (bool (Pythia8::Thrust::*)(const class Pythia8::Event &)) &Pythia8::Thrust::analyze, "C++: Pythia8::Thrust::analyze(const class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("thrust", (double (Pythia8::Thrust::*)() const) &Pythia8::Thrust::thrust, "C++: Pythia8::Thrust::thrust() const --> double");
		cl.def("tMajor", (double (Pythia8::Thrust::*)() const) &Pythia8::Thrust::tMajor, "C++: Pythia8::Thrust::tMajor() const --> double");
		cl.def("tMinor", (double (Pythia8::Thrust::*)() const) &Pythia8::Thrust::tMinor, "C++: Pythia8::Thrust::tMinor() const --> double");
		cl.def("oblateness", (double (Pythia8::Thrust::*)() const) &Pythia8::Thrust::oblateness, "C++: Pythia8::Thrust::oblateness() const --> double");
		cl.def("eventAxis", (class Pythia8::Vec4 (Pythia8::Thrust::*)(int) const) &Pythia8::Thrust::eventAxis, "C++: Pythia8::Thrust::eventAxis(int) const --> class Pythia8::Vec4", pybind11::arg("i"));
		cl.def("list", (void (Pythia8::Thrust::*)() const) &Pythia8::Thrust::list, "C++: Pythia8::Thrust::list() const --> void");
		cl.def("nError", (int (Pythia8::Thrust::*)() const) &Pythia8::Thrust::nError, "C++: Pythia8::Thrust::nError() const --> int");
	}
	{ // Pythia8::SingleClusterJet file:Pythia8/Analysis.h line:128
		pybind11::class_<Pythia8::SingleClusterJet, std::shared_ptr<Pythia8::SingleClusterJet>> cl(M("Pythia8"), "SingleClusterJet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SingleClusterJet(); } ), "doc" );
		cl.def( pybind11::init( [](class Pythia8::Vec4 const & a0){ return new Pythia8::SingleClusterJet(a0); } ), "doc" , pybind11::arg("pJetIn"));
		cl.def( pybind11::init<class Pythia8::Vec4, int>(), pybind11::arg("pJetIn"), pybind11::arg("motherIn") );

		cl.def( pybind11::init( [](Pythia8::SingleClusterJet const &o){ return new Pythia8::SingleClusterJet(o); } ) );
		cl.def_readwrite("pJet", &Pythia8::SingleClusterJet::pJet);
		cl.def_readwrite("mother", &Pythia8::SingleClusterJet::mother);
		cl.def_readwrite("daughter", &Pythia8::SingleClusterJet::daughter);
		cl.def_readwrite("multiplicity", &Pythia8::SingleClusterJet::multiplicity);
		cl.def_readwrite("isAssigned", &Pythia8::SingleClusterJet::isAssigned);
		cl.def_readwrite("pAbs", &Pythia8::SingleClusterJet::pAbs);
		cl.def_readwrite("pTemp", &Pythia8::SingleClusterJet::pTemp);
		cl.def("assign", (class Pythia8::SingleClusterJet & (Pythia8::SingleClusterJet::*)(const class Pythia8::SingleClusterJet &)) &Pythia8::SingleClusterJet::operator=, "C++: Pythia8::SingleClusterJet::operator=(const class Pythia8::SingleClusterJet &) --> class Pythia8::SingleClusterJet &", pybind11::return_value_policy::reference, pybind11::arg("j"));
	}
	// Pythia8::dist2Fun(int, const class Pythia8::SingleClusterJet &, const class Pythia8::SingleClusterJet &) file:Pythia8/Analysis.h line:155
	M("Pythia8").def("dist2Fun", (double (*)(int, const class Pythia8::SingleClusterJet &, const class Pythia8::SingleClusterJet &)) &Pythia8::dist2Fun, "C++: Pythia8::dist2Fun(int, const class Pythia8::SingleClusterJet &, const class Pythia8::SingleClusterJet &) --> double", pybind11::arg("measure"), pybind11::arg("j1"), pybind11::arg("j2"));

	{ // Pythia8::ClusterJet file:Pythia8/Analysis.h line:179
		pybind11::class_<Pythia8::ClusterJet, std::shared_ptr<Pythia8::ClusterJet>> cl(M("Pythia8"), "ClusterJet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ClusterJet(); } ), "doc" );
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0){ return new Pythia8::ClusterJet(a0); } ), "doc" , pybind11::arg("measureIn"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, int const & a1){ return new Pythia8::ClusterJet(a0, a1); } ), "doc" , pybind11::arg("measureIn"), pybind11::arg("selectIn"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, int const & a1, int const & a2){ return new Pythia8::ClusterJet(a0, a1, a2); } ), "doc" , pybind11::arg("measureIn"), pybind11::arg("selectIn"), pybind11::arg("massSetIn"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, int const & a1, int const & a2, bool const & a3){ return new Pythia8::ClusterJet(a0, a1, a2, a3); } ), "doc" , pybind11::arg("measureIn"), pybind11::arg("selectIn"), pybind11::arg("massSetIn"), pybind11::arg("preclusterIn"));
		cl.def( pybind11::init<std::string, int, int, bool, bool>(), pybind11::arg("measureIn"), pybind11::arg("selectIn"), pybind11::arg("massSetIn"), pybind11::arg("preclusterIn"), pybind11::arg("reassignIn") );

		cl.def("analyze", [](Pythia8::ClusterJet &o, const class Pythia8::Event & a0, double const & a1, double const & a2) -> bool { return o.analyze(a0, a1, a2); }, "", pybind11::arg("event"), pybind11::arg("yScaleIn"), pybind11::arg("pTscaleIn"));
		cl.def("analyze", [](Pythia8::ClusterJet &o, const class Pythia8::Event & a0, double const & a1, double const & a2, int const & a3) -> bool { return o.analyze(a0, a1, a2, a3); }, "", pybind11::arg("event"), pybind11::arg("yScaleIn"), pybind11::arg("pTscaleIn"), pybind11::arg("nJetMinIn"));
		cl.def("analyze", (bool (Pythia8::ClusterJet::*)(const class Pythia8::Event &, double, double, int, int)) &Pythia8::ClusterJet::analyze, "C++: Pythia8::ClusterJet::analyze(const class Pythia8::Event &, double, double, int, int) --> bool", pybind11::arg("event"), pybind11::arg("yScaleIn"), pybind11::arg("pTscaleIn"), pybind11::arg("nJetMinIn"), pybind11::arg("nJetMaxIn"));
		cl.def("size", (int (Pythia8::ClusterJet::*)() const) &Pythia8::ClusterJet::size, "C++: Pythia8::ClusterJet::size() const --> int");
		cl.def("p", (class Pythia8::Vec4 (Pythia8::ClusterJet::*)(int) const) &Pythia8::ClusterJet::p, "C++: Pythia8::ClusterJet::p(int) const --> class Pythia8::Vec4", pybind11::arg("i"));
		cl.def("mult", (int (Pythia8::ClusterJet::*)(int) const) &Pythia8::ClusterJet::mult, "C++: Pythia8::ClusterJet::mult(int) const --> int", pybind11::arg("i"));
		cl.def("jetAssignment", (int (Pythia8::ClusterJet::*)(int) const) &Pythia8::ClusterJet::jetAssignment, "C++: Pythia8::ClusterJet::jetAssignment(int) const --> int", pybind11::arg("i"));
		cl.def("list", (void (Pythia8::ClusterJet::*)() const) &Pythia8::ClusterJet::list, "C++: Pythia8::ClusterJet::list() const --> void");
		cl.def("distanceSize", (int (Pythia8::ClusterJet::*)() const) &Pythia8::ClusterJet::distanceSize, "C++: Pythia8::ClusterJet::distanceSize() const --> int");
		cl.def("distance", (double (Pythia8::ClusterJet::*)(int) const) &Pythia8::ClusterJet::distance, "C++: Pythia8::ClusterJet::distance(int) const --> double", pybind11::arg("i"));
		cl.def("nError", (int (Pythia8::ClusterJet::*)() const) &Pythia8::ClusterJet::nError, "C++: Pythia8::ClusterJet::nError() const --> int");
	}
	{ // Pythia8::SingleCell file:Pythia8/Analysis.h line:258
		pybind11::class_<Pythia8::SingleCell, std::shared_ptr<Pythia8::SingleCell>> cl(M("Pythia8"), "SingleCell", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SingleCell(); } ), "doc" );
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::SingleCell(a0); } ), "doc" , pybind11::arg("iCellIn"));
		cl.def( pybind11::init( [](int const & a0, double const & a1){ return new Pythia8::SingleCell(a0, a1); } ), "doc" , pybind11::arg("iCellIn"), pybind11::arg("etaCellIn"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2){ return new Pythia8::SingleCell(a0, a1, a2); } ), "doc" , pybind11::arg("iCellIn"), pybind11::arg("etaCellIn"), pybind11::arg("phiCellIn"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3){ return new Pythia8::SingleCell(a0, a1, a2, a3); } ), "doc" , pybind11::arg("iCellIn"), pybind11::arg("etaCellIn"), pybind11::arg("phiCellIn"), pybind11::arg("eTcellIn"));
		cl.def( pybind11::init<int, double, double, double, int>(), pybind11::arg("iCellIn"), pybind11::arg("etaCellIn"), pybind11::arg("phiCellIn"), pybind11::arg("eTcellIn"), pybind11::arg("multiplicityIn") );

		cl.def_readwrite("iCell", &Pythia8::SingleCell::iCell);
		cl.def_readwrite("etaCell", &Pythia8::SingleCell::etaCell);
		cl.def_readwrite("phiCell", &Pythia8::SingleCell::phiCell);
		cl.def_readwrite("eTcell", &Pythia8::SingleCell::eTcell);
		cl.def_readwrite("multiplicity", &Pythia8::SingleCell::multiplicity);
		cl.def_readwrite("canBeSeed", &Pythia8::SingleCell::canBeSeed);
		cl.def_readwrite("isUsed", &Pythia8::SingleCell::isUsed);
		cl.def_readwrite("isAssigned", &Pythia8::SingleCell::isAssigned);
	}
	{ // Pythia8::SingleCellJet file:Pythia8/Analysis.h line:282
		pybind11::class_<Pythia8::SingleCellJet, std::shared_ptr<Pythia8::SingleCellJet>> cl(M("Pythia8"), "SingleCellJet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SingleCellJet(); } ), "doc" );
		cl.def( pybind11::init( [](double const & a0){ return new Pythia8::SingleCellJet(a0); } ), "doc" , pybind11::arg("eTjetIn"));
		cl.def( pybind11::init( [](double const & a0, double const & a1){ return new Pythia8::SingleCellJet(a0, a1); } ), "doc" , pybind11::arg("eTjetIn"), pybind11::arg("etaCenterIn"));
		cl.def( pybind11::init( [](double const & a0, double const & a1, double const & a2){ return new Pythia8::SingleCellJet(a0, a1, a2); } ), "doc" , pybind11::arg("eTjetIn"), pybind11::arg("etaCenterIn"), pybind11::arg("phiCenterIn"));
		cl.def( pybind11::init( [](double const & a0, double const & a1, double const & a2, double const & a3){ return new Pythia8::SingleCellJet(a0, a1, a2, a3); } ), "doc" , pybind11::arg("eTjetIn"), pybind11::arg("etaCenterIn"), pybind11::arg("phiCenterIn"), pybind11::arg("etaWeightedIn"));
		cl.def( pybind11::init( [](double const & a0, double const & a1, double const & a2, double const & a3, double const & a4){ return new Pythia8::SingleCellJet(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("eTjetIn"), pybind11::arg("etaCenterIn"), pybind11::arg("phiCenterIn"), pybind11::arg("etaWeightedIn"), pybind11::arg("phiWeightedIn"));
		cl.def( pybind11::init( [](double const & a0, double const & a1, double const & a2, double const & a3, double const & a4, int const & a5){ return new Pythia8::SingleCellJet(a0, a1, a2, a3, a4, a5); } ), "doc" , pybind11::arg("eTjetIn"), pybind11::arg("etaCenterIn"), pybind11::arg("phiCenterIn"), pybind11::arg("etaWeightedIn"), pybind11::arg("phiWeightedIn"), pybind11::arg("multiplicityIn"));
		cl.def( pybind11::init<double, double, double, double, double, int, class Pythia8::Vec4>(), pybind11::arg("eTjetIn"), pybind11::arg("etaCenterIn"), pybind11::arg("phiCenterIn"), pybind11::arg("etaWeightedIn"), pybind11::arg("phiWeightedIn"), pybind11::arg("multiplicityIn"), pybind11::arg("pMassiveIn") );

		cl.def_readwrite("eTjet", &Pythia8::SingleCellJet::eTjet);
		cl.def_readwrite("etaCenter", &Pythia8::SingleCellJet::etaCenter);
		cl.def_readwrite("phiCenter", &Pythia8::SingleCellJet::phiCenter);
		cl.def_readwrite("etaWeighted", &Pythia8::SingleCellJet::etaWeighted);
		cl.def_readwrite("phiWeighted", &Pythia8::SingleCellJet::phiWeighted);
		cl.def_readwrite("multiplicity", &Pythia8::SingleCellJet::multiplicity);
		cl.def_readwrite("pMassive", &Pythia8::SingleCellJet::pMassive);
	}
	{ // Pythia8::CellJet file:Pythia8/Analysis.h line:307
		pybind11::class_<Pythia8::CellJet, std::shared_ptr<Pythia8::CellJet>> cl(M("Pythia8"), "CellJet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::CellJet(); } ), "doc" );
		cl.def( pybind11::init( [](double const & a0){ return new Pythia8::CellJet(a0); } ), "doc" , pybind11::arg("etaMaxIn"));
		cl.def( pybind11::init( [](double const & a0, int const & a1){ return new Pythia8::CellJet(a0, a1); } ), "doc" , pybind11::arg("etaMaxIn"), pybind11::arg("nEtaIn"));
		cl.def( pybind11::init( [](double const & a0, int const & a1, int const & a2){ return new Pythia8::CellJet(a0, a1, a2); } ), "doc" , pybind11::arg("etaMaxIn"), pybind11::arg("nEtaIn"), pybind11::arg("nPhiIn"));
		cl.def( pybind11::init( [](double const & a0, int const & a1, int const & a2, int const & a3){ return new Pythia8::CellJet(a0, a1, a2, a3); } ), "doc" , pybind11::arg("etaMaxIn"), pybind11::arg("nEtaIn"), pybind11::arg("nPhiIn"), pybind11::arg("selectIn"));
		cl.def( pybind11::init( [](double const & a0, int const & a1, int const & a2, int const & a3, int const & a4){ return new Pythia8::CellJet(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("etaMaxIn"), pybind11::arg("nEtaIn"), pybind11::arg("nPhiIn"), pybind11::arg("selectIn"), pybind11::arg("smearIn"));
		cl.def( pybind11::init( [](double const & a0, int const & a1, int const & a2, int const & a3, int const & a4, double const & a5){ return new Pythia8::CellJet(a0, a1, a2, a3, a4, a5); } ), "doc" , pybind11::arg("etaMaxIn"), pybind11::arg("nEtaIn"), pybind11::arg("nPhiIn"), pybind11::arg("selectIn"), pybind11::arg("smearIn"), pybind11::arg("resolutionIn"));
		cl.def( pybind11::init( [](double const & a0, int const & a1, int const & a2, int const & a3, int const & a4, double const & a5, double const & a6){ return new Pythia8::CellJet(a0, a1, a2, a3, a4, a5, a6); } ), "doc" , pybind11::arg("etaMaxIn"), pybind11::arg("nEtaIn"), pybind11::arg("nPhiIn"), pybind11::arg("selectIn"), pybind11::arg("smearIn"), pybind11::arg("resolutionIn"), pybind11::arg("upperCutIn"));
		cl.def( pybind11::init( [](double const & a0, int const & a1, int const & a2, int const & a3, int const & a4, double const & a5, double const & a6, double const & a7){ return new Pythia8::CellJet(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc" , pybind11::arg("etaMaxIn"), pybind11::arg("nEtaIn"), pybind11::arg("nPhiIn"), pybind11::arg("selectIn"), pybind11::arg("smearIn"), pybind11::arg("resolutionIn"), pybind11::arg("upperCutIn"), pybind11::arg("thresholdIn"));
		cl.def( pybind11::init<double, int, int, int, int, double, double, double, class Pythia8::Rndm *>(), pybind11::arg("etaMaxIn"), pybind11::arg("nEtaIn"), pybind11::arg("nPhiIn"), pybind11::arg("selectIn"), pybind11::arg("smearIn"), pybind11::arg("resolutionIn"), pybind11::arg("upperCutIn"), pybind11::arg("thresholdIn"), pybind11::arg("rndmPtrIn") );

		cl.def( pybind11::init( [](Pythia8::CellJet const &o){ return new Pythia8::CellJet(o); } ) );
		cl.def("analyze", [](Pythia8::CellJet &o, const class Pythia8::Event & a0) -> bool { return o.analyze(a0); }, "", pybind11::arg("event"));
		cl.def("analyze", [](Pythia8::CellJet &o, const class Pythia8::Event & a0, double const & a1) -> bool { return o.analyze(a0, a1); }, "", pybind11::arg("event"), pybind11::arg("eTjetMinIn"));
		cl.def("analyze", [](Pythia8::CellJet &o, const class Pythia8::Event & a0, double const & a1, double const & a2) -> bool { return o.analyze(a0, a1, a2); }, "", pybind11::arg("event"), pybind11::arg("eTjetMinIn"), pybind11::arg("coneRadiusIn"));
		cl.def("analyze", (bool (Pythia8::CellJet::*)(const class Pythia8::Event &, double, double, double)) &Pythia8::CellJet::analyze, "C++: Pythia8::CellJet::analyze(const class Pythia8::Event &, double, double, double) --> bool", pybind11::arg("event"), pybind11::arg("eTjetMinIn"), pybind11::arg("coneRadiusIn"), pybind11::arg("eTseedIn"));
		cl.def("size", (int (Pythia8::CellJet::*)() const) &Pythia8::CellJet::size, "C++: Pythia8::CellJet::size() const --> int");
		cl.def("eT", (double (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::eT, "C++: Pythia8::CellJet::eT(int) const --> double", pybind11::arg("i"));
		cl.def("etaCenter", (double (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::etaCenter, "C++: Pythia8::CellJet::etaCenter(int) const --> double", pybind11::arg("i"));
		cl.def("phiCenter", (double (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::phiCenter, "C++: Pythia8::CellJet::phiCenter(int) const --> double", pybind11::arg("i"));
		cl.def("etaWeighted", (double (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::etaWeighted, "C++: Pythia8::CellJet::etaWeighted(int) const --> double", pybind11::arg("i"));
		cl.def("phiWeighted", (double (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::phiWeighted, "C++: Pythia8::CellJet::phiWeighted(int) const --> double", pybind11::arg("i"));
		cl.def("multiplicity", (int (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::multiplicity, "C++: Pythia8::CellJet::multiplicity(int) const --> int", pybind11::arg("i"));
		cl.def("pMassless", (class Pythia8::Vec4 (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::pMassless, "C++: Pythia8::CellJet::pMassless(int) const --> class Pythia8::Vec4", pybind11::arg("i"));
		cl.def("pMassive", (class Pythia8::Vec4 (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::pMassive, "C++: Pythia8::CellJet::pMassive(int) const --> class Pythia8::Vec4", pybind11::arg("i"));
		cl.def("m", (double (Pythia8::CellJet::*)(int) const) &Pythia8::CellJet::m, "C++: Pythia8::CellJet::m(int) const --> double", pybind11::arg("i"));
		cl.def("list", (void (Pythia8::CellJet::*)() const) &Pythia8::CellJet::list, "C++: Pythia8::CellJet::list() const --> void");
		cl.def("nError", (int (Pythia8::CellJet::*)() const) &Pythia8::CellJet::nError, "C++: Pythia8::CellJet::nError() const --> int");
	}
}
