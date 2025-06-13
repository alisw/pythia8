#include <Pythia8/Analysis.h>
#include <Pythia8/Basics.h>
#include <Pythia8/Event.h>
#include <Pythia8/Info.h>
#include <Pythia8/MathTools.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ResonanceWidths.h>
#include <functional>
#include <istream>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <string>
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

// Pythia8::SlowJetHook file:Pythia8/Analysis.h line:371
struct PyCallBack_Pythia8_SlowJetHook : public Pythia8::SlowJetHook {
	using Pythia8::SlowJetHook::SlowJetHook;

	bool include(int a0, const class Pythia8::Event & a1, class Pythia8::Vec4 & a2, double & a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SlowJetHook *>(this), "include");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SlowJetHook::include\"");
	}
};

// Pythia8::SlowJet file:Pythia8/Analysis.h line:422
struct PyCallBack_Pythia8_SlowJet : public Pythia8::SlowJet {
	using Pythia8::SlowJet::SlowJet;

	bool doStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SlowJet *>(this), "doStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SlowJet::doStep();
	}
	void findNext() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SlowJet *>(this), "findNext");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SlowJet::findNext();
	}
};

void bind_Pythia8_Analysis_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::SlowJetHook file:Pythia8/Analysis.h line:371
		pybind11::class_<Pythia8::SlowJetHook, std::shared_ptr<Pythia8::SlowJetHook>, PyCallBack_Pythia8_SlowJetHook> cl(M("Pythia8"), "SlowJetHook", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_SlowJetHook(); } ) );
		cl.def("include", (bool (Pythia8::SlowJetHook::*)(int, const class Pythia8::Event &, class Pythia8::Vec4 &, double &)) &Pythia8::SlowJetHook::include, "C++: Pythia8::SlowJetHook::include(int, const class Pythia8::Event &, class Pythia8::Vec4 &, double &) --> bool", pybind11::arg("iSel"), pybind11::arg("event"), pybind11::arg("pSel"), pybind11::arg("mSel"));
		cl.def("assign", (class Pythia8::SlowJetHook & (Pythia8::SlowJetHook::*)(const class Pythia8::SlowJetHook &)) &Pythia8::SlowJetHook::operator=, "C++: Pythia8::SlowJetHook::operator=(const class Pythia8::SlowJetHook &) --> class Pythia8::SlowJetHook &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SingleSlowJet file:Pythia8/Analysis.h line:395
		pybind11::class_<Pythia8::SingleSlowJet, std::shared_ptr<Pythia8::SingleSlowJet>> cl(M("Pythia8"), "SingleSlowJet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SingleSlowJet(); } ), "doc" );
		cl.def( pybind11::init( [](class Pythia8::Vec4 const & a0){ return new Pythia8::SingleSlowJet(a0); } ), "doc" , pybind11::arg("pIn"));
		cl.def( pybind11::init( [](class Pythia8::Vec4 const & a0, double const & a1){ return new Pythia8::SingleSlowJet(a0, a1); } ), "doc" , pybind11::arg("pIn"), pybind11::arg("pT2In"));
		cl.def( pybind11::init( [](class Pythia8::Vec4 const & a0, double const & a1, double const & a2){ return new Pythia8::SingleSlowJet(a0, a1, a2); } ), "doc" , pybind11::arg("pIn"), pybind11::arg("pT2In"), pybind11::arg("yIn"));
		cl.def( pybind11::init( [](class Pythia8::Vec4 const & a0, double const & a1, double const & a2, double const & a3){ return new Pythia8::SingleSlowJet(a0, a1, a2, a3); } ), "doc" , pybind11::arg("pIn"), pybind11::arg("pT2In"), pybind11::arg("yIn"), pybind11::arg("phiIn"));
		cl.def( pybind11::init<class Pythia8::Vec4, double, double, double, int>(), pybind11::arg("pIn"), pybind11::arg("pT2In"), pybind11::arg("yIn"), pybind11::arg("phiIn"), pybind11::arg("idxIn") );

		cl.def( pybind11::init( [](Pythia8::SingleSlowJet const &o){ return new Pythia8::SingleSlowJet(o); } ) );
		cl.def_readwrite("p", &Pythia8::SingleSlowJet::p);
		cl.def_readwrite("pT2", &Pythia8::SingleSlowJet::pT2);
		cl.def_readwrite("y", &Pythia8::SingleSlowJet::y);
		cl.def_readwrite("phi", &Pythia8::SingleSlowJet::phi);
		cl.def_readwrite("mult", &Pythia8::SingleSlowJet::mult);
		cl.def_readwrite("idx", &Pythia8::SingleSlowJet::idx);
		cl.def("assign", (class Pythia8::SingleSlowJet & (Pythia8::SingleSlowJet::*)(const class Pythia8::SingleSlowJet &)) &Pythia8::SingleSlowJet::operator=, "C++: Pythia8::SingleSlowJet::operator=(const class Pythia8::SingleSlowJet &) --> class Pythia8::SingleSlowJet &", pybind11::return_value_policy::reference, pybind11::arg("ssj"));
	}
	{ // Pythia8::SlowJet file:Pythia8/Analysis.h line:422
		pybind11::class_<Pythia8::SlowJet, std::shared_ptr<Pythia8::SlowJet>, PyCallBack_Pythia8_SlowJet> cl(M("Pythia8"), "SlowJet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](int const & a0, double const & a1){ return new Pythia8::SlowJet(a0, a1); }, [](int const & a0, double const & a1){ return new PyCallBack_Pythia8_SlowJet(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2){ return new Pythia8::SlowJet(a0, a1, a2); }, [](int const & a0, double const & a1, double const & a2){ return new PyCallBack_Pythia8_SlowJet(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3){ return new Pythia8::SlowJet(a0, a1, a2, a3); }, [](int const & a0, double const & a1, double const & a2, double const & a3){ return new PyCallBack_Pythia8_SlowJet(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4){ return new Pythia8::SlowJet(a0, a1, a2, a3, a4); }, [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4){ return new PyCallBack_Pythia8_SlowJet(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5){ return new Pythia8::SlowJet(a0, a1, a2, a3, a4, a5); }, [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5){ return new PyCallBack_Pythia8_SlowJet(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5, class Pythia8::SlowJetHook * a6){ return new Pythia8::SlowJet(a0, a1, a2, a3, a4, a5, a6); }, [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5, class Pythia8::SlowJetHook * a6){ return new PyCallBack_Pythia8_SlowJet(a0, a1, a2, a3, a4, a5, a6); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5, class Pythia8::SlowJetHook * a6, bool const & a7){ return new Pythia8::SlowJet(a0, a1, a2, a3, a4, a5, a6, a7); }, [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5, class Pythia8::SlowJetHook * a6, bool const & a7){ return new PyCallBack_Pythia8_SlowJet(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc");
		cl.def( pybind11::init<int, double, double, double, int, int, class Pythia8::SlowJetHook *, bool, bool>(), pybind11::arg("powerIn"), pybind11::arg("Rin"), pybind11::arg("pTjetMinIn"), pybind11::arg("etaMaxIn"), pybind11::arg("selectIn"), pybind11::arg("massSetIn"), pybind11::arg("sjHookPtrIn"), pybind11::arg("useFJcoreIn"), pybind11::arg("useStandardRin") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_SlowJet const &o){ return new PyCallBack_Pythia8_SlowJet(o); } ) );
		cl.def( pybind11::init( [](Pythia8::SlowJet const &o){ return new Pythia8::SlowJet(o); } ) );
		cl.def_readwrite("power", &Pythia8::SlowJet::power);
		cl.def_readwrite("R", &Pythia8::SlowJet::R);
		cl.def_readwrite("pTjetMin", &Pythia8::SlowJet::pTjetMin);
		cl.def_readwrite("etaMax", &Pythia8::SlowJet::etaMax);
		cl.def_readwrite("R2", &Pythia8::SlowJet::R2);
		cl.def_readwrite("pT2jetMin", &Pythia8::SlowJet::pT2jetMin);
		cl.def_readwrite("select", &Pythia8::SlowJet::select);
		cl.def_readwrite("massSet", &Pythia8::SlowJet::massSet);
		cl.def_readwrite("useFJcore", &Pythia8::SlowJet::useFJcore);
		cl.def_readwrite("useStandardR", &Pythia8::SlowJet::useStandardR);
		cl.def_readwrite("isAnti", &Pythia8::SlowJet::isAnti);
		cl.def_readwrite("isKT", &Pythia8::SlowJet::isKT);
		cl.def_readwrite("cutInEta", &Pythia8::SlowJet::cutInEta);
		cl.def_readwrite("chargedOnly", &Pythia8::SlowJet::chargedOnly);
		cl.def_readwrite("visibleOnly", &Pythia8::SlowJet::visibleOnly);
		cl.def_readwrite("modifyMass", &Pythia8::SlowJet::modifyMass);
		cl.def_readwrite("noHook", &Pythia8::SlowJet::noHook);
		cl.def_readwrite("clusters", &Pythia8::SlowJet::clusters);
		cl.def_readwrite("jets", &Pythia8::SlowJet::jets);
		cl.def_readwrite("diB", &Pythia8::SlowJet::diB);
		cl.def_readwrite("dij", &Pythia8::SlowJet::dij);
		cl.def_readwrite("origSize", &Pythia8::SlowJet::origSize);
		cl.def_readwrite("clSize", &Pythia8::SlowJet::clSize);
		cl.def_readwrite("clLast", &Pythia8::SlowJet::clLast);
		cl.def_readwrite("jtSize", &Pythia8::SlowJet::jtSize);
		cl.def_readwrite("iMin", &Pythia8::SlowJet::iMin);
		cl.def_readwrite("jMin", &Pythia8::SlowJet::jMin);
		cl.def_readwrite("dPhi", &Pythia8::SlowJet::dPhi);
		cl.def_readwrite("dijTemp", &Pythia8::SlowJet::dijTemp);
		cl.def_readwrite("dMin", &Pythia8::SlowJet::dMin);
		cl.def("analyze", (bool (Pythia8::SlowJet::*)(const class Pythia8::Event &)) &Pythia8::SlowJet::analyze, "C++: Pythia8::SlowJet::analyze(const class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("setup", (bool (Pythia8::SlowJet::*)(const class Pythia8::Event &)) &Pythia8::SlowJet::setup, "C++: Pythia8::SlowJet::setup(const class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("doStep", (bool (Pythia8::SlowJet::*)()) &Pythia8::SlowJet::doStep, "C++: Pythia8::SlowJet::doStep() --> bool");
		cl.def("doNSteps", (bool (Pythia8::SlowJet::*)(int)) &Pythia8::SlowJet::doNSteps, "C++: Pythia8::SlowJet::doNSteps(int) --> bool", pybind11::arg("nStep"));
		cl.def("stopAtN", (bool (Pythia8::SlowJet::*)(int)) &Pythia8::SlowJet::stopAtN, "C++: Pythia8::SlowJet::stopAtN(int) --> bool", pybind11::arg("nStop"));
		cl.def("sizeOrig", (int (Pythia8::SlowJet::*)() const) &Pythia8::SlowJet::sizeOrig, "C++: Pythia8::SlowJet::sizeOrig() const --> int");
		cl.def("sizeJet", (int (Pythia8::SlowJet::*)() const) &Pythia8::SlowJet::sizeJet, "C++: Pythia8::SlowJet::sizeJet() const --> int");
		cl.def("sizeAll", (int (Pythia8::SlowJet::*)() const) &Pythia8::SlowJet::sizeAll, "C++: Pythia8::SlowJet::sizeAll() const --> int");
		cl.def("pT", (double (Pythia8::SlowJet::*)(int) const) &Pythia8::SlowJet::pT, "C++: Pythia8::SlowJet::pT(int) const --> double", pybind11::arg("i"));
		cl.def("y", (double (Pythia8::SlowJet::*)(int) const) &Pythia8::SlowJet::y, "C++: Pythia8::SlowJet::y(int) const --> double", pybind11::arg("i"));
		cl.def("phi", (double (Pythia8::SlowJet::*)(int) const) &Pythia8::SlowJet::phi, "C++: Pythia8::SlowJet::phi(int) const --> double", pybind11::arg("i"));
		cl.def("p", (class Pythia8::Vec4 (Pythia8::SlowJet::*)(int) const) &Pythia8::SlowJet::p, "C++: Pythia8::SlowJet::p(int) const --> class Pythia8::Vec4", pybind11::arg("i"));
		cl.def("m", (double (Pythia8::SlowJet::*)(int) const) &Pythia8::SlowJet::m, "C++: Pythia8::SlowJet::m(int) const --> double", pybind11::arg("i"));
		cl.def("multiplicity", (int (Pythia8::SlowJet::*)(int) const) &Pythia8::SlowJet::multiplicity, "C++: Pythia8::SlowJet::multiplicity(int) const --> int", pybind11::arg("i"));
		cl.def("iNext", (int (Pythia8::SlowJet::*)() const) &Pythia8::SlowJet::iNext, "C++: Pythia8::SlowJet::iNext() const --> int");
		cl.def("jNext", (int (Pythia8::SlowJet::*)() const) &Pythia8::SlowJet::jNext, "C++: Pythia8::SlowJet::jNext() const --> int");
		cl.def("dNext", (double (Pythia8::SlowJet::*)() const) &Pythia8::SlowJet::dNext, "C++: Pythia8::SlowJet::dNext() const --> double");
		cl.def("list", [](Pythia8::SlowJet const &o) -> void { return o.list(); }, "");
		cl.def("list", (void (Pythia8::SlowJet::*)(bool) const) &Pythia8::SlowJet::list, "C++: Pythia8::SlowJet::list(bool) const --> void", pybind11::arg("listAll"));
		cl.def("constituents", (class std::vector<int, class std::allocator<int> > (Pythia8::SlowJet::*)(int)) &Pythia8::SlowJet::constituents, "C++: Pythia8::SlowJet::constituents(int) --> class std::vector<int, class std::allocator<int> >", pybind11::arg("j"));
		cl.def("clusConstituents", (class std::vector<int, class std::allocator<int> > (Pythia8::SlowJet::*)(int)) &Pythia8::SlowJet::clusConstituents, "C++: Pythia8::SlowJet::clusConstituents(int) --> class std::vector<int, class std::allocator<int> >", pybind11::arg("j"));
		cl.def("jetAssignment", (int (Pythia8::SlowJet::*)(int)) &Pythia8::SlowJet::jetAssignment, "C++: Pythia8::SlowJet::jetAssignment(int) --> int", pybind11::arg("i"));
		cl.def("removeJet", (void (Pythia8::SlowJet::*)(int)) &Pythia8::SlowJet::removeJet, "C++: Pythia8::SlowJet::removeJet(int) --> void", pybind11::arg("i"));
		cl.def("findNext", (void (Pythia8::SlowJet::*)()) &Pythia8::SlowJet::findNext, "C++: Pythia8::SlowJet::findNext() --> void");
		cl.def("clusterFJ", (bool (Pythia8::SlowJet::*)()) &Pythia8::SlowJet::clusterFJ, "C++: Pythia8::SlowJet::clusterFJ() --> bool");
		cl.def("assign", (class Pythia8::SlowJet & (Pythia8::SlowJet::*)(const class Pythia8::SlowJet &)) &Pythia8::SlowJet::operator=, "C++: Pythia8::SlowJet::operator=(const class Pythia8::SlowJet &) --> class Pythia8::SlowJet &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	// Pythia8::gammaReal(double) file:Pythia8/MathTools.h line:19
	M("Pythia8").def("gammaReal", (double (*)(double)) &Pythia8::gammaReal, "C++: Pythia8::gammaReal(double) --> double", pybind11::arg("x"));

	// Pythia8::besselI0(double) file:Pythia8/MathTools.h line:22
	M("Pythia8").def("besselI0", (double (*)(double)) &Pythia8::besselI0, "C++: Pythia8::besselI0(double) --> double", pybind11::arg("x"));

	// Pythia8::besselI1(double) file:Pythia8/MathTools.h line:23
	M("Pythia8").def("besselI1", (double (*)(double)) &Pythia8::besselI1, "C++: Pythia8::besselI1(double) --> double", pybind11::arg("x"));

	// Pythia8::besselK0(double) file:Pythia8/MathTools.h line:24
	M("Pythia8").def("besselK0", (double (*)(double)) &Pythia8::besselK0, "C++: Pythia8::besselK0(double) --> double", pybind11::arg("x"));

	// Pythia8::besselK1(double) file:Pythia8/MathTools.h line:25
	M("Pythia8").def("besselK1", (double (*)(double)) &Pythia8::besselK1, "C++: Pythia8::besselK1(double) --> double", pybind11::arg("x"));

	// Pythia8::integrateGauss(double &, class std::function<double (double)>, double, double, double) file:Pythia8/MathTools.h line:28
	M("Pythia8").def("integrateGauss", [](double & a0, class std::function<double (double)> const & a1, double const & a2, double const & a3) -> bool { return Pythia8::integrateGauss(a0, a1, a2, a3); }, "", pybind11::arg("resultOut"), pybind11::arg("f"), pybind11::arg("xLo"), pybind11::arg("xHi"));
	M("Pythia8").def("integrateGauss", (bool (*)(double &, class std::function<double (double)>, double, double, double)) &Pythia8::integrateGauss, "C++: Pythia8::integrateGauss(double &, class std::function<double (double)>, double, double, double) --> bool", pybind11::arg("resultOut"), pybind11::arg("f"), pybind11::arg("xLo"), pybind11::arg("xHi"), pybind11::arg("tol"));

	// Pythia8::brent(double &, class std::function<double (double)>, double, double, double, double, int) file:Pythia8/MathTools.h line:32
	M("Pythia8").def("brent", [](double & a0, class std::function<double (double)> const & a1, double const & a2, double const & a3, double const & a4) -> bool { return Pythia8::brent(a0, a1, a2, a3, a4); }, "", pybind11::arg("solutionOut"), pybind11::arg("f"), pybind11::arg("target"), pybind11::arg("xLo"), pybind11::arg("xHi"));
	M("Pythia8").def("brent", [](double & a0, class std::function<double (double)> const & a1, double const & a2, double const & a3, double const & a4, double const & a5) -> bool { return Pythia8::brent(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("solutionOut"), pybind11::arg("f"), pybind11::arg("target"), pybind11::arg("xLo"), pybind11::arg("xHi"), pybind11::arg("tol"));
	M("Pythia8").def("brent", (bool (*)(double &, class std::function<double (double)>, double, double, double, double, int)) &Pythia8::brent, "C++: Pythia8::brent(double &, class std::function<double (double)>, double, double, double, double, int) --> bool", pybind11::arg("solutionOut"), pybind11::arg("f"), pybind11::arg("target"), pybind11::arg("xLo"), pybind11::arg("xHi"), pybind11::arg("tol"), pybind11::arg("maxIter"));

	// Pythia8::gramDet(double, double, double, double, double, double) file:Pythia8/MathTools.h line:36
	M("Pythia8").def("gramDet", (double (*)(double, double, double, double, double, double)) &Pythia8::gramDet, "C++: Pythia8::gramDet(double, double, double, double, double, double) --> double", pybind11::arg("s01tilde"), pybind11::arg("s12tilde"), pybind11::arg("s02tilde"), pybind11::arg("m0"), pybind11::arg("m1"), pybind11::arg("m2"));

	// Pythia8::gramDet(class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4) file:Pythia8/MathTools.h line:38
	M("Pythia8").def("gramDet", (double (*)(class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4)) &Pythia8::gramDet, "C++: Pythia8::gramDet(class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4) --> double", pybind11::arg("p0"), pybind11::arg("p1"), pybind11::arg("p2"));

	// Pythia8::Li2(const double, const double, const double) file:Pythia8/MathTools.h line:41
	M("Pythia8").def("Li2", [](const double & a0) -> double { return Pythia8::Li2(a0); }, "", pybind11::arg(""));
	M("Pythia8").def("Li2", [](const double & a0, const double & a1) -> double { return Pythia8::Li2(a0, a1); }, "", pybind11::arg(""), pybind11::arg("kmax"));
	M("Pythia8").def("Li2", (double (*)(const double, const double, const double)) &Pythia8::Li2, "C++: Pythia8::Li2(const double, const double, const double) --> double", pybind11::arg(""), pybind11::arg("kmax"), pybind11::arg("xerr"));

}
