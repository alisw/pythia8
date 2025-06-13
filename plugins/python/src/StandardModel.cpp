#include <Pythia8/Basics.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/Settings.h>
#include <Pythia8/StandardModel.h>
#include <functional>
#include <istream>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
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

// Pythia8::AlphaStrong file:Pythia8/StandardModel.h line:23
struct PyCallBack_Pythia8_AlphaStrong : public Pythia8::AlphaStrong {
	using Pythia8::AlphaStrong::AlphaStrong;

	void init(double a0, int a1, int a2, bool a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlphaStrong *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AlphaStrong::init(a0, a1, a2, a3);
	}
	void setThresholds(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlphaStrong *>(this), "setThresholds");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AlphaStrong::setThresholds(a0, a1, a2);
	}
};

// Pythia8::AlphaSUN file:Pythia8/StandardModel.h line:226
struct PyCallBack_Pythia8_AlphaSUN : public Pythia8::AlphaSUN {
	using Pythia8::AlphaSUN::AlphaSUN;

	void initAlpha(int a0, int a1, int a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlphaSUN *>(this), "initAlpha");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AlphaSUN::initAlpha(a0, a1, a2, a3, a4);
	}
	void initLambda(int a0, int a1, int a2, double a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlphaSUN *>(this), "initLambda");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AlphaSUN::initLambda(a0, a1, a2, a3);
	}
};

void bind_Pythia8_StandardModel(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::AlphaStrong file:Pythia8/StandardModel.h line:23
		pybind11::class_<Pythia8::AlphaStrong, std::shared_ptr<Pythia8::AlphaStrong>, PyCallBack_Pythia8_AlphaStrong> cl(M("Pythia8"), "AlphaStrong", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::AlphaStrong(); }, [](){ return new PyCallBack_Pythia8_AlphaStrong(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_AlphaStrong const &o){ return new PyCallBack_Pythia8_AlphaStrong(o); } ) );
		cl.def( pybind11::init( [](Pythia8::AlphaStrong const &o){ return new Pythia8::AlphaStrong(o); } ) );
		cl.def_readwrite("isInit", &Pythia8::AlphaStrong::isInit);
		cl.def_readwrite("order", &Pythia8::AlphaStrong::order);
		cl.def_readwrite("nfmax", &Pythia8::AlphaStrong::nfmax);
		cl.def_readwrite("Lambda3Save", &Pythia8::AlphaStrong::Lambda3Save);
		cl.def_readwrite("Lambda4Save", &Pythia8::AlphaStrong::Lambda4Save);
		cl.def_readwrite("Lambda5Save", &Pythia8::AlphaStrong::Lambda5Save);
		cl.def_readwrite("Lambda6Save", &Pythia8::AlphaStrong::Lambda6Save);
		cl.def_readwrite("Lambda3Save2", &Pythia8::AlphaStrong::Lambda3Save2);
		cl.def_readwrite("Lambda4Save2", &Pythia8::AlphaStrong::Lambda4Save2);
		cl.def_readwrite("Lambda5Save2", &Pythia8::AlphaStrong::Lambda5Save2);
		cl.def_readwrite("Lambda6Save2", &Pythia8::AlphaStrong::Lambda6Save2);
		cl.def_readwrite("scale2Min", &Pythia8::AlphaStrong::scale2Min);
		cl.def_readwrite("mc", &Pythia8::AlphaStrong::mc);
		cl.def_readwrite("mb", &Pythia8::AlphaStrong::mb);
		cl.def_readwrite("mt", &Pythia8::AlphaStrong::mt);
		cl.def_readwrite("mc2", &Pythia8::AlphaStrong::mc2);
		cl.def_readwrite("mb2", &Pythia8::AlphaStrong::mb2);
		cl.def_readwrite("mt2", &Pythia8::AlphaStrong::mt2);
		cl.def_readwrite("useCMW", &Pythia8::AlphaStrong::useCMW);
		cl.def("init", [](Pythia8::AlphaStrong &o) -> void { return o.init(); }, "");
		cl.def("init", [](Pythia8::AlphaStrong &o, double const & a0) -> void { return o.init(a0); }, "", pybind11::arg("valueIn"));
		cl.def("init", [](Pythia8::AlphaStrong &o, double const & a0, int const & a1) -> void { return o.init(a0, a1); }, "", pybind11::arg("valueIn"), pybind11::arg("orderIn"));
		cl.def("init", [](Pythia8::AlphaStrong &o, double const & a0, int const & a1, int const & a2) -> void { return o.init(a0, a1, a2); }, "", pybind11::arg("valueIn"), pybind11::arg("orderIn"), pybind11::arg("nfmaxIn"));
		cl.def("init", (void (Pythia8::AlphaStrong::*)(double, int, int, bool)) &Pythia8::AlphaStrong::init, "C++: Pythia8::AlphaStrong::init(double, int, int, bool) --> void", pybind11::arg("valueIn"), pybind11::arg("orderIn"), pybind11::arg("nfmaxIn"), pybind11::arg("useCMWIn"));
		cl.def("setThresholds", (void (Pythia8::AlphaStrong::*)(double, double, double)) &Pythia8::AlphaStrong::setThresholds, "C++: Pythia8::AlphaStrong::setThresholds(double, double, double) --> void", pybind11::arg("mcIn"), pybind11::arg("mbIn"), pybind11::arg("mtIn"));
		cl.def("alphaS", (double (Pythia8::AlphaStrong::*)(double)) &Pythia8::AlphaStrong::alphaS, "C++: Pythia8::AlphaStrong::alphaS(double) --> double", pybind11::arg("scale2"));
		cl.def("alphaS1Ord", (double (Pythia8::AlphaStrong::*)(double)) &Pythia8::AlphaStrong::alphaS1Ord, "C++: Pythia8::AlphaStrong::alphaS1Ord(double) --> double", pybind11::arg("scale2"));
		cl.def("alphaS2OrdCorr", (double (Pythia8::AlphaStrong::*)(double)) &Pythia8::AlphaStrong::alphaS2OrdCorr, "C++: Pythia8::AlphaStrong::alphaS2OrdCorr(double) --> double", pybind11::arg("scale2"));
		cl.def("Lambda3", (double (Pythia8::AlphaStrong::*)() const) &Pythia8::AlphaStrong::Lambda3, "C++: Pythia8::AlphaStrong::Lambda3() const --> double");
		cl.def("Lambda4", (double (Pythia8::AlphaStrong::*)() const) &Pythia8::AlphaStrong::Lambda4, "C++: Pythia8::AlphaStrong::Lambda4() const --> double");
		cl.def("Lambda5", (double (Pythia8::AlphaStrong::*)() const) &Pythia8::AlphaStrong::Lambda5, "C++: Pythia8::AlphaStrong::Lambda5() const --> double");
		cl.def("Lambda6", (double (Pythia8::AlphaStrong::*)() const) &Pythia8::AlphaStrong::Lambda6, "C++: Pythia8::AlphaStrong::Lambda6() const --> double");
		cl.def("muThres", (double (Pythia8::AlphaStrong::*)(int)) &Pythia8::AlphaStrong::muThres, "C++: Pythia8::AlphaStrong::muThres(int) --> double", pybind11::arg("idQ"));
		cl.def("muThres2", (double (Pythia8::AlphaStrong::*)(int)) &Pythia8::AlphaStrong::muThres2, "C++: Pythia8::AlphaStrong::muThres2(int) --> double", pybind11::arg("idQ"));
		cl.def("facCMW", (double (Pythia8::AlphaStrong::*)(int)) &Pythia8::AlphaStrong::facCMW, "C++: Pythia8::AlphaStrong::facCMW(int) --> double", pybind11::arg("nFin"));
		cl.def("assign", (class Pythia8::AlphaStrong & (Pythia8::AlphaStrong::*)(const class Pythia8::AlphaStrong &)) &Pythia8::AlphaStrong::operator=, "C++: Pythia8::AlphaStrong::operator=(const class Pythia8::AlphaStrong &) --> class Pythia8::AlphaStrong &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::AlphaEM file:Pythia8/StandardModel.h line:106
		pybind11::class_<Pythia8::AlphaEM, std::shared_ptr<Pythia8::AlphaEM>> cl(M("Pythia8"), "AlphaEM", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::AlphaEM(); } ) );
		cl.def( pybind11::init( [](Pythia8::AlphaEM const &o){ return new Pythia8::AlphaEM(o); } ) );
		cl.def("init", (void (Pythia8::AlphaEM::*)(int, class Pythia8::Settings *)) &Pythia8::AlphaEM::init, "C++: Pythia8::AlphaEM::init(int, class Pythia8::Settings *) --> void", pybind11::arg("orderIn"), pybind11::arg("settingsPtr"));
		cl.def("alphaEM", (double (Pythia8::AlphaEM::*)(double)) &Pythia8::AlphaEM::alphaEM, "C++: Pythia8::AlphaEM::alphaEM(double) --> double", pybind11::arg("scale2"));
		cl.def("assign", (class Pythia8::AlphaEM & (Pythia8::AlphaEM::*)(const class Pythia8::AlphaEM &)) &Pythia8::AlphaEM::operator=, "C++: Pythia8::AlphaEM::operator=(const class Pythia8::AlphaEM &) --> class Pythia8::AlphaEM &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::CoupSM file:Pythia8/StandardModel.h line:135
		pybind11::class_<Pythia8::CoupSM, std::shared_ptr<Pythia8::CoupSM>> cl(M("Pythia8"), "CoupSM", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::CoupSM(); } ) );
		cl.def( pybind11::init( [](Pythia8::CoupSM const &o){ return new Pythia8::CoupSM(o); } ) );
		cl.def_readwrite("s2tW", &Pythia8::CoupSM::s2tW);
		cl.def_readwrite("c2tW", &Pythia8::CoupSM::c2tW);
		cl.def_readwrite("s2tWbar", &Pythia8::CoupSM::s2tWbar);
		cl.def_readwrite("GFermi", &Pythia8::CoupSM::GFermi);
		cl.def_readwrite("alphaSlocal", &Pythia8::CoupSM::alphaSlocal);
		cl.def_readwrite("alphaEMlocal", &Pythia8::CoupSM::alphaEMlocal);
		cl.def("init", (void (Pythia8::CoupSM::*)(class Pythia8::Settings &, class Pythia8::Rndm *)) &Pythia8::CoupSM::init, "C++: Pythia8::CoupSM::init(class Pythia8::Settings &, class Pythia8::Rndm *) --> void", pybind11::arg("settings"), pybind11::arg("rndmPtrIn"));
		cl.def("alphaS", (double (Pythia8::CoupSM::*)(double)) &Pythia8::CoupSM::alphaS, "C++: Pythia8::CoupSM::alphaS(double) --> double", pybind11::arg("scale2"));
		cl.def("alphaS1Ord", (double (Pythia8::CoupSM::*)(double)) &Pythia8::CoupSM::alphaS1Ord, "C++: Pythia8::CoupSM::alphaS1Ord(double) --> double", pybind11::arg("scale2"));
		cl.def("alphaS2OrdCorr", (double (Pythia8::CoupSM::*)(double)) &Pythia8::CoupSM::alphaS2OrdCorr, "C++: Pythia8::CoupSM::alphaS2OrdCorr(double) --> double", pybind11::arg("scale2"));
		cl.def("Lambda3", (double (Pythia8::CoupSM::*)() const) &Pythia8::CoupSM::Lambda3, "C++: Pythia8::CoupSM::Lambda3() const --> double");
		cl.def("Lambda4", (double (Pythia8::CoupSM::*)() const) &Pythia8::CoupSM::Lambda4, "C++: Pythia8::CoupSM::Lambda4() const --> double");
		cl.def("Lambda5", (double (Pythia8::CoupSM::*)() const) &Pythia8::CoupSM::Lambda5, "C++: Pythia8::CoupSM::Lambda5() const --> double");
		cl.def("alphaEM", (double (Pythia8::CoupSM::*)(double)) &Pythia8::CoupSM::alphaEM, "C++: Pythia8::CoupSM::alphaEM(double) --> double", pybind11::arg("scale2"));
		cl.def("sin2thetaW", (double (Pythia8::CoupSM::*)()) &Pythia8::CoupSM::sin2thetaW, "C++: Pythia8::CoupSM::sin2thetaW() --> double");
		cl.def("cos2thetaW", (double (Pythia8::CoupSM::*)()) &Pythia8::CoupSM::cos2thetaW, "C++: Pythia8::CoupSM::cos2thetaW() --> double");
		cl.def("sin2thetaWbar", (double (Pythia8::CoupSM::*)()) &Pythia8::CoupSM::sin2thetaWbar, "C++: Pythia8::CoupSM::sin2thetaWbar() --> double");
		cl.def("GF", (double (Pythia8::CoupSM::*)()) &Pythia8::CoupSM::GF, "C++: Pythia8::CoupSM::GF() --> double");
		cl.def("ef", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::ef, "C++: Pythia8::CoupSM::ef(int) --> double", pybind11::arg("idAbs"));
		cl.def("vf", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::vf, "C++: Pythia8::CoupSM::vf(int) --> double", pybind11::arg("idAbs"));
		cl.def("af", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::af, "C++: Pythia8::CoupSM::af(int) --> double", pybind11::arg("idAbs"));
		cl.def("t3f", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::t3f, "C++: Pythia8::CoupSM::t3f(int) --> double", pybind11::arg("idAbs"));
		cl.def("lf", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::lf, "C++: Pythia8::CoupSM::lf(int) --> double", pybind11::arg("idAbs"));
		cl.def("rf", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::rf, "C++: Pythia8::CoupSM::rf(int) --> double", pybind11::arg("idAbs"));
		cl.def("ef2", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::ef2, "C++: Pythia8::CoupSM::ef2(int) --> double", pybind11::arg("idAbs"));
		cl.def("vf2", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::vf2, "C++: Pythia8::CoupSM::vf2(int) --> double", pybind11::arg("idAbs"));
		cl.def("af2", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::af2, "C++: Pythia8::CoupSM::af2(int) --> double", pybind11::arg("idAbs"));
		cl.def("efvf", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::efvf, "C++: Pythia8::CoupSM::efvf(int) --> double", pybind11::arg("idAbs"));
		cl.def("vf2af2", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::vf2af2, "C++: Pythia8::CoupSM::vf2af2(int) --> double", pybind11::arg("idAbs"));
		cl.def("VCKMgen", (double (Pythia8::CoupSM::*)(int, int)) &Pythia8::CoupSM::VCKMgen, "C++: Pythia8::CoupSM::VCKMgen(int, int) --> double", pybind11::arg("genU"), pybind11::arg("genD"));
		cl.def("V2CKMgen", (double (Pythia8::CoupSM::*)(int, int)) &Pythia8::CoupSM::V2CKMgen, "C++: Pythia8::CoupSM::V2CKMgen(int, int) --> double", pybind11::arg("genU"), pybind11::arg("genD"));
		cl.def("VCKMid", (double (Pythia8::CoupSM::*)(int, int)) &Pythia8::CoupSM::VCKMid, "C++: Pythia8::CoupSM::VCKMid(int, int) --> double", pybind11::arg("id1"), pybind11::arg("id2"));
		cl.def("V2CKMid", (double (Pythia8::CoupSM::*)(int, int)) &Pythia8::CoupSM::V2CKMid, "C++: Pythia8::CoupSM::V2CKMid(int, int) --> double", pybind11::arg("id1"), pybind11::arg("id2"));
		cl.def("V2CKMsum", (double (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::V2CKMsum, "C++: Pythia8::CoupSM::V2CKMsum(int) --> double", pybind11::arg("id"));
		cl.def("V2CKMpick", (int (Pythia8::CoupSM::*)(int)) &Pythia8::CoupSM::V2CKMpick, "C++: Pythia8::CoupSM::V2CKMpick(int) --> int", pybind11::arg("id"));
		cl.def("assign", (class Pythia8::CoupSM & (Pythia8::CoupSM::*)(const class Pythia8::CoupSM &)) &Pythia8::CoupSM::operator=, "C++: Pythia8::CoupSM::operator=(const class Pythia8::CoupSM &) --> class Pythia8::CoupSM &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::AlphaSUN file:Pythia8/StandardModel.h line:226
		pybind11::class_<Pythia8::AlphaSUN, std::shared_ptr<Pythia8::AlphaSUN>, PyCallBack_Pythia8_AlphaSUN> cl(M("Pythia8"), "AlphaSUN", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::AlphaSUN(); }, [](){ return new PyCallBack_Pythia8_AlphaSUN(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_AlphaSUN const &o){ return new PyCallBack_Pythia8_AlphaSUN(o); } ) );
		cl.def( pybind11::init( [](Pythia8::AlphaSUN const &o){ return new Pythia8::AlphaSUN(o); } ) );
		cl.def("initAlpha", [](Pythia8::AlphaSUN &o, int const & a0, int const & a1) -> void { return o.initAlpha(a0, a1); }, "", pybind11::arg("nCin"), pybind11::arg("nFin"));
		cl.def("initAlpha", [](Pythia8::AlphaSUN &o, int const & a0, int const & a1, int const & a2) -> void { return o.initAlpha(a0, a1, a2); }, "", pybind11::arg("nCin"), pybind11::arg("nFin"), pybind11::arg("orderIn"));
		cl.def("initAlpha", [](Pythia8::AlphaSUN &o, int const & a0, int const & a1, int const & a2, double const & a3) -> void { return o.initAlpha(a0, a1, a2, a3); }, "", pybind11::arg("nCin"), pybind11::arg("nFin"), pybind11::arg("orderIn"), pybind11::arg("alphaIn"));
		cl.def("initAlpha", (void (Pythia8::AlphaSUN::*)(int, int, int, double, double)) &Pythia8::AlphaSUN::initAlpha, "C++: Pythia8::AlphaSUN::initAlpha(int, int, int, double, double) --> void", pybind11::arg("nCin"), pybind11::arg("nFin"), pybind11::arg("orderIn"), pybind11::arg("alphaIn"), pybind11::arg("scaleIn"));
		cl.def("initLambda", [](Pythia8::AlphaSUN &o, int const & a0, int const & a1) -> void { return o.initLambda(a0, a1); }, "", pybind11::arg("nCin"), pybind11::arg("nFin"));
		cl.def("initLambda", [](Pythia8::AlphaSUN &o, int const & a0, int const & a1, int const & a2) -> void { return o.initLambda(a0, a1, a2); }, "", pybind11::arg("nCin"), pybind11::arg("nFin"), pybind11::arg("orderIn"));
		cl.def("initLambda", (void (Pythia8::AlphaSUN::*)(int, int, int, double)) &Pythia8::AlphaSUN::initLambda, "C++: Pythia8::AlphaSUN::initLambda(int, int, int, double) --> void", pybind11::arg("nCin"), pybind11::arg("nFin"), pybind11::arg("orderIn"), pybind11::arg("LambdaIn"));
		cl.def("alpha", (double (Pythia8::AlphaSUN::*)(double)) &Pythia8::AlphaSUN::alpha, "C++: Pythia8::AlphaSUN::alpha(double) --> double", pybind11::arg("scale2in"));
		cl.def("alpha1Ord", (double (Pythia8::AlphaSUN::*)(double)) &Pythia8::AlphaSUN::alpha1Ord, "C++: Pythia8::AlphaSUN::alpha1Ord(double) --> double", pybind11::arg("scale2in"));
		cl.def("alpha2OrdCorr", (double (Pythia8::AlphaSUN::*)(double)) &Pythia8::AlphaSUN::alpha2OrdCorr, "C++: Pythia8::AlphaSUN::alpha2OrdCorr(double) --> double", pybind11::arg("scale2in"));
		cl.def("Lambda", (double (Pythia8::AlphaSUN::*)() const) &Pythia8::AlphaSUN::Lambda, "C++: Pythia8::AlphaSUN::Lambda() const --> double");
		cl.def("assign", (class Pythia8::AlphaSUN & (Pythia8::AlphaSUN::*)(const class Pythia8::AlphaSUN &)) &Pythia8::AlphaSUN::operator=, "C++: Pythia8::AlphaSUN::operator=(const class Pythia8::AlphaSUN &) --> class Pythia8::AlphaSUN &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::DecayChannel file:Pythia8/ParticleData.h line:35
		pybind11::class_<Pythia8::DecayChannel, std::shared_ptr<Pythia8::DecayChannel>> cl(M("Pythia8"), "DecayChannel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::DecayChannel(); } ), "doc" );
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::DecayChannel(a0); } ), "doc" , pybind11::arg("onModeIn"));
		cl.def( pybind11::init( [](int const & a0, double const & a1){ return new Pythia8::DecayChannel(a0, a1); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, int const & a2){ return new Pythia8::DecayChannel(a0, a1, a2); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, int const & a2, int const & a3){ return new Pythia8::DecayChannel(a0, a1, a2, a3); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"), pybind11::arg("prod0"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, int const & a2, int const & a3, int const & a4){ return new Pythia8::DecayChannel(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"), pybind11::arg("prod0"), pybind11::arg("prod1"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, int const & a2, int const & a3, int const & a4, int const & a5){ return new Pythia8::DecayChannel(a0, a1, a2, a3, a4, a5); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"), pybind11::arg("prod0"), pybind11::arg("prod1"), pybind11::arg("prod2"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6){ return new Pythia8::DecayChannel(a0, a1, a2, a3, a4, a5, a6); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"), pybind11::arg("prod0"), pybind11::arg("prod1"), pybind11::arg("prod2"), pybind11::arg("prod3"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7){ return new Pythia8::DecayChannel(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"), pybind11::arg("prod0"), pybind11::arg("prod1"), pybind11::arg("prod2"), pybind11::arg("prod3"), pybind11::arg("prod4"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8){ return new Pythia8::DecayChannel(a0, a1, a2, a3, a4, a5, a6, a7, a8); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"), pybind11::arg("prod0"), pybind11::arg("prod1"), pybind11::arg("prod2"), pybind11::arg("prod3"), pybind11::arg("prod4"), pybind11::arg("prod5"));
		cl.def( pybind11::init( [](int const & a0, double const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9){ return new Pythia8::DecayChannel(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); } ), "doc" , pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"), pybind11::arg("prod0"), pybind11::arg("prod1"), pybind11::arg("prod2"), pybind11::arg("prod3"), pybind11::arg("prod4"), pybind11::arg("prod5"), pybind11::arg("prod6"));
		cl.def( pybind11::init<int, double, int, int, int, int, int, int, int, int, int>(), pybind11::arg("onModeIn"), pybind11::arg("bRatioIn"), pybind11::arg("meModeIn"), pybind11::arg("prod0"), pybind11::arg("prod1"), pybind11::arg("prod2"), pybind11::arg("prod3"), pybind11::arg("prod4"), pybind11::arg("prod5"), pybind11::arg("prod6"), pybind11::arg("prod7") );

		cl.def( pybind11::init( [](Pythia8::DecayChannel const &o){ return new Pythia8::DecayChannel(o); } ) );
		cl.def("assign", (class Pythia8::DecayChannel & (Pythia8::DecayChannel::*)(const class Pythia8::DecayChannel &)) &Pythia8::DecayChannel::operator=, "C++: Pythia8::DecayChannel::operator=(const class Pythia8::DecayChannel &) --> class Pythia8::DecayChannel &", pybind11::return_value_policy::reference, pybind11::arg("oldDC"));
		cl.def("onMode", (void (Pythia8::DecayChannel::*)(int)) &Pythia8::DecayChannel::onMode, "C++: Pythia8::DecayChannel::onMode(int) --> void", pybind11::arg("onModeIn"));
		cl.def("bRatio", [](Pythia8::DecayChannel &o, double const & a0) -> void { return o.bRatio(a0); }, "", pybind11::arg("bRatioIn"));
		cl.def("bRatio", (void (Pythia8::DecayChannel::*)(double, bool)) &Pythia8::DecayChannel::bRatio, "C++: Pythia8::DecayChannel::bRatio(double, bool) --> void", pybind11::arg("bRatioIn"), pybind11::arg("countAsChanged"));
		cl.def("rescaleBR", (void (Pythia8::DecayChannel::*)(double)) &Pythia8::DecayChannel::rescaleBR, "C++: Pythia8::DecayChannel::rescaleBR(double) --> void", pybind11::arg("fac"));
		cl.def("meMode", (void (Pythia8::DecayChannel::*)(int)) &Pythia8::DecayChannel::meMode, "C++: Pythia8::DecayChannel::meMode(int) --> void", pybind11::arg("meModeIn"));
		cl.def("multiplicity", (void (Pythia8::DecayChannel::*)(int)) &Pythia8::DecayChannel::multiplicity, "C++: Pythia8::DecayChannel::multiplicity(int) --> void", pybind11::arg("multIn"));
		cl.def("product", (void (Pythia8::DecayChannel::*)(int, int)) &Pythia8::DecayChannel::product, "C++: Pythia8::DecayChannel::product(int, int) --> void", pybind11::arg("i"), pybind11::arg("prodIn"));
		cl.def("setHasChanged", (void (Pythia8::DecayChannel::*)(bool)) &Pythia8::DecayChannel::setHasChanged, "C++: Pythia8::DecayChannel::setHasChanged(bool) --> void", pybind11::arg("hasChangedIn"));
		cl.def("onMode", (int (Pythia8::DecayChannel::*)() const) &Pythia8::DecayChannel::onMode, "C++: Pythia8::DecayChannel::onMode() const --> int");
		cl.def("bRatio", (double (Pythia8::DecayChannel::*)() const) &Pythia8::DecayChannel::bRatio, "C++: Pythia8::DecayChannel::bRatio() const --> double");
		cl.def("meMode", (int (Pythia8::DecayChannel::*)() const) &Pythia8::DecayChannel::meMode, "C++: Pythia8::DecayChannel::meMode() const --> int");
		cl.def("multiplicity", (int (Pythia8::DecayChannel::*)() const) &Pythia8::DecayChannel::multiplicity, "C++: Pythia8::DecayChannel::multiplicity() const --> int");
		cl.def("product", (int (Pythia8::DecayChannel::*)(int) const) &Pythia8::DecayChannel::product, "C++: Pythia8::DecayChannel::product(int) const --> int", pybind11::arg("i"));
		cl.def("hasChanged", (bool (Pythia8::DecayChannel::*)() const) &Pythia8::DecayChannel::hasChanged, "C++: Pythia8::DecayChannel::hasChanged() const --> bool");
		cl.def("contains", (bool (Pythia8::DecayChannel::*)(int) const) &Pythia8::DecayChannel::contains, "C++: Pythia8::DecayChannel::contains(int) const --> bool", pybind11::arg("id1"));
		cl.def("contains", (bool (Pythia8::DecayChannel::*)(int, int) const) &Pythia8::DecayChannel::contains, "C++: Pythia8::DecayChannel::contains(int, int) const --> bool", pybind11::arg("id1"), pybind11::arg("id2"));
		cl.def("contains", (bool (Pythia8::DecayChannel::*)(int, int, int) const) &Pythia8::DecayChannel::contains, "C++: Pythia8::DecayChannel::contains(int, int, int) const --> bool", pybind11::arg("id1"), pybind11::arg("id2"), pybind11::arg("id3"));
		cl.def("currentBR", (void (Pythia8::DecayChannel::*)(double)) &Pythia8::DecayChannel::currentBR, "C++: Pythia8::DecayChannel::currentBR(double) --> void", pybind11::arg("currentBRIn"));
		cl.def("currentBR", (double (Pythia8::DecayChannel::*)() const) &Pythia8::DecayChannel::currentBR, "C++: Pythia8::DecayChannel::currentBR() const --> double");
		cl.def("onShellWidth", (void (Pythia8::DecayChannel::*)(double)) &Pythia8::DecayChannel::onShellWidth, "C++: Pythia8::DecayChannel::onShellWidth(double) --> void", pybind11::arg("onShellWidthIn"));
		cl.def("onShellWidth", (double (Pythia8::DecayChannel::*)() const) &Pythia8::DecayChannel::onShellWidth, "C++: Pythia8::DecayChannel::onShellWidth() const --> double");
		cl.def("onShellWidthFactor", (void (Pythia8::DecayChannel::*)(double)) &Pythia8::DecayChannel::onShellWidthFactor, "C++: Pythia8::DecayChannel::onShellWidthFactor(double) --> void", pybind11::arg("factor"));
		cl.def("openSec", (void (Pythia8::DecayChannel::*)(int, double)) &Pythia8::DecayChannel::openSec, "C++: Pythia8::DecayChannel::openSec(int, double) --> void", pybind11::arg("idSgn"), pybind11::arg("openSecIn"));
		cl.def("openSec", (double (Pythia8::DecayChannel::*)(int) const) &Pythia8::DecayChannel::openSec, "C++: Pythia8::DecayChannel::openSec(int) const --> double", pybind11::arg("idSgn"));
	}
}
