#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/Weights.h>
#include <cwchar>
#include <fstream>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <sstream> // __str__
#include <streambuf>
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

// Pythia8::LHAup file:Pythia8/LesHouches.h line:78
struct PyCallBack_Pythia8_LHAup : public Pythia8::LHAup {
	using Pythia8::LHAup::LHAup;

	void newEventFile(const char * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAup *>(this), "newEventFile");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LHAup::newEventFile(a0);
	}
	bool fileFound() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAup *>(this), "fileFound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::fileFound();
	}
	bool useExternal() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAup *>(this), "useExternal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::useExternal();
	}
	bool setInit() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAup *>(this), "setInit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"LHAup::setInit\"");
	}
	bool setEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAup *>(this), "setEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"LHAup::setEvent\"");
	}
	bool skipEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAup *>(this), "skipEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::skipEvent(a0);
	}
	bool openLHEF(class std::basic_string<char> a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAup *>(this), "openLHEF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::openLHEF(a0);
	}
	bool closeLHEF(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAup *>(this), "closeLHEF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::closeLHEF(a0);
	}
};

void bind_Pythia8_LesHouches(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::LHAProcess file:Pythia8/LesHouches.h line:28
		pybind11::class_<Pythia8::LHAProcess, std::shared_ptr<Pythia8::LHAProcess>> cl(M("Pythia8"), "LHAProcess", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAProcess(); } ) );
		cl.def( pybind11::init<int, double, double, double>(), pybind11::arg("idProcIn"), pybind11::arg("xSecIn"), pybind11::arg("xErrIn"), pybind11::arg("xMaxIn") );

		cl.def( pybind11::init( [](Pythia8::LHAProcess const &o){ return new Pythia8::LHAProcess(o); } ) );
		cl.def_readwrite("idProc", &Pythia8::LHAProcess::idProc);
		cl.def_readwrite("xSecProc", &Pythia8::LHAProcess::xSecProc);
		cl.def_readwrite("xErrProc", &Pythia8::LHAProcess::xErrProc);
		cl.def_readwrite("xMaxProc", &Pythia8::LHAProcess::xMaxProc);
		cl.def("assign", (class Pythia8::LHAProcess & (Pythia8::LHAProcess::*)(const class Pythia8::LHAProcess &)) &Pythia8::LHAProcess::operator=, "C++: Pythia8::LHAProcess::operator=(const class Pythia8::LHAProcess &) --> class Pythia8::LHAProcess &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LHAParticle file:Pythia8/LesHouches.h line:48
		pybind11::class_<Pythia8::LHAParticle, std::shared_ptr<Pythia8::LHAParticle>> cl(M("Pythia8"), "LHAParticle", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAParticle(); } ) );
		cl.def( pybind11::init<int, int, int, int, int, int, double, double, double, double, double, double, double, double>(), pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("pzIn"), pybind11::arg("eIn"), pybind11::arg("mIn"), pybind11::arg("tauIn"), pybind11::arg("spinIn"), pybind11::arg("scaleIn") );

		cl.def( pybind11::init( [](Pythia8::LHAParticle const &o){ return new Pythia8::LHAParticle(o); } ) );
		cl.def_readwrite("idPart", &Pythia8::LHAParticle::idPart);
		cl.def_readwrite("statusPart", &Pythia8::LHAParticle::statusPart);
		cl.def_readwrite("mother1Part", &Pythia8::LHAParticle::mother1Part);
		cl.def_readwrite("mother2Part", &Pythia8::LHAParticle::mother2Part);
		cl.def_readwrite("col1Part", &Pythia8::LHAParticle::col1Part);
		cl.def_readwrite("col2Part", &Pythia8::LHAParticle::col2Part);
		cl.def_readwrite("pxPart", &Pythia8::LHAParticle::pxPart);
		cl.def_readwrite("pyPart", &Pythia8::LHAParticle::pyPart);
		cl.def_readwrite("pzPart", &Pythia8::LHAParticle::pzPart);
		cl.def_readwrite("ePart", &Pythia8::LHAParticle::ePart);
		cl.def_readwrite("mPart", &Pythia8::LHAParticle::mPart);
		cl.def_readwrite("tauPart", &Pythia8::LHAParticle::tauPart);
		cl.def_readwrite("spinPart", &Pythia8::LHAParticle::spinPart);
		cl.def_readwrite("scalePart", &Pythia8::LHAParticle::scalePart);
		cl.def("assign", (class Pythia8::LHAParticle & (Pythia8::LHAParticle::*)(const class Pythia8::LHAParticle &)) &Pythia8::LHAParticle::operator=, "C++: Pythia8::LHAParticle::operator=(const class Pythia8::LHAParticle &) --> class Pythia8::LHAParticle &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LHAup file:Pythia8/LesHouches.h line:78
		pybind11::class_<Pythia8::LHAup, std::shared_ptr<Pythia8::LHAup>, PyCallBack_Pythia8_LHAup> cl(M("Pythia8"), "LHAup", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_LHAup(); } ), "doc");
		cl.def( pybind11::init<int>(), pybind11::arg("strategyIn") );

		cl.def_readwrite("nupSave", &Pythia8::LHAup::nupSave);
		cl.def_readwrite("idprupSave", &Pythia8::LHAup::idprupSave);
		cl.def_readwrite("xwgtupSave", &Pythia8::LHAup::xwgtupSave);
		cl.def_readwrite("scalupSave", &Pythia8::LHAup::scalupSave);
		cl.def_readwrite("aqedupSave", &Pythia8::LHAup::aqedupSave);
		cl.def_readwrite("aqcdupSave", &Pythia8::LHAup::aqcdupSave);
		cl.def_readwrite("xSecSumSave", &Pythia8::LHAup::xSecSumSave);
		cl.def_readwrite("xErrSumSave", &Pythia8::LHAup::xErrSumSave);
		cl.def_readwrite("particlesSave", &Pythia8::LHAup::particlesSave);
		cl.def_readwrite("getPDFSave", &Pythia8::LHAup::getPDFSave);
		cl.def_readwrite("getScale", &Pythia8::LHAup::getScale);
		cl.def_readwrite("getScaleShowers", &Pythia8::LHAup::getScaleShowers);
		cl.def_readwrite("id1InSave", &Pythia8::LHAup::id1InSave);
		cl.def_readwrite("id2InSave", &Pythia8::LHAup::id2InSave);
		cl.def_readwrite("id1pdfInSave", &Pythia8::LHAup::id1pdfInSave);
		cl.def_readwrite("id2pdfInSave", &Pythia8::LHAup::id2pdfInSave);
		cl.def_readwrite("x1InSave", &Pythia8::LHAup::x1InSave);
		cl.def_readwrite("x2InSave", &Pythia8::LHAup::x2InSave);
		cl.def_readwrite("x1pdfInSave", &Pythia8::LHAup::x1pdfInSave);
		cl.def_readwrite("x2pdfInSave", &Pythia8::LHAup::x2pdfInSave);
		cl.def_readwrite("scalePDFInSave", &Pythia8::LHAup::scalePDFInSave);
		cl.def_readwrite("pdf1InSave", &Pythia8::LHAup::pdf1InSave);
		cl.def_readwrite("pdf2InSave", &Pythia8::LHAup::pdf2InSave);
		cl.def_readwrite("fileName", &Pythia8::LHAup::fileName);
		cl.def("setPtr", (void (Pythia8::LHAup::*)(class Pythia8::Info *)) &Pythia8::LHAup::setPtr, "C++: Pythia8::LHAup::setPtr(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("newEventFile", (void (Pythia8::LHAup::*)(const char *)) &Pythia8::LHAup::newEventFile, "C++: Pythia8::LHAup::newEventFile(const char *) --> void", pybind11::arg(""));
		cl.def("fileFound", (bool (Pythia8::LHAup::*)()) &Pythia8::LHAup::fileFound, "C++: Pythia8::LHAup::fileFound() --> bool");
		cl.def("useExternal", (bool (Pythia8::LHAup::*)()) &Pythia8::LHAup::useExternal, "C++: Pythia8::LHAup::useExternal() --> bool");
		cl.def("setInit", (bool (Pythia8::LHAup::*)()) &Pythia8::LHAup::setInit, "C++: Pythia8::LHAup::setInit() --> bool");
		cl.def("idBeamA", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::idBeamA, "C++: Pythia8::LHAup::idBeamA() const --> int");
		cl.def("idBeamB", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::idBeamB, "C++: Pythia8::LHAup::idBeamB() const --> int");
		cl.def("eBeamA", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::eBeamA, "C++: Pythia8::LHAup::eBeamA() const --> double");
		cl.def("eBeamB", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::eBeamB, "C++: Pythia8::LHAup::eBeamB() const --> double");
		cl.def("pdfGroupBeamA", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::pdfGroupBeamA, "C++: Pythia8::LHAup::pdfGroupBeamA() const --> int");
		cl.def("pdfGroupBeamB", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::pdfGroupBeamB, "C++: Pythia8::LHAup::pdfGroupBeamB() const --> int");
		cl.def("pdfSetBeamA", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::pdfSetBeamA, "C++: Pythia8::LHAup::pdfSetBeamA() const --> int");
		cl.def("pdfSetBeamB", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::pdfSetBeamB, "C++: Pythia8::LHAup::pdfSetBeamB() const --> int");
		cl.def("strategy", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::strategy, "C++: Pythia8::LHAup::strategy() const --> int");
		cl.def("sizeProc", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::sizeProc, "C++: Pythia8::LHAup::sizeProc() const --> int");
		cl.def("idProcess", (int (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::idProcess, "C++: Pythia8::LHAup::idProcess(int) const --> int", pybind11::arg("proc"));
		cl.def("xSec", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::xSec, "C++: Pythia8::LHAup::xSec(int) const --> double", pybind11::arg("proc"));
		cl.def("xErr", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::xErr, "C++: Pythia8::LHAup::xErr(int) const --> double", pybind11::arg("proc"));
		cl.def("xMax", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::xMax, "C++: Pythia8::LHAup::xMax(int) const --> double", pybind11::arg("proc"));
		cl.def("xSecSum", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::xSecSum, "C++: Pythia8::LHAup::xSecSum() const --> double");
		cl.def("xErrSum", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::xErrSum, "C++: Pythia8::LHAup::xErrSum() const --> double");
		cl.def("listInit", (void (Pythia8::LHAup::*)()) &Pythia8::LHAup::listInit, "C++: Pythia8::LHAup::listInit() --> void");
		cl.def("setEvent", [](Pythia8::LHAup &o) -> bool { return o.setEvent(); }, "");
		cl.def("setEvent", (bool (Pythia8::LHAup::*)(int)) &Pythia8::LHAup::setEvent, "C++: Pythia8::LHAup::setEvent(int) --> bool", pybind11::arg("idProcIn"));
		cl.def("idProcess", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::idProcess, "C++: Pythia8::LHAup::idProcess() const --> int");
		cl.def("weight", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::weight, "C++: Pythia8::LHAup::weight() const --> double");
		cl.def("scale", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::scale, "C++: Pythia8::LHAup::scale() const --> double");
		cl.def("alphaQED", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::alphaQED, "C++: Pythia8::LHAup::alphaQED() const --> double");
		cl.def("alphaQCD", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::alphaQCD, "C++: Pythia8::LHAup::alphaQCD() const --> double");
		cl.def("sizePart", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::sizePart, "C++: Pythia8::LHAup::sizePart() const --> int");
		cl.def("id", (int (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::id, "C++: Pythia8::LHAup::id(int) const --> int", pybind11::arg("part"));
		cl.def("status", (int (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::status, "C++: Pythia8::LHAup::status(int) const --> int", pybind11::arg("part"));
		cl.def("mother1", (int (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::mother1, "C++: Pythia8::LHAup::mother1(int) const --> int", pybind11::arg("part"));
		cl.def("mother2", (int (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::mother2, "C++: Pythia8::LHAup::mother2(int) const --> int", pybind11::arg("part"));
		cl.def("col1", (int (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::col1, "C++: Pythia8::LHAup::col1(int) const --> int", pybind11::arg("part"));
		cl.def("col2", (int (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::col2, "C++: Pythia8::LHAup::col2(int) const --> int", pybind11::arg("part"));
		cl.def("px", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::px, "C++: Pythia8::LHAup::px(int) const --> double", pybind11::arg("part"));
		cl.def("py", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::py, "C++: Pythia8::LHAup::py(int) const --> double", pybind11::arg("part"));
		cl.def("pz", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::pz, "C++: Pythia8::LHAup::pz(int) const --> double", pybind11::arg("part"));
		cl.def("e", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::e, "C++: Pythia8::LHAup::e(int) const --> double", pybind11::arg("part"));
		cl.def("m", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::m, "C++: Pythia8::LHAup::m(int) const --> double", pybind11::arg("part"));
		cl.def("tau", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::tau, "C++: Pythia8::LHAup::tau(int) const --> double", pybind11::arg("part"));
		cl.def("spin", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::spin, "C++: Pythia8::LHAup::spin(int) const --> double", pybind11::arg("part"));
		cl.def("scale", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::scale, "C++: Pythia8::LHAup::scale(int) const --> double", pybind11::arg("part"));
		cl.def("id1", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::id1, "C++: Pythia8::LHAup::id1() const --> int");
		cl.def("id2", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::id2, "C++: Pythia8::LHAup::id2() const --> int");
		cl.def("x1", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::x1, "C++: Pythia8::LHAup::x1() const --> double");
		cl.def("x2", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::x2, "C++: Pythia8::LHAup::x2() const --> double");
		cl.def("pdfIsSet", (bool (Pythia8::LHAup::*)() const) &Pythia8::LHAup::pdfIsSet, "C++: Pythia8::LHAup::pdfIsSet() const --> bool");
		cl.def("id1pdf", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::id1pdf, "C++: Pythia8::LHAup::id1pdf() const --> int");
		cl.def("id2pdf", (int (Pythia8::LHAup::*)() const) &Pythia8::LHAup::id2pdf, "C++: Pythia8::LHAup::id2pdf() const --> int");
		cl.def("x1pdf", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::x1pdf, "C++: Pythia8::LHAup::x1pdf() const --> double");
		cl.def("x2pdf", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::x2pdf, "C++: Pythia8::LHAup::x2pdf() const --> double");
		cl.def("scalePDF", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::scalePDF, "C++: Pythia8::LHAup::scalePDF() const --> double");
		cl.def("pdf1", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::pdf1, "C++: Pythia8::LHAup::pdf1() const --> double");
		cl.def("pdf2", (double (Pythia8::LHAup::*)() const) &Pythia8::LHAup::pdf2, "C++: Pythia8::LHAup::pdf2() const --> double");
		cl.def("scaleShowersIsSet", (bool (Pythia8::LHAup::*)() const) &Pythia8::LHAup::scaleShowersIsSet, "C++: Pythia8::LHAup::scaleShowersIsSet() const --> bool");
		cl.def("scaleShowers", (double (Pythia8::LHAup::*)(int) const) &Pythia8::LHAup::scaleShowers, "C++: Pythia8::LHAup::scaleShowers(int) const --> double", pybind11::arg("i"));
		cl.def("listEvent", (void (Pythia8::LHAup::*)()) &Pythia8::LHAup::listEvent, "C++: Pythia8::LHAup::listEvent() --> void");
		cl.def("skipEvent", (bool (Pythia8::LHAup::*)(int)) &Pythia8::LHAup::skipEvent, "C++: Pythia8::LHAup::skipEvent(int) --> bool", pybind11::arg("nSkip"));
		cl.def("openLHEF", (bool (Pythia8::LHAup::*)(std::string)) &Pythia8::LHAup::openLHEF, "C++: Pythia8::LHAup::openLHEF(std::string) --> bool", pybind11::arg("fileNameIn"));
		cl.def("closeLHEF", [](Pythia8::LHAup &o) -> bool { return o.closeLHEF(); }, "");
		cl.def("closeLHEF", (bool (Pythia8::LHAup::*)(bool)) &Pythia8::LHAup::closeLHEF, "C++: Pythia8::LHAup::closeLHEF(bool) --> bool", pybind11::arg("updateInit"));
		cl.def("initLHEF", (bool (Pythia8::LHAup::*)()) &Pythia8::LHAup::initLHEF, "C++: Pythia8::LHAup::initLHEF() --> bool");
		cl.def("eventLHEF", [](Pythia8::LHAup &o) -> bool { return o.eventLHEF(); }, "");
		cl.def("eventLHEF", (bool (Pythia8::LHAup::*)(bool)) &Pythia8::LHAup::eventLHEF, "C++: Pythia8::LHAup::eventLHEF(bool) --> bool", pybind11::arg("verbose"));
		cl.def("getFileName", (std::string (Pythia8::LHAup::*)() const) &Pythia8::LHAup::getFileName, "C++: Pythia8::LHAup::getFileName() const --> std::string");
		cl.def("setBeamA", [](Pythia8::LHAup &o, int const & a0, double const & a1) -> void { return o.setBeamA(a0, a1); }, "", pybind11::arg("idIn"), pybind11::arg("eIn"));
		cl.def("setBeamA", [](Pythia8::LHAup &o, int const & a0, double const & a1, int const & a2) -> void { return o.setBeamA(a0, a1, a2); }, "", pybind11::arg("idIn"), pybind11::arg("eIn"), pybind11::arg("pdfGroupIn"));
		cl.def("setBeamA", (void (Pythia8::LHAup::*)(int, double, int, int)) &Pythia8::LHAup::setBeamA, "C++: Pythia8::LHAup::setBeamA(int, double, int, int) --> void", pybind11::arg("idIn"), pybind11::arg("eIn"), pybind11::arg("pdfGroupIn"), pybind11::arg("pdfSetIn"));
		cl.def("setBeamB", [](Pythia8::LHAup &o, int const & a0, double const & a1) -> void { return o.setBeamB(a0, a1); }, "", pybind11::arg("idIn"), pybind11::arg("eIn"));
		cl.def("setBeamB", [](Pythia8::LHAup &o, int const & a0, double const & a1, int const & a2) -> void { return o.setBeamB(a0, a1, a2); }, "", pybind11::arg("idIn"), pybind11::arg("eIn"), pybind11::arg("pdfGroupIn"));
		cl.def("setBeamB", (void (Pythia8::LHAup::*)(int, double, int, int)) &Pythia8::LHAup::setBeamB, "C++: Pythia8::LHAup::setBeamB(int, double, int, int) --> void", pybind11::arg("idIn"), pybind11::arg("eIn"), pybind11::arg("pdfGroupIn"), pybind11::arg("pdfSetIn"));
		cl.def("setStrategy", (void (Pythia8::LHAup::*)(int)) &Pythia8::LHAup::setStrategy, "C++: Pythia8::LHAup::setStrategy(int) --> void", pybind11::arg("strategyIn"));
		cl.def("addProcess", [](Pythia8::LHAup &o, int const & a0) -> void { return o.addProcess(a0); }, "", pybind11::arg("idProcIn"));
		cl.def("addProcess", [](Pythia8::LHAup &o, int const & a0, double const & a1) -> void { return o.addProcess(a0, a1); }, "", pybind11::arg("idProcIn"), pybind11::arg("xSecIn"));
		cl.def("addProcess", [](Pythia8::LHAup &o, int const & a0, double const & a1, double const & a2) -> void { return o.addProcess(a0, a1, a2); }, "", pybind11::arg("idProcIn"), pybind11::arg("xSecIn"), pybind11::arg("xErrIn"));
		cl.def("addProcess", (void (Pythia8::LHAup::*)(int, double, double, double)) &Pythia8::LHAup::addProcess, "C++: Pythia8::LHAup::addProcess(int, double, double, double) --> void", pybind11::arg("idProcIn"), pybind11::arg("xSecIn"), pybind11::arg("xErrIn"), pybind11::arg("xMaxIn"));
		cl.def("setXSec", (void (Pythia8::LHAup::*)(int, double)) &Pythia8::LHAup::setXSec, "C++: Pythia8::LHAup::setXSec(int, double) --> void", pybind11::arg("iP"), pybind11::arg("xSecIn"));
		cl.def("setXErr", (void (Pythia8::LHAup::*)(int, double)) &Pythia8::LHAup::setXErr, "C++: Pythia8::LHAup::setXErr(int, double) --> void", pybind11::arg("iP"), pybind11::arg("xErrIn"));
		cl.def("setXMax", (void (Pythia8::LHAup::*)(int, double)) &Pythia8::LHAup::setXMax, "C++: Pythia8::LHAup::setXMax(int, double) --> void", pybind11::arg("iP"), pybind11::arg("xMaxIn"));
		cl.def("setProcess", [](Pythia8::LHAup &o) -> void { return o.setProcess(); }, "");
		cl.def("setProcess", [](Pythia8::LHAup &o, int const & a0) -> void { return o.setProcess(a0); }, "", pybind11::arg("idProcIn"));
		cl.def("setProcess", [](Pythia8::LHAup &o, int const & a0, double const & a1) -> void { return o.setProcess(a0, a1); }, "", pybind11::arg("idProcIn"), pybind11::arg("weightIn"));
		cl.def("setProcess", [](Pythia8::LHAup &o, int const & a0, double const & a1, double const & a2) -> void { return o.setProcess(a0, a1, a2); }, "", pybind11::arg("idProcIn"), pybind11::arg("weightIn"), pybind11::arg("scaleIn"));
		cl.def("setProcess", [](Pythia8::LHAup &o, int const & a0, double const & a1, double const & a2, double const & a3) -> void { return o.setProcess(a0, a1, a2, a3); }, "", pybind11::arg("idProcIn"), pybind11::arg("weightIn"), pybind11::arg("scaleIn"), pybind11::arg("alphaQEDIn"));
		cl.def("setProcess", (void (Pythia8::LHAup::*)(int, double, double, double, double)) &Pythia8::LHAup::setProcess, "C++: Pythia8::LHAup::setProcess(int, double, double, double, double) --> void", pybind11::arg("idProcIn"), pybind11::arg("weightIn"), pybind11::arg("scaleIn"), pybind11::arg("alphaQEDIn"), pybind11::arg("alphaQCDIn"));
		cl.def("addParticle", (void (Pythia8::LHAup::*)(class Pythia8::LHAParticle)) &Pythia8::LHAup::addParticle, "C++: Pythia8::LHAup::addParticle(class Pythia8::LHAParticle) --> void", pybind11::arg("particleIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0) -> void { return o.addParticle(a0); }, "", pybind11::arg("idIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1) -> void { return o.addParticle(a0, a1); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2) -> void { return o.addParticle(a0, a1, a2); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3) -> void { return o.addParticle(a0, a1, a2, a3); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4) -> void { return o.addParticle(a0, a1, a2, a3, a4); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5) -> void { return o.addParticle(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, double const & a6) -> void { return o.addParticle(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, double const & a6, double const & a7) -> void { return o.addParticle(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"), pybind11::arg("pyIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, double const & a6, double const & a7, double const & a8) -> void { return o.addParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("pzIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, double const & a6, double const & a7, double const & a8, double const & a9) -> void { return o.addParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("pzIn"), pybind11::arg("eIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, double const & a6, double const & a7, double const & a8, double const & a9, double const & a10) -> void { return o.addParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("pzIn"), pybind11::arg("eIn"), pybind11::arg("mIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, double const & a6, double const & a7, double const & a8, double const & a9, double const & a10, double const & a11) -> void { return o.addParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("pzIn"), pybind11::arg("eIn"), pybind11::arg("mIn"), pybind11::arg("tauIn"));
		cl.def("addParticle", [](Pythia8::LHAup &o, int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, double const & a6, double const & a7, double const & a8, double const & a9, double const & a10, double const & a11, double const & a12) -> void { return o.addParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("pzIn"), pybind11::arg("eIn"), pybind11::arg("mIn"), pybind11::arg("tauIn"), pybind11::arg("spinIn"));
		cl.def("addParticle", (void (Pythia8::LHAup::*)(int, int, int, int, int, int, double, double, double, double, double, double, double, double)) &Pythia8::LHAup::addParticle, "C++: Pythia8::LHAup::addParticle(int, int, int, int, int, int, double, double, double, double, double, double, double, double) --> void", pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("col1In"), pybind11::arg("col2In"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("pzIn"), pybind11::arg("eIn"), pybind11::arg("mIn"), pybind11::arg("tauIn"), pybind11::arg("spinIn"), pybind11::arg("scaleIn"));
		cl.def("setIdX", (void (Pythia8::LHAup::*)(int, int, double, double)) &Pythia8::LHAup::setIdX, "C++: Pythia8::LHAup::setIdX(int, int, double, double) --> void", pybind11::arg("id1In"), pybind11::arg("id2In"), pybind11::arg("x1In"), pybind11::arg("x2In"));
		cl.def("setPdf", (void (Pythia8::LHAup::*)(int, int, double, double, double, double, double, bool)) &Pythia8::LHAup::setPdf, "C++: Pythia8::LHAup::setPdf(int, int, double, double, double, double, double, bool) --> void", pybind11::arg("id1pdfIn"), pybind11::arg("id2pdfIn"), pybind11::arg("x1pdfIn"), pybind11::arg("x2pdfIn"), pybind11::arg("scalePDFIn"), pybind11::arg("pdf1In"), pybind11::arg("pdf2In"), pybind11::arg("pdfIsSetIn"));
		cl.def("setScaleShowers", [](Pythia8::LHAup &o, double const & a0) -> void { return o.setScaleShowers(a0); }, "", pybind11::arg("scaleIn1"));
		cl.def("setScaleShowers", (void (Pythia8::LHAup::*)(double, double)) &Pythia8::LHAup::setScaleShowers, "C++: Pythia8::LHAup::setScaleShowers(double, double) --> void", pybind11::arg("scaleIn1"), pybind11::arg("scaleIn2"));
		cl.def("setInitLHEF", [](Pythia8::LHAup &o, class std::basic_istream<char> & a0) -> bool { return o.setInitLHEF(a0); }, "", pybind11::arg("is"));
		cl.def("setInitLHEF", (bool (Pythia8::LHAup::*)(class std::basic_istream<char> &, bool)) &Pythia8::LHAup::setInitLHEF, "C++: Pythia8::LHAup::setInitLHEF(class std::basic_istream<char> &, bool) --> bool", pybind11::arg("is"), pybind11::arg("readHeaders"));
		cl.def("setNewEventLHEF", (bool (Pythia8::LHAup::*)(class std::basic_istream<char> &)) &Pythia8::LHAup::setNewEventLHEF, "C++: Pythia8::LHAup::setNewEventLHEF(class std::basic_istream<char> &) --> bool", pybind11::arg("is"));
		cl.def("setOldEventLHEF", (bool (Pythia8::LHAup::*)()) &Pythia8::LHAup::setOldEventLHEF, "C++: Pythia8::LHAup::setOldEventLHEF() --> bool");
		cl.def("openFile", (class std::basic_istream<char> * (Pythia8::LHAup::*)(const char *, class std::basic_ifstream<char> &)) &Pythia8::LHAup::openFile, "C++: Pythia8::LHAup::openFile(const char *, class std::basic_ifstream<char> &) --> class std::basic_istream<char> *", pybind11::return_value_policy::automatic, pybind11::arg("fn"), pybind11::arg("ifs"));
		cl.def("setInfoHeader", (void (Pythia8::LHAup::*)(const std::string &, const std::string &)) &Pythia8::LHAup::setInfoHeader, "C++: Pythia8::LHAup::setInfoHeader(const std::string &, const std::string &) --> void", pybind11::arg("key"), pybind11::arg("val"));
	}
}
