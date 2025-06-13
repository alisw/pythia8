#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PhysicsBase.h>
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

// Pythia8::BeamParticle file:Pythia8/BeamParticle.h line:133
struct PyCallBack_Pythia8_BeamParticle : public Pythia8::BeamParticle {
	using Pythia8::BeamParticle::BeamParticle;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamParticle *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamParticle *>(this), "onBeginEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onBeginEvent();
	}
	void onEndEvent(enum Pythia8::PhysicsBase::Status a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamParticle *>(this), "onEndEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onEndEvent(a0);
	}
	void onStat() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamParticle *>(this), "onStat");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onStat();
	}
};

void bind_Pythia8_BeamParticle(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::BeamParticle file:Pythia8/BeamParticle.h line:133
		pybind11::class_<Pythia8::BeamParticle, std::shared_ptr<Pythia8::BeamParticle>, PyCallBack_Pythia8_BeamParticle, Pythia8::PhysicsBase> cl(M("Pythia8"), "BeamParticle", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::BeamParticle(); }, [](){ return new PyCallBack_Pythia8_BeamParticle(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_BeamParticle const &o){ return new PyCallBack_Pythia8_BeamParticle(o); } ) );
		cl.def( pybind11::init( [](Pythia8::BeamParticle const &o){ return new Pythia8::BeamParticle(o); } ) );
		cl.def("init", (void (Pythia8::BeamParticle::*)(int, double, double, double, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, bool, class Pythia8::StringFlav *)) &Pythia8::BeamParticle::init, "C++: Pythia8::BeamParticle::init(int, double, double, double, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, bool, class Pythia8::StringFlav *) --> void", pybind11::arg("idIn"), pybind11::arg("pzIn"), pybind11::arg("eIn"), pybind11::arg("mIn"), pybind11::arg("pdfInPtr"), pybind11::arg("pdfHardInPtr"), pybind11::arg("isUnresolvedIn"), pybind11::arg("flavSelPtrIn"));
		cl.def("initID", (void (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::initID, "C++: Pythia8::BeamParticle::initID(int) --> void", pybind11::arg("idIn"));
		cl.def("initPDFPtr", (void (Pythia8::BeamParticle::*)(class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>)) &Pythia8::BeamParticle::initPDFPtr, "C++: Pythia8::BeamParticle::initPDFPtr(class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>) --> void", pybind11::arg("pdfInPtr"), pybind11::arg("pdfHardInPtr"));
		cl.def("initUnres", (void (Pythia8::BeamParticle::*)(class std::shared_ptr<class Pythia8::PDF>)) &Pythia8::BeamParticle::initUnres, "C++: Pythia8::BeamParticle::initUnres(class std::shared_ptr<class Pythia8::PDF>) --> void", pybind11::arg("pdfUnresInPtr"));
		cl.def("initSwitchID", (void (Pythia8::BeamParticle::*)(const class std::vector<class std::shared_ptr<class Pythia8::PDF>, class std::allocator<class std::shared_ptr<class Pythia8::PDF> > > &)) &Pythia8::BeamParticle::initSwitchID, "C++: Pythia8::BeamParticle::initSwitchID(const class std::vector<class std::shared_ptr<class Pythia8::PDF>, class std::allocator<class std::shared_ptr<class Pythia8::PDF> > > &) --> void", pybind11::arg("pdfSavePtrsIn"));
		cl.def("newValenceContent", (void (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::newValenceContent, "C++: Pythia8::BeamParticle::newValenceContent() --> void");
		cl.def("setValenceContent", [](Pythia8::BeamParticle &o, int const & a0) -> void { return o.setValenceContent(a0); }, "", pybind11::arg("idq1"));
		cl.def("setValenceContent", [](Pythia8::BeamParticle &o, int const & a0, int const & a1) -> void { return o.setValenceContent(a0, a1); }, "", pybind11::arg("idq1"), pybind11::arg("idq2"));
		cl.def("setValenceContent", (void (Pythia8::BeamParticle::*)(int, int, int)) &Pythia8::BeamParticle::setValenceContent, "C++: Pythia8::BeamParticle::setValenceContent(int, int, int) --> void", pybind11::arg("idq1"), pybind11::arg("idq2"), pybind11::arg("idq3"));
		cl.def("setBeamID", [](Pythia8::BeamParticle &o, int const & a0) -> void { return o.setBeamID(a0); }, "", pybind11::arg("idIn"));
		cl.def("setBeamID", (void (Pythia8::BeamParticle::*)(int, int)) &Pythia8::BeamParticle::setBeamID, "C++: Pythia8::BeamParticle::setBeamID(int, int) --> void", pybind11::arg("idIn"), pybind11::arg("iPDFin"));
		cl.def("newPzE", (void (Pythia8::BeamParticle::*)(double, double)) &Pythia8::BeamParticle::newPzE, "C++: Pythia8::BeamParticle::newPzE(double, double) --> void", pybind11::arg("pzIn"), pybind11::arg("eIn"));
		cl.def("newM", (void (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::newM, "C++: Pythia8::BeamParticle::newM(double) --> void", pybind11::arg("mIn"));
		cl.def("id", (int (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::id, "C++: Pythia8::BeamParticle::id() const --> int");
		cl.def("idVMD", (int (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::idVMD, "C++: Pythia8::BeamParticle::idVMD() const --> int");
		cl.def("p", (class Pythia8::Vec4 (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::p, "C++: Pythia8::BeamParticle::p() const --> class Pythia8::Vec4");
		cl.def("px", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::px, "C++: Pythia8::BeamParticle::px() const --> double");
		cl.def("py", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::py, "C++: Pythia8::BeamParticle::py() const --> double");
		cl.def("pz", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::pz, "C++: Pythia8::BeamParticle::pz() const --> double");
		cl.def("e", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::e, "C++: Pythia8::BeamParticle::e() const --> double");
		cl.def("m", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::m, "C++: Pythia8::BeamParticle::m() const --> double");
		cl.def("mVMD", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::mVMD, "C++: Pythia8::BeamParticle::mVMD() const --> double");
		cl.def("scaleVMD", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::scaleVMD, "C++: Pythia8::BeamParticle::scaleVMD() const --> double");
		cl.def("isLepton", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::isLepton, "C++: Pythia8::BeamParticle::isLepton() const --> bool");
		cl.def("isUnresolved", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::isUnresolved, "C++: Pythia8::BeamParticle::isUnresolved() const --> bool");
		cl.def("isHadron", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::isHadron, "C++: Pythia8::BeamParticle::isHadron() const --> bool");
		cl.def("isMeson", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::isMeson, "C++: Pythia8::BeamParticle::isMeson() const --> bool");
		cl.def("isBaryon", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::isBaryon, "C++: Pythia8::BeamParticle::isBaryon() const --> bool");
		cl.def("isGamma", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::isGamma, "C++: Pythia8::BeamParticle::isGamma() const --> bool");
		cl.def("hasResGamma", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::hasResGamma, "C++: Pythia8::BeamParticle::hasResGamma() const --> bool");
		cl.def("hasVMDstate", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::hasVMDstate, "C++: Pythia8::BeamParticle::hasVMDstate() const --> bool");
		cl.def("xMax", [](Pythia8::BeamParticle &o) -> double { return o.xMax(); }, "");
		cl.def("xMax", (double (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::xMax, "C++: Pythia8::BeamParticle::xMax(int) --> double", pybind11::arg("iSkip"));
		cl.def("xfHard", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfHard, "C++: Pythia8::BeamParticle::xfHard(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfMax", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfMax, "C++: Pythia8::BeamParticle::xfMax(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfFlux", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfFlux, "C++: Pythia8::BeamParticle::xfFlux(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfApprox", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfApprox, "C++: Pythia8::BeamParticle::xfApprox(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfGamma", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfGamma, "C++: Pythia8::BeamParticle::xfGamma(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfSame", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfSame, "C++: Pythia8::BeamParticle::xfSame(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xf", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xf, "C++: Pythia8::BeamParticle::xf(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfVal", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfVal, "C++: Pythia8::BeamParticle::xfVal(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfSea", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfSea, "C++: Pythia8::BeamParticle::xfSea(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfMPI", (double (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::xfMPI, "C++: Pythia8::BeamParticle::xfMPI(int, double, double) --> double", pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfISR", (double (Pythia8::BeamParticle::*)(int, int, double, double)) &Pythia8::BeamParticle::xfISR, "C++: Pythia8::BeamParticle::xfISR(int, int, double, double) --> double", pybind11::arg("indexMPI"), pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("insideBounds", (bool (Pythia8::BeamParticle::*)(double, double)) &Pythia8::BeamParticle::insideBounds, "C++: Pythia8::BeamParticle::insideBounds(double, double) --> bool", pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("alphaS", (double (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::alphaS, "C++: Pythia8::BeamParticle::alphaS(double) --> double", pybind11::arg("Q2"));
		cl.def("mQuarkPDF", (double (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::mQuarkPDF, "C++: Pythia8::BeamParticle::mQuarkPDF(int) --> double", pybind11::arg("idIn"));
		cl.def("nMembers", (int (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::nMembers, "C++: Pythia8::BeamParticle::nMembers() --> int");
		cl.def("calcPDFEnvelope", (void (Pythia8::BeamParticle::*)(int, double, double, int)) &Pythia8::BeamParticle::calcPDFEnvelope, "C++: Pythia8::BeamParticle::calcPDFEnvelope(int, double, double, int) --> void", pybind11::arg("idNow"), pybind11::arg("xNow"), pybind11::arg("Q2Now"), pybind11::arg("valSea"));
		cl.def("calcPDFEnvelope", (void (Pythia8::BeamParticle::*)(struct std::pair<int, int>, struct std::pair<double, double>, double, int)) &Pythia8::BeamParticle::calcPDFEnvelope, "C++: Pythia8::BeamParticle::calcPDFEnvelope(struct std::pair<int, int>, struct std::pair<double, double>, double, int) --> void", pybind11::arg("idNows"), pybind11::arg("xNows"), pybind11::arg("Q2Now"), pybind11::arg("valSea"));
		cl.def("getPDFEnvelope", (struct Pythia8::PDF::PDFEnvelope (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::getPDFEnvelope, "C++: Pythia8::BeamParticle::getPDFEnvelope() --> struct Pythia8::PDF::PDFEnvelope");
		cl.def("pickValSeaComp", (int (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::pickValSeaComp, "C++: Pythia8::BeamParticle::pickValSeaComp() --> int");
		cl.def("initBeamKind", (void (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::initBeamKind, "C++: Pythia8::BeamParticle::initBeamKind() --> void");
		cl.def("__getitem__", (class Pythia8::ResolvedParton & (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::operator[], "C++: Pythia8::BeamParticle::operator[](int) --> class Pythia8::ResolvedParton &", pybind11::return_value_policy::reference, pybind11::arg("i"));
		cl.def("at", (class Pythia8::ResolvedParton & (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::at, "C++: Pythia8::BeamParticle::at(int) --> class Pythia8::ResolvedParton &", pybind11::return_value_policy::reference, pybind11::arg("i"));
		cl.def("size", (int (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::size, "C++: Pythia8::BeamParticle::size() const --> int");
		cl.def("sizeInit", (int (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::sizeInit, "C++: Pythia8::BeamParticle::sizeInit() const --> int");
		cl.def("clear", (void (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::clear, "C++: Pythia8::BeamParticle::clear() --> void");
		cl.def("resetGamma", (void (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::resetGamma, "C++: Pythia8::BeamParticle::resetGamma() --> void");
		cl.def("resetGammaInLepton", (void (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::resetGammaInLepton, "C++: Pythia8::BeamParticle::resetGammaInLepton() --> void");
		cl.def("append", [](Pythia8::BeamParticle &o, int const & a0, int const & a1, double const & a2) -> int { return o.append(a0, a1, a2); }, "", pybind11::arg("iPos"), pybind11::arg("idIn"), pybind11::arg("x"));
		cl.def("append", (int (Pythia8::BeamParticle::*)(int, int, double, int)) &Pythia8::BeamParticle::append, "C++: Pythia8::BeamParticle::append(int, int, double, int) --> int", pybind11::arg("iPos"), pybind11::arg("idIn"), pybind11::arg("x"), pybind11::arg("companion"));
		cl.def("popBack", (void (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::popBack, "C++: Pythia8::BeamParticle::popBack() --> void");
		cl.def("list", (void (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::list, "C++: Pythia8::BeamParticle::list() const --> void");
		cl.def("nValenceKinds", (int (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::nValenceKinds, "C++: Pythia8::BeamParticle::nValenceKinds() const --> int");
		cl.def("nValence", (int (Pythia8::BeamParticle::*)(int) const) &Pythia8::BeamParticle::nValence, "C++: Pythia8::BeamParticle::nValence(int) const --> int", pybind11::arg("idIn"));
		cl.def("isUnresolvedLepton", (bool (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::isUnresolvedLepton, "C++: Pythia8::BeamParticle::isUnresolvedLepton() --> bool");
		cl.def("remnantFlavours", [](Pythia8::BeamParticle &o, class Pythia8::Event & a0) -> bool { return o.remnantFlavours(a0); }, "", pybind11::arg("event"));
		cl.def("remnantFlavours", (bool (Pythia8::BeamParticle::*)(class Pythia8::Event &, bool)) &Pythia8::BeamParticle::remnantFlavours, "C++: Pythia8::BeamParticle::remnantFlavours(class Pythia8::Event &, bool) --> bool", pybind11::arg("event"), pybind11::arg("isDIS"));
		cl.def("remnantColours", (bool (Pythia8::BeamParticle::*)(class Pythia8::Event &, class std::vector<int, class std::allocator<int> > &, class std::vector<int, class std::allocator<int> > &)) &Pythia8::BeamParticle::remnantColours, "C++: Pythia8::BeamParticle::remnantColours(class Pythia8::Event &, class std::vector<int, class std::allocator<int> > &, class std::vector<int, class std::allocator<int> > &) --> bool", pybind11::arg("event"), pybind11::arg("colFrom"), pybind11::arg("colTo"));
		cl.def("xRemnant", (double (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::xRemnant, "C++: Pythia8::BeamParticle::xRemnant(int) --> double", pybind11::arg("i"));
		cl.def("hasJunction", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::hasJunction, "C++: Pythia8::BeamParticle::hasJunction() const --> bool");
		cl.def("junctionCol", (int (Pythia8::BeamParticle::*)(int) const) &Pythia8::BeamParticle::junctionCol, "C++: Pythia8::BeamParticle::junctionCol(int) const --> int", pybind11::arg("i"));
		cl.def("junctionCol", (void (Pythia8::BeamParticle::*)(int, int)) &Pythia8::BeamParticle::junctionCol, "C++: Pythia8::BeamParticle::junctionCol(int, int) --> void", pybind11::arg("i"), pybind11::arg("col"));
		cl.def("pickGluon", (bool (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::pickGluon, "C++: Pythia8::BeamParticle::pickGluon(double) --> bool", pybind11::arg("mDiff"));
		cl.def("pickValence", (int (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::pickValence, "C++: Pythia8::BeamParticle::pickValence() --> int");
		cl.def("pickRemnant", (int (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::pickRemnant, "C++: Pythia8::BeamParticle::pickRemnant() const --> int");
		cl.def("zShare", (double (Pythia8::BeamParticle::*)(double, double, double)) &Pythia8::BeamParticle::zShare, "C++: Pythia8::BeamParticle::zShare(double, double, double) --> double", pybind11::arg("mDiff"), pybind11::arg("m1"), pybind11::arg("m2"));
		cl.def("pxShare", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::pxShare, "C++: Pythia8::BeamParticle::pxShare() const --> double");
		cl.def("pyShare", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::pyShare, "C++: Pythia8::BeamParticle::pyShare() const --> double");
		cl.def("remnantFlavoursNew", (bool (Pythia8::BeamParticle::*)(class Pythia8::Event &)) &Pythia8::BeamParticle::remnantFlavoursNew, "C++: Pythia8::BeamParticle::remnantFlavoursNew(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("findColSetup", (void (Pythia8::BeamParticle::*)(class Pythia8::Event &)) &Pythia8::BeamParticle::findColSetup, "C++: Pythia8::BeamParticle::findColSetup(class Pythia8::Event &) --> void", pybind11::arg("event"));
		cl.def("setInitialCol", (void (Pythia8::BeamParticle::*)(class Pythia8::Event &)) &Pythia8::BeamParticle::setInitialCol, "C++: Pythia8::BeamParticle::setInitialCol(class Pythia8::Event &) --> void", pybind11::arg("event"));
		cl.def("updateCol", (void (Pythia8::BeamParticle::*)(class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > >)) &Pythia8::BeamParticle::updateCol, "C++: Pythia8::BeamParticle::updateCol(class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > >) --> void", pybind11::arg("colourChanges"));
		cl.def("getColUpdates", (class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::getColUpdates, "C++: Pythia8::BeamParticle::getColUpdates() --> class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > >");
		cl.def("gammaInitiatorIsVal", (bool (Pythia8::BeamParticle::*)(int, int, double, double)) &Pythia8::BeamParticle::gammaInitiatorIsVal, "C++: Pythia8::BeamParticle::gammaInitiatorIsVal(int, int, double, double) --> bool", pybind11::arg("iResolved"), pybind11::arg("id"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("gammaInitiatorIsVal", (bool (Pythia8::BeamParticle::*)(int, double)) &Pythia8::BeamParticle::gammaInitiatorIsVal, "C++: Pythia8::BeamParticle::gammaInitiatorIsVal(int, double) --> bool", pybind11::arg("iResolved"), pybind11::arg("Q2"));
		cl.def("getGammaValFlavour", (int (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::getGammaValFlavour, "C++: Pythia8::BeamParticle::getGammaValFlavour() --> int");
		cl.def("gammaValSeaComp", (int (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::gammaValSeaComp, "C++: Pythia8::BeamParticle::gammaValSeaComp(int) --> int", pybind11::arg("iResolved"));
		cl.def("posVal", (void (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::posVal, "C++: Pythia8::BeamParticle::posVal(int) --> void", pybind11::arg("iPosValIn"));
		cl.def("gamVal", (void (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::gamVal, "C++: Pythia8::BeamParticle::gamVal(int) --> void", pybind11::arg("iGamValIn"));
		cl.def("gamVal", (int (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::gamVal, "C++: Pythia8::BeamParticle::gamVal() --> int");
		cl.def("resolvedGamma", (void (Pythia8::BeamParticle::*)(bool)) &Pythia8::BeamParticle::resolvedGamma, "C++: Pythia8::BeamParticle::resolvedGamma(bool) --> void", pybind11::arg("isResolved"));
		cl.def("resolvedGamma", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::resolvedGamma, "C++: Pythia8::BeamParticle::resolvedGamma() const --> bool");
		cl.def("setGammaMode", (void (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::setGammaMode, "C++: Pythia8::BeamParticle::setGammaMode(int) --> void", pybind11::arg("gammaModeIn"));
		cl.def("getGammaMode", (int (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::getGammaMode, "C++: Pythia8::BeamParticle::getGammaMode() const --> int");
		cl.def("isResolvedUnresolved", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::isResolvedUnresolved, "C++: Pythia8::BeamParticle::isResolvedUnresolved() const --> bool");
		cl.def("initGammaInBeam", (void (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::initGammaInBeam, "C++: Pythia8::BeamParticle::initGammaInBeam() --> void");
		cl.def("gammaInBeam", (bool (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::gammaInBeam, "C++: Pythia8::BeamParticle::gammaInBeam() const --> bool");
		cl.def("setVMDstate", [](Pythia8::BeamParticle &o, bool const & a0, int const & a1, double const & a2, double const & a3) -> void { return o.setVMDstate(a0, a1, a2, a3); }, "", pybind11::arg("isVMDIn"), pybind11::arg("idIn"), pybind11::arg("mIn"), pybind11::arg("scaleIn"));
		cl.def("setVMDstate", (void (Pythia8::BeamParticle::*)(bool, int, double, double, bool)) &Pythia8::BeamParticle::setVMDstate, "C++: Pythia8::BeamParticle::setVMDstate(bool, int, double, double, bool) --> void", pybind11::arg("isVMDIn"), pybind11::arg("idIn"), pybind11::arg("mIn"), pybind11::arg("scaleIn"), pybind11::arg("reassignState"));
		cl.def("pT2gamma2qqbar", (void (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::pT2gamma2qqbar, "C++: Pythia8::BeamParticle::pT2gamma2qqbar(double) --> void", pybind11::arg("pT2in"));
		cl.def("pT2gamma2qqbar", (double (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::pT2gamma2qqbar, "C++: Pythia8::BeamParticle::pT2gamma2qqbar() --> double");
		cl.def("pTMPI", (void (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::pTMPI, "C++: Pythia8::BeamParticle::pTMPI(double) --> void", pybind11::arg("pTminMPIin"));
		cl.def("roomFor1Remnant", (bool (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::roomFor1Remnant, "C++: Pythia8::BeamParticle::roomFor1Remnant(double) --> bool", pybind11::arg("eCM"));
		cl.def("roomFor1Remnant", (bool (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::roomFor1Remnant, "C++: Pythia8::BeamParticle::roomFor1Remnant(int, double, double) --> bool", pybind11::arg("id1"), pybind11::arg("x1"), pybind11::arg("eCM"));
		cl.def("roomFor2Remnants", (bool (Pythia8::BeamParticle::*)(int, double, double)) &Pythia8::BeamParticle::roomFor2Remnants, "C++: Pythia8::BeamParticle::roomFor2Remnants(int, double, double) --> bool", pybind11::arg("id1"), pybind11::arg("x1"), pybind11::arg("eCM"));
		cl.def("roomForRemnants", (bool (Pythia8::BeamParticle::*)(class Pythia8::BeamParticle)) &Pythia8::BeamParticle::roomForRemnants, "C++: Pythia8::BeamParticle::roomForRemnants(class Pythia8::BeamParticle) --> bool", pybind11::arg("beamOther"));
		cl.def("remnantMass", (double (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::remnantMass, "C++: Pythia8::BeamParticle::remnantMass(int) --> double", pybind11::arg("idIn"));
		cl.def("gammaPDFxDependence", (double (Pythia8::BeamParticle::*)(int, double)) &Pythia8::BeamParticle::gammaPDFxDependence, "C++: Pythia8::BeamParticle::gammaPDFxDependence(int, double) --> double", pybind11::arg("flavour"), pybind11::arg("x"));
		cl.def("gammaPDFRefScale", (double (Pythia8::BeamParticle::*)(int)) &Pythia8::BeamParticle::gammaPDFRefScale, "C++: Pythia8::BeamParticle::gammaPDFRefScale(int) --> double", pybind11::arg("flavour"));
		cl.def("xIntegratedPDFs", (double (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::xIntegratedPDFs, "C++: Pythia8::BeamParticle::xIntegratedPDFs(double) --> double", pybind11::arg("Q2"));
		cl.def("xGammaPDF", (void (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::xGammaPDF, "C++: Pythia8::BeamParticle::xGammaPDF() --> void");
		cl.def("xGamma", (void (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::xGamma, "C++: Pythia8::BeamParticle::xGamma(double) --> void", pybind11::arg("xGmIn"));
		cl.def("Q2Gamma", (void (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::Q2Gamma, "C++: Pythia8::BeamParticle::Q2Gamma(double) --> void", pybind11::arg("Q2GmIn"));
		cl.def("newGammaKTPhi", (void (Pythia8::BeamParticle::*)(double, double)) &Pythia8::BeamParticle::newGammaKTPhi, "C++: Pythia8::BeamParticle::newGammaKTPhi(double, double) --> void", pybind11::arg("kTIn"), pybind11::arg("phiIn"));
		cl.def("xGammaMin", (double (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::xGammaMin, "C++: Pythia8::BeamParticle::xGammaMin() --> double");
		cl.def("xGammaHadr", (double (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::xGammaHadr, "C++: Pythia8::BeamParticle::xGammaHadr() --> double");
		cl.def("gammaFluxIntApprox", (double (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::gammaFluxIntApprox, "C++: Pythia8::BeamParticle::gammaFluxIntApprox() --> double");
		cl.def("hasApproxGammaFlux", (bool (Pythia8::BeamParticle::*)()) &Pythia8::BeamParticle::hasApproxGammaFlux, "C++: Pythia8::BeamParticle::hasApproxGammaFlux() --> bool");
		cl.def("xGamma", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::xGamma, "C++: Pythia8::BeamParticle::xGamma() const --> double");
		cl.def("Q2Gamma", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::Q2Gamma, "C++: Pythia8::BeamParticle::Q2Gamma() const --> double");
		cl.def("gammaKTx", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::gammaKTx, "C++: Pythia8::BeamParticle::gammaKTx() const --> double");
		cl.def("gammaKTy", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::gammaKTy, "C++: Pythia8::BeamParticle::gammaKTy() const --> double");
		cl.def("gammaKT", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::gammaKT, "C++: Pythia8::BeamParticle::gammaKT() const --> double");
		cl.def("gammaPhi", (double (Pythia8::BeamParticle::*)() const) &Pythia8::BeamParticle::gammaPhi, "C++: Pythia8::BeamParticle::gammaPhi() const --> double");
		cl.def("xPom", [](Pythia8::BeamParticle &o) -> void { return o.xPom(); }, "");
		cl.def("xPom", (void (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::xPom, "C++: Pythia8::BeamParticle::xPom(double) --> void", pybind11::arg("xpom"));
		cl.def("sampleXgamma", (double (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::sampleXgamma, "C++: Pythia8::BeamParticle::sampleXgamma(double) --> double", pybind11::arg("xMinIn"));
		cl.def("sampleQ2gamma", (double (Pythia8::BeamParticle::*)(double)) &Pythia8::BeamParticle::sampleQ2gamma, "C++: Pythia8::BeamParticle::sampleQ2gamma(double) --> double", pybind11::arg("Q2min"));
		cl.def("assign", (class Pythia8::BeamParticle & (Pythia8::BeamParticle::*)(const class Pythia8::BeamParticle &)) &Pythia8::BeamParticle::operator=, "C++: Pythia8::BeamParticle::operator=(const class Pythia8::BeamParticle &) --> class Pythia8::BeamParticle &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
