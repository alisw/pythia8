#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/LesHouches.h>
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

// Pythia8::BeamSetup file:Pythia8/BeamSetup.h line:33
struct PyCallBack_Pythia8_BeamSetup : public Pythia8::BeamSetup {
	using Pythia8::BeamSetup::BeamSetup;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamSetup *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BeamSetup::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamSetup *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamSetup *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamSetup *>(this), "onStat");
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

void bind_Pythia8_BeamSetup(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::BeamSetup file:Pythia8/BeamSetup.h line:33
		pybind11::class_<Pythia8::BeamSetup, std::shared_ptr<Pythia8::BeamSetup>, PyCallBack_Pythia8_BeamSetup, Pythia8::PhysicsBase> cl(M("Pythia8"), "BeamSetup", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::BeamSetup(); }, [](){ return new PyCallBack_Pythia8_BeamSetup(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_BeamSetup const &o){ return new PyCallBack_Pythia8_BeamSetup(o); } ) );
		cl.def( pybind11::init( [](Pythia8::BeamSetup const &o){ return new Pythia8::BeamSetup(o); } ) );
		cl.def_readwrite("doLHA", &Pythia8::BeamSetup::doLHA);
		cl.def_readwrite("useNewLHA", &Pythia8::BeamSetup::useNewLHA);
		cl.def_readwrite("skipInit", &Pythia8::BeamSetup::skipInit);
		cl.def_readwrite("doMomentumSpread", &Pythia8::BeamSetup::doMomentumSpread);
		cl.def_readwrite("doVertexSpread", &Pythia8::BeamSetup::doVertexSpread);
		cl.def_readwrite("doVarEcm", &Pythia8::BeamSetup::doVarEcm);
		cl.def_readwrite("allowIDAswitch", &Pythia8::BeamSetup::allowIDAswitch);
		cl.def_readwrite("hasSwitchedIDs", &Pythia8::BeamSetup::hasSwitchedIDs);
		cl.def_readwrite("beamA2gamma", &Pythia8::BeamSetup::beamA2gamma);
		cl.def_readwrite("beamB2gamma", &Pythia8::BeamSetup::beamB2gamma);
		cl.def_readwrite("idA", &Pythia8::BeamSetup::idA);
		cl.def_readwrite("idB", &Pythia8::BeamSetup::idB);
		cl.def_readwrite("frameType", &Pythia8::BeamSetup::frameType);
		cl.def_readwrite("boostType", &Pythia8::BeamSetup::boostType);
		cl.def_readwrite("iPDFAsave", &Pythia8::BeamSetup::iPDFAsave);
		cl.def_readwrite("gammaMode", &Pythia8::BeamSetup::gammaMode);
		cl.def_readwrite("mA", &Pythia8::BeamSetup::mA);
		cl.def_readwrite("mB", &Pythia8::BeamSetup::mB);
		cl.def_readwrite("pxA", &Pythia8::BeamSetup::pxA);
		cl.def_readwrite("pxB", &Pythia8::BeamSetup::pxB);
		cl.def_readwrite("pyA", &Pythia8::BeamSetup::pyA);
		cl.def_readwrite("pyB", &Pythia8::BeamSetup::pyB);
		cl.def_readwrite("pzA", &Pythia8::BeamSetup::pzA);
		cl.def_readwrite("pzB", &Pythia8::BeamSetup::pzB);
		cl.def_readwrite("eA", &Pythia8::BeamSetup::eA);
		cl.def_readwrite("eB", &Pythia8::BeamSetup::eB);
		cl.def_readwrite("pzAcm", &Pythia8::BeamSetup::pzAcm);
		cl.def_readwrite("pzBcm", &Pythia8::BeamSetup::pzBcm);
		cl.def_readwrite("eCM", &Pythia8::BeamSetup::eCM);
		cl.def_readwrite("betaZ", &Pythia8::BeamSetup::betaZ);
		cl.def_readwrite("gammaZ", &Pythia8::BeamSetup::gammaZ);
		cl.def_readwrite("pAinit", &Pythia8::BeamSetup::pAinit);
		cl.def_readwrite("pBinit", &Pythia8::BeamSetup::pBinit);
		cl.def_readwrite("pAnow", &Pythia8::BeamSetup::pAnow);
		cl.def_readwrite("pBnow", &Pythia8::BeamSetup::pBnow);
		cl.def_readwrite("MfromCM", &Pythia8::BeamSetup::MfromCM);
		cl.def_readwrite("MtoCM", &Pythia8::BeamSetup::MtoCM);
		cl.def_readwrite("lhaUpPtr", &Pythia8::BeamSetup::lhaUpPtr);
		cl.def_readwrite("beamA", &Pythia8::BeamSetup::beamA);
		cl.def_readwrite("beamB", &Pythia8::BeamSetup::beamB);
		cl.def_readwrite("beamPomA", &Pythia8::BeamSetup::beamPomA);
		cl.def_readwrite("beamPomB", &Pythia8::BeamSetup::beamPomB);
		cl.def_readwrite("beamGamA", &Pythia8::BeamSetup::beamGamA);
		cl.def_readwrite("beamGamB", &Pythia8::BeamSetup::beamGamB);
		cl.def_readwrite("beamVMDA", &Pythia8::BeamSetup::beamVMDA);
		cl.def_readwrite("beamVMDB", &Pythia8::BeamSetup::beamVMDB);
		cl.def_readwrite("idAList", &Pythia8::BeamSetup::idAList);
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1) -> bool { return o.setPDFPtr(a0, a1); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2) -> bool { return o.setPDFPtr(a0, a1, a2); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3) -> bool { return o.setPDFPtr(a0, a1, a2, a3); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6, class std::shared_ptr<class Pythia8::PDF> const & a7) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6, class std::shared_ptr<class Pythia8::PDF> const & a7, class std::shared_ptr<class Pythia8::PDF> const & a8) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"), pybind11::arg("pdfHardGamAPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6, class std::shared_ptr<class Pythia8::PDF> const & a7, class std::shared_ptr<class Pythia8::PDF> const & a8, class std::shared_ptr<class Pythia8::PDF> const & a9) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"), pybind11::arg("pdfHardGamAPtrIn"), pybind11::arg("pdfHardGamBPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6, class std::shared_ptr<class Pythia8::PDF> const & a7, class std::shared_ptr<class Pythia8::PDF> const & a8, class std::shared_ptr<class Pythia8::PDF> const & a9, class std::shared_ptr<class Pythia8::PDF> const & a10) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"), pybind11::arg("pdfHardGamAPtrIn"), pybind11::arg("pdfHardGamBPtrIn"), pybind11::arg("pdfUnresAPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6, class std::shared_ptr<class Pythia8::PDF> const & a7, class std::shared_ptr<class Pythia8::PDF> const & a8, class std::shared_ptr<class Pythia8::PDF> const & a9, class std::shared_ptr<class Pythia8::PDF> const & a10, class std::shared_ptr<class Pythia8::PDF> const & a11) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"), pybind11::arg("pdfHardGamAPtrIn"), pybind11::arg("pdfHardGamBPtrIn"), pybind11::arg("pdfUnresAPtrIn"), pybind11::arg("pdfUnresBPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6, class std::shared_ptr<class Pythia8::PDF> const & a7, class std::shared_ptr<class Pythia8::PDF> const & a8, class std::shared_ptr<class Pythia8::PDF> const & a9, class std::shared_ptr<class Pythia8::PDF> const & a10, class std::shared_ptr<class Pythia8::PDF> const & a11, class std::shared_ptr<class Pythia8::PDF> const & a12) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"), pybind11::arg("pdfHardGamAPtrIn"), pybind11::arg("pdfHardGamBPtrIn"), pybind11::arg("pdfUnresAPtrIn"), pybind11::arg("pdfUnresBPtrIn"), pybind11::arg("pdfUnresGamAPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6, class std::shared_ptr<class Pythia8::PDF> const & a7, class std::shared_ptr<class Pythia8::PDF> const & a8, class std::shared_ptr<class Pythia8::PDF> const & a9, class std::shared_ptr<class Pythia8::PDF> const & a10, class std::shared_ptr<class Pythia8::PDF> const & a11, class std::shared_ptr<class Pythia8::PDF> const & a12, class std::shared_ptr<class Pythia8::PDF> const & a13) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"), pybind11::arg("pdfHardGamAPtrIn"), pybind11::arg("pdfHardGamBPtrIn"), pybind11::arg("pdfUnresAPtrIn"), pybind11::arg("pdfUnresBPtrIn"), pybind11::arg("pdfUnresGamAPtrIn"), pybind11::arg("pdfUnresGamBPtrIn"));
		cl.def("setPDFPtr", [](Pythia8::BeamSetup &o, class std::shared_ptr<class Pythia8::PDF> const & a0, class std::shared_ptr<class Pythia8::PDF> const & a1, class std::shared_ptr<class Pythia8::PDF> const & a2, class std::shared_ptr<class Pythia8::PDF> const & a3, class std::shared_ptr<class Pythia8::PDF> const & a4, class std::shared_ptr<class Pythia8::PDF> const & a5, class std::shared_ptr<class Pythia8::PDF> const & a6, class std::shared_ptr<class Pythia8::PDF> const & a7, class std::shared_ptr<class Pythia8::PDF> const & a8, class std::shared_ptr<class Pythia8::PDF> const & a9, class std::shared_ptr<class Pythia8::PDF> const & a10, class std::shared_ptr<class Pythia8::PDF> const & a11, class std::shared_ptr<class Pythia8::PDF> const & a12, class std::shared_ptr<class Pythia8::PDF> const & a13, class std::shared_ptr<class Pythia8::PDF> const & a14) -> bool { return o.setPDFPtr(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14); }, "", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"), pybind11::arg("pdfHardGamAPtrIn"), pybind11::arg("pdfHardGamBPtrIn"), pybind11::arg("pdfUnresAPtrIn"), pybind11::arg("pdfUnresBPtrIn"), pybind11::arg("pdfUnresGamAPtrIn"), pybind11::arg("pdfUnresGamBPtrIn"), pybind11::arg("pdfVMDAPtrIn"));
		cl.def("setPDFPtr", (bool (Pythia8::BeamSetup::*)(class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>)) &Pythia8::BeamSetup::setPDFPtr, "C++: Pythia8::BeamSetup::setPDFPtr(class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>) --> bool", pybind11::arg("pdfAPtrIn"), pybind11::arg("pdfBPtrIn"), pybind11::arg("pdfHardAPtrIn"), pybind11::arg("pdfHardBPtrIn"), pybind11::arg("pdfPomAPtrIn"), pybind11::arg("pdfPomBPtrIn"), pybind11::arg("pdfGamAPtrIn"), pybind11::arg("pdfGamBPtrIn"), pybind11::arg("pdfHardGamAPtrIn"), pybind11::arg("pdfHardGamBPtrIn"), pybind11::arg("pdfUnresAPtrIn"), pybind11::arg("pdfUnresBPtrIn"), pybind11::arg("pdfUnresGamAPtrIn"), pybind11::arg("pdfUnresGamBPtrIn"), pybind11::arg("pdfVMDAPtrIn"), pybind11::arg("pdfVMDBPtrIn"));
		cl.def("setPDFAPtr", (bool (Pythia8::BeamSetup::*)(class std::shared_ptr<class Pythia8::PDF>)) &Pythia8::BeamSetup::setPDFAPtr, "C++: Pythia8::BeamSetup::setPDFAPtr(class std::shared_ptr<class Pythia8::PDF>) --> bool", pybind11::arg("pdfAPtrIn"));
		cl.def("setPDFBPtr", (bool (Pythia8::BeamSetup::*)(class std::shared_ptr<class Pythia8::PDF>)) &Pythia8::BeamSetup::setPDFBPtr, "C++: Pythia8::BeamSetup::setPDFBPtr(class std::shared_ptr<class Pythia8::PDF>) --> bool", pybind11::arg("pdfBPtrIn"));
		cl.def("setPhotonFluxPtr", (bool (Pythia8::BeamSetup::*)(class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>)) &Pythia8::BeamSetup::setPhotonFluxPtr, "C++: Pythia8::BeamSetup::setPhotonFluxPtr(class std::shared_ptr<class Pythia8::PDF>, class std::shared_ptr<class Pythia8::PDF>) --> bool", pybind11::arg("photonFluxAIn"), pybind11::arg("photonFluxBIn"));
		cl.def("setLHAupPtr", (bool (Pythia8::BeamSetup::*)(class std::shared_ptr<class Pythia8::LHAup>)) &Pythia8::BeamSetup::setLHAupPtr, "C++: Pythia8::BeamSetup::setLHAupPtr(class std::shared_ptr<class Pythia8::LHAup>) --> bool", pybind11::arg("lhaUpPtrIn"));
		cl.def("represent", (int (Pythia8::BeamSetup::*)(int) const) &Pythia8::BeamSetup::represent, "C++: Pythia8::BeamSetup::represent(int) const --> int", pybind11::arg("idIn"));
		cl.def("setBeamIDs", [](Pythia8::BeamSetup &o, int const & a0) -> bool { return o.setBeamIDs(a0); }, "", pybind11::arg("idAin"));
		cl.def("setBeamIDs", (bool (Pythia8::BeamSetup::*)(int, int)) &Pythia8::BeamSetup::setBeamIDs, "C++: Pythia8::BeamSetup::setBeamIDs(int, int) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"));
		cl.def("setKinematics", (bool (Pythia8::BeamSetup::*)(double)) &Pythia8::BeamSetup::setKinematics, "C++: Pythia8::BeamSetup::setKinematics(double) --> bool", pybind11::arg("eCMIn"));
		cl.def("setKinematics", (bool (Pythia8::BeamSetup::*)(double, double)) &Pythia8::BeamSetup::setKinematics, "C++: Pythia8::BeamSetup::setKinematics(double, double) --> bool", pybind11::arg("eAIn"), pybind11::arg("eBIn"));
		cl.def("setKinematics", (bool (Pythia8::BeamSetup::*)(double, double, double, double, double, double)) &Pythia8::BeamSetup::setKinematics, "C++: Pythia8::BeamSetup::setKinematics(double, double, double, double, double, double) --> bool", pybind11::arg("pxAIn"), pybind11::arg("pyAIn"), pybind11::arg("pzAIn"), pybind11::arg("pxBIn"), pybind11::arg("pyBIn"), pybind11::arg("pzBIn"));
		cl.def("setKinematics", (bool (Pythia8::BeamSetup::*)(class Pythia8::Vec4, class Pythia8::Vec4)) &Pythia8::BeamSetup::setKinematics, "C++: Pythia8::BeamSetup::setKinematics(class Pythia8::Vec4, class Pythia8::Vec4) --> bool", pybind11::arg("pAIn"), pybind11::arg("pBIn"));
		cl.def("setBeamShapePtr", (bool (Pythia8::BeamSetup::*)(class std::shared_ptr<class Pythia8::BeamShape>)) &Pythia8::BeamSetup::setBeamShapePtr, "C++: Pythia8::BeamSetup::setBeamShapePtr(class std::shared_ptr<class Pythia8::BeamShape>) --> bool", pybind11::arg("beamShapePtrIn"));
		cl.def("getBeamShapePtr", (class std::shared_ptr<class Pythia8::BeamShape> (Pythia8::BeamSetup::*)()) &Pythia8::BeamSetup::getBeamShapePtr, "C++: Pythia8::BeamSetup::getBeamShapePtr() --> class std::shared_ptr<class Pythia8::BeamShape>");
		cl.def("getPDFPtr", [](Pythia8::BeamSetup &o, int const & a0) -> std::shared_ptr<class Pythia8::PDF> { return o.getPDFPtr(a0); }, "", pybind11::arg("idIn"));
		cl.def("getPDFPtr", [](Pythia8::BeamSetup &o, int const & a0, int const & a1) -> std::shared_ptr<class Pythia8::PDF> { return o.getPDFPtr(a0, a1); }, "", pybind11::arg("idIn"), pybind11::arg("sequence"));
		cl.def("getPDFPtr", [](Pythia8::BeamSetup &o, int const & a0, int const & a1, class std::basic_string<char> const & a2) -> std::shared_ptr<class Pythia8::PDF> { return o.getPDFPtr(a0, a1, a2); }, "", pybind11::arg("idIn"), pybind11::arg("sequence"), pybind11::arg("beam"));
		cl.def("getPDFPtr", (class std::shared_ptr<class Pythia8::PDF> (Pythia8::BeamSetup::*)(int, int, std::string, bool)) &Pythia8::BeamSetup::getPDFPtr, "C++: Pythia8::BeamSetup::getPDFPtr(int, int, std::string, bool) --> class std::shared_ptr<class Pythia8::PDF>", pybind11::arg("idIn"), pybind11::arg("sequence"), pybind11::arg("beam"), pybind11::arg("resolved"));
		cl.def("initFrame", (bool (Pythia8::BeamSetup::*)()) &Pythia8::BeamSetup::initFrame, "C++: Pythia8::BeamSetup::initFrame() --> bool");
		cl.def("initBeams", (bool (Pythia8::BeamSetup::*)(bool, class Pythia8::StringFlav *)) &Pythia8::BeamSetup::initBeams, "C++: Pythia8::BeamSetup::initBeams(bool, class Pythia8::StringFlav *) --> bool", pybind11::arg("doNonPertIn"), pybind11::arg("flavSelPtr"));
		cl.def("getVMDsideA", (bool (Pythia8::BeamSetup::*)()) &Pythia8::BeamSetup::getVMDsideA, "C++: Pythia8::BeamSetup::getVMDsideA() --> bool");
		cl.def("getVMDsideB", (bool (Pythia8::BeamSetup::*)()) &Pythia8::BeamSetup::getVMDsideB, "C++: Pythia8::BeamSetup::getVMDsideB() --> bool");
		cl.def("clear", (void (Pythia8::BeamSetup::*)()) &Pythia8::BeamSetup::clear, "C++: Pythia8::BeamSetup::clear() --> void");
		cl.def("newValenceContent", (void (Pythia8::BeamSetup::*)()) &Pythia8::BeamSetup::newValenceContent, "C++: Pythia8::BeamSetup::newValenceContent() --> void");
		cl.def("nextKinematics", (void (Pythia8::BeamSetup::*)()) &Pythia8::BeamSetup::nextKinematics, "C++: Pythia8::BeamSetup::nextKinematics() --> void");
		cl.def("boostAndVertex", (void (Pythia8::BeamSetup::*)(class Pythia8::Event &, class Pythia8::Event &, bool, bool)) &Pythia8::BeamSetup::boostAndVertex, "C++: Pythia8::BeamSetup::boostAndVertex(class Pythia8::Event &, class Pythia8::Event &, bool, bool) --> void", pybind11::arg("process"), pybind11::arg("event"), pybind11::arg("toLab"), pybind11::arg("setVertex"));
		cl.def("list", (void (Pythia8::BeamSetup::*)() const) &Pythia8::BeamSetup::list, "C++: Pythia8::BeamSetup::list() const --> void");
		cl.def("onInitInfoPtr", (void (Pythia8::BeamSetup::*)()) &Pythia8::BeamSetup::onInitInfoPtr, "C++: Pythia8::BeamSetup::onInitInfoPtr() --> void");
		cl.def("assign", (class Pythia8::BeamSetup & (Pythia8::BeamSetup::*)(const class Pythia8::BeamSetup &)) &Pythia8::BeamSetup::operator=, "C++: Pythia8::BeamSetup::operator=(const class Pythia8::BeamSetup &) --> class Pythia8::BeamSetup &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
