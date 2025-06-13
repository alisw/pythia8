#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/GammaKinematics.h>
#include <Pythia8/HardDiffraction.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonLevel.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/RHadrons.h>
#include <Pythia8/ResonanceDecays.h>
#include <Pythia8/SigmaProcess.h>
#include <Pythia8/SpaceShower.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/TimeShower.h>
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

// Pythia8::HardDiffraction file:Pythia8/HardDiffraction.h line:31
struct PyCallBack_Pythia8_HardDiffraction : public Pythia8::HardDiffraction {
	using Pythia8::HardDiffraction::HardDiffraction;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardDiffraction *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardDiffraction *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardDiffraction *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardDiffraction *>(this), "onStat");
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

// Pythia8::ResonanceDecays file:Pythia8/ResonanceDecays.h line:28
struct PyCallBack_Pythia8_ResonanceDecays : public Pythia8::ResonanceDecays {
	using Pythia8::ResonanceDecays::ResonanceDecays;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecays *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecays *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecays *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecays *>(this), "onStat");
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

// Pythia8::PartonLevel file:Pythia8/PartonLevel.h line:45
struct PyCallBack_Pythia8_PartonLevel : public Pythia8::PartonLevel {
	using Pythia8::PartonLevel::PartonLevel;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonLevel *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartonLevel::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonLevel *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonLevel *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonLevel *>(this), "onStat");
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

// Pythia8::GammaKinematics file:Pythia8/GammaKinematics.h line:23
struct PyCallBack_Pythia8_GammaKinematics : public Pythia8::GammaKinematics {
	using Pythia8::GammaKinematics::GammaKinematics;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GammaKinematics *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GammaKinematics *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GammaKinematics *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GammaKinematics *>(this), "onStat");
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

// Pythia8::PhaseSpace file:Pythia8/PhaseSpace.h line:42
struct PyCallBack_Pythia8_PhaseSpace : public Pythia8::PhaseSpace {
	using Pythia8::PhaseSpace::PhaseSpace;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"PhaseSpace::setupSampling\"");
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"PhaseSpace::trialKin\"");
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"PhaseSpace::finalKin\"");
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "sigmaSumSigned");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PhaseSpace::sigmaSumSigned();
	}
	bool isResolved() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "isResolved");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace::isResolved();
	}
	void rescaleSigma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "rescaleSigma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhaseSpace::rescaleSigma(a0);
	}
	void rescaleMomenta(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "rescaleMomenta");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhaseSpace::rescaleMomenta(a0);
	}
	double weightGammaPDFApprox() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "weightGammaPDFApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PhaseSpace::weightGammaPDFApprox();
	}
	void setGammaKinPtr(class Pythia8::GammaKinematics * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "setGammaKinPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhaseSpace::setGammaKinPtr(a0);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace *>(this), "onStat");
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

void bind_Pythia8_HardDiffraction(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::HardDiffraction file:Pythia8/HardDiffraction.h line:31
		pybind11::class_<Pythia8::HardDiffraction, std::shared_ptr<Pythia8::HardDiffraction>, PyCallBack_Pythia8_HardDiffraction, Pythia8::PhysicsBase> cl(M("Pythia8"), "HardDiffraction", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HardDiffraction(); }, [](){ return new PyCallBack_Pythia8_HardDiffraction(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HardDiffraction const &o){ return new PyCallBack_Pythia8_HardDiffraction(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HardDiffraction const &o){ return new Pythia8::HardDiffraction(o); } ) );
		cl.def("init", (void (Pythia8::HardDiffraction::*)(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *)) &Pythia8::HardDiffraction::init, "C++: Pythia8::HardDiffraction::init(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *) --> void", pybind11::arg("beamAPtrIn"), pybind11::arg("beamBPtrIn"));
		cl.def("isDiffractive", [](Pythia8::HardDiffraction &o) -> bool { return o.isDiffractive(); }, "");
		cl.def("isDiffractive", [](Pythia8::HardDiffraction &o, int const & a0) -> bool { return o.isDiffractive(a0); }, "", pybind11::arg("iBeamIn"));
		cl.def("isDiffractive", [](Pythia8::HardDiffraction &o, int const & a0, int const & a1) -> bool { return o.isDiffractive(a0, a1); }, "", pybind11::arg("iBeamIn"), pybind11::arg("partonIn"));
		cl.def("isDiffractive", [](Pythia8::HardDiffraction &o, int const & a0, int const & a1, double const & a2) -> bool { return o.isDiffractive(a0, a1, a2); }, "", pybind11::arg("iBeamIn"), pybind11::arg("partonIn"), pybind11::arg("xIn"));
		cl.def("isDiffractive", [](Pythia8::HardDiffraction &o, int const & a0, int const & a1, double const & a2, double const & a3) -> bool { return o.isDiffractive(a0, a1, a2, a3); }, "", pybind11::arg("iBeamIn"), pybind11::arg("partonIn"), pybind11::arg("xIn"), pybind11::arg("Q2In"));
		cl.def("isDiffractive", (bool (Pythia8::HardDiffraction::*)(int, int, double, double, double)) &Pythia8::HardDiffraction::isDiffractive, "C++: Pythia8::HardDiffraction::isDiffractive(int, int, double, double, double) --> bool", pybind11::arg("iBeamIn"), pybind11::arg("partonIn"), pybind11::arg("xIn"), pybind11::arg("Q2In"), pybind11::arg("xfIncIn"));
		cl.def("getXPomeronA", (double (Pythia8::HardDiffraction::*)()) &Pythia8::HardDiffraction::getXPomeronA, "C++: Pythia8::HardDiffraction::getXPomeronA() --> double");
		cl.def("getXPomeronB", (double (Pythia8::HardDiffraction::*)()) &Pythia8::HardDiffraction::getXPomeronB, "C++: Pythia8::HardDiffraction::getXPomeronB() --> double");
		cl.def("getTPomeronA", (double (Pythia8::HardDiffraction::*)()) &Pythia8::HardDiffraction::getTPomeronA, "C++: Pythia8::HardDiffraction::getTPomeronA() --> double");
		cl.def("getTPomeronB", (double (Pythia8::HardDiffraction::*)()) &Pythia8::HardDiffraction::getTPomeronB, "C++: Pythia8::HardDiffraction::getTPomeronB() --> double");
		cl.def("getThetaPomeronA", (double (Pythia8::HardDiffraction::*)()) &Pythia8::HardDiffraction::getThetaPomeronA, "C++: Pythia8::HardDiffraction::getThetaPomeronA() --> double");
		cl.def("getThetaPomeronB", (double (Pythia8::HardDiffraction::*)()) &Pythia8::HardDiffraction::getThetaPomeronB, "C++: Pythia8::HardDiffraction::getThetaPomeronB() --> double");
		cl.def("assign", (class Pythia8::HardDiffraction & (Pythia8::HardDiffraction::*)(const class Pythia8::HardDiffraction &)) &Pythia8::HardDiffraction::operator=, "C++: Pythia8::HardDiffraction::operator=(const class Pythia8::HardDiffraction &) --> class Pythia8::HardDiffraction &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceDecays file:Pythia8/ResonanceDecays.h line:28
		pybind11::class_<Pythia8::ResonanceDecays, std::shared_ptr<Pythia8::ResonanceDecays>, PyCallBack_Pythia8_ResonanceDecays, Pythia8::PhysicsBase> cl(M("Pythia8"), "ResonanceDecays", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ResonanceDecays(); }, [](){ return new PyCallBack_Pythia8_ResonanceDecays(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_ResonanceDecays const &o){ return new PyCallBack_Pythia8_ResonanceDecays(o); } ) );
		cl.def( pybind11::init( [](Pythia8::ResonanceDecays const &o){ return new Pythia8::ResonanceDecays(o); } ) );
		cl.def("init", (void (Pythia8::ResonanceDecays::*)()) &Pythia8::ResonanceDecays::init, "C++: Pythia8::ResonanceDecays::init() --> void");
		cl.def("next", [](Pythia8::ResonanceDecays &o, class Pythia8::Event & a0) -> bool { return o.next(a0); }, "", pybind11::arg("process"));
		cl.def("next", (bool (Pythia8::ResonanceDecays::*)(class Pythia8::Event &, int)) &Pythia8::ResonanceDecays::next, "C++: Pythia8::ResonanceDecays::next(class Pythia8::Event &, int) --> bool", pybind11::arg("process"), pybind11::arg("iDecNow"));
		cl.def("assign", (class Pythia8::ResonanceDecays & (Pythia8::ResonanceDecays::*)(const class Pythia8::ResonanceDecays &)) &Pythia8::ResonanceDecays::operator=, "C++: Pythia8::ResonanceDecays::operator=(const class Pythia8::ResonanceDecays &) --> class Pythia8::ResonanceDecays &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PartonLevel file:Pythia8/PartonLevel.h line:45
		pybind11::class_<Pythia8::PartonLevel, std::shared_ptr<Pythia8::PartonLevel>, PyCallBack_Pythia8_PartonLevel, Pythia8::PhysicsBase> cl(M("Pythia8"), "PartonLevel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PartonLevel(); }, [](){ return new PyCallBack_Pythia8_PartonLevel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_PartonLevel const &o){ return new PyCallBack_Pythia8_PartonLevel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::PartonLevel const &o){ return new Pythia8::PartonLevel(o); } ) );
		cl.def_readwrite("timesDecPtr", &Pythia8::PartonLevel::timesDecPtr);
		cl.def_readwrite("timesPtr", &Pythia8::PartonLevel::timesPtr);
		cl.def_readwrite("spacePtr", &Pythia8::PartonLevel::spacePtr);
		cl.def("init", (bool (Pythia8::PartonLevel::*)(class std::shared_ptr<class Pythia8::TimeShower>, class std::shared_ptr<class Pythia8::TimeShower>, class std::shared_ptr<class Pythia8::SpaceShower>, class std::shared_ptr<class Pythia8::RHadrons>, class std::shared_ptr<class Pythia8::MergingHooks>, class std::shared_ptr<class Pythia8::PartonVertex>, class std::shared_ptr<class Pythia8::StringInteractions>, bool)) &Pythia8::PartonLevel::init, "C++: Pythia8::PartonLevel::init(class std::shared_ptr<class Pythia8::TimeShower>, class std::shared_ptr<class Pythia8::TimeShower>, class std::shared_ptr<class Pythia8::SpaceShower>, class std::shared_ptr<class Pythia8::RHadrons>, class std::shared_ptr<class Pythia8::MergingHooks>, class std::shared_ptr<class Pythia8::PartonVertex>, class std::shared_ptr<class Pythia8::StringInteractions>, bool) --> bool", pybind11::arg("timesDecPtrIn"), pybind11::arg("timesPtrIn"), pybind11::arg("spacePtrIn"), pybind11::arg("rHadronsPtrIn"), pybind11::arg("mergingHooksPtr"), pybind11::arg("partonVertexPtrIn"), pybind11::arg("stringInteractionPtrIn"), pybind11::arg("useAsTrial"));
		cl.def("initSwitchID", (void (Pythia8::PartonLevel::*)(const class std::vector<int, class std::allocator<int> > &)) &Pythia8::PartonLevel::initSwitchID, "C++: Pythia8::PartonLevel::initSwitchID(const class std::vector<int, class std::allocator<int> > &) --> void", pybind11::arg("idAList"));
		cl.def("setBeamID", [](Pythia8::PartonLevel &o) -> void { return o.setBeamID(); }, "");
		cl.def("setBeamID", (void (Pythia8::PartonLevel::*)(int)) &Pythia8::PartonLevel::setBeamID, "C++: Pythia8::PartonLevel::setBeamID(int) --> void", pybind11::arg("iPDFA"));
		cl.def("next", (bool (Pythia8::PartonLevel::*)(class Pythia8::Event &, class Pythia8::Event &)) &Pythia8::PartonLevel::next, "C++: Pythia8::PartonLevel::next(class Pythia8::Event &, class Pythia8::Event &) --> bool", pybind11::arg("process"), pybind11::arg("event"));
		cl.def("setupShowerSys", (void (Pythia8::PartonLevel::*)(class Pythia8::Event &, class Pythia8::Event &)) &Pythia8::PartonLevel::setupShowerSys, "C++: Pythia8::PartonLevel::setupShowerSys(class Pythia8::Event &, class Pythia8::Event &) --> void", pybind11::arg("process"), pybind11::arg("event"));
		cl.def("resonanceShowers", (bool (Pythia8::PartonLevel::*)(class Pythia8::Event &, class Pythia8::Event &, bool)) &Pythia8::PartonLevel::resonanceShowers, "C++: Pythia8::PartonLevel::resonanceShowers(class Pythia8::Event &, class Pythia8::Event &, bool) --> bool", pybind11::arg("process"), pybind11::arg("event"), pybind11::arg("skipForR"));
		cl.def("wzDecayShowers", (bool (Pythia8::PartonLevel::*)(class Pythia8::Event &)) &Pythia8::PartonLevel::wzDecayShowers, "C++: Pythia8::PartonLevel::wzDecayShowers(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("hasVetoed", (bool (Pythia8::PartonLevel::*)() const) &Pythia8::PartonLevel::hasVetoed, "C++: Pythia8::PartonLevel::hasVetoed() const --> bool");
		cl.def("hasVetoedDiff", (bool (Pythia8::PartonLevel::*)() const) &Pythia8::PartonLevel::hasVetoedDiff, "C++: Pythia8::PartonLevel::hasVetoedDiff() const --> bool");
		cl.def("hasVetoedMerging", (bool (Pythia8::PartonLevel::*)() const) &Pythia8::PartonLevel::hasVetoedMerging, "C++: Pythia8::PartonLevel::hasVetoedMerging() const --> bool");
		cl.def("accumulate", (void (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::accumulate, "C++: Pythia8::PartonLevel::accumulate() --> void");
		cl.def("statistics", [](Pythia8::PartonLevel &o) -> void { return o.statistics(); }, "");
		cl.def("statistics", (void (Pythia8::PartonLevel::*)(bool)) &Pythia8::PartonLevel::statistics, "C++: Pythia8::PartonLevel::statistics(bool) --> void", pybind11::arg("reset"));
		cl.def("resetStatistics", (void (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::resetStatistics, "C++: Pythia8::PartonLevel::resetStatistics() --> void");
		cl.def("resetTrial", (void (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::resetTrial, "C++: Pythia8::PartonLevel::resetTrial() --> void");
		cl.def("pTLastInShower", (double (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::pTLastInShower, "C++: Pythia8::PartonLevel::pTLastInShower() --> double");
		cl.def("typeLastInShower", (int (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::typeLastInShower, "C++: Pythia8::PartonLevel::typeLastInShower() --> int");
		cl.def("canEnhanceTrial", (bool (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::canEnhanceTrial, "C++: Pythia8::PartonLevel::canEnhanceTrial() --> bool");
		cl.def("getEnhancedTrialPT", (double (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::getEnhancedTrialPT, "C++: Pythia8::PartonLevel::getEnhancedTrialPT() --> double");
		cl.def("getEnhancedTrialWeight", (double (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::getEnhancedTrialWeight, "C++: Pythia8::PartonLevel::getEnhancedTrialWeight() --> double");
		cl.def("onInitInfoPtr", (void (Pythia8::PartonLevel::*)()) &Pythia8::PartonLevel::onInitInfoPtr, "C++: Pythia8::PartonLevel::onInitInfoPtr() --> void");
		cl.def("assign", (class Pythia8::PartonLevel & (Pythia8::PartonLevel::*)(const class Pythia8::PartonLevel &)) &Pythia8::PartonLevel::operator=, "C++: Pythia8::PartonLevel::operator=(const class Pythia8::PartonLevel &) --> class Pythia8::PartonLevel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::GammaKinematics file:Pythia8/GammaKinematics.h line:23
		pybind11::class_<Pythia8::GammaKinematics, std::shared_ptr<Pythia8::GammaKinematics>, PyCallBack_Pythia8_GammaKinematics, Pythia8::PhysicsBase> cl(M("Pythia8"), "GammaKinematics", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::GammaKinematics(); }, [](){ return new PyCallBack_Pythia8_GammaKinematics(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_GammaKinematics const &o){ return new PyCallBack_Pythia8_GammaKinematics(o); } ) );
		cl.def( pybind11::init( [](Pythia8::GammaKinematics const &o){ return new Pythia8::GammaKinematics(o); } ) );
		cl.def("init", (bool (Pythia8::GammaKinematics::*)()) &Pythia8::GammaKinematics::init, "C++: Pythia8::GammaKinematics::init() --> bool");
		cl.def("sampleKTgamma", [](Pythia8::GammaKinematics &o) -> bool { return o.sampleKTgamma(); }, "");
		cl.def("sampleKTgamma", (bool (Pythia8::GammaKinematics::*)(bool)) &Pythia8::GammaKinematics::sampleKTgamma, "C++: Pythia8::GammaKinematics::sampleKTgamma(bool) --> bool", pybind11::arg("nonDiff"));
		cl.def("deriveKin", (bool (Pythia8::GammaKinematics::*)(double, double, double, double)) &Pythia8::GammaKinematics::deriveKin, "C++: Pythia8::GammaKinematics::deriveKin(double, double, double, double) --> bool", pybind11::arg("xGamma"), pybind11::arg("Q2gamma"), pybind11::arg("m2beam"), pybind11::arg("eCM2"));
		cl.def("finalize", (bool (Pythia8::GammaKinematics::*)()) &Pythia8::GammaKinematics::finalize, "C++: Pythia8::GammaKinematics::finalize() --> bool");
		cl.def("trialKinSoftPhaseSpaceSampling", (bool (Pythia8::GammaKinematics::*)()) &Pythia8::GammaKinematics::trialKinSoftPhaseSpaceSampling, "C++: Pythia8::GammaKinematics::trialKinSoftPhaseSpaceSampling() --> bool");
		cl.def("fluxWeight", (double (Pythia8::GammaKinematics::*)()) &Pythia8::GammaKinematics::fluxWeight, "C++: Pythia8::GammaKinematics::fluxWeight() --> double");
		cl.def("setupSoftPhaseSpaceSampling", (double (Pythia8::GammaKinematics::*)(double)) &Pythia8::GammaKinematics::setupSoftPhaseSpaceSampling, "C++: Pythia8::GammaKinematics::setupSoftPhaseSpaceSampling(double) --> double", pybind11::arg("sigmaMax"));
		cl.def("calcNewSHat", (double (Pythia8::GammaKinematics::*)(double)) &Pythia8::GammaKinematics::calcNewSHat, "C++: Pythia8::GammaKinematics::calcNewSHat(double) --> double", pybind11::arg("sHatOld"));
		cl.def("getQ2gamma1", (double (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::getQ2gamma1, "C++: Pythia8::GammaKinematics::getQ2gamma1() const --> double");
		cl.def("getQ2gamma2", (double (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::getQ2gamma2, "C++: Pythia8::GammaKinematics::getQ2gamma2() const --> double");
		cl.def("getPhi1", (double (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::getPhi1, "C++: Pythia8::GammaKinematics::getPhi1() const --> double");
		cl.def("getPhi2", (double (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::getPhi2, "C++: Pythia8::GammaKinematics::getPhi2() const --> double");
		cl.def("getKT1", (double (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::getKT1, "C++: Pythia8::GammaKinematics::getKT1() const --> double");
		cl.def("getKT2", (double (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::getKT2, "C++: Pythia8::GammaKinematics::getKT2() const --> double");
		cl.def("eCMsub", (double (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::eCMsub, "C++: Pythia8::GammaKinematics::eCMsub() const --> double");
		cl.def("weight", (double (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::weight, "C++: Pythia8::GammaKinematics::weight() const --> double");
		cl.def("idInA", (int (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::idInA, "C++: Pythia8::GammaKinematics::idInA() const --> int");
		cl.def("idInB", (int (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::idInB, "C++: Pythia8::GammaKinematics::idInB() const --> int");
		cl.def("hasNewSHat", (bool (Pythia8::GammaKinematics::*)() const) &Pythia8::GammaKinematics::hasNewSHat, "C++: Pythia8::GammaKinematics::hasNewSHat() const --> bool");
		cl.def("assign", (class Pythia8::GammaKinematics & (Pythia8::GammaKinematics::*)(const class Pythia8::GammaKinematics &)) &Pythia8::GammaKinematics::operator=, "C++: Pythia8::GammaKinematics::operator=(const class Pythia8::GammaKinematics &) --> class Pythia8::GammaKinematics &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpace file:Pythia8/PhaseSpace.h line:42
		pybind11::class_<Pythia8::PhaseSpace, std::shared_ptr<Pythia8::PhaseSpace>, PyCallBack_Pythia8_PhaseSpace, Pythia8::PhysicsBase> cl(M("Pythia8"), "PhaseSpace", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_PhaseSpace(); } ) );
		cl.def(pybind11::init<PyCallBack_Pythia8_PhaseSpace const &>());
		cl.def_readwrite("sigmaProcessPtr", &Pythia8::PhaseSpace::sigmaProcessPtr);
		cl.def_readwrite("lhaUpPtr", &Pythia8::PhaseSpace::lhaUpPtr);
		cl.def_readwrite("useBreitWigners", &Pythia8::PhaseSpace::useBreitWigners);
		cl.def_readwrite("doEnergySpread", &Pythia8::PhaseSpace::doEnergySpread);
		cl.def_readwrite("showSearch", &Pythia8::PhaseSpace::showSearch);
		cl.def_readwrite("showViolation", &Pythia8::PhaseSpace::showViolation);
		cl.def_readwrite("increaseMaximum", &Pythia8::PhaseSpace::increaseMaximum);
		cl.def_readwrite("hasQ2Min", &Pythia8::PhaseSpace::hasQ2Min);
		cl.def_readwrite("gmZmodeGlobal", &Pythia8::PhaseSpace::gmZmodeGlobal);
		cl.def_readwrite("mHatGlobalMin", &Pythia8::PhaseSpace::mHatGlobalMin);
		cl.def_readwrite("mHatGlobalMax", &Pythia8::PhaseSpace::mHatGlobalMax);
		cl.def_readwrite("pTHatGlobalMin", &Pythia8::PhaseSpace::pTHatGlobalMin);
		cl.def_readwrite("pTHatGlobalMax", &Pythia8::PhaseSpace::pTHatGlobalMax);
		cl.def_readwrite("Q2GlobalMin", &Pythia8::PhaseSpace::Q2GlobalMin);
		cl.def_readwrite("pTHatMinDiverge", &Pythia8::PhaseSpace::pTHatMinDiverge);
		cl.def_readwrite("minWidthBreitWigners", &Pythia8::PhaseSpace::minWidthBreitWigners);
		cl.def_readwrite("minWidthNarrowBW", &Pythia8::PhaseSpace::minWidthNarrowBW);
		cl.def_readwrite("idA", &Pythia8::PhaseSpace::idA);
		cl.def_readwrite("idB", &Pythia8::PhaseSpace::idB);
		cl.def_readwrite("idAold", &Pythia8::PhaseSpace::idAold);
		cl.def_readwrite("idBold", &Pythia8::PhaseSpace::idBold);
		cl.def_readwrite("idAgm", &Pythia8::PhaseSpace::idAgm);
		cl.def_readwrite("idBgm", &Pythia8::PhaseSpace::idBgm);
		cl.def_readwrite("mA", &Pythia8::PhaseSpace::mA);
		cl.def_readwrite("mB", &Pythia8::PhaseSpace::mB);
		cl.def_readwrite("eCM", &Pythia8::PhaseSpace::eCM);
		cl.def_readwrite("s", &Pythia8::PhaseSpace::s);
		cl.def_readwrite("sigmaMxGm", &Pythia8::PhaseSpace::sigmaMxGm);
		cl.def_readwrite("hasLeptonBeamA", &Pythia8::PhaseSpace::hasLeptonBeamA);
		cl.def_readwrite("hasLeptonBeamB", &Pythia8::PhaseSpace::hasLeptonBeamB);
		cl.def_readwrite("hasOneLeptonBeam", &Pythia8::PhaseSpace::hasOneLeptonBeam);
		cl.def_readwrite("hasTwoLeptonBeams", &Pythia8::PhaseSpace::hasTwoLeptonBeams);
		cl.def_readwrite("hasPointGammaA", &Pythia8::PhaseSpace::hasPointGammaA);
		cl.def_readwrite("hasPointGammaB", &Pythia8::PhaseSpace::hasPointGammaB);
		cl.def_readwrite("hasOnePointParticle", &Pythia8::PhaseSpace::hasOnePointParticle);
		cl.def_readwrite("hasTwoPointParticles", &Pythia8::PhaseSpace::hasTwoPointParticles);
		cl.def_readwrite("hasGamma", &Pythia8::PhaseSpace::hasGamma);
		cl.def_readwrite("hasVMD", &Pythia8::PhaseSpace::hasVMD);
		cl.def_readwrite("newSigmaMx", &Pythia8::PhaseSpace::newSigmaMx);
		cl.def_readwrite("canModifySigma", &Pythia8::PhaseSpace::canModifySigma);
		cl.def_readwrite("canBiasSelection", &Pythia8::PhaseSpace::canBiasSelection);
		cl.def_readwrite("canBias2Sel", &Pythia8::PhaseSpace::canBias2Sel);
		cl.def_readwrite("gmZmode", &Pythia8::PhaseSpace::gmZmode);
		cl.def_readwrite("bias2SelPow", &Pythia8::PhaseSpace::bias2SelPow);
		cl.def_readwrite("bias2SelRef", &Pythia8::PhaseSpace::bias2SelRef);
		cl.def_readwrite("wtBW", &Pythia8::PhaseSpace::wtBW);
		cl.def_readwrite("sigmaNw", &Pythia8::PhaseSpace::sigmaNw);
		cl.def_readwrite("sigmaMx", &Pythia8::PhaseSpace::sigmaMx);
		cl.def_readwrite("sigmaPos", &Pythia8::PhaseSpace::sigmaPos);
		cl.def_readwrite("sigmaNeg", &Pythia8::PhaseSpace::sigmaNeg);
		cl.def_readwrite("biasWt", &Pythia8::PhaseSpace::biasWt);
		cl.def_readwrite("mHatMin", &Pythia8::PhaseSpace::mHatMin);
		cl.def_readwrite("mHatMax", &Pythia8::PhaseSpace::mHatMax);
		cl.def_readwrite("sHatMin", &Pythia8::PhaseSpace::sHatMin);
		cl.def_readwrite("sHatMax", &Pythia8::PhaseSpace::sHatMax);
		cl.def_readwrite("pTHatMin", &Pythia8::PhaseSpace::pTHatMin);
		cl.def_readwrite("pTHatMax", &Pythia8::PhaseSpace::pTHatMax);
		cl.def_readwrite("pT2HatMin", &Pythia8::PhaseSpace::pT2HatMin);
		cl.def_readwrite("pT2HatMax", &Pythia8::PhaseSpace::pT2HatMax);
		cl.def_readwrite("x1H", &Pythia8::PhaseSpace::x1H);
		cl.def_readwrite("x2H", &Pythia8::PhaseSpace::x2H);
		cl.def_readwrite("m3", &Pythia8::PhaseSpace::m3);
		cl.def_readwrite("m4", &Pythia8::PhaseSpace::m4);
		cl.def_readwrite("m5", &Pythia8::PhaseSpace::m5);
		cl.def_readwrite("s3", &Pythia8::PhaseSpace::s3);
		cl.def_readwrite("s4", &Pythia8::PhaseSpace::s4);
		cl.def_readwrite("s5", &Pythia8::PhaseSpace::s5);
		cl.def_readwrite("mHat", &Pythia8::PhaseSpace::mHat);
		cl.def_readwrite("sH", &Pythia8::PhaseSpace::sH);
		cl.def_readwrite("tH", &Pythia8::PhaseSpace::tH);
		cl.def_readwrite("uH", &Pythia8::PhaseSpace::uH);
		cl.def_readwrite("pAbs", &Pythia8::PhaseSpace::pAbs);
		cl.def_readwrite("p2Abs", &Pythia8::PhaseSpace::p2Abs);
		cl.def_readwrite("pTH", &Pythia8::PhaseSpace::pTH);
		cl.def_readwrite("theta", &Pythia8::PhaseSpace::theta);
		cl.def_readwrite("phi", &Pythia8::PhaseSpace::phi);
		cl.def_readwrite("betaZ", &Pythia8::PhaseSpace::betaZ);
		cl.def_readwrite("idResA", &Pythia8::PhaseSpace::idResA);
		cl.def_readwrite("idResB", &Pythia8::PhaseSpace::idResB);
		cl.def_readwrite("mResA", &Pythia8::PhaseSpace::mResA);
		cl.def_readwrite("mResB", &Pythia8::PhaseSpace::mResB);
		cl.def_readwrite("GammaResA", &Pythia8::PhaseSpace::GammaResA);
		cl.def_readwrite("GammaResB", &Pythia8::PhaseSpace::GammaResB);
		cl.def_readwrite("tauResA", &Pythia8::PhaseSpace::tauResA);
		cl.def_readwrite("tauResB", &Pythia8::PhaseSpace::tauResB);
		cl.def_readwrite("widResA", &Pythia8::PhaseSpace::widResA);
		cl.def_readwrite("widResB", &Pythia8::PhaseSpace::widResB);
		cl.def_readwrite("sameResMass", &Pythia8::PhaseSpace::sameResMass);
		cl.def_readwrite("useMirrorWeight", &Pythia8::PhaseSpace::useMirrorWeight);
		cl.def_readwrite("hasNegZ", &Pythia8::PhaseSpace::hasNegZ);
		cl.def_readwrite("hasPosZ", &Pythia8::PhaseSpace::hasPosZ);
		cl.def_readwrite("tau", &Pythia8::PhaseSpace::tau);
		cl.def_readwrite("y", &Pythia8::PhaseSpace::y);
		cl.def_readwrite("z", &Pythia8::PhaseSpace::z);
		cl.def_readwrite("tauMin", &Pythia8::PhaseSpace::tauMin);
		cl.def_readwrite("tauMax", &Pythia8::PhaseSpace::tauMax);
		cl.def_readwrite("yMax", &Pythia8::PhaseSpace::yMax);
		cl.def_readwrite("zMin", &Pythia8::PhaseSpace::zMin);
		cl.def_readwrite("zMax", &Pythia8::PhaseSpace::zMax);
		cl.def_readwrite("ratio34", &Pythia8::PhaseSpace::ratio34);
		cl.def_readwrite("unity34", &Pythia8::PhaseSpace::unity34);
		cl.def_readwrite("zNeg", &Pythia8::PhaseSpace::zNeg);
		cl.def_readwrite("zPos", &Pythia8::PhaseSpace::zPos);
		cl.def_readwrite("wtTau", &Pythia8::PhaseSpace::wtTau);
		cl.def_readwrite("wtY", &Pythia8::PhaseSpace::wtY);
		cl.def_readwrite("wtZ", &Pythia8::PhaseSpace::wtZ);
		cl.def_readwrite("wt3Body", &Pythia8::PhaseSpace::wt3Body);
		cl.def_readwrite("runBW3H", &Pythia8::PhaseSpace::runBW3H);
		cl.def_readwrite("runBW4H", &Pythia8::PhaseSpace::runBW4H);
		cl.def_readwrite("runBW5H", &Pythia8::PhaseSpace::runBW5H);
		cl.def_readwrite("intTau0", &Pythia8::PhaseSpace::intTau0);
		cl.def_readwrite("intTau1", &Pythia8::PhaseSpace::intTau1);
		cl.def_readwrite("intTau2", &Pythia8::PhaseSpace::intTau2);
		cl.def_readwrite("intTau3", &Pythia8::PhaseSpace::intTau3);
		cl.def_readwrite("intTau4", &Pythia8::PhaseSpace::intTau4);
		cl.def_readwrite("intTau5", &Pythia8::PhaseSpace::intTau5);
		cl.def_readwrite("intTau6", &Pythia8::PhaseSpace::intTau6);
		cl.def_readwrite("intY0", &Pythia8::PhaseSpace::intY0);
		cl.def_readwrite("intY12", &Pythia8::PhaseSpace::intY12);
		cl.def_readwrite("intY34", &Pythia8::PhaseSpace::intY34);
		cl.def_readwrite("intY56", &Pythia8::PhaseSpace::intY56);
		cl.def_readwrite("mTchan1", &Pythia8::PhaseSpace::mTchan1);
		cl.def_readwrite("sTchan1", &Pythia8::PhaseSpace::sTchan1);
		cl.def_readwrite("mTchan2", &Pythia8::PhaseSpace::mTchan2);
		cl.def_readwrite("sTchan2", &Pythia8::PhaseSpace::sTchan2);
		cl.def_readwrite("frac3Flat", &Pythia8::PhaseSpace::frac3Flat);
		cl.def_readwrite("frac3Pow1", &Pythia8::PhaseSpace::frac3Pow1);
		cl.def_readwrite("frac3Pow2", &Pythia8::PhaseSpace::frac3Pow2);
		cl.def_readwrite("zNegMin", &Pythia8::PhaseSpace::zNegMin);
		cl.def_readwrite("zNegMax", &Pythia8::PhaseSpace::zNegMax);
		cl.def_readwrite("zPosMin", &Pythia8::PhaseSpace::zPosMin);
		cl.def_readwrite("zPosMax", &Pythia8::PhaseSpace::zPosMax);
		cl.def_readwrite("p3cm", &Pythia8::PhaseSpace::p3cm);
		cl.def_readwrite("p4cm", &Pythia8::PhaseSpace::p4cm);
		cl.def_readwrite("p5cm", &Pythia8::PhaseSpace::p5cm);
		cl.def_readwrite("nTau", &Pythia8::PhaseSpace::nTau);
		cl.def_readwrite("nY", &Pythia8::PhaseSpace::nY);
		cl.def_readwrite("nZ", &Pythia8::PhaseSpace::nZ);
		cl.def("init", (void (Pythia8::PhaseSpace::*)(bool, class std::shared_ptr<class Pythia8::SigmaProcess>)) &Pythia8::PhaseSpace::init, "C++: Pythia8::PhaseSpace::init(bool, class std::shared_ptr<class Pythia8::SigmaProcess>) --> void", pybind11::arg("isFirst"), pybind11::arg("sigmaProcessPtrIn"));
		cl.def("updateBeamIDs", (void (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::updateBeamIDs, "C++: Pythia8::PhaseSpace::updateBeamIDs() --> void");
		cl.def("newECM", (void (Pythia8::PhaseSpace::*)(double)) &Pythia8::PhaseSpace::newECM, "C++: Pythia8::PhaseSpace::newECM(double) --> void", pybind11::arg("eCMin"));
		cl.def("setLHAPtr", (void (Pythia8::PhaseSpace::*)(class std::shared_ptr<class Pythia8::LHAup>)) &Pythia8::PhaseSpace::setLHAPtr, "C++: Pythia8::PhaseSpace::setLHAPtr(class std::shared_ptr<class Pythia8::LHAup>) --> void", pybind11::arg("lhaUpPtrIn"));
		cl.def("setupSampling", (bool (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::setupSampling, "C++: Pythia8::PhaseSpace::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace &o) -> bool { return o.trialKin(); }, "");
		cl.def("trialKin", [](Pythia8::PhaseSpace &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg("inEvent"));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace::*)(bool, bool)) &Pythia8::PhaseSpace::trialKin, "C++: Pythia8::PhaseSpace::trialKin(bool, bool) --> bool", pybind11::arg("inEvent"), pybind11::arg("repeatSame"));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::finalKin, "C++: Pythia8::PhaseSpace::finalKin() --> bool");
		cl.def("decayKinematics", (void (Pythia8::PhaseSpace::*)(class Pythia8::Event &)) &Pythia8::PhaseSpace::decayKinematics, "C++: Pythia8::PhaseSpace::decayKinematics(class Pythia8::Event &) --> void", pybind11::arg("process"));
		cl.def("sigmaNow", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::sigmaNow, "C++: Pythia8::PhaseSpace::sigmaNow() const --> double");
		cl.def("sigmaMax", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::sigmaMax, "C++: Pythia8::PhaseSpace::sigmaMax() const --> double");
		cl.def("biasSelectionWeight", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::biasSelectionWeight, "C++: Pythia8::PhaseSpace::biasSelectionWeight() const --> double");
		cl.def("newSigmaMax", (bool (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::newSigmaMax, "C++: Pythia8::PhaseSpace::newSigmaMax() const --> bool");
		cl.def("setSigmaMax", (void (Pythia8::PhaseSpace::*)(double)) &Pythia8::PhaseSpace::setSigmaMax, "C++: Pythia8::PhaseSpace::setSigmaMax(double) --> void", pybind11::arg("sigmaMaxIn"));
		cl.def("sigmaMaxSwitch", (double (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::sigmaMaxSwitch, "C++: Pythia8::PhaseSpace::sigmaMaxSwitch() --> double");
		cl.def("sigmaSumSigned", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::sigmaSumSigned, "C++: Pythia8::PhaseSpace::sigmaSumSigned() const --> double");
		cl.def("p", (class Pythia8::Vec4 (Pythia8::PhaseSpace::*)(int) const) &Pythia8::PhaseSpace::p, "C++: Pythia8::PhaseSpace::p(int) const --> class Pythia8::Vec4", pybind11::arg("i"));
		cl.def("m", (double (Pythia8::PhaseSpace::*)(int) const) &Pythia8::PhaseSpace::m, "C++: Pythia8::PhaseSpace::m(int) const --> double", pybind11::arg("i"));
		cl.def("setP", (void (Pythia8::PhaseSpace::*)(int, class Pythia8::Vec4)) &Pythia8::PhaseSpace::setP, "C++: Pythia8::PhaseSpace::setP(int, class Pythia8::Vec4) --> void", pybind11::arg("i"), pybind11::arg("pNew"));
		cl.def("ecm", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::ecm, "C++: Pythia8::PhaseSpace::ecm() const --> double");
		cl.def("x1", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::x1, "C++: Pythia8::PhaseSpace::x1() const --> double");
		cl.def("x2", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::x2, "C++: Pythia8::PhaseSpace::x2() const --> double");
		cl.def("sHat", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::sHat, "C++: Pythia8::PhaseSpace::sHat() const --> double");
		cl.def("tHat", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::tHat, "C++: Pythia8::PhaseSpace::tHat() const --> double");
		cl.def("uHat", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::uHat, "C++: Pythia8::PhaseSpace::uHat() const --> double");
		cl.def("pTHat", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::pTHat, "C++: Pythia8::PhaseSpace::pTHat() const --> double");
		cl.def("thetaHat", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::thetaHat, "C++: Pythia8::PhaseSpace::thetaHat() const --> double");
		cl.def("phiHat", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::phiHat, "C++: Pythia8::PhaseSpace::phiHat() const --> double");
		cl.def("runBW3", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::runBW3, "C++: Pythia8::PhaseSpace::runBW3() const --> double");
		cl.def("runBW4", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::runBW4, "C++: Pythia8::PhaseSpace::runBW4() const --> double");
		cl.def("runBW5", (double (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::runBW5, "C++: Pythia8::PhaseSpace::runBW5() const --> double");
		cl.def("isResolved", (bool (Pythia8::PhaseSpace::*)() const) &Pythia8::PhaseSpace::isResolved, "C++: Pythia8::PhaseSpace::isResolved() const --> bool");
		cl.def("rescaleSigma", (void (Pythia8::PhaseSpace::*)(double)) &Pythia8::PhaseSpace::rescaleSigma, "C++: Pythia8::PhaseSpace::rescaleSigma(double) --> void", pybind11::arg(""));
		cl.def("rescaleMomenta", (void (Pythia8::PhaseSpace::*)(double)) &Pythia8::PhaseSpace::rescaleMomenta, "C++: Pythia8::PhaseSpace::rescaleMomenta(double) --> void", pybind11::arg(""));
		cl.def("weightGammaPDFApprox", (double (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::weightGammaPDFApprox, "C++: Pythia8::PhaseSpace::weightGammaPDFApprox() --> double");
		cl.def("setGammaKinPtr", (void (Pythia8::PhaseSpace::*)(class Pythia8::GammaKinematics *)) &Pythia8::PhaseSpace::setGammaKinPtr, "C++: Pythia8::PhaseSpace::setGammaKinPtr(class Pythia8::GammaKinematics *) --> void", pybind11::arg("gammaKinPtrIn"));
		cl.def("decayKinematicsStep", (void (Pythia8::PhaseSpace::*)(class Pythia8::Event &, int)) &Pythia8::PhaseSpace::decayKinematicsStep, "C++: Pythia8::PhaseSpace::decayKinematicsStep(class Pythia8::Event &, int) --> void", pybind11::arg("process"), pybind11::arg("iRes"));
		cl.def("setup3Body", (void (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::setup3Body, "C++: Pythia8::PhaseSpace::setup3Body() --> void");
		cl.def("setupSampling123", (bool (Pythia8::PhaseSpace::*)(bool, bool)) &Pythia8::PhaseSpace::setupSampling123, "C++: Pythia8::PhaseSpace::setupSampling123(bool, bool) --> bool", pybind11::arg("is2"), pybind11::arg("is3"));
		cl.def("trialKin123", [](Pythia8::PhaseSpace &o, bool const & a0, bool const & a1) -> bool { return o.trialKin123(a0, a1); }, "", pybind11::arg("is2"), pybind11::arg("is3"));
		cl.def("trialKin123", (bool (Pythia8::PhaseSpace::*)(bool, bool, bool)) &Pythia8::PhaseSpace::trialKin123, "C++: Pythia8::PhaseSpace::trialKin123(bool, bool, bool) --> bool", pybind11::arg("is2"), pybind11::arg("is3"), pybind11::arg("inEvent"));
		cl.def("limitTau", (bool (Pythia8::PhaseSpace::*)(bool, bool)) &Pythia8::PhaseSpace::limitTau, "C++: Pythia8::PhaseSpace::limitTau(bool, bool) --> bool", pybind11::arg("is2"), pybind11::arg("is3"));
		cl.def("limitY", (bool (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::limitY, "C++: Pythia8::PhaseSpace::limitY() --> bool");
		cl.def("limitZ", (bool (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::limitZ, "C++: Pythia8::PhaseSpace::limitZ() --> bool");
		cl.def("selectTau", (void (Pythia8::PhaseSpace::*)(int, double, bool)) &Pythia8::PhaseSpace::selectTau, "C++: Pythia8::PhaseSpace::selectTau(int, double, bool) --> void", pybind11::arg("iTau"), pybind11::arg("tauVal"), pybind11::arg("is2"));
		cl.def("selectY", (void (Pythia8::PhaseSpace::*)(int, double)) &Pythia8::PhaseSpace::selectY, "C++: Pythia8::PhaseSpace::selectY(int, double) --> void", pybind11::arg("iY"), pybind11::arg("yVal"));
		cl.def("selectZ", (void (Pythia8::PhaseSpace::*)(int, double)) &Pythia8::PhaseSpace::selectZ, "C++: Pythia8::PhaseSpace::selectZ(int, double) --> void", pybind11::arg("iZ"), pybind11::arg("zVal"));
		cl.def("select3Body", (bool (Pythia8::PhaseSpace::*)()) &Pythia8::PhaseSpace::select3Body, "C++: Pythia8::PhaseSpace::select3Body() --> bool");
		cl.def("setupMass1", (void (Pythia8::PhaseSpace::*)(int)) &Pythia8::PhaseSpace::setupMass1, "C++: Pythia8::PhaseSpace::setupMass1(int) --> void", pybind11::arg("iM"));
		cl.def("setupMass2", (void (Pythia8::PhaseSpace::*)(int, double)) &Pythia8::PhaseSpace::setupMass2, "C++: Pythia8::PhaseSpace::setupMass2(int, double) --> void", pybind11::arg("iM"), pybind11::arg("distToThresh"));
		cl.def("trialMass", (void (Pythia8::PhaseSpace::*)(int)) &Pythia8::PhaseSpace::trialMass, "C++: Pythia8::PhaseSpace::trialMass(int) --> void", pybind11::arg("iM"));
		cl.def("weightMass", (double (Pythia8::PhaseSpace::*)(int)) &Pythia8::PhaseSpace::weightMass, "C++: Pythia8::PhaseSpace::weightMass(int) --> double", pybind11::arg("iM"));
		cl.def("tRange", (struct std::pair<double, double> (Pythia8::PhaseSpace::*)(double, double, double, double, double)) &Pythia8::PhaseSpace::tRange, "C++: Pythia8::PhaseSpace::tRange(double, double, double, double, double) --> struct std::pair<double, double>", pybind11::arg("sIn"), pybind11::arg("s1In"), pybind11::arg("s2In"), pybind11::arg("s3In"), pybind11::arg("s4In"));
		cl.def("tInRange", (bool (Pythia8::PhaseSpace::*)(double, double, double, double, double, double)) &Pythia8::PhaseSpace::tInRange, "C++: Pythia8::PhaseSpace::tInRange(double, double, double, double, double, double) --> bool", pybind11::arg("tIn"), pybind11::arg("sIn"), pybind11::arg("s1In"), pybind11::arg("s2In"), pybind11::arg("s3In"), pybind11::arg("s4In"));
		cl.def("assign", (class Pythia8::PhaseSpace & (Pythia8::PhaseSpace::*)(const class Pythia8::PhaseSpace &)) &Pythia8::PhaseSpace::operator=, "C++: Pythia8::PhaseSpace::operator=(const class Pythia8::PhaseSpace &) --> class Pythia8::PhaseSpace &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
