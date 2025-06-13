#include <Pythia8/GammaKinematics.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/PhysicsBase.h>
#include <sstream> // __str__

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

// Pythia8::PhaseSpace2to1tauy file:Pythia8/PhaseSpace.h line:299
struct PyCallBack_Pythia8_PhaseSpace2to1tauy : public Pythia8::PhaseSpace2to1tauy {
	using Pythia8::PhaseSpace2to1tauy::PhaseSpace2to1tauy;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to1tauy::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to1tauy::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to1tauy::finalKin();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "sigmaSumSigned");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "isResolved");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "rescaleSigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "rescaleMomenta");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "weightGammaPDFApprox");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to1tauy *>(this), "onStat");
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

// Pythia8::PhaseSpace2to2tauyz file:Pythia8/PhaseSpace.h line:328
struct PyCallBack_Pythia8_PhaseSpace2to2tauyz : public Pythia8::PhaseSpace2to2tauyz {
	using Pythia8::PhaseSpace2to2tauyz::PhaseSpace2to2tauyz;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2tauyz::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2tauyz::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2tauyz::finalKin();
	}
	void rescaleMomenta(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "rescaleMomenta");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhaseSpace2to2tauyz::rescaleMomenta(a0);
	}
	void rescaleSigma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "rescaleSigma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhaseSpace2to2tauyz::rescaleSigma(a0);
	}
	double weightGammaPDFApprox() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "weightGammaPDFApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PhaseSpace2to2tauyz::weightGammaPDFApprox();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "sigmaSumSigned");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "isResolved");
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
	void setGammaKinPtr(class Pythia8::GammaKinematics * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2tauyz *>(this), "onStat");
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

// Pythia8::PhaseSpace2to2elastic file:Pythia8/PhaseSpace.h line:376
struct PyCallBack_Pythia8_PhaseSpace2to2elastic : public Pythia8::PhaseSpace2to2elastic {
	using Pythia8::PhaseSpace2to2elastic::PhaseSpace2to2elastic;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2elastic::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2elastic::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2elastic::finalKin();
	}
	bool isResolved() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "isResolved");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2elastic::isResolved();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "sigmaSumSigned");
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
	void rescaleSigma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "rescaleSigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "rescaleMomenta");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "weightGammaPDFApprox");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2elastic *>(this), "onStat");
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

// Pythia8::PhaseSpace2to2diffractive file:Pythia8/PhaseSpace.h line:412
struct PyCallBack_Pythia8_PhaseSpace2to2diffractive : public Pythia8::PhaseSpace2to2diffractive {
	using Pythia8::PhaseSpace2to2diffractive::PhaseSpace2to2diffractive;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2diffractive::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2diffractive::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2diffractive::finalKin();
	}
	bool isResolved() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "isResolved");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2diffractive::isResolved();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "sigmaSumSigned");
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
	void rescaleSigma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "rescaleSigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "rescaleMomenta");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "weightGammaPDFApprox");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2diffractive *>(this), "onStat");
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

// Pythia8::PhaseSpace2to3diffractive file:Pythia8/PhaseSpace.h line:459
struct PyCallBack_Pythia8_PhaseSpace2to3diffractive : public Pythia8::PhaseSpace2to3diffractive {
	using Pythia8::PhaseSpace2to3diffractive::PhaseSpace2to3diffractive;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3diffractive::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3diffractive::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3diffractive::finalKin();
	}
	bool isResolved() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "isResolved");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3diffractive::isResolved();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "sigmaSumSigned");
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
	void rescaleSigma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "rescaleSigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "rescaleMomenta");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "weightGammaPDFApprox");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3diffractive *>(this), "onStat");
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

// Pythia8::PhaseSpace2to2nondiffractive file:Pythia8/PhaseSpace.h line:495
struct PyCallBack_Pythia8_PhaseSpace2to2nondiffractive : public Pythia8::PhaseSpace2to2nondiffractive {
	using Pythia8::PhaseSpace2to2nondiffractive::PhaseSpace2to2nondiffractive;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2nondiffractive::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2nondiffractive::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to2nondiffractive::finalKin();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "sigmaSumSigned");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "isResolved");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "rescaleSigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "rescaleMomenta");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "weightGammaPDFApprox");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to2nondiffractive *>(this), "onStat");
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

// Pythia8::PhaseSpace2to3tauycyl file:Pythia8/PhaseSpace.h line:520
struct PyCallBack_Pythia8_PhaseSpace2to3tauycyl : public Pythia8::PhaseSpace2to3tauycyl {
	using Pythia8::PhaseSpace2to3tauycyl::PhaseSpace2to3tauycyl;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3tauycyl::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3tauycyl::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3tauycyl::finalKin();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "sigmaSumSigned");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "isResolved");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "rescaleSigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "rescaleMomenta");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "weightGammaPDFApprox");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3tauycyl *>(this), "onStat");
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

// Pythia8::PhaseSpace2to3yyycyl file:Pythia8/PhaseSpace.h line:558
struct PyCallBack_Pythia8_PhaseSpace2to3yyycyl : public Pythia8::PhaseSpace2to3yyycyl {
	using Pythia8::PhaseSpace2to3yyycyl::PhaseSpace2to3yyycyl;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3yyycyl::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3yyycyl::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpace2to3yyycyl::finalKin();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "sigmaSumSigned");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "isResolved");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "rescaleSigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "rescaleMomenta");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "weightGammaPDFApprox");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpace2to3yyycyl *>(this), "onStat");
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

// Pythia8::PhaseSpaceLHA file:Pythia8/PhaseSpace.h line:594
struct PyCallBack_Pythia8_PhaseSpaceLHA : public Pythia8::PhaseSpaceLHA {
	using Pythia8::PhaseSpaceLHA::PhaseSpaceLHA;

	bool setupSampling() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "setupSampling");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpaceLHA::setupSampling();
	}
	bool trialKin(bool a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "trialKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpaceLHA::trialKin(a0, a1);
	}
	bool finalKin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "finalKin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PhaseSpaceLHA::finalKin();
	}
	double sigmaSumSigned() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "sigmaSumSigned");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PhaseSpaceLHA::sigmaSumSigned();
	}
	bool isResolved() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "isResolved");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "rescaleSigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "rescaleMomenta");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "weightGammaPDFApprox");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "setGammaKinPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhaseSpaceLHA *>(this), "onStat");
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

void bind_Pythia8_PhaseSpace(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::PhaseSpace2to1tauy file:Pythia8/PhaseSpace.h line:299
		pybind11::class_<Pythia8::PhaseSpace2to1tauy, std::shared_ptr<Pythia8::PhaseSpace2to1tauy>, PyCallBack_Pythia8_PhaseSpace2to1tauy, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpace2to1tauy", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpace2to1tauy(); }, [](){ return new PyCallBack_Pythia8_PhaseSpace2to1tauy(); } ) );
		cl.def("setupSampling", (bool (Pythia8::PhaseSpace2to1tauy::*)()) &Pythia8::PhaseSpace2to1tauy::setupSampling, "C++: Pythia8::PhaseSpace2to1tauy::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to1tauy &o) -> bool { return o.trialKin(); }, "");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to1tauy &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg("inEvent"));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace2to1tauy::*)(bool, bool)) &Pythia8::PhaseSpace2to1tauy::trialKin, "C++: Pythia8::PhaseSpace2to1tauy::trialKin(bool, bool) --> bool", pybind11::arg("inEvent"), pybind11::arg(""));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace2to1tauy::*)()) &Pythia8::PhaseSpace2to1tauy::finalKin, "C++: Pythia8::PhaseSpace2to1tauy::finalKin() --> bool");
		cl.def("assign", (class Pythia8::PhaseSpace2to1tauy & (Pythia8::PhaseSpace2to1tauy::*)(const class Pythia8::PhaseSpace2to1tauy &)) &Pythia8::PhaseSpace2to1tauy::operator=, "C++: Pythia8::PhaseSpace2to1tauy::operator=(const class Pythia8::PhaseSpace2to1tauy &) --> class Pythia8::PhaseSpace2to1tauy &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpace2to2tauyz file:Pythia8/PhaseSpace.h line:328
		pybind11::class_<Pythia8::PhaseSpace2to2tauyz, std::shared_ptr<Pythia8::PhaseSpace2to2tauyz>, PyCallBack_Pythia8_PhaseSpace2to2tauyz, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpace2to2tauyz", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpace2to2tauyz(); }, [](){ return new PyCallBack_Pythia8_PhaseSpace2to2tauyz(); } ) );
		cl.def("setupSampling", (bool (Pythia8::PhaseSpace2to2tauyz::*)()) &Pythia8::PhaseSpace2to2tauyz::setupSampling, "C++: Pythia8::PhaseSpace2to2tauyz::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to2tauyz &o) -> bool { return o.trialKin(); }, "");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to2tauyz &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg("inEvent"));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace2to2tauyz::*)(bool, bool)) &Pythia8::PhaseSpace2to2tauyz::trialKin, "C++: Pythia8::PhaseSpace2to2tauyz::trialKin(bool, bool) --> bool", pybind11::arg("inEvent"), pybind11::arg(""));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace2to2tauyz::*)()) &Pythia8::PhaseSpace2to2tauyz::finalKin, "C++: Pythia8::PhaseSpace2to2tauyz::finalKin() --> bool");
		cl.def("rescaleMomenta", (void (Pythia8::PhaseSpace2to2tauyz::*)(double)) &Pythia8::PhaseSpace2to2tauyz::rescaleMomenta, "C++: Pythia8::PhaseSpace2to2tauyz::rescaleMomenta(double) --> void", pybind11::arg("sHatNew"));
		cl.def("rescaleSigma", (void (Pythia8::PhaseSpace2to2tauyz::*)(double)) &Pythia8::PhaseSpace2to2tauyz::rescaleSigma, "C++: Pythia8::PhaseSpace2to2tauyz::rescaleSigma(double) --> void", pybind11::arg("sHatNew"));
		cl.def("weightGammaPDFApprox", (double (Pythia8::PhaseSpace2to2tauyz::*)()) &Pythia8::PhaseSpace2to2tauyz::weightGammaPDFApprox, "C++: Pythia8::PhaseSpace2to2tauyz::weightGammaPDFApprox() --> double");
		cl.def("assign", (class Pythia8::PhaseSpace2to2tauyz & (Pythia8::PhaseSpace2to2tauyz::*)(const class Pythia8::PhaseSpace2to2tauyz &)) &Pythia8::PhaseSpace2to2tauyz::operator=, "C++: Pythia8::PhaseSpace2to2tauyz::operator=(const class Pythia8::PhaseSpace2to2tauyz &) --> class Pythia8::PhaseSpace2to2tauyz &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpace2to2elastic file:Pythia8/PhaseSpace.h line:376
		pybind11::class_<Pythia8::PhaseSpace2to2elastic, std::shared_ptr<Pythia8::PhaseSpace2to2elastic>, PyCallBack_Pythia8_PhaseSpace2to2elastic, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpace2to2elastic", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpace2to2elastic(); }, [](){ return new PyCallBack_Pythia8_PhaseSpace2to2elastic(); } ) );
		cl.def("setupSampling", (bool (Pythia8::PhaseSpace2to2elastic::*)()) &Pythia8::PhaseSpace2to2elastic::setupSampling, "C++: Pythia8::PhaseSpace2to2elastic::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to2elastic &o) -> bool { return o.trialKin(); }, "");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to2elastic &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg("inEvent"));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace2to2elastic::*)(bool, bool)) &Pythia8::PhaseSpace2to2elastic::trialKin, "C++: Pythia8::PhaseSpace2to2elastic::trialKin(bool, bool) --> bool", pybind11::arg("inEvent"), pybind11::arg(""));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace2to2elastic::*)()) &Pythia8::PhaseSpace2to2elastic::finalKin, "C++: Pythia8::PhaseSpace2to2elastic::finalKin() --> bool");
		cl.def("isResolved", (bool (Pythia8::PhaseSpace2to2elastic::*)() const) &Pythia8::PhaseSpace2to2elastic::isResolved, "C++: Pythia8::PhaseSpace2to2elastic::isResolved() const --> bool");
		cl.def("assign", (class Pythia8::PhaseSpace2to2elastic & (Pythia8::PhaseSpace2to2elastic::*)(const class Pythia8::PhaseSpace2to2elastic &)) &Pythia8::PhaseSpace2to2elastic::operator=, "C++: Pythia8::PhaseSpace2to2elastic::operator=(const class Pythia8::PhaseSpace2to2elastic &) --> class Pythia8::PhaseSpace2to2elastic &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpace2to2diffractive file:Pythia8/PhaseSpace.h line:412
		pybind11::class_<Pythia8::PhaseSpace2to2diffractive, std::shared_ptr<Pythia8::PhaseSpace2to2diffractive>, PyCallBack_Pythia8_PhaseSpace2to2diffractive, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpace2to2diffractive", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpace2to2diffractive(); }, [](){ return new PyCallBack_Pythia8_PhaseSpace2to2diffractive(); } ), "doc");
		cl.def( pybind11::init( [](bool const & a0){ return new Pythia8::PhaseSpace2to2diffractive(a0); }, [](bool const & a0){ return new PyCallBack_Pythia8_PhaseSpace2to2diffractive(a0); } ), "doc");
		cl.def( pybind11::init<bool, bool>(), pybind11::arg("isDiffAin"), pybind11::arg("isDiffBin") );

		cl.def("setupSampling", (bool (Pythia8::PhaseSpace2to2diffractive::*)()) &Pythia8::PhaseSpace2to2diffractive::setupSampling, "C++: Pythia8::PhaseSpace2to2diffractive::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to2diffractive &o) -> bool { return o.trialKin(); }, "");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to2diffractive &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg("inEvent"));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace2to2diffractive::*)(bool, bool)) &Pythia8::PhaseSpace2to2diffractive::trialKin, "C++: Pythia8::PhaseSpace2to2diffractive::trialKin(bool, bool) --> bool", pybind11::arg("inEvent"), pybind11::arg(""));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace2to2diffractive::*)()) &Pythia8::PhaseSpace2to2diffractive::finalKin, "C++: Pythia8::PhaseSpace2to2diffractive::finalKin() --> bool");
		cl.def("isResolved", (bool (Pythia8::PhaseSpace2to2diffractive::*)() const) &Pythia8::PhaseSpace2to2diffractive::isResolved, "C++: Pythia8::PhaseSpace2to2diffractive::isResolved() const --> bool");
		cl.def("assign", (class Pythia8::PhaseSpace2to2diffractive & (Pythia8::PhaseSpace2to2diffractive::*)(const class Pythia8::PhaseSpace2to2diffractive &)) &Pythia8::PhaseSpace2to2diffractive::operator=, "C++: Pythia8::PhaseSpace2to2diffractive::operator=(const class Pythia8::PhaseSpace2to2diffractive &) --> class Pythia8::PhaseSpace2to2diffractive &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpace2to3diffractive file:Pythia8/PhaseSpace.h line:459
		pybind11::class_<Pythia8::PhaseSpace2to3diffractive, std::shared_ptr<Pythia8::PhaseSpace2to3diffractive>, PyCallBack_Pythia8_PhaseSpace2to3diffractive, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpace2to3diffractive", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpace2to3diffractive(); }, [](){ return new PyCallBack_Pythia8_PhaseSpace2to3diffractive(); } ) );
		cl.def("setupSampling", (bool (Pythia8::PhaseSpace2to3diffractive::*)()) &Pythia8::PhaseSpace2to3diffractive::setupSampling, "C++: Pythia8::PhaseSpace2to3diffractive::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to3diffractive &o) -> bool { return o.trialKin(); }, "");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to3diffractive &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg("inEvent"));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace2to3diffractive::*)(bool, bool)) &Pythia8::PhaseSpace2to3diffractive::trialKin, "C++: Pythia8::PhaseSpace2to3diffractive::trialKin(bool, bool) --> bool", pybind11::arg("inEvent"), pybind11::arg(""));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace2to3diffractive::*)()) &Pythia8::PhaseSpace2to3diffractive::finalKin, "C++: Pythia8::PhaseSpace2to3diffractive::finalKin() --> bool");
		cl.def("isResolved", (bool (Pythia8::PhaseSpace2to3diffractive::*)() const) &Pythia8::PhaseSpace2to3diffractive::isResolved, "C++: Pythia8::PhaseSpace2to3diffractive::isResolved() const --> bool");
		cl.def("assign", (class Pythia8::PhaseSpace2to3diffractive & (Pythia8::PhaseSpace2to3diffractive::*)(const class Pythia8::PhaseSpace2to3diffractive &)) &Pythia8::PhaseSpace2to3diffractive::operator=, "C++: Pythia8::PhaseSpace2to3diffractive::operator=(const class Pythia8::PhaseSpace2to3diffractive &) --> class Pythia8::PhaseSpace2to3diffractive &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpace2to2nondiffractive file:Pythia8/PhaseSpace.h line:495
		pybind11::class_<Pythia8::PhaseSpace2to2nondiffractive, std::shared_ptr<Pythia8::PhaseSpace2to2nondiffractive>, PyCallBack_Pythia8_PhaseSpace2to2nondiffractive, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpace2to2nondiffractive", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpace2to2nondiffractive(); }, [](){ return new PyCallBack_Pythia8_PhaseSpace2to2nondiffractive(); } ) );
		cl.def("setupSampling", (bool (Pythia8::PhaseSpace2to2nondiffractive::*)()) &Pythia8::PhaseSpace2to2nondiffractive::setupSampling, "C++: Pythia8::PhaseSpace2to2nondiffractive::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to2nondiffractive &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg(""));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace2to2nondiffractive::*)(bool, bool)) &Pythia8::PhaseSpace2to2nondiffractive::trialKin, "C++: Pythia8::PhaseSpace2to2nondiffractive::trialKin(bool, bool) --> bool", pybind11::arg(""), pybind11::arg(""));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace2to2nondiffractive::*)()) &Pythia8::PhaseSpace2to2nondiffractive::finalKin, "C++: Pythia8::PhaseSpace2to2nondiffractive::finalKin() --> bool");
		cl.def("assign", (class Pythia8::PhaseSpace2to2nondiffractive & (Pythia8::PhaseSpace2to2nondiffractive::*)(const class Pythia8::PhaseSpace2to2nondiffractive &)) &Pythia8::PhaseSpace2to2nondiffractive::operator=, "C++: Pythia8::PhaseSpace2to2nondiffractive::operator=(const class Pythia8::PhaseSpace2to2nondiffractive &) --> class Pythia8::PhaseSpace2to2nondiffractive &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpace2to3tauycyl file:Pythia8/PhaseSpace.h line:520
		pybind11::class_<Pythia8::PhaseSpace2to3tauycyl, std::shared_ptr<Pythia8::PhaseSpace2to3tauycyl>, PyCallBack_Pythia8_PhaseSpace2to3tauycyl, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpace2to3tauycyl", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpace2to3tauycyl(); }, [](){ return new PyCallBack_Pythia8_PhaseSpace2to3tauycyl(); } ) );
		cl.def("setupSampling", (bool (Pythia8::PhaseSpace2to3tauycyl::*)()) &Pythia8::PhaseSpace2to3tauycyl::setupSampling, "C++: Pythia8::PhaseSpace2to3tauycyl::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to3tauycyl &o) -> bool { return o.trialKin(); }, "");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to3tauycyl &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg("inEvent"));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace2to3tauycyl::*)(bool, bool)) &Pythia8::PhaseSpace2to3tauycyl::trialKin, "C++: Pythia8::PhaseSpace2to3tauycyl::trialKin(bool, bool) --> bool", pybind11::arg("inEvent"), pybind11::arg(""));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace2to3tauycyl::*)()) &Pythia8::PhaseSpace2to3tauycyl::finalKin, "C++: Pythia8::PhaseSpace2to3tauycyl::finalKin() --> bool");
		cl.def("assign", (class Pythia8::PhaseSpace2to3tauycyl & (Pythia8::PhaseSpace2to3tauycyl::*)(const class Pythia8::PhaseSpace2to3tauycyl &)) &Pythia8::PhaseSpace2to3tauycyl::operator=, "C++: Pythia8::PhaseSpace2to3tauycyl::operator=(const class Pythia8::PhaseSpace2to3tauycyl &) --> class Pythia8::PhaseSpace2to3tauycyl &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpace2to3yyycyl file:Pythia8/PhaseSpace.h line:558
		pybind11::class_<Pythia8::PhaseSpace2to3yyycyl, std::shared_ptr<Pythia8::PhaseSpace2to3yyycyl>, PyCallBack_Pythia8_PhaseSpace2to3yyycyl, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpace2to3yyycyl", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpace2to3yyycyl(); }, [](){ return new PyCallBack_Pythia8_PhaseSpace2to3yyycyl(); } ) );
		cl.def("setupSampling", (bool (Pythia8::PhaseSpace2to3yyycyl::*)()) &Pythia8::PhaseSpace2to3yyycyl::setupSampling, "C++: Pythia8::PhaseSpace2to3yyycyl::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to3yyycyl &o) -> bool { return o.trialKin(); }, "");
		cl.def("trialKin", [](Pythia8::PhaseSpace2to3yyycyl &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg("inEvent"));
		cl.def("trialKin", (bool (Pythia8::PhaseSpace2to3yyycyl::*)(bool, bool)) &Pythia8::PhaseSpace2to3yyycyl::trialKin, "C++: Pythia8::PhaseSpace2to3yyycyl::trialKin(bool, bool) --> bool", pybind11::arg("inEvent"), pybind11::arg(""));
		cl.def("finalKin", (bool (Pythia8::PhaseSpace2to3yyycyl::*)()) &Pythia8::PhaseSpace2to3yyycyl::finalKin, "C++: Pythia8::PhaseSpace2to3yyycyl::finalKin() --> bool");
		cl.def("assign", (class Pythia8::PhaseSpace2to3yyycyl & (Pythia8::PhaseSpace2to3yyycyl::*)(const class Pythia8::PhaseSpace2to3yyycyl &)) &Pythia8::PhaseSpace2to3yyycyl::operator=, "C++: Pythia8::PhaseSpace2to3yyycyl::operator=(const class Pythia8::PhaseSpace2to3yyycyl &) --> class Pythia8::PhaseSpace2to3yyycyl &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PhaseSpaceLHA file:Pythia8/PhaseSpace.h line:594
		pybind11::class_<Pythia8::PhaseSpaceLHA, std::shared_ptr<Pythia8::PhaseSpaceLHA>, PyCallBack_Pythia8_PhaseSpaceLHA, Pythia8::PhaseSpace> cl(M("Pythia8"), "PhaseSpaceLHA", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhaseSpaceLHA(); }, [](){ return new PyCallBack_Pythia8_PhaseSpaceLHA(); } ) );
		cl.def("setupSampling", (bool (Pythia8::PhaseSpaceLHA::*)()) &Pythia8::PhaseSpaceLHA::setupSampling, "C++: Pythia8::PhaseSpaceLHA::setupSampling() --> bool");
		cl.def("trialKin", [](Pythia8::PhaseSpaceLHA &o, bool const & a0) -> bool { return o.trialKin(a0); }, "", pybind11::arg(""));
		cl.def("trialKin", (bool (Pythia8::PhaseSpaceLHA::*)(bool, bool)) &Pythia8::PhaseSpaceLHA::trialKin, "C++: Pythia8::PhaseSpaceLHA::trialKin(bool, bool) --> bool", pybind11::arg(""), pybind11::arg("repeatSame"));
		cl.def("finalKin", (bool (Pythia8::PhaseSpaceLHA::*)()) &Pythia8::PhaseSpaceLHA::finalKin, "C++: Pythia8::PhaseSpaceLHA::finalKin() --> bool");
		cl.def("sigmaSumSigned", (double (Pythia8::PhaseSpaceLHA::*)() const) &Pythia8::PhaseSpaceLHA::sigmaSumSigned, "C++: Pythia8::PhaseSpaceLHA::sigmaSumSigned() const --> double");
		cl.def("assign", (class Pythia8::PhaseSpaceLHA & (Pythia8::PhaseSpaceLHA::*)(const class Pythia8::PhaseSpaceLHA &)) &Pythia8::PhaseSpaceLHA::operator=, "C++: Pythia8::PhaseSpaceLHA::operator=(const class Pythia8::PhaseSpaceLHA &) --> class Pythia8::PhaseSpaceLHA &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
