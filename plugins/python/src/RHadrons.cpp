#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/HadronLevel.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/NucleonExcitations.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ParticleDecays.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/RHadrons.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
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

// Pythia8::RHadrons file:Pythia8/RHadrons.h line:21
struct PyCallBack_Pythia8_RHadrons : public Pythia8::RHadrons {
	using Pythia8::RHadrons::RHadrons;

	bool init(class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2, class std::shared_ptr<class Pythia8::FragmentationModifierBase> a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RHadrons *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return RHadrons::init(a0, a1, a2, a3);
	}
	bool fragment(int a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RHadrons *>(this), "fragment");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return RHadrons::fragment(a0, a1, a2, a3, a4);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RHadrons *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RHadrons *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RHadrons *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RHadrons *>(this), "onStat");
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

// Pythia8::HadronLevel file:Pythia8/HadronLevel.h line:44
struct PyCallBack_Pythia8_HadronLevel : public Pythia8::HadronLevel {
	using Pythia8::HadronLevel::HadronLevel;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HadronLevel *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HadronLevel::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HadronLevel *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HadronLevel *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HadronLevel *>(this), "onStat");
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

// Pythia8::StringInteractions file:Pythia8/StringInteractions.h line:28
struct PyCallBack_Pythia8_StringInteractions : public Pythia8::StringInteractions {
	using Pythia8::StringInteractions::StringInteractions;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringInteractions *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringInteractions::init();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringInteractions *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringInteractions *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringInteractions *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringInteractions *>(this), "onStat");
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

// Pythia8::ColourReconnectionBase file:Pythia8/StringInteractions.h line:78
struct PyCallBack_Pythia8_ColourReconnectionBase : public Pythia8::ColourReconnectionBase {
	using Pythia8::ColourReconnectionBase::ColourReconnectionBase;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnectionBase *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ColourReconnectionBase::init();
	}
	void reassignBeamPtrs(class Pythia8::BeamParticle * a0, class Pythia8::BeamParticle * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnectionBase *>(this), "reassignBeamPtrs");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColourReconnectionBase::reassignBeamPtrs(a0, a1);
	}
	bool next(class Pythia8::Event & a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnectionBase *>(this), "next");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ColourReconnectionBase::next\"");
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnectionBase *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnectionBase *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnectionBase *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnectionBase *>(this), "onStat");
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

// Pythia8::DipoleSwingBase file:Pythia8/StringInteractions.h line:106
struct PyCallBack_Pythia8_DipoleSwingBase : public Pythia8::DipoleSwingBase {
	using Pythia8::DipoleSwingBase::DipoleSwingBase;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return DipoleSwingBase::init();
	}
	void reassignBeamPtrs(class Pythia8::BeamParticle * a0, class Pythia8::BeamParticle * a1, int a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "reassignBeamPtrs");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DipoleSwingBase::reassignBeamPtrs(a0, a1, a2);
	}
	void prepare(int a0, class Pythia8::Event & a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "prepare");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"DipoleSwingBase::prepare\"");
	}
	void rescatterUpdate(int a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "rescatterUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"DipoleSwingBase::rescatterUpdate\"");
	}
	void update(int a0, class Pythia8::Event & a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"DipoleSwingBase::update\"");
	}
	double pTnext(class Pythia8::Event & a0, double a1, double a2, bool a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "pTnext");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"DipoleSwingBase::pTnext\"");
	}
	bool swing(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "swing");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"DipoleSwingBase::swing\"");
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DipoleSwingBase *>(this), "onStat");
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

void bind_Pythia8_RHadrons(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::RHadrons file:Pythia8/RHadrons.h line:21
		pybind11::class_<Pythia8::RHadrons, std::shared_ptr<Pythia8::RHadrons>, PyCallBack_Pythia8_RHadrons, Pythia8::FragmentationModel> cl(M("Pythia8"), "RHadrons", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::RHadrons(); }, [](){ return new PyCallBack_Pythia8_RHadrons(); } ) );
		cl.def("init", [](Pythia8::RHadrons &o) -> bool { return o.init(); }, "");
		cl.def("init", [](Pythia8::RHadrons &o, class Pythia8::StringFlav * a0) -> bool { return o.init(a0); }, "", pybind11::arg("flavSelPtrIn"));
		cl.def("init", [](Pythia8::RHadrons &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1) -> bool { return o.init(a0, a1); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"));
		cl.def("init", [](Pythia8::RHadrons &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2) -> bool { return o.init(a0, a1, a2); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"));
		cl.def("init", (bool (Pythia8::RHadrons::*)(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>)) &Pythia8::RHadrons::init, "C++: Pythia8::RHadrons::init(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>) --> bool", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"), pybind11::arg("fragModPtrIn"));
		cl.def("fragment", [](Pythia8::RHadrons &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2) -> bool { return o.fragment(a0, a1, a2); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"));
		cl.def("fragment", [](Pythia8::RHadrons &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool const & a3) -> bool { return o.fragment(a0, a1, a2, a3); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"));
		cl.def("fragment", (bool (Pythia8::RHadrons::*)(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool)) &Pythia8::RHadrons::fragment, "C++: Pythia8::RHadrons::fragment(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool) --> bool", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"), pybind11::arg("systemRecoil"));
		cl.def("decay", (bool (Pythia8::RHadrons::*)(class Pythia8::Event &)) &Pythia8::RHadrons::decay, "C++: Pythia8::RHadrons::decay(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("givesRHadron", (bool (Pythia8::RHadrons::*)(int)) &Pythia8::RHadrons::givesRHadron, "C++: Pythia8::RHadrons::givesRHadron(int) --> bool", pybind11::arg("id"));
		cl.def("exist", (bool (Pythia8::RHadrons::*)()) &Pythia8::RHadrons::exist, "C++: Pythia8::RHadrons::exist() --> bool");
		cl.def("trace", (int (Pythia8::RHadrons::*)(int)) &Pythia8::RHadrons::trace, "C++: Pythia8::RHadrons::trace(int) --> int", pybind11::arg("i"));
		cl.def("assign", (class Pythia8::RHadrons & (Pythia8::RHadrons::*)(const class Pythia8::RHadrons &)) &Pythia8::RHadrons::operator=, "C++: Pythia8::RHadrons::operator=(const class Pythia8::RHadrons &) --> class Pythia8::RHadrons &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HadronLevel file:Pythia8/HadronLevel.h line:44
		pybind11::class_<Pythia8::HadronLevel, std::shared_ptr<Pythia8::HadronLevel>, PyCallBack_Pythia8_HadronLevel, Pythia8::PhysicsBase> cl(M("Pythia8"), "HadronLevel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HadronLevel(); }, [](){ return new PyCallBack_Pythia8_HadronLevel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HadronLevel const &o){ return new PyCallBack_Pythia8_HadronLevel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HadronLevel const &o){ return new Pythia8::HadronLevel(o); } ) );
		cl.def("init", (bool (Pythia8::HadronLevel::*)(class std::shared_ptr<class Pythia8::TimeShower>, class std::shared_ptr<class Pythia8::RHadrons>, class std::shared_ptr<class Pythia8::LundFragmentation>, class std::vector<class std::shared_ptr<class Pythia8::FragmentationModel>, class std::allocator<class std::shared_ptr<class Pythia8::FragmentationModel> > > *, class std::shared_ptr<class Pythia8::DecayHandler>, class std::vector<int, class std::allocator<int> >, class std::shared_ptr<class Pythia8::StringInteractions>, class std::shared_ptr<class Pythia8::PartonVertex>, class Pythia8::SigmaLowEnergy &, class Pythia8::NucleonExcitations &)) &Pythia8::HadronLevel::init, "C++: Pythia8::HadronLevel::init(class std::shared_ptr<class Pythia8::TimeShower>, class std::shared_ptr<class Pythia8::RHadrons>, class std::shared_ptr<class Pythia8::LundFragmentation>, class std::vector<class std::shared_ptr<class Pythia8::FragmentationModel>, class std::allocator<class std::shared_ptr<class Pythia8::FragmentationModel> > > *, class std::shared_ptr<class Pythia8::DecayHandler>, class std::vector<int, class std::allocator<int> >, class std::shared_ptr<class Pythia8::StringInteractions>, class std::shared_ptr<class Pythia8::PartonVertex>, class Pythia8::SigmaLowEnergy &, class Pythia8::NucleonExcitations &) --> bool", pybind11::arg("timesDecPtrIn"), pybind11::arg("rHadronsPtrIn"), pybind11::arg("fragPtrIn"), pybind11::arg("fragPtrsIn"), pybind11::arg("decayHandlePtr"), pybind11::arg("handledParticles"), pybind11::arg("stringInteractionsPtrIn"), pybind11::arg("partonVertexPtrIn"), pybind11::arg("sigmaLowEnergyIn"), pybind11::arg("nucleonExcitationsIn"));
		cl.def("getStringFlavPtr", (class Pythia8::StringFlav * (Pythia8::HadronLevel::*)()) &Pythia8::HadronLevel::getStringFlavPtr, "C++: Pythia8::HadronLevel::getStringFlavPtr() --> class Pythia8::StringFlav *", pybind11::return_value_policy::automatic);
		cl.def("next", (bool (Pythia8::HadronLevel::*)(class Pythia8::Event &)) &Pythia8::HadronLevel::next, "C++: Pythia8::HadronLevel::next(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("decay", (bool (Pythia8::HadronLevel::*)(int, class Pythia8::Event &)) &Pythia8::HadronLevel::decay, "C++: Pythia8::HadronLevel::decay(int, class Pythia8::Event &) --> bool", pybind11::arg("iDec"), pybind11::arg("event"));
		cl.def("moreDecays", (bool (Pythia8::HadronLevel::*)(class Pythia8::Event &)) &Pythia8::HadronLevel::moreDecays, "C++: Pythia8::HadronLevel::moreDecays(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("rescatter", (bool (Pythia8::HadronLevel::*)(class Pythia8::Event &)) &Pythia8::HadronLevel::rescatter, "C++: Pythia8::HadronLevel::rescatter(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("initLowEnergyProcesses", (bool (Pythia8::HadronLevel::*)()) &Pythia8::HadronLevel::initLowEnergyProcesses, "C++: Pythia8::HadronLevel::initLowEnergyProcesses() --> bool");
		cl.def("pickLowEnergyProcess", (int (Pythia8::HadronLevel::*)(int, int, double, double, double)) &Pythia8::HadronLevel::pickLowEnergyProcess, "C++: Pythia8::HadronLevel::pickLowEnergyProcess(int, int, double, double, double) --> int", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"), pybind11::arg("mA"), pybind11::arg("mB"));
		cl.def("doLowEnergyProcess", (bool (Pythia8::HadronLevel::*)(int, int, int, class Pythia8::Event &)) &Pythia8::HadronLevel::doLowEnergyProcess, "C++: Pythia8::HadronLevel::doLowEnergyProcess(int, int, int, class Pythia8::Event &) --> bool", pybind11::arg("i1"), pybind11::arg("i2"), pybind11::arg("procTypeIn"), pybind11::arg("event"));
		cl.def("hasVetoedHadronize", (bool (Pythia8::HadronLevel::*)() const) &Pythia8::HadronLevel::hasVetoedHadronize, "C++: Pythia8::HadronLevel::hasVetoedHadronize() const --> bool");
		cl.def("onInitInfoPtr", (void (Pythia8::HadronLevel::*)()) &Pythia8::HadronLevel::onInitInfoPtr, "C++: Pythia8::HadronLevel::onInitInfoPtr() --> void");
		cl.def("assign", (class Pythia8::HadronLevel & (Pythia8::HadronLevel::*)(const class Pythia8::HadronLevel &)) &Pythia8::HadronLevel::operator=, "C++: Pythia8::HadronLevel::operator=(const class Pythia8::HadronLevel &) --> class Pythia8::HadronLevel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::StringInteractions file:Pythia8/StringInteractions.h line:28
		pybind11::class_<Pythia8::StringInteractions, std::shared_ptr<Pythia8::StringInteractions>, PyCallBack_Pythia8_StringInteractions, Pythia8::PhysicsBase> cl(M("Pythia8"), "StringInteractions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::StringInteractions(); }, [](){ return new PyCallBack_Pythia8_StringInteractions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_StringInteractions const &o){ return new PyCallBack_Pythia8_StringInteractions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::StringInteractions const &o){ return new Pythia8::StringInteractions(o); } ) );
		cl.def_readwrite("colrecPtr", &Pythia8::StringInteractions::colrecPtr);
		cl.def_readwrite("dipswingPtr", &Pythia8::StringInteractions::dipswingPtr);
		cl.def_readwrite("stringrepPtr", &Pythia8::StringInteractions::stringrepPtr);
		cl.def_readwrite("fragmodPtr", &Pythia8::StringInteractions::fragmodPtr);
		cl.def("init", (bool (Pythia8::StringInteractions::*)()) &Pythia8::StringInteractions::init, "C++: Pythia8::StringInteractions::init() --> bool");
		cl.def("getColourReconnections", (class std::shared_ptr<class Pythia8::ColourReconnectionBase> (Pythia8::StringInteractions::*)() const) &Pythia8::StringInteractions::getColourReconnections, "C++: Pythia8::StringInteractions::getColourReconnections() const --> class std::shared_ptr<class Pythia8::ColourReconnectionBase>");
		cl.def("getDipoleSwing", (class std::shared_ptr<class Pythia8::DipoleSwingBase> (Pythia8::StringInteractions::*)() const) &Pythia8::StringInteractions::getDipoleSwing, "C++: Pythia8::StringInteractions::getDipoleSwing() const --> class std::shared_ptr<class Pythia8::DipoleSwingBase>");
		cl.def("getStringRepulsion", (class std::shared_ptr<class Pythia8::StringRepulsionBase> (Pythia8::StringInteractions::*)() const) &Pythia8::StringInteractions::getStringRepulsion, "C++: Pythia8::StringInteractions::getStringRepulsion() const --> class std::shared_ptr<class Pythia8::StringRepulsionBase>");
		cl.def("getFragmentationModifier", (class std::shared_ptr<class Pythia8::FragmentationModifierBase> (Pythia8::StringInteractions::*)() const) &Pythia8::StringInteractions::getFragmentationModifier, "C++: Pythia8::StringInteractions::getFragmentationModifier() const --> class std::shared_ptr<class Pythia8::FragmentationModifierBase>");
		cl.def("assign", (class Pythia8::StringInteractions & (Pythia8::StringInteractions::*)(const class Pythia8::StringInteractions &)) &Pythia8::StringInteractions::operator=, "C++: Pythia8::StringInteractions::operator=(const class Pythia8::StringInteractions &) --> class Pythia8::StringInteractions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ColourReconnectionBase file:Pythia8/StringInteractions.h line:78
		pybind11::class_<Pythia8::ColourReconnectionBase, std::shared_ptr<Pythia8::ColourReconnectionBase>, PyCallBack_Pythia8_ColourReconnectionBase, Pythia8::PhysicsBase> cl(M("Pythia8"), "ColourReconnectionBase", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_ColourReconnectionBase(); } ) );
		cl.def(pybind11::init<PyCallBack_Pythia8_ColourReconnectionBase const &>());
		cl.def("init", (bool (Pythia8::ColourReconnectionBase::*)()) &Pythia8::ColourReconnectionBase::init, "C++: Pythia8::ColourReconnectionBase::init() --> bool");
		cl.def("reassignBeamPtrs", (void (Pythia8::ColourReconnectionBase::*)(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *)) &Pythia8::ColourReconnectionBase::reassignBeamPtrs, "C++: Pythia8::ColourReconnectionBase::reassignBeamPtrs(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *) --> void", pybind11::arg("beamAPtrIn"), pybind11::arg("beamBPtrIn"));
		cl.def("next", (bool (Pythia8::ColourReconnectionBase::*)(class Pythia8::Event &, int)) &Pythia8::ColourReconnectionBase::next, "C++: Pythia8::ColourReconnectionBase::next(class Pythia8::Event &, int) --> bool", pybind11::arg("event"), pybind11::arg("oldSize"));
		cl.def("assign", (class Pythia8::ColourReconnectionBase & (Pythia8::ColourReconnectionBase::*)(const class Pythia8::ColourReconnectionBase &)) &Pythia8::ColourReconnectionBase::operator=, "C++: Pythia8::ColourReconnectionBase::operator=(const class Pythia8::ColourReconnectionBase &) --> class Pythia8::ColourReconnectionBase &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::DipoleSwingBase file:Pythia8/StringInteractions.h line:106
		pybind11::class_<Pythia8::DipoleSwingBase, std::shared_ptr<Pythia8::DipoleSwingBase>, PyCallBack_Pythia8_DipoleSwingBase, Pythia8::PhysicsBase> cl(M("Pythia8"), "DipoleSwingBase", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_DipoleSwingBase(); } ) );
		cl.def_readwrite("beamOffset", &Pythia8::DipoleSwingBase::beamOffset);
		cl.def("init", (bool (Pythia8::DipoleSwingBase::*)()) &Pythia8::DipoleSwingBase::init, "C++: Pythia8::DipoleSwingBase::init() --> bool");
		cl.def("reassignBeamPtrs", [](Pythia8::DipoleSwingBase &o, class Pythia8::BeamParticle * a0, class Pythia8::BeamParticle * a1) -> void { return o.reassignBeamPtrs(a0, a1); }, "", pybind11::arg("beamAPtrIn"), pybind11::arg("beamBPtrIn"));
		cl.def("reassignBeamPtrs", (void (Pythia8::DipoleSwingBase::*)(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *, int)) &Pythia8::DipoleSwingBase::reassignBeamPtrs, "C++: Pythia8::DipoleSwingBase::reassignBeamPtrs(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *, int) --> void", pybind11::arg("beamAPtrIn"), pybind11::arg("beamBPtrIn"), pybind11::arg("beamOffsetIn"));
		cl.def("prepare", [](Pythia8::DipoleSwingBase &o, int const & a0, class Pythia8::Event & a1) -> void { return o.prepare(a0, a1); }, "", pybind11::arg(""), pybind11::arg(""));
		cl.def("prepare", (void (Pythia8::DipoleSwingBase::*)(int, class Pythia8::Event &, bool)) &Pythia8::DipoleSwingBase::prepare, "C++: Pythia8::DipoleSwingBase::prepare(int, class Pythia8::Event &, bool) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("rescatterUpdate", (void (Pythia8::DipoleSwingBase::*)(int, class Pythia8::Event &)) &Pythia8::DipoleSwingBase::rescatterUpdate, "C++: Pythia8::DipoleSwingBase::rescatterUpdate(int, class Pythia8::Event &) --> void", pybind11::arg(""), pybind11::arg(""));
		cl.def("update", [](Pythia8::DipoleSwingBase &o, int const & a0, class Pythia8::Event & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg(""), pybind11::arg(""));
		cl.def("update", (void (Pythia8::DipoleSwingBase::*)(int, class Pythia8::Event &, bool)) &Pythia8::DipoleSwingBase::update, "C++: Pythia8::DipoleSwingBase::update(int, class Pythia8::Event &, bool) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("pTnext", [](Pythia8::DipoleSwingBase &o, class Pythia8::Event & a0, double const & a1, double const & a2) -> double { return o.pTnext(a0, a1, a2); }, "", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("pTnext", [](Pythia8::DipoleSwingBase &o, class Pythia8::Event & a0, double const & a1, double const & a2, bool const & a3) -> double { return o.pTnext(a0, a1, a2, a3); }, "", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("pTnext", (double (Pythia8::DipoleSwingBase::*)(class Pythia8::Event &, double, double, bool, bool)) &Pythia8::DipoleSwingBase::pTnext, "C++: Pythia8::DipoleSwingBase::pTnext(class Pythia8::Event &, double, double, bool, bool) --> double", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("swing", (bool (Pythia8::DipoleSwingBase::*)(class Pythia8::Event &)) &Pythia8::DipoleSwingBase::swing, "C++: Pythia8::DipoleSwingBase::swing(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("assign", (class Pythia8::DipoleSwingBase & (Pythia8::DipoleSwingBase::*)(const class Pythia8::DipoleSwingBase &)) &Pythia8::DipoleSwingBase::operator=, "C++: Pythia8::DipoleSwingBase::operator=(const class Pythia8::DipoleSwingBase &) --> class Pythia8::DipoleSwingBase &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
