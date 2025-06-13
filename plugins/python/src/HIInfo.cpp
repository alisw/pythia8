#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/HIBasics.h>
#include <Pythia8/HIInfo.h>
#include <Pythia8/HINucleusModel.h>
#include <Pythia8/HISubCollisionModel.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/HeavyIons.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/Merging.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ParticleDecays.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/ShowerModel.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaProcess.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/SusyLesHouches.h>
#include <Pythia8/UserHooks.h>
#include <Pythia8/Weights.h>
#include <complex>
#include <cwchar>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
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

// Pythia8::HIUserHooks file:Pythia8/HIInfo.h line:313
struct PyCallBack_Pythia8_HIUserHooks : public Pythia8::HIUserHooks {
	using Pythia8::HIUserHooks::HIUserHooks;

	void init(int a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HIUserHooks::init(a0, a1);
	}
	bool hasImpactParameterGenerator() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "hasImpactParameterGenerator");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::hasImpactParameterGenerator();
	}
	class std::shared_ptr<class Pythia8::ImpactParameterGenerator> impactParameterGenerator() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "impactParameterGenerator");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::ImpactParameterGenerator>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::ImpactParameterGenerator>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::ImpactParameterGenerator>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::ImpactParameterGenerator>>(std::move(o));
		}
		return HIUserHooks::impactParameterGenerator();
	}
	bool hasProjectileModel() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "hasProjectileModel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::hasProjectileModel();
	}
	class std::shared_ptr<class Pythia8::NucleusModel> projectileModel() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "projectileModel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::NucleusModel>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::NucleusModel>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::NucleusModel>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::NucleusModel>>(std::move(o));
		}
		return HIUserHooks::projectileModel();
	}
	bool hasTargetModel() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "hasTargetModel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::hasTargetModel();
	}
	class std::shared_ptr<class Pythia8::NucleusModel> targetModel() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "targetModel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::NucleusModel>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::NucleusModel>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::NucleusModel>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::NucleusModel>>(std::move(o));
		}
		return HIUserHooks::targetModel();
	}
	bool hasSubCollisionModel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "hasSubCollisionModel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::hasSubCollisionModel();
	}
	class std::shared_ptr<class Pythia8::SubCollisionModel> subCollisionModel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "subCollisionModel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::SubCollisionModel>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::SubCollisionModel>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::SubCollisionModel>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::SubCollisionModel>>(std::move(o));
		}
		return HIUserHooks::subCollisionModel();
	}
	bool hasEventOrdering() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "hasEventOrdering");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::hasEventOrdering();
	}
	double eventOrdering(const class Pythia8::Event & a0, const class Pythia8::Info & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "eventOrdering");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HIUserHooks::eventOrdering(a0, a1);
	}
	bool canFixIsoSpin() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "canFixIsoSpin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::canFixIsoSpin();
	}
	bool fixIsoSpin(class Pythia8::EventInfo & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "fixIsoSpin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::fixIsoSpin(a0);
	}
	bool canShiftEvent() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "canShiftEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::canShiftEvent();
	}
	class Pythia8::EventInfo & shiftEvent(class Pythia8::EventInfo & a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "shiftEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::EventInfo &>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::EventInfo &> caster;
				return pybind11::detail::cast_ref<class Pythia8::EventInfo &>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::EventInfo &>(std::move(o));
		}
		return HIUserHooks::shiftEvent(a0);
	}
	bool canForceHadronLevel() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "canForceHadronLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::canForceHadronLevel();
	}
	bool forceHadronLevel(class Pythia8::Pythia & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "forceHadronLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::forceHadronLevel(a0);
	}
	bool canFindRecoilers() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "canFindRecoilers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HIUserHooks::canFindRecoilers();
	}
	using _binder_ret_0 = class std::vector<int, class std::allocator<int> >;
	_binder_ret_0 findRecoilers(const class Pythia8::Event & a0, bool a1, int a2, int a3, const class Pythia8::Vec4 & a4, const class Pythia8::Vec4 & a5) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HIUserHooks *>(this), "findRecoilers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return HIUserHooks::findRecoilers(a0, a1, a2, a3, a4, a5);
	}
};

// Pythia8::HeavyIons file:Pythia8/HeavyIons.h line:31
struct PyCallBack_Pythia8_HeavyIons : public Pythia8::HeavyIons {
	using Pythia8::HeavyIons::HeavyIons;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"HeavyIons::init\"");
	}
	bool next() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "next");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"HeavyIons::next\"");
	}
	bool setKinematics(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "setKinematics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HeavyIons::setKinematics(a0);
	}
	bool setKinematics(double a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "setKinematics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HeavyIons::setKinematics(a0, a1);
	}
	bool setKinematics(double a0, double a1, double a2, double a3, double a4, double a5) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "setKinematics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HeavyIons::setKinematics(a0, a1, a2, a3, a4, a5);
	}
	bool setKinematics(class Pythia8::Vec4 a0, class Pythia8::Vec4 a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "setKinematics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HeavyIons::setKinematics(a0, a1);
	}
	bool setBeamIDs(int a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "setBeamIDs");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HeavyIons::setBeamIDs(a0, a1);
	}
	void stat() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "stat");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HeavyIons::stat();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons *>(this), "onStat");
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

// Pythia8::HeavyIons::InfoGrabber file:Pythia8/HeavyIons.h line:147
struct PyCallBack_Pythia8_HeavyIons_InfoGrabber : public Pythia8::HeavyIons::InfoGrabber {
	using Pythia8::HeavyIons::InfoGrabber::InfoGrabber;

	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::initAfterBeams();
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canModifySigma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canModifySigma();
	}
	double multiplySigmaBy(const class Pythia8::SigmaProcess * a0, const class Pythia8::PhaseSpace * a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "multiplySigmaBy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return UserHooks::multiplySigmaBy(a0, a1, a2);
	}
	bool canBiasSelection() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canBiasSelection");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canBiasSelection();
	}
	double biasSelectionBy(const class Pythia8::SigmaProcess * a0, const class Pythia8::PhaseSpace * a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "biasSelectionBy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return UserHooks::biasSelectionBy(a0, a1, a2);
	}
	double biasedSelectionWeight() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "biasedSelectionWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return UserHooks::biasedSelectionWeight();
	}
	bool canVetoProcessLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoProcessLevel();
	}
	bool doVetoProcessLevel(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoProcessLevel(a0);
	}
	bool canSetLowEnergySigma(int a0, int a1) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canSetLowEnergySigma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canSetLowEnergySigma(a0, a1);
	}
	double doSetLowEnergySigma(int a0, int a1, double a2, double a3, double a4) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doSetLowEnergySigma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return UserHooks::doSetLowEnergySigma(a0, a1, a2, a3, a4);
	}
	bool canVetoResonanceDecays() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoResonanceDecays");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoResonanceDecays();
	}
	bool doVetoResonanceDecays(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoResonanceDecays");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoResonanceDecays(a0);
	}
	bool canVetoPT() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoPT");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoPT();
	}
	double scaleVetoPT() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "scaleVetoPT");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return UserHooks::scaleVetoPT();
	}
	bool doVetoPT(int a0, const class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoPT");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoPT(a0, a1);
	}
	bool canVetoStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoStep();
	}
	int numberVetoStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "numberVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return UserHooks::numberVetoStep();
	}
	bool doVetoStep(int a0, int a1, int a2, const class Pythia8::Event & a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoStep(a0, a1, a2, a3);
	}
	bool canVetoMPIStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoMPIStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoMPIStep();
	}
	int numberVetoMPIStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "numberVetoMPIStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return UserHooks::numberVetoMPIStep();
	}
	bool doVetoMPIStep(int a0, const class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoMPIStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoMPIStep(a0, a1);
	}
	bool canVetoPartonLevelEarly() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoPartonLevelEarly();
	}
	bool doVetoPartonLevelEarly(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoPartonLevelEarly(a0);
	}
	bool retryPartonLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "retryPartonLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::retryPartonLevel();
	}
	bool canVetoPartonLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoPartonLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoPartonLevel();
	}
	bool doVetoPartonLevel(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoPartonLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoPartonLevel(a0);
	}
	bool canSetResonanceScale() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canSetResonanceScale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canSetResonanceScale();
	}
	double scaleResonance(int a0, const class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "scaleResonance");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return UserHooks::scaleResonance(a0, a1);
	}
	bool canVetoISREmission() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoISREmission");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoISREmission();
	}
	bool doVetoISREmission(int a0, const class Pythia8::Event & a1, int a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoISREmission");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoISREmission(a0, a1, a2);
	}
	bool canVetoFSREmission() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoFSREmission");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoFSREmission();
	}
	bool doVetoFSREmission(int a0, const class Pythia8::Event & a1, int a2, bool a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoFSREmission");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoFSREmission(a0, a1, a2, a3);
	}
	bool canVetoMPIEmission() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoMPIEmission");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoMPIEmission();
	}
	bool doVetoMPIEmission(int a0, const class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoMPIEmission");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoMPIEmission(a0, a1);
	}
	bool canReconnectResonanceSystems() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canReconnectResonanceSystems");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canReconnectResonanceSystems();
	}
	bool doReconnectResonanceSystems(int a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doReconnectResonanceSystems");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doReconnectResonanceSystems(a0, a1);
	}
	bool canChangeFragPar() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canChangeFragPar");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canChangeFragPar();
	}
	void setStringEnds(const class Pythia8::StringEnd * a0, const class Pythia8::StringEnd * a1, class std::vector<int, class std::allocator<int> > a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "setStringEnds");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return UserHooks::setStringEnds(a0, a1, a2);
	}
	bool doChangeFragPar(class Pythia8::StringFlav * a0, class Pythia8::StringZ * a1, class Pythia8::StringPT * a2, int a3, double a4, class std::vector<int, class std::allocator<int> > a5, const class Pythia8::StringEnd * a6) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doChangeFragPar");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doChangeFragPar(a0, a1, a2, a3, a4, a5, a6);
	}
	bool canVetoFragmentation() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoFragmentation");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoFragmentation();
	}
	bool doVetoFragmentation(class Pythia8::Particle a0, const class Pythia8::StringEnd * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoFragmentation");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoFragmentation(a0, a1);
	}
	bool doVetoFragmentation(class Pythia8::Particle a0, class Pythia8::Particle a1, const class Pythia8::StringEnd * a2, const class Pythia8::StringEnd * a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoFragmentation");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoFragmentation(a0, a1, a2, a3);
	}
	bool doVetoFinalTwo(class Pythia8::Particle a0, class Pythia8::Particle a1, const class Pythia8::StringEnd * a2, const class Pythia8::StringEnd * a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoFinalTwo");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoFinalTwo(a0, a1, a2, a3);
	}
	bool canVetoAfterHadronization() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canVetoAfterHadronization");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canVetoAfterHadronization();
	}
	bool doVetoAfterHadronization(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doVetoAfterHadronization");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::doVetoAfterHadronization(a0);
	}
	bool canSetImpactParameter() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "canSetImpactParameter");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::canSetImpactParameter();
	}
	double doSetImpactParameter() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "doSetImpactParameter");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return UserHooks::doSetImpactParameter();
	}
	bool onEndHadronLevel(class Pythia8::HadronLevel & a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "onEndHadronLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return UserHooks::onEndHadronLevel(a0, a1);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return UserHooks::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HeavyIons::InfoGrabber *>(this), "onStat");
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

// Pythia8::Angantyr file:Pythia8/HeavyIons.h line:162
struct PyCallBack_Pythia8_Angantyr : public Pythia8::Angantyr {
	using Pythia8::Angantyr::Angantyr;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Angantyr::init();
	}
	bool next() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "next");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Angantyr::next();
	}
	bool setKinematics(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "setKinematics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Angantyr::setKinematics(a0);
	}
	bool setKinematics(double a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "setKinematics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Angantyr::setKinematics(a0, a1);
	}
	bool setKinematics(double a0, double a1, double a2, double a3, double a4, double a5) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "setKinematics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Angantyr::setKinematics(a0, a1, a2, a3, a4, a5);
	}
	bool setKinematics(class Pythia8::Vec4 a0, class Pythia8::Vec4 a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "setKinematics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Angantyr::setKinematics(a0, a1);
	}
	bool setBeamIDs(int a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "setBeamIDs");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Angantyr::setBeamIDs(a0, a1);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Angantyr::onInitInfoPtr();
	}
	void stat() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "stat");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HeavyIons::stat();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Angantyr *>(this), "onStat");
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

void bind_Pythia8_HIInfo(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::HIUserHooks file:Pythia8/HIInfo.h line:313
		pybind11::class_<Pythia8::HIUserHooks, std::shared_ptr<Pythia8::HIUserHooks>, PyCallBack_Pythia8_HIUserHooks> cl(M("Pythia8"), "HIUserHooks", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HIUserHooks(); }, [](){ return new PyCallBack_Pythia8_HIUserHooks(); } ) );
		cl.def_readwrite("idProjSave", &Pythia8::HIUserHooks::idProjSave);
		cl.def_readwrite("idTargSave", &Pythia8::HIUserHooks::idTargSave);
		cl.def("init", (void (Pythia8::HIUserHooks::*)(int, int)) &Pythia8::HIUserHooks::init, "C++: Pythia8::HIUserHooks::init(int, int) --> void", pybind11::arg("idProjIn"), pybind11::arg("idTargIn"));
		cl.def("hasImpactParameterGenerator", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::hasImpactParameterGenerator, "C++: Pythia8::HIUserHooks::hasImpactParameterGenerator() const --> bool");
		cl.def("impactParameterGenerator", (class std::shared_ptr<class Pythia8::ImpactParameterGenerator> (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::impactParameterGenerator, "C++: Pythia8::HIUserHooks::impactParameterGenerator() const --> class std::shared_ptr<class Pythia8::ImpactParameterGenerator>");
		cl.def("hasProjectileModel", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::hasProjectileModel, "C++: Pythia8::HIUserHooks::hasProjectileModel() const --> bool");
		cl.def("projectileModel", (class std::shared_ptr<class Pythia8::NucleusModel> (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::projectileModel, "C++: Pythia8::HIUserHooks::projectileModel() const --> class std::shared_ptr<class Pythia8::NucleusModel>");
		cl.def("hasTargetModel", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::hasTargetModel, "C++: Pythia8::HIUserHooks::hasTargetModel() const --> bool");
		cl.def("targetModel", (class std::shared_ptr<class Pythia8::NucleusModel> (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::targetModel, "C++: Pythia8::HIUserHooks::targetModel() const --> class std::shared_ptr<class Pythia8::NucleusModel>");
		cl.def("hasSubCollisionModel", (bool (Pythia8::HIUserHooks::*)()) &Pythia8::HIUserHooks::hasSubCollisionModel, "C++: Pythia8::HIUserHooks::hasSubCollisionModel() --> bool");
		cl.def("subCollisionModel", (class std::shared_ptr<class Pythia8::SubCollisionModel> (Pythia8::HIUserHooks::*)()) &Pythia8::HIUserHooks::subCollisionModel, "C++: Pythia8::HIUserHooks::subCollisionModel() --> class std::shared_ptr<class Pythia8::SubCollisionModel>");
		cl.def("hasEventOrdering", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::hasEventOrdering, "C++: Pythia8::HIUserHooks::hasEventOrdering() const --> bool");
		cl.def("eventOrdering", (double (Pythia8::HIUserHooks::*)(const class Pythia8::Event &, const class Pythia8::Info &)) &Pythia8::HIUserHooks::eventOrdering, "C++: Pythia8::HIUserHooks::eventOrdering(const class Pythia8::Event &, const class Pythia8::Info &) --> double", pybind11::arg(""), pybind11::arg(""));
		cl.def("canFixIsoSpin", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::canFixIsoSpin, "C++: Pythia8::HIUserHooks::canFixIsoSpin() const --> bool");
		cl.def("fixIsoSpin", (bool (Pythia8::HIUserHooks::*)(class Pythia8::EventInfo &)) &Pythia8::HIUserHooks::fixIsoSpin, "C++: Pythia8::HIUserHooks::fixIsoSpin(class Pythia8::EventInfo &) --> bool", pybind11::arg(""));
		cl.def("canShiftEvent", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::canShiftEvent, "C++: Pythia8::HIUserHooks::canShiftEvent() const --> bool");
		cl.def("shiftEvent", (class Pythia8::EventInfo & (Pythia8::HIUserHooks::*)(class Pythia8::EventInfo &) const) &Pythia8::HIUserHooks::shiftEvent, "C++: Pythia8::HIUserHooks::shiftEvent(class Pythia8::EventInfo &) const --> class Pythia8::EventInfo &", pybind11::return_value_policy::reference, pybind11::arg("ei"));
		cl.def("canAddNucleonExcitation", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::canAddNucleonExcitation, "C++: Pythia8::HIUserHooks::canAddNucleonExcitation() const --> bool");
		cl.def("addNucleonExcitation", (bool (Pythia8::HIUserHooks::*)(class Pythia8::EventInfo &, class Pythia8::EventInfo &, bool) const) &Pythia8::HIUserHooks::addNucleonExcitation, "C++: Pythia8::HIUserHooks::addNucleonExcitation(class Pythia8::EventInfo &, class Pythia8::EventInfo &, bool) const --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("canForceHadronLevel", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::canForceHadronLevel, "C++: Pythia8::HIUserHooks::canForceHadronLevel() const --> bool");
		cl.def("forceHadronLevel", (bool (Pythia8::HIUserHooks::*)(class Pythia8::Pythia &)) &Pythia8::HIUserHooks::forceHadronLevel, "C++: Pythia8::HIUserHooks::forceHadronLevel(class Pythia8::Pythia &) --> bool", pybind11::arg(""));
		cl.def("canFindRecoilers", (bool (Pythia8::HIUserHooks::*)() const) &Pythia8::HIUserHooks::canFindRecoilers, "C++: Pythia8::HIUserHooks::canFindRecoilers() const --> bool");
		cl.def("findRecoilers", (class std::vector<int, class std::allocator<int> > (Pythia8::HIUserHooks::*)(const class Pythia8::Event &, bool, int, int, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) const) &Pythia8::HIUserHooks::findRecoilers, "C++: Pythia8::HIUserHooks::findRecoilers(const class Pythia8::Event &, bool, int, int, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) const --> class std::vector<int, class std::allocator<int> >", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("assign", (class Pythia8::HIUserHooks & (Pythia8::HIUserHooks::*)(const class Pythia8::HIUserHooks &)) &Pythia8::HIUserHooks::operator=, "C++: Pythia8::HIUserHooks::operator=(const class Pythia8::HIUserHooks &) --> class Pythia8::HIUserHooks &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HeavyIons file:Pythia8/HeavyIons.h line:31
		pybind11::class_<Pythia8::HeavyIons, std::shared_ptr<Pythia8::HeavyIons>, PyCallBack_Pythia8_HeavyIons, Pythia8::PhysicsBase> cl(M("Pythia8"), "HeavyIons", "");
		pybind11::handle cl_type = cl;

		{ // Pythia8::HeavyIons::InfoGrabber file:Pythia8/HeavyIons.h line:147
			auto & enclosing_class = cl;
			pybind11::class_<Pythia8::HeavyIons::InfoGrabber, std::shared_ptr<Pythia8::HeavyIons::InfoGrabber>, PyCallBack_Pythia8_HeavyIons_InfoGrabber, Pythia8::UserHooks> cl(enclosing_class, "InfoGrabber", "");
			pybind11::handle cl_type = cl;

			cl.def( pybind11::init( [](){ return new Pythia8::HeavyIons::InfoGrabber(); }, [](){ return new PyCallBack_Pythia8_HeavyIons_InfoGrabber(); } ) );
			cl.def("getInfo", (class Pythia8::Info * (Pythia8::HeavyIons::InfoGrabber::*)()) &Pythia8::HeavyIons::InfoGrabber::getInfo, "C++: Pythia8::HeavyIons::InfoGrabber::getInfo() --> class Pythia8::Info *", pybind11::return_value_policy::automatic);
			cl.def("assign", (struct Pythia8::HeavyIons::InfoGrabber & (Pythia8::HeavyIons::InfoGrabber::*)(const struct Pythia8::HeavyIons::InfoGrabber &)) &Pythia8::HeavyIons::InfoGrabber::operator=, "C++: Pythia8::HeavyIons::InfoGrabber::operator=(const struct Pythia8::HeavyIons::InfoGrabber &) --> struct Pythia8::HeavyIons::InfoGrabber &", pybind11::return_value_policy::reference, pybind11::arg(""));
		}

		cl.def( pybind11::init<class Pythia8::Pythia &>(), pybind11::arg("mainPythiaIn") );

		cl.def_readwrite("hiInfo", &Pythia8::HeavyIons::hiInfo);
		cl.def_readwrite("idProj", &Pythia8::HeavyIons::idProj);
		cl.def_readwrite("idTarg", &Pythia8::HeavyIons::idTarg);
		cl.def_readwrite("sigTotNN", &Pythia8::HeavyIons::sigTotNN);
		cl.def_readwrite("HIHooksPtr", &Pythia8::HeavyIons::HIHooksPtr);
		cl.def_readwrite("pythia", &Pythia8::HeavyIons::pythia);
		cl.def_readwrite("pythiaNames", &Pythia8::HeavyIons::pythiaNames);
		cl.def_readwrite("info", &Pythia8::HeavyIons::info);
		cl.def("init", (bool (Pythia8::HeavyIons::*)()) &Pythia8::HeavyIons::init, "C++: Pythia8::HeavyIons::init() --> bool");
		cl.def("next", (bool (Pythia8::HeavyIons::*)()) &Pythia8::HeavyIons::next, "C++: Pythia8::HeavyIons::next() --> bool");
		cl.def_static("addSpecialSettings", (void (*)(class Pythia8::Settings &)) &Pythia8::HeavyIons::addSpecialSettings, "C++: Pythia8::HeavyIons::addSpecialSettings(class Pythia8::Settings &) --> void", pybind11::arg("settings"));
		cl.def_static("isHeavyIon", (bool (*)(class Pythia8::Settings &)) &Pythia8::HeavyIons::isHeavyIon, "C++: Pythia8::HeavyIons::isHeavyIon(class Pythia8::Settings &) --> bool", pybind11::arg("settings"));
		cl.def("setHIUserHooksPtr", (bool (Pythia8::HeavyIons::*)(class std::shared_ptr<class Pythia8::HIUserHooks>)) &Pythia8::HeavyIons::setHIUserHooksPtr, "C++: Pythia8::HeavyIons::setHIUserHooksPtr(class std::shared_ptr<class Pythia8::HIUserHooks>) --> bool", pybind11::arg("userHooksPtrIn"));
		cl.def("setKinematics", (bool (Pythia8::HeavyIons::*)(double)) &Pythia8::HeavyIons::setKinematics, "C++: Pythia8::HeavyIons::setKinematics(double) --> bool", pybind11::arg(""));
		cl.def("setKinematics", (bool (Pythia8::HeavyIons::*)(double, double)) &Pythia8::HeavyIons::setKinematics, "C++: Pythia8::HeavyIons::setKinematics(double, double) --> bool", pybind11::arg(""), pybind11::arg(""));
		cl.def("setKinematics", (bool (Pythia8::HeavyIons::*)(double, double, double, double, double, double)) &Pythia8::HeavyIons::setKinematics, "C++: Pythia8::HeavyIons::setKinematics(double, double, double, double, double, double) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("setKinematics", (bool (Pythia8::HeavyIons::*)(class Pythia8::Vec4, class Pythia8::Vec4)) &Pythia8::HeavyIons::setKinematics, "C++: Pythia8::HeavyIons::setKinematics(class Pythia8::Vec4, class Pythia8::Vec4) --> bool", pybind11::arg(""), pybind11::arg(""));
		cl.def("setBeamIDs", [](Pythia8::HeavyIons &o, int const & a0) -> bool { return o.setBeamIDs(a0); }, "", pybind11::arg(""));
		cl.def("setBeamIDs", (bool (Pythia8::HeavyIons::*)(int, int)) &Pythia8::HeavyIons::setBeamIDs, "C++: Pythia8::HeavyIons::setBeamIDs(int, int) --> bool", pybind11::arg(""), pybind11::arg(""));
		cl.def("stat", (void (Pythia8::HeavyIons::*)()) &Pythia8::HeavyIons::stat, "C++: Pythia8::HeavyIons::stat() --> void");
		cl.def("updateInfo", (void (Pythia8::HeavyIons::*)()) &Pythia8::HeavyIons::updateInfo, "C++: Pythia8::HeavyIons::updateInfo() --> void");
		cl.def("clearProcessLevel", (void (Pythia8::HeavyIons::*)(class Pythia8::Pythia &)) &Pythia8::HeavyIons::clearProcessLevel, "C++: Pythia8::HeavyIons::clearProcessLevel(class Pythia8::Pythia &) --> void", pybind11::arg("pyt"));
		cl.def_static("setupSpecials", (void (*)(class Pythia8::Settings &, std::string)) &Pythia8::HeavyIons::setupSpecials, "C++: Pythia8::HeavyIons::setupSpecials(class Pythia8::Settings &, std::string) --> void", pybind11::arg("settings"), pybind11::arg("match"));
		cl.def_static("setupSpecials", (void (*)(class Pythia8::Pythia &, std::string)) &Pythia8::HeavyIons::setupSpecials, "C++: Pythia8::HeavyIons::setupSpecials(class Pythia8::Pythia &, std::string) --> void", pybind11::arg("p"), pybind11::arg("match"));
		cl.def("assign", (class Pythia8::HeavyIons & (Pythia8::HeavyIons::*)(const class Pythia8::HeavyIons &)) &Pythia8::HeavyIons::operator=, "C++: Pythia8::HeavyIons::operator=(const class Pythia8::HeavyIons &) --> class Pythia8::HeavyIons &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::Angantyr file:Pythia8/HeavyIons.h line:162
		pybind11::class_<Pythia8::Angantyr, std::shared_ptr<Pythia8::Angantyr>, PyCallBack_Pythia8_Angantyr, Pythia8::HeavyIons> cl(M("Pythia8"), "Angantyr", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Pythia &>(), pybind11::arg("mainPythiaIn") );


		pybind11::enum_<Pythia8::Angantyr::PythiaObject>(cl, "PythiaObject", pybind11::arithmetic(), "")
			.value("HADRON", Pythia8::Angantyr::PythiaObject::HADRON)
			.value("MBIAS", Pythia8::Angantyr::PythiaObject::MBIAS)
			.value("SASD", Pythia8::Angantyr::PythiaObject::SASD)
			.value("SIGPP", Pythia8::Angantyr::PythiaObject::SIGPP)
			.value("SIGPN", Pythia8::Angantyr::PythiaObject::SIGPN)
			.value("SIGNP", Pythia8::Angantyr::PythiaObject::SIGNP)
			.value("SIGNN", Pythia8::Angantyr::PythiaObject::SIGNN)
			.value("ALL", Pythia8::Angantyr::PythiaObject::ALL)
			.export_values();

		cl.def("init", (bool (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::init, "C++: Pythia8::Angantyr::init() --> bool");
		cl.def("next", (bool (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::next, "C++: Pythia8::Angantyr::next() --> bool");
		cl.def("setUserHooksPtr", (bool (Pythia8::Angantyr::*)(enum Pythia8::Angantyr::PythiaObject, class std::shared_ptr<class Pythia8::UserHooks>)) &Pythia8::Angantyr::setUserHooksPtr, "C++: Pythia8::Angantyr::setUserHooksPtr(enum Pythia8::Angantyr::PythiaObject, class std::shared_ptr<class Pythia8::UserHooks>) --> bool", pybind11::arg("sel"), pybind11::arg("userHooksPtrIn"));
		cl.def("setKinematicsCM", (bool (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::setKinematicsCM, "C++: Pythia8::Angantyr::setKinematicsCM() --> bool");
		cl.def("setKinematics", (bool (Pythia8::Angantyr::*)(double)) &Pythia8::Angantyr::setKinematics, "C++: Pythia8::Angantyr::setKinematics(double) --> bool", pybind11::arg("eCMIn"));
		cl.def("setKinematics", (bool (Pythia8::Angantyr::*)(double, double)) &Pythia8::Angantyr::setKinematics, "C++: Pythia8::Angantyr::setKinematics(double, double) --> bool", pybind11::arg("eAIn"), pybind11::arg("eBIn"));
		cl.def("setKinematics", (bool (Pythia8::Angantyr::*)(double, double, double, double, double, double)) &Pythia8::Angantyr::setKinematics, "C++: Pythia8::Angantyr::setKinematics(double, double, double, double, double, double) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("setKinematics", (bool (Pythia8::Angantyr::*)(class Pythia8::Vec4, class Pythia8::Vec4)) &Pythia8::Angantyr::setKinematics, "C++: Pythia8::Angantyr::setKinematics(class Pythia8::Vec4, class Pythia8::Vec4) --> bool", pybind11::arg(""), pybind11::arg(""));
		cl.def("setKinematics", (bool (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::setKinematics, "C++: Pythia8::Angantyr::setKinematics() --> bool");
		cl.def("setBeamIDs", [](Pythia8::Angantyr &o, int const & a0) -> bool { return o.setBeamIDs(a0); }, "", pybind11::arg("idAIn"));
		cl.def("setBeamIDs", (bool (Pythia8::Angantyr::*)(int, int)) &Pythia8::Angantyr::setBeamIDs, "C++: Pythia8::Angantyr::setBeamIDs(int, int) --> bool", pybind11::arg("idAIn"), pybind11::arg("idBIn"));
		cl.def("unifyFrames", (void (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::unifyFrames, "C++: Pythia8::Angantyr::unifyFrames() --> void");
		cl.def("banner", (void (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::banner, "C++: Pythia8::Angantyr::banner() const --> void");
		cl.def("subCollisions", (const class Pythia8::SubCollisionSet & (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::subCollisions, "C++: Pythia8::Angantyr::subCollisions() const --> const class Pythia8::SubCollisionSet &", pybind11::return_value_policy::reference);
		cl.def("subCollisionModel", (const class Pythia8::SubCollisionModel & (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::subCollisionModel, "C++: Pythia8::Angantyr::subCollisionModel() const --> const class Pythia8::SubCollisionModel &", pybind11::return_value_policy::reference);
		cl.def("subCollPtr", (class Pythia8::SubCollisionModel * (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::subCollPtr, "C++: Pythia8::Angantyr::subCollPtr() --> class Pythia8::SubCollisionModel *", pybind11::return_value_policy::automatic);
		cl.def("impactParameterGenerator", (const class Pythia8::ImpactParameterGenerator (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::impactParameterGenerator, "C++: Pythia8::Angantyr::impactParameterGenerator() const --> const class Pythia8::ImpactParameterGenerator");
		cl.def("projectile", (const class Pythia8::Nucleus & (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::projectile, "C++: Pythia8::Angantyr::projectile() const --> const class Pythia8::Nucleus &", pybind11::return_value_policy::reference);
		cl.def("target", (const class Pythia8::Nucleus & (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::target, "C++: Pythia8::Angantyr::target() const --> const class Pythia8::Nucleus &", pybind11::return_value_policy::reference);
		cl.def("projectileModel", (const class Pythia8::NucleusModel & (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::projectileModel, "C++: Pythia8::Angantyr::projectileModel() const --> const class Pythia8::NucleusModel &", pybind11::return_value_policy::reference);
		cl.def("targetModel", (const class Pythia8::NucleusModel & (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::targetModel, "C++: Pythia8::Angantyr::targetModel() const --> const class Pythia8::NucleusModel &", pybind11::return_value_policy::reference);
		cl.def("sigmaNN", (const class Pythia8::SigmaTotal (Pythia8::Angantyr::*)() const) &Pythia8::Angantyr::sigmaNN, "C++: Pythia8::Angantyr::sigmaNN() const --> const class Pythia8::SigmaTotal");
		cl.def("onInitInfoPtr", (void (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::onInitInfoPtr, "C++: Pythia8::Angantyr::onInitInfoPtr() --> void");
		cl.def("setBeamKinematics", (void (Pythia8::Angantyr::*)(int, int)) &Pythia8::Angantyr::setBeamKinematics, "C++: Pythia8::Angantyr::setBeamKinematics(int, int) --> void", pybind11::arg("idA"), pybind11::arg("idB"));
		cl.def("init", [](Pythia8::Angantyr &o, enum Pythia8::Angantyr::PythiaObject const & a0, class std::basic_string<char> const & a1) -> bool { return o.init(a0, a1); }, "", pybind11::arg("sel"), pybind11::arg("name"));
		cl.def("init", (bool (Pythia8::Angantyr::*)(enum Pythia8::Angantyr::PythiaObject, std::string, int)) &Pythia8::Angantyr::init, "C++: Pythia8::Angantyr::init(enum Pythia8::Angantyr::PythiaObject, std::string, int) --> bool", pybind11::arg("sel"), pybind11::arg("name"), pybind11::arg("n"));
		cl.def("mkEventInfo", [](Pythia8::Angantyr &o, class Pythia8::Pythia & a0, class Pythia8::Info & a1) -> Pythia8::EventInfo { return o.mkEventInfo(a0, a1); }, "", pybind11::arg(""), pybind11::arg(""));
		cl.def("mkEventInfo", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(class Pythia8::Pythia &, class Pythia8::Info &, const class Pythia8::SubCollision *)) &Pythia8::Angantyr::mkEventInfo, "C++: Pythia8::Angantyr::mkEventInfo(class Pythia8::Pythia &, class Pythia8::Info &, const class Pythia8::SubCollision *) --> class Pythia8::EventInfo", pybind11::arg(""), pybind11::arg(""), pybind11::arg("coll"));
		cl.def("getSignal", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getSignal, "C++: Pythia8::Angantyr::getSignal(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getND", (class Pythia8::EventInfo (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::getND, "C++: Pythia8::Angantyr::getND() --> class Pythia8::EventInfo");
		cl.def("getND", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getND, "C++: Pythia8::Angantyr::getND(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getEl", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getEl, "C++: Pythia8::Angantyr::getEl(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getSDP", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getSDP, "C++: Pythia8::Angantyr::getSDP(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getSDT", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getSDT, "C++: Pythia8::Angantyr::getSDT(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getDD", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getDD, "C++: Pythia8::Angantyr::getDD(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getCD", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getCD, "C++: Pythia8::Angantyr::getCD(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getSDabsP", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getSDabsP, "C++: Pythia8::Angantyr::getSDabsP(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getSDabsT", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision &)) &Pythia8::Angantyr::getSDabsT, "C++: Pythia8::Angantyr::getSDabsT(const class Pythia8::SubCollision &) --> class Pythia8::EventInfo", pybind11::arg("coll"));
		cl.def("getMBIAS", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision *, int)) &Pythia8::Angantyr::getMBIAS, "C++: Pythia8::Angantyr::getMBIAS(const class Pythia8::SubCollision *, int) --> class Pythia8::EventInfo", pybind11::arg("coll"), pybind11::arg("procid"));
		cl.def("getSASD", (class Pythia8::EventInfo (Pythia8::Angantyr::*)(const class Pythia8::SubCollision *, int)) &Pythia8::Angantyr::getSASD, "C++: Pythia8::Angantyr::getSASD(const class Pythia8::SubCollision *, int) --> class Pythia8::EventInfo", pybind11::arg("coll"), pybind11::arg("procid"));
		cl.def("addSASD", (void (Pythia8::Angantyr::*)(const class Pythia8::SubCollisionSet &)) &Pythia8::Angantyr::addSASD, "C++: Pythia8::Angantyr::addSASD(const class Pythia8::SubCollisionSet &) --> void", pybind11::arg("subCollsIn"));
		cl.def("addSDsecond", (void (Pythia8::Angantyr::*)(const class Pythia8::SubCollisionSet &)) &Pythia8::Angantyr::addSDsecond, "C++: Pythia8::Angantyr::addSDsecond(const class Pythia8::SubCollisionSet &) --> void", pybind11::arg("subCollsIn"));
		cl.def("addCDsecond", (void (Pythia8::Angantyr::*)(const class Pythia8::SubCollisionSet &)) &Pythia8::Angantyr::addCDsecond, "C++: Pythia8::Angantyr::addCDsecond(const class Pythia8::SubCollisionSet &) --> void", pybind11::arg("subCollsIn"));
		cl.def("addELsecond", (void (Pythia8::Angantyr::*)(const class Pythia8::SubCollisionSet &)) &Pythia8::Angantyr::addELsecond, "C++: Pythia8::Angantyr::addELsecond(const class Pythia8::SubCollisionSet &) --> void", pybind11::arg("subCollsIn"));
		cl.def("resetEvent", (void (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::resetEvent, "C++: Pythia8::Angantyr::resetEvent() --> void");
		cl.def("setupFullCollision", (bool (Pythia8::Angantyr::*)(class Pythia8::EventInfo &, const class Pythia8::SubCollision &, enum Pythia8::Nucleon::Status, enum Pythia8::Nucleon::Status)) &Pythia8::Angantyr::setupFullCollision, "C++: Pythia8::Angantyr::setupFullCollision(class Pythia8::EventInfo &, const class Pythia8::SubCollision &, enum Pythia8::Nucleon::Status, enum Pythia8::Nucleon::Status) --> bool", pybind11::arg("ei"), pybind11::arg("coll"), pybind11::arg("projStatus"), pybind11::arg("targStatus"));
		cl.def("isRemnant", [](Pythia8::Angantyr const &o, const class Pythia8::EventInfo & a0, int const & a1) -> bool { return o.isRemnant(a0, a1); }, "", pybind11::arg("ei"), pybind11::arg("i"));
		cl.def("isRemnant", (bool (Pythia8::Angantyr::*)(const class Pythia8::EventInfo &, int, int) const) &Pythia8::Angantyr::isRemnant, "C++: Pythia8::Angantyr::isRemnant(const class Pythia8::EventInfo &, int, int) const --> bool", pybind11::arg("ei"), pybind11::arg("i"), pybind11::arg("past"));
		cl.def("fixIsoSpin", (bool (Pythia8::Angantyr::*)(class Pythia8::EventInfo &)) &Pythia8::Angantyr::fixIsoSpin, "C++: Pythia8::Angantyr::fixIsoSpin(class Pythia8::EventInfo &) --> bool", pybind11::arg("ei"));
		cl.def("shiftEvent", (class Pythia8::EventInfo & (Pythia8::Angantyr::*)(class Pythia8::EventInfo &)) &Pythia8::Angantyr::shiftEvent, "C++: Pythia8::Angantyr::shiftEvent(class Pythia8::EventInfo &) --> class Pythia8::EventInfo &", pybind11::return_value_policy::reference, pybind11::arg("ei"));
		cl.def_static("getBeam", (int (*)(class Pythia8::Event &, int)) &Pythia8::Angantyr::getBeam, "C++: Pythia8::Angantyr::getBeam(class Pythia8::Event &, int) --> int", pybind11::arg("ev"), pybind11::arg("i"));
		cl.def("nextSASD", (bool (Pythia8::Angantyr::*)(int)) &Pythia8::Angantyr::nextSASD, "C++: Pythia8::Angantyr::nextSASD(int) --> bool", pybind11::arg("proc"));
		cl.def("addNucleonExcitation", [](Pythia8::Angantyr &o, class Pythia8::EventInfo & a0, class Pythia8::EventInfo & a1) -> bool { return o.addNucleonExcitation(a0, a1); }, "", pybind11::arg("orig"), pybind11::arg("add"));
		cl.def("addNucleonExcitation", (bool (Pythia8::Angantyr::*)(class Pythia8::EventInfo &, class Pythia8::EventInfo &, bool)) &Pythia8::Angantyr::addNucleonExcitation, "C++: Pythia8::Angantyr::addNucleonExcitation(class Pythia8::EventInfo &, class Pythia8::EventInfo &, bool) --> bool", pybind11::arg("orig"), pybind11::arg("add"), pybind11::arg("colConnect"));
		cl.def("findRecoilers", (class std::vector<int, class std::allocator<int> > (Pythia8::Angantyr::*)(const class Pythia8::Event &, bool, int, int, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &)) &Pythia8::Angantyr::findRecoilers, "C++: Pythia8::Angantyr::findRecoilers(const class Pythia8::Event &, bool, int, int, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) --> class std::vector<int, class std::allocator<int> >", pybind11::arg("e"), pybind11::arg("tside"), pybind11::arg("beam"), pybind11::arg("end"), pybind11::arg("pdiff"), pybind11::arg("pbeam"));
		cl.def("addSubEvent", (void (Pythia8::Angantyr::*)(class Pythia8::Event &, class Pythia8::Event &)) &Pythia8::Angantyr::addSubEvent, "C++: Pythia8::Angantyr::addSubEvent(class Pythia8::Event &, class Pythia8::Event &) --> void", pybind11::arg("evnt"), pybind11::arg("sub"));
		cl.def_static("addJunctions", (void (*)(class Pythia8::Event &, class Pythia8::Event &, int)) &Pythia8::Angantyr::addJunctions, "C++: Pythia8::Angantyr::addJunctions(class Pythia8::Event &, class Pythia8::Event &, int) --> void", pybind11::arg("evnt"), pybind11::arg("sub"), pybind11::arg("coloff"));
		cl.def("addNucleusRemnants", (bool (Pythia8::Angantyr::*)()) &Pythia8::Angantyr::addNucleusRemnants, "C++: Pythia8::Angantyr::addNucleusRemnants() --> bool");
		cl.def_static("mT2", (double (*)(const class Pythia8::Vec4 &)) &Pythia8::Angantyr::mT2, "C++: Pythia8::Angantyr::mT2(const class Pythia8::Vec4 &) --> double", pybind11::arg("p"));
		cl.def_static("mT", (double (*)(const class Pythia8::Vec4 &)) &Pythia8::Angantyr::mT, "C++: Pythia8::Angantyr::mT(const class Pythia8::Vec4 &) --> double", pybind11::arg("p"));
		cl.def("assign", (class Pythia8::Angantyr & (Pythia8::Angantyr::*)(const class Pythia8::Angantyr &)) &Pythia8::Angantyr::operator=, "C++: Pythia8::Angantyr::operator=(const class Pythia8::Angantyr &) --> class Pythia8::Angantyr &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
