#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/Merging.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonLevel.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Ropewalk.h>
#include <Pythia8/Settings.h>
#include <Pythia8/ShowerModel.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/SimpleSpaceShower.h>
#include <Pythia8/SpaceShower.h>
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

// Pythia8::RopeFragPars file:Pythia8/Ropewalk.h line:304
struct PyCallBack_Pythia8_RopeFragPars : public Pythia8::RopeFragPars {
	using Pythia8::RopeFragPars::RopeFragPars;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopeFragPars *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopeFragPars *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopeFragPars *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopeFragPars *>(this), "onStat");
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

// Pythia8::FlavourRope file:Pythia8/Ropewalk.h line:374
struct PyCallBack_Pythia8_FlavourRope : public Pythia8::FlavourRope {
	using Pythia8::FlavourRope::FlavourRope;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FlavourRope *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return FlavourRope::init();
	}
	bool doChangeFragPar(class Pythia8::StringFlav * a0, class Pythia8::StringZ * a1, class Pythia8::StringPT * a2, double a3, class std::vector<int, class std::allocator<int> > a4, int a5) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FlavourRope *>(this), "doChangeFragPar");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return FlavourRope::doChangeFragPar(a0, a1, a2, a3, a4, a5);
	}
	bool initEvent(class Pythia8::Event & a0, class Pythia8::ColConfig & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FlavourRope *>(this), "initEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return FlavourRope::initEvent(a0, a1);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FlavourRope *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return FlavourRope::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FlavourRope *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FlavourRope *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FlavourRope *>(this), "onStat");
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

// Pythia8::RopewalkShover file:Pythia8/Ropewalk.h line:457
struct PyCallBack_Pythia8_RopewalkShover : public Pythia8::RopewalkShover {
	using Pythia8::RopewalkShover::RopewalkShover;

	bool stringRepulsion(class Pythia8::Event & a0, class Pythia8::ColConfig & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopewalkShover *>(this), "stringRepulsion");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return RopewalkShover::stringRepulsion(a0, a1);
	}
	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopewalkShover *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringRepulsionBase::init();
	}
	bool hadronRepulsion(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopewalkShover *>(this), "hadronRepulsion");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringRepulsionBase::hadronRepulsion(a0);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopewalkShover *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopewalkShover *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopewalkShover *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::RopewalkShover *>(this), "onStat");
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

// Pythia8::ShowerModel file:Pythia8/ShowerModel.h line:28
struct PyCallBack_Pythia8_ShowerModel : public Pythia8::ShowerModel {
	using Pythia8::ShowerModel::ShowerModel;

	bool init(class std::shared_ptr<class Pythia8::Merging> a0, class std::shared_ptr<class Pythia8::MergingHooks> a1, class std::shared_ptr<class Pythia8::PartonVertex> a2, class Pythia8::WeightContainer * a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ShowerModel::init\"");
	}
	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ShowerModel::initAfterBeams\"");
	}
	class std::shared_ptr<class Pythia8::TimeShower> getTimeShower() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "getTimeShower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::TimeShower>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::TimeShower>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::TimeShower>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::TimeShower>>(std::move(o));
		}
		return ShowerModel::getTimeShower();
	}
	class std::shared_ptr<class Pythia8::TimeShower> getTimeDecShower() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "getTimeDecShower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::TimeShower>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::TimeShower>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::TimeShower>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::TimeShower>>(std::move(o));
		}
		return ShowerModel::getTimeDecShower();
	}
	class std::shared_ptr<class Pythia8::SpaceShower> getSpaceShower() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "getSpaceShower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::SpaceShower>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::SpaceShower>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::SpaceShower>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::SpaceShower>>(std::move(o));
		}
		return ShowerModel::getSpaceShower();
	}
	class std::shared_ptr<class Pythia8::Merging> getMerging() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "getMerging");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::Merging>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::Merging>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::Merging>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::Merging>>(std::move(o));
		}
		return ShowerModel::getMerging();
	}
	class std::shared_ptr<class Pythia8::MergingHooks> getMergingHooks() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "getMergingHooks");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::MergingHooks>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::MergingHooks>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::MergingHooks>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::MergingHooks>>(std::move(o));
		}
		return ShowerModel::getMergingHooks();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ShowerModel *>(this), "onStat");
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

// Pythia8::SimpleShowerModel file:Pythia8/ShowerModel.h line:81
struct PyCallBack_Pythia8_SimpleShowerModel : public Pythia8::SimpleShowerModel {
	using Pythia8::SimpleShowerModel::SimpleShowerModel;

	bool init(class std::shared_ptr<class Pythia8::Merging> a0, class std::shared_ptr<class Pythia8::MergingHooks> a1, class std::shared_ptr<class Pythia8::PartonVertex> a2, class Pythia8::WeightContainer * a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleShowerModel::init(a0, a1, a2, a3);
	}
	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleShowerModel::initAfterBeams();
	}
	class std::shared_ptr<class Pythia8::TimeShower> getTimeShower() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "getTimeShower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::TimeShower>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::TimeShower>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::TimeShower>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::TimeShower>>(std::move(o));
		}
		return ShowerModel::getTimeShower();
	}
	class std::shared_ptr<class Pythia8::TimeShower> getTimeDecShower() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "getTimeDecShower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::TimeShower>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::TimeShower>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::TimeShower>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::TimeShower>>(std::move(o));
		}
		return ShowerModel::getTimeDecShower();
	}
	class std::shared_ptr<class Pythia8::SpaceShower> getSpaceShower() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "getSpaceShower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::SpaceShower>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::SpaceShower>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::SpaceShower>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::SpaceShower>>(std::move(o));
		}
		return ShowerModel::getSpaceShower();
	}
	class std::shared_ptr<class Pythia8::Merging> getMerging() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "getMerging");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::Merging>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::Merging>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::Merging>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::Merging>>(std::move(o));
		}
		return ShowerModel::getMerging();
	}
	class std::shared_ptr<class Pythia8::MergingHooks> getMergingHooks() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "getMergingHooks");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class std::shared_ptr<class Pythia8::MergingHooks>>::value) {
				static pybind11::detail::override_caster_t<class std::shared_ptr<class Pythia8::MergingHooks>> caster;
				return pybind11::detail::cast_ref<class std::shared_ptr<class Pythia8::MergingHooks>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::shared_ptr<class Pythia8::MergingHooks>>(std::move(o));
		}
		return ShowerModel::getMergingHooks();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleShowerModel *>(this), "onStat");
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

void bind_Pythia8_Ropewalk(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::RopeFragPars file:Pythia8/Ropewalk.h line:304
		pybind11::class_<Pythia8::RopeFragPars, std::shared_ptr<Pythia8::RopeFragPars>, PyCallBack_Pythia8_RopeFragPars, Pythia8::PhysicsBase> cl(M("Pythia8"), "RopeFragPars", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::RopeFragPars(); }, [](){ return new PyCallBack_Pythia8_RopeFragPars(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_RopeFragPars const &o){ return new PyCallBack_Pythia8_RopeFragPars(o); } ) );
		cl.def( pybind11::init( [](Pythia8::RopeFragPars const &o){ return new Pythia8::RopeFragPars(o); } ) );
		cl.def("init", (bool (Pythia8::RopeFragPars::*)()) &Pythia8::RopeFragPars::init, "C++: Pythia8::RopeFragPars::init() --> bool");
		cl.def("getEffectiveParameters", (class std::map<std::string, double, struct std::less<std::string >, class std::allocator<struct std::pair<const std::string, double> > > (Pythia8::RopeFragPars::*)(double)) &Pythia8::RopeFragPars::getEffectiveParameters, "C++: Pythia8::RopeFragPars::getEffectiveParameters(double) --> class std::map<std::string, double, struct std::less<std::string >, class std::allocator<struct std::pair<const std::string, double> > >", pybind11::arg("h"));
		cl.def("assign", (class Pythia8::RopeFragPars & (Pythia8::RopeFragPars::*)(const class Pythia8::RopeFragPars &)) &Pythia8::RopeFragPars::operator=, "C++: Pythia8::RopeFragPars::operator=(const class Pythia8::RopeFragPars &) --> class Pythia8::RopeFragPars &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::FlavourRope file:Pythia8/Ropewalk.h line:374
		pybind11::class_<Pythia8::FlavourRope, std::shared_ptr<Pythia8::FlavourRope>, PyCallBack_Pythia8_FlavourRope, Pythia8::FragmentationModifierBase> cl(M("Pythia8"), "FlavourRope", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Ropewalk &>(), pybind11::arg("rwIn") );

		cl.def("init", (bool (Pythia8::FlavourRope::*)()) &Pythia8::FlavourRope::init, "C++: Pythia8::FlavourRope::init() --> bool");
		cl.def("doChangeFragPar", (bool (Pythia8::FlavourRope::*)(class Pythia8::StringFlav *, class Pythia8::StringZ *, class Pythia8::StringPT *, double, class std::vector<int, class std::allocator<int> >, int)) &Pythia8::FlavourRope::doChangeFragPar, "C++: Pythia8::FlavourRope::doChangeFragPar(class Pythia8::StringFlav *, class Pythia8::StringZ *, class Pythia8::StringPT *, double, class std::vector<int, class std::allocator<int> >, int) --> bool", pybind11::arg("flavPtr"), pybind11::arg("zPtr"), pybind11::arg("pTPtr"), pybind11::arg("m2Had"), pybind11::arg("iParton"), pybind11::arg("endId"));
		cl.def("setEnhancement", (void (Pythia8::FlavourRope::*)(double)) &Pythia8::FlavourRope::setEnhancement, "C++: Pythia8::FlavourRope::setEnhancement(double) --> void", pybind11::arg("hIn"));
		cl.def("setEventPtr", (void (Pythia8::FlavourRope::*)(class Pythia8::Event &)) &Pythia8::FlavourRope::setEventPtr, "C++: Pythia8::FlavourRope::setEventPtr(class Pythia8::Event &) --> void", pybind11::arg("event"));
		cl.def("initEvent", (bool (Pythia8::FlavourRope::*)(class Pythia8::Event &, class Pythia8::ColConfig &)) &Pythia8::FlavourRope::initEvent, "C++: Pythia8::FlavourRope::initEvent(class Pythia8::Event &, class Pythia8::ColConfig &) --> bool", pybind11::arg("event"), pybind11::arg("colConfig"));
		cl.def("onInitInfoPtr", (void (Pythia8::FlavourRope::*)()) &Pythia8::FlavourRope::onInitInfoPtr, "C++: Pythia8::FlavourRope::onInitInfoPtr() --> void");
		cl.def("assign", (class Pythia8::FlavourRope & (Pythia8::FlavourRope::*)(const class Pythia8::FlavourRope &)) &Pythia8::FlavourRope::operator=, "C++: Pythia8::FlavourRope::operator=(const class Pythia8::FlavourRope &) --> class Pythia8::FlavourRope &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::RopewalkShover file:Pythia8/Ropewalk.h line:457
		pybind11::class_<Pythia8::RopewalkShover, std::shared_ptr<Pythia8::RopewalkShover>, PyCallBack_Pythia8_RopewalkShover, Pythia8::StringRepulsionBase> cl(M("Pythia8"), "RopewalkShover", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Ropewalk &>(), pybind11::arg("rwIn") );

		cl.def("stringRepulsion", (bool (Pythia8::RopewalkShover::*)(class Pythia8::Event &, class Pythia8::ColConfig &)) &Pythia8::RopewalkShover::stringRepulsion, "C++: Pythia8::RopewalkShover::stringRepulsion(class Pythia8::Event &, class Pythia8::ColConfig &) --> bool", pybind11::arg("event"), pybind11::arg("colConfig"));
		cl.def("assign", (class Pythia8::RopewalkShover & (Pythia8::RopewalkShover::*)(const class Pythia8::RopewalkShover &)) &Pythia8::RopewalkShover::operator=, "C++: Pythia8::RopewalkShover::operator=(const class Pythia8::RopewalkShover &) --> class Pythia8::RopewalkShover &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ShowerModel file:Pythia8/ShowerModel.h line:28
		pybind11::class_<Pythia8::ShowerModel, std::shared_ptr<Pythia8::ShowerModel>, PyCallBack_Pythia8_ShowerModel, Pythia8::PhysicsBase> cl(M("Pythia8"), "ShowerModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_ShowerModel(); } ) );
		cl.def_readwrite("timesPtr", &Pythia8::ShowerModel::timesPtr);
		cl.def_readwrite("timesDecPtr", &Pythia8::ShowerModel::timesDecPtr);
		cl.def_readwrite("spacePtr", &Pythia8::ShowerModel::spacePtr);
		cl.def_readwrite("mergingPtr", &Pythia8::ShowerModel::mergingPtr);
		cl.def_readwrite("mergingHooksPtr", &Pythia8::ShowerModel::mergingHooksPtr);
		cl.def("init", (bool (Pythia8::ShowerModel::*)(class std::shared_ptr<class Pythia8::Merging>, class std::shared_ptr<class Pythia8::MergingHooks>, class std::shared_ptr<class Pythia8::PartonVertex>, class Pythia8::WeightContainer *)) &Pythia8::ShowerModel::init, "C++: Pythia8::ShowerModel::init(class std::shared_ptr<class Pythia8::Merging>, class std::shared_ptr<class Pythia8::MergingHooks>, class std::shared_ptr<class Pythia8::PartonVertex>, class Pythia8::WeightContainer *) --> bool", pybind11::arg("mergPtrIn"), pybind11::arg("mergHooksPtrIn"), pybind11::arg("partonVertexPtrIn"), pybind11::arg("weightContainerPtrIn"));
		cl.def("initAfterBeams", (bool (Pythia8::ShowerModel::*)()) &Pythia8::ShowerModel::initAfterBeams, "C++: Pythia8::ShowerModel::initAfterBeams() --> bool");
		cl.def("getTimeShower", (class std::shared_ptr<class Pythia8::TimeShower> (Pythia8::ShowerModel::*)() const) &Pythia8::ShowerModel::getTimeShower, "C++: Pythia8::ShowerModel::getTimeShower() const --> class std::shared_ptr<class Pythia8::TimeShower>");
		cl.def("getTimeDecShower", (class std::shared_ptr<class Pythia8::TimeShower> (Pythia8::ShowerModel::*)() const) &Pythia8::ShowerModel::getTimeDecShower, "C++: Pythia8::ShowerModel::getTimeDecShower() const --> class std::shared_ptr<class Pythia8::TimeShower>");
		cl.def("getSpaceShower", (class std::shared_ptr<class Pythia8::SpaceShower> (Pythia8::ShowerModel::*)() const) &Pythia8::ShowerModel::getSpaceShower, "C++: Pythia8::ShowerModel::getSpaceShower() const --> class std::shared_ptr<class Pythia8::SpaceShower>");
		cl.def("getMerging", (class std::shared_ptr<class Pythia8::Merging> (Pythia8::ShowerModel::*)() const) &Pythia8::ShowerModel::getMerging, "C++: Pythia8::ShowerModel::getMerging() const --> class std::shared_ptr<class Pythia8::Merging>");
		cl.def("getMergingHooks", (class std::shared_ptr<class Pythia8::MergingHooks> (Pythia8::ShowerModel::*)() const) &Pythia8::ShowerModel::getMergingHooks, "C++: Pythia8::ShowerModel::getMergingHooks() const --> class std::shared_ptr<class Pythia8::MergingHooks>");
		cl.def("assign", (class Pythia8::ShowerModel & (Pythia8::ShowerModel::*)(const class Pythia8::ShowerModel &)) &Pythia8::ShowerModel::operator=, "C++: Pythia8::ShowerModel::operator=(const class Pythia8::ShowerModel &) --> class Pythia8::ShowerModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SimpleShowerModel file:Pythia8/ShowerModel.h line:81
		pybind11::class_<Pythia8::SimpleShowerModel, std::shared_ptr<Pythia8::SimpleShowerModel>, PyCallBack_Pythia8_SimpleShowerModel, Pythia8::ShowerModel> cl(M("Pythia8"), "SimpleShowerModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SimpleShowerModel(); }, [](){ return new PyCallBack_Pythia8_SimpleShowerModel(); } ) );
		cl.def("init", (bool (Pythia8::SimpleShowerModel::*)(class std::shared_ptr<class Pythia8::Merging>, class std::shared_ptr<class Pythia8::MergingHooks>, class std::shared_ptr<class Pythia8::PartonVertex>, class Pythia8::WeightContainer *)) &Pythia8::SimpleShowerModel::init, "C++: Pythia8::SimpleShowerModel::init(class std::shared_ptr<class Pythia8::Merging>, class std::shared_ptr<class Pythia8::MergingHooks>, class std::shared_ptr<class Pythia8::PartonVertex>, class Pythia8::WeightContainer *) --> bool", pybind11::arg("mergPtrIn"), pybind11::arg("mergHooksPtrIn"), pybind11::arg("partonVertexPtrIn"), pybind11::arg("weightContainerPtrIn"));
		cl.def("initAfterBeams", (bool (Pythia8::SimpleShowerModel::*)()) &Pythia8::SimpleShowerModel::initAfterBeams, "C++: Pythia8::SimpleShowerModel::initAfterBeams() --> bool");
		cl.def("assign", (class Pythia8::SimpleShowerModel & (Pythia8::SimpleShowerModel::*)(const class Pythia8::SimpleShowerModel &)) &Pythia8::SimpleShowerModel::operator=, "C++: Pythia8::SimpleShowerModel::operator=(const class Pythia8::SimpleShowerModel &) --> class Pythia8::SimpleShowerModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SpaceDipoleEnd file:Pythia8/SimpleSpaceShower.h line:22
		pybind11::class_<Pythia8::SpaceDipoleEnd, std::shared_ptr<Pythia8::SpaceDipoleEnd>> cl(M("Pythia8"), "SpaceDipoleEnd", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SpaceDipoleEnd(); } ), "doc" );
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::SpaceDipoleEnd(a0); } ), "doc" , pybind11::arg("systemIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Pythia8::SpaceDipoleEnd(a0, a1); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3, a4, a5); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colTypeIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3, a4, a5, a6); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colTypeIn"), pybind11::arg("chgTypeIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colTypeIn"), pybind11::arg("chgTypeIn"), pybind11::arg("weakTypeIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7, int const & a8){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colTypeIn"), pybind11::arg("chgTypeIn"), pybind11::arg("weakTypeIn"), pybind11::arg("MEtypeIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7, int const & a8, bool const & a9){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colTypeIn"), pybind11::arg("chgTypeIn"), pybind11::arg("weakTypeIn"), pybind11::arg("MEtypeIn"), pybind11::arg("normalRecoilIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7, int const & a8, bool const & a9, int const & a10){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colTypeIn"), pybind11::arg("chgTypeIn"), pybind11::arg("weakTypeIn"), pybind11::arg("MEtypeIn"), pybind11::arg("normalRecoilIn"), pybind11::arg("weakPolIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, double const & a4, int const & a5, int const & a6, int const & a7, int const & a8, bool const & a9, int const & a10, int const & a11){ return new Pythia8::SpaceDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); } ), "doc" , pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colTypeIn"), pybind11::arg("chgTypeIn"), pybind11::arg("weakTypeIn"), pybind11::arg("MEtypeIn"), pybind11::arg("normalRecoilIn"), pybind11::arg("weakPolIn"), pybind11::arg("iColPartnerIn"));
		cl.def( pybind11::init<int, int, int, int, double, int, int, int, int, bool, int, int, int>(), pybind11::arg("systemIn"), pybind11::arg("sideIn"), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colTypeIn"), pybind11::arg("chgTypeIn"), pybind11::arg("weakTypeIn"), pybind11::arg("MEtypeIn"), pybind11::arg("normalRecoilIn"), pybind11::arg("weakPolIn"), pybind11::arg("iColPartnerIn"), pybind11::arg("idColPartnerIn") );

		cl.def_readwrite("system", &Pythia8::SpaceDipoleEnd::system);
		cl.def_readwrite("side", &Pythia8::SpaceDipoleEnd::side);
		cl.def_readwrite("iRadiator", &Pythia8::SpaceDipoleEnd::iRadiator);
		cl.def_readwrite("iRecoiler", &Pythia8::SpaceDipoleEnd::iRecoiler);
		cl.def_readwrite("pTmax", &Pythia8::SpaceDipoleEnd::pTmax);
		cl.def_readwrite("colType", &Pythia8::SpaceDipoleEnd::colType);
		cl.def_readwrite("chgType", &Pythia8::SpaceDipoleEnd::chgType);
		cl.def_readwrite("weakType", &Pythia8::SpaceDipoleEnd::weakType);
		cl.def_readwrite("MEtype", &Pythia8::SpaceDipoleEnd::MEtype);
		cl.def_readwrite("normalRecoil", &Pythia8::SpaceDipoleEnd::normalRecoil);
		cl.def_readwrite("weakPol", &Pythia8::SpaceDipoleEnd::weakPol);
		cl.def_readwrite("iColPartner", &Pythia8::SpaceDipoleEnd::iColPartner);
		cl.def_readwrite("idColPartner", &Pythia8::SpaceDipoleEnd::idColPartner);
		cl.def_readwrite("nBranch", &Pythia8::SpaceDipoleEnd::nBranch);
		cl.def_readwrite("idDaughter", &Pythia8::SpaceDipoleEnd::idDaughter);
		cl.def_readwrite("idMother", &Pythia8::SpaceDipoleEnd::idMother);
		cl.def_readwrite("idSister", &Pythia8::SpaceDipoleEnd::idSister);
		cl.def_readwrite("iFinPol", &Pythia8::SpaceDipoleEnd::iFinPol);
		cl.def_readwrite("x1", &Pythia8::SpaceDipoleEnd::x1);
		cl.def_readwrite("x2", &Pythia8::SpaceDipoleEnd::x2);
		cl.def_readwrite("m2Dip", &Pythia8::SpaceDipoleEnd::m2Dip);
		cl.def_readwrite("pT2", &Pythia8::SpaceDipoleEnd::pT2);
		cl.def_readwrite("z", &Pythia8::SpaceDipoleEnd::z);
		cl.def_readwrite("xMo", &Pythia8::SpaceDipoleEnd::xMo);
		cl.def_readwrite("Q2", &Pythia8::SpaceDipoleEnd::Q2);
		cl.def_readwrite("mSister", &Pythia8::SpaceDipoleEnd::mSister);
		cl.def_readwrite("m2Sister", &Pythia8::SpaceDipoleEnd::m2Sister);
		cl.def_readwrite("pT2corr", &Pythia8::SpaceDipoleEnd::pT2corr);
		cl.def_readwrite("pT2Old", &Pythia8::SpaceDipoleEnd::pT2Old);
		cl.def_readwrite("zOld", &Pythia8::SpaceDipoleEnd::zOld);
		cl.def_readwrite("asymPol", &Pythia8::SpaceDipoleEnd::asymPol);
		cl.def_readwrite("m2IF", &Pythia8::SpaceDipoleEnd::m2IF);
		cl.def_readwrite("mColPartner", &Pythia8::SpaceDipoleEnd::mColPartner);
		cl.def_readwrite("pAccept", &Pythia8::SpaceDipoleEnd::pAccept);
		cl.def("store", (void (Pythia8::SpaceDipoleEnd::*)(int, int, int, double, double, double, double, double, double, double, double, double, double, int, double, double)) &Pythia8::SpaceDipoleEnd::store, "C++: Pythia8::SpaceDipoleEnd::store(int, int, int, double, double, double, double, double, double, double, double, double, double, int, double, double) --> void", pybind11::arg("idDaughterIn"), pybind11::arg("idMotherIn"), pybind11::arg("idSisterIn"), pybind11::arg("x1In"), pybind11::arg("x2In"), pybind11::arg("m2DipIn"), pybind11::arg("pT2In"), pybind11::arg("zIn"), pybind11::arg("xMoIn"), pybind11::arg("Q2In"), pybind11::arg("mSisterIn"), pybind11::arg("m2SisterIn"), pybind11::arg("pT2corrIn"), pybind11::arg("iColPartnerIn"), pybind11::arg("m2IFIn"), pybind11::arg("mColPartnerIn"));
	}
}
