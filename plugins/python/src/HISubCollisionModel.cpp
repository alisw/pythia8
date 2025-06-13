#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/HIInfo.h>
#include <Pythia8/HINucleusModel.h>
#include <Pythia8/HISubCollisionModel.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/Weights.h>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <set>
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

// Pythia8::BlackSubCollisionModel file:Pythia8/HISubCollisionModel.h line:323
struct PyCallBack_Pythia8_BlackSubCollisionModel : public Pythia8::BlackSubCollisionModel {
	using Pythia8::BlackSubCollisionModel::BlackSubCollisionModel;

	using _binder_ret_0 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_0 minParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "minParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return BlackSubCollisionModel::minParm();
	}
	using _binder_ret_1 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_1 defParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "defParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		return BlackSubCollisionModel::defParm();
	}
	using _binder_ret_2 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_2 maxParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "maxParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		return BlackSubCollisionModel::maxParm();
	}
	struct Pythia8::SubCollisionModel::SigEst getSig() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "getSig");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::SubCollisionModel::SigEst>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::SubCollisionModel::SigEst> caster;
				return pybind11::detail::cast_ref<struct Pythia8::SubCollisionModel::SigEst>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::SubCollisionModel::SigEst>(std::move(o));
		}
		return BlackSubCollisionModel::getSig();
	}
	class Pythia8::SubCollisionSet getCollisions(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "getCollisions");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::SubCollisionSet>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::SubCollisionSet> caster;
				return pybind11::detail::cast_ref<class Pythia8::SubCollisionSet>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::SubCollisionSet>(std::move(o));
		}
		return BlackSubCollisionModel::getCollisions(a0, a1);
	}
	bool init(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SubCollisionModel::init(a0, a1, a2);
	}
};

// Pythia8::NaiveSubCollisionModel file:Pythia8/HISubCollisionModel.h line:359
struct PyCallBack_Pythia8_NaiveSubCollisionModel : public Pythia8::NaiveSubCollisionModel {
	using Pythia8::NaiveSubCollisionModel::NaiveSubCollisionModel;

	using _binder_ret_0 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_0 minParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "minParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return NaiveSubCollisionModel::minParm();
	}
	using _binder_ret_1 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_1 defParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "defParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		return NaiveSubCollisionModel::defParm();
	}
	using _binder_ret_2 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_2 maxParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "maxParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		return NaiveSubCollisionModel::maxParm();
	}
	struct Pythia8::SubCollisionModel::SigEst getSig() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "getSig");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::SubCollisionModel::SigEst>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::SubCollisionModel::SigEst> caster;
				return pybind11::detail::cast_ref<struct Pythia8::SubCollisionModel::SigEst>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::SubCollisionModel::SigEst>(std::move(o));
		}
		return NaiveSubCollisionModel::getSig();
	}
	class Pythia8::SubCollisionSet getCollisions(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "getCollisions");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::SubCollisionSet>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::SubCollisionSet> caster;
				return pybind11::detail::cast_ref<class Pythia8::SubCollisionSet>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::SubCollisionSet>(std::move(o));
		}
		return NaiveSubCollisionModel::getCollisions(a0, a1);
	}
	bool init(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SubCollisionModel::init(a0, a1, a2);
	}
};

// Pythia8::FluctuatingSubCollisionModel file:Pythia8/HISubCollisionModel.h line:399
struct PyCallBack_Pythia8_FluctuatingSubCollisionModel : public Pythia8::FluctuatingSubCollisionModel {
	using Pythia8::FluctuatingSubCollisionModel::FluctuatingSubCollisionModel;

	class Pythia8::SubCollisionSet getCollisions(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FluctuatingSubCollisionModel *>(this), "getCollisions");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::SubCollisionSet>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::SubCollisionSet> caster;
				return pybind11::detail::cast_ref<class Pythia8::SubCollisionSet>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::SubCollisionSet>(std::move(o));
		}
		return FluctuatingSubCollisionModel::getCollisions(a0, a1);
	}
	struct Pythia8::SubCollisionModel::SigEst getSig() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FluctuatingSubCollisionModel *>(this), "getSig");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::SubCollisionModel::SigEst>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::SubCollisionModel::SigEst> caster;
				return pybind11::detail::cast_ref<struct Pythia8::SubCollisionModel::SigEst>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::SubCollisionModel::SigEst>(std::move(o));
		}
		return FluctuatingSubCollisionModel::getSig();
	}
	double pickRadiusProj() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FluctuatingSubCollisionModel *>(this), "pickRadiusProj");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"FluctuatingSubCollisionModel::pickRadiusProj\"");
	}
	double pickRadiusTarg() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FluctuatingSubCollisionModel *>(this), "pickRadiusTarg");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"FluctuatingSubCollisionModel::pickRadiusTarg\"");
	}
	bool init(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FluctuatingSubCollisionModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SubCollisionModel::init(a0, a1, a2);
	}
	using _binder_ret_0 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_0 minParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FluctuatingSubCollisionModel *>(this), "minParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SubCollisionModel::minParm\"");
	}
	using _binder_ret_1 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_1 defParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FluctuatingSubCollisionModel *>(this), "defParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SubCollisionModel::defParm\"");
	}
	using _binder_ret_2 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_2 maxParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FluctuatingSubCollisionModel *>(this), "maxParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SubCollisionModel::maxParm\"");
	}
};

// Pythia8::DoubleStrikmanSubCollisionModel file:Pythia8/HISubCollisionModel.h line:461
struct PyCallBack_Pythia8_DoubleStrikmanSubCollisionModel : public Pythia8::DoubleStrikmanSubCollisionModel {
	using Pythia8::DoubleStrikmanSubCollisionModel::DoubleStrikmanSubCollisionModel;

	using _binder_ret_0 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_0 minParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DoubleStrikmanSubCollisionModel *>(this), "minParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return DoubleStrikmanSubCollisionModel::minParm();
	}
	using _binder_ret_1 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_1 defParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DoubleStrikmanSubCollisionModel *>(this), "defParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		return DoubleStrikmanSubCollisionModel::defParm();
	}
	using _binder_ret_2 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_2 maxParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DoubleStrikmanSubCollisionModel *>(this), "maxParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		return DoubleStrikmanSubCollisionModel::maxParm();
	}
	double pickRadiusProj() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DoubleStrikmanSubCollisionModel *>(this), "pickRadiusProj");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return DoubleStrikmanSubCollisionModel::pickRadiusProj();
	}
	double pickRadiusTarg() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DoubleStrikmanSubCollisionModel *>(this), "pickRadiusTarg");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return DoubleStrikmanSubCollisionModel::pickRadiusTarg();
	}
	class Pythia8::SubCollisionSet getCollisions(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DoubleStrikmanSubCollisionModel *>(this), "getCollisions");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::SubCollisionSet>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::SubCollisionSet> caster;
				return pybind11::detail::cast_ref<class Pythia8::SubCollisionSet>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::SubCollisionSet>(std::move(o));
		}
		return FluctuatingSubCollisionModel::getCollisions(a0, a1);
	}
	struct Pythia8::SubCollisionModel::SigEst getSig() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DoubleStrikmanSubCollisionModel *>(this), "getSig");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::SubCollisionModel::SigEst>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::SubCollisionModel::SigEst> caster;
				return pybind11::detail::cast_ref<struct Pythia8::SubCollisionModel::SigEst>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::SubCollisionModel::SigEst>(std::move(o));
		}
		return FluctuatingSubCollisionModel::getSig();
	}
	bool init(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DoubleStrikmanSubCollisionModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SubCollisionModel::init(a0, a1, a2);
	}
};

// Pythia8::ImpactParameterGenerator file:Pythia8/HISubCollisionModel.h line:511
struct PyCallBack_Pythia8_ImpactParameterGenerator : public Pythia8::ImpactParameterGenerator {
	using Pythia8::ImpactParameterGenerator::ImpactParameterGenerator;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ImpactParameterGenerator *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ImpactParameterGenerator::init();
	}
	class Pythia8::Vec4 generate(double & a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ImpactParameterGenerator *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Vec4>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Vec4> caster;
				return pybind11::detail::cast_ref<class Pythia8::Vec4>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Vec4>(std::move(o));
		}
		return ImpactParameterGenerator::generate(a0);
	}
	double xSecScale() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ImpactParameterGenerator *>(this), "xSecScale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ImpactParameterGenerator::xSecScale();
	}
};

// Pythia8::LogNormalSubCollisionModel file:Pythia8/HISubCollisionModel.h line:578
struct PyCallBack_Pythia8_LogNormalSubCollisionModel : public Pythia8::LogNormalSubCollisionModel {
	using Pythia8::LogNormalSubCollisionModel::LogNormalSubCollisionModel;

	using _binder_ret_0 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_0 minParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LogNormalSubCollisionModel *>(this), "minParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return LogNormalSubCollisionModel::minParm();
	}
	using _binder_ret_1 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_1 defParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LogNormalSubCollisionModel *>(this), "defParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		return LogNormalSubCollisionModel::defParm();
	}
	using _binder_ret_2 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_2 maxParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LogNormalSubCollisionModel *>(this), "maxParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		return LogNormalSubCollisionModel::maxParm();
	}
	double pickRadiusProj() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LogNormalSubCollisionModel *>(this), "pickRadiusProj");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return LogNormalSubCollisionModel::pickRadiusProj();
	}
	double pickRadiusTarg() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LogNormalSubCollisionModel *>(this), "pickRadiusTarg");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return LogNormalSubCollisionModel::pickRadiusTarg();
	}
	class Pythia8::SubCollisionSet getCollisions(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LogNormalSubCollisionModel *>(this), "getCollisions");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::SubCollisionSet>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::SubCollisionSet> caster;
				return pybind11::detail::cast_ref<class Pythia8::SubCollisionSet>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::SubCollisionSet>(std::move(o));
		}
		return FluctuatingSubCollisionModel::getCollisions(a0, a1);
	}
	struct Pythia8::SubCollisionModel::SigEst getSig() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LogNormalSubCollisionModel *>(this), "getSig");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::SubCollisionModel::SigEst>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::SubCollisionModel::SigEst> caster;
				return pybind11::detail::cast_ref<struct Pythia8::SubCollisionModel::SigEst>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::SubCollisionModel::SigEst>(std::move(o));
		}
		return FluctuatingSubCollisionModel::getSig();
	}
	bool init(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LogNormalSubCollisionModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SubCollisionModel::init(a0, a1, a2);
	}
};

void bind_Pythia8_HISubCollisionModel(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::BlackSubCollisionModel file:Pythia8/HISubCollisionModel.h line:323
		pybind11::class_<Pythia8::BlackSubCollisionModel, std::shared_ptr<Pythia8::BlackSubCollisionModel>, PyCallBack_Pythia8_BlackSubCollisionModel, Pythia8::SubCollisionModel> cl(M("Pythia8"), "BlackSubCollisionModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::BlackSubCollisionModel(); }, [](){ return new PyCallBack_Pythia8_BlackSubCollisionModel(); } ) );
		cl.def("minParm", (class std::vector<double, class std::allocator<double> > (Pythia8::BlackSubCollisionModel::*)() const) &Pythia8::BlackSubCollisionModel::minParm, "C++: Pythia8::BlackSubCollisionModel::minParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("defParm", (class std::vector<double, class std::allocator<double> > (Pythia8::BlackSubCollisionModel::*)() const) &Pythia8::BlackSubCollisionModel::defParm, "C++: Pythia8::BlackSubCollisionModel::defParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("maxParm", (class std::vector<double, class std::allocator<double> > (Pythia8::BlackSubCollisionModel::*)() const) &Pythia8::BlackSubCollisionModel::maxParm, "C++: Pythia8::BlackSubCollisionModel::maxParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("getSig", (struct Pythia8::SubCollisionModel::SigEst (Pythia8::BlackSubCollisionModel::*)() const) &Pythia8::BlackSubCollisionModel::getSig, "C++: Pythia8::BlackSubCollisionModel::getSig() const --> struct Pythia8::SubCollisionModel::SigEst");
		cl.def("getCollisions", (class Pythia8::SubCollisionSet (Pythia8::BlackSubCollisionModel::*)(class Pythia8::Nucleus &, class Pythia8::Nucleus &)) &Pythia8::BlackSubCollisionModel::getCollisions, "C++: Pythia8::BlackSubCollisionModel::getCollisions(class Pythia8::Nucleus &, class Pythia8::Nucleus &) --> class Pythia8::SubCollisionSet", pybind11::arg("proj"), pybind11::arg("targ"));
		cl.def("assign", (class Pythia8::BlackSubCollisionModel & (Pythia8::BlackSubCollisionModel::*)(const class Pythia8::BlackSubCollisionModel &)) &Pythia8::BlackSubCollisionModel::operator=, "C++: Pythia8::BlackSubCollisionModel::operator=(const class Pythia8::BlackSubCollisionModel &) --> class Pythia8::BlackSubCollisionModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::NaiveSubCollisionModel file:Pythia8/HISubCollisionModel.h line:359
		pybind11::class_<Pythia8::NaiveSubCollisionModel, std::shared_ptr<Pythia8::NaiveSubCollisionModel>, PyCallBack_Pythia8_NaiveSubCollisionModel, Pythia8::SubCollisionModel> cl(M("Pythia8"), "NaiveSubCollisionModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::NaiveSubCollisionModel(); }, [](){ return new PyCallBack_Pythia8_NaiveSubCollisionModel(); } ) );
		cl.def("minParm", (class std::vector<double, class std::allocator<double> > (Pythia8::NaiveSubCollisionModel::*)() const) &Pythia8::NaiveSubCollisionModel::minParm, "C++: Pythia8::NaiveSubCollisionModel::minParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("defParm", (class std::vector<double, class std::allocator<double> > (Pythia8::NaiveSubCollisionModel::*)() const) &Pythia8::NaiveSubCollisionModel::defParm, "C++: Pythia8::NaiveSubCollisionModel::defParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("maxParm", (class std::vector<double, class std::allocator<double> > (Pythia8::NaiveSubCollisionModel::*)() const) &Pythia8::NaiveSubCollisionModel::maxParm, "C++: Pythia8::NaiveSubCollisionModel::maxParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("getSig", (struct Pythia8::SubCollisionModel::SigEst (Pythia8::NaiveSubCollisionModel::*)() const) &Pythia8::NaiveSubCollisionModel::getSig, "C++: Pythia8::NaiveSubCollisionModel::getSig() const --> struct Pythia8::SubCollisionModel::SigEst");
		cl.def("getCollisions", (class Pythia8::SubCollisionSet (Pythia8::NaiveSubCollisionModel::*)(class Pythia8::Nucleus &, class Pythia8::Nucleus &)) &Pythia8::NaiveSubCollisionModel::getCollisions, "C++: Pythia8::NaiveSubCollisionModel::getCollisions(class Pythia8::Nucleus &, class Pythia8::Nucleus &) --> class Pythia8::SubCollisionSet", pybind11::arg("proj"), pybind11::arg("targ"));
		cl.def("assign", (class Pythia8::NaiveSubCollisionModel & (Pythia8::NaiveSubCollisionModel::*)(const class Pythia8::NaiveSubCollisionModel &)) &Pythia8::NaiveSubCollisionModel::operator=, "C++: Pythia8::NaiveSubCollisionModel::operator=(const class Pythia8::NaiveSubCollisionModel &) --> class Pythia8::NaiveSubCollisionModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::FluctuatingSubCollisionModel file:Pythia8/HISubCollisionModel.h line:399
		pybind11::class_<Pythia8::FluctuatingSubCollisionModel, std::shared_ptr<Pythia8::FluctuatingSubCollisionModel>, PyCallBack_Pythia8_FluctuatingSubCollisionModel, Pythia8::SubCollisionModel> cl(M("Pythia8"), "FluctuatingSubCollisionModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int, int>(), pybind11::arg("nParmIn"), pybind11::arg("modein") );

		cl.def(pybind11::init<PyCallBack_Pythia8_FluctuatingSubCollisionModel const &>());
		cl.def_readwrite("opacityMode", &Pythia8::FluctuatingSubCollisionModel::opacityMode);
		cl.def("getCollisions", (class Pythia8::SubCollisionSet (Pythia8::FluctuatingSubCollisionModel::*)(class Pythia8::Nucleus &, class Pythia8::Nucleus &)) &Pythia8::FluctuatingSubCollisionModel::getCollisions, "C++: Pythia8::FluctuatingSubCollisionModel::getCollisions(class Pythia8::Nucleus &, class Pythia8::Nucleus &) --> class Pythia8::SubCollisionSet", pybind11::arg("proj"), pybind11::arg("targ"));
		cl.def("getSig", (struct Pythia8::SubCollisionModel::SigEst (Pythia8::FluctuatingSubCollisionModel::*)() const) &Pythia8::FluctuatingSubCollisionModel::getSig, "C++: Pythia8::FluctuatingSubCollisionModel::getSig() const --> struct Pythia8::SubCollisionModel::SigEst");
		cl.def("pickRadiusProj", (double (Pythia8::FluctuatingSubCollisionModel::*)() const) &Pythia8::FluctuatingSubCollisionModel::pickRadiusProj, "C++: Pythia8::FluctuatingSubCollisionModel::pickRadiusProj() const --> double");
		cl.def("pickRadiusTarg", (double (Pythia8::FluctuatingSubCollisionModel::*)() const) &Pythia8::FluctuatingSubCollisionModel::pickRadiusTarg, "C++: Pythia8::FluctuatingSubCollisionModel::pickRadiusTarg() const --> double");
	}
	{ // Pythia8::DoubleStrikmanSubCollisionModel file:Pythia8/HISubCollisionModel.h line:461
		pybind11::class_<Pythia8::DoubleStrikmanSubCollisionModel, std::shared_ptr<Pythia8::DoubleStrikmanSubCollisionModel>, PyCallBack_Pythia8_DoubleStrikmanSubCollisionModel, Pythia8::FluctuatingSubCollisionModel> cl(M("Pythia8"), "DoubleStrikmanSubCollisionModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::DoubleStrikmanSubCollisionModel(); }, [](){ return new PyCallBack_Pythia8_DoubleStrikmanSubCollisionModel(); } ), "doc");
		cl.def( pybind11::init<int>(), pybind11::arg("modeIn") );

		cl.def("minParm", (class std::vector<double, class std::allocator<double> > (Pythia8::DoubleStrikmanSubCollisionModel::*)() const) &Pythia8::DoubleStrikmanSubCollisionModel::minParm, "C++: Pythia8::DoubleStrikmanSubCollisionModel::minParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("defParm", (class std::vector<double, class std::allocator<double> > (Pythia8::DoubleStrikmanSubCollisionModel::*)() const) &Pythia8::DoubleStrikmanSubCollisionModel::defParm, "C++: Pythia8::DoubleStrikmanSubCollisionModel::defParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("maxParm", (class std::vector<double, class std::allocator<double> > (Pythia8::DoubleStrikmanSubCollisionModel::*)() const) &Pythia8::DoubleStrikmanSubCollisionModel::maxParm, "C++: Pythia8::DoubleStrikmanSubCollisionModel::maxParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("pickRadiusProj", (double (Pythia8::DoubleStrikmanSubCollisionModel::*)() const) &Pythia8::DoubleStrikmanSubCollisionModel::pickRadiusProj, "C++: Pythia8::DoubleStrikmanSubCollisionModel::pickRadiusProj() const --> double");
		cl.def("pickRadiusTarg", (double (Pythia8::DoubleStrikmanSubCollisionModel::*)() const) &Pythia8::DoubleStrikmanSubCollisionModel::pickRadiusTarg, "C++: Pythia8::DoubleStrikmanSubCollisionModel::pickRadiusTarg() const --> double");
	}
	{ // Pythia8::ImpactParameterGenerator file:Pythia8/HISubCollisionModel.h line:511
		pybind11::class_<Pythia8::ImpactParameterGenerator, std::shared_ptr<Pythia8::ImpactParameterGenerator>, PyCallBack_Pythia8_ImpactParameterGenerator> cl(M("Pythia8"), "ImpactParameterGenerator", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ImpactParameterGenerator(); }, [](){ return new PyCallBack_Pythia8_ImpactParameterGenerator(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_ImpactParameterGenerator const &o){ return new PyCallBack_Pythia8_ImpactParameterGenerator(o); } ) );
		cl.def( pybind11::init( [](Pythia8::ImpactParameterGenerator const &o){ return new Pythia8::ImpactParameterGenerator(o); } ) );
		cl.def("init", (bool (Pythia8::ImpactParameterGenerator::*)()) &Pythia8::ImpactParameterGenerator::init, "C++: Pythia8::ImpactParameterGenerator::init() --> bool");
		cl.def("initPtr", (void (Pythia8::ImpactParameterGenerator::*)(class Pythia8::Info &, class Pythia8::SubCollisionModel &, class Pythia8::NucleusModel &, class Pythia8::NucleusModel &)) &Pythia8::ImpactParameterGenerator::initPtr, "C++: Pythia8::ImpactParameterGenerator::initPtr(class Pythia8::Info &, class Pythia8::SubCollisionModel &, class Pythia8::NucleusModel &, class Pythia8::NucleusModel &) --> void", pybind11::arg("infoIn"), pybind11::arg("collIn"), pybind11::arg("projIn"), pybind11::arg("targIn"));
		cl.def("generate", (class Pythia8::Vec4 (Pythia8::ImpactParameterGenerator::*)(double &) const) &Pythia8::ImpactParameterGenerator::generate, "C++: Pythia8::ImpactParameterGenerator::generate(double &) const --> class Pythia8::Vec4", pybind11::arg("weight"));
		cl.def("xSecScale", (double (Pythia8::ImpactParameterGenerator::*)() const) &Pythia8::ImpactParameterGenerator::xSecScale, "C++: Pythia8::ImpactParameterGenerator::xSecScale() const --> double");
		cl.def("width", (void (Pythia8::ImpactParameterGenerator::*)(double)) &Pythia8::ImpactParameterGenerator::width, "C++: Pythia8::ImpactParameterGenerator::width(double) --> void", pybind11::arg("widthIn"));
		cl.def("width", (double (Pythia8::ImpactParameterGenerator::*)() const) &Pythia8::ImpactParameterGenerator::width, "C++: Pythia8::ImpactParameterGenerator::width() const --> double");
		cl.def("updateWidth", (void (Pythia8::ImpactParameterGenerator::*)()) &Pythia8::ImpactParameterGenerator::updateWidth, "C++: Pythia8::ImpactParameterGenerator::updateWidth() --> void");
		cl.def("assign", (class Pythia8::ImpactParameterGenerator & (Pythia8::ImpactParameterGenerator::*)(const class Pythia8::ImpactParameterGenerator &)) &Pythia8::ImpactParameterGenerator::operator=, "C++: Pythia8::ImpactParameterGenerator::operator=(const class Pythia8::ImpactParameterGenerator &) --> class Pythia8::ImpactParameterGenerator &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LogNormalSubCollisionModel file:Pythia8/HISubCollisionModel.h line:578
		pybind11::class_<Pythia8::LogNormalSubCollisionModel, std::shared_ptr<Pythia8::LogNormalSubCollisionModel>, PyCallBack_Pythia8_LogNormalSubCollisionModel, Pythia8::FluctuatingSubCollisionModel> cl(M("Pythia8"), "LogNormalSubCollisionModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LogNormalSubCollisionModel(); }, [](){ return new PyCallBack_Pythia8_LogNormalSubCollisionModel(); } ), "doc");
		cl.def( pybind11::init<int>(), pybind11::arg("modeIn") );

		cl.def("minParm", (class std::vector<double, class std::allocator<double> > (Pythia8::LogNormalSubCollisionModel::*)() const) &Pythia8::LogNormalSubCollisionModel::minParm, "C++: Pythia8::LogNormalSubCollisionModel::minParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("defParm", (class std::vector<double, class std::allocator<double> > (Pythia8::LogNormalSubCollisionModel::*)() const) &Pythia8::LogNormalSubCollisionModel::defParm, "C++: Pythia8::LogNormalSubCollisionModel::defParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("maxParm", (class std::vector<double, class std::allocator<double> > (Pythia8::LogNormalSubCollisionModel::*)() const) &Pythia8::LogNormalSubCollisionModel::maxParm, "C++: Pythia8::LogNormalSubCollisionModel::maxParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("pickRadiusProj", (double (Pythia8::LogNormalSubCollisionModel::*)() const) &Pythia8::LogNormalSubCollisionModel::pickRadiusProj, "C++: Pythia8::LogNormalSubCollisionModel::pickRadiusProj() const --> double");
		cl.def("pickRadiusTarg", (double (Pythia8::LogNormalSubCollisionModel::*)() const) &Pythia8::LogNormalSubCollisionModel::pickRadiusTarg, "C++: Pythia8::LogNormalSubCollisionModel::pickRadiusTarg() const --> double");
	}
	{ // Pythia8::HIInfo file:Pythia8/HIInfo.h line:27
		pybind11::class_<Pythia8::HIInfo, std::shared_ptr<Pythia8::HIInfo>> cl(M("Pythia8"), "HIInfo", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HIInfo(); } ) );
		cl.def( pybind11::init( [](Pythia8::HIInfo const &o){ return new Pythia8::HIInfo(o); } ) );
		cl.def("b", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::b, "C++: Pythia8::HIInfo::b() const --> double");
		cl.def("phi", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::phi, "C++: Pythia8::HIInfo::phi() const --> double");
		cl.def("T", (double (Pythia8::HIInfo::*)()) &Pythia8::HIInfo::T, "C++: Pythia8::HIInfo::T() --> double");
		cl.def("avNDb", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::avNDb, "C++: Pythia8::HIInfo::avNDb() const --> double");
		cl.def("nAttempts", (long (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nAttempts, "C++: Pythia8::HIInfo::nAttempts() const --> long");
		cl.def("nAccepted", (long (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nAccepted, "C++: Pythia8::HIInfo::nAccepted() const --> long");
		cl.def("nCollTot", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nCollTot, "C++: Pythia8::HIInfo::nCollTot() const --> int");
		cl.def("nCollND", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nCollND, "C++: Pythia8::HIInfo::nCollND() const --> int");
		cl.def("nCollSDP", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nCollSDP, "C++: Pythia8::HIInfo::nCollSDP() const --> int");
		cl.def("nCollSDT", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nCollSDT, "C++: Pythia8::HIInfo::nCollSDT() const --> int");
		cl.def("nCollDD", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nCollDD, "C++: Pythia8::HIInfo::nCollDD() const --> int");
		cl.def("nCollCD", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nCollCD, "C++: Pythia8::HIInfo::nCollCD() const --> int");
		cl.def("nCollEl", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nCollEl, "C++: Pythia8::HIInfo::nCollEl() const --> int");
		cl.def("nPartProj", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nPartProj, "C++: Pythia8::HIInfo::nPartProj() const --> int");
		cl.def("nAbsProj", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nAbsProj, "C++: Pythia8::HIInfo::nAbsProj() const --> int");
		cl.def("nDiffProj", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nDiffProj, "C++: Pythia8::HIInfo::nDiffProj() const --> int");
		cl.def("nElProj", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nElProj, "C++: Pythia8::HIInfo::nElProj() const --> int");
		cl.def("nPartTarg", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nPartTarg, "C++: Pythia8::HIInfo::nPartTarg() const --> int");
		cl.def("nAbsTarg", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nAbsTarg, "C++: Pythia8::HIInfo::nAbsTarg() const --> int");
		cl.def("nDiffTarg", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nDiffTarg, "C++: Pythia8::HIInfo::nDiffTarg() const --> int");
		cl.def("nElTarg", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nElTarg, "C++: Pythia8::HIInfo::nElTarg() const --> int");
		cl.def("weight", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::weight, "C++: Pythia8::HIInfo::weight() const --> double");
		cl.def("weightSum", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::weightSum, "C++: Pythia8::HIInfo::weightSum() const --> double");
		cl.def("nFail", (int (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::nFail, "C++: Pythia8::HIInfo::nFail() const --> int");
		cl.def("failedExcitation", (void (Pythia8::HIInfo::*)(const class Pythia8::SubCollision &)) &Pythia8::HIInfo::failedExcitation, "C++: Pythia8::HIInfo::failedExcitation(const class Pythia8::SubCollision &) --> void", pybind11::arg("subColl"));
		cl.def("glauberReset", (void (Pythia8::HIInfo::*)()) &Pythia8::HIInfo::glauberReset, "C++: Pythia8::HIInfo::glauberReset() --> void");
		cl.def("glauberTot", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberTot, "C++: Pythia8::HIInfo::glauberTot() const --> double");
		cl.def("glauberTotErr", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberTotErr, "C++: Pythia8::HIInfo::glauberTotErr() const --> double");
		cl.def("glauberND", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberND, "C++: Pythia8::HIInfo::glauberND() const --> double");
		cl.def("glauberNDErr", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberNDErr, "C++: Pythia8::HIInfo::glauberNDErr() const --> double");
		cl.def("glauberINEL", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberINEL, "C++: Pythia8::HIInfo::glauberINEL() const --> double");
		cl.def("glauberINELErr", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberINELErr, "C++: Pythia8::HIInfo::glauberINELErr() const --> double");
		cl.def("glauberEL", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberEL, "C++: Pythia8::HIInfo::glauberEL() const --> double");
		cl.def("glauberELErr", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberELErr, "C++: Pythia8::HIInfo::glauberELErr() const --> double");
		cl.def("glauberDiffT", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberDiffT, "C++: Pythia8::HIInfo::glauberDiffT() const --> double");
		cl.def("glauberDiffTErr", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberDiffTErr, "C++: Pythia8::HIInfo::glauberDiffTErr() const --> double");
		cl.def("glauberDiffP", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberDiffP, "C++: Pythia8::HIInfo::glauberDiffP() const --> double");
		cl.def("glauberDiffPErr", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberDiffPErr, "C++: Pythia8::HIInfo::glauberDiffPErr() const --> double");
		cl.def("glauberDDiff", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberDDiff, "C++: Pythia8::HIInfo::glauberDDiff() const --> double");
		cl.def("glauberDDiffErr", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberDDiffErr, "C++: Pythia8::HIInfo::glauberDDiffErr() const --> double");
		cl.def("glauberBSlope", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberBSlope, "C++: Pythia8::HIInfo::glauberBSlope() const --> double");
		cl.def("glauberBSlopeErr", (double (Pythia8::HIInfo::*)() const) &Pythia8::HIInfo::glauberBSlopeErr, "C++: Pythia8::HIInfo::glauberBSlopeErr() const --> double");
		cl.def("subCollisionsPtr", (const class Pythia8::SubCollisionSet * (Pythia8::HIInfo::*)()) &Pythia8::HIInfo::subCollisionsPtr, "C++: Pythia8::HIInfo::subCollisionsPtr() --> const class Pythia8::SubCollisionSet *", pybind11::return_value_policy::automatic);
		cl.def("assign", (class Pythia8::HIInfo & (Pythia8::HIInfo::*)(const class Pythia8::HIInfo &)) &Pythia8::HIInfo::operator=, "C++: Pythia8::HIInfo::operator=(const class Pythia8::HIInfo &) --> class Pythia8::HIInfo &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
