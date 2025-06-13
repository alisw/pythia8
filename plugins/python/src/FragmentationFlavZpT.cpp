#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/PhysicsBase.h>
#include <iterator>
#include <memory>
#include <sstream> // __str__
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

// Pythia8::StringFlav file:Pythia8/FragmentationFlavZpT.h line:76
struct PyCallBack_Pythia8_StringFlav : public Pythia8::StringFlav {
	using Pythia8::StringFlav::StringFlav;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StringFlav::init();
	}
	void init(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StringFlav::init(a0, a1, a2);
	}
	class Pythia8::FlavContainer pick(class Pythia8::FlavContainer & a0, double a1, double a2, bool a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "pick");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::FlavContainer>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::FlavContainer> caster;
				return pybind11::detail::cast_ref<class Pythia8::FlavContainer>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::FlavContainer>(std::move(o));
		}
		return StringFlav::pick(a0, a1, a2, a3);
	}
	class Pythia8::FlavContainer pickGauss(class Pythia8::FlavContainer & a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "pickGauss");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::FlavContainer>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::FlavContainer> caster;
				return pybind11::detail::cast_ref<class Pythia8::FlavContainer>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::FlavContainer>(std::move(o));
		}
		return StringFlav::pickGauss(a0, a1);
	}
	class Pythia8::FlavContainer pickThermal(class Pythia8::FlavContainer & a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "pickThermal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::FlavContainer>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::FlavContainer> caster;
				return pybind11::detail::cast_ref<class Pythia8::FlavContainer>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::FlavContainer>(std::move(o));
		}
		return StringFlav::pickThermal(a0, a1, a2);
	}
	int combine(class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "combine");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::combine(a0, a1);
	}
	int combineId(int a0, int a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "combineId");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::combineId(a0, a1, a2);
	}
	using _binder_ret_0 = struct std::pair<int, int>;
	_binder_ret_0 combineDiquarkJunction(int a0, int a1, int a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "combineDiquarkJunction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return StringFlav::combineDiquarkJunction(a0, a1, a2);
	}
	int combineToLightest(int a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "combineToLightest");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::combineToLightest(a0, a1);
	}
	int idLightestNeutralMeson() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "idLightestNeutralMeson");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::idLightestNeutralMeson();
	}
	int getHadronIDwin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "getHadronIDwin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::getHadronIDwin();
	}
	int combineLastThermal(class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1, double a2, double a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "combineLastThermal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::combineLastThermal(a0, a1, a2, a3);
	}
	int getHadronID(class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1, double a2, double a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "getHadronID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::getHadronID(a0, a1, a2, a3, a4);
	}
	double getHadronMassWin(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "getHadronMassWin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringFlav::getHadronMassWin(a0);
	}
	void initDerived() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "initDerived");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StringFlav::initDerived();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFlav *>(this), "onStat");
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

// Pythia8::StringZ file:Pythia8/FragmentationFlavZpT.h line:268
struct PyCallBack_Pythia8_StringZ : public Pythia8::StringZ {
	using Pythia8::StringZ::StringZ;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringZ::init();
	}
	double zFrag(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "zFrag");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::zFrag(a0, a1, a2);
	}
	double zLund(double a0, double a1, double a2, double a3, double a4, int a5, bool a6, bool a7, bool a8, bool a9) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "zLund");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::zLund(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	double zPeterson(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "zPeterson");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::zPeterson(a0);
	}
	double zLundMax(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "zLundMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::zLundMax(a0, a1, a2);
	}
	double stopMass() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "stopMass");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::stopMass();
	}
	double stopNewFlav() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "stopNewFlav");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::stopNewFlav();
	}
	double stopSmear() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "stopSmear");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::stopSmear();
	}
	double aAreaLund() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "aAreaLund");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::aAreaLund();
	}
	double bAreaLund() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "bAreaLund");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return StringZ::bAreaLund();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringZ *>(this), "onStat");
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

void bind_Pythia8_FragmentationFlavZpT(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::StringFlav file:Pythia8/FragmentationFlavZpT.h line:76
		pybind11::class_<Pythia8::StringFlav, std::shared_ptr<Pythia8::StringFlav>, PyCallBack_Pythia8_StringFlav, Pythia8::PhysicsBase> cl(M("Pythia8"), "StringFlav", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::StringFlav(); }, [](){ return new PyCallBack_Pythia8_StringFlav(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_StringFlav const &o){ return new PyCallBack_Pythia8_StringFlav(o); } ) );
		cl.def( pybind11::init( [](Pythia8::StringFlav const &o){ return new Pythia8::StringFlav(o); } ) );
		cl.def_readwrite("suppressLeadingB", &Pythia8::StringFlav::suppressLeadingB);
		cl.def_readwrite("mT2suppression", &Pythia8::StringFlav::mT2suppression);
		cl.def_readwrite("useWidthPre", &Pythia8::StringFlav::useWidthPre);
		cl.def_readwrite("probQQtoQ", &Pythia8::StringFlav::probQQtoQ);
		cl.def_readwrite("probStoUD", &Pythia8::StringFlav::probStoUD);
		cl.def_readwrite("probSQtoQQ", &Pythia8::StringFlav::probSQtoQQ);
		cl.def_readwrite("probQQ1toQQ0", &Pythia8::StringFlav::probQQ1toQQ0);
		cl.def_readwrite("probQandQQ", &Pythia8::StringFlav::probQandQQ);
		cl.def_readwrite("probQandS", &Pythia8::StringFlav::probQandS);
		cl.def_readwrite("probQandSinQQ", &Pythia8::StringFlav::probQandSinQQ);
		cl.def_readwrite("probQQ1corr", &Pythia8::StringFlav::probQQ1corr);
		cl.def_readwrite("probQQ1corrInv", &Pythia8::StringFlav::probQQ1corrInv);
		cl.def_readwrite("probQQ1norm", &Pythia8::StringFlav::probQQ1norm);
		cl.def_readwrite("etaSup", &Pythia8::StringFlav::etaSup);
		cl.def_readwrite("etaPrimeSup", &Pythia8::StringFlav::etaPrimeSup);
		cl.def_readwrite("decupletSup", &Pythia8::StringFlav::decupletSup);
		cl.def_readwrite("popcornRate", &Pythia8::StringFlav::popcornRate);
		cl.def_readwrite("popcornSpair", &Pythia8::StringFlav::popcornSpair);
		cl.def_readwrite("popcornSmeson", &Pythia8::StringFlav::popcornSmeson);
		cl.def_readwrite("popFrac", &Pythia8::StringFlav::popFrac);
		cl.def_readwrite("lightLeadingBSup", &Pythia8::StringFlav::lightLeadingBSup);
		cl.def_readwrite("heavyLeadingBSup", &Pythia8::StringFlav::heavyLeadingBSup);
		cl.def_readwrite("qqKappa", &Pythia8::StringFlav::qqKappa);
		cl.def_readwrite("probStoUDSav", &Pythia8::StringFlav::probStoUDSav);
		cl.def_readwrite("probQQtoQSav", &Pythia8::StringFlav::probQQtoQSav);
		cl.def_readwrite("probSQtoQQSav", &Pythia8::StringFlav::probSQtoQQSav);
		cl.def_readwrite("probQQ1toQQ0Sav", &Pythia8::StringFlav::probQQ1toQQ0Sav);
		cl.def_readwrite("alphaQQSav", &Pythia8::StringFlav::alphaQQSav);
		cl.def_readwrite("sigmaHad", &Pythia8::StringFlav::sigmaHad);
		cl.def_readwrite("widthPreStrange", &Pythia8::StringFlav::widthPreStrange);
		cl.def_readwrite("widthPreDiquark", &Pythia8::StringFlav::widthPreDiquark);
		cl.def_readwrite("thermalModel", &Pythia8::StringFlav::thermalModel);
		cl.def_readwrite("mesonNonetL1", &Pythia8::StringFlav::mesonNonetL1);
		cl.def_readwrite("temperature", &Pythia8::StringFlav::temperature);
		cl.def_readwrite("tempPreFactor", &Pythia8::StringFlav::tempPreFactor);
		cl.def_readwrite("nNewQuark", &Pythia8::StringFlav::nNewQuark);
		cl.def_readwrite("closePacking", &Pythia8::StringFlav::closePacking);
		cl.def_readwrite("doEnhanceDiquark", &Pythia8::StringFlav::doEnhanceDiquark);
		cl.def_readwrite("enhanceStrange", &Pythia8::StringFlav::enhanceStrange);
		cl.def_readwrite("enhancePT", &Pythia8::StringFlav::enhancePT);
		cl.def_readwrite("enhanceDiquark", &Pythia8::StringFlav::enhanceDiquark);
		cl.def_readwrite("exponentMPI", &Pythia8::StringFlav::exponentMPI);
		cl.def_readwrite("exponentNSP", &Pythia8::StringFlav::exponentNSP);
		cl.def_readwrite("hadronConstIDs", &Pythia8::StringFlav::hadronConstIDs);
		cl.def_readwrite("possibleHadrons", &Pythia8::StringFlav::possibleHadrons);
		cl.def_readwrite("possibleRatePrefacs", &Pythia8::StringFlav::possibleRatePrefacs);
		cl.def_readwrite("possibleHadronsLast", &Pythia8::StringFlav::possibleHadronsLast);
		cl.def_readwrite("possibleRatePrefacsLast", &Pythia8::StringFlav::possibleRatePrefacsLast);
		cl.def_readwrite("hadronIDwin", &Pythia8::StringFlav::hadronIDwin);
		cl.def_readwrite("idNewWin", &Pythia8::StringFlav::idNewWin);
		cl.def_readwrite("hadronMassWin", &Pythia8::StringFlav::hadronMassWin);
		cl.def("init", (void (Pythia8::StringFlav::*)()) &Pythia8::StringFlav::init, "C++: Pythia8::StringFlav::init() --> void");
		cl.def("init", (void (Pythia8::StringFlav::*)(double, double, double)) &Pythia8::StringFlav::init, "C++: Pythia8::StringFlav::init(double, double, double) --> void", pybind11::arg("kappaModifier"), pybind11::arg("strangeJunc"), pybind11::arg("probQQmod"));
		cl.def("pickLightQ", (int (Pythia8::StringFlav::*)()) &Pythia8::StringFlav::pickLightQ, "C++: Pythia8::StringFlav::pickLightQ() --> int");
		cl.def("pick", [](Pythia8::StringFlav &o, class Pythia8::FlavContainer & a0) -> Pythia8::FlavContainer { return o.pick(a0); }, "", pybind11::arg("flavOld"));
		cl.def("pick", [](Pythia8::StringFlav &o, class Pythia8::FlavContainer & a0, double const & a1) -> Pythia8::FlavContainer { return o.pick(a0, a1); }, "", pybind11::arg("flavOld"), pybind11::arg("pT"));
		cl.def("pick", [](Pythia8::StringFlav &o, class Pythia8::FlavContainer & a0, double const & a1, double const & a2) -> Pythia8::FlavContainer { return o.pick(a0, a1, a2); }, "", pybind11::arg("flavOld"), pybind11::arg("pT"), pybind11::arg("kappaModifier"));
		cl.def("pick", (class Pythia8::FlavContainer (Pythia8::StringFlav::*)(class Pythia8::FlavContainer &, double, double, bool)) &Pythia8::StringFlav::pick, "C++: Pythia8::StringFlav::pick(class Pythia8::FlavContainer &, double, double, bool) --> class Pythia8::FlavContainer", pybind11::arg("flavOld"), pybind11::arg("pT"), pybind11::arg("kappaModifier"), pybind11::arg("allowPop"));
		cl.def("pickGauss", [](Pythia8::StringFlav &o, class Pythia8::FlavContainer & a0) -> Pythia8::FlavContainer { return o.pickGauss(a0); }, "", pybind11::arg("flavOld"));
		cl.def("pickGauss", (class Pythia8::FlavContainer (Pythia8::StringFlav::*)(class Pythia8::FlavContainer &, bool)) &Pythia8::StringFlav::pickGauss, "C++: Pythia8::StringFlav::pickGauss(class Pythia8::FlavContainer &, bool) --> class Pythia8::FlavContainer", pybind11::arg("flavOld"), pybind11::arg("allowPop"));
		cl.def("pickThermal", (class Pythia8::FlavContainer (Pythia8::StringFlav::*)(class Pythia8::FlavContainer &, double, double)) &Pythia8::StringFlav::pickThermal, "C++: Pythia8::StringFlav::pickThermal(class Pythia8::FlavContainer &, double, double) --> class Pythia8::FlavContainer", pybind11::arg("flavOld"), pybind11::arg("pT"), pybind11::arg("kappaModifier"));
		cl.def("combine", (int (Pythia8::StringFlav::*)(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &)) &Pythia8::StringFlav::combine, "C++: Pythia8::StringFlav::combine(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &) --> int", pybind11::arg("flav1"), pybind11::arg("flav2"));
		cl.def("combineId", [](Pythia8::StringFlav &o, int const & a0, int const & a1) -> int { return o.combineId(a0, a1); }, "", pybind11::arg("id1"), pybind11::arg("id2"));
		cl.def("combineId", (int (Pythia8::StringFlav::*)(int, int, bool)) &Pythia8::StringFlav::combineId, "C++: Pythia8::StringFlav::combineId(int, int, bool) --> int", pybind11::arg("id1"), pybind11::arg("id2"), pybind11::arg("keepTrying"));
		cl.def("combineDiquarkJunction", (struct std::pair<int, int> (Pythia8::StringFlav::*)(int, int, int)) &Pythia8::StringFlav::combineDiquarkJunction, "C++: Pythia8::StringFlav::combineDiquarkJunction(int, int, int) --> struct std::pair<int, int>", pybind11::arg("id1"), pybind11::arg("id2"), pybind11::arg("id3"));
		cl.def("combineToLightest", (int (Pythia8::StringFlav::*)(int, int)) &Pythia8::StringFlav::combineToLightest, "C++: Pythia8::StringFlav::combineToLightest(int, int) --> int", pybind11::arg("id1"), pybind11::arg("id2"));
		cl.def("idLightestNeutralMeson", (int (Pythia8::StringFlav::*)()) &Pythia8::StringFlav::idLightestNeutralMeson, "C++: Pythia8::StringFlav::idLightestNeutralMeson() --> int");
		cl.def("getHadronIDwin", (int (Pythia8::StringFlav::*)()) &Pythia8::StringFlav::getHadronIDwin, "C++: Pythia8::StringFlav::getHadronIDwin() --> int");
		cl.def("combineLastThermal", (int (Pythia8::StringFlav::*)(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &, double, double)) &Pythia8::StringFlav::combineLastThermal, "C++: Pythia8::StringFlav::combineLastThermal(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &, double, double) --> int", pybind11::arg("flav1"), pybind11::arg("flav2"), pybind11::arg("pT"), pybind11::arg("kappaModifier"));
		cl.def("getHadronID", [](Pythia8::StringFlav &o, class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1) -> int { return o.getHadronID(a0, a1); }, "", pybind11::arg("flav1"), pybind11::arg("flav2"));
		cl.def("getHadronID", [](Pythia8::StringFlav &o, class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1, double const & a2) -> int { return o.getHadronID(a0, a1, a2); }, "", pybind11::arg("flav1"), pybind11::arg("flav2"), pybind11::arg("pT"));
		cl.def("getHadronID", [](Pythia8::StringFlav &o, class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1, double const & a2, double const & a3) -> int { return o.getHadronID(a0, a1, a2, a3); }, "", pybind11::arg("flav1"), pybind11::arg("flav2"), pybind11::arg("pT"), pybind11::arg("kappaModifier"));
		cl.def("getHadronID", (int (Pythia8::StringFlav::*)(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &, double, double, bool)) &Pythia8::StringFlav::getHadronID, "C++: Pythia8::StringFlav::getHadronID(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &, double, double, bool) --> int", pybind11::arg("flav1"), pybind11::arg("flav2"), pybind11::arg("pT"), pybind11::arg("kappaModifier"), pybind11::arg("finalTwo"));
		cl.def("getHadronMassWin", (double (Pythia8::StringFlav::*)(int)) &Pythia8::StringFlav::getHadronMassWin, "C++: Pythia8::StringFlav::getHadronMassWin(int) --> double", pybind11::arg("idHad"));
		cl.def("assignPopQ", (void (Pythia8::StringFlav::*)(class Pythia8::FlavContainer &)) &Pythia8::StringFlav::assignPopQ, "C++: Pythia8::StringFlav::assignPopQ(class Pythia8::FlavContainer &) --> void", pybind11::arg("flav"));
		cl.def("makeDiquark", [](Pythia8::StringFlav &o, int const & a0, int const & a1) -> int { return o.makeDiquark(a0, a1); }, "", pybind11::arg("id1"), pybind11::arg("id2"));
		cl.def("makeDiquark", (int (Pythia8::StringFlav::*)(int, int, int)) &Pythia8::StringFlav::makeDiquark, "C++: Pythia8::StringFlav::makeDiquark(int, int, int) --> int", pybind11::arg("id1"), pybind11::arg("id2"), pybind11::arg("idHad"));
		cl.def("addQuarkDiquark", (void (Pythia8::StringFlav::*)(class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > &, int, int, int)) &Pythia8::StringFlav::addQuarkDiquark, "C++: Pythia8::StringFlav::addQuarkDiquark(class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > &, int, int, int) --> void", pybind11::arg("quarkCombis"), pybind11::arg("qID"), pybind11::arg("diqID"), pybind11::arg("hadronID"));
		cl.def("getMesonSpinCounter", (int (Pythia8::StringFlav::*)(int)) &Pythia8::StringFlav::getMesonSpinCounter, "C++: Pythia8::StringFlav::getMesonSpinCounter(int) --> int", pybind11::arg("hadronID"));
		cl.def("getFlavourSpinRatios", (double (Pythia8::StringFlav::*)(int, int)) &Pythia8::StringFlav::getFlavourSpinRatios, "C++: Pythia8::StringFlav::getFlavourSpinRatios(int, int) --> double", pybind11::arg("i"), pybind11::arg("j"));
		cl.def("initDerived", (void (Pythia8::StringFlav::*)()) &Pythia8::StringFlav::initDerived, "C++: Pythia8::StringFlav::initDerived() --> void");
		cl.def("assign", (class Pythia8::StringFlav & (Pythia8::StringFlav::*)(const class Pythia8::StringFlav &)) &Pythia8::StringFlav::operator=, "C++: Pythia8::StringFlav::operator=(const class Pythia8::StringFlav &) --> class Pythia8::StringFlav &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	// Pythia8::LundFFRaw(double, double, double, double, double) file:Pythia8/FragmentationFlavZpT.h line:258
	M("Pythia8").def("LundFFRaw", (double (*)(double, double, double, double, double)) &Pythia8::LundFFRaw, "C++: Pythia8::LundFFRaw(double, double, double, double, double) --> double", pybind11::arg("z"), pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("mT2"));

	// Pythia8::LundFFAvg(double, double, double, double) file:Pythia8/FragmentationFlavZpT.h line:260
	M("Pythia8").def("LundFFAvg", (double (*)(double, double, double, double)) &Pythia8::LundFFAvg, "C++: Pythia8::LundFFAvg(double, double, double, double) --> double", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("mT2"), pybind11::arg("tol"));

	// Pythia8::LundFFRms(double, double, double, double) file:Pythia8/FragmentationFlavZpT.h line:262
	M("Pythia8").def("LundFFRms", (double (*)(double, double, double, double)) &Pythia8::LundFFRms, "C++: Pythia8::LundFFRms(double, double, double, double) --> double", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("mT2"), pybind11::arg("tol"));

	{ // Pythia8::StringZ file:Pythia8/FragmentationFlavZpT.h line:268
		pybind11::class_<Pythia8::StringZ, std::shared_ptr<Pythia8::StringZ>, PyCallBack_Pythia8_StringZ, Pythia8::PhysicsBase> cl(M("Pythia8"), "StringZ", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::StringZ(); }, [](){ return new PyCallBack_Pythia8_StringZ(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_StringZ const &o){ return new PyCallBack_Pythia8_StringZ(o); } ) );
		cl.def( pybind11::init( [](Pythia8::StringZ const &o){ return new Pythia8::StringZ(o); } ) );
		cl.def_readwrite("useNonStandC", &Pythia8::StringZ::useNonStandC);
		cl.def_readwrite("useNonStandB", &Pythia8::StringZ::useNonStandB);
		cl.def_readwrite("useNonStandH", &Pythia8::StringZ::useNonStandH);
		cl.def_readwrite("usePetersonC", &Pythia8::StringZ::usePetersonC);
		cl.def_readwrite("usePetersonB", &Pythia8::StringZ::usePetersonB);
		cl.def_readwrite("usePetersonH", &Pythia8::StringZ::usePetersonH);
		cl.def_readwrite("useOldAExtra", &Pythia8::StringZ::useOldAExtra);
		cl.def_readwrite("mc2", &Pythia8::StringZ::mc2);
		cl.def_readwrite("mb2", &Pythia8::StringZ::mb2);
		cl.def_readwrite("aLund", &Pythia8::StringZ::aLund);
		cl.def_readwrite("bLund", &Pythia8::StringZ::bLund);
		cl.def_readwrite("aExtraSQuark", &Pythia8::StringZ::aExtraSQuark);
		cl.def_readwrite("aExtraDiquark", &Pythia8::StringZ::aExtraDiquark);
		cl.def_readwrite("rFactC", &Pythia8::StringZ::rFactC);
		cl.def_readwrite("rFactB", &Pythia8::StringZ::rFactB);
		cl.def_readwrite("rFactH", &Pythia8::StringZ::rFactH);
		cl.def_readwrite("aNonC", &Pythia8::StringZ::aNonC);
		cl.def_readwrite("aNonB", &Pythia8::StringZ::aNonB);
		cl.def_readwrite("aNonH", &Pythia8::StringZ::aNonH);
		cl.def_readwrite("bNonC", &Pythia8::StringZ::bNonC);
		cl.def_readwrite("bNonB", &Pythia8::StringZ::bNonB);
		cl.def_readwrite("bNonH", &Pythia8::StringZ::bNonH);
		cl.def_readwrite("epsilonC", &Pythia8::StringZ::epsilonC);
		cl.def_readwrite("epsilonB", &Pythia8::StringZ::epsilonB);
		cl.def_readwrite("epsilonH", &Pythia8::StringZ::epsilonH);
		cl.def_readwrite("stopM", &Pythia8::StringZ::stopM);
		cl.def_readwrite("stopNF", &Pythia8::StringZ::stopNF);
		cl.def_readwrite("stopS", &Pythia8::StringZ::stopS);
		cl.def("init", (bool (Pythia8::StringZ::*)()) &Pythia8::StringZ::init, "C++: Pythia8::StringZ::init() --> bool");
		cl.def("zFrag", [](Pythia8::StringZ &o, int const & a0) -> double { return o.zFrag(a0); }, "", pybind11::arg("idOld"));
		cl.def("zFrag", [](Pythia8::StringZ &o, int const & a0, int const & a1) -> double { return o.zFrag(a0, a1); }, "", pybind11::arg("idOld"), pybind11::arg("idNew"));
		cl.def("zFrag", (double (Pythia8::StringZ::*)(int, int, double)) &Pythia8::StringZ::zFrag, "C++: Pythia8::StringZ::zFrag(int, int, double) --> double", pybind11::arg("idOld"), pybind11::arg("idNew"), pybind11::arg("mT2"));
		cl.def("zLund", [](Pythia8::StringZ &o, double const & a0, double const & a1) -> double { return o.zLund(a0, a1); }, "", pybind11::arg("a"), pybind11::arg("b"));
		cl.def("zLund", [](Pythia8::StringZ &o, double const & a0, double const & a1, double const & a2) -> double { return o.zLund(a0, a1, a2); }, "", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"));
		cl.def("zLund", [](Pythia8::StringZ &o, double const & a0, double const & a1, double const & a2, double const & a3) -> double { return o.zLund(a0, a1, a2, a3); }, "", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("head"));
		cl.def("zLund", [](Pythia8::StringZ &o, double const & a0, double const & a1, double const & a2, double const & a3, double const & a4) -> double { return o.zLund(a0, a1, a2, a3, a4); }, "", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("head"), pybind11::arg("bNow"));
		cl.def("zLund", [](Pythia8::StringZ &o, double const & a0, double const & a1, double const & a2, double const & a3, double const & a4, int const & a5) -> double { return o.zLund(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("head"), pybind11::arg("bNow"), pybind11::arg("idFrag"));
		cl.def("zLund", [](Pythia8::StringZ &o, double const & a0, double const & a1, double const & a2, double const & a3, double const & a4, int const & a5, bool const & a6) -> double { return o.zLund(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("head"), pybind11::arg("bNow"), pybind11::arg("idFrag"), pybind11::arg("isOldSQuark"));
		cl.def("zLund", [](Pythia8::StringZ &o, double const & a0, double const & a1, double const & a2, double const & a3, double const & a4, int const & a5, bool const & a6, bool const & a7) -> double { return o.zLund(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("head"), pybind11::arg("bNow"), pybind11::arg("idFrag"), pybind11::arg("isOldSQuark"), pybind11::arg("isNewSQuark"));
		cl.def("zLund", [](Pythia8::StringZ &o, double const & a0, double const & a1, double const & a2, double const & a3, double const & a4, int const & a5, bool const & a6, bool const & a7, bool const & a8) -> double { return o.zLund(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("head"), pybind11::arg("bNow"), pybind11::arg("idFrag"), pybind11::arg("isOldSQuark"), pybind11::arg("isNewSQuark"), pybind11::arg("isOldDiquark"));
		cl.def("zLund", (double (Pythia8::StringZ::*)(double, double, double, double, double, int, bool, bool, bool, bool)) &Pythia8::StringZ::zLund, "C++: Pythia8::StringZ::zLund(double, double, double, double, double, int, bool, bool, bool, bool) --> double", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"), pybind11::arg("head"), pybind11::arg("bNow"), pybind11::arg("idFrag"), pybind11::arg("isOldSQuark"), pybind11::arg("isNewSQuark"), pybind11::arg("isOldDiquark"), pybind11::arg("isNewDiquark"));
		cl.def("zPeterson", (double (Pythia8::StringZ::*)(double)) &Pythia8::StringZ::zPeterson, "C++: Pythia8::StringZ::zPeterson(double) --> double", pybind11::arg("epsilon"));
		cl.def("zLundMax", [](Pythia8::StringZ &o, double const & a0, double const & a1) -> double { return o.zLundMax(a0, a1); }, "", pybind11::arg("a"), pybind11::arg("b"));
		cl.def("zLundMax", (double (Pythia8::StringZ::*)(double, double, double)) &Pythia8::StringZ::zLundMax, "C++: Pythia8::StringZ::zLundMax(double, double, double) --> double", pybind11::arg("a"), pybind11::arg("b"), pybind11::arg("c"));
		cl.def("stopMass", (double (Pythia8::StringZ::*)()) &Pythia8::StringZ::stopMass, "C++: Pythia8::StringZ::stopMass() --> double");
		cl.def("stopNewFlav", (double (Pythia8::StringZ::*)()) &Pythia8::StringZ::stopNewFlav, "C++: Pythia8::StringZ::stopNewFlav() --> double");
		cl.def("stopSmear", (double (Pythia8::StringZ::*)()) &Pythia8::StringZ::stopSmear, "C++: Pythia8::StringZ::stopSmear() --> double");
		cl.def("aAreaLund", (double (Pythia8::StringZ::*)()) &Pythia8::StringZ::aAreaLund, "C++: Pythia8::StringZ::aAreaLund() --> double");
		cl.def("bAreaLund", (double (Pythia8::StringZ::*)()) &Pythia8::StringZ::bAreaLund, "C++: Pythia8::StringZ::bAreaLund() --> double");
		cl.def("deriveABLund", [](Pythia8::StringZ &o) -> bool { return o.deriveABLund(); }, "");
		cl.def("deriveABLund", [](Pythia8::StringZ &o, bool const & a0) -> bool { return o.deriveABLund(a0); }, "", pybind11::arg("derivaA"));
		cl.def("deriveABLund", [](Pythia8::StringZ &o, bool const & a0, bool const & a1) -> bool { return o.deriveABLund(a0, a1); }, "", pybind11::arg("derivaA"), pybind11::arg("deriveAExtraDiquark"));
		cl.def("deriveABLund", (bool (Pythia8::StringZ::*)(bool, bool, bool)) &Pythia8::StringZ::deriveABLund, "C++: Pythia8::StringZ::deriveABLund(bool, bool, bool) --> bool", pybind11::arg("derivaA"), pybind11::arg("deriveAExtraDiquark"), pybind11::arg("deriveAExtraSQuark"));
		cl.def("deriveBLund", (double (Pythia8::StringZ::*)(double, double, double)) &Pythia8::StringZ::deriveBLund, "C++: Pythia8::StringZ::deriveBLund(double, double, double) --> double", pybind11::arg("avgZ"), pybind11::arg("a"), pybind11::arg("mT2ref"));
		cl.def("assign", (class Pythia8::StringZ & (Pythia8::StringZ::*)(const class Pythia8::StringZ &)) &Pythia8::StringZ::operator=, "C++: Pythia8::StringZ::operator=(const class Pythia8::StringZ &) --> class Pythia8::StringZ &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
