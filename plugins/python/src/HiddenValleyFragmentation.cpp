#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/HiddenValleyFragmentation.h>
#include <Pythia8/Info.h>
#include <Pythia8/JunctionSplitting.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/NucleonExcitations.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/StringLength.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/Weights.h>
#include <cwchar>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
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

// Pythia8::HVStringFlav file:Pythia8/HiddenValleyFragmentation.h line:21
struct PyCallBack_Pythia8_HVStringFlav : public Pythia8::HVStringFlav {
	using Pythia8::HVStringFlav::HVStringFlav;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HVStringFlav::init();
	}
	class Pythia8::FlavContainer pick(class Pythia8::FlavContainer & a0, double a1, double a2, bool a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "pick");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::FlavContainer>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::FlavContainer> caster;
				return pybind11::detail::cast_ref<class Pythia8::FlavContainer>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::FlavContainer>(std::move(o));
		}
		return HVStringFlav::pick(a0, a1, a2, a3);
	}
	int combine(class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "combine");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return HVStringFlav::combine(a0, a1);
	}
	int idLightestNeutralMeson() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "idLightestNeutralMeson");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return HVStringFlav::idLightestNeutralMeson();
	}
	void init(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "init");
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
	class Pythia8::FlavContainer pickGauss(class Pythia8::FlavContainer & a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "pickGauss");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "pickThermal");
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
	int combineId(int a0, int a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "combineId");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "combineDiquarkJunction");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "combineToLightest");
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
	int getHadronIDwin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "getHadronIDwin");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "combineLastThermal");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "getHadronID");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "getHadronMassWin");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "initDerived");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringFlav *>(this), "onStat");
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

// Pythia8::HVStringPT file:Pythia8/HiddenValleyFragmentation.h line:60
struct PyCallBack_Pythia8_HVStringPT : public Pythia8::HVStringPT {
	using Pythia8::HVStringPT::HVStringPT;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringPT *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HVStringPT::init();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringPT *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringPT *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringPT *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringPT *>(this), "onStat");
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

// Pythia8::HVStringZ file:Pythia8/HiddenValleyFragmentation.h line:88
struct PyCallBack_Pythia8_HVStringZ : public Pythia8::HVStringZ {
	using Pythia8::HVStringZ::HVStringZ;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HVStringZ::init();
	}
	double zFrag(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "zFrag");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HVStringZ::zFrag(a0, a1, a2);
	}
	double stopMass() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "stopMass");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HVStringZ::stopMass();
	}
	double stopNewFlav() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "stopNewFlav");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HVStringZ::stopNewFlav();
	}
	double stopSmear() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "stopSmear");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HVStringZ::stopSmear();
	}
	double zLund(double a0, double a1, double a2, double a3, double a4, int a5, bool a6, bool a7, bool a8, bool a9) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "zLund");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "zPeterson");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "zLundMax");
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
	double aAreaLund() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "aAreaLund");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "bAreaLund");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HVStringZ *>(this), "onStat");
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

// Pythia8::HiddenValleyFragmentation file:Pythia8/HiddenValleyFragmentation.h line:127
struct PyCallBack_Pythia8_HiddenValleyFragmentation : public Pythia8::HiddenValleyFragmentation {
	using Pythia8::HiddenValleyFragmentation::HiddenValleyFragmentation;

	bool init(class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2, class std::shared_ptr<class Pythia8::FragmentationModifierBase> a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HiddenValleyFragmentation *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HiddenValleyFragmentation::init(a0, a1, a2, a3);
	}
	bool fragment(int a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HiddenValleyFragmentation *>(this), "fragment");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HiddenValleyFragmentation::fragment(a0, a1, a2, a3, a4);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HiddenValleyFragmentation *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HiddenValleyFragmentation::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HiddenValleyFragmentation *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HiddenValleyFragmentation *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HiddenValleyFragmentation *>(this), "onStat");
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

// Pythia8::JunctionSplitting file:Pythia8/JunctionSplitting.h line:30
struct PyCallBack_Pythia8_JunctionSplitting : public Pythia8::JunctionSplitting {
	using Pythia8::JunctionSplitting::JunctionSplitting;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JunctionSplitting *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return JunctionSplitting::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JunctionSplitting *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JunctionSplitting *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JunctionSplitting *>(this), "onStat");
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

// Pythia8::NucleonExcitations file:Pythia8/NucleonExcitations.h line:23
struct PyCallBack_Pythia8_NucleonExcitations : public Pythia8::NucleonExcitations {
	using Pythia8::NucleonExcitations::NucleonExcitations;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleonExcitations *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleonExcitations *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleonExcitations *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleonExcitations *>(this), "onStat");
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

void bind_Pythia8_HiddenValleyFragmentation(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::HVStringFlav file:Pythia8/HiddenValleyFragmentation.h line:21
		pybind11::class_<Pythia8::HVStringFlav, std::shared_ptr<Pythia8::HVStringFlav>, PyCallBack_Pythia8_HVStringFlav, Pythia8::StringFlav> cl(M("Pythia8"), "HVStringFlav", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HVStringFlav(); }, [](){ return new PyCallBack_Pythia8_HVStringFlav(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HVStringFlav const &o){ return new PyCallBack_Pythia8_HVStringFlav(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HVStringFlav const &o){ return new Pythia8::HVStringFlav(o); } ) );
		cl.def("init", (void (Pythia8::HVStringFlav::*)()) &Pythia8::HVStringFlav::init, "C++: Pythia8::HVStringFlav::init() --> void");
		cl.def("pick", (class Pythia8::FlavContainer (Pythia8::HVStringFlav::*)(class Pythia8::FlavContainer &, double, double, bool)) &Pythia8::HVStringFlav::pick, "C++: Pythia8::HVStringFlav::pick(class Pythia8::FlavContainer &, double, double, bool) --> class Pythia8::FlavContainer", pybind11::arg("flavOld"), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("combine", (int (Pythia8::HVStringFlav::*)(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &)) &Pythia8::HVStringFlav::combine, "C++: Pythia8::HVStringFlav::combine(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &) --> int", pybind11::arg("flav1"), pybind11::arg("flav2"));
		cl.def("idLightestNeutralMeson", (int (Pythia8::HVStringFlav::*)()) &Pythia8::HVStringFlav::idLightestNeutralMeson, "C++: Pythia8::HVStringFlav::idLightestNeutralMeson() --> int");
		cl.def("assign", (class Pythia8::HVStringFlav & (Pythia8::HVStringFlav::*)(const class Pythia8::HVStringFlav &)) &Pythia8::HVStringFlav::operator=, "C++: Pythia8::HVStringFlav::operator=(const class Pythia8::HVStringFlav &) --> class Pythia8::HVStringFlav &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HVStringPT file:Pythia8/HiddenValleyFragmentation.h line:60
		pybind11::class_<Pythia8::HVStringPT, std::shared_ptr<Pythia8::HVStringPT>, PyCallBack_Pythia8_HVStringPT, Pythia8::StringPT> cl(M("Pythia8"), "HVStringPT", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HVStringPT(); }, [](){ return new PyCallBack_Pythia8_HVStringPT(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HVStringPT const &o){ return new PyCallBack_Pythia8_HVStringPT(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HVStringPT const &o){ return new Pythia8::HVStringPT(o); } ) );
		cl.def("preinit", (void (Pythia8::HVStringPT::*)(int, double)) &Pythia8::HVStringPT::preinit, "C++: Pythia8::HVStringPT::preinit(int, double) --> void", pybind11::arg("setabsigmaIn"), pybind11::arg("rescalebsigmaIn"));
		cl.def("init", (void (Pythia8::HVStringPT::*)()) &Pythia8::HVStringPT::init, "C++: Pythia8::HVStringPT::init() --> void");
		cl.def("assign", (class Pythia8::HVStringPT & (Pythia8::HVStringPT::*)(const class Pythia8::HVStringPT &)) &Pythia8::HVStringPT::operator=, "C++: Pythia8::HVStringPT::operator=(const class Pythia8::HVStringPT &) --> class Pythia8::HVStringPT &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HVStringZ file:Pythia8/HiddenValleyFragmentation.h line:88
		pybind11::class_<Pythia8::HVStringZ, std::shared_ptr<Pythia8::HVStringZ>, PyCallBack_Pythia8_HVStringZ, Pythia8::StringZ> cl(M("Pythia8"), "HVStringZ", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HVStringZ(); }, [](){ return new PyCallBack_Pythia8_HVStringZ(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HVStringZ const &o){ return new PyCallBack_Pythia8_HVStringZ(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HVStringZ const &o){ return new Pythia8::HVStringZ(o); } ) );
		cl.def("preinit", (void (Pythia8::HVStringZ::*)(int, double, double)) &Pythia8::HVStringZ::preinit, "C++: Pythia8::HVStringZ::preinit(int, double, double) --> void", pybind11::arg("setabsigmaIn"), pybind11::arg("rescalebsigmaIn"), pybind11::arg("mVecRatioIn"));
		cl.def("init", (bool (Pythia8::HVStringZ::*)()) &Pythia8::HVStringZ::init, "C++: Pythia8::HVStringZ::init() --> bool");
		cl.def("zFrag", [](Pythia8::HVStringZ &o, int const & a0) -> double { return o.zFrag(a0); }, "", pybind11::arg("idOld"));
		cl.def("zFrag", [](Pythia8::HVStringZ &o, int const & a0, int const & a1) -> double { return o.zFrag(a0, a1); }, "", pybind11::arg("idOld"), pybind11::arg("idNew"));
		cl.def("zFrag", (double (Pythia8::HVStringZ::*)(int, int, double)) &Pythia8::HVStringZ::zFrag, "C++: Pythia8::HVStringZ::zFrag(int, int, double) --> double", pybind11::arg("idOld"), pybind11::arg("idNew"), pybind11::arg("mT2"));
		cl.def("stopMass", (double (Pythia8::HVStringZ::*)()) &Pythia8::HVStringZ::stopMass, "C++: Pythia8::HVStringZ::stopMass() --> double");
		cl.def("stopNewFlav", (double (Pythia8::HVStringZ::*)()) &Pythia8::HVStringZ::stopNewFlav, "C++: Pythia8::HVStringZ::stopNewFlav() --> double");
		cl.def("stopSmear", (double (Pythia8::HVStringZ::*)()) &Pythia8::HVStringZ::stopSmear, "C++: Pythia8::HVStringZ::stopSmear() --> double");
		cl.def("assign", (class Pythia8::HVStringZ & (Pythia8::HVStringZ::*)(const class Pythia8::HVStringZ &)) &Pythia8::HVStringZ::operator=, "C++: Pythia8::HVStringZ::operator=(const class Pythia8::HVStringZ &) --> class Pythia8::HVStringZ &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HiddenValleyFragmentation file:Pythia8/HiddenValleyFragmentation.h line:127
		pybind11::class_<Pythia8::HiddenValleyFragmentation, std::shared_ptr<Pythia8::HiddenValleyFragmentation>, PyCallBack_Pythia8_HiddenValleyFragmentation, Pythia8::FragmentationModel> cl(M("Pythia8"), "HiddenValleyFragmentation", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HiddenValleyFragmentation(); }, [](){ return new PyCallBack_Pythia8_HiddenValleyFragmentation(); } ) );
		cl.def("init", [](Pythia8::HiddenValleyFragmentation &o) -> bool { return o.init(); }, "");
		cl.def("init", [](Pythia8::HiddenValleyFragmentation &o, class Pythia8::StringFlav * a0) -> bool { return o.init(a0); }, "", pybind11::arg("flavSelPtrIn"));
		cl.def("init", [](Pythia8::HiddenValleyFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1) -> bool { return o.init(a0, a1); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"));
		cl.def("init", [](Pythia8::HiddenValleyFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2) -> bool { return o.init(a0, a1, a2); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"));
		cl.def("init", (bool (Pythia8::HiddenValleyFragmentation::*)(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>)) &Pythia8::HiddenValleyFragmentation::init, "C++: Pythia8::HiddenValleyFragmentation::init(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>) --> bool", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"), pybind11::arg("fragModPtrIn"));
		cl.def("fragment", [](Pythia8::HiddenValleyFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2) -> bool { return o.fragment(a0, a1, a2); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"));
		cl.def("fragment", [](Pythia8::HiddenValleyFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool const & a3) -> bool { return o.fragment(a0, a1, a2, a3); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"));
		cl.def("fragment", (bool (Pythia8::HiddenValleyFragmentation::*)(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool)) &Pythia8::HiddenValleyFragmentation::fragment, "C++: Pythia8::HiddenValleyFragmentation::fragment(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool) --> bool", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"), pybind11::arg("systemRecoil"));
		cl.def("onInitInfoPtr", (void (Pythia8::HiddenValleyFragmentation::*)()) &Pythia8::HiddenValleyFragmentation::onInitInfoPtr, "C++: Pythia8::HiddenValleyFragmentation::onInitInfoPtr() --> void");
		cl.def("assign", (class Pythia8::HiddenValleyFragmentation & (Pythia8::HiddenValleyFragmentation::*)(const class Pythia8::HiddenValleyFragmentation &)) &Pythia8::HiddenValleyFragmentation::operator=, "C++: Pythia8::HiddenValleyFragmentation::operator=(const class Pythia8::HiddenValleyFragmentation &) --> class Pythia8::HiddenValleyFragmentation &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::StringLength file:Pythia8/StringLength.h line:23
		pybind11::class_<Pythia8::StringLength, std::shared_ptr<Pythia8::StringLength>> cl(M("Pythia8"), "StringLength", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::StringLength(); } ) );
		cl.def( pybind11::init( [](Pythia8::StringLength const &o){ return new Pythia8::StringLength(o); } ) );
		cl.def("init", (void (Pythia8::StringLength::*)(class Pythia8::Info *, class Pythia8::Settings &)) &Pythia8::StringLength::init, "C++: Pythia8::StringLength::init(class Pythia8::Info *, class Pythia8::Settings &) --> void", pybind11::arg("infoPtrIn"), pybind11::arg("settings"));
		cl.def("getLength", [](Pythia8::StringLength const &o, const class Pythia8::Vec4 & a0, const class Pythia8::Vec4 & a1) -> double { return o.getLength(a0, a1); }, "", pybind11::arg("p"), pybind11::arg("v"));
		cl.def("getLength", (double (Pythia8::StringLength::*)(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, bool) const) &Pythia8::StringLength::getLength, "C++: Pythia8::StringLength::getLength(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, bool) const --> double", pybind11::arg("p"), pybind11::arg("v"), pybind11::arg("isJunc"));
		cl.def("getStringLength", (double (Pythia8::StringLength::*)(class Pythia8::Event &, int, int)) &Pythia8::StringLength::getStringLength, "C++: Pythia8::StringLength::getStringLength(class Pythia8::Event &, int, int) --> double", pybind11::arg("event"), pybind11::arg("i"), pybind11::arg("j"));
		cl.def("getStringLength", (double (Pythia8::StringLength::*)(class Pythia8::Vec4, class Pythia8::Vec4) const) &Pythia8::StringLength::getStringLength, "C++: Pythia8::StringLength::getStringLength(class Pythia8::Vec4, class Pythia8::Vec4) const --> double", pybind11::arg("p1"), pybind11::arg("p2"));
		cl.def("getJuncLength", (double (Pythia8::StringLength::*)(class Pythia8::Event &, int, int, int)) &Pythia8::StringLength::getJuncLength, "C++: Pythia8::StringLength::getJuncLength(class Pythia8::Event &, int, int, int) --> double", pybind11::arg("event"), pybind11::arg("i"), pybind11::arg("j"), pybind11::arg("k"));
		cl.def("getJuncLength", (double (Pythia8::StringLength::*)(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) const) &Pythia8::StringLength::getJuncLength, "C++: Pythia8::StringLength::getJuncLength(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) const --> double", pybind11::arg("p1"), pybind11::arg("p2"), pybind11::arg("p3"));
		cl.def("getJuncLength", (double (Pythia8::StringLength::*)(class Pythia8::Event &, int, int, int, int)) &Pythia8::StringLength::getJuncLength, "C++: Pythia8::StringLength::getJuncLength(class Pythia8::Event &, int, int, int, int) --> double", pybind11::arg("event"), pybind11::arg("i"), pybind11::arg("j"), pybind11::arg("k"), pybind11::arg("l"));
		cl.def("getJuncLength", (double (Pythia8::StringLength::*)(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) const) &Pythia8::StringLength::getJuncLength, "C++: Pythia8::StringLength::getJuncLength(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) const --> double", pybind11::arg("p1"), pybind11::arg("p2"), pybind11::arg("p3"), pybind11::arg("p4"));
		cl.def("assign", (class Pythia8::StringLength & (Pythia8::StringLength::*)(const class Pythia8::StringLength &)) &Pythia8::StringLength::operator=, "C++: Pythia8::StringLength::operator=(const class Pythia8::StringLength &) --> class Pythia8::StringLength &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::JunctionSplitting file:Pythia8/JunctionSplitting.h line:30
		pybind11::class_<Pythia8::JunctionSplitting, std::shared_ptr<Pythia8::JunctionSplitting>, PyCallBack_Pythia8_JunctionSplitting, Pythia8::PhysicsBase> cl(M("Pythia8"), "JunctionSplitting", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::JunctionSplitting(); }, [](){ return new PyCallBack_Pythia8_JunctionSplitting(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_JunctionSplitting const &o){ return new PyCallBack_Pythia8_JunctionSplitting(o); } ) );
		cl.def( pybind11::init( [](Pythia8::JunctionSplitting const &o){ return new Pythia8::JunctionSplitting(o); } ) );
		cl.def("init", (void (Pythia8::JunctionSplitting::*)()) &Pythia8::JunctionSplitting::init, "C++: Pythia8::JunctionSplitting::init() --> void");
		cl.def("checkColours", (bool (Pythia8::JunctionSplitting::*)(class Pythia8::Event &)) &Pythia8::JunctionSplitting::checkColours, "C++: Pythia8::JunctionSplitting::checkColours(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("onInitInfoPtr", (void (Pythia8::JunctionSplitting::*)()) &Pythia8::JunctionSplitting::onInitInfoPtr, "C++: Pythia8::JunctionSplitting::onInitInfoPtr() --> void");
		cl.def("assign", (class Pythia8::JunctionSplitting & (Pythia8::JunctionSplitting::*)(const class Pythia8::JunctionSplitting &)) &Pythia8::JunctionSplitting::operator=, "C++: Pythia8::JunctionSplitting::operator=(const class Pythia8::JunctionSplitting &) --> class Pythia8::JunctionSplitting &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::NucleonExcitations file:Pythia8/NucleonExcitations.h line:23
		pybind11::class_<Pythia8::NucleonExcitations, std::shared_ptr<Pythia8::NucleonExcitations>, PyCallBack_Pythia8_NucleonExcitations, Pythia8::PhysicsBase> cl(M("Pythia8"), "NucleonExcitations", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::NucleonExcitations(); }, [](){ return new PyCallBack_Pythia8_NucleonExcitations(); } ) );
		cl.def("init", (bool (Pythia8::NucleonExcitations::*)(std::string)) &Pythia8::NucleonExcitations::init, "C++: Pythia8::NucleonExcitations::init(std::string) --> bool", pybind11::arg("path"));
		cl.def("init", (bool (Pythia8::NucleonExcitations::*)(class std::basic_istream<char> &)) &Pythia8::NucleonExcitations::init, "C++: Pythia8::NucleonExcitations::init(class std::basic_istream<char> &) --> bool", pybind11::arg("stream"));
		cl.def("check", (bool (Pythia8::NucleonExcitations::*)()) &Pythia8::NucleonExcitations::check, "C++: Pythia8::NucleonExcitations::check() --> bool");
		cl.def("getExcitationMasks", (class std::vector<int, class std::allocator<int> > (Pythia8::NucleonExcitations::*)() const) &Pythia8::NucleonExcitations::getExcitationMasks, "C++: Pythia8::NucleonExcitations::getExcitationMasks() const --> class std::vector<int, class std::allocator<int> >");
		cl.def("getChannels", (class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > (Pythia8::NucleonExcitations::*)() const) &Pythia8::NucleonExcitations::getChannels, "C++: Pythia8::NucleonExcitations::getChannels() const --> class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > >");
		cl.def("sigmaExTotal", (double (Pythia8::NucleonExcitations::*)(double) const) &Pythia8::NucleonExcitations::sigmaExTotal, "C++: Pythia8::NucleonExcitations::sigmaExTotal(double) const --> double", pybind11::arg("eCM"));
		cl.def("sigmaExPartial", (double (Pythia8::NucleonExcitations::*)(double, int, int) const) &Pythia8::NucleonExcitations::sigmaExPartial, "C++: Pythia8::NucleonExcitations::sigmaExPartial(double, int, int) const --> double", pybind11::arg("eCM"), pybind11::arg("maskC"), pybind11::arg("maskD"));
		cl.def("pickExcitation", (bool (Pythia8::NucleonExcitations::*)(int, int, double, int &, double &, int &, double &)) &Pythia8::NucleonExcitations::pickExcitation, "C++: Pythia8::NucleonExcitations::pickExcitation(int, int, double, int &, double &, int &, double &) --> bool", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"), pybind11::arg("idCOut"), pybind11::arg("mCOut"), pybind11::arg("idDOut"), pybind11::arg("mDOut"));
		cl.def("sigmaCalc", (double (Pythia8::NucleonExcitations::*)(double) const) &Pythia8::NucleonExcitations::sigmaCalc, "C++: Pythia8::NucleonExcitations::sigmaCalc(double) const --> double", pybind11::arg("eCM"));
		cl.def("sigmaCalc", (double (Pythia8::NucleonExcitations::*)(double, int, int) const) &Pythia8::NucleonExcitations::sigmaCalc, "C++: Pythia8::NucleonExcitations::sigmaCalc(double, int, int) const --> double", pybind11::arg("eCM"), pybind11::arg("maskC"), pybind11::arg("maskD"));
		cl.def("parameterizeAll", [](Pythia8::NucleonExcitations &o, int const & a0) -> bool { return o.parameterizeAll(a0); }, "", pybind11::arg("precision"));
		cl.def("parameterizeAll", (bool (Pythia8::NucleonExcitations::*)(int, double)) &Pythia8::NucleonExcitations::parameterizeAll, "C++: Pythia8::NucleonExcitations::parameterizeAll(int, double) --> bool", pybind11::arg("precision"), pybind11::arg("threshold"));
		cl.def("save", (bool (Pythia8::NucleonExcitations::*)(std::ostream &) const) &Pythia8::NucleonExcitations::save, "C++: Pythia8::NucleonExcitations::save(std::ostream &) const --> bool", pybind11::arg("stream"));
		cl.def("save", [](Pythia8::NucleonExcitations const &o) -> bool { return o.save(); }, "");
		cl.def("save", (bool (Pythia8::NucleonExcitations::*)(std::string) const) &Pythia8::NucleonExcitations::save, "C++: Pythia8::NucleonExcitations::save(std::string) const --> bool", pybind11::arg("file"));
	}
}
