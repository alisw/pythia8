#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
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
#include <unordered_map>
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

// Pythia8::WeightsBase file:Pythia8/Weights.h line:36
struct PyCallBack_Pythia8_WeightsBase : public Pythia8::WeightsBase {
	using Pythia8::WeightsBase::WeightsBase;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::init();
	}
	void init(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::init(a0);
	}
	void clear() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "clear");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::clear();
	}
	void bookVectors(class std::vector<double, class std::allocator<double> > a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "bookVectors");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::bookVectors(a0, a1);
	}
	void collectWeightValues(class std::vector<double, class std::allocator<double> > & a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "collectWeightValues");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::collectWeightValues(a0, a1);
	}
	void collectWeightNames(class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "collectWeightNames");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::collectWeightNames(a0);
	}
	double getWeightsValue(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "getWeightsValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return WeightsBase::getWeightsValue(a0);
	}
	void reweightValueByIndex(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "reweightValueByIndex");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::reweightValueByIndex(a0, a1);
	}
	void reweightValueByName(class std::basic_string<char> a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsBase *>(this), "reweightValueByName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::reweightValueByName(a0, a1);
	}
};

// Pythia8::WeightsShower file:Pythia8/Weights.h line:145
struct PyCallBack_Pythia8_WeightsShower : public Pythia8::WeightsShower {
	using Pythia8::WeightsShower::WeightsShower;

	void init(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsShower::init(a0);
	}
	void initWeightGroups(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "initWeightGroups");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsShower::initWeightGroups(a0);
	}
	int nWeightGroups() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "nWeightGroups");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return WeightsShower::nWeightGroups();
	}
	class std::basic_string<char> getGroupName(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "getGroupName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::basic_string<char>>::value) {
				static pybind11::detail::override_caster_t<class std::basic_string<char>> caster;
				return pybind11::detail::cast_ref<class std::basic_string<char>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::basic_string<char>>(std::move(o));
		}
		return WeightsShower::getGroupName(a0);
	}
	double getGroupWeight(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "getGroupWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return WeightsShower::getGroupWeight(a0);
	}
	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::init();
	}
	void clear() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "clear");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::clear();
	}
	void bookVectors(class std::vector<double, class std::allocator<double> > a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "bookVectors");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::bookVectors(a0, a1);
	}
	void collectWeightValues(class std::vector<double, class std::allocator<double> > & a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "collectWeightValues");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::collectWeightValues(a0, a1);
	}
	void collectWeightNames(class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "collectWeightNames");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::collectWeightNames(a0);
	}
	double getWeightsValue(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "getWeightsValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return WeightsBase::getWeightsValue(a0);
	}
	void reweightValueByIndex(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "reweightValueByIndex");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::reweightValueByIndex(a0, a1);
	}
	void reweightValueByName(class std::basic_string<char> a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsShower *>(this), "reweightValueByName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::reweightValueByName(a0, a1);
	}
};

// Pythia8::WeightsSimpleShower file:Pythia8/Weights.h line:164
struct PyCallBack_Pythia8_WeightsSimpleShower : public Pythia8::WeightsSimpleShower {
	using Pythia8::WeightsSimpleShower::WeightsSimpleShower;

	void init(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsSimpleShower::init(a0);
	}
	void collectWeightNames(class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "collectWeightNames");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsSimpleShower::collectWeightNames(a0);
	}
	void collectWeightValues(class std::vector<double, class std::allocator<double> > & a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "collectWeightValues");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsSimpleShower::collectWeightValues(a0, a1);
	}
	void initWeightGroups(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "initWeightGroups");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsSimpleShower::initWeightGroups(a0);
	}
	class std::basic_string<char> getGroupName(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "getGroupName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class std::basic_string<char>>::value) {
				static pybind11::detail::override_caster_t<class std::basic_string<char>> caster;
				return pybind11::detail::cast_ref<class std::basic_string<char>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class std::basic_string<char>>(std::move(o));
		}
		return WeightsSimpleShower::getGroupName(a0);
	}
	double getGroupWeight(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "getGroupWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return WeightsSimpleShower::getGroupWeight(a0);
	}
	int nWeightGroups() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "nWeightGroups");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return WeightsSimpleShower::nWeightGroups();
	}
	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::init();
	}
	void clear() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "clear");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::clear();
	}
	void bookVectors(class std::vector<double, class std::allocator<double> > a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "bookVectors");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::bookVectors(a0, a1);
	}
	double getWeightsValue(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "getWeightsValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return WeightsBase::getWeightsValue(a0);
	}
	void reweightValueByIndex(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "reweightValueByIndex");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::reweightValueByIndex(a0, a1);
	}
	void reweightValueByName(class std::basic_string<char> a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsSimpleShower *>(this), "reweightValueByName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsBase::reweightValueByName(a0, a1);
	}
};

void bind_Pythia8_LHEF3(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::Reader file:Pythia8/LHEF3.h line:814
		pybind11::class_<Pythia8::Reader, std::shared_ptr<Pythia8::Reader>> cl(M("Pythia8"), "Reader", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<std::string>(), pybind11::arg("filenameIn") );

		cl.def( pybind11::init<class std::basic_istream<char> *>(), pybind11::arg("is") );

		cl.def_readwrite("filename", &Pythia8::Reader::filename);
		cl.def_readwrite("currentLine", &Pythia8::Reader::currentLine);
		cl.def_readwrite("isGood", &Pythia8::Reader::isGood);
		cl.def_readwrite("version", &Pythia8::Reader::version);
		cl.def_readwrite("outsideBlock", &Pythia8::Reader::outsideBlock);
		cl.def_readwrite("headerBlock", &Pythia8::Reader::headerBlock);
		cl.def_readwrite("headerComments", &Pythia8::Reader::headerComments);
		cl.def_readwrite("heprup", &Pythia8::Reader::heprup);
		cl.def_readwrite("initComments", &Pythia8::Reader::initComments);
		cl.def_readwrite("hepeup", &Pythia8::Reader::hepeup);
		cl.def_readwrite("eventComments", &Pythia8::Reader::eventComments);
		cl.def_readwrite("weights_detailed_vec", &Pythia8::Reader::weights_detailed_vec);
		cl.def_readwrite("weightnames_detailed_vec", &Pythia8::Reader::weightnames_detailed_vec);
		cl.def("setup", (bool (Pythia8::Reader::*)(std::string)) &Pythia8::Reader::setup, "C++: Pythia8::Reader::setup(std::string) --> bool", pybind11::arg("filenameIn"));
		cl.def("readEvent", [](Pythia8::Reader &o) -> bool { return o.readEvent(); }, "");
		cl.def("readEvent", (bool (Pythia8::Reader::*)(class Pythia8::HEPEUP *)) &Pythia8::Reader::readEvent, "C++: Pythia8::Reader::readEvent(class Pythia8::HEPEUP *) --> bool", pybind11::arg("peup"));
		cl.def("clearEvent", (void (Pythia8::Reader::*)()) &Pythia8::Reader::clearEvent, "C++: Pythia8::Reader::clearEvent() --> void");
		cl.def("getLine", (bool (Pythia8::Reader::*)()) &Pythia8::Reader::getLine, "C++: Pythia8::Reader::getLine() --> bool");
		cl.def("weights_detailed_vector", (class std::vector<double, class std::allocator<double> > (Pythia8::Reader::*)()) &Pythia8::Reader::weights_detailed_vector, "C++: Pythia8::Reader::weights_detailed_vector() --> class std::vector<double, class std::allocator<double> >");
		cl.def("weightnames_detailed_vector", (class std::vector<std::string, class std::allocator<std::string > > (Pythia8::Reader::*)()) &Pythia8::Reader::weightnames_detailed_vector, "C++: Pythia8::Reader::weightnames_detailed_vector() --> class std::vector<std::string, class std::allocator<std::string > >");
	}
	{ // Pythia8::WeightsBase file:Pythia8/Weights.h line:36
		pybind11::class_<Pythia8::WeightsBase, std::shared_ptr<Pythia8::WeightsBase>, PyCallBack_Pythia8_WeightsBase> cl(M("Pythia8"), "WeightsBase", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](PyCallBack_Pythia8_WeightsBase const &o){ return new PyCallBack_Pythia8_WeightsBase(o); } ) );
		cl.def( pybind11::init( [](Pythia8::WeightsBase const &o){ return new Pythia8::WeightsBase(o); } ) );
		cl.def( pybind11::init( [](){ return new Pythia8::WeightsBase(); }, [](){ return new PyCallBack_Pythia8_WeightsBase(); } ) );
		cl.def_readwrite("weightValues", &Pythia8::WeightsBase::weightValues);
		cl.def_readwrite("weightNames", &Pythia8::WeightsBase::weightNames);
		cl.def("init", (void (Pythia8::WeightsBase::*)()) &Pythia8::WeightsBase::init, "C++: Pythia8::WeightsBase::init() --> void");
		cl.def("init", (void (Pythia8::WeightsBase::*)(bool)) &Pythia8::WeightsBase::init, "C++: Pythia8::WeightsBase::init(bool) --> void", pybind11::arg(""));
		cl.def("clear", (void (Pythia8::WeightsBase::*)()) &Pythia8::WeightsBase::clear, "C++: Pythia8::WeightsBase::clear() --> void");
		cl.def("bookVectors", (void (Pythia8::WeightsBase::*)(class std::vector<double, class std::allocator<double> >, class std::vector<std::string, class std::allocator<std::string > >)) &Pythia8::WeightsBase::bookVectors, "C++: Pythia8::WeightsBase::bookVectors(class std::vector<double, class std::allocator<double> >, class std::vector<std::string, class std::allocator<std::string > >) --> void", pybind11::arg("weights"), pybind11::arg("names"));
		cl.def("collectWeightValues", [](Pythia8::WeightsBase &o, class std::vector<double, class std::allocator<double> > & a0) -> void { return o.collectWeightValues(a0); }, "", pybind11::arg("outputWeights"));
		cl.def("collectWeightValues", (void (Pythia8::WeightsBase::*)(class std::vector<double, class std::allocator<double> > &, double)) &Pythia8::WeightsBase::collectWeightValues, "C++: Pythia8::WeightsBase::collectWeightValues(class std::vector<double, class std::allocator<double> > &, double) --> void", pybind11::arg("outputWeights"), pybind11::arg("norm"));
		cl.def("collectWeightNames", (void (Pythia8::WeightsBase::*)(class std::vector<std::string, class std::allocator<std::string > > &)) &Pythia8::WeightsBase::collectWeightNames, "C++: Pythia8::WeightsBase::collectWeightNames(class std::vector<std::string, class std::allocator<std::string > > &) --> void", pybind11::arg("outputNames"));
		cl.def("getWeightsName", (std::string (Pythia8::WeightsBase::*)(int) const) &Pythia8::WeightsBase::getWeightsName, "C++: Pythia8::WeightsBase::getWeightsName(int) const --> std::string", pybind11::arg("iPos"));
		cl.def("getWeightsValue", (double (Pythia8::WeightsBase::*)(int) const) &Pythia8::WeightsBase::getWeightsValue, "C++: Pythia8::WeightsBase::getWeightsValue(int) const --> double", pybind11::arg("iPos"));
		cl.def("getWeightsSize", (int (Pythia8::WeightsBase::*)() const) &Pythia8::WeightsBase::getWeightsSize, "C++: Pythia8::WeightsBase::getWeightsSize() const --> int");
		cl.def("bookWeight", [](Pythia8::WeightsBase &o, class std::basic_string<char> const & a0) -> void { return o.bookWeight(a0); }, "", pybind11::arg("name"));
		cl.def("bookWeight", (void (Pythia8::WeightsBase::*)(std::string, double)) &Pythia8::WeightsBase::bookWeight, "C++: Pythia8::WeightsBase::bookWeight(std::string, double) --> void", pybind11::arg("name"), pybind11::arg("defaultValue"));
		cl.def("setValueByIndex", (void (Pythia8::WeightsBase::*)(int, double)) &Pythia8::WeightsBase::setValueByIndex, "C++: Pythia8::WeightsBase::setValueByIndex(int, double) --> void", pybind11::arg("iPos"), pybind11::arg("val"));
		cl.def("setValueByName", (void (Pythia8::WeightsBase::*)(std::string, double)) &Pythia8::WeightsBase::setValueByName, "C++: Pythia8::WeightsBase::setValueByName(std::string, double) --> void", pybind11::arg("name"), pybind11::arg("val"));
		cl.def("reweightValueByIndex", (void (Pythia8::WeightsBase::*)(int, double)) &Pythia8::WeightsBase::reweightValueByIndex, "C++: Pythia8::WeightsBase::reweightValueByIndex(int, double) --> void", pybind11::arg("iPos"), pybind11::arg("val"));
		cl.def("reweightValueByName", (void (Pythia8::WeightsBase::*)(std::string, double)) &Pythia8::WeightsBase::reweightValueByName, "C++: Pythia8::WeightsBase::reweightValueByName(std::string, double) --> void", pybind11::arg("name"), pybind11::arg("val"));
		cl.def("findIndexOfName", (int (Pythia8::WeightsBase::*)(std::string)) &Pythia8::WeightsBase::findIndexOfName, "C++: Pythia8::WeightsBase::findIndexOfName(std::string) --> int", pybind11::arg("name"));
		cl.def("setPtrs", (void (Pythia8::WeightsBase::*)(class Pythia8::Info *)) &Pythia8::WeightsBase::setPtrs, "C++: Pythia8::WeightsBase::setPtrs(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("assign", (class Pythia8::WeightsBase & (Pythia8::WeightsBase::*)(const class Pythia8::WeightsBase &)) &Pythia8::WeightsBase::operator=, "C++: Pythia8::WeightsBase::operator=(const class Pythia8::WeightsBase &) --> class Pythia8::WeightsBase &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::WeightsShower file:Pythia8/Weights.h line:145
		pybind11::class_<Pythia8::WeightsShower, std::shared_ptr<Pythia8::WeightsShower>, PyCallBack_Pythia8_WeightsShower, Pythia8::WeightsBase> cl(M("Pythia8"), "WeightsShower", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](PyCallBack_Pythia8_WeightsShower const &o){ return new PyCallBack_Pythia8_WeightsShower(o); } ) );
		cl.def( pybind11::init( [](Pythia8::WeightsShower const &o){ return new Pythia8::WeightsShower(o); } ) );
		cl.def( pybind11::init( [](){ return new Pythia8::WeightsShower(); }, [](){ return new PyCallBack_Pythia8_WeightsShower(); } ) );
		cl.def("init", (void (Pythia8::WeightsShower::*)(bool)) &Pythia8::WeightsShower::init, "C++: Pythia8::WeightsShower::init(bool) --> void", pybind11::arg(""));
		cl.def("initWeightGroups", [](Pythia8::WeightsShower &o) -> void { return o.initWeightGroups(); }, "");
		cl.def("initWeightGroups", (void (Pythia8::WeightsShower::*)(bool)) &Pythia8::WeightsShower::initWeightGroups, "C++: Pythia8::WeightsShower::initWeightGroups(bool) --> void", pybind11::arg(""));
		cl.def("nWeightGroups", (int (Pythia8::WeightsShower::*)() const) &Pythia8::WeightsShower::nWeightGroups, "C++: Pythia8::WeightsShower::nWeightGroups() const --> int");
		cl.def("getGroupName", (std::string (Pythia8::WeightsShower::*)(int) const) &Pythia8::WeightsShower::getGroupName, "C++: Pythia8::WeightsShower::getGroupName(int) const --> std::string", pybind11::arg(""));
		cl.def("getGroupWeight", (double (Pythia8::WeightsShower::*)(int) const) &Pythia8::WeightsShower::getGroupWeight, "C++: Pythia8::WeightsShower::getGroupWeight(int) const --> double", pybind11::arg(""));
		cl.def("assign", (class Pythia8::WeightsShower & (Pythia8::WeightsShower::*)(const class Pythia8::WeightsShower &)) &Pythia8::WeightsShower::operator=, "C++: Pythia8::WeightsShower::operator=(const class Pythia8::WeightsShower &) --> class Pythia8::WeightsShower &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::WeightsSimpleShower file:Pythia8/Weights.h line:164
		pybind11::class_<Pythia8::WeightsSimpleShower, std::shared_ptr<Pythia8::WeightsSimpleShower>, PyCallBack_Pythia8_WeightsSimpleShower, Pythia8::WeightsShower> cl(M("Pythia8"), "WeightsSimpleShower", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](PyCallBack_Pythia8_WeightsSimpleShower const &o){ return new PyCallBack_Pythia8_WeightsSimpleShower(o); } ) );
		cl.def( pybind11::init( [](Pythia8::WeightsSimpleShower const &o){ return new Pythia8::WeightsSimpleShower(o); } ) );
		cl.def( pybind11::init( [](){ return new Pythia8::WeightsSimpleShower(); }, [](){ return new PyCallBack_Pythia8_WeightsSimpleShower(); } ) );
		cl.def_readwrite("varPDFplus", &Pythia8::WeightsSimpleShower::varPDFplus);
		cl.def_readwrite("varPDFminus", &Pythia8::WeightsSimpleShower::varPDFminus);
		cl.def_readwrite("varPDFmember", &Pythia8::WeightsSimpleShower::varPDFmember);
		cl.def_readwrite("uniqueShowerVars", &Pythia8::WeightsSimpleShower::uniqueShowerVars);
		cl.def_readwrite("enhanceFactors", &Pythia8::WeightsSimpleShower::enhanceFactors);
		cl.def_readwrite("externalVariations", &Pythia8::WeightsSimpleShower::externalVariations);
		cl.def_readwrite("externalVarNames", &Pythia8::WeightsSimpleShower::externalVarNames);
		cl.def_readwrite("externalGroupNames", &Pythia8::WeightsSimpleShower::externalGroupNames);
		cl.def_readwrite("initialNameSave", &Pythia8::WeightsSimpleShower::initialNameSave);
		cl.def_readwrite("externalMap", &Pythia8::WeightsSimpleShower::externalMap);
		cl.def_readwrite("externalVariationsSize", &Pythia8::WeightsSimpleShower::externalVariationsSize);
		cl.def_readwrite("mergingVarNames", &Pythia8::WeightsSimpleShower::mergingVarNames);
		cl.def_readwrite("pTEnhanced", &Pythia8::WeightsSimpleShower::pTEnhanced);
		cl.def_readwrite("wtEnhanced", &Pythia8::WeightsSimpleShower::wtEnhanced);
		cl.def("init", (void (Pythia8::WeightsSimpleShower::*)(bool)) &Pythia8::WeightsSimpleShower::init, "C++: Pythia8::WeightsSimpleShower::init(bool) --> void", pybind11::arg("doMerging"));
		cl.def("collectWeightNames", (void (Pythia8::WeightsSimpleShower::*)(class std::vector<std::string, class std::allocator<std::string > > &)) &Pythia8::WeightsSimpleShower::collectWeightNames, "C++: Pythia8::WeightsSimpleShower::collectWeightNames(class std::vector<std::string, class std::allocator<std::string > > &) --> void", pybind11::arg("outputNames"));
		cl.def("collectWeightValues", [](Pythia8::WeightsSimpleShower &o, class std::vector<double, class std::allocator<double> > & a0) -> void { return o.collectWeightValues(a0); }, "", pybind11::arg("outputWeights"));
		cl.def("collectWeightValues", (void (Pythia8::WeightsSimpleShower::*)(class std::vector<double, class std::allocator<double> > &, double)) &Pythia8::WeightsSimpleShower::collectWeightValues, "C++: Pythia8::WeightsSimpleShower::collectWeightValues(class std::vector<double, class std::allocator<double> > &, double) --> void", pybind11::arg("outputWeights"), pybind11::arg("norm"));
		cl.def("initWeightGroups", [](Pythia8::WeightsSimpleShower &o) -> void { return o.initWeightGroups(); }, "");
		cl.def("initWeightGroups", (void (Pythia8::WeightsSimpleShower::*)(bool)) &Pythia8::WeightsSimpleShower::initWeightGroups, "C++: Pythia8::WeightsSimpleShower::initWeightGroups(bool) --> void", pybind11::arg(""));
		cl.def("getGroupName", (std::string (Pythia8::WeightsSimpleShower::*)(int) const) &Pythia8::WeightsSimpleShower::getGroupName, "C++: Pythia8::WeightsSimpleShower::getGroupName(int) const --> std::string", pybind11::arg("iGN"));
		cl.def("getGroupWeight", (double (Pythia8::WeightsSimpleShower::*)(int) const) &Pythia8::WeightsSimpleShower::getGroupWeight, "C++: Pythia8::WeightsSimpleShower::getGroupWeight(int) const --> double", pybind11::arg("iGW"));
		cl.def("nWeightGroups", (int (Pythia8::WeightsSimpleShower::*)() const) &Pythia8::WeightsSimpleShower::nWeightGroups, "C++: Pythia8::WeightsSimpleShower::nWeightGroups() const --> int");
		cl.def("initUniqueShowerVars", (bool (Pythia8::WeightsSimpleShower::*)()) &Pythia8::WeightsSimpleShower::initUniqueShowerVars, "C++: Pythia8::WeightsSimpleShower::initUniqueShowerVars() --> bool");
		cl.def("setEnhancedTrial", (void (Pythia8::WeightsSimpleShower::*)(double, double)) &Pythia8::WeightsSimpleShower::setEnhancedTrial, "C++: Pythia8::WeightsSimpleShower::setEnhancedTrial(double, double) --> void", pybind11::arg("pTIn"), pybind11::arg("wtIn"));
		cl.def("getEnhancedTrialPT", (double (Pythia8::WeightsSimpleShower::*)()) &Pythia8::WeightsSimpleShower::getEnhancedTrialPT, "C++: Pythia8::WeightsSimpleShower::getEnhancedTrialPT() --> double");
		cl.def("getEnhancedTrialWeight", (double (Pythia8::WeightsSimpleShower::*)()) &Pythia8::WeightsSimpleShower::getEnhancedTrialWeight, "C++: Pythia8::WeightsSimpleShower::getEnhancedTrialWeight() --> double");
		cl.def("getUniqueShowerVars", (class std::vector<std::string, class std::allocator<std::string > > (Pythia8::WeightsSimpleShower::*)()) &Pythia8::WeightsSimpleShower::getUniqueShowerVars, "C++: Pythia8::WeightsSimpleShower::getUniqueShowerVars() --> class std::vector<std::string, class std::allocator<std::string > >");
		cl.def("getUniqueShowerVars", (class std::vector<std::string, class std::allocator<std::string > > (Pythia8::WeightsSimpleShower::*)(class std::vector<std::string, class std::allocator<std::string > >)) &Pythia8::WeightsSimpleShower::getUniqueShowerVars, "C++: Pythia8::WeightsSimpleShower::getUniqueShowerVars(class std::vector<std::string, class std::allocator<std::string > >) --> class std::vector<std::string, class std::allocator<std::string > >", pybind11::arg("keys"));
		cl.def("initEnhanceFactors", (bool (Pythia8::WeightsSimpleShower::*)()) &Pythia8::WeightsSimpleShower::initEnhanceFactors, "C++: Pythia8::WeightsSimpleShower::initEnhanceFactors() --> bool");
		cl.def("getEnhanceFactors", (class std::unordered_map<std::string, double, struct std::hash<std::string>, struct std::equal_to<std::string >, class std::allocator<struct std::pair<const std::string, double> > > (Pythia8::WeightsSimpleShower::*)()) &Pythia8::WeightsSimpleShower::getEnhanceFactors, "C++: Pythia8::WeightsSimpleShower::getEnhanceFactors() --> class std::unordered_map<std::string, double, struct std::hash<std::string>, struct std::equal_to<std::string >, class std::allocator<struct std::pair<const std::string, double> > >");
		cl.def("getInitialName", (std::string (Pythia8::WeightsSimpleShower::*)(int) const) &Pythia8::WeightsSimpleShower::getInitialName, "C++: Pythia8::WeightsSimpleShower::getInitialName(int) const --> std::string", pybind11::arg("iG"));
		cl.def("getMuRWeightVector", (class std::vector<double, class std::allocator<double> > (Pythia8::WeightsSimpleShower::*)()) &Pythia8::WeightsSimpleShower::getMuRWeightVector, "C++: Pythia8::WeightsSimpleShower::getMuRWeightVector() --> class std::vector<double, class std::allocator<double> >");
		cl.def("assign", (class Pythia8::WeightsSimpleShower & (Pythia8::WeightsSimpleShower::*)(const class Pythia8::WeightsSimpleShower &)) &Pythia8::WeightsSimpleShower::operator=, "C++: Pythia8::WeightsSimpleShower::operator=(const class Pythia8::WeightsSimpleShower &) --> class Pythia8::WeightsSimpleShower &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
