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
#include <functional>
#include <iterator>
#include <map>
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

// Pythia8::WeightsMerging file:Pythia8/Weights.h line:236
struct PyCallBack_Pythia8_WeightsMerging : public Pythia8::WeightsMerging {
	using Pythia8::WeightsMerging::WeightsMerging;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsMerging::init();
	}
	void clear() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "clear");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsMerging::clear();
	}
	void bookVectors(class std::vector<double, class std::allocator<double> > a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "bookVectors");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsMerging::bookVectors(a0, a1);
	}
	double getWeightsValue(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "getWeightsValue");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return WeightsMerging::getWeightsValue(a0);
	}
	void reweightValueByIndex(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "reweightValueByIndex");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsMerging::reweightValueByIndex(a0, a1);
	}
	void reweightValueByName(class std::basic_string<char> a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "reweightValueByName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsMerging::reweightValueByName(a0, a1);
	}
	void collectWeightNames(class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "collectWeightNames");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsMerging::collectWeightNames(a0);
	}
	void collectWeightValues(class std::vector<double, class std::allocator<double> > & a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "collectWeightValues");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsMerging::collectWeightValues(a0, a1);
	}
	void init(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsMerging *>(this), "init");
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
};

// Pythia8::WeightsLHEF file:Pythia8/Weights.h line:314
struct PyCallBack_Pythia8_WeightsLHEF : public Pythia8::WeightsLHEF {
	using Pythia8::WeightsLHEF::WeightsLHEF;

	void clear() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "clear");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsLHEF::clear();
	}
	void bookVectors(class std::vector<double, class std::allocator<double> > a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "bookVectors");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsLHEF::bookVectors(a0, a1);
	}
	void collectWeightValues(class std::vector<double, class std::allocator<double> > & a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "collectWeightValues");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsLHEF::collectWeightValues(a0, a1);
	}
	void collectWeightNames(class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "collectWeightNames");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsLHEF::collectWeightNames(a0);
	}
	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "init");
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
	double getWeightsValue(int a0) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "getWeightsValue");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "reweightValueByIndex");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsLHEF *>(this), "reweightValueByName");
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

// Pythia8::WeightsFragmentation file:Pythia8/Weights.h line:354
struct PyCallBack_Pythia8_WeightsFragmentation : public Pythia8::WeightsFragmentation {
	using Pythia8::WeightsFragmentation::WeightsFragmentation;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsFragmentation::init();
	}
	void clear() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "clear");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsFragmentation::clear();
	}
	void collectWeightNames(class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "collectWeightNames");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsFragmentation::collectWeightNames(a0);
	}
	void collectWeightValues(class std::vector<double, class std::allocator<double> > & a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "collectWeightValues");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return WeightsFragmentation::collectWeightValues(a0, a1);
	}
	void init(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "init");
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
	void bookVectors(class std::vector<double, class std::allocator<double> > a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "bookVectors");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "getWeightsValue");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "reweightValueByIndex");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WeightsFragmentation *>(this), "reweightValueByName");
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

void bind_Pythia8_Weights(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::WeightsMerging file:Pythia8/Weights.h line:236
		pybind11::class_<Pythia8::WeightsMerging, std::shared_ptr<Pythia8::WeightsMerging>, PyCallBack_Pythia8_WeightsMerging, Pythia8::WeightsBase> cl(M("Pythia8"), "WeightsMerging", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](PyCallBack_Pythia8_WeightsMerging const &o){ return new PyCallBack_Pythia8_WeightsMerging(o); } ) );
		cl.def( pybind11::init( [](Pythia8::WeightsMerging const &o){ return new Pythia8::WeightsMerging(o); } ) );
		cl.def( pybind11::init( [](){ return new Pythia8::WeightsMerging(); }, [](){ return new PyCallBack_Pythia8_WeightsMerging(); } ) );
		cl.def_readwrite("muRVarLHEFindex", &Pythia8::WeightsMerging::muRVarLHEFindex);
		cl.def_readwrite("weightValuesFirst", &Pythia8::WeightsMerging::weightValuesFirst);
		cl.def_readwrite("weightValuesP", &Pythia8::WeightsMerging::weightValuesP);
		cl.def_readwrite("weightValuesPC", &Pythia8::WeightsMerging::weightValuesPC);
		cl.def_readwrite("weightValuesFirstP", &Pythia8::WeightsMerging::weightValuesFirstP);
		cl.def_readwrite("weightValuesFirstPC", &Pythia8::WeightsMerging::weightValuesFirstPC);
		cl.def_readwrite("isNLO", &Pythia8::WeightsMerging::isNLO);
		cl.def("init", (void (Pythia8::WeightsMerging::*)()) &Pythia8::WeightsMerging::init, "C++: Pythia8::WeightsMerging::init() --> void");
		cl.def("clear", (void (Pythia8::WeightsMerging::*)()) &Pythia8::WeightsMerging::clear, "C++: Pythia8::WeightsMerging::clear() --> void");
		cl.def("bookWeight", (void (Pythia8::WeightsMerging::*)(std::string, double, double)) &Pythia8::WeightsMerging::bookWeight, "C++: Pythia8::WeightsMerging::bookWeight(std::string, double, double) --> void", pybind11::arg("name"), pybind11::arg("value"), pybind11::arg("valueFirst"));
		cl.def("bookVectors", (void (Pythia8::WeightsMerging::*)(class std::vector<double, class std::allocator<double> >, class std::vector<std::string, class std::allocator<std::string > >)) &Pythia8::WeightsMerging::bookVectors, "C++: Pythia8::WeightsMerging::bookVectors(class std::vector<double, class std::allocator<double> >, class std::vector<std::string, class std::allocator<std::string > >) --> void", pybind11::arg("weights"), pybind11::arg("names"));
		cl.def("bookVectors", (void (Pythia8::WeightsMerging::*)(class std::vector<double, class std::allocator<double> >, class std::vector<double, class std::allocator<double> >, class std::vector<std::string, class std::allocator<std::string > >)) &Pythia8::WeightsMerging::bookVectors, "C++: Pythia8::WeightsMerging::bookVectors(class std::vector<double, class std::allocator<double> >, class std::vector<double, class std::allocator<double> >, class std::vector<std::string, class std::allocator<std::string > >) --> void", pybind11::arg("weights"), pybind11::arg("weightsFirst"), pybind11::arg("names"));
		cl.def("getWeightsValue", (double (Pythia8::WeightsMerging::*)(int) const) &Pythia8::WeightsMerging::getWeightsValue, "C++: Pythia8::WeightsMerging::getWeightsValue(int) const --> double", pybind11::arg("iPos"));
		cl.def("getWeightsValueP", (double (Pythia8::WeightsMerging::*)(int) const) &Pythia8::WeightsMerging::getWeightsValueP, "C++: Pythia8::WeightsMerging::getWeightsValueP(int) const --> double", pybind11::arg("iPos"));
		cl.def("getWeightsValuePC", (double (Pythia8::WeightsMerging::*)(int) const) &Pythia8::WeightsMerging::getWeightsValuePC, "C++: Pythia8::WeightsMerging::getWeightsValuePC(int) const --> double", pybind11::arg("iPos"));
		cl.def("reweightValueByIndex", (void (Pythia8::WeightsMerging::*)(int, double)) &Pythia8::WeightsMerging::reweightValueByIndex, "C++: Pythia8::WeightsMerging::reweightValueByIndex(int, double) --> void", pybind11::arg("iPos"), pybind11::arg("val"));
		cl.def("reweightValueByName", (void (Pythia8::WeightsMerging::*)(std::string, double)) &Pythia8::WeightsMerging::reweightValueByName, "C++: Pythia8::WeightsMerging::reweightValueByName(std::string, double) --> void", pybind11::arg("name"), pybind11::arg("val"));
		cl.def("setValueFirstByIndex", (void (Pythia8::WeightsMerging::*)(int, double)) &Pythia8::WeightsMerging::setValueFirstByIndex, "C++: Pythia8::WeightsMerging::setValueFirstByIndex(int, double) --> void", pybind11::arg("iPos"), pybind11::arg("val"));
		cl.def("setValueFirstByName", (void (Pythia8::WeightsMerging::*)(std::string, double)) &Pythia8::WeightsMerging::setValueFirstByName, "C++: Pythia8::WeightsMerging::setValueFirstByName(std::string, double) --> void", pybind11::arg("name"), pybind11::arg("val"));
		cl.def("setValueVector", (void (Pythia8::WeightsMerging::*)(class std::vector<double, class std::allocator<double> >)) &Pythia8::WeightsMerging::setValueVector, "C++: Pythia8::WeightsMerging::setValueVector(class std::vector<double, class std::allocator<double> >) --> void", pybind11::arg("ValueVector"));
		cl.def("setValueFirstVector", (void (Pythia8::WeightsMerging::*)(class std::vector<double, class std::allocator<double> >)) &Pythia8::WeightsMerging::setValueFirstVector, "C++: Pythia8::WeightsMerging::setValueFirstVector(class std::vector<double, class std::allocator<double> >) --> void", pybind11::arg("ValueFirstVector"));
		cl.def("getMuRVarFactors", (class std::vector<double, class std::allocator<double> > (Pythia8::WeightsMerging::*)()) &Pythia8::WeightsMerging::getMuRVarFactors, "C++: Pythia8::WeightsMerging::getMuRVarFactors() --> class std::vector<double, class std::allocator<double> >");
		cl.def("setLHEFvariationMapping", (void (Pythia8::WeightsMerging::*)()) &Pythia8::WeightsMerging::setLHEFvariationMapping, "C++: Pythia8::WeightsMerging::setLHEFvariationMapping() --> void");
		cl.def("collectWeightNames", (void (Pythia8::WeightsMerging::*)(class std::vector<std::string, class std::allocator<std::string > > &)) &Pythia8::WeightsMerging::collectWeightNames, "C++: Pythia8::WeightsMerging::collectWeightNames(class std::vector<std::string, class std::allocator<std::string > > &) --> void", pybind11::arg("outputNames"));
		cl.def("collectWeightValues", [](Pythia8::WeightsMerging &o, class std::vector<double, class std::allocator<double> > & a0) -> void { return o.collectWeightValues(a0); }, "", pybind11::arg("outputWeights"));
		cl.def("collectWeightValues", (void (Pythia8::WeightsMerging::*)(class std::vector<double, class std::allocator<double> > &, double)) &Pythia8::WeightsMerging::collectWeightValues, "C++: Pythia8::WeightsMerging::collectWeightValues(class std::vector<double, class std::allocator<double> > &, double) --> void", pybind11::arg("outputWeights"), pybind11::arg("norm"));
		cl.def("assign", (class Pythia8::WeightsMerging & (Pythia8::WeightsMerging::*)(const class Pythia8::WeightsMerging &)) &Pythia8::WeightsMerging::operator=, "C++: Pythia8::WeightsMerging::operator=(const class Pythia8::WeightsMerging &) --> class Pythia8::WeightsMerging &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::WeightsLHEF file:Pythia8/Weights.h line:314
		pybind11::class_<Pythia8::WeightsLHEF, std::shared_ptr<Pythia8::WeightsLHEF>, PyCallBack_Pythia8_WeightsLHEF, Pythia8::WeightsBase> cl(M("Pythia8"), "WeightsLHEF", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](PyCallBack_Pythia8_WeightsLHEF const &o){ return new PyCallBack_Pythia8_WeightsLHEF(o); } ) );
		cl.def( pybind11::init( [](Pythia8::WeightsLHEF const &o){ return new Pythia8::WeightsLHEF(o); } ) );
		cl.def( pybind11::init( [](){ return new Pythia8::WeightsLHEF(); }, [](){ return new PyCallBack_Pythia8_WeightsLHEF(); } ) );
		cl.def_readwrite("centralWeight", &Pythia8::WeightsLHEF::centralWeight);
		cl.def_readwrite("muRvars", &Pythia8::WeightsLHEF::muRvars);
		cl.def("clear", (void (Pythia8::WeightsLHEF::*)()) &Pythia8::WeightsLHEF::clear, "C++: Pythia8::WeightsLHEF::clear() --> void");
		cl.def("bookVectors", (void (Pythia8::WeightsLHEF::*)(class std::vector<double, class std::allocator<double> >, class std::vector<std::string, class std::allocator<std::string > >)) &Pythia8::WeightsLHEF::bookVectors, "C++: Pythia8::WeightsLHEF::bookVectors(class std::vector<double, class std::allocator<double> >, class std::vector<std::string, class std::allocator<std::string > >) --> void", pybind11::arg("weights"), pybind11::arg("names"));
		cl.def("collectWeightValues", [](Pythia8::WeightsLHEF &o, class std::vector<double, class std::allocator<double> > & a0) -> void { return o.collectWeightValues(a0); }, "", pybind11::arg("outputWeights"));
		cl.def("collectWeightValues", (void (Pythia8::WeightsLHEF::*)(class std::vector<double, class std::allocator<double> > &, double)) &Pythia8::WeightsLHEF::collectWeightValues, "C++: Pythia8::WeightsLHEF::collectWeightValues(class std::vector<double, class std::allocator<double> > &, double) --> void", pybind11::arg("outputWeights"), pybind11::arg("norm"));
		cl.def("collectWeightNames", (void (Pythia8::WeightsLHEF::*)(class std::vector<std::string, class std::allocator<std::string > > &)) &Pythia8::WeightsLHEF::collectWeightNames, "C++: Pythia8::WeightsLHEF::collectWeightNames(class std::vector<std::string, class std::allocator<std::string > > &) --> void", pybind11::arg("outputNames"));
		cl.def("convertNames", (class std::vector<std::string, class std::allocator<std::string > > (Pythia8::WeightsLHEF::*)(class std::vector<std::string, class std::allocator<std::string > >)) &Pythia8::WeightsLHEF::convertNames, "C++: Pythia8::WeightsLHEF::convertNames(class std::vector<std::string, class std::allocator<std::string > >) --> class std::vector<std::string, class std::allocator<std::string > >", pybind11::arg("names"));
		cl.def("identifyVariationsFromLHAinit", (void (Pythia8::WeightsLHEF::*)(class std::map<std::string, struct Pythia8::LHAweight, struct std::less<std::string >, class std::allocator<struct std::pair<const std::string, struct Pythia8::LHAweight> > > *)) &Pythia8::WeightsLHEF::identifyVariationsFromLHAinit, "C++: Pythia8::WeightsLHEF::identifyVariationsFromLHAinit(class std::map<std::string, struct Pythia8::LHAweight, struct std::less<std::string >, class std::allocator<struct std::pair<const std::string, struct Pythia8::LHAweight> > > *) --> void", pybind11::arg("weights"));
		cl.def("assign", (class Pythia8::WeightsLHEF & (Pythia8::WeightsLHEF::*)(const class Pythia8::WeightsLHEF &)) &Pythia8::WeightsLHEF::operator=, "C++: Pythia8::WeightsLHEF::operator=(const class Pythia8::WeightsLHEF &) --> class Pythia8::WeightsLHEF &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::WeightsFragmentation file:Pythia8/Weights.h line:354
		pybind11::class_<Pythia8::WeightsFragmentation, std::shared_ptr<Pythia8::WeightsFragmentation>, PyCallBack_Pythia8_WeightsFragmentation, Pythia8::WeightsBase> cl(M("Pythia8"), "WeightsFragmentation", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](PyCallBack_Pythia8_WeightsFragmentation const &o){ return new PyCallBack_Pythia8_WeightsFragmentation(o); } ) );
		cl.def( pybind11::init( [](Pythia8::WeightsFragmentation const &o){ return new Pythia8::WeightsFragmentation(o); } ) );
		cl.def( pybind11::init( [](){ return new Pythia8::WeightsFragmentation(); }, [](){ return new PyCallBack_Pythia8_WeightsFragmentation(); } ) );

		pybind11::enum_<Pythia8::WeightsFragmentation::FactIndex>(cl, "FactIndex", pybind11::arithmetic(), "")
			.value("Z", Pythia8::WeightsFragmentation::FactIndex::Z)
			.value("Flav", Pythia8::WeightsFragmentation::FactIndex::Flav)
			.value("PT", Pythia8::WeightsFragmentation::FactIndex::PT)
			.export_values();

		cl.def_readwrite("weightParms", &Pythia8::WeightsFragmentation::weightParms);
		cl.def_readwrite("externalGroupNames", &Pythia8::WeightsFragmentation::externalGroupNames);
		cl.def_readwrite("externalMap", &Pythia8::WeightsFragmentation::externalMap);
		cl.def_readwrite("flavBreaks", &Pythia8::WeightsFragmentation::flavBreaks);
		cl.def("init", (void (Pythia8::WeightsFragmentation::*)()) &Pythia8::WeightsFragmentation::init, "C++: Pythia8::WeightsFragmentation::init() --> void");
		cl.def("clear", (void (Pythia8::WeightsFragmentation::*)()) &Pythia8::WeightsFragmentation::clear, "C++: Pythia8::WeightsFragmentation::clear() --> void");
		cl.def("nWeightGroups", (int (Pythia8::WeightsFragmentation::*)() const) &Pythia8::WeightsFragmentation::nWeightGroups, "C++: Pythia8::WeightsFragmentation::nWeightGroups() const --> int");
		cl.def("getGroupName", (std::string (Pythia8::WeightsFragmentation::*)(int) const) &Pythia8::WeightsFragmentation::getGroupName, "C++: Pythia8::WeightsFragmentation::getGroupName(int) const --> std::string", pybind11::arg("iGN"));
		cl.def("getGroupWeight", (double (Pythia8::WeightsFragmentation::*)(int) const) &Pythia8::WeightsFragmentation::getGroupWeight, "C++: Pythia8::WeightsFragmentation::getGroupWeight(int) const --> double", pybind11::arg("iGW"));
		cl.def("collectWeightNames", (void (Pythia8::WeightsFragmentation::*)(class std::vector<std::string, class std::allocator<std::string > > &)) &Pythia8::WeightsFragmentation::collectWeightNames, "C++: Pythia8::WeightsFragmentation::collectWeightNames(class std::vector<std::string, class std::allocator<std::string > > &) --> void", pybind11::arg("outputNames"));
		cl.def("collectWeightValues", [](Pythia8::WeightsFragmentation &o, class std::vector<double, class std::allocator<double> > & a0) -> void { return o.collectWeightValues(a0); }, "", pybind11::arg("outputWeights"));
		cl.def("collectWeightValues", (void (Pythia8::WeightsFragmentation::*)(class std::vector<double, class std::allocator<double> > &, double)) &Pythia8::WeightsFragmentation::collectWeightValues, "C++: Pythia8::WeightsFragmentation::collectWeightValues(class std::vector<double, class std::allocator<double> > &, double) --> void", pybind11::arg("outputWeights"), pybind11::arg("norm"));
		cl.def("flavParms", (class std::vector<double, class std::allocator<double> > (Pythia8::WeightsFragmentation::*)(double, double, double, double)) &Pythia8::WeightsFragmentation::flavParms, "C++: Pythia8::WeightsFragmentation::flavParms(double, double, double, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("xi"), pybind11::arg("rho"), pybind11::arg("x"), pybind11::arg("y"));
		cl.def("flavWeight", (double (Pythia8::WeightsFragmentation::*)(const class std::vector<double, class std::allocator<double> > &)) &Pythia8::WeightsFragmentation::flavWeight, "C++: Pythia8::WeightsFragmentation::flavWeight(const class std::vector<double, class std::allocator<double> > &) --> double", pybind11::arg("parms"));
		cl.def("flavWeight", (double (Pythia8::WeightsFragmentation::*)(const class std::vector<double, class std::allocator<double> > &, const class std::vector<int, class std::allocator<int> > &)) &Pythia8::WeightsFragmentation::flavWeight, "C++: Pythia8::WeightsFragmentation::flavWeight(const class std::vector<double, class std::allocator<double> > &, const class std::vector<int, class std::allocator<int> > &) --> double", pybind11::arg("parms"), pybind11::arg("breaks"));
	}
	{ // Pythia8::WeightContainer file:Pythia8/Weights.h line:433
		pybind11::class_<Pythia8::WeightContainer, std::shared_ptr<Pythia8::WeightContainer>> cl(M("Pythia8"), "WeightContainer", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::WeightContainer(); } ) );
		cl.def( pybind11::init( [](Pythia8::WeightContainer const &o){ return new Pythia8::WeightContainer(o); } ) );
		cl.def_readwrite("weightNominal", &Pythia8::WeightContainer::weightNominal);
		cl.def_readwrite("weightsLHEF", &Pythia8::WeightContainer::weightsLHEF);
		cl.def_readwrite("weightsSimpleShower", &Pythia8::WeightContainer::weightsSimpleShower);
		cl.def_readwrite("weightsMerging", &Pythia8::WeightContainer::weightsMerging);
		cl.def_readwrite("weightsUserHooks", &Pythia8::WeightContainer::weightsUserHooks);
		cl.def("setWeightNominal", (void (Pythia8::WeightContainer::*)(double)) &Pythia8::WeightContainer::setWeightNominal, "C++: Pythia8::WeightContainer::setWeightNominal(double) --> void", pybind11::arg("weightNow"));
		cl.def("collectWeightNominal", (double (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::collectWeightNominal, "C++: Pythia8::WeightContainer::collectWeightNominal() --> double");
		cl.def("numberOfWeights", (int (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::numberOfWeights, "C++: Pythia8::WeightContainer::numberOfWeights() --> int");
		cl.def("weightValueByIndex", [](Pythia8::WeightContainer &o) -> double { return o.weightValueByIndex(); }, "");
		cl.def("weightValueByIndex", (double (Pythia8::WeightContainer::*)(int)) &Pythia8::WeightContainer::weightValueByIndex, "C++: Pythia8::WeightContainer::weightValueByIndex(int) --> double", pybind11::arg("key"));
		cl.def("weightNameByIndex", [](Pythia8::WeightContainer &o) -> std::string { return o.weightNameByIndex(); }, "");
		cl.def("weightNameByIndex", (std::string (Pythia8::WeightContainer::*)(int)) &Pythia8::WeightContainer::weightNameByIndex, "C++: Pythia8::WeightContainer::weightNameByIndex(int) --> std::string", pybind11::arg("key"));
		cl.def("weightValueVector", (class std::vector<double, class std::allocator<double> > (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::weightValueVector, "C++: Pythia8::WeightContainer::weightValueVector() --> class std::vector<double, class std::allocator<double> >");
		cl.def("weightNameVector", (class std::vector<std::string, class std::allocator<std::string > > (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::weightNameVector, "C++: Pythia8::WeightContainer::weightNameVector() --> class std::vector<std::string, class std::allocator<std::string > >");
		cl.def("clear", (void (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::clear, "C++: Pythia8::WeightContainer::clear() --> void");
		cl.def("clearTotal", (void (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::clearTotal, "C++: Pythia8::WeightContainer::clearTotal() --> void");
		cl.def("init", (void (Pythia8::WeightContainer::*)(bool)) &Pythia8::WeightContainer::init, "C++: Pythia8::WeightContainer::init(bool) --> void", pybind11::arg("doMerging"));
		cl.def("initPtrs", (void (Pythia8::WeightContainer::*)(class Pythia8::Info *)) &Pythia8::WeightContainer::initPtrs, "C++: Pythia8::WeightContainer::initPtrs(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("initXsecVec", (void (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::initXsecVec, "C++: Pythia8::WeightContainer::initXsecVec() --> void");
		cl.def("getSampleXsec", (class std::vector<double, class std::allocator<double> > (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::getSampleXsec, "C++: Pythia8::WeightContainer::getSampleXsec() --> class std::vector<double, class std::allocator<double> >");
		cl.def("getTotalXsec", (class std::vector<double, class std::allocator<double> > (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::getTotalXsec, "C++: Pythia8::WeightContainer::getTotalXsec() --> class std::vector<double, class std::allocator<double> >");
		cl.def("getSampleXsecErr", (class std::vector<double, class std::allocator<double> > (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::getSampleXsecErr, "C++: Pythia8::WeightContainer::getSampleXsecErr() --> class std::vector<double, class std::allocator<double> >");
		cl.def("getTotalXsecErr", (class std::vector<double, class std::allocator<double> > (Pythia8::WeightContainer::*)()) &Pythia8::WeightContainer::getTotalXsecErr, "C++: Pythia8::WeightContainer::getTotalXsecErr() --> class std::vector<double, class std::allocator<double> >");
		cl.def("accumulateXsec", [](Pythia8::WeightContainer &o) -> void { return o.accumulateXsec(); }, "");
		cl.def("accumulateXsec", (void (Pythia8::WeightContainer::*)(double)) &Pythia8::WeightContainer::accumulateXsec, "C++: Pythia8::WeightContainer::accumulateXsec(double) --> void", pybind11::arg("norm"));
	}
}
