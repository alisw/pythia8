#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/SusyLesHouches.h>
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

// Pythia8::SigmaTotAux file:Pythia8/SigmaTotal.h line:33
struct PyCallBack_Pythia8_SigmaTotAux : public Pythia8::SigmaTotAux {
	using Pythia8::SigmaTotAux::SigmaTotAux;

	void init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SigmaTotAux::init\"");
	}
	bool calcTot(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "calcTot");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::calcTot(a0, a1, a2);
	}
	bool calcTotEl(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "calcTotEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::calcTotEl(a0, a1, a2, a3, a4);
	}
	double dsigmaEl(double a0, bool a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "dsigmaEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotAux::dsigmaEl(a0, a1, a2);
	}
	bool calcDiff(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "calcDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::calcDiff(a0, a1, a2, a3, a4);
	}
	double dsigmaSD(double a0, double a1, bool a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "dsigmaSD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotAux::dsigmaSD(a0, a1, a2, a3);
	}
	bool splitDiff() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "splitDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::splitDiff();
	}
	double dsigmaDD(double a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "dsigmaDD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotAux::dsigmaDD(a0, a1, a2, a3);
	}
	double dsigmaCD(double a0, double a1, double a2, double a3, int a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "dsigmaCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotAux::dsigmaCD(a0, a1, a2, a3, a4);
	}
	double mMinCD() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "mMinCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotAux::mMinCD();
	}
	bool initCoulomb(class Pythia8::Settings & a0, class Pythia8::ParticleData * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "initCoulomb");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::initCoulomb(a0, a1);
	}
	bool addCoulomb() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "addCoulomb");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::addCoulomb();
	}
	double dsigmaElCoulomb(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotAux *>(this), "dsigmaElCoulomb");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotAux::dsigmaElCoulomb(a0);
	}
};

// Pythia8::SigmaTotal file:Pythia8/SigmaTotal.h line:141
struct PyCallBack_Pythia8_SigmaTotal : public Pythia8::SigmaTotal {
	using Pythia8::SigmaTotal::SigmaTotal;

	bool splitDiff() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotal *>(this), "splitDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotal::splitDiff();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotal *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotal *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotal *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotal *>(this), "onStat");
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

// Pythia8::SigmaTotOwn file:Pythia8/SigmaTotal.h line:242
struct PyCallBack_Pythia8_SigmaTotOwn : public Pythia8::SigmaTotOwn {
	using Pythia8::SigmaTotOwn::SigmaTotOwn;

	void init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SigmaTotOwn::init(a0);
	}
	bool calcTotEl(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "calcTotEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotOwn::calcTotEl(a0, a1, a2, a3, a4);
	}
	double dsigmaEl(double a0, bool a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "dsigmaEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotOwn::dsigmaEl(a0, a1, a2);
	}
	bool calcDiff(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "calcDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotOwn::calcDiff(a0, a1, a2, a3, a4);
	}
	double dsigmaSD(double a0, double a1, bool a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "dsigmaSD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotOwn::dsigmaSD(a0, a1, a2, a3);
	}
	double dsigmaDD(double a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "dsigmaDD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotOwn::dsigmaDD(a0, a1, a2, a3);
	}
	double dsigmaCD(double a0, double a1, double a2, double a3, int a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "dsigmaCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotOwn::dsigmaCD(a0, a1, a2, a3, a4);
	}
	double mMinCD() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "mMinCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotOwn::mMinCD();
	}
	bool calcTot(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "calcTot");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::calcTot(a0, a1, a2);
	}
	bool splitDiff() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "splitDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::splitDiff();
	}
	bool initCoulomb(class Pythia8::Settings & a0, class Pythia8::ParticleData * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "initCoulomb");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::initCoulomb(a0, a1);
	}
	bool addCoulomb() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "addCoulomb");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaTotAux::addCoulomb();
	}
	double dsigmaElCoulomb(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaTotOwn *>(this), "dsigmaElCoulomb");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaTotAux::dsigmaElCoulomb(a0);
	}
};

void bind_Pythia8_SigmaTotal(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::SigmaTotAux file:Pythia8/SigmaTotal.h line:33
		pybind11::class_<Pythia8::SigmaTotAux, std::shared_ptr<Pythia8::SigmaTotAux>, PyCallBack_Pythia8_SigmaTotAux> cl(M("Pythia8"), "SigmaTotAux", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_SigmaTotAux(); } ) );
		cl.def(pybind11::init<PyCallBack_Pythia8_SigmaTotAux const &>());
		cl.def_readwrite("isExpEl", &Pythia8::SigmaTotAux::isExpEl);
		cl.def_readwrite("hasCou", &Pythia8::SigmaTotAux::hasCou);
		cl.def_readwrite("sigTot", &Pythia8::SigmaTotAux::sigTot);
		cl.def_readwrite("rhoOwn", &Pythia8::SigmaTotAux::rhoOwn);
		cl.def_readwrite("sigEl", &Pythia8::SigmaTotAux::sigEl);
		cl.def_readwrite("bEl", &Pythia8::SigmaTotAux::bEl);
		cl.def_readwrite("sigTotCou", &Pythia8::SigmaTotAux::sigTotCou);
		cl.def_readwrite("sigElCou", &Pythia8::SigmaTotAux::sigElCou);
		cl.def_readwrite("sigXB", &Pythia8::SigmaTotAux::sigXB);
		cl.def_readwrite("sigAX", &Pythia8::SigmaTotAux::sigAX);
		cl.def_readwrite("sigXX", &Pythia8::SigmaTotAux::sigXX);
		cl.def_readwrite("sigAXB", &Pythia8::SigmaTotAux::sigAXB);
		cl.def_readwrite("sigNDtmp", &Pythia8::SigmaTotAux::sigNDtmp);
		cl.def_readwrite("idA", &Pythia8::SigmaTotAux::idA);
		cl.def_readwrite("idB", &Pythia8::SigmaTotAux::idB);
		cl.def_readwrite("tryCoulomb", &Pythia8::SigmaTotAux::tryCoulomb);
		cl.def_readwrite("chgSgn", &Pythia8::SigmaTotAux::chgSgn);
		cl.def_readwrite("tAbsMin", &Pythia8::SigmaTotAux::tAbsMin);
		cl.def_readwrite("lambda", &Pythia8::SigmaTotAux::lambda);
		cl.def_readwrite("phaseCst", &Pythia8::SigmaTotAux::phaseCst);
		cl.def("init", (void (Pythia8::SigmaTotAux::*)(class Pythia8::Info *)) &Pythia8::SigmaTotAux::init, "C++: Pythia8::SigmaTotAux::init(class Pythia8::Info *) --> void", pybind11::arg(""));
		cl.def("calcTot", (bool (Pythia8::SigmaTotAux::*)(int, int, double)) &Pythia8::SigmaTotAux::calcTot, "C++: Pythia8::SigmaTotAux::calcTot(int, int, double) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("calcTotEl", (bool (Pythia8::SigmaTotAux::*)(int, int, double, double, double)) &Pythia8::SigmaTotAux::calcTotEl, "C++: Pythia8::SigmaTotAux::calcTotEl(int, int, double, double, double) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaEl", [](Pythia8::SigmaTotAux &o, double const & a0) -> double { return o.dsigmaEl(a0); }, "", pybind11::arg(""));
		cl.def("dsigmaEl", [](Pythia8::SigmaTotAux &o, double const & a0, bool const & a1) -> double { return o.dsigmaEl(a0, a1); }, "", pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaEl", (double (Pythia8::SigmaTotAux::*)(double, bool, bool)) &Pythia8::SigmaTotAux::dsigmaEl, "C++: Pythia8::SigmaTotAux::dsigmaEl(double, bool, bool) --> double", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("calcDiff", (bool (Pythia8::SigmaTotAux::*)(int, int, double, double, double)) &Pythia8::SigmaTotAux::calcDiff, "C++: Pythia8::SigmaTotAux::calcDiff(int, int, double, double, double) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaSD", [](Pythia8::SigmaTotAux &o, double const & a0, double const & a1) -> double { return o.dsigmaSD(a0, a1); }, "", pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaSD", [](Pythia8::SigmaTotAux &o, double const & a0, double const & a1, bool const & a2) -> double { return o.dsigmaSD(a0, a1, a2); }, "", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaSD", (double (Pythia8::SigmaTotAux::*)(double, double, bool, int)) &Pythia8::SigmaTotAux::dsigmaSD, "C++: Pythia8::SigmaTotAux::dsigmaSD(double, double, bool, int) --> double", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("splitDiff", (bool (Pythia8::SigmaTotAux::*)()) &Pythia8::SigmaTotAux::splitDiff, "C++: Pythia8::SigmaTotAux::splitDiff() --> bool");
		cl.def("dsigmaDD", [](Pythia8::SigmaTotAux &o, double const & a0, double const & a1, double const & a2) -> double { return o.dsigmaDD(a0, a1, a2); }, "", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaDD", (double (Pythia8::SigmaTotAux::*)(double, double, double, int)) &Pythia8::SigmaTotAux::dsigmaDD, "C++: Pythia8::SigmaTotAux::dsigmaDD(double, double, double, int) --> double", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaCD", (double (Pythia8::SigmaTotAux::*)(double, double, double, double, int)) &Pythia8::SigmaTotAux::dsigmaCD, "C++: Pythia8::SigmaTotAux::dsigmaCD(double, double, double, double, int) --> double", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("mMinCD", (double (Pythia8::SigmaTotAux::*)()) &Pythia8::SigmaTotAux::mMinCD, "C++: Pythia8::SigmaTotAux::mMinCD() --> double");
		cl.def("tRange", (struct std::pair<double, double> (Pythia8::SigmaTotAux::*)(double, double, double, double, double)) &Pythia8::SigmaTotAux::tRange, "C++: Pythia8::SigmaTotAux::tRange(double, double, double, double, double) --> struct std::pair<double, double>", pybind11::arg("sIn"), pybind11::arg("s1In"), pybind11::arg("s2In"), pybind11::arg("s3In"), pybind11::arg("s4In"));
		cl.def("tInRange", (bool (Pythia8::SigmaTotAux::*)(double, double, double, double, double, double)) &Pythia8::SigmaTotAux::tInRange, "C++: Pythia8::SigmaTotAux::tInRange(double, double, double, double, double, double) --> bool", pybind11::arg("tIn"), pybind11::arg("sIn"), pybind11::arg("s1In"), pybind11::arg("s2In"), pybind11::arg("s3In"), pybind11::arg("s4In"));
		cl.def("pFormFac", (double (Pythia8::SigmaTotAux::*)(double)) &Pythia8::SigmaTotAux::pFormFac, "C++: Pythia8::SigmaTotAux::pFormFac(double) --> double", pybind11::arg("tIn"));
		cl.def("initCoulomb", (bool (Pythia8::SigmaTotAux::*)(class Pythia8::Settings &, class Pythia8::ParticleData *)) &Pythia8::SigmaTotAux::initCoulomb, "C++: Pythia8::SigmaTotAux::initCoulomb(class Pythia8::Settings &, class Pythia8::ParticleData *) --> bool", pybind11::arg("settings"), pybind11::arg("particleDataPtrIn"));
		cl.def("addCoulomb", (bool (Pythia8::SigmaTotAux::*)()) &Pythia8::SigmaTotAux::addCoulomb, "C++: Pythia8::SigmaTotAux::addCoulomb() --> bool");
		cl.def("dsigmaElCoulomb", (double (Pythia8::SigmaTotAux::*)(double)) &Pythia8::SigmaTotAux::dsigmaElCoulomb, "C++: Pythia8::SigmaTotAux::dsigmaElCoulomb(double) --> double", pybind11::arg("t"));
		cl.def("assign", (class Pythia8::SigmaTotAux & (Pythia8::SigmaTotAux::*)(const class Pythia8::SigmaTotAux &)) &Pythia8::SigmaTotAux::operator=, "C++: Pythia8::SigmaTotAux::operator=(const class Pythia8::SigmaTotAux &) --> class Pythia8::SigmaTotAux &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SigmaTotal file:Pythia8/SigmaTotal.h line:141
		pybind11::class_<Pythia8::SigmaTotal, std::shared_ptr<Pythia8::SigmaTotal>, PyCallBack_Pythia8_SigmaTotal, Pythia8::PhysicsBase> cl(M("Pythia8"), "SigmaTotal", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SigmaTotal(); }, [](){ return new PyCallBack_Pythia8_SigmaTotal(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_SigmaTotal const &o){ return new PyCallBack_Pythia8_SigmaTotal(o); } ) );
		cl.def( pybind11::init( [](Pythia8::SigmaTotal const &o){ return new Pythia8::SigmaTotal(o); } ) );
		cl.def("init", (void (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::init, "C++: Pythia8::SigmaTotal::init() --> void");
		cl.def("calc", (bool (Pythia8::SigmaTotal::*)(int, int, double)) &Pythia8::SigmaTotal::calc, "C++: Pythia8::SigmaTotal::calc(int, int, double) --> bool", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"));
		cl.def("hasSigmaTot", (bool (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::hasSigmaTot, "C++: Pythia8::SigmaTotal::hasSigmaTot() --> bool");
		cl.def("sigmaTot", (double (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::sigmaTot, "C++: Pythia8::SigmaTotal::sigmaTot() --> double");
		cl.def("rho", (double (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::rho, "C++: Pythia8::SigmaTotal::rho() --> double");
		cl.def("sigmaEl", (double (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::sigmaEl, "C++: Pythia8::SigmaTotal::sigmaEl() --> double");
		cl.def("bElIsExp", (bool (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::bElIsExp, "C++: Pythia8::SigmaTotal::bElIsExp() --> bool");
		cl.def("bSlopeEl", (double (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::bSlopeEl, "C++: Pythia8::SigmaTotal::bSlopeEl() --> double");
		cl.def("hasCoulomb", (bool (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::hasCoulomb, "C++: Pythia8::SigmaTotal::hasCoulomb() --> bool");
		cl.def("calcTotEl", (bool (Pythia8::SigmaTotal::*)(int, int, double, double, double)) &Pythia8::SigmaTotal::calcTotEl, "C++: Pythia8::SigmaTotal::calcTotEl(int, int, double, double, double) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg("sIn"), pybind11::arg("mAin"), pybind11::arg("mBin"));
		cl.def("dsigmaEl", [](Pythia8::SigmaTotal &o, double const & a0) -> double { return o.dsigmaEl(a0); }, "", pybind11::arg("t"));
		cl.def("dsigmaEl", [](Pythia8::SigmaTotal &o, double const & a0, bool const & a1) -> double { return o.dsigmaEl(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("useCoulomb"));
		cl.def("dsigmaEl", (double (Pythia8::SigmaTotal::*)(double, bool, bool)) &Pythia8::SigmaTotal::dsigmaEl, "C++: Pythia8::SigmaTotal::dsigmaEl(double, bool, bool) --> double", pybind11::arg("t"), pybind11::arg("useCoulomb"), pybind11::arg("onlyPomerons"));
		cl.def("sigmaXB", (double (Pythia8::SigmaTotal::*)() const) &Pythia8::SigmaTotal::sigmaXB, "C++: Pythia8::SigmaTotal::sigmaXB() const --> double");
		cl.def("sigmaAX", (double (Pythia8::SigmaTotal::*)() const) &Pythia8::SigmaTotal::sigmaAX, "C++: Pythia8::SigmaTotal::sigmaAX() const --> double");
		cl.def("sigmaXX", (double (Pythia8::SigmaTotal::*)() const) &Pythia8::SigmaTotal::sigmaXX, "C++: Pythia8::SigmaTotal::sigmaXX() const --> double");
		cl.def("sigmaAXB", (double (Pythia8::SigmaTotal::*)() const) &Pythia8::SigmaTotal::sigmaAXB, "C++: Pythia8::SigmaTotal::sigmaAXB() const --> double");
		cl.def("sigmaND", (double (Pythia8::SigmaTotal::*)() const) &Pythia8::SigmaTotal::sigmaND, "C++: Pythia8::SigmaTotal::sigmaND() const --> double");
		cl.def("dsigmaSD", [](Pythia8::SigmaTotal &o, double const & a0, double const & a1) -> double { return o.dsigmaSD(a0, a1); }, "", pybind11::arg("xi"), pybind11::arg("t"));
		cl.def("dsigmaSD", [](Pythia8::SigmaTotal &o, double const & a0, double const & a1, bool const & a2) -> double { return o.dsigmaSD(a0, a1, a2); }, "", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg("isXB"));
		cl.def("dsigmaSD", (double (Pythia8::SigmaTotal::*)(double, double, bool, int)) &Pythia8::SigmaTotal::dsigmaSD, "C++: Pythia8::SigmaTotal::dsigmaSD(double, double, bool, int) --> double", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg("isXB"), pybind11::arg("step"));
		cl.def("splitDiff", (bool (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::splitDiff, "C++: Pythia8::SigmaTotal::splitDiff() --> bool");
		cl.def("dsigmaDD", [](Pythia8::SigmaTotal &o, double const & a0, double const & a1, double const & a2) -> double { return o.dsigmaDD(a0, a1, a2); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"));
		cl.def("dsigmaDD", (double (Pythia8::SigmaTotal::*)(double, double, double, int)) &Pythia8::SigmaTotal::dsigmaDD, "C++: Pythia8::SigmaTotal::dsigmaDD(double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"), pybind11::arg("step"));
		cl.def("dsigmaCD", [](Pythia8::SigmaTotal &o, double const & a0, double const & a1, double const & a2, double const & a3) -> double { return o.dsigmaCD(a0, a1, a2, a3); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"));
		cl.def("dsigmaCD", (double (Pythia8::SigmaTotal::*)(double, double, double, double, int)) &Pythia8::SigmaTotal::dsigmaCD, "C++: Pythia8::SigmaTotal::dsigmaCD(double, double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"), pybind11::arg("step"));
		cl.def("mMinCD", (double (Pythia8::SigmaTotal::*)()) &Pythia8::SigmaTotal::mMinCD, "C++: Pythia8::SigmaTotal::mMinCD() --> double");
		cl.def("chooseVMDstates", (void (Pythia8::SigmaTotal::*)(int, int, double, int)) &Pythia8::SigmaTotal::chooseVMDstates, "C++: Pythia8::SigmaTotal::chooseVMDstates(int, int, double, int) --> void", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"), pybind11::arg("processCode"));
		cl.def("tRange", (struct std::pair<double, double> (Pythia8::SigmaTotal::*)(double, double, double, double, double)) &Pythia8::SigmaTotal::tRange, "C++: Pythia8::SigmaTotal::tRange(double, double, double, double, double) --> struct std::pair<double, double>", pybind11::arg("sIn"), pybind11::arg("s1In"), pybind11::arg("s2In"), pybind11::arg("s3In"), pybind11::arg("s4In"));
		cl.def("tInRange", (bool (Pythia8::SigmaTotal::*)(double, double, double, double, double, double)) &Pythia8::SigmaTotal::tInRange, "C++: Pythia8::SigmaTotal::tInRange(double, double, double, double, double, double) --> bool", pybind11::arg("tIn"), pybind11::arg("sIn"), pybind11::arg("s1In"), pybind11::arg("s2In"), pybind11::arg("s3In"), pybind11::arg("s4In"));
		cl.def("assign", (class Pythia8::SigmaTotal & (Pythia8::SigmaTotal::*)(const class Pythia8::SigmaTotal &)) &Pythia8::SigmaTotal::operator=, "C++: Pythia8::SigmaTotal::operator=(const class Pythia8::SigmaTotal &) --> class Pythia8::SigmaTotal &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SigmaTotOwn file:Pythia8/SigmaTotal.h line:242
		pybind11::class_<Pythia8::SigmaTotOwn, std::shared_ptr<Pythia8::SigmaTotOwn>, PyCallBack_Pythia8_SigmaTotOwn, Pythia8::SigmaTotAux> cl(M("Pythia8"), "SigmaTotOwn", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SigmaTotOwn(); }, [](){ return new PyCallBack_Pythia8_SigmaTotOwn(); } ) );
		cl.def("init", (void (Pythia8::SigmaTotOwn::*)(class Pythia8::Info *)) &Pythia8::SigmaTotOwn::init, "C++: Pythia8::SigmaTotOwn::init(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("calcTotEl", (bool (Pythia8::SigmaTotOwn::*)(int, int, double, double, double)) &Pythia8::SigmaTotOwn::calcTotEl, "C++: Pythia8::SigmaTotOwn::calcTotEl(int, int, double, double, double) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaEl", [](Pythia8::SigmaTotOwn &o, double const & a0) -> double { return o.dsigmaEl(a0); }, "", pybind11::arg("t"));
		cl.def("dsigmaEl", [](Pythia8::SigmaTotOwn &o, double const & a0, bool const & a1) -> double { return o.dsigmaEl(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("useCoulomb"));
		cl.def("dsigmaEl", (double (Pythia8::SigmaTotOwn::*)(double, bool, bool)) &Pythia8::SigmaTotOwn::dsigmaEl, "C++: Pythia8::SigmaTotOwn::dsigmaEl(double, bool, bool) --> double", pybind11::arg("t"), pybind11::arg("useCoulomb"), pybind11::arg(""));
		cl.def("calcDiff", (bool (Pythia8::SigmaTotOwn::*)(int, int, double, double, double)) &Pythia8::SigmaTotOwn::calcDiff, "C++: Pythia8::SigmaTotOwn::calcDiff(int, int, double, double, double) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg("sIn"), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaSD", [](Pythia8::SigmaTotOwn &o, double const & a0, double const & a1) -> double { return o.dsigmaSD(a0, a1); }, "", pybind11::arg("xi"), pybind11::arg("t"));
		cl.def("dsigmaSD", [](Pythia8::SigmaTotOwn &o, double const & a0, double const & a1, bool const & a2) -> double { return o.dsigmaSD(a0, a1, a2); }, "", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg(""));
		cl.def("dsigmaSD", (double (Pythia8::SigmaTotOwn::*)(double, double, bool, int)) &Pythia8::SigmaTotOwn::dsigmaSD, "C++: Pythia8::SigmaTotOwn::dsigmaSD(double, double, bool, int) --> double", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaDD", [](Pythia8::SigmaTotOwn &o, double const & a0, double const & a1, double const & a2) -> double { return o.dsigmaDD(a0, a1, a2); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"));
		cl.def("dsigmaDD", (double (Pythia8::SigmaTotOwn::*)(double, double, double, int)) &Pythia8::SigmaTotOwn::dsigmaDD, "C++: Pythia8::SigmaTotOwn::dsigmaDD(double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"), pybind11::arg(""));
		cl.def("dsigmaCD", [](Pythia8::SigmaTotOwn &o, double const & a0, double const & a1, double const & a2, double const & a3) -> double { return o.dsigmaCD(a0, a1, a2, a3); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"));
		cl.def("dsigmaCD", (double (Pythia8::SigmaTotOwn::*)(double, double, double, double, int)) &Pythia8::SigmaTotOwn::dsigmaCD, "C++: Pythia8::SigmaTotOwn::dsigmaCD(double, double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"), pybind11::arg(""));
		cl.def("mMinCD", (double (Pythia8::SigmaTotOwn::*)()) &Pythia8::SigmaTotOwn::mMinCD, "C++: Pythia8::SigmaTotOwn::mMinCD() --> double");
		cl.def("assign", (class Pythia8::SigmaTotOwn & (Pythia8::SigmaTotOwn::*)(const class Pythia8::SigmaTotOwn &)) &Pythia8::SigmaTotOwn::operator=, "C++: Pythia8::SigmaTotOwn::operator=(const class Pythia8::SigmaTotOwn &) --> class Pythia8::SigmaTotOwn &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
