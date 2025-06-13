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
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/SusyLesHouches.h>
#include <Pythia8/Weights.h>
#include <complex>
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

// Pythia8::ResonanceGeneric file:Pythia8/ResonanceWidths.h line:152
struct PyCallBack_Pythia8_ResonanceGeneric : public Pythia8::ResonanceGeneric {
	using Pythia8::ResonanceGeneric::ResonanceGeneric;

	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGeneric *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceGeneric::allowCalc();
	}
	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGeneric *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGeneric *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGeneric *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGeneric *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGeneric *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceGmZ file:Pythia8/ResonanceWidths.h line:168
struct PyCallBack_Pythia8_ResonanceGmZ : public Pythia8::ResonanceGmZ {
	using Pythia8::ResonanceGmZ::ResonanceGmZ;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGmZ *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGmZ *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGmZ *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGmZ *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGmZ *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGmZ *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceW file:Pythia8/ResonanceWidths.h line:197
struct PyCallBack_Pythia8_ResonanceW : public Pythia8::ResonanceW {
	using Pythia8::ResonanceW::ResonanceW;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceW *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceW *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceW *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceW *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceW *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceW *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceTop file:Pythia8/ResonanceWidths.h line:224
struct PyCallBack_Pythia8_ResonanceTop : public Pythia8::ResonanceTop {
	using Pythia8::ResonanceTop::ResonanceTop;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceTop *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceTop *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceTop *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceTop *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceTop *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceTop *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceFour file:Pythia8/ResonanceWidths.h line:252
struct PyCallBack_Pythia8_ResonanceFour : public Pythia8::ResonanceFour {
	using Pythia8::ResonanceFour::ResonanceFour;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceFour *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceFour *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceFour *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceFour *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceFour *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceFour *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceH file:Pythia8/ResonanceWidths.h line:280
struct PyCallBack_Pythia8_ResonanceH : public Pythia8::ResonanceH {
	using Pythia8::ResonanceH::ResonanceH;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceH *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceH *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceH *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceH *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceH *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceH *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceHchg file:Pythia8/ResonanceWidths.h line:333
struct PyCallBack_Pythia8_ResonanceHchg : public Pythia8::ResonanceHchg {
	using Pythia8::ResonanceHchg::ResonanceHchg;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchg *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchg *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchg *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchg *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchg *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchg *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceZprime file:Pythia8/ResonanceWidths.h line:362
struct PyCallBack_Pythia8_ResonanceZprime : public Pythia8::ResonanceZprime {
	using Pythia8::ResonanceZprime::ResonanceZprime;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZprime *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZprime *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZprime *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZprime *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZprime *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZprime *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceWprime file:Pythia8/ResonanceWidths.h line:395
struct PyCallBack_Pythia8_ResonanceWprime : public Pythia8::ResonanceWprime {
	using Pythia8::ResonanceWprime::ResonanceWprime;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWprime *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWprime *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWprime *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWprime *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWprime *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWprime *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceRhorizontal file:Pythia8/ResonanceWidths.h line:423
struct PyCallBack_Pythia8_ResonanceRhorizontal : public Pythia8::ResonanceRhorizontal {
	using Pythia8::ResonanceRhorizontal::ResonanceRhorizontal;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceRhorizontal *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceRhorizontal *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceRhorizontal *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceRhorizontal *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceRhorizontal *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceRhorizontal *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceExcited file:Pythia8/ResonanceWidths.h line:450
struct PyCallBack_Pythia8_ResonanceExcited : public Pythia8::ResonanceExcited {
	using Pythia8::ResonanceExcited::ResonanceExcited;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceExcited *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceExcited *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceExcited *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceExcited *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceExcited *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceExcited *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceGraviton file:Pythia8/ResonanceWidths.h line:478
struct PyCallBack_Pythia8_ResonanceGraviton : public Pythia8::ResonanceGraviton {
	using Pythia8::ResonanceGraviton::ResonanceGraviton;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGraviton *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGraviton *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGraviton *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGraviton *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGraviton *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceGraviton *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceKKgluon file:Pythia8/ResonanceWidths.h line:510
struct PyCallBack_Pythia8_ResonanceKKgluon : public Pythia8::ResonanceKKgluon {
	using Pythia8::ResonanceKKgluon::ResonanceKKgluon;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceKKgluon *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceKKgluon *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceKKgluon *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceKKgluon *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceKKgluon *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceKKgluon *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceLeptoquark file:Pythia8/ResonanceWidths.h line:546
struct PyCallBack_Pythia8_ResonanceLeptoquark : public Pythia8::ResonanceLeptoquark {
	using Pythia8::ResonanceLeptoquark::ResonanceLeptoquark;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceLeptoquark *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceLeptoquark *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceLeptoquark *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceLeptoquark *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceLeptoquark *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceLeptoquark *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceNuRight file:Pythia8/ResonanceWidths.h line:573
struct PyCallBack_Pythia8_ResonanceNuRight : public Pythia8::ResonanceNuRight {
	using Pythia8::ResonanceNuRight::ResonanceNuRight;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceNuRight *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceNuRight *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceNuRight *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceNuRight *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceNuRight *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceNuRight *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceZRight file:Pythia8/ResonanceWidths.h line:600
struct PyCallBack_Pythia8_ResonanceZRight : public Pythia8::ResonanceZRight {
	using Pythia8::ResonanceZRight::ResonanceZRight;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZRight *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZRight *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZRight *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZRight *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZRight *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceZRight *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceWRight file:Pythia8/ResonanceWidths.h line:627
struct PyCallBack_Pythia8_ResonanceWRight : public Pythia8::ResonanceWRight {
	using Pythia8::ResonanceWRight::ResonanceWRight;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWRight *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWRight *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWRight *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWRight *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWRight *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWRight *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceHchgchgLeft file:Pythia8/ResonanceWidths.h line:654
struct PyCallBack_Pythia8_ResonanceHchgchgLeft : public Pythia8::ResonanceHchgchgLeft {
	using Pythia8::ResonanceHchgchgLeft::ResonanceHchgchgLeft;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgLeft *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgLeft *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgLeft *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgLeft *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgLeft *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgLeft *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

// Pythia8::ResonanceHchgchgRight file:Pythia8/ResonanceWidths.h line:682
struct PyCallBack_Pythia8_ResonanceHchgchgRight : public Pythia8::ResonanceHchgchgRight {
	using Pythia8::ResonanceHchgchgRight::ResonanceHchgchgRight;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgRight *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::init(a0);
	}
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgRight *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::initConstants();
	}
	bool initBSM() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgRight *>(this), "initBSM");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::initBSM();
	}
	bool allowCalc() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgRight *>(this), "allowCalc");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceWidths::allowCalc();
	}
	void calcPreFac(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgRight *>(this), "calcPreFac");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcPreFac(a0);
	}
	void calcWidth(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceHchgchgRight *>(this), "calcWidth");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ResonanceWidths::calcWidth(a0);
	}
};

void bind_Pythia8_ResonanceWidths(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::ResonanceGeneric file:Pythia8/ResonanceWidths.h line:152
		pybind11::class_<Pythia8::ResonanceGeneric, std::shared_ptr<Pythia8::ResonanceGeneric>, PyCallBack_Pythia8_ResonanceGeneric, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceGeneric", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("allowCalc", (bool (Pythia8::ResonanceGeneric::*)()) &Pythia8::ResonanceGeneric::allowCalc, "C++: Pythia8::ResonanceGeneric::allowCalc() --> bool");
		cl.def("assign", (class Pythia8::ResonanceGeneric & (Pythia8::ResonanceGeneric::*)(const class Pythia8::ResonanceGeneric &)) &Pythia8::ResonanceGeneric::operator=, "C++: Pythia8::ResonanceGeneric::operator=(const class Pythia8::ResonanceGeneric &) --> class Pythia8::ResonanceGeneric &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceGmZ file:Pythia8/ResonanceWidths.h line:168
		pybind11::class_<Pythia8::ResonanceGmZ, std::shared_ptr<Pythia8::ResonanceGmZ>, PyCallBack_Pythia8_ResonanceGmZ, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceGmZ", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceGmZ & (Pythia8::ResonanceGmZ::*)(const class Pythia8::ResonanceGmZ &)) &Pythia8::ResonanceGmZ::operator=, "C++: Pythia8::ResonanceGmZ::operator=(const class Pythia8::ResonanceGmZ &) --> class Pythia8::ResonanceGmZ &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceW file:Pythia8/ResonanceWidths.h line:197
		pybind11::class_<Pythia8::ResonanceW, std::shared_ptr<Pythia8::ResonanceW>, PyCallBack_Pythia8_ResonanceW, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceW", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceW & (Pythia8::ResonanceW::*)(const class Pythia8::ResonanceW &)) &Pythia8::ResonanceW::operator=, "C++: Pythia8::ResonanceW::operator=(const class Pythia8::ResonanceW &) --> class Pythia8::ResonanceW &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceTop file:Pythia8/ResonanceWidths.h line:224
		pybind11::class_<Pythia8::ResonanceTop, std::shared_ptr<Pythia8::ResonanceTop>, PyCallBack_Pythia8_ResonanceTop, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceTop", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceTop & (Pythia8::ResonanceTop::*)(const class Pythia8::ResonanceTop &)) &Pythia8::ResonanceTop::operator=, "C++: Pythia8::ResonanceTop::operator=(const class Pythia8::ResonanceTop &) --> class Pythia8::ResonanceTop &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceFour file:Pythia8/ResonanceWidths.h line:252
		pybind11::class_<Pythia8::ResonanceFour, std::shared_ptr<Pythia8::ResonanceFour>, PyCallBack_Pythia8_ResonanceFour, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceFour", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceFour & (Pythia8::ResonanceFour::*)(const class Pythia8::ResonanceFour &)) &Pythia8::ResonanceFour::operator=, "C++: Pythia8::ResonanceFour::operator=(const class Pythia8::ResonanceFour &) --> class Pythia8::ResonanceFour &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceH file:Pythia8/ResonanceWidths.h line:280
		pybind11::class_<Pythia8::ResonanceH, std::shared_ptr<Pythia8::ResonanceH>, PyCallBack_Pythia8_ResonanceH, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceH", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int, int>(), pybind11::arg("higgsTypeIn"), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceH & (Pythia8::ResonanceH::*)(const class Pythia8::ResonanceH &)) &Pythia8::ResonanceH::operator=, "C++: Pythia8::ResonanceH::operator=(const class Pythia8::ResonanceH &) --> class Pythia8::ResonanceH &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceHchg file:Pythia8/ResonanceWidths.h line:333
		pybind11::class_<Pythia8::ResonanceHchg, std::shared_ptr<Pythia8::ResonanceHchg>, PyCallBack_Pythia8_ResonanceHchg, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceHchg", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceHchg & (Pythia8::ResonanceHchg::*)(const class Pythia8::ResonanceHchg &)) &Pythia8::ResonanceHchg::operator=, "C++: Pythia8::ResonanceHchg::operator=(const class Pythia8::ResonanceHchg &) --> class Pythia8::ResonanceHchg &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceZprime file:Pythia8/ResonanceWidths.h line:362
		pybind11::class_<Pythia8::ResonanceZprime, std::shared_ptr<Pythia8::ResonanceZprime>, PyCallBack_Pythia8_ResonanceZprime, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceZprime", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceZprime & (Pythia8::ResonanceZprime::*)(const class Pythia8::ResonanceZprime &)) &Pythia8::ResonanceZprime::operator=, "C++: Pythia8::ResonanceZprime::operator=(const class Pythia8::ResonanceZprime &) --> class Pythia8::ResonanceZprime &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceWprime file:Pythia8/ResonanceWidths.h line:395
		pybind11::class_<Pythia8::ResonanceWprime, std::shared_ptr<Pythia8::ResonanceWprime>, PyCallBack_Pythia8_ResonanceWprime, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceWprime", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceWprime & (Pythia8::ResonanceWprime::*)(const class Pythia8::ResonanceWprime &)) &Pythia8::ResonanceWprime::operator=, "C++: Pythia8::ResonanceWprime::operator=(const class Pythia8::ResonanceWprime &) --> class Pythia8::ResonanceWprime &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceRhorizontal file:Pythia8/ResonanceWidths.h line:423
		pybind11::class_<Pythia8::ResonanceRhorizontal, std::shared_ptr<Pythia8::ResonanceRhorizontal>, PyCallBack_Pythia8_ResonanceRhorizontal, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceRhorizontal", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceRhorizontal & (Pythia8::ResonanceRhorizontal::*)(const class Pythia8::ResonanceRhorizontal &)) &Pythia8::ResonanceRhorizontal::operator=, "C++: Pythia8::ResonanceRhorizontal::operator=(const class Pythia8::ResonanceRhorizontal &) --> class Pythia8::ResonanceRhorizontal &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceExcited file:Pythia8/ResonanceWidths.h line:450
		pybind11::class_<Pythia8::ResonanceExcited, std::shared_ptr<Pythia8::ResonanceExcited>, PyCallBack_Pythia8_ResonanceExcited, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceExcited", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceExcited & (Pythia8::ResonanceExcited::*)(const class Pythia8::ResonanceExcited &)) &Pythia8::ResonanceExcited::operator=, "C++: Pythia8::ResonanceExcited::operator=(const class Pythia8::ResonanceExcited &) --> class Pythia8::ResonanceExcited &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceGraviton file:Pythia8/ResonanceWidths.h line:478
		pybind11::class_<Pythia8::ResonanceGraviton, std::shared_ptr<Pythia8::ResonanceGraviton>, PyCallBack_Pythia8_ResonanceGraviton, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceGraviton", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceGraviton & (Pythia8::ResonanceGraviton::*)(const class Pythia8::ResonanceGraviton &)) &Pythia8::ResonanceGraviton::operator=, "C++: Pythia8::ResonanceGraviton::operator=(const class Pythia8::ResonanceGraviton &) --> class Pythia8::ResonanceGraviton &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceKKgluon file:Pythia8/ResonanceWidths.h line:510
		pybind11::class_<Pythia8::ResonanceKKgluon, std::shared_ptr<Pythia8::ResonanceKKgluon>, PyCallBack_Pythia8_ResonanceKKgluon, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceKKgluon", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceKKgluon & (Pythia8::ResonanceKKgluon::*)(const class Pythia8::ResonanceKKgluon &)) &Pythia8::ResonanceKKgluon::operator=, "C++: Pythia8::ResonanceKKgluon::operator=(const class Pythia8::ResonanceKKgluon &) --> class Pythia8::ResonanceKKgluon &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceLeptoquark file:Pythia8/ResonanceWidths.h line:546
		pybind11::class_<Pythia8::ResonanceLeptoquark, std::shared_ptr<Pythia8::ResonanceLeptoquark>, PyCallBack_Pythia8_ResonanceLeptoquark, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceLeptoquark", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceLeptoquark & (Pythia8::ResonanceLeptoquark::*)(const class Pythia8::ResonanceLeptoquark &)) &Pythia8::ResonanceLeptoquark::operator=, "C++: Pythia8::ResonanceLeptoquark::operator=(const class Pythia8::ResonanceLeptoquark &) --> class Pythia8::ResonanceLeptoquark &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceNuRight file:Pythia8/ResonanceWidths.h line:573
		pybind11::class_<Pythia8::ResonanceNuRight, std::shared_ptr<Pythia8::ResonanceNuRight>, PyCallBack_Pythia8_ResonanceNuRight, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceNuRight", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceNuRight & (Pythia8::ResonanceNuRight::*)(const class Pythia8::ResonanceNuRight &)) &Pythia8::ResonanceNuRight::operator=, "C++: Pythia8::ResonanceNuRight::operator=(const class Pythia8::ResonanceNuRight &) --> class Pythia8::ResonanceNuRight &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceZRight file:Pythia8/ResonanceWidths.h line:600
		pybind11::class_<Pythia8::ResonanceZRight, std::shared_ptr<Pythia8::ResonanceZRight>, PyCallBack_Pythia8_ResonanceZRight, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceZRight", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceZRight & (Pythia8::ResonanceZRight::*)(const class Pythia8::ResonanceZRight &)) &Pythia8::ResonanceZRight::operator=, "C++: Pythia8::ResonanceZRight::operator=(const class Pythia8::ResonanceZRight &) --> class Pythia8::ResonanceZRight &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceWRight file:Pythia8/ResonanceWidths.h line:627
		pybind11::class_<Pythia8::ResonanceWRight, std::shared_ptr<Pythia8::ResonanceWRight>, PyCallBack_Pythia8_ResonanceWRight, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceWRight", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceWRight & (Pythia8::ResonanceWRight::*)(const class Pythia8::ResonanceWRight &)) &Pythia8::ResonanceWRight::operator=, "C++: Pythia8::ResonanceWRight::operator=(const class Pythia8::ResonanceWRight &) --> class Pythia8::ResonanceWRight &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceHchgchgLeft file:Pythia8/ResonanceWidths.h line:654
		pybind11::class_<Pythia8::ResonanceHchgchgLeft, std::shared_ptr<Pythia8::ResonanceHchgchgLeft>, PyCallBack_Pythia8_ResonanceHchgchgLeft, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceHchgchgLeft", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceHchgchgLeft & (Pythia8::ResonanceHchgchgLeft::*)(const class Pythia8::ResonanceHchgchgLeft &)) &Pythia8::ResonanceHchgchgLeft::operator=, "C++: Pythia8::ResonanceHchgchgLeft::operator=(const class Pythia8::ResonanceHchgchgLeft &) --> class Pythia8::ResonanceHchgchgLeft &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResonanceHchgchgRight file:Pythia8/ResonanceWidths.h line:682
		pybind11::class_<Pythia8::ResonanceHchgchgRight, std::shared_ptr<Pythia8::ResonanceHchgchgRight>, PyCallBack_Pythia8_ResonanceHchgchgRight, Pythia8::ResonanceWidths> cl(M("Pythia8"), "ResonanceHchgchgRight", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int>(), pybind11::arg("idResIn") );

		cl.def("assign", (class Pythia8::ResonanceHchgchgRight & (Pythia8::ResonanceHchgchgRight::*)(const class Pythia8::ResonanceHchgchgRight &)) &Pythia8::ResonanceHchgchgRight::operator=, "C++: Pythia8::ResonanceHchgchgRight::operator=(const class Pythia8::ResonanceHchgchgRight &) --> class Pythia8::ResonanceHchgchgRight &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LHblock file:Pythia8/SusyLesHouches.h line:26
		pybind11::class_<Pythia8::LHblock<std::string>, std::shared_ptr<Pythia8::LHblock<std::string>>> cl(M("Pythia8"), "LHblock_std_string_t", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHblock<std::string>(); } ) );
		cl.def( pybind11::init( [](Pythia8::LHblock<std::string> const &o){ return new Pythia8::LHblock<std::string>(o); } ) );
		cl.def_readwrite("entry", &Pythia8::LHblock<std::string>::entry);
		cl.def("exists", (bool (Pythia8::LHblock<std::string>::*)()) &Pythia8::LHblock<std::string >::exists, "C++: Pythia8::LHblock<std::basic_string<char> >::exists() --> bool");
		cl.def("clear", (void (Pythia8::LHblock<std::string>::*)()) &Pythia8::LHblock<std::string >::clear, "C++: Pythia8::LHblock<std::basic_string<char> >::clear() --> void");
		cl.def("set", (int (Pythia8::LHblock<std::string>::*)(int, std::string)) &Pythia8::LHblock<std::string >::set, "C++: Pythia8::LHblock<std::basic_string<char> >::set(int, std::string) --> int", pybind11::arg("iIn"), pybind11::arg("valIn"));
		cl.def("set", [](Pythia8::LHblock<std::string> &o, class std::basic_istringstream<char> & a0) -> int { return o.set(a0); }, "", pybind11::arg("linestream"));
		cl.def("set", (int (Pythia8::LHblock<std::string>::*)(class std::basic_istringstream<char> &, bool)) &Pythia8::LHblock<std::string >::set, "C++: Pythia8::LHblock<std::basic_string<char> >::set(class std::basic_istringstream<char> &, bool) --> int", pybind11::arg("linestream"), pybind11::arg("indexed"));
		cl.def("set", (int (Pythia8::LHblock<std::string>::*)(int, class std::basic_istringstream<char> &)) &Pythia8::LHblock<std::string >::set, "C++: Pythia8::LHblock<std::basic_string<char> >::set(int, class std::basic_istringstream<char> &) --> int", pybind11::arg("iIn"), pybind11::arg("linestream"));
		cl.def("set", (void (Pythia8::LHblock<std::string>::*)(std::string)) &Pythia8::LHblock<std::string >::set, "C++: Pythia8::LHblock<std::basic_string<char> >::set(std::string) --> void", pybind11::arg("valIn"));
		cl.def("exists", (bool (Pythia8::LHblock<std::string>::*)(int)) &Pythia8::LHblock<std::string >::exists, "C++: Pythia8::LHblock<std::basic_string<char> >::exists(int) --> bool", pybind11::arg("iIn"));
		cl.def("__call__", (std::string (Pythia8::LHblock<std::string>::*)()) &Pythia8::LHblock<std::string >::operator(), "C++: Pythia8::LHblock<std::basic_string<char> >::operator()() --> std::string");
		cl.def("__call__", (std::string (Pythia8::LHblock<std::string>::*)(int)) &Pythia8::LHblock<std::string >::operator(), "C++: Pythia8::LHblock<std::basic_string<char> >::operator()(int) --> std::string", pybind11::arg("iIn"));
		cl.def("size", (int (Pythia8::LHblock<std::string>::*)()) &Pythia8::LHblock<std::string >::size, "C++: Pythia8::LHblock<std::basic_string<char> >::size() --> int");
		cl.def("first", (int (Pythia8::LHblock<std::string>::*)()) &Pythia8::LHblock<std::string >::first, "C++: Pythia8::LHblock<std::basic_string<char> >::first() --> int");
		cl.def("next", (int (Pythia8::LHblock<std::string>::*)()) &Pythia8::LHblock<std::string >::next, "C++: Pythia8::LHblock<std::basic_string<char> >::next() --> int");
		cl.def("list", (void (Pythia8::LHblock<std::string>::*)()) &Pythia8::LHblock<std::string >::list, "C++: Pythia8::LHblock<std::basic_string<char> >::list() --> void");
		cl.def("setq", (void (Pythia8::LHblock<std::string>::*)(double)) &Pythia8::LHblock<std::string >::setq, "C++: Pythia8::LHblock<std::basic_string<char> >::setq(double) --> void", pybind11::arg("qIn"));
		cl.def("q", (double (Pythia8::LHblock<std::string>::*)()) &Pythia8::LHblock<std::string >::q, "C++: Pythia8::LHblock<std::basic_string<char> >::q() --> double");
		cl.def("assign", (class Pythia8::LHblock<std::string > & (Pythia8::LHblock<std::string>::*)(const class Pythia8::LHblock<std::string > &)) &Pythia8::LHblock<std::string >::operator=, "C++: Pythia8::LHblock<std::basic_string<char> >::operator=(const class Pythia8::LHblock<std::string > &) --> class Pythia8::LHblock<std::string > &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
