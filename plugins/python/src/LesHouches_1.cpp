#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/HIInfo.h>
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

// Pythia8::LHAupLHEF file:Pythia8/LesHouches.h line:346
struct PyCallBack_Pythia8_LHAupLHEF : public Pythia8::LHAupLHEF {
	using Pythia8::LHAupLHEF::LHAupLHEF;

	void newEventFile(const char * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupLHEF *>(this), "newEventFile");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LHAupLHEF::newEventFile(a0);
	}
	bool fileFound() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupLHEF *>(this), "fileFound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupLHEF::fileFound();
	}
	bool useExternal() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupLHEF *>(this), "useExternal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupLHEF::useExternal();
	}
	bool setInit() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupLHEF *>(this), "setInit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupLHEF::setInit();
	}
	bool setEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupLHEF *>(this), "setEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupLHEF::setEvent(a0);
	}
	bool skipEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupLHEF *>(this), "skipEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupLHEF::skipEvent(a0);
	}
	bool openLHEF(class std::basic_string<char> a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupLHEF *>(this), "openLHEF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::openLHEF(a0);
	}
	bool closeLHEF(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupLHEF *>(this), "closeLHEF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::closeLHEF(a0);
	}
};

// Pythia8::LHAupFromPYTHIA8 file:Pythia8/LesHouches.h line:484
struct PyCallBack_Pythia8_LHAupFromPYTHIA8 : public Pythia8::LHAupFromPYTHIA8 {
	using Pythia8::LHAupFromPYTHIA8::LHAupFromPYTHIA8;

	bool setInit() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupFromPYTHIA8 *>(this), "setInit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupFromPYTHIA8::setInit();
	}
	bool setEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupFromPYTHIA8 *>(this), "setEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupFromPYTHIA8::setEvent(a0);
	}
	void newEventFile(const char * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupFromPYTHIA8 *>(this), "newEventFile");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LHAup::newEventFile(a0);
	}
	bool fileFound() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupFromPYTHIA8 *>(this), "fileFound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::fileFound();
	}
	bool useExternal() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupFromPYTHIA8 *>(this), "useExternal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::useExternal();
	}
	bool skipEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupFromPYTHIA8 *>(this), "skipEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::skipEvent(a0);
	}
	bool openLHEF(class std::basic_string<char> a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupFromPYTHIA8 *>(this), "openLHEF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::openLHEF(a0);
	}
	bool closeLHEF(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupFromPYTHIA8 *>(this), "closeLHEF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::closeLHEF(a0);
	}
};

// Pythia8::LHEF3FromPythia8 file:Pythia8/LesHouches.h line:524
struct PyCallBack_Pythia8_LHEF3FromPythia8 : public Pythia8::LHEF3FromPythia8 {
	using Pythia8::LHEF3FromPythia8::LHEF3FromPythia8;

	bool setInit() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHEF3FromPythia8 *>(this), "setInit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHEF3FromPythia8::setInit();
	}
	bool setEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHEF3FromPythia8 *>(this), "setEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHEF3FromPythia8::setEvent(a0);
	}
	bool openLHEF(class std::basic_string<char> a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHEF3FromPythia8 *>(this), "openLHEF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHEF3FromPythia8::openLHEF(a0);
	}
	bool closeLHEF(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHEF3FromPythia8 *>(this), "closeLHEF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHEF3FromPythia8::closeLHEF(a0);
	}
	void newEventFile(const char * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHEF3FromPythia8 *>(this), "newEventFile");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LHAup::newEventFile(a0);
	}
	bool fileFound() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHEF3FromPythia8 *>(this), "fileFound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::fileFound();
	}
	bool useExternal() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHEF3FromPythia8 *>(this), "useExternal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::useExternal();
	}
	bool skipEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHEF3FromPythia8 *>(this), "skipEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAup::skipEvent(a0);
	}
};

// Pythia8::ResonanceWidths file:Pythia8/ResonanceWidths.h line:34
struct PyCallBack_Pythia8_ResonanceWidths : public Pythia8::ResonanceWidths {
	using Pythia8::ResonanceWidths::ResonanceWidths;

	bool init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWidths *>(this), "init");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWidths *>(this), "initConstants");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWidths *>(this), "initBSM");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWidths *>(this), "allowCalc");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWidths *>(this), "calcPreFac");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceWidths *>(this), "calcWidth");
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

void bind_Pythia8_LesHouches_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::LHAupLHEF file:Pythia8/LesHouches.h line:346
		pybind11::class_<Pythia8::LHAupLHEF, std::shared_ptr<Pythia8::LHAupLHEF>, PyCallBack_Pythia8_LHAupLHEF, Pythia8::LHAup> cl(M("Pythia8"), "LHAupLHEF", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](class Pythia8::Info * a0, class std::basic_istream<char> * a1, class std::basic_istream<char> * a2){ return new Pythia8::LHAupLHEF(a0, a1, a2); }, [](class Pythia8::Info * a0, class std::basic_istream<char> * a1, class std::basic_istream<char> * a2){ return new PyCallBack_Pythia8_LHAupLHEF(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init( [](class Pythia8::Info * a0, class std::basic_istream<char> * a1, class std::basic_istream<char> * a2, bool const & a3){ return new Pythia8::LHAupLHEF(a0, a1, a2, a3); }, [](class Pythia8::Info * a0, class std::basic_istream<char> * a1, class std::basic_istream<char> * a2, bool const & a3){ return new PyCallBack_Pythia8_LHAupLHEF(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init<class Pythia8::Info *, class std::basic_istream<char> *, class std::basic_istream<char> *, bool, bool>(), pybind11::arg("infoPtrIn"), pybind11::arg("isIn"), pybind11::arg("isHeadIn"), pybind11::arg("readHeadersIn"), pybind11::arg("setScalesFromLHEFIn") );

		cl.def( pybind11::init( [](class Pythia8::Info * a0, const char * a1){ return new Pythia8::LHAupLHEF(a0, a1); }, [](class Pythia8::Info * a0, const char * a1){ return new PyCallBack_Pythia8_LHAupLHEF(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](class Pythia8::Info * a0, const char * a1, const char * a2){ return new Pythia8::LHAupLHEF(a0, a1, a2); }, [](class Pythia8::Info * a0, const char * a1, const char * a2){ return new PyCallBack_Pythia8_LHAupLHEF(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init( [](class Pythia8::Info * a0, const char * a1, const char * a2, bool const & a3){ return new Pythia8::LHAupLHEF(a0, a1, a2, a3); }, [](class Pythia8::Info * a0, const char * a1, const char * a2, bool const & a3){ return new PyCallBack_Pythia8_LHAupLHEF(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init<class Pythia8::Info *, const char *, const char *, bool, bool>(), pybind11::arg("infoPtrIn"), pybind11::arg("filenameIn"), pybind11::arg("headerIn"), pybind11::arg("readHeadersIn"), pybind11::arg("setScalesFromLHEFIn") );

		cl.def("closeAllFiles", (void (Pythia8::LHAupLHEF::*)()) &Pythia8::LHAupLHEF::closeAllFiles, "C++: Pythia8::LHAupLHEF::closeAllFiles() --> void");
		cl.def("newEventFile", (void (Pythia8::LHAupLHEF::*)(const char *)) &Pythia8::LHAupLHEF::newEventFile, "C++: Pythia8::LHAupLHEF::newEventFile(const char *) --> void", pybind11::arg("filenameIn"));
		cl.def("fileFound", (bool (Pythia8::LHAupLHEF::*)()) &Pythia8::LHAupLHEF::fileFound, "C++: Pythia8::LHAupLHEF::fileFound() --> bool");
		cl.def("useExternal", (bool (Pythia8::LHAupLHEF::*)()) &Pythia8::LHAupLHEF::useExternal, "C++: Pythia8::LHAupLHEF::useExternal() --> bool");
		cl.def("setInit", (bool (Pythia8::LHAupLHEF::*)()) &Pythia8::LHAupLHEF::setInit, "C++: Pythia8::LHAupLHEF::setInit() --> bool");
		cl.def("setInitLHEF", (bool (Pythia8::LHAupLHEF::*)(class std::basic_istream<char> &, bool)) &Pythia8::LHAupLHEF::setInitLHEF, "C++: Pythia8::LHAupLHEF::setInitLHEF(class std::basic_istream<char> &, bool) --> bool", pybind11::arg("isIn"), pybind11::arg("readHead"));
		cl.def("setEvent", [](Pythia8::LHAupLHEF &o) -> bool { return o.setEvent(); }, "");
		cl.def("setEvent", (bool (Pythia8::LHAupLHEF::*)(int)) &Pythia8::LHAupLHEF::setEvent, "C++: Pythia8::LHAupLHEF::setEvent(int) --> bool", pybind11::arg(""));
		cl.def("skipEvent", (bool (Pythia8::LHAupLHEF::*)(int)) &Pythia8::LHAupLHEF::skipEvent, "C++: Pythia8::LHAupLHEF::skipEvent(int) --> bool", pybind11::arg("nSkip"));
		cl.def("setNewEventLHEF", (bool (Pythia8::LHAupLHEF::*)()) &Pythia8::LHAupLHEF::setNewEventLHEF, "C++: Pythia8::LHAupLHEF::setNewEventLHEF() --> bool");
		cl.def("updateSigma", (bool (Pythia8::LHAupLHEF::*)()) &Pythia8::LHAupLHEF::updateSigma, "C++: Pythia8::LHAupLHEF::updateSigma() --> bool");
		cl.def("getLine", [](Pythia8::LHAupLHEF &o, class std::basic_string<char> & a0) -> bool { return o.getLine(a0); }, "", pybind11::arg("line"));
		cl.def("getLine", (bool (Pythia8::LHAupLHEF::*)(std::string &, bool)) &Pythia8::LHAupLHEF::getLine, "C++: Pythia8::LHAupLHEF::getLine(std::string &, bool) --> bool", pybind11::arg("line"), pybind11::arg("header"));
	}
	{ // Pythia8::LHAupFromPYTHIA8 file:Pythia8/LesHouches.h line:484
		pybind11::class_<Pythia8::LHAupFromPYTHIA8, std::shared_ptr<Pythia8::LHAupFromPYTHIA8>, PyCallBack_Pythia8_LHAupFromPYTHIA8, Pythia8::LHAup> cl(M("Pythia8"), "LHAupFromPYTHIA8", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Event *, const class Pythia8::Info *>(), pybind11::arg("processPtrIn"), pybind11::arg("infoPtrIn") );

		cl.def("setInit", (bool (Pythia8::LHAupFromPYTHIA8::*)()) &Pythia8::LHAupFromPYTHIA8::setInit, "C++: Pythia8::LHAupFromPYTHIA8::setInit() --> bool");
		cl.def("setEvent", [](Pythia8::LHAupFromPYTHIA8 &o) -> bool { return o.setEvent(); }, "");
		cl.def("setEvent", (bool (Pythia8::LHAupFromPYTHIA8::*)(int)) &Pythia8::LHAupFromPYTHIA8::setEvent, "C++: Pythia8::LHAupFromPYTHIA8::setEvent(int) --> bool", pybind11::arg(""));
		cl.def("updateSigma", (bool (Pythia8::LHAupFromPYTHIA8::*)()) &Pythia8::LHAupFromPYTHIA8::updateSigma, "C++: Pythia8::LHAupFromPYTHIA8::updateSigma() --> bool");
	}
	{ // Pythia8::LHEF3FromPythia8 file:Pythia8/LesHouches.h line:524
		pybind11::class_<Pythia8::LHEF3FromPythia8, std::shared_ptr<Pythia8::LHEF3FromPythia8>, PyCallBack_Pythia8_LHEF3FromPythia8, Pythia8::LHAup> cl(M("Pythia8"), "LHEF3FromPythia8", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](class Pythia8::Pythia * a0){ return new Pythia8::LHEF3FromPythia8(a0); }, [](class Pythia8::Pythia * a0){ return new PyCallBack_Pythia8_LHEF3FromPythia8(a0); } ), "doc");
		cl.def( pybind11::init( [](class Pythia8::Pythia * a0, int const & a1){ return new Pythia8::LHEF3FromPythia8(a0, a1); }, [](class Pythia8::Pythia * a0, int const & a1){ return new PyCallBack_Pythia8_LHEF3FromPythia8(a0, a1); } ), "doc");
		cl.def( pybind11::init<class Pythia8::Pythia *, int, bool>(), pybind11::arg("pythiaPtrIn"), pybind11::arg("pDigitsIn"), pybind11::arg("writeToFileIn") );

		cl.def( pybind11::init( [](class Pythia8::Event * a0, const class Pythia8::Info * a1){ return new Pythia8::LHEF3FromPythia8(a0, a1); }, [](class Pythia8::Event * a0, const class Pythia8::Info * a1){ return new PyCallBack_Pythia8_LHEF3FromPythia8(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](class Pythia8::Event * a0, const class Pythia8::Info * a1, int const & a2){ return new Pythia8::LHEF3FromPythia8(a0, a1, a2); }, [](class Pythia8::Event * a0, const class Pythia8::Info * a1, int const & a2){ return new PyCallBack_Pythia8_LHEF3FromPythia8(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<class Pythia8::Event *, const class Pythia8::Info *, int, bool>(), pybind11::arg("eventPtrIn"), pybind11::arg("infoPtrIn"), pybind11::arg("pDigitsIn"), pybind11::arg("writeToFileIn") );

		cl.def_readwrite("heprup", &Pythia8::LHEF3FromPythia8::heprup);
		cl.def_readwrite("hepeup", &Pythia8::LHEF3FromPythia8::hepeup);
		cl.def("setInit", (bool (Pythia8::LHEF3FromPythia8::*)()) &Pythia8::LHEF3FromPythia8::setInit, "C++: Pythia8::LHEF3FromPythia8::setInit() --> bool");
		cl.def("setEventPtr", (void (Pythia8::LHEF3FromPythia8::*)(class Pythia8::Event *)) &Pythia8::LHEF3FromPythia8::setEventPtr, "C++: Pythia8::LHEF3FromPythia8::setEventPtr(class Pythia8::Event *) --> void", pybind11::arg("evPtr"));
		cl.def("setEvent", [](Pythia8::LHEF3FromPythia8 &o) -> bool { return o.setEvent(); }, "");
		cl.def("setEvent", (bool (Pythia8::LHEF3FromPythia8::*)(int)) &Pythia8::LHEF3FromPythia8::setEvent, "C++: Pythia8::LHEF3FromPythia8::setEvent(int) --> bool", pybind11::arg(""));
		cl.def("getEventString", (std::string (Pythia8::LHEF3FromPythia8::*)()) &Pythia8::LHEF3FromPythia8::getEventString, "C++: Pythia8::LHEF3FromPythia8::getEventString() --> std::string");
		cl.def("openLHEF", (bool (Pythia8::LHEF3FromPythia8::*)(std::string)) &Pythia8::LHEF3FromPythia8::openLHEF, "C++: Pythia8::LHEF3FromPythia8::openLHEF(std::string) --> bool", pybind11::arg("fileNameIn"));
		cl.def("closeLHEF", [](Pythia8::LHEF3FromPythia8 &o) -> bool { return o.closeLHEF(); }, "");
		cl.def("closeLHEF", (bool (Pythia8::LHEF3FromPythia8::*)(bool)) &Pythia8::LHEF3FromPythia8::closeLHEF, "C++: Pythia8::LHEF3FromPythia8::closeLHEF(bool) --> bool", pybind11::arg("updateInit"));
	}
	{ // Pythia8::ResonanceWidths file:Pythia8/ResonanceWidths.h line:34
		pybind11::class_<Pythia8::ResonanceWidths, std::shared_ptr<Pythia8::ResonanceWidths>, PyCallBack_Pythia8_ResonanceWidths> cl(M("Pythia8"), "ResonanceWidths", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ResonanceWidths(); }, [](){ return new PyCallBack_Pythia8_ResonanceWidths(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_ResonanceWidths const &o){ return new PyCallBack_Pythia8_ResonanceWidths(o); } ) );
		cl.def( pybind11::init( [](Pythia8::ResonanceWidths const &o){ return new Pythia8::ResonanceWidths(o); } ) );
		cl.def_readwrite("idRes", &Pythia8::ResonanceWidths::idRes);
		cl.def_readwrite("hasAntiRes", &Pythia8::ResonanceWidths::hasAntiRes);
		cl.def_readwrite("doForceWidth", &Pythia8::ResonanceWidths::doForceWidth);
		cl.def_readwrite("isGeneric", &Pythia8::ResonanceWidths::isGeneric);
		cl.def_readwrite("allowCalcWidth", &Pythia8::ResonanceWidths::allowCalcWidth);
		cl.def_readwrite("minWidth", &Pythia8::ResonanceWidths::minWidth);
		cl.def_readwrite("minThreshold", &Pythia8::ResonanceWidths::minThreshold);
		cl.def_readwrite("mRes", &Pythia8::ResonanceWidths::mRes);
		cl.def_readwrite("GammaRes", &Pythia8::ResonanceWidths::GammaRes);
		cl.def_readwrite("m2Res", &Pythia8::ResonanceWidths::m2Res);
		cl.def_readwrite("GamMRat", &Pythia8::ResonanceWidths::GamMRat);
		cl.def_readwrite("openPos", &Pythia8::ResonanceWidths::openPos);
		cl.def_readwrite("openNeg", &Pythia8::ResonanceWidths::openNeg);
		cl.def_readwrite("forceFactor", &Pythia8::ResonanceWidths::forceFactor);
		cl.def_readwrite("iChannel", &Pythia8::ResonanceWidths::iChannel);
		cl.def_readwrite("onMode", &Pythia8::ResonanceWidths::onMode);
		cl.def_readwrite("meMode", &Pythia8::ResonanceWidths::meMode);
		cl.def_readwrite("mult", &Pythia8::ResonanceWidths::mult);
		cl.def_readwrite("id1", &Pythia8::ResonanceWidths::id1);
		cl.def_readwrite("id2", &Pythia8::ResonanceWidths::id2);
		cl.def_readwrite("id3", &Pythia8::ResonanceWidths::id3);
		cl.def_readwrite("id1Abs", &Pythia8::ResonanceWidths::id1Abs);
		cl.def_readwrite("id2Abs", &Pythia8::ResonanceWidths::id2Abs);
		cl.def_readwrite("id3Abs", &Pythia8::ResonanceWidths::id3Abs);
		cl.def_readwrite("idInFlav", &Pythia8::ResonanceWidths::idInFlav);
		cl.def_readwrite("widNow", &Pythia8::ResonanceWidths::widNow);
		cl.def_readwrite("mHat", &Pythia8::ResonanceWidths::mHat);
		cl.def_readwrite("mf1", &Pythia8::ResonanceWidths::mf1);
		cl.def_readwrite("mf2", &Pythia8::ResonanceWidths::mf2);
		cl.def_readwrite("mf3", &Pythia8::ResonanceWidths::mf3);
		cl.def_readwrite("mr1", &Pythia8::ResonanceWidths::mr1);
		cl.def_readwrite("mr2", &Pythia8::ResonanceWidths::mr2);
		cl.def_readwrite("mr3", &Pythia8::ResonanceWidths::mr3);
		cl.def_readwrite("ps", &Pythia8::ResonanceWidths::ps);
		cl.def_readwrite("kinFac", &Pythia8::ResonanceWidths::kinFac);
		cl.def_readwrite("alpEM", &Pythia8::ResonanceWidths::alpEM);
		cl.def_readwrite("alpS", &Pythia8::ResonanceWidths::alpS);
		cl.def_readwrite("colQ", &Pythia8::ResonanceWidths::colQ);
		cl.def_readwrite("preFac", &Pythia8::ResonanceWidths::preFac);
		cl.def_readwrite("particlePtr", &Pythia8::ResonanceWidths::particlePtr);
		cl.def("initBasic", [](Pythia8::ResonanceWidths &o, int const & a0) -> void { return o.initBasic(a0); }, "", pybind11::arg("idResIn"));
		cl.def("initBasic", (void (Pythia8::ResonanceWidths::*)(int, bool)) &Pythia8::ResonanceWidths::initBasic, "C++: Pythia8::ResonanceWidths::initBasic(int, bool) --> void", pybind11::arg("idResIn"), pybind11::arg("isGenericIn"));
		cl.def("init", (bool (Pythia8::ResonanceWidths::*)(class Pythia8::Info *)) &Pythia8::ResonanceWidths::init, "C++: Pythia8::ResonanceWidths::init(class Pythia8::Info *) --> bool", pybind11::arg("infoPtrIn"));
		cl.def("id", (int (Pythia8::ResonanceWidths::*)() const) &Pythia8::ResonanceWidths::id, "C++: Pythia8::ResonanceWidths::id() const --> int");
		cl.def("width", [](Pythia8::ResonanceWidths &o, int const & a0, double const & a1) -> double { return o.width(a0, a1); }, "", pybind11::arg("idSgn"), pybind11::arg("mHatIn"));
		cl.def("width", [](Pythia8::ResonanceWidths &o, int const & a0, double const & a1, int const & a2) -> double { return o.width(a0, a1, a2); }, "", pybind11::arg("idSgn"), pybind11::arg("mHatIn"), pybind11::arg("idInFlavIn"));
		cl.def("width", [](Pythia8::ResonanceWidths &o, int const & a0, double const & a1, int const & a2, bool const & a3) -> double { return o.width(a0, a1, a2, a3); }, "", pybind11::arg("idSgn"), pybind11::arg("mHatIn"), pybind11::arg("idInFlavIn"), pybind11::arg("openOnly"));
		cl.def("width", [](Pythia8::ResonanceWidths &o, int const & a0, double const & a1, int const & a2, bool const & a3, bool const & a4) -> double { return o.width(a0, a1, a2, a3, a4); }, "", pybind11::arg("idSgn"), pybind11::arg("mHatIn"), pybind11::arg("idInFlavIn"), pybind11::arg("openOnly"), pybind11::arg("setBR"));
		cl.def("width", [](Pythia8::ResonanceWidths &o, int const & a0, double const & a1, int const & a2, bool const & a3, bool const & a4, int const & a5) -> double { return o.width(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("idSgn"), pybind11::arg("mHatIn"), pybind11::arg("idInFlavIn"), pybind11::arg("openOnly"), pybind11::arg("setBR"), pybind11::arg("idOutFlav1"));
		cl.def("width", (double (Pythia8::ResonanceWidths::*)(int, double, int, bool, bool, int, int)) &Pythia8::ResonanceWidths::width, "C++: Pythia8::ResonanceWidths::width(int, double, int, bool, bool, int, int) --> double", pybind11::arg("idSgn"), pybind11::arg("mHatIn"), pybind11::arg("idInFlavIn"), pybind11::arg("openOnly"), pybind11::arg("setBR"), pybind11::arg("idOutFlav1"), pybind11::arg("idOutFlav2"));
		cl.def("widthOpen", [](Pythia8::ResonanceWidths &o, int const & a0, double const & a1) -> double { return o.widthOpen(a0, a1); }, "", pybind11::arg("idSgn"), pybind11::arg("mHatIn"));
		cl.def("widthOpen", (double (Pythia8::ResonanceWidths::*)(int, double, int)) &Pythia8::ResonanceWidths::widthOpen, "C++: Pythia8::ResonanceWidths::widthOpen(int, double, int) --> double", pybind11::arg("idSgn"), pybind11::arg("mHatIn"), pybind11::arg("idIn"));
		cl.def("widthStore", [](Pythia8::ResonanceWidths &o, int const & a0, double const & a1) -> double { return o.widthStore(a0, a1); }, "", pybind11::arg("idSgn"), pybind11::arg("mHatIn"));
		cl.def("widthStore", (double (Pythia8::ResonanceWidths::*)(int, double, int)) &Pythia8::ResonanceWidths::widthStore, "C++: Pythia8::ResonanceWidths::widthStore(int, double, int) --> double", pybind11::arg("idSgn"), pybind11::arg("mHatIn"), pybind11::arg("idIn"));
		cl.def("openFrac", (double (Pythia8::ResonanceWidths::*)(int)) &Pythia8::ResonanceWidths::openFrac, "C++: Pythia8::ResonanceWidths::openFrac(int) --> double", pybind11::arg("idSgn"));
		cl.def("widthRescaleFactor", (double (Pythia8::ResonanceWidths::*)()) &Pythia8::ResonanceWidths::widthRescaleFactor, "C++: Pythia8::ResonanceWidths::widthRescaleFactor() --> double");
		cl.def("widthChan", (double (Pythia8::ResonanceWidths::*)(double, int, int)) &Pythia8::ResonanceWidths::widthChan, "C++: Pythia8::ResonanceWidths::widthChan(double, int, int) --> double", pybind11::arg("mHatIn"), pybind11::arg("idOutFlav1"), pybind11::arg("idOutFlav2"));
		cl.def("initConstants", (void (Pythia8::ResonanceWidths::*)()) &Pythia8::ResonanceWidths::initConstants, "C++: Pythia8::ResonanceWidths::initConstants() --> void");
		cl.def("initBSM", (bool (Pythia8::ResonanceWidths::*)()) &Pythia8::ResonanceWidths::initBSM, "C++: Pythia8::ResonanceWidths::initBSM() --> bool");
		cl.def("allowCalc", (bool (Pythia8::ResonanceWidths::*)()) &Pythia8::ResonanceWidths::allowCalc, "C++: Pythia8::ResonanceWidths::allowCalc() --> bool");
		cl.def("calcPreFac", [](Pythia8::ResonanceWidths &o) -> void { return o.calcPreFac(); }, "");
		cl.def("calcPreFac", (void (Pythia8::ResonanceWidths::*)(bool)) &Pythia8::ResonanceWidths::calcPreFac, "C++: Pythia8::ResonanceWidths::calcPreFac(bool) --> void", pybind11::arg(""));
		cl.def("calcWidth", [](Pythia8::ResonanceWidths &o) -> void { return o.calcWidth(); }, "");
		cl.def("calcWidth", (void (Pythia8::ResonanceWidths::*)(bool)) &Pythia8::ResonanceWidths::calcWidth, "C++: Pythia8::ResonanceWidths::calcWidth(bool) --> void", pybind11::arg(""));
		cl.def("numInt1BW", [](Pythia8::ResonanceWidths &o, double const & a0, double const & a1, double const & a2, double const & a3, double const & a4) -> double { return o.numInt1BW(a0, a1, a2, a3, a4); }, "", pybind11::arg("mHatIn"), pybind11::arg("m1"), pybind11::arg("Gamma1"), pybind11::arg("mMin1"), pybind11::arg("m2"));
		cl.def("numInt1BW", (double (Pythia8::ResonanceWidths::*)(double, double, double, double, double, int)) &Pythia8::ResonanceWidths::numInt1BW, "C++: Pythia8::ResonanceWidths::numInt1BW(double, double, double, double, double, int) --> double", pybind11::arg("mHatIn"), pybind11::arg("m1"), pybind11::arg("Gamma1"), pybind11::arg("mMin1"), pybind11::arg("m2"), pybind11::arg("psMode"));
		cl.def("numInt2BW", [](Pythia8::ResonanceWidths &o, double const & a0, double const & a1, double const & a2, double const & a3, double const & a4, double const & a5, double const & a6) -> double { return o.numInt2BW(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("mHatIn"), pybind11::arg("m1"), pybind11::arg("Gamma1"), pybind11::arg("mMin1"), pybind11::arg("m2"), pybind11::arg("Gamma2"), pybind11::arg("mMin2"));
		cl.def("numInt2BW", (double (Pythia8::ResonanceWidths::*)(double, double, double, double, double, double, double, int)) &Pythia8::ResonanceWidths::numInt2BW, "C++: Pythia8::ResonanceWidths::numInt2BW(double, double, double, double, double, double, double, int) --> double", pybind11::arg("mHatIn"), pybind11::arg("m1"), pybind11::arg("Gamma1"), pybind11::arg("mMin1"), pybind11::arg("m2"), pybind11::arg("Gamma2"), pybind11::arg("mMin2"), pybind11::arg("psMode"));
		cl.def("assign", (class Pythia8::ResonanceWidths & (Pythia8::ResonanceWidths::*)(const class Pythia8::ResonanceWidths &)) &Pythia8::ResonanceWidths::operator=, "C++: Pythia8::ResonanceWidths::operator=(const class Pythia8::ResonanceWidths &) --> class Pythia8::ResonanceWidths &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
