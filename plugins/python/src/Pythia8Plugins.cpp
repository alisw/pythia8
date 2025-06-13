#include <Pythia8/Analysis.h>
#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/GammaKinematics.h>
#include <Pythia8/HIInfo.h>
#include <Pythia8/HadronLevel.h>
#include <Pythia8/HeavyIons.h>
#include <Pythia8/Info.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/Merging.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/NucleonExcitations.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ParticleDecays.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/RHadrons.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/SLHAinterface.h>
#include <Pythia8/Settings.h>
#include <Pythia8/ShowerModel.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaProcess.h>
#include <Pythia8/StringFragmentation.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/TimeShower.h>
#include <Pythia8/UserHooks.h>
#include <Pythia8Plugins/ColourReconnectionHooks.h>
#include <Pythia8Plugins/CombineMatchingInput.h>
#include <Pythia8Plugins/GeneratorInput.h>
#include <Pythia8Plugins/JetMatching.h>
#include <Pythia8Plugins/aMCatNLOHooks.h>
#include <functional>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
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

// Pythia8::amcnlo_unitarised_interface file:Pythia8Plugins/aMCatNLOHooks.h line:28
struct PyCallBack_Pythia8_amcnlo_unitarised_interface : public Pythia8::amcnlo_unitarised_interface {
	using Pythia8::amcnlo_unitarised_interface::amcnlo_unitarised_interface;

	bool canVetoProcessLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return amcnlo_unitarised_interface::canVetoProcessLevel();
	}
	bool doVetoProcessLevel(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return amcnlo_unitarised_interface::doVetoProcessLevel(a0);
	}
	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "initAfterBeams");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "biasedSelectionWeight");
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
	bool canSetLowEnergySigma(int a0, int a1) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "numberVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::amcnlo_unitarised_interface *>(this), "onStat");
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

// Pythia8::MBReconUserHooks file:Pythia8Plugins/ColourReconnectionHooks.h line:33
struct PyCallBack_Pythia8_MBReconUserHooks : public Pythia8::MBReconUserHooks {
	using Pythia8::MBReconUserHooks::MBReconUserHooks;

	bool canReconnectResonanceSystems() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canReconnectResonanceSystems");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return MBReconUserHooks::canReconnectResonanceSystems();
	}
	bool doReconnectResonanceSystems(int a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doReconnectResonanceSystems");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return MBReconUserHooks::doReconnectResonanceSystems(a0, a1);
	}
	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "initAfterBeams");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "biasedSelectionWeight");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoProcessLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoProcessLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "numberVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoMPIEmission");
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
	bool canChangeFragPar() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MBReconUserHooks *>(this), "onStat");
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

// Pythia8::TopReconUserHooks file:Pythia8Plugins/ColourReconnectionHooks.h line:619
struct PyCallBack_Pythia8_TopReconUserHooks : public Pythia8::TopReconUserHooks {
	using Pythia8::TopReconUserHooks::TopReconUserHooks;

	bool canReconnectResonanceSystems() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canReconnectResonanceSystems");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return TopReconUserHooks::canReconnectResonanceSystems();
	}
	bool doReconnectResonanceSystems(int a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doReconnectResonanceSystems");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return TopReconUserHooks::doReconnectResonanceSystems(a0, a1);
	}
	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "initAfterBeams");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "biasedSelectionWeight");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoProcessLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoProcessLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "numberVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoMPIEmission");
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
	bool canChangeFragPar() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::TopReconUserHooks *>(this), "onStat");
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

// Pythia8::LHAupAlpgen file:Pythia8Plugins/GeneratorInput.h line:83
struct PyCallBack_Pythia8_LHAupAlpgen : public Pythia8::LHAupAlpgen {
	using Pythia8::LHAupAlpgen::LHAupAlpgen;

	bool fileFound() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupAlpgen *>(this), "fileFound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupAlpgen::fileFound();
	}
	bool setInit() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupAlpgen *>(this), "setInit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupAlpgen::setInit();
	}
	bool setEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupAlpgen *>(this), "setEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupAlpgen::setEvent(a0);
	}
	void newEventFile(const char * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupAlpgen *>(this), "newEventFile");
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
	bool useExternal() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupAlpgen *>(this), "useExternal");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupAlpgen *>(this), "skipEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupAlpgen *>(this), "openLHEF");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupAlpgen *>(this), "closeLHEF");
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

// Pythia8::AlpgenHooks file:Pythia8Plugins/GeneratorInput.h line:133
struct PyCallBack_Pythia8_AlpgenHooks : public Pythia8::AlpgenHooks {
	using Pythia8::AlpgenHooks::AlpgenHooks;

	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return AlpgenHooks::initAfterBeams();
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "biasedSelectionWeight");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoProcessLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoProcessLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "numberVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::AlpgenHooks *>(this), "onStat");
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

// Pythia8::HJSlowJet file:Pythia8Plugins/JetMatching.h line:29
struct PyCallBack_Pythia8_HJSlowJet : public Pythia8::HJSlowJet {
	using Pythia8::HJSlowJet::HJSlowJet;

	void findNext() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HJSlowJet *>(this), "findNext");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HJSlowJet::findNext();
	}
	bool doStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HJSlowJet *>(this), "doStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SlowJet::doStep();
	}
};

// Pythia8::JetMatching file:Pythia8Plugins/JetMatching.h line:85
struct PyCallBack_Pythia8_JetMatching : public Pythia8::JetMatching {
	using Pythia8::JetMatching::JetMatching;

	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::initAfterBeams\"");
	}
	bool canVetoProcessLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::canVetoProcessLevel();
	}
	bool doVetoProcessLevel(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::doVetoProcessLevel(a0);
	}
	bool canVetoPartonLevelEarly() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::canVetoPartonLevelEarly();
	}
	bool doVetoPartonLevelEarly(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::doVetoPartonLevelEarly(a0);
	}
	int numberVetoStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "numberVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return JetMatching::numberVetoStep();
	}
	bool canVetoStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::canVetoStep();
	}
	bool doVetoStep(int a0, int a1, int a2, const class Pythia8::Event & a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::doVetoStep(a0, a1, a2, a3);
	}
	void sortIncomingProcess(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "sortIncomingProcess");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::sortIncomingProcess\"");
	}
	void jetAlgorithmInput(const class Pythia8::Event & a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "jetAlgorithmInput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::jetAlgorithmInput\"");
	}
	void runJetAlgorithm() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "runJetAlgorithm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::runJetAlgorithm\"");
	}
	bool matchPartonsToJets(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "matchPartonsToJets");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJets\"");
	}
	int matchPartonsToJetsLight() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "matchPartonsToJetsLight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJetsLight\"");
	}
	int matchPartonsToJetsHeavy() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "matchPartonsToJetsHeavy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJetsHeavy\"");
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "biasedSelectionWeight");
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
	bool canSetLowEnergySigma(int a0, int a1) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoPT");
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
	bool canVetoMPIStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoMPIStep");
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
	bool retryPartonLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatching *>(this), "onStat");
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

void bind_Pythia8Plugins_aMCatNLOHooks(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::amcnlo_unitarised_interface file:Pythia8Plugins/aMCatNLOHooks.h line:28
		pybind11::class_<Pythia8::amcnlo_unitarised_interface, std::shared_ptr<Pythia8::amcnlo_unitarised_interface>, PyCallBack_Pythia8_amcnlo_unitarised_interface, Pythia8::UserHooks> cl(M("Pythia8"), "amcnlo_unitarised_interface", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::amcnlo_unitarised_interface(); }, [](){ return new PyCallBack_Pythia8_amcnlo_unitarised_interface(); } ) );
		cl.def( pybind11::init<int>(), pybind11::arg("mergingSchemeIn") );

		cl.def("getNormFactor", (double (Pythia8::amcnlo_unitarised_interface::*)()) &Pythia8::amcnlo_unitarised_interface::getNormFactor, "C++: Pythia8::amcnlo_unitarised_interface::getNormFactor() --> double");
		cl.def("canVetoProcessLevel", (bool (Pythia8::amcnlo_unitarised_interface::*)()) &Pythia8::amcnlo_unitarised_interface::canVetoProcessLevel, "C++: Pythia8::amcnlo_unitarised_interface::canVetoProcessLevel() --> bool");
		cl.def("doVetoProcessLevel", (bool (Pythia8::amcnlo_unitarised_interface::*)(class Pythia8::Event &)) &Pythia8::amcnlo_unitarised_interface::doVetoProcessLevel, "C++: Pythia8::amcnlo_unitarised_interface::doVetoProcessLevel(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("assign", (class Pythia8::amcnlo_unitarised_interface & (Pythia8::amcnlo_unitarised_interface::*)(const class Pythia8::amcnlo_unitarised_interface &)) &Pythia8::amcnlo_unitarised_interface::operator=, "C++: Pythia8::amcnlo_unitarised_interface::operator=(const class Pythia8::amcnlo_unitarised_interface &) --> class Pythia8::amcnlo_unitarised_interface &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::MBReconUserHooks file:Pythia8Plugins/ColourReconnectionHooks.h line:33
		pybind11::class_<Pythia8::MBReconUserHooks, std::shared_ptr<Pythia8::MBReconUserHooks>, PyCallBack_Pythia8_MBReconUserHooks, Pythia8::UserHooks> cl(M("Pythia8"), "MBReconUserHooks", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::MBReconUserHooks(); }, [](){ return new PyCallBack_Pythia8_MBReconUserHooks(); } ), "doc");
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::MBReconUserHooks(a0); }, [](int const & a0){ return new PyCallBack_Pythia8_MBReconUserHooks(a0); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Pythia8::MBReconUserHooks(a0, a1); }, [](int const & a0, int const & a1){ return new PyCallBack_Pythia8_MBReconUserHooks(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2){ return new Pythia8::MBReconUserHooks(a0, a1, a2); }, [](int const & a0, int const & a1, double const & a2){ return new PyCallBack_Pythia8_MBReconUserHooks(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<int, int, double, double>(), pybind11::arg("modeIn"), pybind11::arg("flipIn"), pybind11::arg("dLamCutIn"), pybind11::arg("fracGluonIn") );

		cl.def("canReconnectResonanceSystems", (bool (Pythia8::MBReconUserHooks::*)()) &Pythia8::MBReconUserHooks::canReconnectResonanceSystems, "C++: Pythia8::MBReconUserHooks::canReconnectResonanceSystems() --> bool");
		cl.def("doReconnectResonanceSystems", (bool (Pythia8::MBReconUserHooks::*)(int, class Pythia8::Event &)) &Pythia8::MBReconUserHooks::doReconnectResonanceSystems, "C++: Pythia8::MBReconUserHooks::doReconnectResonanceSystems(int, class Pythia8::Event &) --> bool", pybind11::arg(""), pybind11::arg("event"));
		cl.def("assign", (class Pythia8::MBReconUserHooks & (Pythia8::MBReconUserHooks::*)(const class Pythia8::MBReconUserHooks &)) &Pythia8::MBReconUserHooks::operator=, "C++: Pythia8::MBReconUserHooks::operator=(const class Pythia8::MBReconUserHooks &) --> class Pythia8::MBReconUserHooks &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::TopReconUserHooks file:Pythia8Plugins/ColourReconnectionHooks.h line:619
		pybind11::class_<Pythia8::TopReconUserHooks, std::shared_ptr<Pythia8::TopReconUserHooks>, PyCallBack_Pythia8_TopReconUserHooks, Pythia8::UserHooks> cl(M("Pythia8"), "TopReconUserHooks", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::TopReconUserHooks(); }, [](){ return new PyCallBack_Pythia8_TopReconUserHooks(); } ), "doc");
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::TopReconUserHooks(a0); }, [](int const & a0){ return new PyCallBack_Pythia8_TopReconUserHooks(a0); } ), "doc");
		cl.def( pybind11::init<int, double>(), pybind11::arg("modeIn"), pybind11::arg("strengthIn") );

		cl.def("canReconnectResonanceSystems", (bool (Pythia8::TopReconUserHooks::*)()) &Pythia8::TopReconUserHooks::canReconnectResonanceSystems, "C++: Pythia8::TopReconUserHooks::canReconnectResonanceSystems() --> bool");
		cl.def("doReconnectResonanceSystems", (bool (Pythia8::TopReconUserHooks::*)(int, class Pythia8::Event &)) &Pythia8::TopReconUserHooks::doReconnectResonanceSystems, "C++: Pythia8::TopReconUserHooks::doReconnectResonanceSystems(int, class Pythia8::Event &) --> bool", pybind11::arg(""), pybind11::arg("event"));
		cl.def("assign", (class Pythia8::TopReconUserHooks & (Pythia8::TopReconUserHooks::*)(const class Pythia8::TopReconUserHooks &)) &Pythia8::TopReconUserHooks::operator=, "C++: Pythia8::TopReconUserHooks::operator=(const class Pythia8::TopReconUserHooks &) --> class Pythia8::TopReconUserHooks &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::AlpgenPar file:Pythia8Plugins/GeneratorInput.h line:35
		pybind11::class_<Pythia8::AlpgenPar, std::shared_ptr<Pythia8::AlpgenPar>> cl(M("Pythia8"), "AlpgenPar", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::AlpgenPar(); } ) );
		cl.def( pybind11::init( [](Pythia8::AlpgenPar const &o){ return new Pythia8::AlpgenPar(o); } ) );
		cl.def("parse", (bool (Pythia8::AlpgenPar::*)(const std::string)) &Pythia8::AlpgenPar::parse, "C++: Pythia8::AlpgenPar::parse(const std::string) --> bool", pybind11::arg("paramStr"));
		cl.def("extractRunParam", (void (Pythia8::AlpgenPar::*)(std::string)) &Pythia8::AlpgenPar::extractRunParam, "C++: Pythia8::AlpgenPar::extractRunParam(std::string) --> void", pybind11::arg("line"));
		cl.def("haveParam", (bool (Pythia8::AlpgenPar::*)(const std::string &)) &Pythia8::AlpgenPar::haveParam, "C++: Pythia8::AlpgenPar::haveParam(const std::string &) --> bool", pybind11::arg("paramIn"));
		cl.def("getParam", (double (Pythia8::AlpgenPar::*)(const std::string &)) &Pythia8::AlpgenPar::getParam, "C++: Pythia8::AlpgenPar::getParam(const std::string &) --> double", pybind11::arg("paramIn"));
		cl.def("getParamAsInt", (int (Pythia8::AlpgenPar::*)(const std::string &)) &Pythia8::AlpgenPar::getParamAsInt, "C++: Pythia8::AlpgenPar::getParamAsInt(const std::string &) --> int", pybind11::arg("paramIn"));
		cl.def("printParams", (void (Pythia8::AlpgenPar::*)()) &Pythia8::AlpgenPar::printParams, "C++: Pythia8::AlpgenPar::printParams() --> void");
	}
	{ // Pythia8::LHAupAlpgen file:Pythia8Plugins/GeneratorInput.h line:83
		pybind11::class_<Pythia8::LHAupAlpgen, std::shared_ptr<Pythia8::LHAupAlpgen>, PyCallBack_Pythia8_LHAupAlpgen, Pythia8::LHAup> cl(M("Pythia8"), "LHAupAlpgen", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<const char *>(), pybind11::arg("baseFNin") );

		cl.def("fileFound", (bool (Pythia8::LHAupAlpgen::*)()) &Pythia8::LHAupAlpgen::fileFound, "C++: Pythia8::LHAupAlpgen::fileFound() --> bool");
		cl.def("setInit", (bool (Pythia8::LHAupAlpgen::*)()) &Pythia8::LHAupAlpgen::setInit, "C++: Pythia8::LHAupAlpgen::setInit() --> bool");
		cl.def("setEvent", (bool (Pythia8::LHAupAlpgen::*)(int)) &Pythia8::LHAupAlpgen::setEvent, "C++: Pythia8::LHAupAlpgen::setEvent(int) --> bool", pybind11::arg(""));
		cl.def("printParticles", (void (Pythia8::LHAupAlpgen::*)()) &Pythia8::LHAupAlpgen::printParticles, "C++: Pythia8::LHAupAlpgen::printParticles() --> void");
	}
	{ // Pythia8::AlpgenHooks file:Pythia8Plugins/GeneratorInput.h line:133
		pybind11::class_<Pythia8::AlpgenHooks, std::shared_ptr<Pythia8::AlpgenHooks>, PyCallBack_Pythia8_AlpgenHooks, Pythia8::UserHooks> cl(M("Pythia8"), "AlpgenHooks", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Pythia &>(), pybind11::arg("pythia") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_AlpgenHooks const &o){ return new PyCallBack_Pythia8_AlpgenHooks(o); } ) );
		cl.def( pybind11::init( [](Pythia8::AlpgenHooks const &o){ return new Pythia8::AlpgenHooks(o); } ) );
		cl.def("initAfterBeams", (bool (Pythia8::AlpgenHooks::*)()) &Pythia8::AlpgenHooks::initAfterBeams, "C++: Pythia8::AlpgenHooks::initAfterBeams() --> bool");
		cl.def("assign", (class Pythia8::AlpgenHooks & (Pythia8::AlpgenHooks::*)(const class Pythia8::AlpgenHooks &)) &Pythia8::AlpgenHooks::operator=, "C++: Pythia8::AlpgenHooks::operator=(const class Pythia8::AlpgenHooks &) --> class Pythia8::AlpgenHooks &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::MadgraphPar file:Pythia8Plugins/GeneratorInput.h line:156
		pybind11::class_<Pythia8::MadgraphPar, std::shared_ptr<Pythia8::MadgraphPar>> cl(M("Pythia8"), "MadgraphPar", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::MadgraphPar(); } ) );
		cl.def( pybind11::init( [](Pythia8::MadgraphPar const &o){ return new Pythia8::MadgraphPar(o); } ) );
		cl.def("parse", (bool (Pythia8::MadgraphPar::*)(const std::string)) &Pythia8::MadgraphPar::parse, "C++: Pythia8::MadgraphPar::parse(const std::string) --> bool", pybind11::arg("paramStr"));
		cl.def("extractRunParam", (void (Pythia8::MadgraphPar::*)(std::string)) &Pythia8::MadgraphPar::extractRunParam, "C++: Pythia8::MadgraphPar::extractRunParam(std::string) --> void", pybind11::arg("line"));
		cl.def("haveParam", (bool (Pythia8::MadgraphPar::*)(const std::string &)) &Pythia8::MadgraphPar::haveParam, "C++: Pythia8::MadgraphPar::haveParam(const std::string &) --> bool", pybind11::arg("paramIn"));
		cl.def("getParam", (double (Pythia8::MadgraphPar::*)(const std::string &)) &Pythia8::MadgraphPar::getParam, "C++: Pythia8::MadgraphPar::getParam(const std::string &) --> double", pybind11::arg("paramIn"));
		cl.def("getParamAsInt", (int (Pythia8::MadgraphPar::*)(const std::string &)) &Pythia8::MadgraphPar::getParamAsInt, "C++: Pythia8::MadgraphPar::getParamAsInt(const std::string &) --> int", pybind11::arg("paramIn"));
		cl.def("printParams", (void (Pythia8::MadgraphPar::*)()) &Pythia8::MadgraphPar::printParams, "C++: Pythia8::MadgraphPar::printParams() --> void");
	}
	{ // Pythia8::HJSlowJet file:Pythia8Plugins/JetMatching.h line:29
		pybind11::class_<Pythia8::HJSlowJet, std::shared_ptr<Pythia8::HJSlowJet>, PyCallBack_Pythia8_HJSlowJet, Pythia8::SlowJet> cl(M("Pythia8"), "HJSlowJet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](int const & a0, double const & a1){ return new Pythia8::HJSlowJet(a0, a1); }, [](int const & a0, double const & a1){ return new PyCallBack_Pythia8_HJSlowJet(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2){ return new Pythia8::HJSlowJet(a0, a1, a2); }, [](int const & a0, double const & a1, double const & a2){ return new PyCallBack_Pythia8_HJSlowJet(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3){ return new Pythia8::HJSlowJet(a0, a1, a2, a3); }, [](int const & a0, double const & a1, double const & a2, double const & a3){ return new PyCallBack_Pythia8_HJSlowJet(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4){ return new Pythia8::HJSlowJet(a0, a1, a2, a3, a4); }, [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4){ return new PyCallBack_Pythia8_HJSlowJet(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5){ return new Pythia8::HJSlowJet(a0, a1, a2, a3, a4, a5); }, [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5){ return new PyCallBack_Pythia8_HJSlowJet(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5, class Pythia8::SlowJetHook * a6){ return new Pythia8::HJSlowJet(a0, a1, a2, a3, a4, a5, a6); }, [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5, class Pythia8::SlowJetHook * a6){ return new PyCallBack_Pythia8_HJSlowJet(a0, a1, a2, a3, a4, a5, a6); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5, class Pythia8::SlowJetHook * a6, bool const & a7){ return new Pythia8::HJSlowJet(a0, a1, a2, a3, a4, a5, a6, a7); }, [](int const & a0, double const & a1, double const & a2, double const & a3, int const & a4, int const & a5, class Pythia8::SlowJetHook * a6, bool const & a7){ return new PyCallBack_Pythia8_HJSlowJet(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc");
		cl.def( pybind11::init<int, double, double, double, int, int, class Pythia8::SlowJetHook *, bool, bool>(), pybind11::arg("powerIn"), pybind11::arg("Rin"), pybind11::arg("pTjetMinIn"), pybind11::arg("etaMaxIn"), pybind11::arg("selectIn"), pybind11::arg("massSetIn"), pybind11::arg("sjHookPtrIn"), pybind11::arg("useFJcoreIn"), pybind11::arg("useStandardRin") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_HJSlowJet const &o){ return new PyCallBack_Pythia8_HJSlowJet(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HJSlowJet const &o){ return new Pythia8::HJSlowJet(o); } ) );
		cl.def("findNext", (void (Pythia8::HJSlowJet::*)()) &Pythia8::HJSlowJet::findNext, "C++: Pythia8::HJSlowJet::findNext() --> void");
		cl.def("assign", (class Pythia8::HJSlowJet & (Pythia8::HJSlowJet::*)(const class Pythia8::HJSlowJet &)) &Pythia8::HJSlowJet::operator=, "C++: Pythia8::HJSlowJet::operator=(const class Pythia8::HJSlowJet &) --> class Pythia8::HJSlowJet &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::JetMatching file:Pythia8Plugins/JetMatching.h line:85
		pybind11::class_<Pythia8::JetMatching, std::shared_ptr<Pythia8::JetMatching>, PyCallBack_Pythia8_JetMatching, Pythia8::UserHooks> cl(M("Pythia8"), "JetMatching", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_JetMatching(); } ) );
		cl.def(pybind11::init<PyCallBack_Pythia8_JetMatching const &>());

		pybind11::enum_<Pythia8::JetMatching::vetoStatus>(cl, "vetoStatus", pybind11::arithmetic(), "")
			.value("NONE", Pythia8::JetMatching::vetoStatus::NONE)
			.value("LESS_JETS", Pythia8::JetMatching::vetoStatus::LESS_JETS)
			.value("MORE_JETS", Pythia8::JetMatching::vetoStatus::MORE_JETS)
			.value("HARD_JET", Pythia8::JetMatching::vetoStatus::HARD_JET)
			.value("UNMATCHED_PARTON", Pythia8::JetMatching::vetoStatus::UNMATCHED_PARTON)
			.value("INCLUSIVE_VETO", Pythia8::JetMatching::vetoStatus::INCLUSIVE_VETO)
			.export_values();


		pybind11::enum_<Pythia8::JetMatching::partonTypes>(cl, "partonTypes", pybind11::arithmetic(), "")
			.value("ID_CHARM", Pythia8::JetMatching::partonTypes::ID_CHARM)
			.value("ID_BOT", Pythia8::JetMatching::partonTypes::ID_BOT)
			.value("ID_TOP", Pythia8::JetMatching::partonTypes::ID_TOP)
			.value("ID_LEPMIN", Pythia8::JetMatching::partonTypes::ID_LEPMIN)
			.value("ID_LEPMAX", Pythia8::JetMatching::partonTypes::ID_LEPMAX)
			.value("ID_GLUON", Pythia8::JetMatching::partonTypes::ID_GLUON)
			.value("ID_PHOTON", Pythia8::JetMatching::partonTypes::ID_PHOTON)
			.export_values();

		cl.def_readwrite("doMerge", &Pythia8::JetMatching::doMerge);
		cl.def_readwrite("doShowerKt", &Pythia8::JetMatching::doShowerKt);
		cl.def_readwrite("nJetMax", &Pythia8::JetMatching::nJetMax);
		cl.def_readwrite("nJet", &Pythia8::JetMatching::nJet);
		cl.def_readwrite("jetAlgorithm", &Pythia8::JetMatching::jetAlgorithm);
		cl.def_readwrite("eTjetMin", &Pythia8::JetMatching::eTjetMin);
		cl.def_readwrite("coneRadius", &Pythia8::JetMatching::coneRadius);
		cl.def_readwrite("etaJetMax", &Pythia8::JetMatching::etaJetMax);
		cl.def_readwrite("etaJetMaxAlgo", &Pythia8::JetMatching::etaJetMaxAlgo);
		cl.def_readwrite("slowJetPower", &Pythia8::JetMatching::slowJetPower);
		cl.def_readwrite("eventProcessOrig", &Pythia8::JetMatching::eventProcessOrig);
		cl.def_readwrite("eventProcess", &Pythia8::JetMatching::eventProcess);
		cl.def_readwrite("workEventJet", &Pythia8::JetMatching::workEventJet);
		cl.def_readwrite("jetMomenta", &Pythia8::JetMatching::jetMomenta);
		cl.def_readwrite("nEta", &Pythia8::JetMatching::nEta);
		cl.def_readwrite("nPhi", &Pythia8::JetMatching::nPhi);
		cl.def_readwrite("eTseed", &Pythia8::JetMatching::eTseed);
		cl.def_readwrite("eTthreshold", &Pythia8::JetMatching::eTthreshold);
		cl.def_readwrite("jetAllow", &Pythia8::JetMatching::jetAllow);
		cl.def_readwrite("jetMatch", &Pythia8::JetMatching::jetMatch);
		cl.def_readwrite("exclusiveMode", &Pythia8::JetMatching::exclusiveMode);
		cl.def_readwrite("coneMatchLight", &Pythia8::JetMatching::coneMatchLight);
		cl.def_readwrite("coneRadiusHeavy", &Pythia8::JetMatching::coneRadiusHeavy);
		cl.def_readwrite("coneMatchHeavy", &Pythia8::JetMatching::coneMatchHeavy);
		cl.def_readwrite("exclusive", &Pythia8::JetMatching::exclusive);
		cl.def_readwrite("eTpTlightMin", &Pythia8::JetMatching::eTpTlightMin);
		cl.def_readwrite("messages", &Pythia8::JetMatching::messages);
		cl.def("initAfterBeams", (bool (Pythia8::JetMatching::*)()) &Pythia8::JetMatching::initAfterBeams, "C++: Pythia8::JetMatching::initAfterBeams() --> bool");
		cl.def("canVetoProcessLevel", (bool (Pythia8::JetMatching::*)()) &Pythia8::JetMatching::canVetoProcessLevel, "C++: Pythia8::JetMatching::canVetoProcessLevel() --> bool");
		cl.def("doVetoProcessLevel", (bool (Pythia8::JetMatching::*)(class Pythia8::Event &)) &Pythia8::JetMatching::doVetoProcessLevel, "C++: Pythia8::JetMatching::doVetoProcessLevel(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("canVetoPartonLevelEarly", (bool (Pythia8::JetMatching::*)()) &Pythia8::JetMatching::canVetoPartonLevelEarly, "C++: Pythia8::JetMatching::canVetoPartonLevelEarly() --> bool");
		cl.def("doVetoPartonLevelEarly", (bool (Pythia8::JetMatching::*)(const class Pythia8::Event &)) &Pythia8::JetMatching::doVetoPartonLevelEarly, "C++: Pythia8::JetMatching::doVetoPartonLevelEarly(const class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("numberVetoStep", (int (Pythia8::JetMatching::*)()) &Pythia8::JetMatching::numberVetoStep, "C++: Pythia8::JetMatching::numberVetoStep() --> int");
		cl.def("canVetoStep", (bool (Pythia8::JetMatching::*)()) &Pythia8::JetMatching::canVetoStep, "C++: Pythia8::JetMatching::canVetoStep() --> bool");
		cl.def("doVetoStep", (bool (Pythia8::JetMatching::*)(int, int, int, const class Pythia8::Event &)) &Pythia8::JetMatching::doVetoStep, "C++: Pythia8::JetMatching::doVetoStep(int, int, int, const class Pythia8::Event &) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("sortIncomingProcess", (void (Pythia8::JetMatching::*)(const class Pythia8::Event &)) &Pythia8::JetMatching::sortIncomingProcess, "C++: Pythia8::JetMatching::sortIncomingProcess(const class Pythia8::Event &) --> void", pybind11::arg(""));
		cl.def("jetAlgorithmInput", (void (Pythia8::JetMatching::*)(const class Pythia8::Event &, int)) &Pythia8::JetMatching::jetAlgorithmInput, "C++: Pythia8::JetMatching::jetAlgorithmInput(const class Pythia8::Event &, int) --> void", pybind11::arg(""), pybind11::arg(""));
		cl.def("runJetAlgorithm", (void (Pythia8::JetMatching::*)()) &Pythia8::JetMatching::runJetAlgorithm, "C++: Pythia8::JetMatching::runJetAlgorithm() --> void");
		cl.def("matchPartonsToJets", (bool (Pythia8::JetMatching::*)(int)) &Pythia8::JetMatching::matchPartonsToJets, "C++: Pythia8::JetMatching::matchPartonsToJets(int) --> bool", pybind11::arg(""));
		cl.def("matchPartonsToJetsLight", (int (Pythia8::JetMatching::*)()) &Pythia8::JetMatching::matchPartonsToJetsLight, "C++: Pythia8::JetMatching::matchPartonsToJetsLight() --> int");
		cl.def("matchPartonsToJetsHeavy", (int (Pythia8::JetMatching::*)()) &Pythia8::JetMatching::matchPartonsToJetsHeavy, "C++: Pythia8::JetMatching::matchPartonsToJetsHeavy() --> int");
		cl.def("errorMsg", (void (Pythia8::JetMatching::*)(std::string)) &Pythia8::JetMatching::errorMsg, "C++: Pythia8::JetMatching::errorMsg(std::string) --> void", pybind11::arg("messageIn"));
		cl.def("assign", (class Pythia8::JetMatching & (Pythia8::JetMatching::*)(const class Pythia8::JetMatching &)) &Pythia8::JetMatching::operator=, "C++: Pythia8::JetMatching::operator=(const class Pythia8::JetMatching &) --> class Pythia8::JetMatching &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/GammaKinematics.h>
#include <Pythia8/HIInfo.h>
#include <Pythia8/HadronLevel.h>
#include <Pythia8/HeavyIons.h>
#include <Pythia8/Info.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/Merging.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/NucleonExcitations.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ParticleDecays.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/RHadrons.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/SLHAinterface.h>
#include <Pythia8/Settings.h>
#include <Pythia8/ShowerModel.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaProcess.h>
#include <Pythia8/StringFragmentation.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/TimeShower.h>
#include <Pythia8/UserHooks.h>
#include <Pythia8Plugins/CombineMatchingInput.h>
#include <Pythia8Plugins/FlavorVariations.h>
#include <Pythia8Plugins/InputParser.h>
#include <Pythia8Plugins/JetMatching.h>
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

// Pythia8::JetMatchingAlpgen file:Pythia8Plugins/JetMatching.h line:241
struct PyCallBack_Pythia8_JetMatchingAlpgen : public Pythia8::JetMatchingAlpgen {
	using Pythia8::JetMatchingAlpgen::JetMatchingAlpgen;

	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingAlpgen::initAfterBeams();
	}
	bool canVetoProcessLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::canVetoProcessLevel();
	}
	bool doVetoProcessLevel(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::doVetoProcessLevel(a0);
	}
	bool canVetoPartonLevelEarly() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::canVetoPartonLevelEarly();
	}
	bool doVetoPartonLevelEarly(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::doVetoPartonLevelEarly(a0);
	}
	int numberVetoStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "numberVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return JetMatching::numberVetoStep();
	}
	bool canVetoStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::canVetoStep();
	}
	bool doVetoStep(int a0, int a1, int a2, const class Pythia8::Event & a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::doVetoStep(a0, a1, a2, a3);
	}
	void sortIncomingProcess(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "sortIncomingProcess");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::sortIncomingProcess\"");
	}
	void jetAlgorithmInput(const class Pythia8::Event & a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "jetAlgorithmInput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::jetAlgorithmInput\"");
	}
	void runJetAlgorithm() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "runJetAlgorithm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::runJetAlgorithm\"");
	}
	bool matchPartonsToJets(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "matchPartonsToJets");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJets\"");
	}
	int matchPartonsToJetsLight() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "matchPartonsToJetsLight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJetsLight\"");
	}
	int matchPartonsToJetsHeavy() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "matchPartonsToJetsHeavy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJetsHeavy\"");
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "biasedSelectionWeight");
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
	bool canSetLowEnergySigma(int a0, int a1) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoPT");
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
	bool canVetoMPIStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoMPIStep");
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
	bool retryPartonLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgen *>(this), "onStat");
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

// Pythia8::JetMatchingMadgraph file:Pythia8Plugins/JetMatching.h line:274
struct PyCallBack_Pythia8_JetMatchingMadgraph : public Pythia8::JetMatchingMadgraph {
	using Pythia8::JetMatchingMadgraph::JetMatchingMadgraph;

	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraph::initAfterBeams();
	}
	bool canVetoProcessLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraph::canVetoProcessLevel();
	}
	bool doVetoProcessLevel(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraph::doVetoProcessLevel(a0);
	}
	int numberVetoStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "numberVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return JetMatchingMadgraph::numberVetoStep();
	}
	bool canVetoStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraph::canVetoStep();
	}
	bool doVetoStep(int a0, int a1, int a2, const class Pythia8::Event & a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoStep");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraph::doVetoStep(a0, a1, a2, a3);
	}
	void sortIncomingProcess(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "sortIncomingProcess");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return JetMatchingMadgraph::sortIncomingProcess(a0);
	}
	void jetAlgorithmInput(const class Pythia8::Event & a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "jetAlgorithmInput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return JetMatchingMadgraph::jetAlgorithmInput(a0, a1);
	}
	void runJetAlgorithm() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "runJetAlgorithm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return JetMatchingMadgraph::runJetAlgorithm();
	}
	bool matchPartonsToJets(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "matchPartonsToJets");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraph::matchPartonsToJets(a0);
	}
	int matchPartonsToJetsLight() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "matchPartonsToJetsLight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return JetMatchingMadgraph::matchPartonsToJetsLight();
	}
	int matchPartonsToJetsHeavy() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "matchPartonsToJetsHeavy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return JetMatchingMadgraph::matchPartonsToJetsHeavy();
	}
	bool canVetoPartonLevelEarly() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::canVetoPartonLevelEarly();
	}
	bool doVetoPartonLevelEarly(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatching::doVetoPartonLevelEarly(a0);
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "biasedSelectionWeight");
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
	bool canSetLowEnergySigma(int a0, int a1) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoPT");
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
	bool canVetoMPIStep() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoMPIStep");
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
	bool retryPartonLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraph *>(this), "onStat");
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

// Pythia8::JetMatchingAlpgenInputAlpgen file:Pythia8Plugins/CombineMatchingInput.h line:31
struct PyCallBack_Pythia8_JetMatchingAlpgenInputAlpgen : public Pythia8::JetMatchingAlpgenInputAlpgen {
	using Pythia8::JetMatchingAlpgenInputAlpgen::JetMatchingAlpgenInputAlpgen;

	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingAlpgenInputAlpgen::initAfterBeams();
	}
	bool canVetoProcessLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingAlpgenInputAlpgen::canVetoProcessLevel();
	}
	bool doVetoProcessLevel(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingAlpgenInputAlpgen::doVetoProcessLevel(a0);
	}
	bool canVetoPartonLevelEarly() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingAlpgenInputAlpgen::canVetoPartonLevelEarly();
	}
	bool doVetoPartonLevelEarly(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingAlpgenInputAlpgen::doVetoPartonLevelEarly(a0);
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "biasedSelectionWeight");
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
	bool canSetLowEnergySigma(int a0, int a1) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "numberVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoMPIStep");
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
	bool retryPartonLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "onStat");
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
	void sortIncomingProcess(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "sortIncomingProcess");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::sortIncomingProcess\"");
	}
	void jetAlgorithmInput(const class Pythia8::Event & a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "jetAlgorithmInput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::jetAlgorithmInput\"");
	}
	void runJetAlgorithm() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "runJetAlgorithm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::runJetAlgorithm\"");
	}
	bool matchPartonsToJets(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "matchPartonsToJets");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJets\"");
	}
	int matchPartonsToJetsLight() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "matchPartonsToJetsLight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJetsLight\"");
	}
	int matchPartonsToJetsHeavy() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingAlpgenInputAlpgen *>(this), "matchPartonsToJetsHeavy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"JetMatching::matchPartonsToJetsHeavy\"");
	}
};

// Pythia8::JetMatchingMadgraphInputAlpgen file:Pythia8Plugins/CombineMatchingInput.h line:73
struct PyCallBack_Pythia8_JetMatchingMadgraphInputAlpgen : public Pythia8::JetMatchingMadgraphInputAlpgen {
	using Pythia8::JetMatchingMadgraphInputAlpgen::JetMatchingMadgraphInputAlpgen;

	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraphInputAlpgen::initAfterBeams();
	}
	bool canVetoProcessLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraphInputAlpgen::canVetoProcessLevel();
	}
	bool doVetoProcessLevel(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraphInputAlpgen::doVetoProcessLevel(a0);
	}
	bool canVetoPartonLevelEarly() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraphInputAlpgen::canVetoPartonLevelEarly();
	}
	bool doVetoPartonLevelEarly(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoPartonLevelEarly");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraphInputAlpgen::doVetoPartonLevelEarly(a0);
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "biasedSelectionWeight");
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
	bool canSetLowEnergySigma(int a0, int a1) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "numberVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoMPIStep");
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
	bool retryPartonLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "onStat");
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
	void sortIncomingProcess(const class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "sortIncomingProcess");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return JetMatchingMadgraph::sortIncomingProcess(a0);
	}
	void jetAlgorithmInput(const class Pythia8::Event & a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "jetAlgorithmInput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return JetMatchingMadgraph::jetAlgorithmInput(a0, a1);
	}
	void runJetAlgorithm() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "runJetAlgorithm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return JetMatchingMadgraph::runJetAlgorithm();
	}
	bool matchPartonsToJets(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "matchPartonsToJets");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return JetMatchingMadgraph::matchPartonsToJets(a0);
	}
	int matchPartonsToJetsLight() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "matchPartonsToJetsLight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return JetMatchingMadgraph::matchPartonsToJetsLight();
	}
	int matchPartonsToJetsHeavy() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::JetMatchingMadgraphInputAlpgen *>(this), "matchPartonsToJetsHeavy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return JetMatchingMadgraph::matchPartonsToJetsHeavy();
	}
};

void bind_Pythia8Plugins_JetMatching(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::JetMatchingAlpgen file:Pythia8Plugins/JetMatching.h line:241
		pybind11::class_<Pythia8::JetMatchingAlpgen, std::shared_ptr<Pythia8::JetMatchingAlpgen>, PyCallBack_Pythia8_JetMatchingAlpgen, Pythia8::JetMatching> cl(M("Pythia8"), "JetMatchingAlpgen", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::JetMatchingAlpgen(); }, [](){ return new PyCallBack_Pythia8_JetMatchingAlpgen(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_JetMatchingAlpgen const &o){ return new PyCallBack_Pythia8_JetMatchingAlpgen(o); } ) );
		cl.def( pybind11::init( [](Pythia8::JetMatchingAlpgen const &o){ return new Pythia8::JetMatchingAlpgen(o); } ) );
		cl.def("initAfterBeams", (bool (Pythia8::JetMatchingAlpgen::*)()) &Pythia8::JetMatchingAlpgen::initAfterBeams, "C++: Pythia8::JetMatchingAlpgen::initAfterBeams() --> bool");
		cl.def("assign", (class Pythia8::JetMatchingAlpgen & (Pythia8::JetMatchingAlpgen::*)(const class Pythia8::JetMatchingAlpgen &)) &Pythia8::JetMatchingAlpgen::operator=, "C++: Pythia8::JetMatchingAlpgen::operator=(const class Pythia8::JetMatchingAlpgen &) --> class Pythia8::JetMatchingAlpgen &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::JetMatchingMadgraph file:Pythia8Plugins/JetMatching.h line:274
		pybind11::class_<Pythia8::JetMatchingMadgraph, std::shared_ptr<Pythia8::JetMatchingMadgraph>, PyCallBack_Pythia8_JetMatchingMadgraph, Pythia8::JetMatching> cl(M("Pythia8"), "JetMatchingMadgraph", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::JetMatchingMadgraph(); }, [](){ return new PyCallBack_Pythia8_JetMatchingMadgraph(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_JetMatchingMadgraph const &o){ return new PyCallBack_Pythia8_JetMatchingMadgraph(o); } ) );
		cl.def( pybind11::init( [](Pythia8::JetMatchingMadgraph const &o){ return new Pythia8::JetMatchingMadgraph(o); } ) );
		cl.def("initAfterBeams", (bool (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::initAfterBeams, "C++: Pythia8::JetMatchingMadgraph::initAfterBeams() --> bool");
		cl.def("canVetoProcessLevel", (bool (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::canVetoProcessLevel, "C++: Pythia8::JetMatchingMadgraph::canVetoProcessLevel() --> bool");
		cl.def("doVetoProcessLevel", (bool (Pythia8::JetMatchingMadgraph::*)(class Pythia8::Event &)) &Pythia8::JetMatchingMadgraph::doVetoProcessLevel, "C++: Pythia8::JetMatchingMadgraph::doVetoProcessLevel(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("numberVetoStep", (int (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::numberVetoStep, "C++: Pythia8::JetMatchingMadgraph::numberVetoStep() --> int");
		cl.def("canVetoStep", (bool (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::canVetoStep, "C++: Pythia8::JetMatchingMadgraph::canVetoStep() --> bool");
		cl.def("doVetoStep", (bool (Pythia8::JetMatchingMadgraph::*)(int, int, int, const class Pythia8::Event &)) &Pythia8::JetMatchingMadgraph::doVetoStep, "C++: Pythia8::JetMatchingMadgraph::doVetoStep(int, int, int, const class Pythia8::Event &) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("getDJR", (class std::vector<double, class std::allocator<double> > (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::getDJR, "C++: Pythia8::JetMatchingMadgraph::getDJR() --> class std::vector<double, class std::allocator<double> >");
		cl.def("nMEpartons", (struct std::pair<int, int> (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::nMEpartons, "C++: Pythia8::JetMatchingMadgraph::nMEpartons() --> struct std::pair<int, int>");
		cl.def("getWorkEventJet", (class Pythia8::Event (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::getWorkEventJet, "C++: Pythia8::JetMatchingMadgraph::getWorkEventJet() --> class Pythia8::Event");
		cl.def("getProcessSubset", (class Pythia8::Event (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::getProcessSubset, "C++: Pythia8::JetMatchingMadgraph::getProcessSubset() --> class Pythia8::Event");
		cl.def("getExclusive", (bool (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::getExclusive, "C++: Pythia8::JetMatchingMadgraph::getExclusive() --> bool");
		cl.def("getPTfirst", (double (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::getPTfirst, "C++: Pythia8::JetMatchingMadgraph::getPTfirst() --> double");
		cl.def("sortIncomingProcess", (void (Pythia8::JetMatchingMadgraph::*)(const class Pythia8::Event &)) &Pythia8::JetMatchingMadgraph::sortIncomingProcess, "C++: Pythia8::JetMatchingMadgraph::sortIncomingProcess(const class Pythia8::Event &) --> void", pybind11::arg(""));
		cl.def("jetAlgorithmInput", (void (Pythia8::JetMatchingMadgraph::*)(const class Pythia8::Event &, int)) &Pythia8::JetMatchingMadgraph::jetAlgorithmInput, "C++: Pythia8::JetMatchingMadgraph::jetAlgorithmInput(const class Pythia8::Event &, int) --> void", pybind11::arg(""), pybind11::arg(""));
		cl.def("runJetAlgorithm", (void (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::runJetAlgorithm, "C++: Pythia8::JetMatchingMadgraph::runJetAlgorithm() --> void");
		cl.def("matchPartonsToJets", (bool (Pythia8::JetMatchingMadgraph::*)(int)) &Pythia8::JetMatchingMadgraph::matchPartonsToJets, "C++: Pythia8::JetMatchingMadgraph::matchPartonsToJets(int) --> bool", pybind11::arg(""));
		cl.def("matchPartonsToJetsLight", (int (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::matchPartonsToJetsLight, "C++: Pythia8::JetMatchingMadgraph::matchPartonsToJetsLight() --> int");
		cl.def("matchPartonsToJetsHeavy", (int (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::matchPartonsToJetsHeavy, "C++: Pythia8::JetMatchingMadgraph::matchPartonsToJetsHeavy() --> int");
		cl.def("matchPartonsToJetsOther", (int (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::matchPartonsToJetsOther, "C++: Pythia8::JetMatchingMadgraph::matchPartonsToJetsOther() --> int");
		cl.def("doShowerKtVeto", (bool (Pythia8::JetMatchingMadgraph::*)(double)) &Pythia8::JetMatchingMadgraph::doShowerKtVeto, "C++: Pythia8::JetMatchingMadgraph::doShowerKtVeto(double) --> bool", pybind11::arg("pTfirst"));
		cl.def("clearDJR", (void (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::clearDJR, "C++: Pythia8::JetMatchingMadgraph::clearDJR() --> void");
		cl.def("setDJR", (void (Pythia8::JetMatchingMadgraph::*)(const class Pythia8::Event &)) &Pythia8::JetMatchingMadgraph::setDJR, "C++: Pythia8::JetMatchingMadgraph::setDJR(const class Pythia8::Event &) --> void", pybind11::arg("event"));
		cl.def("clear_nMEpartons", (void (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::clear_nMEpartons, "C++: Pythia8::JetMatchingMadgraph::clear_nMEpartons() --> void");
		cl.def("set_nMEpartons", (void (Pythia8::JetMatchingMadgraph::*)(const int, const int)) &Pythia8::JetMatchingMadgraph::set_nMEpartons, "C++: Pythia8::JetMatchingMadgraph::set_nMEpartons(const int, const int) --> void", pybind11::arg("nOrig"), pybind11::arg("nMatch"));
		cl.def("npNLO", (int (Pythia8::JetMatchingMadgraph::*)()) &Pythia8::JetMatchingMadgraph::npNLO, "C++: Pythia8::JetMatchingMadgraph::npNLO() --> int");
		cl.def("assign", (class Pythia8::JetMatchingMadgraph & (Pythia8::JetMatchingMadgraph::*)(const class Pythia8::JetMatchingMadgraph &)) &Pythia8::JetMatchingMadgraph::operator=, "C++: Pythia8::JetMatchingMadgraph::operator=(const class Pythia8::JetMatchingMadgraph &) --> class Pythia8::JetMatchingMadgraph &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::JetMatchingAlpgenInputAlpgen file:Pythia8Plugins/CombineMatchingInput.h line:31
		pybind11::class_<Pythia8::JetMatchingAlpgenInputAlpgen, std::shared_ptr<Pythia8::JetMatchingAlpgenInputAlpgen>, PyCallBack_Pythia8_JetMatchingAlpgenInputAlpgen, Pythia8::AlpgenHooks, Pythia8::JetMatchingAlpgen> cl(M("Pythia8"), "JetMatchingAlpgenInputAlpgen", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Pythia &>(), pybind11::arg("pythia") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_JetMatchingAlpgenInputAlpgen const &o){ return new PyCallBack_Pythia8_JetMatchingAlpgenInputAlpgen(o); } ) );
		cl.def( pybind11::init( [](Pythia8::JetMatchingAlpgenInputAlpgen const &o){ return new Pythia8::JetMatchingAlpgenInputAlpgen(o); } ) );
		cl.def("initAfterBeams", (bool (Pythia8::JetMatchingAlpgenInputAlpgen::*)()) &Pythia8::JetMatchingAlpgenInputAlpgen::initAfterBeams, "C++: Pythia8::JetMatchingAlpgenInputAlpgen::initAfterBeams() --> bool");
		cl.def("canVetoProcessLevel", (bool (Pythia8::JetMatchingAlpgenInputAlpgen::*)()) &Pythia8::JetMatchingAlpgenInputAlpgen::canVetoProcessLevel, "C++: Pythia8::JetMatchingAlpgenInputAlpgen::canVetoProcessLevel() --> bool");
		cl.def("doVetoProcessLevel", (bool (Pythia8::JetMatchingAlpgenInputAlpgen::*)(class Pythia8::Event &)) &Pythia8::JetMatchingAlpgenInputAlpgen::doVetoProcessLevel, "C++: Pythia8::JetMatchingAlpgenInputAlpgen::doVetoProcessLevel(class Pythia8::Event &) --> bool", pybind11::arg("proc"));
		cl.def("canVetoPartonLevelEarly", (bool (Pythia8::JetMatchingAlpgenInputAlpgen::*)()) &Pythia8::JetMatchingAlpgenInputAlpgen::canVetoPartonLevelEarly, "C++: Pythia8::JetMatchingAlpgenInputAlpgen::canVetoPartonLevelEarly() --> bool");
		cl.def("doVetoPartonLevelEarly", (bool (Pythia8::JetMatchingAlpgenInputAlpgen::*)(const class Pythia8::Event &)) &Pythia8::JetMatchingAlpgenInputAlpgen::doVetoPartonLevelEarly, "C++: Pythia8::JetMatchingAlpgenInputAlpgen::doVetoPartonLevelEarly(const class Pythia8::Event &) --> bool", pybind11::arg("proc"));
		cl.def("assign", (class Pythia8::JetMatchingAlpgenInputAlpgen & (Pythia8::JetMatchingAlpgenInputAlpgen::*)(const class Pythia8::JetMatchingAlpgenInputAlpgen &)) &Pythia8::JetMatchingAlpgenInputAlpgen::operator=, "C++: Pythia8::JetMatchingAlpgenInputAlpgen::operator=(const class Pythia8::JetMatchingAlpgenInputAlpgen &) --> class Pythia8::JetMatchingAlpgenInputAlpgen &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::JetMatchingMadgraphInputAlpgen file:Pythia8Plugins/CombineMatchingInput.h line:73
		pybind11::class_<Pythia8::JetMatchingMadgraphInputAlpgen, std::shared_ptr<Pythia8::JetMatchingMadgraphInputAlpgen>, PyCallBack_Pythia8_JetMatchingMadgraphInputAlpgen, Pythia8::AlpgenHooks, Pythia8::JetMatchingMadgraph> cl(M("Pythia8"), "JetMatchingMadgraphInputAlpgen", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Pythia &>(), pybind11::arg("pythia") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_JetMatchingMadgraphInputAlpgen const &o){ return new PyCallBack_Pythia8_JetMatchingMadgraphInputAlpgen(o); } ) );
		cl.def( pybind11::init( [](Pythia8::JetMatchingMadgraphInputAlpgen const &o){ return new Pythia8::JetMatchingMadgraphInputAlpgen(o); } ) );
		cl.def("initAfterBeams", (bool (Pythia8::JetMatchingMadgraphInputAlpgen::*)()) &Pythia8::JetMatchingMadgraphInputAlpgen::initAfterBeams, "C++: Pythia8::JetMatchingMadgraphInputAlpgen::initAfterBeams() --> bool");
		cl.def("canVetoProcessLevel", (bool (Pythia8::JetMatchingMadgraphInputAlpgen::*)()) &Pythia8::JetMatchingMadgraphInputAlpgen::canVetoProcessLevel, "C++: Pythia8::JetMatchingMadgraphInputAlpgen::canVetoProcessLevel() --> bool");
		cl.def("doVetoProcessLevel", (bool (Pythia8::JetMatchingMadgraphInputAlpgen::*)(class Pythia8::Event &)) &Pythia8::JetMatchingMadgraphInputAlpgen::doVetoProcessLevel, "C++: Pythia8::JetMatchingMadgraphInputAlpgen::doVetoProcessLevel(class Pythia8::Event &) --> bool", pybind11::arg("proc"));
		cl.def("canVetoPartonLevelEarly", (bool (Pythia8::JetMatchingMadgraphInputAlpgen::*)()) &Pythia8::JetMatchingMadgraphInputAlpgen::canVetoPartonLevelEarly, "C++: Pythia8::JetMatchingMadgraphInputAlpgen::canVetoPartonLevelEarly() --> bool");
		cl.def("doVetoPartonLevelEarly", (bool (Pythia8::JetMatchingMadgraphInputAlpgen::*)(const class Pythia8::Event &)) &Pythia8::JetMatchingMadgraphInputAlpgen::doVetoPartonLevelEarly, "C++: Pythia8::JetMatchingMadgraphInputAlpgen::doVetoPartonLevelEarly(const class Pythia8::Event &) --> bool", pybind11::arg("proc"));
		cl.def("assign", (class Pythia8::JetMatchingMadgraphInputAlpgen & (Pythia8::JetMatchingMadgraphInputAlpgen::*)(const class Pythia8::JetMatchingMadgraphInputAlpgen &)) &Pythia8::JetMatchingMadgraphInputAlpgen::operator=, "C++: Pythia8::JetMatchingMadgraphInputAlpgen::operator=(const class Pythia8::JetMatchingMadgraphInputAlpgen &) --> class Pythia8::JetMatchingMadgraphInputAlpgen &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::CombineMatchingInput file:Pythia8Plugins/CombineMatchingInput.h line:112
		pybind11::class_<Pythia8::CombineMatchingInput, std::shared_ptr<Pythia8::CombineMatchingInput>> cl(M("Pythia8"), "CombineMatchingInput", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::CombineMatchingInput(); } ) );
		cl.def_readwrite("hook", &Pythia8::CombineMatchingInput::hook);
		cl.def("setHook", (void (Pythia8::CombineMatchingInput::*)(class Pythia8::Pythia &)) &Pythia8::CombineMatchingInput::setHook, "C++: Pythia8::CombineMatchingInput::setHook(class Pythia8::Pythia &) --> void", pybind11::arg("pythia"));
	}
	{ // Pythia8::FlavorVariations file:Pythia8Plugins/FlavorVariations.h line:20
		pybind11::class_<Pythia8::FlavorVariations, std::shared_ptr<Pythia8::FlavorVariations>> cl(M("Pythia8"), "FlavorVariations", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Settings &>(), pybind11::arg("settings") );

		cl.def( pybind11::init<double, double, double, double>(), pybind11::arg("xi"), pybind11::arg("rho"), pybind11::arg("x"), pybind11::arg("y") );

		cl.def("read", (class std::vector<int, class std::allocator<int> > (Pythia8::FlavorVariations::*)(std::string)) &Pythia8::FlavorVariations::read, "C++: Pythia8::FlavorVariations::read(std::string) --> class std::vector<int, class std::allocator<int> >", pybind11::arg("breaks"));
		cl.def("write", (std::string (Pythia8::FlavorVariations::*)(const class std::vector<int, class std::allocator<int> > &)) &Pythia8::FlavorVariations::write, "C++: Pythia8::FlavorVariations::write(const class std::vector<int, class std::allocator<int> > &) --> std::string", pybind11::arg("breaks"));
		cl.def("parms", (class std::vector<double, class std::allocator<double> > (Pythia8::FlavorVariations::*)(double, double, double, double)) &Pythia8::FlavorVariations::parms, "C++: Pythia8::FlavorVariations::parms(double, double, double, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("xi"), pybind11::arg("rho"), pybind11::arg("x"), pybind11::arg("y"));
		cl.def("weight", (double (Pythia8::FlavorVariations::*)(const class std::vector<double, class std::allocator<double> > &, const class std::vector<int, class std::allocator<int> > &)) &Pythia8::FlavorVariations::weight, "C++: Pythia8::FlavorVariations::weight(const class std::vector<double, class std::allocator<double> > &, const class std::vector<int, class std::allocator<int> > &) --> double", pybind11::arg("parms"), pybind11::arg("breaks"));
	}
	{ // Pythia8::InputParser file:Pythia8Plugins/InputParser.h line:25
		pybind11::class_<Pythia8::InputParser, std::shared_ptr<Pythia8::InputParser>> cl(M("Pythia8"), "InputParser", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::InputParser(); } ), "doc" );
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0){ return new Pythia8::InputParser(a0); } ), "doc" , pybind11::arg("usage"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > const & a1){ return new Pythia8::InputParser(a0, a1); } ), "doc" , pybind11::arg("usage"), pybind11::arg("examples"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > const & a1, class std::basic_string<char> const & a2){ return new Pythia8::InputParser(a0, a1, a2); } ), "doc" , pybind11::arg("usage"), pybind11::arg("examples"), pybind11::arg("extra"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > const & a1, class std::basic_string<char> const & a2, class std::basic_ostream<char> * a3){ return new Pythia8::InputParser(a0, a1, a2, a3); } ), "doc" , pybind11::arg("usage"), pybind11::arg("examples"), pybind11::arg("extra"), pybind11::arg("stream"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > > const & a1, class std::basic_string<char> const & a2, class std::basic_ostream<char> * a3, class std::basic_string<char> const & a4){ return new Pythia8::InputParser(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("usage"), pybind11::arg("examples"), pybind11::arg("extra"), pybind11::arg("stream"), pybind11::arg("optName"));
		cl.def( pybind11::init<std::string, class std::vector<std::string, class std::allocator<std::string > >, std::string, std::ostream *, std::string, class std::set<std::string, struct std::less<std::string >, class std::allocator<std::string > >>(), pybind11::arg("usage"), pybind11::arg("examples"), pybind11::arg("extra"), pybind11::arg("stream"), pybind11::arg("optName"), pybind11::arg("aliases") );


		pybind11::enum_<Pythia8::InputParser::Status>(cl, "Status", pybind11::arithmetic(), "")
			.value("Valid", Pythia8::InputParser::Status::Valid)
			.value("Invalid", Pythia8::InputParser::Status::Invalid)
			.value("Help", Pythia8::InputParser::Status::Help)
			.export_values();

		cl.def("get", (bool (Pythia8::InputParser::*)(const std::string &)) &Pythia8::InputParser::get<bool>, "C++: Pythia8::InputParser::get(const std::string &) --> bool", pybind11::arg("optName"));
		cl.def("add", [](Pythia8::InputParser &o, const class std::basic_string<char> & a0, const class std::basic_string<char> & a1) -> bool { return o.add(a0, a1); }, "", pybind11::arg("optName"), pybind11::arg("defString"));
		cl.def("add", [](Pythia8::InputParser &o, const class std::basic_string<char> & a0, const class std::basic_string<char> & a1, const class std::basic_string<char> & a2) -> bool { return o.add(a0, a1, a2); }, "", pybind11::arg("optName"), pybind11::arg("defString"), pybind11::arg("helpText"));
		cl.def("add", (bool (Pythia8::InputParser::*)(const std::string &, const std::string &, const std::string &, class std::set<std::string, struct std::less<std::string >, class std::allocator<std::string > >)) &Pythia8::InputParser::add, "C++: Pythia8::InputParser::add(const std::string &, const std::string &, const std::string &, class std::set<std::string, struct std::less<std::string >, class std::allocator<std::string > >) --> bool", pybind11::arg("optName"), pybind11::arg("defString"), pybind11::arg("helpText"), pybind11::arg("aliases"));
		cl.def("require", [](Pythia8::InputParser &o, const class std::basic_string<char> & a0) -> bool { return o.require(a0); }, "", pybind11::arg("optName"));
		cl.def("require", [](Pythia8::InputParser &o, const class std::basic_string<char> & a0, const class std::basic_string<char> & a1) -> bool { return o.require(a0, a1); }, "", pybind11::arg("optName"), pybind11::arg("helpText"));
		cl.def("require", (bool (Pythia8::InputParser::*)(const std::string &, const std::string &, class std::set<std::string, struct std::less<std::string >, class std::allocator<std::string > >)) &Pythia8::InputParser::require, "C++: Pythia8::InputParser::require(const std::string &, const std::string &, class std::set<std::string, struct std::less<std::string >, class std::allocator<std::string > >) --> bool", pybind11::arg("optName"), pybind11::arg("helpText"), pybind11::arg("aliases"));
		cl.def("has", (bool (Pythia8::InputParser::*)(const std::string &) const) &Pythia8::InputParser::has, "C++: Pythia8::InputParser::has(const std::string &) const --> bool", pybind11::arg("optName"));
		cl.def("help", (const std::string (Pythia8::InputParser::*)() const) &Pythia8::InputParser::help, "C++: Pythia8::InputParser::help() const --> const std::string");
	}
}
#include <Pythia8/Basics.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/HIInfo.h>
#include <Pythia8/HeavyIons.h>
#include <Pythia8/Info.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Merging.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ParticleDecays.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/ShowerModel.h>
#include <Pythia8/SigmaProcess.h>
#include <Pythia8/UserHooks.h>
#include <Pythia8Plugins/CombineMatchingInput.h>
#include <Pythia8Plugins/JetMatching.h>
#include <Pythia8Plugins/LHAHelaconia.h>
#include <Pythia8Plugins/LHAMadgraph.h>
#include <istream>
#include <iterator>
#include <memory>
#include <sstream> // __str__
#include <string>
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

// Pythia8::LHAupHelaconia file:Pythia8Plugins/LHAHelaconia.h line:35
struct PyCallBack_Pythia8_LHAupHelaconia : public Pythia8::LHAupHelaconia {
	using Pythia8::LHAupHelaconia::LHAupHelaconia;

	bool setInit() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupHelaconia *>(this), "setInit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupHelaconia::setInit();
	}
	bool setEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupHelaconia *>(this), "setEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupHelaconia::setEvent(a0);
	}
	void newEventFile(const char * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupHelaconia *>(this), "newEventFile");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupHelaconia *>(this), "fileFound");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupHelaconia *>(this), "useExternal");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupHelaconia *>(this), "skipEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupHelaconia *>(this), "openLHEF");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupHelaconia *>(this), "closeLHEF");
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

// Pythia8::LHAupMadgraph file:Pythia8Plugins/LHAMadgraph.h line:63
struct PyCallBack_Pythia8_LHAupMadgraph : public Pythia8::LHAupMadgraph {
	using Pythia8::LHAupMadgraph::LHAupMadgraph;

	bool setInit() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupMadgraph *>(this), "setInit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupMadgraph::setInit();
	}
	bool setEvent(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupMadgraph *>(this), "setEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LHAupMadgraph::setEvent(a0);
	}
	void newEventFile(const char * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupMadgraph *>(this), "newEventFile");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupMadgraph *>(this), "fileFound");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupMadgraph *>(this), "useExternal");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupMadgraph *>(this), "skipEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupMadgraph *>(this), "openLHEF");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LHAupMadgraph *>(this), "closeLHEF");
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

void bind_Pythia8Plugins_LHAHelaconia(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::LHAupHelaconia file:Pythia8Plugins/LHAHelaconia.h line:35
		pybind11::class_<Pythia8::LHAupHelaconia, std::shared_ptr<Pythia8::LHAupHelaconia>, PyCallBack_Pythia8_LHAupHelaconia, Pythia8::LHAup> cl(M("Pythia8"), "LHAupHelaconia", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](class Pythia8::Pythia * a0){ return new Pythia8::LHAupHelaconia(a0); }, [](class Pythia8::Pythia * a0){ return new PyCallBack_Pythia8_LHAupHelaconia(a0); } ), "doc");
		cl.def( pybind11::init( [](class Pythia8::Pythia * a0, class std::basic_string<char> const & a1){ return new Pythia8::LHAupHelaconia(a0, a1); }, [](class Pythia8::Pythia * a0, class std::basic_string<char> const & a1){ return new PyCallBack_Pythia8_LHAupHelaconia(a0, a1); } ), "doc");
		cl.def( pybind11::init<class Pythia8::Pythia *, std::string, std::string>(), pybind11::arg("pythiaIn"), pybind11::arg("dirIn"), pybind11::arg("exeIn") );

		cl.def("readString", (bool (Pythia8::LHAupHelaconia::*)(std::string)) &Pythia8::LHAupHelaconia::readString, "C++: Pythia8::LHAupHelaconia::readString(std::string) --> bool", pybind11::arg("line"));
		cl.def("setEvents", (void (Pythia8::LHAupHelaconia::*)(int)) &Pythia8::LHAupHelaconia::setEvents, "C++: Pythia8::LHAupHelaconia::setEvents(int) --> void", pybind11::arg("eventsIn"));
		cl.def("setSeed", [](Pythia8::LHAupHelaconia &o, int const & a0) -> bool { return o.setSeed(a0); }, "", pybind11::arg("seedIn"));
		cl.def("setSeed", (bool (Pythia8::LHAupHelaconia::*)(int, int)) &Pythia8::LHAupHelaconia::setSeed, "C++: Pythia8::LHAupHelaconia::setSeed(int, int) --> bool", pybind11::arg("seedIn"), pybind11::arg("runsIn"));
		cl.def("setInit", (bool (Pythia8::LHAupHelaconia::*)()) &Pythia8::LHAupHelaconia::setInit, "C++: Pythia8::LHAupHelaconia::setInit() --> bool");
		cl.def("setEvent", [](Pythia8::LHAupHelaconia &o) -> bool { return o.setEvent(); }, "");
		cl.def("setEvent", (bool (Pythia8::LHAupHelaconia::*)(int)) &Pythia8::LHAupHelaconia::setEvent, "C++: Pythia8::LHAupHelaconia::setEvent(int) --> bool", pybind11::arg(""));
		cl.def("execute", (bool (Pythia8::LHAupHelaconia::*)(std::string)) &Pythia8::LHAupHelaconia::execute, "C++: Pythia8::LHAupHelaconia::execute(std::string) --> bool", pybind11::arg("line"));
		cl.def("run", [](Pythia8::LHAupHelaconia &o, int const & a0) -> bool { return o.run(a0); }, "", pybind11::arg("eventsIn"));
		cl.def("run", (bool (Pythia8::LHAupHelaconia::*)(int, int)) &Pythia8::LHAupHelaconia::run, "C++: Pythia8::LHAupHelaconia::run(int, int) --> bool", pybind11::arg("eventsIn"), pybind11::arg("seedIn"));
		cl.def("reader", (bool (Pythia8::LHAupHelaconia::*)(bool)) &Pythia8::LHAupHelaconia::reader, "C++: Pythia8::LHAupHelaconia::reader(bool) --> bool", pybind11::arg("init"));
		cl.def("convert", (int (Pythia8::LHAupHelaconia::*)(int)) &Pythia8::LHAupHelaconia::convert, "C++: Pythia8::LHAupHelaconia::convert(int) --> int", pybind11::arg("idIn"));
		cl.def("errorMsg", (void (Pythia8::LHAupHelaconia::*)(std::string)) &Pythia8::LHAupHelaconia::errorMsg, "C++: Pythia8::LHAupHelaconia::errorMsg(std::string) --> void", pybind11::arg("messageIn"));
	}
	{ // Pythia8::LHAupMadgraph file:Pythia8Plugins/LHAMadgraph.h line:63
		pybind11::class_<Pythia8::LHAupMadgraph, std::shared_ptr<Pythia8::LHAupMadgraph>, PyCallBack_Pythia8_LHAupMadgraph, Pythia8::LHAup> cl(M("Pythia8"), "LHAupMadgraph", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](class Pythia8::Pythia * a0){ return new Pythia8::LHAupMadgraph(a0); }, [](class Pythia8::Pythia * a0){ return new PyCallBack_Pythia8_LHAupMadgraph(a0); } ), "doc");
		cl.def( pybind11::init( [](class Pythia8::Pythia * a0, bool const & a1){ return new Pythia8::LHAupMadgraph(a0, a1); }, [](class Pythia8::Pythia * a0, bool const & a1){ return new PyCallBack_Pythia8_LHAupMadgraph(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](class Pythia8::Pythia * a0, bool const & a1, class std::basic_string<char> const & a2){ return new Pythia8::LHAupMadgraph(a0, a1, a2); }, [](class Pythia8::Pythia * a0, bool const & a1, class std::basic_string<char> const & a2){ return new PyCallBack_Pythia8_LHAupMadgraph(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<class Pythia8::Pythia *, bool, std::string, std::string>(), pybind11::arg("pythiaIn"), pybind11::arg("matchIn"), pybind11::arg("dirIn"), pybind11::arg("exeIn") );


		pybind11::enum_<Pythia8::LHAupMadgraph::Stage>(cl, "Stage", pybind11::arithmetic(), "")
			.value("Auto", Pythia8::LHAupMadgraph::Stage::Auto)
			.value("Configure", Pythia8::LHAupMadgraph::Stage::Configure)
			.value("Generate", Pythia8::LHAupMadgraph::Stage::Generate)
			.value("Launch", Pythia8::LHAupMadgraph::Stage::Launch)
			.export_values();

		cl.def("readString", [](Pythia8::LHAupMadgraph &o, class std::basic_string<char> const & a0) -> bool { return o.readString(a0); }, "", pybind11::arg("line"));
		cl.def("readString", (bool (Pythia8::LHAupMadgraph::*)(std::string, enum Pythia8::LHAupMadgraph::Stage)) &Pythia8::LHAupMadgraph::readString, "C++: Pythia8::LHAupMadgraph::readString(std::string, enum Pythia8::LHAupMadgraph::Stage) --> bool", pybind11::arg("line"), pybind11::arg("stage"));
		cl.def("addCard", (void (Pythia8::LHAupMadgraph::*)(std::string, std::string)) &Pythia8::LHAupMadgraph::addCard, "C++: Pythia8::LHAupMadgraph::addCard(std::string, std::string) --> void", pybind11::arg("src"), pybind11::arg("dst"));
		cl.def("setEvents", (void (Pythia8::LHAupMadgraph::*)(int)) &Pythia8::LHAupMadgraph::setEvents, "C++: Pythia8::LHAupMadgraph::setEvents(int) --> void", pybind11::arg("eventsIn"));
		cl.def("setSeed", [](Pythia8::LHAupMadgraph &o, int const & a0) -> bool { return o.setSeed(a0); }, "", pybind11::arg("seedIn"));
		cl.def("setSeed", (bool (Pythia8::LHAupMadgraph::*)(int, int)) &Pythia8::LHAupMadgraph::setSeed, "C++: Pythia8::LHAupMadgraph::setSeed(int, int) --> bool", pybind11::arg("seedIn"), pybind11::arg("runsIn"));
		cl.def("setJets", (void (Pythia8::LHAupMadgraph::*)(int)) &Pythia8::LHAupMadgraph::setJets, "C++: Pythia8::LHAupMadgraph::setJets(int) --> void", pybind11::arg("jetsIn"));
		cl.def("setInit", (bool (Pythia8::LHAupMadgraph::*)()) &Pythia8::LHAupMadgraph::setInit, "C++: Pythia8::LHAupMadgraph::setInit() --> bool");
		cl.def("setEvent", [](Pythia8::LHAupMadgraph &o) -> bool { return o.setEvent(); }, "");
		cl.def("setEvent", (bool (Pythia8::LHAupMadgraph::*)(int)) &Pythia8::LHAupMadgraph::setEvent, "C++: Pythia8::LHAupMadgraph::setEvent(int) --> bool", pybind11::arg(""));
		cl.def("execute", (bool (Pythia8::LHAupMadgraph::*)(std::string)) &Pythia8::LHAupMadgraph::execute, "C++: Pythia8::LHAupMadgraph::execute(std::string) --> bool", pybind11::arg("line"));
		cl.def("configure", (bool (Pythia8::LHAupMadgraph::*)()) &Pythia8::LHAupMadgraph::configure, "C++: Pythia8::LHAupMadgraph::configure() --> bool");
		cl.def("generate", (bool (Pythia8::LHAupMadgraph::*)()) &Pythia8::LHAupMadgraph::generate, "C++: Pythia8::LHAupMadgraph::generate() --> bool");
		cl.def("launch", (bool (Pythia8::LHAupMadgraph::*)()) &Pythia8::LHAupMadgraph::launch, "C++: Pythia8::LHAupMadgraph::launch() --> bool");
		cl.def("run", [](Pythia8::LHAupMadgraph &o, int const & a0) -> bool { return o.run(a0); }, "", pybind11::arg("eventsIn"));
		cl.def("run", (bool (Pythia8::LHAupMadgraph::*)(int, int)) &Pythia8::LHAupMadgraph::run, "C++: Pythia8::LHAupMadgraph::run(int, int) --> bool", pybind11::arg("eventsIn"), pybind11::arg("seedIn"));
		cl.def("reader", (bool (Pythia8::LHAupMadgraph::*)(bool)) &Pythia8::LHAupMadgraph::reader, "C++: Pythia8::LHAupMadgraph::reader(bool) --> bool", pybind11::arg("init"));
		cl.def("errorMsg", (void (Pythia8::LHAupMadgraph::*)(std::string)) &Pythia8::LHAupMadgraph::errorMsg, "C++: Pythia8::LHAupMadgraph::errorMsg(std::string) --> void", pybind11::arg("messageIn"));
	}
}
#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/GammaKinematics.h>
#include <Pythia8/HadronLevel.h>
#include <Pythia8/Info.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/NucleonExcitations.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ParticleDecays.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/RHadrons.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/SLHAinterface.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaProcess.h>
#include <Pythia8/StringFragmentation.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/TimeShower.h>
#include <Pythia8/UserHooks.h>
#include <Pythia8Plugins/ProgressLog.h>
#include <Pythia8Plugins/PythiaCascade.h>
#include <Pythia8Plugins/ResonanceDecayFilterHook.h>
#include <Pythia8Plugins/SetLHEDecayProductHook.h>
#include <Pythia8Plugins/Visualisation.h>
#include <functional>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
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

// Pythia8::ResonanceDecayFilterHook file:Pythia8Plugins/ResonanceDecayFilterHook.h line:49
struct PyCallBack_Pythia8_ResonanceDecayFilterHook : public Pythia8::ResonanceDecayFilterHook {
	using Pythia8::ResonanceDecayFilterHook::ResonanceDecayFilterHook;

	bool canVetoResonanceDecays() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoResonanceDecays");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceDecayFilterHook::canVetoResonanceDecays();
	}
	bool doVetoResonanceDecays(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoResonanceDecays");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceDecayFilterHook::doVetoResonanceDecays(a0);
	}
	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ResonanceDecayFilterHook::initAfterBeams();
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "biasedSelectionWeight");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoProcessLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoProcessLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doSetLowEnergySigma");
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
	bool canVetoPT() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "numberVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ResonanceDecayFilterHook *>(this), "onStat");
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

// Pythia8::SetLHEDecayProductHook file:Pythia8Plugins/SetLHEDecayProductHook.h line:24
struct PyCallBack_Pythia8_SetLHEDecayProductHook : public Pythia8::SetLHEDecayProductHook {
	using Pythia8::SetLHEDecayProductHook::SetLHEDecayProductHook;

	bool canVetoProcessLevel() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SetLHEDecayProductHook::canVetoProcessLevel();
	}
	bool doVetoProcessLevel(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoProcessLevel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SetLHEDecayProductHook::doVetoProcessLevel(a0);
	}
	bool initAfterBeams() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "initAfterBeams");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SetLHEDecayProductHook::initAfterBeams();
	}
	bool canModifySigma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canModifySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "multiplySigmaBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canBiasSelection");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "biasSelectionBy");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "biasedSelectionWeight");
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
	bool canSetLowEnergySigma(int a0, int a1) const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doSetLowEnergySigma");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoResonanceDecays");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "scaleVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoPT");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "numberVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "numberVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoMPIStep");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoPartonLevelEarly");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "retryPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoPartonLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canSetResonanceScale");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "scaleResonance");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoISREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoFSREmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoMPIEmission");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doReconnectResonanceSystems");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "setStringEnds");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doChangeFragPar");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoFragmentation");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoFinalTwo");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doVetoAfterHadronization");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "canSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "doSetImpactParameter");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "onEndHadronLevel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SetLHEDecayProductHook *>(this), "onStat");
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

void bind_Pythia8Plugins_ProgressLog(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::ProgressLog file:Pythia8Plugins/ProgressLog.h line:27
		pybind11::class_<Pythia8::ProgressLog, std::shared_ptr<Pythia8::ProgressLog>> cl(M("Pythia8"), "ProgressLog", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](long const & a0){ return new Pythia8::ProgressLog(a0); } ), "doc" , pybind11::arg("NIn"));
		cl.def( pybind11::init<long, int>(), pybind11::arg("NIn"), pybind11::arg("maxsec") );

		cl.def("__call__", (void (Pythia8::ProgressLog::*)(int)) &Pythia8::ProgressLog::operator(), "C++: Pythia8::ProgressLog::operator()(int) --> void", pybind11::arg("cnt"));
		cl.def("__call__", (void (Pythia8::ProgressLog::*)()) &Pythia8::ProgressLog::operator(), "C++: Pythia8::ProgressLog::operator()() --> void");
		cl.def("tick", (void (Pythia8::ProgressLog::*)(long, long)) &Pythia8::ProgressLog::tick, "C++: Pythia8::ProgressLog::tick(long, long) --> void", pybind11::arg("i"), pybind11::arg("n"));
		cl.def("fclock", (double (Pythia8::ProgressLog::*)()) &Pythia8::ProgressLog::fclock, "C++: Pythia8::ProgressLog::fclock() --> double");
		cl.def("statusTime", (bool (Pythia8::ProgressLog::*)(long, long) const) &Pythia8::ProgressLog::statusTime, "C++: Pythia8::ProgressLog::statusTime(long, long) const --> bool", pybind11::arg("i"), pybind11::arg("n"));
		cl.def("init", (void (Pythia8::ProgressLog::*)(long)) &Pythia8::ProgressLog::init, "C++: Pythia8::ProgressLog::init(long) --> void", pybind11::arg("n"));
	}
	{ // Pythia8::PythiaCascade file:Pythia8Plugins/PythiaCascade.h line:61
		pybind11::class_<Pythia8::PythiaCascade, std::shared_ptr<Pythia8::PythiaCascade>> cl(M("Pythia8"), "PythiaCascade", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PythiaCascade(); } ) );
		cl.def("init", [](Pythia8::PythiaCascade &o) -> bool { return o.init(); }, "");
		cl.def("init", [](Pythia8::PythiaCascade &o, double const & a0) -> bool { return o.init(a0); }, "", pybind11::arg("eMaxIn"));
		cl.def("init", [](Pythia8::PythiaCascade &o, double const & a0, bool const & a1) -> bool { return o.init(a0, a1); }, "", pybind11::arg("eMaxIn"), pybind11::arg("listFinalIn"));
		cl.def("init", [](Pythia8::PythiaCascade &o, double const & a0, bool const & a1, bool const & a2) -> bool { return o.init(a0, a1, a2); }, "", pybind11::arg("eMaxIn"), pybind11::arg("listFinalIn"), pybind11::arg("rapidDecaysIn"));
		cl.def("init", [](Pythia8::PythiaCascade &o, double const & a0, bool const & a1, bool const & a2, double const & a3) -> bool { return o.init(a0, a1, a2, a3); }, "", pybind11::arg("eMaxIn"), pybind11::arg("listFinalIn"), pybind11::arg("rapidDecaysIn"), pybind11::arg("smallTau0In"));
		cl.def("init", [](Pythia8::PythiaCascade &o, double const & a0, bool const & a1, bool const & a2, double const & a3, int const & a4) -> bool { return o.init(a0, a1, a2, a3, a4); }, "", pybind11::arg("eMaxIn"), pybind11::arg("listFinalIn"), pybind11::arg("rapidDecaysIn"), pybind11::arg("smallTau0In"), pybind11::arg("reuseMPI"));
		cl.def("init", (bool (Pythia8::PythiaCascade::*)(double, bool, bool, double, int, std::string)) &Pythia8::PythiaCascade::init, "C++: Pythia8::PythiaCascade::init(double, bool, bool, double, int, std::string) --> bool", pybind11::arg("eMaxIn"), pybind11::arg("listFinalIn"), pybind11::arg("rapidDecaysIn"), pybind11::arg("smallTau0In"), pybind11::arg("reuseMPI"), pybind11::arg("initFile"));
		cl.def("nCollAvg", (double (Pythia8::PythiaCascade::*)(int)) &Pythia8::PythiaCascade::nCollAvg, "C++: Pythia8::PythiaCascade::nCollAvg(int) --> double", pybind11::arg("A"));
		cl.def("sigmaSetuphN", (bool (Pythia8::PythiaCascade::*)(int, class Pythia8::Vec4, double)) &Pythia8::PythiaCascade::sigmaSetuphN, "C++: Pythia8::PythiaCascade::sigmaSetuphN(int, class Pythia8::Vec4, double) --> bool", pybind11::arg("idNowIn"), pybind11::arg("pNowIn"), pybind11::arg("mNowIn"));
		cl.def("sigmahA", (double (Pythia8::PythiaCascade::*)(int)) &Pythia8::PythiaCascade::sigmahA, "C++: Pythia8::PythiaCascade::sigmahA(int) --> double", pybind11::arg("A"));
		cl.def("nextColl", [](Pythia8::PythiaCascade &o, int const & a0, int const & a1) -> Pythia8::Event & { return o.nextColl(a0, a1); }, "", pybind11::return_value_policy::reference, pybind11::arg("Znow"), pybind11::arg("Anow"));
		cl.def("nextColl", (class Pythia8::Event & (Pythia8::PythiaCascade::*)(int, int, class Pythia8::Vec4)) &Pythia8::PythiaCascade::nextColl, "C++: Pythia8::PythiaCascade::nextColl(int, int, class Pythia8::Vec4) --> class Pythia8::Event &", pybind11::return_value_policy::reference, pybind11::arg("Znow"), pybind11::arg("Anow"), pybind11::arg("vNow"));
		cl.def("nextDecay", [](Pythia8::PythiaCascade &o, int const & a0, class Pythia8::Vec4 const & a1, double const & a2) -> Pythia8::Event & { return o.nextDecay(a0, a1, a2); }, "", pybind11::return_value_policy::reference, pybind11::arg("idNowIn"), pybind11::arg("pNowIn"), pybind11::arg("mNowIn"));
		cl.def("nextDecay", (class Pythia8::Event & (Pythia8::PythiaCascade::*)(int, class Pythia8::Vec4, double, class Pythia8::Vec4)) &Pythia8::PythiaCascade::nextDecay, "C++: Pythia8::PythiaCascade::nextDecay(int, class Pythia8::Vec4, double, class Pythia8::Vec4) --> class Pythia8::Event &", pybind11::return_value_policy::reference, pybind11::arg("idNowIn"), pybind11::arg("pNowIn"), pybind11::arg("mNowIn"), pybind11::arg("vNow"));
		cl.def("compress", (void (Pythia8::PythiaCascade::*)()) &Pythia8::PythiaCascade::compress, "C++: Pythia8::PythiaCascade::compress() --> void");
		cl.def("stat", (void (Pythia8::PythiaCascade::*)()) &Pythia8::PythiaCascade::stat, "C++: Pythia8::PythiaCascade::stat() --> void");
		cl.def("particleData", (class Pythia8::ParticleData & (Pythia8::PythiaCascade::*)()) &Pythia8::PythiaCascade::particleData, "C++: Pythia8::PythiaCascade::particleData() --> class Pythia8::ParticleData &", pybind11::return_value_policy::reference);
		cl.def("rndm", (class Pythia8::Rndm & (Pythia8::PythiaCascade::*)()) &Pythia8::PythiaCascade::rndm, "C++: Pythia8::PythiaCascade::rndm() --> class Pythia8::Rndm &", pybind11::return_value_policy::reference);
	}
	{ // Pythia8::ResonanceDecayFilterHook file:Pythia8Plugins/ResonanceDecayFilterHook.h line:49
		pybind11::class_<Pythia8::ResonanceDecayFilterHook, std::shared_ptr<Pythia8::ResonanceDecayFilterHook>, PyCallBack_Pythia8_ResonanceDecayFilterHook, Pythia8::UserHooks> cl(M("Pythia8"), "ResonanceDecayFilterHook", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Settings &>(), pybind11::arg("settings") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_ResonanceDecayFilterHook const &o){ return new PyCallBack_Pythia8_ResonanceDecayFilterHook(o); } ) );
		cl.def( pybind11::init( [](Pythia8::ResonanceDecayFilterHook const &o){ return new Pythia8::ResonanceDecayFilterHook(o); } ) );
		cl.def("canVetoResonanceDecays", (bool (Pythia8::ResonanceDecayFilterHook::*)()) &Pythia8::ResonanceDecayFilterHook::canVetoResonanceDecays, "C++: Pythia8::ResonanceDecayFilterHook::canVetoResonanceDecays() --> bool");
		cl.def("doVetoResonanceDecays", (bool (Pythia8::ResonanceDecayFilterHook::*)(class Pythia8::Event &)) &Pythia8::ResonanceDecayFilterHook::doVetoResonanceDecays, "C++: Pythia8::ResonanceDecayFilterHook::doVetoResonanceDecays(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("initAfterBeams", (bool (Pythia8::ResonanceDecayFilterHook::*)()) &Pythia8::ResonanceDecayFilterHook::initAfterBeams, "C++: Pythia8::ResonanceDecayFilterHook::initAfterBeams() --> bool");
		cl.def("checkVetoResonanceDecays", (bool (Pythia8::ResonanceDecayFilterHook::*)(const class Pythia8::Event &)) &Pythia8::ResonanceDecayFilterHook::checkVetoResonanceDecays, "C++: Pythia8::ResonanceDecayFilterHook::checkVetoResonanceDecays(const class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("returnCounter", (unsigned long (Pythia8::ResonanceDecayFilterHook::*)()) &Pythia8::ResonanceDecayFilterHook::returnCounter, "C++: Pythia8::ResonanceDecayFilterHook::returnCounter() --> unsigned long");
		cl.def("assign", (class Pythia8::ResonanceDecayFilterHook & (Pythia8::ResonanceDecayFilterHook::*)(const class Pythia8::ResonanceDecayFilterHook &)) &Pythia8::ResonanceDecayFilterHook::operator=, "C++: Pythia8::ResonanceDecayFilterHook::operator=(const class Pythia8::ResonanceDecayFilterHook &) --> class Pythia8::ResonanceDecayFilterHook &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SetLHEDecayProductHook file:Pythia8Plugins/SetLHEDecayProductHook.h line:24
		pybind11::class_<Pythia8::SetLHEDecayProductHook, std::shared_ptr<Pythia8::SetLHEDecayProductHook>, PyCallBack_Pythia8_SetLHEDecayProductHook, Pythia8::UserHooks> cl(M("Pythia8"), "SetLHEDecayProductHook", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Settings &, const class Pythia8::ParticleData *>(), pybind11::arg("settings"), pybind11::arg("pdtPtrIn") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_SetLHEDecayProductHook const &o){ return new PyCallBack_Pythia8_SetLHEDecayProductHook(o); } ) );
		cl.def( pybind11::init( [](Pythia8::SetLHEDecayProductHook const &o){ return new Pythia8::SetLHEDecayProductHook(o); } ) );
		cl.def("canVetoProcessLevel", (bool (Pythia8::SetLHEDecayProductHook::*)()) &Pythia8::SetLHEDecayProductHook::canVetoProcessLevel, "C++: Pythia8::SetLHEDecayProductHook::canVetoProcessLevel() --> bool");
		cl.def("doVetoProcessLevel", (bool (Pythia8::SetLHEDecayProductHook::*)(class Pythia8::Event &)) &Pythia8::SetLHEDecayProductHook::doVetoProcessLevel, "C++: Pythia8::SetLHEDecayProductHook::doVetoProcessLevel(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("initAfterBeams", (bool (Pythia8::SetLHEDecayProductHook::*)()) &Pythia8::SetLHEDecayProductHook::initAfterBeams, "C++: Pythia8::SetLHEDecayProductHook::initAfterBeams() --> bool");
		cl.def("checkVetoProcessLevel", (bool (Pythia8::SetLHEDecayProductHook::*)(class Pythia8::Event &)) &Pythia8::SetLHEDecayProductHook::checkVetoProcessLevel, "C++: Pythia8::SetLHEDecayProductHook::checkVetoProcessLevel(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("returnCounter", (unsigned long (Pythia8::SetLHEDecayProductHook::*)()) &Pythia8::SetLHEDecayProductHook::returnCounter, "C++: Pythia8::SetLHEDecayProductHook::returnCounter() --> unsigned long");
		cl.def("assign", (class Pythia8::SetLHEDecayProductHook & (Pythia8::SetLHEDecayProductHook::*)(const class Pythia8::SetLHEDecayProductHook &)) &Pythia8::SetLHEDecayProductHook::operator=, "C++: Pythia8::SetLHEDecayProductHook::operator=(const class Pythia8::SetLHEDecayProductHook &) --> class Pythia8::SetLHEDecayProductHook &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	// Pythia8::py_status(int) file:Pythia8Plugins/Visualisation.h line:21
	M("Pythia8").def("py_status", (std::string (*)(int)) &Pythia8::py_status, "C++: Pythia8::py_status(int) --> std::string", pybind11::arg("stAbs"));

	// Pythia8::makeArrow(class std::map<struct std::pair<std::string, std::string >, std::string, struct std::less<struct std::pair<std::string, std::string > >, class std::allocator<struct std::pair<const struct std::pair<std::string, std::string >, std::string > > > *, std::string, std::string) file:Pythia8Plugins/Visualisation.h line:37
	M("Pythia8").def("makeArrow", (void (*)(class std::map<struct std::pair<std::string, std::string >, std::string, struct std::less<struct std::pair<std::string, std::string > >, class std::allocator<struct std::pair<const struct std::pair<std::string, std::string >, std::string > > > *, std::string, std::string)) &Pythia8::makeArrow, "C++: Pythia8::makeArrow(class std::map<struct std::pair<std::string, std::string >, std::string, struct std::less<struct std::pair<std::string, std::string > >, class std::allocator<struct std::pair<const struct std::pair<std::string, std::string >, std::string > > > *, std::string, std::string) --> void", pybind11::arg("arrows"), pybind11::arg("identParent"), pybind11::arg("identChild"));

	// Pythia8::printEvent(class Pythia8::Event &, std::string) file:Pythia8Plugins/Visualisation.h line:52
	M("Pythia8").def("printEvent", [](class Pythia8::Event & a0) -> void { return Pythia8::printEvent(a0); }, "", pybind11::arg("evt"));
	M("Pythia8").def("printEvent", (void (*)(class Pythia8::Event &, std::string)) &Pythia8::printEvent, "C++: Pythia8::printEvent(class Pythia8::Event &, std::string) --> void", pybind11::arg("evt"), pybind11::arg("fileName"));

}
