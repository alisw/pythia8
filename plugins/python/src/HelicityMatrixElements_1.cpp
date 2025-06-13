#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/HelicityBasics.h>
#include <Pythia8/HelicityMatrixElements.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/SusyCouplings.h>
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

// Pythia8::HMETau2Meson file:Pythia8/HelicityMatrixElements.h line:344
struct PyCallBack_Pythia8_HMETau2Meson : public Pythia8::HMETau2Meson {
	using Pythia8::HMETau2Meson::HMETau2Meson;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2Meson::initConstants();
	}
	void initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "initHadronicCurrent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2Meson::initHadronicCurrent(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETauDecay::calculateME(a0);
	}
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETauDecay::decayWeightMax(a0);
	}
	void calculateResonanceWeights(class std::vector<double, class std::allocator<double> > & a0, class std::vector<double, class std::allocator<double> > & a1, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "calculateResonanceWeights");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::calculateResonanceWeights(a0, a1, a2);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "initPointers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initPointers(a0, a1, a2);
	}
	class Pythia8::HelicityMatrixElement * initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "initChannel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::HelicityMatrixElement *>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::HelicityMatrixElement *> caster;
				return pybind11::detail::cast_ref<class Pythia8::HelicityMatrixElement *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::HelicityMatrixElement *>(std::move(o));
		}
		return HelicityMatrixElement::initChannel(a0);
	}
	double decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "decayWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeight(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "calculateD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateD(a0);
	}
	void calculateRho(unsigned int a0, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "calculateRho");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateRho(a0, a1);
	}
	struct std::complex<double> breitWigner(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "breitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::breitWigner(a0, a1, a2);
	}
	struct std::complex<double> sBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "sBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::sBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> pBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "pBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::pBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> dBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2Meson *>(this), "dBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::dBreitWigner(a0, a1, a2, a3, a4);
	}
};

// Pythia8::HMETau2TwoLeptons file:Pythia8/HelicityMatrixElements.h line:358
struct PyCallBack_Pythia8_HMETau2TwoLeptons : public Pythia8::HMETau2TwoLeptons {
	using Pythia8::HMETau2TwoLeptons::HMETau2TwoLeptons;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2TwoLeptons::initConstants();
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2TwoLeptons::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2TwoLeptons::calculateME(a0);
	}
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETauDecay::decayWeightMax(a0);
	}
	void initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "initHadronicCurrent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::initHadronicCurrent(a0);
	}
	void calculateResonanceWeights(class std::vector<double, class std::allocator<double> > & a0, class std::vector<double, class std::allocator<double> > & a1, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "calculateResonanceWeights");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::calculateResonanceWeights(a0, a1, a2);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "initPointers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initPointers(a0, a1, a2);
	}
	class Pythia8::HelicityMatrixElement * initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "initChannel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::HelicityMatrixElement *>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::HelicityMatrixElement *> caster;
				return pybind11::detail::cast_ref<class Pythia8::HelicityMatrixElement *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::HelicityMatrixElement *>(std::move(o));
		}
		return HelicityMatrixElement::initChannel(a0);
	}
	double decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "decayWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeight(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "calculateD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateD(a0);
	}
	void calculateRho(unsigned int a0, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "calculateRho");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateRho(a0, a1);
	}
	struct std::complex<double> breitWigner(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "breitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::breitWigner(a0, a1, a2);
	}
	struct std::complex<double> sBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "sBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::sBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> pBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "pBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::pBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> dBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoLeptons *>(this), "dBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::dBreitWigner(a0, a1, a2, a3, a4);
	}
};

// Pythia8::HMETau2TwoMesonsViaVector file:Pythia8/HelicityMatrixElements.h line:375
struct PyCallBack_Pythia8_HMETau2TwoMesonsViaVector : public Pythia8::HMETau2TwoMesonsViaVector {
	using Pythia8::HMETau2TwoMesonsViaVector::HMETau2TwoMesonsViaVector;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2TwoMesonsViaVector::initConstants();
	}
	void initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "initHadronicCurrent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2TwoMesonsViaVector::initHadronicCurrent(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETauDecay::calculateME(a0);
	}
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETauDecay::decayWeightMax(a0);
	}
	void calculateResonanceWeights(class std::vector<double, class std::allocator<double> > & a0, class std::vector<double, class std::allocator<double> > & a1, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "calculateResonanceWeights");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::calculateResonanceWeights(a0, a1, a2);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "initPointers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initPointers(a0, a1, a2);
	}
	class Pythia8::HelicityMatrixElement * initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "initChannel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::HelicityMatrixElement *>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::HelicityMatrixElement *> caster;
				return pybind11::detail::cast_ref<class Pythia8::HelicityMatrixElement *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::HelicityMatrixElement *>(std::move(o));
		}
		return HelicityMatrixElement::initChannel(a0);
	}
	double decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "decayWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeight(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "calculateD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateD(a0);
	}
	void calculateRho(unsigned int a0, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "calculateRho");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateRho(a0, a1);
	}
	struct std::complex<double> breitWigner(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "breitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::breitWigner(a0, a1, a2);
	}
	struct std::complex<double> sBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "sBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::sBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> pBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "pBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::pBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> dBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVector *>(this), "dBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::dBreitWigner(a0, a1, a2, a3, a4);
	}
};

// Pythia8::HMETau2TwoMesonsViaVectorScalar file:Pythia8/HelicityMatrixElements.h line:396
struct PyCallBack_Pythia8_HMETau2TwoMesonsViaVectorScalar : public Pythia8::HMETau2TwoMesonsViaVectorScalar {
	using Pythia8::HMETau2TwoMesonsViaVectorScalar::HMETau2TwoMesonsViaVectorScalar;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2TwoMesonsViaVectorScalar::initConstants();
	}
	void initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "initHadronicCurrent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2TwoMesonsViaVectorScalar::initHadronicCurrent(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETauDecay::calculateME(a0);
	}
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETauDecay::decayWeightMax(a0);
	}
	void calculateResonanceWeights(class std::vector<double, class std::allocator<double> > & a0, class std::vector<double, class std::allocator<double> > & a1, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "calculateResonanceWeights");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::calculateResonanceWeights(a0, a1, a2);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "initPointers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initPointers(a0, a1, a2);
	}
	class Pythia8::HelicityMatrixElement * initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "initChannel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::HelicityMatrixElement *>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::HelicityMatrixElement *> caster;
				return pybind11::detail::cast_ref<class Pythia8::HelicityMatrixElement *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::HelicityMatrixElement *>(std::move(o));
		}
		return HelicityMatrixElement::initChannel(a0);
	}
	double decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "decayWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeight(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "calculateD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateD(a0);
	}
	void calculateRho(unsigned int a0, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "calculateRho");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateRho(a0, a1);
	}
	struct std::complex<double> breitWigner(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "breitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::breitWigner(a0, a1, a2);
	}
	struct std::complex<double> sBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "sBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::sBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> pBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "pBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::pBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> dBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2TwoMesonsViaVectorScalar *>(this), "dBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::dBreitWigner(a0, a1, a2, a3, a4);
	}
};

// Pythia8::HMETau2ThreeMesons file:Pythia8/HelicityMatrixElements.h line:421
struct PyCallBack_Pythia8_HMETau2ThreeMesons : public Pythia8::HMETau2ThreeMesons {
	using Pythia8::HMETau2ThreeMesons::HMETau2ThreeMesons;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initConstants();
	}
	void initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "initHadronicCurrent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initHadronicCurrent(a0);
	}
	void initMode() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "initMode");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initMode();
	}
	void initResonances() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "initResonances");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initResonances();
	}
	void initMomenta(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "initMomenta");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initMomenta(a0);
	}
	struct std::complex<double> F1() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "F1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F1();
	}
	struct std::complex<double> F2() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "F2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F2();
	}
	struct std::complex<double> F3() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "F3");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F3();
	}
	struct std::complex<double> F4() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "F4");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F4();
	}
	double a1PhaseSpace(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "a1PhaseSpace");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETau2ThreeMesons::a1PhaseSpace(a0);
	}
	struct std::complex<double> a1BreitWigner(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "a1BreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::a1BreitWigner(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETauDecay::calculateME(a0);
	}
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETauDecay::decayWeightMax(a0);
	}
	void calculateResonanceWeights(class std::vector<double, class std::allocator<double> > & a0, class std::vector<double, class std::allocator<double> > & a1, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "calculateResonanceWeights");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::calculateResonanceWeights(a0, a1, a2);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "initPointers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initPointers(a0, a1, a2);
	}
	class Pythia8::HelicityMatrixElement * initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "initChannel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::HelicityMatrixElement *>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::HelicityMatrixElement *> caster;
				return pybind11::detail::cast_ref<class Pythia8::HelicityMatrixElement *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::HelicityMatrixElement *>(std::move(o));
		}
		return HelicityMatrixElement::initChannel(a0);
	}
	double decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "decayWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeight(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "calculateD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateD(a0);
	}
	void calculateRho(unsigned int a0, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "calculateRho");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateRho(a0, a1);
	}
	struct std::complex<double> breitWigner(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "breitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::breitWigner(a0, a1, a2);
	}
	struct std::complex<double> sBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "sBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::sBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> pBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "pBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::pBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> dBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesons *>(this), "dBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::dBreitWigner(a0, a1, a2, a3, a4);
	}
};

// Pythia8::HMETau2ThreePions file:Pythia8/HelicityMatrixElements.h line:473
struct PyCallBack_Pythia8_HMETau2ThreePions : public Pythia8::HMETau2ThreePions {
	using Pythia8::HMETau2ThreePions::HMETau2ThreePions;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initConstants();
	}
	void initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "initHadronicCurrent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initHadronicCurrent(a0);
	}
	void initMode() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "initMode");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initMode();
	}
	void initResonances() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "initResonances");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initResonances();
	}
	void initMomenta(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "initMomenta");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initMomenta(a0);
	}
	struct std::complex<double> F1() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "F1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F1();
	}
	struct std::complex<double> F2() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "F2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F2();
	}
	struct std::complex<double> F3() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "F3");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F3();
	}
	struct std::complex<double> F4() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "F4");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F4();
	}
	double a1PhaseSpace(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "a1PhaseSpace");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETau2ThreeMesons::a1PhaseSpace(a0);
	}
	struct std::complex<double> a1BreitWigner(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "a1BreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::a1BreitWigner(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETauDecay::calculateME(a0);
	}
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETauDecay::decayWeightMax(a0);
	}
	void calculateResonanceWeights(class std::vector<double, class std::allocator<double> > & a0, class std::vector<double, class std::allocator<double> > & a1, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "calculateResonanceWeights");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::calculateResonanceWeights(a0, a1, a2);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "initPointers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initPointers(a0, a1, a2);
	}
	class Pythia8::HelicityMatrixElement * initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "initChannel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::HelicityMatrixElement *>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::HelicityMatrixElement *> caster;
				return pybind11::detail::cast_ref<class Pythia8::HelicityMatrixElement *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::HelicityMatrixElement *>(std::move(o));
		}
		return HelicityMatrixElement::initChannel(a0);
	}
	double decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "decayWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeight(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "calculateD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateD(a0);
	}
	void calculateRho(unsigned int a0, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "calculateRho");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateRho(a0, a1);
	}
	struct std::complex<double> breitWigner(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "breitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::breitWigner(a0, a1, a2);
	}
	struct std::complex<double> sBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "sBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::sBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> pBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "pBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::pBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> dBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreePions *>(this), "dBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::dBreitWigner(a0, a1, a2, a3, a4);
	}
};

// Pythia8::HMETau2ThreeMesonsWithKaons file:Pythia8/HelicityMatrixElements.h line:506
struct PyCallBack_Pythia8_HMETau2ThreeMesonsWithKaons : public Pythia8::HMETau2ThreeMesonsWithKaons {
	using Pythia8::HMETau2ThreeMesonsWithKaons::HMETau2ThreeMesonsWithKaons;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initConstants();
	}
	void initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "initHadronicCurrent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initHadronicCurrent(a0);
	}
	void initMode() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "initMode");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initMode();
	}
	void initResonances() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "initResonances");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initResonances();
	}
	void initMomenta(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "initMomenta");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initMomenta(a0);
	}
	struct std::complex<double> F1() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "F1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F1();
	}
	struct std::complex<double> F2() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "F2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F2();
	}
	struct std::complex<double> F3() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "F3");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F3();
	}
	struct std::complex<double> F4() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "F4");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F4();
	}
	double a1PhaseSpace(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "a1PhaseSpace");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETau2ThreeMesons::a1PhaseSpace(a0);
	}
	struct std::complex<double> a1BreitWigner(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "a1BreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::a1BreitWigner(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETauDecay::calculateME(a0);
	}
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETauDecay::decayWeightMax(a0);
	}
	void calculateResonanceWeights(class std::vector<double, class std::allocator<double> > & a0, class std::vector<double, class std::allocator<double> > & a1, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "calculateResonanceWeights");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::calculateResonanceWeights(a0, a1, a2);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "initPointers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initPointers(a0, a1, a2);
	}
	class Pythia8::HelicityMatrixElement * initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "initChannel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::HelicityMatrixElement *>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::HelicityMatrixElement *> caster;
				return pybind11::detail::cast_ref<class Pythia8::HelicityMatrixElement *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::HelicityMatrixElement *>(std::move(o));
		}
		return HelicityMatrixElement::initChannel(a0);
	}
	double decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "decayWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeight(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "calculateD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateD(a0);
	}
	void calculateRho(unsigned int a0, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "calculateRho");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateRho(a0, a1);
	}
	struct std::complex<double> breitWigner(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "breitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::breitWigner(a0, a1, a2);
	}
	struct std::complex<double> sBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "sBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::sBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> pBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "pBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::pBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> dBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsWithKaons *>(this), "dBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::dBreitWigner(a0, a1, a2, a3, a4);
	}
};

// Pythia8::HMETau2ThreeMesonsGeneric file:Pythia8/HelicityMatrixElements.h line:534
struct PyCallBack_Pythia8_HMETau2ThreeMesonsGeneric : public Pythia8::HMETau2ThreeMesonsGeneric {
	using Pythia8::HMETau2ThreeMesonsGeneric::HMETau2ThreeMesonsGeneric;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initConstants();
	}
	void initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "initHadronicCurrent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initHadronicCurrent(a0);
	}
	void initMode() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "initMode");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initMode();
	}
	void initResonances() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "initResonances");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initResonances();
	}
	void initMomenta(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "initMomenta");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETau2ThreeMesons::initMomenta(a0);
	}
	struct std::complex<double> F1() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "F1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F1();
	}
	struct std::complex<double> F2() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "F2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F2();
	}
	struct std::complex<double> F3() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "F3");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F3();
	}
	struct std::complex<double> F4() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "F4");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::F4();
	}
	double a1PhaseSpace(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "a1PhaseSpace");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETau2ThreeMesons::a1PhaseSpace(a0);
	}
	struct std::complex<double> a1BreitWigner(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "a1BreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETau2ThreeMesons::a1BreitWigner(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETauDecay::calculateME(a0);
	}
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HMETauDecay::decayWeightMax(a0);
	}
	void calculateResonanceWeights(class std::vector<double, class std::allocator<double> > & a0, class std::vector<double, class std::allocator<double> > & a1, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "calculateResonanceWeights");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETauDecay::calculateResonanceWeights(a0, a1, a2);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "initPointers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initPointers(a0, a1, a2);
	}
	class Pythia8::HelicityMatrixElement * initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "initChannel");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::HelicityMatrixElement *>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::HelicityMatrixElement *> caster;
				return pybind11::detail::cast_ref<class Pythia8::HelicityMatrixElement *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::HelicityMatrixElement *>(std::move(o));
		}
		return HelicityMatrixElement::initChannel(a0);
	}
	double decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "decayWeight");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeight(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "calculateD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateD(a0);
	}
	void calculateRho(unsigned int a0, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "calculateRho");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::calculateRho(a0, a1);
	}
	struct std::complex<double> breitWigner(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "breitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::breitWigner(a0, a1, a2);
	}
	struct std::complex<double> sBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "sBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::sBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> pBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "pBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::pBreitWigner(a0, a1, a2, a3, a4);
	}
	struct std::complex<double> dBreitWigner(double a0, double a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETau2ThreeMesonsGeneric *>(this), "dBreitWigner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::dBreitWigner(a0, a1, a2, a3, a4);
	}
};

void bind_Pythia8_HelicityMatrixElements_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::HMETau2Meson file:Pythia8/HelicityMatrixElements.h line:344
		pybind11::class_<Pythia8::HMETau2Meson, std::shared_ptr<Pythia8::HMETau2Meson>, PyCallBack_Pythia8_HMETau2Meson, Pythia8::HMETauDecay> cl(M("Pythia8"), "HMETau2Meson", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETau2Meson(); }, [](){ return new PyCallBack_Pythia8_HMETau2Meson(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETau2Meson const &o){ return new PyCallBack_Pythia8_HMETau2Meson(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETau2Meson const &o){ return new Pythia8::HMETau2Meson(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMETau2Meson::*)()) &Pythia8::HMETau2Meson::initConstants, "C++: Pythia8::HMETau2Meson::initConstants() --> void");
		cl.def("initHadronicCurrent", (void (Pythia8::HMETau2Meson::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETau2Meson::initHadronicCurrent, "C++: Pythia8::HMETau2Meson::initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMETau2Meson & (Pythia8::HMETau2Meson::*)(const class Pythia8::HMETau2Meson &)) &Pythia8::HMETau2Meson::operator=, "C++: Pythia8::HMETau2Meson::operator=(const class Pythia8::HMETau2Meson &) --> class Pythia8::HMETau2Meson &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETau2TwoLeptons file:Pythia8/HelicityMatrixElements.h line:358
		pybind11::class_<Pythia8::HMETau2TwoLeptons, std::shared_ptr<Pythia8::HMETau2TwoLeptons>, PyCallBack_Pythia8_HMETau2TwoLeptons, Pythia8::HMETauDecay> cl(M("Pythia8"), "HMETau2TwoLeptons", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETau2TwoLeptons(); }, [](){ return new PyCallBack_Pythia8_HMETau2TwoLeptons(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETau2TwoLeptons const &o){ return new PyCallBack_Pythia8_HMETau2TwoLeptons(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETau2TwoLeptons const &o){ return new Pythia8::HMETau2TwoLeptons(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMETau2TwoLeptons::*)()) &Pythia8::HMETau2TwoLeptons::initConstants, "C++: Pythia8::HMETau2TwoLeptons::initConstants() --> void");
		cl.def("initWaves", (void (Pythia8::HMETau2TwoLeptons::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETau2TwoLeptons::initWaves, "C++: Pythia8::HMETau2TwoLeptons::initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMETau2TwoLeptons::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMETau2TwoLeptons::calculateME, "C++: Pythia8::HMETau2TwoLeptons::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMETau2TwoLeptons & (Pythia8::HMETau2TwoLeptons::*)(const class Pythia8::HMETau2TwoLeptons &)) &Pythia8::HMETau2TwoLeptons::operator=, "C++: Pythia8::HMETau2TwoLeptons::operator=(const class Pythia8::HMETau2TwoLeptons &) --> class Pythia8::HMETau2TwoLeptons &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETau2TwoMesonsViaVector file:Pythia8/HelicityMatrixElements.h line:375
		pybind11::class_<Pythia8::HMETau2TwoMesonsViaVector, std::shared_ptr<Pythia8::HMETau2TwoMesonsViaVector>, PyCallBack_Pythia8_HMETau2TwoMesonsViaVector, Pythia8::HMETauDecay> cl(M("Pythia8"), "HMETau2TwoMesonsViaVector", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETau2TwoMesonsViaVector(); }, [](){ return new PyCallBack_Pythia8_HMETau2TwoMesonsViaVector(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETau2TwoMesonsViaVector const &o){ return new PyCallBack_Pythia8_HMETau2TwoMesonsViaVector(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETau2TwoMesonsViaVector const &o){ return new Pythia8::HMETau2TwoMesonsViaVector(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMETau2TwoMesonsViaVector::*)()) &Pythia8::HMETau2TwoMesonsViaVector::initConstants, "C++: Pythia8::HMETau2TwoMesonsViaVector::initConstants() --> void");
		cl.def("initHadronicCurrent", (void (Pythia8::HMETau2TwoMesonsViaVector::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETau2TwoMesonsViaVector::initHadronicCurrent, "C++: Pythia8::HMETau2TwoMesonsViaVector::initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMETau2TwoMesonsViaVector & (Pythia8::HMETau2TwoMesonsViaVector::*)(const class Pythia8::HMETau2TwoMesonsViaVector &)) &Pythia8::HMETau2TwoMesonsViaVector::operator=, "C++: Pythia8::HMETau2TwoMesonsViaVector::operator=(const class Pythia8::HMETau2TwoMesonsViaVector &) --> class Pythia8::HMETau2TwoMesonsViaVector &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETau2TwoMesonsViaVectorScalar file:Pythia8/HelicityMatrixElements.h line:396
		pybind11::class_<Pythia8::HMETau2TwoMesonsViaVectorScalar, std::shared_ptr<Pythia8::HMETau2TwoMesonsViaVectorScalar>, PyCallBack_Pythia8_HMETau2TwoMesonsViaVectorScalar, Pythia8::HMETauDecay> cl(M("Pythia8"), "HMETau2TwoMesonsViaVectorScalar", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETau2TwoMesonsViaVectorScalar(); }, [](){ return new PyCallBack_Pythia8_HMETau2TwoMesonsViaVectorScalar(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETau2TwoMesonsViaVectorScalar const &o){ return new PyCallBack_Pythia8_HMETau2TwoMesonsViaVectorScalar(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETau2TwoMesonsViaVectorScalar const &o){ return new Pythia8::HMETau2TwoMesonsViaVectorScalar(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMETau2TwoMesonsViaVectorScalar::*)()) &Pythia8::HMETau2TwoMesonsViaVectorScalar::initConstants, "C++: Pythia8::HMETau2TwoMesonsViaVectorScalar::initConstants() --> void");
		cl.def("initHadronicCurrent", (void (Pythia8::HMETau2TwoMesonsViaVectorScalar::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETau2TwoMesonsViaVectorScalar::initHadronicCurrent, "C++: Pythia8::HMETau2TwoMesonsViaVectorScalar::initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMETau2TwoMesonsViaVectorScalar & (Pythia8::HMETau2TwoMesonsViaVectorScalar::*)(const class Pythia8::HMETau2TwoMesonsViaVectorScalar &)) &Pythia8::HMETau2TwoMesonsViaVectorScalar::operator=, "C++: Pythia8::HMETau2TwoMesonsViaVectorScalar::operator=(const class Pythia8::HMETau2TwoMesonsViaVectorScalar &) --> class Pythia8::HMETau2TwoMesonsViaVectorScalar &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETau2ThreeMesons file:Pythia8/HelicityMatrixElements.h line:421
		pybind11::class_<Pythia8::HMETau2ThreeMesons, std::shared_ptr<Pythia8::HMETau2ThreeMesons>, PyCallBack_Pythia8_HMETau2ThreeMesons, Pythia8::HMETauDecay> cl(M("Pythia8"), "HMETau2ThreeMesons", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETau2ThreeMesons(); }, [](){ return new PyCallBack_Pythia8_HMETau2ThreeMesons(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETau2ThreeMesons const &o){ return new PyCallBack_Pythia8_HMETau2ThreeMesons(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETau2ThreeMesons const &o){ return new Pythia8::HMETau2ThreeMesons(o); } ) );

		pybind11::enum_<Pythia8::HMETau2ThreeMesons::Mode>(cl, "Mode", pybind11::arithmetic(), "")
			.value("Pi0Pi0Pim", Pythia8::HMETau2ThreeMesons::Mode::Pi0Pi0Pim)
			.value("PimPimPip", Pythia8::HMETau2ThreeMesons::Mode::PimPimPip)
			.value("Pi0PimK0b", Pythia8::HMETau2ThreeMesons::Mode::Pi0PimK0b)
			.value("PimPipKm", Pythia8::HMETau2ThreeMesons::Mode::PimPipKm)
			.value("Pi0PimEta", Pythia8::HMETau2ThreeMesons::Mode::Pi0PimEta)
			.value("PimKmKp", Pythia8::HMETau2ThreeMesons::Mode::PimKmKp)
			.value("Pi0K0Km", Pythia8::HMETau2ThreeMesons::Mode::Pi0K0Km)
			.value("KlPimKs", Pythia8::HMETau2ThreeMesons::Mode::KlPimKs)
			.value("Pi0Pi0Km", Pythia8::HMETau2ThreeMesons::Mode::Pi0Pi0Km)
			.value("KlKlPim", Pythia8::HMETau2ThreeMesons::Mode::KlKlPim)
			.value("PimKsKs", Pythia8::HMETau2ThreeMesons::Mode::PimKsKs)
			.value("PimK0bK0", Pythia8::HMETau2ThreeMesons::Mode::PimK0bK0)
			.value("Uknown", Pythia8::HMETau2ThreeMesons::Mode::Uknown)
			.export_values();

		cl.def_readwrite("mode", &Pythia8::HMETau2ThreeMesons::mode);
		cl.def_readwrite("s1", &Pythia8::HMETau2ThreeMesons::s1);
		cl.def_readwrite("s2", &Pythia8::HMETau2ThreeMesons::s2);
		cl.def_readwrite("s3", &Pythia8::HMETau2ThreeMesons::s3);
		cl.def_readwrite("s4", &Pythia8::HMETau2ThreeMesons::s4);
		cl.def_readwrite("q", &Pythia8::HMETau2ThreeMesons::q);
		cl.def_readwrite("q2", &Pythia8::HMETau2ThreeMesons::q2);
		cl.def_readwrite("q3", &Pythia8::HMETau2ThreeMesons::q3);
		cl.def_readwrite("q4", &Pythia8::HMETau2ThreeMesons::q4);
		cl.def_readwrite("a1BW", &Pythia8::HMETau2ThreeMesons::a1BW);
		cl.def("initConstants", (void (Pythia8::HMETau2ThreeMesons::*)()) &Pythia8::HMETau2ThreeMesons::initConstants, "C++: Pythia8::HMETau2ThreeMesons::initConstants() --> void");
		cl.def("initHadronicCurrent", (void (Pythia8::HMETau2ThreeMesons::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETau2ThreeMesons::initHadronicCurrent, "C++: Pythia8::HMETau2ThreeMesons::initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("initMode", (void (Pythia8::HMETau2ThreeMesons::*)()) &Pythia8::HMETau2ThreeMesons::initMode, "C++: Pythia8::HMETau2ThreeMesons::initMode() --> void");
		cl.def("initResonances", (void (Pythia8::HMETau2ThreeMesons::*)()) &Pythia8::HMETau2ThreeMesons::initResonances, "C++: Pythia8::HMETau2ThreeMesons::initResonances() --> void");
		cl.def("initMomenta", (void (Pythia8::HMETau2ThreeMesons::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETau2ThreeMesons::initMomenta, "C++: Pythia8::HMETau2ThreeMesons::initMomenta(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("F1", (struct std::complex<double> (Pythia8::HMETau2ThreeMesons::*)()) &Pythia8::HMETau2ThreeMesons::F1, "C++: Pythia8::HMETau2ThreeMesons::F1() --> struct std::complex<double>");
		cl.def("F2", (struct std::complex<double> (Pythia8::HMETau2ThreeMesons::*)()) &Pythia8::HMETau2ThreeMesons::F2, "C++: Pythia8::HMETau2ThreeMesons::F2() --> struct std::complex<double>");
		cl.def("F3", (struct std::complex<double> (Pythia8::HMETau2ThreeMesons::*)()) &Pythia8::HMETau2ThreeMesons::F3, "C++: Pythia8::HMETau2ThreeMesons::F3() --> struct std::complex<double>");
		cl.def("F4", (struct std::complex<double> (Pythia8::HMETau2ThreeMesons::*)()) &Pythia8::HMETau2ThreeMesons::F4, "C++: Pythia8::HMETau2ThreeMesons::F4() --> struct std::complex<double>");
		cl.def("a1PhaseSpace", (double (Pythia8::HMETau2ThreeMesons::*)(double)) &Pythia8::HMETau2ThreeMesons::a1PhaseSpace, "C++: Pythia8::HMETau2ThreeMesons::a1PhaseSpace(double) --> double", pybind11::arg(""));
		cl.def("a1BreitWigner", (struct std::complex<double> (Pythia8::HMETau2ThreeMesons::*)(double)) &Pythia8::HMETau2ThreeMesons::a1BreitWigner, "C++: Pythia8::HMETau2ThreeMesons::a1BreitWigner(double) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("T", (struct std::complex<double> (Pythia8::HMETau2ThreeMesons::*)(double, double, double, class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &)) &Pythia8::HMETau2ThreeMesons::T, "C++: Pythia8::HMETau2ThreeMesons::T(double, double, double, class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &) --> struct std::complex<double>", pybind11::arg("m0"), pybind11::arg("m1"), pybind11::arg("s"), pybind11::arg("M"), pybind11::arg("G"), pybind11::arg("W"));
		cl.def("T", (struct std::complex<double> (Pythia8::HMETau2ThreeMesons::*)(double, class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &)) &Pythia8::HMETau2ThreeMesons::T, "C++: Pythia8::HMETau2ThreeMesons::T(double, class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &) --> struct std::complex<double>", pybind11::arg("s"), pybind11::arg("M"), pybind11::arg("G"), pybind11::arg("W"));
		cl.def("assign", (class Pythia8::HMETau2ThreeMesons & (Pythia8::HMETau2ThreeMesons::*)(const class Pythia8::HMETau2ThreeMesons &)) &Pythia8::HMETau2ThreeMesons::operator=, "C++: Pythia8::HMETau2ThreeMesons::operator=(const class Pythia8::HMETau2ThreeMesons &) --> class Pythia8::HMETau2ThreeMesons &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETau2ThreePions file:Pythia8/HelicityMatrixElements.h line:473
		pybind11::class_<Pythia8::HMETau2ThreePions, std::shared_ptr<Pythia8::HMETau2ThreePions>, PyCallBack_Pythia8_HMETau2ThreePions, Pythia8::HMETau2ThreeMesons> cl(M("Pythia8"), "HMETau2ThreePions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETau2ThreePions(); }, [](){ return new PyCallBack_Pythia8_HMETau2ThreePions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETau2ThreePions const &o){ return new PyCallBack_Pythia8_HMETau2ThreePions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETau2ThreePions const &o){ return new Pythia8::HMETau2ThreePions(o); } ) );
		cl.def("assign", (class Pythia8::HMETau2ThreePions & (Pythia8::HMETau2ThreePions::*)(const class Pythia8::HMETau2ThreePions &)) &Pythia8::HMETau2ThreePions::operator=, "C++: Pythia8::HMETau2ThreePions::operator=(const class Pythia8::HMETau2ThreePions &) --> class Pythia8::HMETau2ThreePions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETau2ThreeMesonsWithKaons file:Pythia8/HelicityMatrixElements.h line:506
		pybind11::class_<Pythia8::HMETau2ThreeMesonsWithKaons, std::shared_ptr<Pythia8::HMETau2ThreeMesonsWithKaons>, PyCallBack_Pythia8_HMETau2ThreeMesonsWithKaons, Pythia8::HMETau2ThreeMesons> cl(M("Pythia8"), "HMETau2ThreeMesonsWithKaons", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETau2ThreeMesonsWithKaons(); }, [](){ return new PyCallBack_Pythia8_HMETau2ThreeMesonsWithKaons(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETau2ThreeMesonsWithKaons const &o){ return new PyCallBack_Pythia8_HMETau2ThreeMesonsWithKaons(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETau2ThreeMesonsWithKaons const &o){ return new Pythia8::HMETau2ThreeMesonsWithKaons(o); } ) );
		cl.def("assign", (class Pythia8::HMETau2ThreeMesonsWithKaons & (Pythia8::HMETau2ThreeMesonsWithKaons::*)(const class Pythia8::HMETau2ThreeMesonsWithKaons &)) &Pythia8::HMETau2ThreeMesonsWithKaons::operator=, "C++: Pythia8::HMETau2ThreeMesonsWithKaons::operator=(const class Pythia8::HMETau2ThreeMesonsWithKaons &) --> class Pythia8::HMETau2ThreeMesonsWithKaons &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETau2ThreeMesonsGeneric file:Pythia8/HelicityMatrixElements.h line:534
		pybind11::class_<Pythia8::HMETau2ThreeMesonsGeneric, std::shared_ptr<Pythia8::HMETau2ThreeMesonsGeneric>, PyCallBack_Pythia8_HMETau2ThreeMesonsGeneric, Pythia8::HMETau2ThreeMesons> cl(M("Pythia8"), "HMETau2ThreeMesonsGeneric", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETau2ThreeMesonsGeneric(); }, [](){ return new PyCallBack_Pythia8_HMETau2ThreeMesonsGeneric(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETau2ThreeMesonsGeneric const &o){ return new PyCallBack_Pythia8_HMETau2ThreeMesonsGeneric(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETau2ThreeMesonsGeneric const &o){ return new Pythia8::HMETau2ThreeMesonsGeneric(o); } ) );
		cl.def("assign", (class Pythia8::HMETau2ThreeMesonsGeneric & (Pythia8::HMETau2ThreeMesonsGeneric::*)(const class Pythia8::HMETau2ThreeMesonsGeneric &)) &Pythia8::HMETau2ThreeMesonsGeneric::operator=, "C++: Pythia8::HMETau2ThreeMesonsGeneric::operator=(const class Pythia8::HMETau2ThreeMesonsGeneric &) --> class Pythia8::HMETau2ThreeMesonsGeneric &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
