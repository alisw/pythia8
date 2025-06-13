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

// Pythia8::HMETwoFermions2W2TwoFermions file:Pythia8/HelicityMatrixElements.h line:131
struct PyCallBack_Pythia8_HMETwoFermions2W2TwoFermions : public Pythia8::HMETwoFermions2W2TwoFermions {
	using Pythia8::HMETwoFermions2W2TwoFermions::HMETwoFermions2W2TwoFermions;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETwoFermions2W2TwoFermions::initConstants();
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETwoFermions2W2TwoFermions::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETwoFermions2W2TwoFermions::calculateME(a0);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "decayWeight");
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
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeightMax(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2W2TwoFermions *>(this), "dBreitWigner");
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

// Pythia8::HMETwoFermions2GammaZ2TwoFermions file:Pythia8/HelicityMatrixElements.h line:155
struct PyCallBack_Pythia8_HMETwoFermions2GammaZ2TwoFermions : public Pythia8::HMETwoFermions2GammaZ2TwoFermions {
	using Pythia8::HMETwoFermions2GammaZ2TwoFermions::HMETwoFermions2GammaZ2TwoFermions;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETwoFermions2GammaZ2TwoFermions::initConstants();
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETwoFermions2GammaZ2TwoFermions::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETwoFermions2GammaZ2TwoFermions::calculateME(a0);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "decayWeight");
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
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeightMax(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoFermions2GammaZ2TwoFermions *>(this), "dBreitWigner");
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

// Pythia8::HMETwoGammas2TwoFermions file:Pythia8/HelicityMatrixElements.h line:208
struct PyCallBack_Pythia8_HMETwoGammas2TwoFermions : public Pythia8::HMETwoGammas2TwoFermions {
	using Pythia8::HMETwoGammas2TwoFermions::HMETwoGammas2TwoFermions;

	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMETwoGammas2TwoFermions::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMETwoGammas2TwoFermions::calculateME(a0);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "decayWeight");
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
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeightMax(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "dBreitWigner");
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
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETwoGammas2TwoFermions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initConstants();
	}
};

// Pythia8::HMEX2TwoFermions file:Pythia8/HelicityMatrixElements.h line:227
struct PyCallBack_Pythia8_HMEX2TwoFermions : public Pythia8::HMEX2TwoFermions {
	using Pythia8::HMEX2TwoFermions::HMEX2TwoFermions;

	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMEX2TwoFermions::initWaves(a0);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "decayWeight");
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
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeightMax(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HelicityMatrixElement::calculateME(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "dBreitWigner");
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
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEX2TwoFermions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initConstants();
	}
};

// Pythia8::HMEW2TwoFermions file:Pythia8/HelicityMatrixElements.h line:239
struct PyCallBack_Pythia8_HMEW2TwoFermions : public Pythia8::HMEW2TwoFermions {
	using Pythia8::HMEW2TwoFermions::HMEW2TwoFermions;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMEW2TwoFermions::initConstants();
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMEW2TwoFermions::calculateME(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMEX2TwoFermions::initWaves(a0);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "decayWeight");
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
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeightMax(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEW2TwoFermions *>(this), "dBreitWigner");
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

// Pythia8::HMEGamma2TwoFermions file:Pythia8/HelicityMatrixElements.h line:260
struct PyCallBack_Pythia8_HMEGamma2TwoFermions : public Pythia8::HMEGamma2TwoFermions {
	using Pythia8::HMEGamma2TwoFermions::HMEGamma2TwoFermions;

	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMEGamma2TwoFermions::calculateME(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMEX2TwoFermions::initWaves(a0);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "decayWeight");
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
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeightMax(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "dBreitWigner");
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
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEGamma2TwoFermions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initConstants();
	}
};

// Pythia8::HMEZ2TwoFermions file:Pythia8/HelicityMatrixElements.h line:272
struct PyCallBack_Pythia8_HMEZ2TwoFermions : public Pythia8::HMEZ2TwoFermions {
	using Pythia8::HMEZ2TwoFermions::HMEZ2TwoFermions;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMEZ2TwoFermions::initConstants();
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMEZ2TwoFermions::calculateME(a0);
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMEX2TwoFermions::initWaves(a0);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "decayWeight");
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
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeightMax(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEZ2TwoFermions *>(this), "dBreitWigner");
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

// Pythia8::HMEHiggs2TwoFermions file:Pythia8/HelicityMatrixElements.h line:300
struct PyCallBack_Pythia8_HMEHiggs2TwoFermions : public Pythia8::HMEHiggs2TwoFermions {
	using Pythia8::HMEHiggs2TwoFermions::HMEHiggs2TwoFermions;

	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMEHiggs2TwoFermions::initConstants();
	}
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HMEHiggs2TwoFermions::initWaves(a0);
	}
	struct std::complex<double> calculateME(class std::vector<int, class std::allocator<int> > a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "calculateME");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<struct std::complex<double>>::value) {
				static pybind11::detail::override_caster_t<struct std::complex<double>> caster;
				return pybind11::detail::cast_ref<struct std::complex<double>>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct std::complex<double>>(std::move(o));
		}
		return HMEHiggs2TwoFermions::calculateME(a0);
	}
	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "decayWeight");
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
	double decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "decayWeightMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return HelicityMatrixElement::decayWeightMax(a0);
	}
	void calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMEHiggs2TwoFermions *>(this), "dBreitWigner");
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

// Pythia8::HMETauDecay file:Pythia8/HelicityMatrixElements.h line:321
struct PyCallBack_Pythia8_HMETauDecay : public Pythia8::HMETauDecay {
	using Pythia8::HMETauDecay::HMETauDecay;

	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "initWaves");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "calculateME");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "decayWeightMax");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "initHadronicCurrent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "calculateResonanceWeights");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "decayWeight");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "dBreitWigner");
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
	void initConstants() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HMETauDecay *>(this), "initConstants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initConstants();
	}
};

void bind_Pythia8_HelicityMatrixElements(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::HMETwoFermions2W2TwoFermions file:Pythia8/HelicityMatrixElements.h line:131
		pybind11::class_<Pythia8::HMETwoFermions2W2TwoFermions, std::shared_ptr<Pythia8::HMETwoFermions2W2TwoFermions>, PyCallBack_Pythia8_HMETwoFermions2W2TwoFermions, Pythia8::HelicityMatrixElement> cl(M("Pythia8"), "HMETwoFermions2W2TwoFermions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETwoFermions2W2TwoFermions(); }, [](){ return new PyCallBack_Pythia8_HMETwoFermions2W2TwoFermions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETwoFermions2W2TwoFermions const &o){ return new PyCallBack_Pythia8_HMETwoFermions2W2TwoFermions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETwoFermions2W2TwoFermions const &o){ return new Pythia8::HMETwoFermions2W2TwoFermions(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMETwoFermions2W2TwoFermions::*)()) &Pythia8::HMETwoFermions2W2TwoFermions::initConstants, "C++: Pythia8::HMETwoFermions2W2TwoFermions::initConstants() --> void");
		cl.def("initWaves", (void (Pythia8::HMETwoFermions2W2TwoFermions::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETwoFermions2W2TwoFermions::initWaves, "C++: Pythia8::HMETwoFermions2W2TwoFermions::initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMETwoFermions2W2TwoFermions::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMETwoFermions2W2TwoFermions::calculateME, "C++: Pythia8::HMETwoFermions2W2TwoFermions::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMETwoFermions2W2TwoFermions & (Pythia8::HMETwoFermions2W2TwoFermions::*)(const class Pythia8::HMETwoFermions2W2TwoFermions &)) &Pythia8::HMETwoFermions2W2TwoFermions::operator=, "C++: Pythia8::HMETwoFermions2W2TwoFermions::operator=(const class Pythia8::HMETwoFermions2W2TwoFermions &) --> class Pythia8::HMETwoFermions2W2TwoFermions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETwoFermions2GammaZ2TwoFermions file:Pythia8/HelicityMatrixElements.h line:155
		pybind11::class_<Pythia8::HMETwoFermions2GammaZ2TwoFermions, std::shared_ptr<Pythia8::HMETwoFermions2GammaZ2TwoFermions>, PyCallBack_Pythia8_HMETwoFermions2GammaZ2TwoFermions, Pythia8::HelicityMatrixElement> cl(M("Pythia8"), "HMETwoFermions2GammaZ2TwoFermions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETwoFermions2GammaZ2TwoFermions(); }, [](){ return new PyCallBack_Pythia8_HMETwoFermions2GammaZ2TwoFermions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETwoFermions2GammaZ2TwoFermions const &o){ return new PyCallBack_Pythia8_HMETwoFermions2GammaZ2TwoFermions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETwoFermions2GammaZ2TwoFermions const &o){ return new Pythia8::HMETwoFermions2GammaZ2TwoFermions(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMETwoFermions2GammaZ2TwoFermions::*)()) &Pythia8::HMETwoFermions2GammaZ2TwoFermions::initConstants, "C++: Pythia8::HMETwoFermions2GammaZ2TwoFermions::initConstants() --> void");
		cl.def("initWaves", (void (Pythia8::HMETwoFermions2GammaZ2TwoFermions::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETwoFermions2GammaZ2TwoFermions::initWaves, "C++: Pythia8::HMETwoFermions2GammaZ2TwoFermions::initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMETwoFermions2GammaZ2TwoFermions::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMETwoFermions2GammaZ2TwoFermions::calculateME, "C++: Pythia8::HMETwoFermions2GammaZ2TwoFermions::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMETwoFermions2GammaZ2TwoFermions & (Pythia8::HMETwoFermions2GammaZ2TwoFermions::*)(const class Pythia8::HMETwoFermions2GammaZ2TwoFermions &)) &Pythia8::HMETwoFermions2GammaZ2TwoFermions::operator=, "C++: Pythia8::HMETwoFermions2GammaZ2TwoFermions::operator=(const class Pythia8::HMETwoFermions2GammaZ2TwoFermions &) --> class Pythia8::HMETwoFermions2GammaZ2TwoFermions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETwoGammas2TwoFermions file:Pythia8/HelicityMatrixElements.h line:208
		pybind11::class_<Pythia8::HMETwoGammas2TwoFermions, std::shared_ptr<Pythia8::HMETwoGammas2TwoFermions>, PyCallBack_Pythia8_HMETwoGammas2TwoFermions, Pythia8::HelicityMatrixElement> cl(M("Pythia8"), "HMETwoGammas2TwoFermions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMETwoGammas2TwoFermions(); }, [](){ return new PyCallBack_Pythia8_HMETwoGammas2TwoFermions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETwoGammas2TwoFermions const &o){ return new PyCallBack_Pythia8_HMETwoGammas2TwoFermions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETwoGammas2TwoFermions const &o){ return new Pythia8::HMETwoGammas2TwoFermions(o); } ) );
		cl.def("initWaves", (void (Pythia8::HMETwoGammas2TwoFermions::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETwoGammas2TwoFermions::initWaves, "C++: Pythia8::HMETwoGammas2TwoFermions::initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMETwoGammas2TwoFermions::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMETwoGammas2TwoFermions::calculateME, "C++: Pythia8::HMETwoGammas2TwoFermions::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMETwoGammas2TwoFermions & (Pythia8::HMETwoGammas2TwoFermions::*)(const class Pythia8::HMETwoGammas2TwoFermions &)) &Pythia8::HMETwoGammas2TwoFermions::operator=, "C++: Pythia8::HMETwoGammas2TwoFermions::operator=(const class Pythia8::HMETwoGammas2TwoFermions &) --> class Pythia8::HMETwoGammas2TwoFermions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMEX2TwoFermions file:Pythia8/HelicityMatrixElements.h line:227
		pybind11::class_<Pythia8::HMEX2TwoFermions, std::shared_ptr<Pythia8::HMEX2TwoFermions>, PyCallBack_Pythia8_HMEX2TwoFermions, Pythia8::HelicityMatrixElement> cl(M("Pythia8"), "HMEX2TwoFermions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMEX2TwoFermions(); }, [](){ return new PyCallBack_Pythia8_HMEX2TwoFermions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMEX2TwoFermions const &o){ return new PyCallBack_Pythia8_HMEX2TwoFermions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMEX2TwoFermions const &o){ return new Pythia8::HMEX2TwoFermions(o); } ) );
		cl.def("initWaves", (void (Pythia8::HMEX2TwoFermions::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMEX2TwoFermions::initWaves, "C++: Pythia8::HMEX2TwoFermions::initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMEX2TwoFermions & (Pythia8::HMEX2TwoFermions::*)(const class Pythia8::HMEX2TwoFermions &)) &Pythia8::HMEX2TwoFermions::operator=, "C++: Pythia8::HMEX2TwoFermions::operator=(const class Pythia8::HMEX2TwoFermions &) --> class Pythia8::HMEX2TwoFermions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMEW2TwoFermions file:Pythia8/HelicityMatrixElements.h line:239
		pybind11::class_<Pythia8::HMEW2TwoFermions, std::shared_ptr<Pythia8::HMEW2TwoFermions>, PyCallBack_Pythia8_HMEW2TwoFermions, Pythia8::HMEX2TwoFermions> cl(M("Pythia8"), "HMEW2TwoFermions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMEW2TwoFermions(); }, [](){ return new PyCallBack_Pythia8_HMEW2TwoFermions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMEW2TwoFermions const &o){ return new PyCallBack_Pythia8_HMEW2TwoFermions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMEW2TwoFermions const &o){ return new Pythia8::HMEW2TwoFermions(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMEW2TwoFermions::*)()) &Pythia8::HMEW2TwoFermions::initConstants, "C++: Pythia8::HMEW2TwoFermions::initConstants() --> void");
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMEW2TwoFermions::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMEW2TwoFermions::calculateME, "C++: Pythia8::HMEW2TwoFermions::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMEW2TwoFermions & (Pythia8::HMEW2TwoFermions::*)(const class Pythia8::HMEW2TwoFermions &)) &Pythia8::HMEW2TwoFermions::operator=, "C++: Pythia8::HMEW2TwoFermions::operator=(const class Pythia8::HMEW2TwoFermions &) --> class Pythia8::HMEW2TwoFermions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMEGamma2TwoFermions file:Pythia8/HelicityMatrixElements.h line:260
		pybind11::class_<Pythia8::HMEGamma2TwoFermions, std::shared_ptr<Pythia8::HMEGamma2TwoFermions>, PyCallBack_Pythia8_HMEGamma2TwoFermions, Pythia8::HMEX2TwoFermions> cl(M("Pythia8"), "HMEGamma2TwoFermions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMEGamma2TwoFermions(); }, [](){ return new PyCallBack_Pythia8_HMEGamma2TwoFermions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMEGamma2TwoFermions const &o){ return new PyCallBack_Pythia8_HMEGamma2TwoFermions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMEGamma2TwoFermions const &o){ return new Pythia8::HMEGamma2TwoFermions(o); } ) );
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMEGamma2TwoFermions::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMEGamma2TwoFermions::calculateME, "C++: Pythia8::HMEGamma2TwoFermions::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMEGamma2TwoFermions & (Pythia8::HMEGamma2TwoFermions::*)(const class Pythia8::HMEGamma2TwoFermions &)) &Pythia8::HMEGamma2TwoFermions::operator=, "C++: Pythia8::HMEGamma2TwoFermions::operator=(const class Pythia8::HMEGamma2TwoFermions &) --> class Pythia8::HMEGamma2TwoFermions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMEZ2TwoFermions file:Pythia8/HelicityMatrixElements.h line:272
		pybind11::class_<Pythia8::HMEZ2TwoFermions, std::shared_ptr<Pythia8::HMEZ2TwoFermions>, PyCallBack_Pythia8_HMEZ2TwoFermions, Pythia8::HMEX2TwoFermions> cl(M("Pythia8"), "HMEZ2TwoFermions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMEZ2TwoFermions(); }, [](){ return new PyCallBack_Pythia8_HMEZ2TwoFermions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMEZ2TwoFermions const &o){ return new PyCallBack_Pythia8_HMEZ2TwoFermions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMEZ2TwoFermions const &o){ return new Pythia8::HMEZ2TwoFermions(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMEZ2TwoFermions::*)()) &Pythia8::HMEZ2TwoFermions::initConstants, "C++: Pythia8::HMEZ2TwoFermions::initConstants() --> void");
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMEZ2TwoFermions::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMEZ2TwoFermions::calculateME, "C++: Pythia8::HMEZ2TwoFermions::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMEZ2TwoFermions & (Pythia8::HMEZ2TwoFermions::*)(const class Pythia8::HMEZ2TwoFermions &)) &Pythia8::HMEZ2TwoFermions::operator=, "C++: Pythia8::HMEZ2TwoFermions::operator=(const class Pythia8::HMEZ2TwoFermions &) --> class Pythia8::HMEZ2TwoFermions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMEHiggs2TwoFermions file:Pythia8/HelicityMatrixElements.h line:300
		pybind11::class_<Pythia8::HMEHiggs2TwoFermions, std::shared_ptr<Pythia8::HMEHiggs2TwoFermions>, PyCallBack_Pythia8_HMEHiggs2TwoFermions, Pythia8::HelicityMatrixElement> cl(M("Pythia8"), "HMEHiggs2TwoFermions", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HMEHiggs2TwoFermions(); }, [](){ return new PyCallBack_Pythia8_HMEHiggs2TwoFermions(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMEHiggs2TwoFermions const &o){ return new PyCallBack_Pythia8_HMEHiggs2TwoFermions(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMEHiggs2TwoFermions const &o){ return new Pythia8::HMEHiggs2TwoFermions(o); } ) );
		cl.def("initConstants", (void (Pythia8::HMEHiggs2TwoFermions::*)()) &Pythia8::HMEHiggs2TwoFermions::initConstants, "C++: Pythia8::HMEHiggs2TwoFermions::initConstants() --> void");
		cl.def("initWaves", (void (Pythia8::HMEHiggs2TwoFermions::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMEHiggs2TwoFermions::initWaves, "C++: Pythia8::HMEHiggs2TwoFermions::initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMEHiggs2TwoFermions::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMEHiggs2TwoFermions::calculateME, "C++: Pythia8::HMEHiggs2TwoFermions::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMEHiggs2TwoFermions & (Pythia8::HMEHiggs2TwoFermions::*)(const class Pythia8::HMEHiggs2TwoFermions &)) &Pythia8::HMEHiggs2TwoFermions::operator=, "C++: Pythia8::HMEHiggs2TwoFermions::operator=(const class Pythia8::HMEHiggs2TwoFermions &) --> class Pythia8::HMEHiggs2TwoFermions &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HMETauDecay file:Pythia8/HelicityMatrixElements.h line:321
		pybind11::class_<Pythia8::HMETauDecay, std::shared_ptr<Pythia8::HMETauDecay>, PyCallBack_Pythia8_HMETauDecay, Pythia8::HelicityMatrixElement> cl(M("Pythia8"), "HMETauDecay", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](PyCallBack_Pythia8_HMETauDecay const &o){ return new PyCallBack_Pythia8_HMETauDecay(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HMETauDecay const &o){ return new Pythia8::HMETauDecay(o); } ) );
		cl.def( pybind11::init( [](){ return new Pythia8::HMETauDecay(); }, [](){ return new PyCallBack_Pythia8_HMETauDecay(); } ) );
		cl.def("initWaves", (void (Pythia8::HMETauDecay::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETauDecay::initWaves, "C++: Pythia8::HMETauDecay::initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HMETauDecay::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HMETauDecay::calculateME, "C++: Pythia8::HMETauDecay::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("decayWeightMax", (double (Pythia8::HMETauDecay::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETauDecay::decayWeightMax, "C++: Pythia8::HMETauDecay::decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> double", pybind11::arg(""));
		cl.def("initHadronicCurrent", (void (Pythia8::HMETauDecay::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HMETauDecay::initHadronicCurrent, "C++: Pythia8::HMETauDecay::initHadronicCurrent(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("calculateResonanceWeights", (void (Pythia8::HMETauDecay::*)(class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > &)) &Pythia8::HMETauDecay::calculateResonanceWeights, "C++: Pythia8::HMETauDecay::calculateResonanceWeights(class std::vector<double, class std::allocator<double> > &, class std::vector<double, class std::allocator<double> > &, class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > &) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("assign", (class Pythia8::HMETauDecay & (Pythia8::HMETauDecay::*)(const class Pythia8::HMETauDecay &)) &Pythia8::HMETauDecay::operator=, "C++: Pythia8::HMETauDecay::operator=(const class Pythia8::HMETauDecay &) --> class Pythia8::HMETauDecay &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
