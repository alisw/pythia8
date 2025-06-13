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
#include <Pythia8/NucleonExcitations.h>
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

// Pythia8::SigmaSaSDL file:Pythia8/SigmaTotal.h line:294
struct PyCallBack_Pythia8_SigmaSaSDL : public Pythia8::SigmaSaSDL {
	using Pythia8::SigmaSaSDL::SigmaSaSDL;

	void init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SigmaSaSDL::init(a0);
	}
	bool calcTotEl(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "calcTotEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaSaSDL::calcTotEl(a0, a1, a2, a3, a4);
	}
	double dsigmaEl(double a0, bool a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "dsigmaEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaSaSDL::dsigmaEl(a0, a1, a2);
	}
	bool calcDiff(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "calcDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaSaSDL::calcDiff(a0, a1, a2, a3, a4);
	}
	double dsigmaSD(double a0, double a1, bool a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "dsigmaSD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaSaSDL::dsigmaSD(a0, a1, a2, a3);
	}
	double dsigmaDD(double a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "dsigmaDD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaSaSDL::dsigmaDD(a0, a1, a2, a3);
	}
	double dsigmaCD(double a0, double a1, double a2, double a3, int a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "dsigmaCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaSaSDL::dsigmaCD(a0, a1, a2, a3, a4);
	}
	double mMinCD() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "mMinCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaSaSDL::mMinCD();
	}
	bool calcTot(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "calcTot");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "splitDiff");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "initCoulomb");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "addCoulomb");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaSaSDL *>(this), "dsigmaElCoulomb");
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

// Pythia8::SigmaMBR file:Pythia8/SigmaTotal.h line:368
struct PyCallBack_Pythia8_SigmaMBR : public Pythia8::SigmaMBR {
	using Pythia8::SigmaMBR::SigmaMBR;

	void init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SigmaMBR::init(a0);
	}
	bool calcTotEl(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "calcTotEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaMBR::calcTotEl(a0, a1, a2, a3, a4);
	}
	double dsigmaEl(double a0, bool a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "dsigmaEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaMBR::dsigmaEl(a0, a1, a2);
	}
	bool calcDiff(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "calcDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaMBR::calcDiff(a0, a1, a2, a3, a4);
	}
	bool splitDiff() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "splitDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaMBR::splitDiff();
	}
	double dsigmaSD(double a0, double a1, bool a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "dsigmaSD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaMBR::dsigmaSD(a0, a1, a2, a3);
	}
	double dsigmaDD(double a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "dsigmaDD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaMBR::dsigmaDD(a0, a1, a2, a3);
	}
	double dsigmaCD(double a0, double a1, double a2, double a3, int a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "dsigmaCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaMBR::dsigmaCD(a0, a1, a2, a3, a4);
	}
	double mMinCD() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "mMinCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaMBR::mMinCD();
	}
	bool calcTot(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "calcTot");
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
	bool initCoulomb(class Pythia8::Settings & a0, class Pythia8::ParticleData * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "initCoulomb");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "addCoulomb");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaMBR *>(this), "dsigmaElCoulomb");
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

// Pythia8::SigmaABMST file:Pythia8/SigmaTotal.h line:432
struct PyCallBack_Pythia8_SigmaABMST : public Pythia8::SigmaABMST {
	using Pythia8::SigmaABMST::SigmaABMST;

	void init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SigmaABMST::init(a0);
	}
	bool calcTotEl(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "calcTotEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaABMST::calcTotEl(a0, a1, a2, a3, a4);
	}
	double dsigmaEl(double a0, bool a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "dsigmaEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaABMST::dsigmaEl(a0, a1, a2);
	}
	bool calcDiff(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "calcDiff");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaABMST::calcDiff(a0, a1, a2, a3, a4);
	}
	double dsigmaSD(double a0, double a1, bool a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "dsigmaSD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaABMST::dsigmaSD(a0, a1, a2, a3);
	}
	double dsigmaDD(double a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "dsigmaDD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaABMST::dsigmaDD(a0, a1, a2, a3);
	}
	double dsigmaCD(double a0, double a1, double a2, double a3, int a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "dsigmaCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaABMST::dsigmaCD(a0, a1, a2, a3, a4);
	}
	double mMinCD() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "mMinCD");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaABMST::mMinCD();
	}
	bool calcTot(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "calcTot");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "splitDiff");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "initCoulomb");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "addCoulomb");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaABMST *>(this), "dsigmaElCoulomb");
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

// Pythia8::SigmaRPP file:Pythia8/SigmaTotal.h line:529
struct PyCallBack_Pythia8_SigmaRPP : public Pythia8::SigmaRPP {
	using Pythia8::SigmaRPP::SigmaRPP;

	void init(class Pythia8::Info * a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SigmaRPP::init(a0);
	}
	bool calcTotEl(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "calcTotEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SigmaRPP::calcTotEl(a0, a1, a2, a3, a4);
	}
	double dsigmaEl(double a0, bool a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "dsigmaEl");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SigmaRPP::dsigmaEl(a0, a1, a2);
	}
	bool calcTot(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "calcTot");
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
	bool calcDiff(int a0, int a1, double a2, double a3, double a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "calcDiff");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "dsigmaSD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "splitDiff");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "dsigmaDD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "dsigmaCD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "mMinCD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "initCoulomb");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "addCoulomb");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaRPP *>(this), "dsigmaElCoulomb");
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

// Pythia8::SigmaLowEnergy file:Pythia8/SigmaLowEnergy.h line:22
struct PyCallBack_Pythia8_SigmaLowEnergy : public Pythia8::SigmaLowEnergy {
	using Pythia8::SigmaLowEnergy::SigmaLowEnergy;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaLowEnergy *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaLowEnergy *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaLowEnergy *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaLowEnergy *>(this), "onStat");
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

void bind_Pythia8_SigmaTotal_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::SigmaSaSDL file:Pythia8/SigmaTotal.h line:294
		pybind11::class_<Pythia8::SigmaSaSDL, std::shared_ptr<Pythia8::SigmaSaSDL>, PyCallBack_Pythia8_SigmaSaSDL, Pythia8::SigmaTotAux> cl(M("Pythia8"), "SigmaSaSDL", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SigmaSaSDL(); }, [](){ return new PyCallBack_Pythia8_SigmaSaSDL(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_SigmaSaSDL const &o){ return new PyCallBack_Pythia8_SigmaSaSDL(o); } ) );
		cl.def( pybind11::init( [](Pythia8::SigmaSaSDL const &o){ return new Pythia8::SigmaSaSDL(o); } ) );
		cl.def("init", (void (Pythia8::SigmaSaSDL::*)(class Pythia8::Info *)) &Pythia8::SigmaSaSDL::init, "C++: Pythia8::SigmaSaSDL::init(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("sigmaTotal", (double (Pythia8::SigmaSaSDL::*)(int, int, double, double, double)) &Pythia8::SigmaSaSDL::sigmaTotal, "C++: Pythia8::SigmaSaSDL::sigmaTotal(int, int, double, double, double) --> double", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg("sIn"), pybind11::arg("mAin"), pybind11::arg("mBin"));
		cl.def("calcTotEl", (bool (Pythia8::SigmaSaSDL::*)(int, int, double, double, double)) &Pythia8::SigmaSaSDL::calcTotEl, "C++: Pythia8::SigmaSaSDL::calcTotEl(int, int, double, double, double) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg("sIn"), pybind11::arg("mAin"), pybind11::arg("mBin"));
		cl.def("dsigmaEl", [](Pythia8::SigmaSaSDL &o, double const & a0) -> double { return o.dsigmaEl(a0); }, "", pybind11::arg("t"));
		cl.def("dsigmaEl", [](Pythia8::SigmaSaSDL &o, double const & a0, bool const & a1) -> double { return o.dsigmaEl(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("useCoulomb"));
		cl.def("dsigmaEl", (double (Pythia8::SigmaSaSDL::*)(double, bool, bool)) &Pythia8::SigmaSaSDL::dsigmaEl, "C++: Pythia8::SigmaSaSDL::dsigmaEl(double, bool, bool) --> double", pybind11::arg("t"), pybind11::arg("useCoulomb"), pybind11::arg(""));
		cl.def("calcDiff", (bool (Pythia8::SigmaSaSDL::*)(int, int, double, double, double)) &Pythia8::SigmaSaSDL::calcDiff, "C++: Pythia8::SigmaSaSDL::calcDiff(int, int, double, double, double) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg("sIn"), pybind11::arg("mAin"), pybind11::arg("mBin"));
		cl.def("dsigmaSD", [](Pythia8::SigmaSaSDL &o, double const & a0, double const & a1, bool const & a2) -> double { return o.dsigmaSD(a0, a1, a2); }, "", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg("isXB"));
		cl.def("dsigmaSD", (double (Pythia8::SigmaSaSDL::*)(double, double, bool, int)) &Pythia8::SigmaSaSDL::dsigmaSD, "C++: Pythia8::SigmaSaSDL::dsigmaSD(double, double, bool, int) --> double", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg("isXB"), pybind11::arg(""));
		cl.def("dsigmaDD", [](Pythia8::SigmaSaSDL &o, double const & a0, double const & a1, double const & a2) -> double { return o.dsigmaDD(a0, a1, a2); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"));
		cl.def("dsigmaDD", (double (Pythia8::SigmaSaSDL::*)(double, double, double, int)) &Pythia8::SigmaSaSDL::dsigmaDD, "C++: Pythia8::SigmaSaSDL::dsigmaDD(double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"), pybind11::arg(""));
		cl.def("dsigmaCD", [](Pythia8::SigmaSaSDL &o, double const & a0, double const & a1, double const & a2, double const & a3) -> double { return o.dsigmaCD(a0, a1, a2, a3); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"));
		cl.def("dsigmaCD", (double (Pythia8::SigmaSaSDL::*)(double, double, double, double, int)) &Pythia8::SigmaSaSDL::dsigmaCD, "C++: Pythia8::SigmaSaSDL::dsigmaCD(double, double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"), pybind11::arg(""));
		cl.def("mMinCD", (double (Pythia8::SigmaSaSDL::*)()) &Pythia8::SigmaSaSDL::mMinCD, "C++: Pythia8::SigmaSaSDL::mMinCD() --> double");
		cl.def("assign", (class Pythia8::SigmaSaSDL & (Pythia8::SigmaSaSDL::*)(const class Pythia8::SigmaSaSDL &)) &Pythia8::SigmaSaSDL::operator=, "C++: Pythia8::SigmaSaSDL::operator=(const class Pythia8::SigmaSaSDL &) --> class Pythia8::SigmaSaSDL &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SigmaMBR file:Pythia8/SigmaTotal.h line:368
		pybind11::class_<Pythia8::SigmaMBR, std::shared_ptr<Pythia8::SigmaMBR>, PyCallBack_Pythia8_SigmaMBR, Pythia8::SigmaTotAux> cl(M("Pythia8"), "SigmaMBR", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SigmaMBR(); }, [](){ return new PyCallBack_Pythia8_SigmaMBR(); } ) );
		cl.def("init", (void (Pythia8::SigmaMBR::*)(class Pythia8::Info *)) &Pythia8::SigmaMBR::init, "C++: Pythia8::SigmaMBR::init(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("calcTotEl", (bool (Pythia8::SigmaMBR::*)(int, int, double, double, double)) &Pythia8::SigmaMBR::calcTotEl, "C++: Pythia8::SigmaMBR::calcTotEl(int, int, double, double, double) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg("sIn"), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaEl", [](Pythia8::SigmaMBR &o, double const & a0) -> double { return o.dsigmaEl(a0); }, "", pybind11::arg("t"));
		cl.def("dsigmaEl", [](Pythia8::SigmaMBR &o, double const & a0, bool const & a1) -> double { return o.dsigmaEl(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("useCoulomb"));
		cl.def("dsigmaEl", (double (Pythia8::SigmaMBR::*)(double, bool, bool)) &Pythia8::SigmaMBR::dsigmaEl, "C++: Pythia8::SigmaMBR::dsigmaEl(double, bool, bool) --> double", pybind11::arg("t"), pybind11::arg("useCoulomb"), pybind11::arg(""));
		cl.def("calcDiff", (bool (Pythia8::SigmaMBR::*)(int, int, double, double, double)) &Pythia8::SigmaMBR::calcDiff, "C++: Pythia8::SigmaMBR::calcDiff(int, int, double, double, double) --> bool", pybind11::arg(""), pybind11::arg(""), pybind11::arg("sIn"), pybind11::arg(""), pybind11::arg(""));
		cl.def("splitDiff", (bool (Pythia8::SigmaMBR::*)()) &Pythia8::SigmaMBR::splitDiff, "C++: Pythia8::SigmaMBR::splitDiff() --> bool");
		cl.def("dsigmaSD", [](Pythia8::SigmaMBR &o, double const & a0, double const & a1) -> double { return o.dsigmaSD(a0, a1); }, "", pybind11::arg("xi"), pybind11::arg("t"));
		cl.def("dsigmaSD", [](Pythia8::SigmaMBR &o, double const & a0, double const & a1, bool const & a2) -> double { return o.dsigmaSD(a0, a1, a2); }, "", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg(""));
		cl.def("dsigmaSD", (double (Pythia8::SigmaMBR::*)(double, double, bool, int)) &Pythia8::SigmaMBR::dsigmaSD, "C++: Pythia8::SigmaMBR::dsigmaSD(double, double, bool, int) --> double", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg(""), pybind11::arg("step"));
		cl.def("dsigmaDD", [](Pythia8::SigmaMBR &o, double const & a0, double const & a1, double const & a2) -> double { return o.dsigmaDD(a0, a1, a2); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"));
		cl.def("dsigmaDD", (double (Pythia8::SigmaMBR::*)(double, double, double, int)) &Pythia8::SigmaMBR::dsigmaDD, "C++: Pythia8::SigmaMBR::dsigmaDD(double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"), pybind11::arg("step"));
		cl.def("dsigmaCD", [](Pythia8::SigmaMBR &o, double const & a0, double const & a1, double const & a2, double const & a3) -> double { return o.dsigmaCD(a0, a1, a2, a3); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"));
		cl.def("dsigmaCD", (double (Pythia8::SigmaMBR::*)(double, double, double, double, int)) &Pythia8::SigmaMBR::dsigmaCD, "C++: Pythia8::SigmaMBR::dsigmaCD(double, double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"), pybind11::arg("step"));
		cl.def("mMinCD", (double (Pythia8::SigmaMBR::*)()) &Pythia8::SigmaMBR::mMinCD, "C++: Pythia8::SigmaMBR::mMinCD() --> double");
		cl.def("assign", (class Pythia8::SigmaMBR & (Pythia8::SigmaMBR::*)(const class Pythia8::SigmaMBR &)) &Pythia8::SigmaMBR::operator=, "C++: Pythia8::SigmaMBR::operator=(const class Pythia8::SigmaMBR &) --> class Pythia8::SigmaMBR &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SigmaABMST file:Pythia8/SigmaTotal.h line:432
		pybind11::class_<Pythia8::SigmaABMST, std::shared_ptr<Pythia8::SigmaABMST>, PyCallBack_Pythia8_SigmaABMST, Pythia8::SigmaTotAux> cl(M("Pythia8"), "SigmaABMST", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SigmaABMST(); }, [](){ return new PyCallBack_Pythia8_SigmaABMST(); } ) );
		cl.def("init", (void (Pythia8::SigmaABMST::*)(class Pythia8::Info *)) &Pythia8::SigmaABMST::init, "C++: Pythia8::SigmaABMST::init(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("calcTotEl", (bool (Pythia8::SigmaABMST::*)(int, int, double, double, double)) &Pythia8::SigmaABMST::calcTotEl, "C++: Pythia8::SigmaABMST::calcTotEl(int, int, double, double, double) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg("sIn"), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaEl", [](Pythia8::SigmaABMST &o, double const & a0) -> double { return o.dsigmaEl(a0); }, "", pybind11::arg("t"));
		cl.def("dsigmaEl", [](Pythia8::SigmaABMST &o, double const & a0, bool const & a1) -> double { return o.dsigmaEl(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("useCoulomb"));
		cl.def("dsigmaEl", (double (Pythia8::SigmaABMST::*)(double, bool, bool)) &Pythia8::SigmaABMST::dsigmaEl, "C++: Pythia8::SigmaABMST::dsigmaEl(double, bool, bool) --> double", pybind11::arg("t"), pybind11::arg("useCoulomb"), pybind11::arg("onlyPomerons"));
		cl.def("calcDiff", (bool (Pythia8::SigmaABMST::*)(int, int, double, double, double)) &Pythia8::SigmaABMST::calcDiff, "C++: Pythia8::SigmaABMST::calcDiff(int, int, double, double, double) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg("sIn"), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaSD", [](Pythia8::SigmaABMST &o, double const & a0, double const & a1) -> double { return o.dsigmaSD(a0, a1); }, "", pybind11::arg("xi"), pybind11::arg("t"));
		cl.def("dsigmaSD", [](Pythia8::SigmaABMST &o, double const & a0, double const & a1, bool const & a2) -> double { return o.dsigmaSD(a0, a1, a2); }, "", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg(""));
		cl.def("dsigmaSD", (double (Pythia8::SigmaABMST::*)(double, double, bool, int)) &Pythia8::SigmaABMST::dsigmaSD, "C++: Pythia8::SigmaABMST::dsigmaSD(double, double, bool, int) --> double", pybind11::arg("xi"), pybind11::arg("t"), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaDD", [](Pythia8::SigmaABMST &o, double const & a0, double const & a1, double const & a2) -> double { return o.dsigmaDD(a0, a1, a2); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"));
		cl.def("dsigmaDD", (double (Pythia8::SigmaABMST::*)(double, double, double, int)) &Pythia8::SigmaABMST::dsigmaDD, "C++: Pythia8::SigmaABMST::dsigmaDD(double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t"), pybind11::arg(""));
		cl.def("dsigmaCD", [](Pythia8::SigmaABMST &o, double const & a0, double const & a1, double const & a2, double const & a3) -> double { return o.dsigmaCD(a0, a1, a2, a3); }, "", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"));
		cl.def("dsigmaCD", (double (Pythia8::SigmaABMST::*)(double, double, double, double, int)) &Pythia8::SigmaABMST::dsigmaCD, "C++: Pythia8::SigmaABMST::dsigmaCD(double, double, double, double, int) --> double", pybind11::arg("xi1"), pybind11::arg("xi2"), pybind11::arg("t1"), pybind11::arg("t2"), pybind11::arg(""));
		cl.def("mMinCD", (double (Pythia8::SigmaABMST::*)()) &Pythia8::SigmaABMST::mMinCD, "C++: Pythia8::SigmaABMST::mMinCD() --> double");
		cl.def("assign", (class Pythia8::SigmaABMST & (Pythia8::SigmaABMST::*)(const class Pythia8::SigmaABMST &)) &Pythia8::SigmaABMST::operator=, "C++: Pythia8::SigmaABMST::operator=(const class Pythia8::SigmaABMST &) --> class Pythia8::SigmaABMST &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SigmaRPP file:Pythia8/SigmaTotal.h line:529
		pybind11::class_<Pythia8::SigmaRPP, std::shared_ptr<Pythia8::SigmaRPP>, PyCallBack_Pythia8_SigmaRPP, Pythia8::SigmaTotAux> cl(M("Pythia8"), "SigmaRPP", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SigmaRPP(); }, [](){ return new PyCallBack_Pythia8_SigmaRPP(); } ) );
		cl.def("init", (void (Pythia8::SigmaRPP::*)(class Pythia8::Info *)) &Pythia8::SigmaRPP::init, "C++: Pythia8::SigmaRPP::init(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("calcTotEl", (bool (Pythia8::SigmaRPP::*)(int, int, double, double, double)) &Pythia8::SigmaRPP::calcTotEl, "C++: Pythia8::SigmaRPP::calcTotEl(int, int, double, double, double) --> bool", pybind11::arg("idAin"), pybind11::arg("idBin"), pybind11::arg("sIn"), pybind11::arg(""), pybind11::arg(""));
		cl.def("dsigmaEl", [](Pythia8::SigmaRPP &o, double const & a0) -> double { return o.dsigmaEl(a0); }, "", pybind11::arg("t"));
		cl.def("dsigmaEl", [](Pythia8::SigmaRPP &o, double const & a0, bool const & a1) -> double { return o.dsigmaEl(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("useCoulomb"));
		cl.def("dsigmaEl", (double (Pythia8::SigmaRPP::*)(double, bool, bool)) &Pythia8::SigmaRPP::dsigmaEl, "C++: Pythia8::SigmaRPP::dsigmaEl(double, bool, bool) --> double", pybind11::arg("t"), pybind11::arg("useCoulomb"), pybind11::arg(""));
		cl.def("assign", (class Pythia8::SigmaRPP & (Pythia8::SigmaRPP::*)(const class Pythia8::SigmaRPP &)) &Pythia8::SigmaRPP::operator=, "C++: Pythia8::SigmaRPP::operator=(const class Pythia8::SigmaRPP &) --> class Pythia8::SigmaRPP &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SigmaLowEnergy file:Pythia8/SigmaLowEnergy.h line:22
		pybind11::class_<Pythia8::SigmaLowEnergy, std::shared_ptr<Pythia8::SigmaLowEnergy>, PyCallBack_Pythia8_SigmaLowEnergy, Pythia8::PhysicsBase> cl(M("Pythia8"), "SigmaLowEnergy", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SigmaLowEnergy(); }, [](){ return new PyCallBack_Pythia8_SigmaLowEnergy(); } ) );
		cl.def("init", (void (Pythia8::SigmaLowEnergy::*)(class Pythia8::NucleonExcitations *)) &Pythia8::SigmaLowEnergy::init, "C++: Pythia8::SigmaLowEnergy::init(class Pythia8::NucleonExcitations *) --> void", pybind11::arg("nucleonExcitationsPtrIn"));
		cl.def("sigmaTotal", (double (Pythia8::SigmaLowEnergy::*)(int, int, double, double, double)) &Pythia8::SigmaLowEnergy::sigmaTotal, "C++: Pythia8::SigmaLowEnergy::sigmaTotal(int, int, double, double, double) --> double", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"), pybind11::arg("mA"), pybind11::arg("mB"));
		cl.def("sigmaTotal", (double (Pythia8::SigmaLowEnergy::*)(int, int, double)) &Pythia8::SigmaLowEnergy::sigmaTotal, "C++: Pythia8::SigmaLowEnergy::sigmaTotal(int, int, double) --> double", pybind11::arg("idAIn"), pybind11::arg("idBIn"), pybind11::arg("eCMIn"));
		cl.def("sigmaPartial", (double (Pythia8::SigmaLowEnergy::*)(int, int, double, double, double, int)) &Pythia8::SigmaLowEnergy::sigmaPartial, "C++: Pythia8::SigmaLowEnergy::sigmaPartial(int, int, double, double, double, int) --> double", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"), pybind11::arg("mA"), pybind11::arg("mB"), pybind11::arg("proc"));
		cl.def("sigmaPartial", (double (Pythia8::SigmaLowEnergy::*)(int, int, double, int)) &Pythia8::SigmaLowEnergy::sigmaPartial, "C++: Pythia8::SigmaLowEnergy::sigmaPartial(int, int, double, int) --> double", pybind11::arg("idAIn"), pybind11::arg("idBIn"), pybind11::arg("eCMIn"), pybind11::arg("proc"));
		cl.def("sigmaPartial", (bool (Pythia8::SigmaLowEnergy::*)(int, int, double, double, double, class std::vector<int, class std::allocator<int> > &, class std::vector<double, class std::allocator<double> > &)) &Pythia8::SigmaLowEnergy::sigmaPartial, "C++: Pythia8::SigmaLowEnergy::sigmaPartial(int, int, double, double, double, class std::vector<int, class std::allocator<int> > &, class std::vector<double, class std::allocator<double> > &) --> bool", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"), pybind11::arg("mA"), pybind11::arg("mB"), pybind11::arg("procsOut"), pybind11::arg("sigmasOut"));
		cl.def("pickProcess", (int (Pythia8::SigmaLowEnergy::*)(int, int, double, double, double)) &Pythia8::SigmaLowEnergy::pickProcess, "C++: Pythia8::SigmaLowEnergy::pickProcess(int, int, double, double, double) --> int", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"), pybind11::arg("mA"), pybind11::arg("mB"));
		cl.def("pickResonance", (int (Pythia8::SigmaLowEnergy::*)(int, int, double)) &Pythia8::SigmaLowEnergy::pickResonance, "C++: Pythia8::SigmaLowEnergy::pickResonance(int, int, double) --> int", pybind11::arg("idA"), pybind11::arg("idB"), pybind11::arg("eCM"));
		cl.def("hasExcitation", (bool (Pythia8::SigmaLowEnergy::*)(int, int) const) &Pythia8::SigmaLowEnergy::hasExcitation, "C++: Pythia8::SigmaLowEnergy::hasExcitation(int, int) const --> bool", pybind11::arg("idAIn"), pybind11::arg("idBIn"));
		cl.def("updateResonances", (void (Pythia8::SigmaLowEnergy::*)()) &Pythia8::SigmaLowEnergy::updateResonances, "C++: Pythia8::SigmaLowEnergy::updateResonances() --> void");
		cl.def("assign", (class Pythia8::SigmaLowEnergy & (Pythia8::SigmaLowEnergy::*)(const class Pythia8::SigmaLowEnergy &)) &Pythia8::SigmaLowEnergy::operator=, "C++: Pythia8::SigmaLowEnergy::operator=(const class Pythia8::SigmaLowEnergy &) --> class Pythia8::SigmaLowEnergy &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
