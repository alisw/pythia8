#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonDistributions.h>
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
#include <ostream>
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

// Pythia8::Nucleus2gamma file:Pythia8/PartonDistributions.h line:957
struct PyCallBack_Pythia8_Nucleus2gamma : public Pythia8::Nucleus2gamma {
	using Pythia8::Nucleus2gamma::Nucleus2gamma;

	void setBeamID(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "setBeamID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setBeamID(a0);
	}
	void setExtrapolate(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "setExtrapolate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setExtrapolate(a0);
	}
	bool insideBounds(double a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "insideBounds");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::insideBounds(a0, a1);
	}
	double alphaS(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "alphaS");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::alphaS(a0);
	}
	double mQuarkPDF(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "mQuarkPDF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::mQuarkPDF(a0);
	}
	int nMembers() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "nMembers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::nMembers();
	}
	void calcPDFEnvelope(int a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	void calcPDFEnvelope(struct std::pair<int, int> a0, struct std::pair<double, double> a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	struct Pythia8::PDF::PDFEnvelope getPDFEnvelope() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "getPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::PDF::PDFEnvelope>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::PDF::PDFEnvelope> caster;
				return pybind11::detail::cast_ref<struct Pythia8::PDF::PDFEnvelope>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::PDF::PDFEnvelope>(std::move(o));
		}
		return PDF::getPDFEnvelope();
	}
	double gammaPDFxDependence(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "gammaPDFxDependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFxDependence(a0, a1);
	}
	double gammaPDFRefScale(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "gammaPDFRefScale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFRefScale(a0);
	}
	int sampleGammaValFlavor(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "sampleGammaValFlavor");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::sampleGammaValFlavor(a0);
	}
	double xfIntegratedTotal(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xfIntegratedTotal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfIntegratedTotal(a0);
	}
	double xGamma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xGamma();
	}
	void xPom(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xPom");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::xPom(a0);
	}
	double xfFlux(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xfFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfFlux(a0, a1, a2);
	}
	double xfApprox(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xfApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfApprox(a0, a1, a2);
	}
	double xfGamma(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xfGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfGamma(a0, a1, a2);
	}
	double intFluxApprox() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "intFluxApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::intFluxApprox();
	}
	bool hasApproxGammaFlux() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "hasApproxGammaFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::hasApproxGammaFlux();
	}
	double getXmin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "getXmin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXmin();
	}
	double getXhadr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "getXhadr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXhadr();
	}
	double sampleXgamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "sampleXgamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleXgamma(a0);
	}
	double sampleQ2gamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "sampleQ2gamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleQ2gamma(a0);
	}
	double fluxQ2dependence(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "fluxQ2dependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::fluxQ2dependence(a0);
	}
	double xfMax(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xfMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfMax(a0, a1, a2);
	}
	double xfSame(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xfSame");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfSame(a0, a1, a2);
	}
	void setVMDscale(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "setVMDscale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setVMDscale(a0);
	}
	void xfUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Nucleus2gamma *>(this), "xfUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"PDF::xfUpdate\"");
	}
};

// Pythia8::EPAexternal file:Pythia8/PartonDistributions.h line:989
struct PyCallBack_Pythia8_EPAexternal : public Pythia8::EPAexternal {
	using Pythia8::EPAexternal::EPAexternal;

	void xfUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xfUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return EPAexternal::xfUpdate(a0, a1, a2);
	}
	double xfFlux(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xfFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return EPAexternal::xfFlux(a0, a1, a2);
	}
	double xfGamma(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xfGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return EPAexternal::xfGamma(a0, a1, a2);
	}
	double xfApprox(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xfApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return EPAexternal::xfApprox(a0, a1, a2);
	}
	double intFluxApprox() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "intFluxApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return EPAexternal::intFluxApprox();
	}
	bool hasApproxGammaFlux() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "hasApproxGammaFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return EPAexternal::hasApproxGammaFlux();
	}
	double getXmin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "getXmin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return EPAexternal::getXmin();
	}
	double getXhadr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "getXhadr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return EPAexternal::getXhadr();
	}
	double sampleXgamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "sampleXgamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return EPAexternal::sampleXgamma(a0);
	}
	double sampleQ2gamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "sampleQ2gamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return EPAexternal::sampleQ2gamma(a0);
	}
	void setBeamID(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "setBeamID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setBeamID(a0);
	}
	void setExtrapolate(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "setExtrapolate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setExtrapolate(a0);
	}
	bool insideBounds(double a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "insideBounds");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::insideBounds(a0, a1);
	}
	double alphaS(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "alphaS");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::alphaS(a0);
	}
	double mQuarkPDF(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "mQuarkPDF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::mQuarkPDF(a0);
	}
	int nMembers() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "nMembers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::nMembers();
	}
	void calcPDFEnvelope(int a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	void calcPDFEnvelope(struct std::pair<int, int> a0, struct std::pair<double, double> a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	struct Pythia8::PDF::PDFEnvelope getPDFEnvelope() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "getPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::PDF::PDFEnvelope>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::PDF::PDFEnvelope> caster;
				return pybind11::detail::cast_ref<struct Pythia8::PDF::PDFEnvelope>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::PDF::PDFEnvelope>(std::move(o));
		}
		return PDF::getPDFEnvelope();
	}
	double gammaPDFxDependence(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "gammaPDFxDependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFxDependence(a0, a1);
	}
	double gammaPDFRefScale(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "gammaPDFRefScale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFRefScale(a0);
	}
	int sampleGammaValFlavor(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "sampleGammaValFlavor");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::sampleGammaValFlavor(a0);
	}
	double xfIntegratedTotal(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xfIntegratedTotal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfIntegratedTotal(a0);
	}
	double xGamma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xGamma();
	}
	void xPom(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xPom");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::xPom(a0);
	}
	double fluxQ2dependence(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "fluxQ2dependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::fluxQ2dependence(a0);
	}
	double xfMax(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xfMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfMax(a0, a1, a2);
	}
	double xfSame(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "xfSame");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfSame(a0, a1, a2);
	}
	void setVMDscale(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPAexternal *>(this), "setVMDscale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setVMDscale(a0);
	}
};

// Pythia8::nPDF file:Pythia8/PartonDistributions.h line:1060
struct PyCallBack_Pythia8_nPDF : public Pythia8::nPDF {
	using Pythia8::nPDF::nPDF;

	void xfUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xfUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return nPDF::xfUpdate(a0, a1, a2);
	}
	void rUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "rUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"nPDF::rUpdate\"");
	}
	void setBeamID(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "setBeamID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setBeamID(a0);
	}
	void setExtrapolate(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "setExtrapolate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setExtrapolate(a0);
	}
	bool insideBounds(double a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "insideBounds");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::insideBounds(a0, a1);
	}
	double alphaS(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "alphaS");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::alphaS(a0);
	}
	double mQuarkPDF(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "mQuarkPDF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::mQuarkPDF(a0);
	}
	int nMembers() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "nMembers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::nMembers();
	}
	void calcPDFEnvelope(int a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	void calcPDFEnvelope(struct std::pair<int, int> a0, struct std::pair<double, double> a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	struct Pythia8::PDF::PDFEnvelope getPDFEnvelope() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "getPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::PDF::PDFEnvelope>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::PDF::PDFEnvelope> caster;
				return pybind11::detail::cast_ref<struct Pythia8::PDF::PDFEnvelope>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::PDF::PDFEnvelope>(std::move(o));
		}
		return PDF::getPDFEnvelope();
	}
	double gammaPDFxDependence(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "gammaPDFxDependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFxDependence(a0, a1);
	}
	double gammaPDFRefScale(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "gammaPDFRefScale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFRefScale(a0);
	}
	int sampleGammaValFlavor(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "sampleGammaValFlavor");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::sampleGammaValFlavor(a0);
	}
	double xfIntegratedTotal(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xfIntegratedTotal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfIntegratedTotal(a0);
	}
	double xGamma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xGamma();
	}
	void xPom(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xPom");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::xPom(a0);
	}
	double xfFlux(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xfFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfFlux(a0, a1, a2);
	}
	double xfApprox(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xfApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfApprox(a0, a1, a2);
	}
	double xfGamma(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xfGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfGamma(a0, a1, a2);
	}
	double intFluxApprox() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "intFluxApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::intFluxApprox();
	}
	bool hasApproxGammaFlux() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "hasApproxGammaFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::hasApproxGammaFlux();
	}
	double getXmin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "getXmin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXmin();
	}
	double getXhadr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "getXhadr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXhadr();
	}
	double sampleXgamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "sampleXgamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleXgamma(a0);
	}
	double sampleQ2gamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "sampleQ2gamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleQ2gamma(a0);
	}
	double fluxQ2dependence(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "fluxQ2dependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::fluxQ2dependence(a0);
	}
	double xfMax(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xfMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfMax(a0, a1, a2);
	}
	double xfSame(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "xfSame");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfSame(a0, a1, a2);
	}
	void setVMDscale(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::nPDF *>(this), "setVMDscale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setVMDscale(a0);
	}
};

// Pythia8::Isospin file:Pythia8/PartonDistributions.h line:1110
struct PyCallBack_Pythia8_Isospin : public Pythia8::Isospin {
	using Pythia8::Isospin::Isospin;

	void rUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "rUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Isospin::rUpdate(a0, a1, a2);
	}
	void xfUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xfUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return nPDF::xfUpdate(a0, a1, a2);
	}
	void setBeamID(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "setBeamID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setBeamID(a0);
	}
	void setExtrapolate(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "setExtrapolate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setExtrapolate(a0);
	}
	bool insideBounds(double a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "insideBounds");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::insideBounds(a0, a1);
	}
	double alphaS(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "alphaS");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::alphaS(a0);
	}
	double mQuarkPDF(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "mQuarkPDF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::mQuarkPDF(a0);
	}
	int nMembers() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "nMembers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::nMembers();
	}
	void calcPDFEnvelope(int a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	void calcPDFEnvelope(struct std::pair<int, int> a0, struct std::pair<double, double> a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	struct Pythia8::PDF::PDFEnvelope getPDFEnvelope() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "getPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::PDF::PDFEnvelope>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::PDF::PDFEnvelope> caster;
				return pybind11::detail::cast_ref<struct Pythia8::PDF::PDFEnvelope>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::PDF::PDFEnvelope>(std::move(o));
		}
		return PDF::getPDFEnvelope();
	}
	double gammaPDFxDependence(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "gammaPDFxDependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFxDependence(a0, a1);
	}
	double gammaPDFRefScale(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "gammaPDFRefScale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFRefScale(a0);
	}
	int sampleGammaValFlavor(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "sampleGammaValFlavor");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::sampleGammaValFlavor(a0);
	}
	double xfIntegratedTotal(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xfIntegratedTotal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfIntegratedTotal(a0);
	}
	double xGamma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xGamma();
	}
	void xPom(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xPom");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::xPom(a0);
	}
	double xfFlux(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xfFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfFlux(a0, a1, a2);
	}
	double xfApprox(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xfApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfApprox(a0, a1, a2);
	}
	double xfGamma(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xfGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfGamma(a0, a1, a2);
	}
	double intFluxApprox() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "intFluxApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::intFluxApprox();
	}
	bool hasApproxGammaFlux() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "hasApproxGammaFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::hasApproxGammaFlux();
	}
	double getXmin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "getXmin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXmin();
	}
	double getXhadr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "getXhadr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXhadr();
	}
	double sampleXgamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "sampleXgamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleXgamma(a0);
	}
	double sampleQ2gamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "sampleQ2gamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleQ2gamma(a0);
	}
	double fluxQ2dependence(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "fluxQ2dependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::fluxQ2dependence(a0);
	}
	double xfMax(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xfMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfMax(a0, a1, a2);
	}
	double xfSame(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "xfSame");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfSame(a0, a1, a2);
	}
	void setVMDscale(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Isospin *>(this), "setVMDscale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setVMDscale(a0);
	}
};

// Pythia8::EPS09 file:Pythia8/PartonDistributions.h line:1126
struct PyCallBack_Pythia8_EPS09 : public Pythia8::EPS09 {
	using Pythia8::EPS09::EPS09;

	void rUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "rUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return EPS09::rUpdate(a0, a1, a2);
	}
	void xfUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xfUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return nPDF::xfUpdate(a0, a1, a2);
	}
	void setBeamID(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "setBeamID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setBeamID(a0);
	}
	void setExtrapolate(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "setExtrapolate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setExtrapolate(a0);
	}
	bool insideBounds(double a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "insideBounds");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::insideBounds(a0, a1);
	}
	double alphaS(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "alphaS");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::alphaS(a0);
	}
	double mQuarkPDF(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "mQuarkPDF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::mQuarkPDF(a0);
	}
	int nMembers() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "nMembers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::nMembers();
	}
	void calcPDFEnvelope(int a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	void calcPDFEnvelope(struct std::pair<int, int> a0, struct std::pair<double, double> a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	struct Pythia8::PDF::PDFEnvelope getPDFEnvelope() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "getPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::PDF::PDFEnvelope>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::PDF::PDFEnvelope> caster;
				return pybind11::detail::cast_ref<struct Pythia8::PDF::PDFEnvelope>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::PDF::PDFEnvelope>(std::move(o));
		}
		return PDF::getPDFEnvelope();
	}
	double gammaPDFxDependence(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "gammaPDFxDependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFxDependence(a0, a1);
	}
	double gammaPDFRefScale(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "gammaPDFRefScale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFRefScale(a0);
	}
	int sampleGammaValFlavor(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "sampleGammaValFlavor");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::sampleGammaValFlavor(a0);
	}
	double xfIntegratedTotal(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xfIntegratedTotal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfIntegratedTotal(a0);
	}
	double xGamma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xGamma();
	}
	void xPom(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xPom");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::xPom(a0);
	}
	double xfFlux(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xfFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfFlux(a0, a1, a2);
	}
	double xfApprox(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xfApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfApprox(a0, a1, a2);
	}
	double xfGamma(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xfGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfGamma(a0, a1, a2);
	}
	double intFluxApprox() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "intFluxApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::intFluxApprox();
	}
	bool hasApproxGammaFlux() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "hasApproxGammaFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::hasApproxGammaFlux();
	}
	double getXmin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "getXmin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXmin();
	}
	double getXhadr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "getXhadr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXhadr();
	}
	double sampleXgamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "sampleXgamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleXgamma(a0);
	}
	double sampleQ2gamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "sampleQ2gamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleQ2gamma(a0);
	}
	double fluxQ2dependence(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "fluxQ2dependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::fluxQ2dependence(a0);
	}
	double xfMax(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xfMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfMax(a0, a1, a2);
	}
	double xfSame(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "xfSame");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfSame(a0, a1, a2);
	}
	void setVMDscale(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPS09 *>(this), "setVMDscale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setVMDscale(a0);
	}
};

// Pythia8::EPPS16 file:Pythia8/PartonDistributions.h line:1170
struct PyCallBack_Pythia8_EPPS16 : public Pythia8::EPPS16 {
	using Pythia8::EPPS16::EPPS16;

	void rUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "rUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return EPPS16::rUpdate(a0, a1, a2);
	}
	void xfUpdate(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xfUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return nPDF::xfUpdate(a0, a1, a2);
	}
	void setBeamID(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "setBeamID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setBeamID(a0);
	}
	void setExtrapolate(bool a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "setExtrapolate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setExtrapolate(a0);
	}
	bool insideBounds(double a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "insideBounds");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::insideBounds(a0, a1);
	}
	double alphaS(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "alphaS");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::alphaS(a0);
	}
	double mQuarkPDF(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "mQuarkPDF");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::mQuarkPDF(a0);
	}
	int nMembers() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "nMembers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::nMembers();
	}
	void calcPDFEnvelope(int a0, double a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	void calcPDFEnvelope(struct std::pair<int, int> a0, struct std::pair<double, double> a1, double a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "calcPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::calcPDFEnvelope(a0, a1, a2, a3);
	}
	struct Pythia8::PDF::PDFEnvelope getPDFEnvelope() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "getPDFEnvelope");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::PDF::PDFEnvelope>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::PDF::PDFEnvelope> caster;
				return pybind11::detail::cast_ref<struct Pythia8::PDF::PDFEnvelope>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::PDF::PDFEnvelope>(std::move(o));
		}
		return PDF::getPDFEnvelope();
	}
	double gammaPDFxDependence(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "gammaPDFxDependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFxDependence(a0, a1);
	}
	double gammaPDFRefScale(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "gammaPDFRefScale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::gammaPDFRefScale(a0);
	}
	int sampleGammaValFlavor(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "sampleGammaValFlavor");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PDF::sampleGammaValFlavor(a0);
	}
	double xfIntegratedTotal(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xfIntegratedTotal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfIntegratedTotal(a0);
	}
	double xGamma() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xGamma();
	}
	void xPom(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xPom");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::xPom(a0);
	}
	double xfFlux(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xfFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfFlux(a0, a1, a2);
	}
	double xfApprox(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xfApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfApprox(a0, a1, a2);
	}
	double xfGamma(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xfGamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfGamma(a0, a1, a2);
	}
	double intFluxApprox() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "intFluxApprox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::intFluxApprox();
	}
	bool hasApproxGammaFlux() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "hasApproxGammaFlux");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return PDF::hasApproxGammaFlux();
	}
	double getXmin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "getXmin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXmin();
	}
	double getXhadr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "getXhadr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::getXhadr();
	}
	double sampleXgamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "sampleXgamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleXgamma(a0);
	}
	double sampleQ2gamma(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "sampleQ2gamma");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::sampleQ2gamma(a0);
	}
	double fluxQ2dependence(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "fluxQ2dependence");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::fluxQ2dependence(a0);
	}
	double xfMax(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xfMax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfMax(a0, a1, a2);
	}
	double xfSame(int a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "xfSame");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PDF::xfSame(a0, a1, a2);
	}
	void setVMDscale(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::EPPS16 *>(this), "setVMDscale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PDF::setVMDscale(a0);
	}
};

void bind_Pythia8_PartonDistributions_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::Nucleus2gamma file:Pythia8/PartonDistributions.h line:957
		pybind11::class_<Pythia8::Nucleus2gamma, std::shared_ptr<Pythia8::Nucleus2gamma>, PyCallBack_Pythia8_Nucleus2gamma, Pythia8::PDF> cl(M("Pythia8"), "Nucleus2gamma", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<int, double, double>(), pybind11::arg("idBeamIn"), pybind11::arg("bMinIn"), pybind11::arg("mNucleonIn") );

		cl.def("assign", (class Pythia8::Nucleus2gamma & (Pythia8::Nucleus2gamma::*)(const class Pythia8::Nucleus2gamma &)) &Pythia8::Nucleus2gamma::operator=, "C++: Pythia8::Nucleus2gamma::operator=(const class Pythia8::Nucleus2gamma &) --> class Pythia8::Nucleus2gamma &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::EPAexternal file:Pythia8/PartonDistributions.h line:989
		pybind11::class_<Pythia8::EPAexternal, std::shared_ptr<Pythia8::EPAexternal>, PyCallBack_Pythia8_EPAexternal, Pythia8::PDF> cl(M("Pythia8"), "EPAexternal", "");
		pybind11::handle cl_type = cl;

		cl.def("xfUpdate", (void (Pythia8::EPAexternal::*)(int, double, double)) &Pythia8::EPAexternal::xfUpdate, "C++: Pythia8::EPAexternal::xfUpdate(int, double, double) --> void", pybind11::arg(""), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfFlux", [](Pythia8::EPAexternal &o, int const & a0, double const & a1) -> double { return o.xfFlux(a0, a1); }, "", pybind11::arg("id"), pybind11::arg("x"));
		cl.def("xfFlux", (double (Pythia8::EPAexternal::*)(int, double, double)) &Pythia8::EPAexternal::xfFlux, "C++: Pythia8::EPAexternal::xfFlux(int, double, double) --> double", pybind11::arg("id"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfGamma", (double (Pythia8::EPAexternal::*)(int, double, double)) &Pythia8::EPAexternal::xfGamma, "C++: Pythia8::EPAexternal::xfGamma(int, double, double) --> double", pybind11::arg("id"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("xfApprox", (double (Pythia8::EPAexternal::*)(int, double, double)) &Pythia8::EPAexternal::xfApprox, "C++: Pythia8::EPAexternal::xfApprox(int, double, double) --> double", pybind11::arg("id"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("intFluxApprox", (double (Pythia8::EPAexternal::*)()) &Pythia8::EPAexternal::intFluxApprox, "C++: Pythia8::EPAexternal::intFluxApprox() --> double");
		cl.def("hasApproxGammaFlux", (bool (Pythia8::EPAexternal::*)()) &Pythia8::EPAexternal::hasApproxGammaFlux, "C++: Pythia8::EPAexternal::hasApproxGammaFlux() --> bool");
		cl.def("getXmin", (double (Pythia8::EPAexternal::*)()) &Pythia8::EPAexternal::getXmin, "C++: Pythia8::EPAexternal::getXmin() --> double");
		cl.def("getXhadr", (double (Pythia8::EPAexternal::*)()) &Pythia8::EPAexternal::getXhadr, "C++: Pythia8::EPAexternal::getXhadr() --> double");
		cl.def("sampleXgamma", (double (Pythia8::EPAexternal::*)(double)) &Pythia8::EPAexternal::sampleXgamma, "C++: Pythia8::EPAexternal::sampleXgamma(double) --> double", pybind11::arg("xMinIn"));
		cl.def("sampleQ2gamma", (double (Pythia8::EPAexternal::*)(double)) &Pythia8::EPAexternal::sampleQ2gamma, "C++: Pythia8::EPAexternal::sampleQ2gamma(double) --> double", pybind11::arg("Q2minIn"));
		cl.def("assign", (class Pythia8::EPAexternal & (Pythia8::EPAexternal::*)(const class Pythia8::EPAexternal &)) &Pythia8::EPAexternal::operator=, "C++: Pythia8::EPAexternal::operator=(const class Pythia8::EPAexternal &) --> class Pythia8::EPAexternal &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::nPDF file:Pythia8/PartonDistributions.h line:1060
		pybind11::class_<Pythia8::nPDF, std::shared_ptr<Pythia8::nPDF>, PyCallBack_Pythia8_nPDF, Pythia8::PDF> cl(M("Pythia8"), "nPDF", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_nPDF(); } ), "doc");
		cl.def( pybind11::init( [](int const & a0){ return new PyCallBack_Pythia8_nPDF(a0); } ), "doc");
		cl.def( pybind11::init<int, class std::shared_ptr<class Pythia8::PDF>>(), pybind11::arg("idBeamIn"), pybind11::arg("protonPDFPtrIn") );

		cl.def(pybind11::init<PyCallBack_Pythia8_nPDF const &>());
		cl.def_readwrite("ruv", &Pythia8::nPDF::ruv);
		cl.def_readwrite("rdv", &Pythia8::nPDF::rdv);
		cl.def_readwrite("ru", &Pythia8::nPDF::ru);
		cl.def_readwrite("rd", &Pythia8::nPDF::rd);
		cl.def_readwrite("rs", &Pythia8::nPDF::rs);
		cl.def_readwrite("rc", &Pythia8::nPDF::rc);
		cl.def_readwrite("rb", &Pythia8::nPDF::rb);
		cl.def_readwrite("rg", &Pythia8::nPDF::rg);
		cl.def("xfUpdate", (void (Pythia8::nPDF::*)(int, double, double)) &Pythia8::nPDF::xfUpdate, "C++: Pythia8::nPDF::xfUpdate(int, double, double) --> void", pybind11::arg("id"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("rUpdate", (void (Pythia8::nPDF::*)(int, double, double)) &Pythia8::nPDF::rUpdate, "C++: Pythia8::nPDF::rUpdate(int, double, double) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("initNPDF", [](Pythia8::nPDF &o, int const & a0) -> void { return o.initNPDF(a0); }, "", pybind11::arg("idBeamIn"));
		cl.def("initNPDF", (void (Pythia8::nPDF::*)(int, class std::shared_ptr<class Pythia8::PDF>)) &Pythia8::nPDF::initNPDF, "C++: Pythia8::nPDF::initNPDF(int, class std::shared_ptr<class Pythia8::PDF>) --> void", pybind11::arg("idBeamIn"), pybind11::arg("protonPDFPtrIn"));
		cl.def("getA", (int (Pythia8::nPDF::*)()) &Pythia8::nPDF::getA, "C++: Pythia8::nPDF::getA() --> int");
		cl.def("getZ", (int (Pythia8::nPDF::*)()) &Pythia8::nPDF::getZ, "C++: Pythia8::nPDF::getZ() --> int");
		cl.def("setMode", (void (Pythia8::nPDF::*)(double)) &Pythia8::nPDF::setMode, "C++: Pythia8::nPDF::setMode(double) --> void", pybind11::arg("zaIn"));
		cl.def("resetMode", (void (Pythia8::nPDF::*)()) &Pythia8::nPDF::resetMode, "C++: Pythia8::nPDF::resetMode() --> void");
		cl.def("assign", (class Pythia8::nPDF & (Pythia8::nPDF::*)(const class Pythia8::nPDF &)) &Pythia8::nPDF::operator=, "C++: Pythia8::nPDF::operator=(const class Pythia8::nPDF &) --> class Pythia8::nPDF &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::Isospin file:Pythia8/PartonDistributions.h line:1110
		pybind11::class_<Pythia8::Isospin, std::shared_ptr<Pythia8::Isospin>, PyCallBack_Pythia8_Isospin, Pythia8::nPDF> cl(M("Pythia8"), "Isospin", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Isospin(); }, [](){ return new PyCallBack_Pythia8_Isospin(); } ), "doc");
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::Isospin(a0); }, [](int const & a0){ return new PyCallBack_Pythia8_Isospin(a0); } ), "doc");
		cl.def( pybind11::init<int, class std::shared_ptr<class Pythia8::PDF>>(), pybind11::arg("idBeamIn"), pybind11::arg("protonPDFPtrIn") );

		cl.def("rUpdate", (void (Pythia8::Isospin::*)(int, double, double)) &Pythia8::Isospin::rUpdate, "C++: Pythia8::Isospin::rUpdate(int, double, double) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("assign", (class Pythia8::Isospin & (Pythia8::Isospin::*)(const class Pythia8::Isospin &)) &Pythia8::Isospin::operator=, "C++: Pythia8::Isospin::operator=(const class Pythia8::Isospin &) --> class Pythia8::Isospin &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::EPS09 file:Pythia8/PartonDistributions.h line:1126
		pybind11::class_<Pythia8::EPS09, std::shared_ptr<Pythia8::EPS09>, PyCallBack_Pythia8_EPS09, Pythia8::nPDF> cl(M("Pythia8"), "EPS09", "");
		pybind11::handle cl_type = cl;

		cl.def("rUpdate", (void (Pythia8::EPS09::*)(int, double, double)) &Pythia8::EPS09::rUpdate, "C++: Pythia8::EPS09::rUpdate(int, double, double) --> void", pybind11::arg("id"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("setErrorSet", (void (Pythia8::EPS09::*)(int)) &Pythia8::EPS09::setErrorSet, "C++: Pythia8::EPS09::setErrorSet(int) --> void", pybind11::arg("iSetIn"));
		cl.def("assign", (class Pythia8::EPS09 & (Pythia8::EPS09::*)(const class Pythia8::EPS09 &)) &Pythia8::EPS09::operator=, "C++: Pythia8::EPS09::operator=(const class Pythia8::EPS09 &) --> class Pythia8::EPS09 &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::EPPS16 file:Pythia8/PartonDistributions.h line:1170
		pybind11::class_<Pythia8::EPPS16, std::shared_ptr<Pythia8::EPPS16>, PyCallBack_Pythia8_EPPS16, Pythia8::nPDF> cl(M("Pythia8"), "EPPS16", "");
		pybind11::handle cl_type = cl;

		cl.def("rUpdate", (void (Pythia8::EPPS16::*)(int, double, double)) &Pythia8::EPPS16::rUpdate, "C++: Pythia8::EPPS16::rUpdate(int, double, double) --> void", pybind11::arg("id"), pybind11::arg("x"), pybind11::arg("Q2"));
		cl.def("setErrorSet", (void (Pythia8::EPPS16::*)(int)) &Pythia8::EPPS16::setErrorSet, "C++: Pythia8::EPPS16::setErrorSet(int) --> void", pybind11::arg("iSetIn"));
		cl.def("assign", (class Pythia8::EPPS16 & (Pythia8::EPPS16::*)(const class Pythia8::EPPS16 &)) &Pythia8::EPPS16::operator=, "C++: Pythia8::EPPS16::operator=(const class Pythia8::EPPS16 &) --> class Pythia8::EPPS16 &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ResolvedParton file:Pythia8/BeamParticle.h line:38
		pybind11::class_<Pythia8::ResolvedParton, std::shared_ptr<Pythia8::ResolvedParton>> cl(M("Pythia8"), "ResolvedParton", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ResolvedParton(); } ), "doc" );
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::ResolvedParton(a0); } ), "doc" , pybind11::arg("iPosIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Pythia8::ResolvedParton(a0, a1); } ), "doc" , pybind11::arg("iPosIn"), pybind11::arg("idIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2){ return new Pythia8::ResolvedParton(a0, a1, a2); } ), "doc" , pybind11::arg("iPosIn"), pybind11::arg("idIn"), pybind11::arg("xIn"));
		cl.def( pybind11::init<int, int, double, int>(), pybind11::arg("iPosIn"), pybind11::arg("idIn"), pybind11::arg("xIn"), pybind11::arg("companionIn") );

		cl.def( pybind11::init( [](Pythia8::ResolvedParton const &o){ return new Pythia8::ResolvedParton(o); } ) );
		cl.def("iPos", (void (Pythia8::ResolvedParton::*)(int)) &Pythia8::ResolvedParton::iPos, "C++: Pythia8::ResolvedParton::iPos(int) --> void", pybind11::arg("iPosIn"));
		cl.def("id", (void (Pythia8::ResolvedParton::*)(int)) &Pythia8::ResolvedParton::id, "C++: Pythia8::ResolvedParton::id(int) --> void", pybind11::arg("idIn"));
		cl.def("x", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::x, "C++: Pythia8::ResolvedParton::x(double) --> void", pybind11::arg("xIn"));
		cl.def("update", (void (Pythia8::ResolvedParton::*)(int, int, double)) &Pythia8::ResolvedParton::update, "C++: Pythia8::ResolvedParton::update(int, int, double) --> void", pybind11::arg("iPosIn"), pybind11::arg("idIn"), pybind11::arg("xIn"));
		cl.def("companion", (void (Pythia8::ResolvedParton::*)(int)) &Pythia8::ResolvedParton::companion, "C++: Pythia8::ResolvedParton::companion(int) --> void", pybind11::arg("companionIn"));
		cl.def("xqCompanion", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::xqCompanion, "C++: Pythia8::ResolvedParton::xqCompanion(double) --> void", pybind11::arg("xqCompIn"));
		cl.def("p", (void (Pythia8::ResolvedParton::*)(class Pythia8::Vec4)) &Pythia8::ResolvedParton::p, "C++: Pythia8::ResolvedParton::p(class Pythia8::Vec4) --> void", pybind11::arg("pIn"));
		cl.def("px", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::px, "C++: Pythia8::ResolvedParton::px(double) --> void", pybind11::arg("pxIn"));
		cl.def("py", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::py, "C++: Pythia8::ResolvedParton::py(double) --> void", pybind11::arg("pyIn"));
		cl.def("pz", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::pz, "C++: Pythia8::ResolvedParton::pz(double) --> void", pybind11::arg("pzIn"));
		cl.def("e", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::e, "C++: Pythia8::ResolvedParton::e(double) --> void", pybind11::arg("eIn"));
		cl.def("m", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::m, "C++: Pythia8::ResolvedParton::m(double) --> void", pybind11::arg("mIn"));
		cl.def("col", (void (Pythia8::ResolvedParton::*)(int)) &Pythia8::ResolvedParton::col, "C++: Pythia8::ResolvedParton::col(int) --> void", pybind11::arg("colIn"));
		cl.def("acol", (void (Pythia8::ResolvedParton::*)(int)) &Pythia8::ResolvedParton::acol, "C++: Pythia8::ResolvedParton::acol(int) --> void", pybind11::arg("acolIn"));
		cl.def("cols", [](Pythia8::ResolvedParton &o) -> void { return o.cols(); }, "");
		cl.def("cols", [](Pythia8::ResolvedParton &o, int const & a0) -> void { return o.cols(a0); }, "", pybind11::arg("colIn"));
		cl.def("cols", (void (Pythia8::ResolvedParton::*)(int, int)) &Pythia8::ResolvedParton::cols, "C++: Pythia8::ResolvedParton::cols(int, int) --> void", pybind11::arg("colIn"), pybind11::arg("acolIn"));
		cl.def("scalePT", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::scalePT, "C++: Pythia8::ResolvedParton::scalePT(double) --> void", pybind11::arg("factorIn"));
		cl.def("scaleX", (void (Pythia8::ResolvedParton::*)(double)) &Pythia8::ResolvedParton::scaleX, "C++: Pythia8::ResolvedParton::scaleX(double) --> void", pybind11::arg("factorIn"));
		cl.def("iPos", (int (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::iPos, "C++: Pythia8::ResolvedParton::iPos() const --> int");
		cl.def("id", (int (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::id, "C++: Pythia8::ResolvedParton::id() const --> int");
		cl.def("x", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::x, "C++: Pythia8::ResolvedParton::x() const --> double");
		cl.def("companion", (int (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::companion, "C++: Pythia8::ResolvedParton::companion() const --> int");
		cl.def("isValence", (bool (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::isValence, "C++: Pythia8::ResolvedParton::isValence() const --> bool");
		cl.def("isUnmatched", (bool (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::isUnmatched, "C++: Pythia8::ResolvedParton::isUnmatched() const --> bool");
		cl.def("isCompanion", (bool (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::isCompanion, "C++: Pythia8::ResolvedParton::isCompanion() const --> bool");
		cl.def("isFromBeam", (bool (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::isFromBeam, "C++: Pythia8::ResolvedParton::isFromBeam() const --> bool");
		cl.def("xqCompanion", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::xqCompanion, "C++: Pythia8::ResolvedParton::xqCompanion() const --> double");
		cl.def("p", (class Pythia8::Vec4 (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::p, "C++: Pythia8::ResolvedParton::p() const --> class Pythia8::Vec4");
		cl.def("px", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::px, "C++: Pythia8::ResolvedParton::px() const --> double");
		cl.def("py", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::py, "C++: Pythia8::ResolvedParton::py() const --> double");
		cl.def("pz", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::pz, "C++: Pythia8::ResolvedParton::pz() const --> double");
		cl.def("e", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::e, "C++: Pythia8::ResolvedParton::e() const --> double");
		cl.def("m", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::m, "C++: Pythia8::ResolvedParton::m() const --> double");
		cl.def("pT", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::pT, "C++: Pythia8::ResolvedParton::pT() const --> double");
		cl.def("mT2", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::mT2, "C++: Pythia8::ResolvedParton::mT2() const --> double");
		cl.def("pPos", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::pPos, "C++: Pythia8::ResolvedParton::pPos() const --> double");
		cl.def("pNeg", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::pNeg, "C++: Pythia8::ResolvedParton::pNeg() const --> double");
		cl.def("col", (int (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::col, "C++: Pythia8::ResolvedParton::col() const --> int");
		cl.def("acol", (int (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::acol, "C++: Pythia8::ResolvedParton::acol() const --> int");
		cl.def("pTfactor", (double (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::pTfactor, "C++: Pythia8::ResolvedParton::pTfactor() const --> double");
		cl.def("hasCol", (bool (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::hasCol, "C++: Pythia8::ResolvedParton::hasCol() const --> bool");
		cl.def("hasAcol", (bool (Pythia8::ResolvedParton::*)() const) &Pythia8::ResolvedParton::hasAcol, "C++: Pythia8::ResolvedParton::hasAcol() const --> bool");
		cl.def("assign", (class Pythia8::ResolvedParton & (Pythia8::ResolvedParton::*)(const class Pythia8::ResolvedParton &)) &Pythia8::ResolvedParton::operator=, "C++: Pythia8::ResolvedParton::operator=(const class Pythia8::ResolvedParton &) --> class Pythia8::ResolvedParton &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
