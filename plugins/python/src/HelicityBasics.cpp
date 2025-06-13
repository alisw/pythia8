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

// Pythia8::HelicityParticle file:Pythia8/HelicityBasics.h line:182
struct PyCallBack_Pythia8_HelicityParticle : public Pythia8::HelicityParticle {
	using Pythia8::HelicityParticle::HelicityParticle;

	int index() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityParticle *>(this), "index");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return HelicityParticle::index();
	}
};

// Pythia8::HelicityMatrixElement file:Pythia8/HelicityMatrixElements.h line:24
struct PyCallBack_Pythia8_HelicityMatrixElement : public Pythia8::HelicityMatrixElement {
	using Pythia8::HelicityMatrixElement::HelicityMatrixElement;

	void initPointers(class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1, class Pythia8::Settings * a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "initPointers");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "initChannel");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "decayWeight");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "decayWeightMax");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "calculateME");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "calculateD");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "calculateRho");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "breitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "sBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "pBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "dBreitWigner");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "initConstants");
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
	void initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HelicityMatrixElement *>(this), "initWaves");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return HelicityMatrixElement::initWaves(a0);
	}
};

void bind_Pythia8_HelicityBasics(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Pythia8::conj(class Pythia8::Wave4) file:Pythia8/HelicityBasics.h line:77
	M("Pythia8").def("conj", (class Pythia8::Wave4 (*)(class Pythia8::Wave4)) &Pythia8::conj, "C++: Pythia8::conj(class Pythia8::Wave4) --> class Pythia8::Wave4", pybind11::arg("w"));

	// Pythia8::epsilon(class Pythia8::Wave4, class Pythia8::Wave4, class Pythia8::Wave4) file:Pythia8/HelicityBasics.h line:80
	M("Pythia8").def("epsilon", (class Pythia8::Wave4 (*)(class Pythia8::Wave4, class Pythia8::Wave4, class Pythia8::Wave4)) &Pythia8::epsilon, "C++: Pythia8::epsilon(class Pythia8::Wave4, class Pythia8::Wave4, class Pythia8::Wave4) --> class Pythia8::Wave4", pybind11::arg("w1"), pybind11::arg("w2"), pybind11::arg("w3"));

	// Pythia8::m2(class Pythia8::Wave4) file:Pythia8/HelicityBasics.h line:83
	M("Pythia8").def("m2", (double (*)(class Pythia8::Wave4)) &Pythia8::m2, "C++: Pythia8::m2(class Pythia8::Wave4) --> double", pybind11::arg("w"));

	// Pythia8::m2(class Pythia8::Wave4, class Pythia8::Wave4) file:Pythia8/HelicityBasics.h line:84
	M("Pythia8").def("m2", (double (*)(class Pythia8::Wave4, class Pythia8::Wave4)) &Pythia8::m2, "C++: Pythia8::m2(class Pythia8::Wave4, class Pythia8::Wave4) --> double", pybind11::arg("w1"), pybind11::arg("w2"));

	{ // Pythia8::GammaMatrix file:Pythia8/HelicityBasics.h line:118
		pybind11::class_<Pythia8::GammaMatrix, std::shared_ptr<Pythia8::GammaMatrix>> cl(M("Pythia8"), "GammaMatrix", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::GammaMatrix(); } ) );
		cl.def( pybind11::init<int>(), pybind11::arg("mu") );

		cl.def( pybind11::init( [](Pythia8::GammaMatrix const &o){ return new Pythia8::GammaMatrix(o); } ) );
		cl.def_readwrite("COMPLEXZERO", &Pythia8::GammaMatrix::COMPLEXZERO);
		cl.def("__call__", (struct std::complex<double> & (Pythia8::GammaMatrix::*)(int, int)) &Pythia8::GammaMatrix::operator(), "C++: Pythia8::GammaMatrix::operator()(int, int) --> struct std::complex<double> &", pybind11::return_value_policy::reference, pybind11::arg("I"), pybind11::arg("J"));
		cl.def("__mul__", (class Pythia8::GammaMatrix (Pythia8::GammaMatrix::*)(struct std::complex<double>)) &Pythia8::GammaMatrix::operator*, "C++: Pythia8::GammaMatrix::operator*(struct std::complex<double>) --> class Pythia8::GammaMatrix", pybind11::arg("s"));
		cl.def("__sub__", (class Pythia8::GammaMatrix (Pythia8::GammaMatrix::*)(struct std::complex<double>)) &Pythia8::GammaMatrix::operator-, "C++: Pythia8::GammaMatrix::operator-(struct std::complex<double>) --> class Pythia8::GammaMatrix", pybind11::arg("s"));
		cl.def("__add__", (class Pythia8::GammaMatrix (Pythia8::GammaMatrix::*)(struct std::complex<double>)) &Pythia8::GammaMatrix::operator+, "C++: Pythia8::GammaMatrix::operator+(struct std::complex<double>) --> class Pythia8::GammaMatrix", pybind11::arg("s"));
	}
	{ // Pythia8::HelicityParticle file:Pythia8/HelicityBasics.h line:182
		pybind11::class_<Pythia8::HelicityParticle, std::shared_ptr<Pythia8::HelicityParticle>, PyCallBack_Pythia8_HelicityParticle, Pythia8::Particle> cl(M("Pythia8"), "HelicityParticle", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HelicityParticle(); }, [](){ return new PyCallBack_Pythia8_HelicityParticle(); } ) );
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::HelicityParticle(a0); }, [](int const & a0){ return new PyCallBack_Pythia8_HelicityParticle(a0); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Pythia8::HelicityParticle(a0, a1); }, [](int const & a0, int const & a1){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2){ return new Pythia8::HelicityParticle(a0, a1, a2); }, [](int const & a0, int const & a1, int const & a2){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3){ return new Pythia8::HelicityParticle(a0, a1, a2, a3); }, [](int const & a0, int const & a1, int const & a2, int const & a3){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9, double const & a10){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9, double const & a10){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9, double const & a10, double const & a11){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9, double const & a10, double const & a11){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9, double const & a10, double const & a11, double const & a12){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9, double const & a10, double const & a11, double const & a12){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9, double const & a10, double const & a11, double const & a12, double const & a13){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, double const & a8, double const & a9, double const & a10, double const & a11, double const & a12, double const & a13){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); } ), "doc");
		cl.def( pybind11::init<int, int, int, int, int, int, int, int, double, double, double, double, double, double, class Pythia8::ParticleData *>(), pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("daughter1In"), pybind11::arg("daughter2In"), pybind11::arg("colIn"), pybind11::arg("acolIn"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("pzIn"), pybind11::arg("eIn"), pybind11::arg("mIn"), pybind11::arg("scaleIn"), pybind11::arg("ptr") );

		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, class Pythia8::Vec4 const & a8){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, class Pythia8::Vec4 const & a8){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, class Pythia8::Vec4 const & a8, double const & a9){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, class Pythia8::Vec4 const & a8, double const & a9){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); } ), "doc");
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, class Pythia8::Vec4 const & a8, double const & a9, double const & a10){ return new Pythia8::HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, [](int const & a0, int const & a1, int const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, class Pythia8::Vec4 const & a8, double const & a9, double const & a10){ return new PyCallBack_Pythia8_HelicityParticle(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); } ), "doc");
		cl.def( pybind11::init<int, int, int, int, int, int, int, int, class Pythia8::Vec4, double, double, class Pythia8::ParticleData *>(), pybind11::arg("idIn"), pybind11::arg("statusIn"), pybind11::arg("mother1In"), pybind11::arg("mother2In"), pybind11::arg("daughter1In"), pybind11::arg("daughter2In"), pybind11::arg("colIn"), pybind11::arg("acolIn"), pybind11::arg("pIn"), pybind11::arg("mIn"), pybind11::arg("scaleIn"), pybind11::arg("ptr") );

		cl.def( pybind11::init( [](const class Pythia8::Particle & a0){ return new Pythia8::HelicityParticle(a0); }, [](const class Pythia8::Particle & a0){ return new PyCallBack_Pythia8_HelicityParticle(a0); } ), "doc");
		cl.def( pybind11::init<const class Pythia8::Particle &, class Pythia8::ParticleData *>(), pybind11::arg("ptIn"), pybind11::arg("ptr") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_HelicityParticle const &o){ return new PyCallBack_Pythia8_HelicityParticle(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HelicityParticle const &o){ return new Pythia8::HelicityParticle(o); } ) );
		cl.def_readwrite("direction", &Pythia8::HelicityParticle::direction);
		cl.def_readwrite("rho", &Pythia8::HelicityParticle::rho);
		cl.def_readwrite("D", &Pythia8::HelicityParticle::D);
		cl.def("wave", (class Pythia8::Wave4 (Pythia8::HelicityParticle::*)(int)) &Pythia8::HelicityParticle::wave, "C++: Pythia8::HelicityParticle::wave(int) --> class Pythia8::Wave4", pybind11::arg("h"));
		cl.def("waveBar", (class Pythia8::Wave4 (Pythia8::HelicityParticle::*)(int)) &Pythia8::HelicityParticle::waveBar, "C++: Pythia8::HelicityParticle::waveBar(int) --> class Pythia8::Wave4", pybind11::arg("h"));
		cl.def("normalize", (void (Pythia8::HelicityParticle::*)(class std::vector<class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > >, class std::allocator<class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > > > &)) &Pythia8::HelicityParticle::normalize, "C++: Pythia8::HelicityParticle::normalize(class std::vector<class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > >, class std::allocator<class std::vector<struct std::complex<double>, class std::allocator<struct std::complex<double> > > > > &) --> void", pybind11::arg("m"));
		cl.def("spinStates", (int (Pythia8::HelicityParticle::*)()) &Pythia8::HelicityParticle::spinStates, "C++: Pythia8::HelicityParticle::spinStates() --> int");
		cl.def("m", (double (Pythia8::HelicityParticle::*)()) &Pythia8::HelicityParticle::m, "C++: Pythia8::HelicityParticle::m() --> double");
		cl.def("m", (void (Pythia8::HelicityParticle::*)(double)) &Pythia8::HelicityParticle::m, "C++: Pythia8::HelicityParticle::m(double) --> void", pybind11::arg("mIn"));
		cl.def("pol", (double (Pythia8::HelicityParticle::*)()) &Pythia8::HelicityParticle::pol, "C++: Pythia8::HelicityParticle::pol() --> double");
		cl.def("pol", (void (Pythia8::HelicityParticle::*)(double)) &Pythia8::HelicityParticle::pol, "C++: Pythia8::HelicityParticle::pol(double) --> void", pybind11::arg("hIn"));
		cl.def("index", (int (Pythia8::HelicityParticle::*)() const) &Pythia8::HelicityParticle::index, "C++: Pythia8::HelicityParticle::index() const --> int");
		cl.def("index", (void (Pythia8::HelicityParticle::*)(int)) &Pythia8::HelicityParticle::index, "C++: Pythia8::HelicityParticle::index(int) --> void", pybind11::arg("indexIn"));
		cl.def("assign", (class Pythia8::HelicityParticle & (Pythia8::HelicityParticle::*)(const class Pythia8::HelicityParticle &)) &Pythia8::HelicityParticle::operator=, "C++: Pythia8::HelicityParticle::operator=(const class Pythia8::HelicityParticle &) --> class Pythia8::HelicityParticle &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HelicityMatrixElement file:Pythia8/HelicityMatrixElements.h line:24
		pybind11::class_<Pythia8::HelicityMatrixElement, std::shared_ptr<Pythia8::HelicityMatrixElement>, PyCallBack_Pythia8_HelicityMatrixElement> cl(M("Pythia8"), "HelicityMatrixElement", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HelicityMatrixElement(); }, [](){ return new PyCallBack_Pythia8_HelicityMatrixElement(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HelicityMatrixElement const &o){ return new PyCallBack_Pythia8_HelicityMatrixElement(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HelicityMatrixElement const &o){ return new Pythia8::HelicityMatrixElement(o); } ) );
		cl.def_readwrite("DECAYWEIGHTMAX", &Pythia8::HelicityMatrixElement::DECAYWEIGHTMAX);
		cl.def_readwrite("gamma", &Pythia8::HelicityMatrixElement::gamma);
		cl.def_readwrite("pMap", &Pythia8::HelicityMatrixElement::pMap);
		cl.def_readwrite("pID", &Pythia8::HelicityMatrixElement::pID);
		cl.def_readwrite("pM", &Pythia8::HelicityMatrixElement::pM);
		cl.def_readwrite("u", &Pythia8::HelicityMatrixElement::u);
		cl.def("initPointers", [](Pythia8::HelicityMatrixElement &o, class Pythia8::ParticleData * a0, class Pythia8::CoupSM * a1) -> void { return o.initPointers(a0, a1); }, "", pybind11::arg(""), pybind11::arg(""));
		cl.def("initPointers", (void (Pythia8::HelicityMatrixElement::*)(class Pythia8::ParticleData *, class Pythia8::CoupSM *, class Pythia8::Settings *)) &Pythia8::HelicityMatrixElement::initPointers, "C++: Pythia8::HelicityMatrixElement::initPointers(class Pythia8::ParticleData *, class Pythia8::CoupSM *, class Pythia8::Settings *) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("initChannel", (class Pythia8::HelicityMatrixElement * (Pythia8::HelicityMatrixElement::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HelicityMatrixElement::initChannel, "C++: Pythia8::HelicityMatrixElement::initChannel(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> class Pythia8::HelicityMatrixElement *", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("decayWeight", (double (Pythia8::HelicityMatrixElement::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HelicityMatrixElement::decayWeight, "C++: Pythia8::HelicityMatrixElement::decayWeight(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> double", pybind11::arg(""));
		cl.def("decayWeightMax", (double (Pythia8::HelicityMatrixElement::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HelicityMatrixElement::decayWeightMax, "C++: Pythia8::HelicityMatrixElement::decayWeightMax(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> double", pybind11::arg(""));
		cl.def("calculateME", (struct std::complex<double> (Pythia8::HelicityMatrixElement::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::HelicityMatrixElement::calculateME, "C++: Pythia8::HelicityMatrixElement::calculateME(class std::vector<int, class std::allocator<int> >) --> struct std::complex<double>", pybind11::arg(""));
		cl.def("calculateD", (void (Pythia8::HelicityMatrixElement::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HelicityMatrixElement::calculateD, "C++: Pythia8::HelicityMatrixElement::calculateD(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("calculateRho", (void (Pythia8::HelicityMatrixElement::*)(unsigned int, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HelicityMatrixElement::calculateRho, "C++: Pythia8::HelicityMatrixElement::calculateRho(unsigned int, class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""), pybind11::arg(""));
		cl.def("setFermionLine", (void (Pythia8::HelicityMatrixElement::*)(int, class Pythia8::HelicityParticle &, class Pythia8::HelicityParticle &)) &Pythia8::HelicityMatrixElement::setFermionLine, "C++: Pythia8::HelicityMatrixElement::setFermionLine(int, class Pythia8::HelicityParticle &, class Pythia8::HelicityParticle &) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("breitWigner", (struct std::complex<double> (Pythia8::HelicityMatrixElement::*)(double, double, double)) &Pythia8::HelicityMatrixElement::breitWigner, "C++: Pythia8::HelicityMatrixElement::breitWigner(double, double, double) --> struct std::complex<double>", pybind11::arg("s"), pybind11::arg("M"), pybind11::arg("G"));
		cl.def("sBreitWigner", (struct std::complex<double> (Pythia8::HelicityMatrixElement::*)(double, double, double, double, double)) &Pythia8::HelicityMatrixElement::sBreitWigner, "C++: Pythia8::HelicityMatrixElement::sBreitWigner(double, double, double, double, double) --> struct std::complex<double>", pybind11::arg("m0"), pybind11::arg("m1"), pybind11::arg("s"), pybind11::arg("M"), pybind11::arg("G"));
		cl.def("pBreitWigner", (struct std::complex<double> (Pythia8::HelicityMatrixElement::*)(double, double, double, double, double)) &Pythia8::HelicityMatrixElement::pBreitWigner, "C++: Pythia8::HelicityMatrixElement::pBreitWigner(double, double, double, double, double) --> struct std::complex<double>", pybind11::arg("m0"), pybind11::arg("m1"), pybind11::arg("s"), pybind11::arg("M"), pybind11::arg("G"));
		cl.def("dBreitWigner", (struct std::complex<double> (Pythia8::HelicityMatrixElement::*)(double, double, double, double, double)) &Pythia8::HelicityMatrixElement::dBreitWigner, "C++: Pythia8::HelicityMatrixElement::dBreitWigner(double, double, double, double, double) --> struct std::complex<double>", pybind11::arg("m0"), pybind11::arg("m1"), pybind11::arg("s"), pybind11::arg("M"), pybind11::arg("G"));
		cl.def("initConstants", (void (Pythia8::HelicityMatrixElement::*)()) &Pythia8::HelicityMatrixElement::initConstants, "C++: Pythia8::HelicityMatrixElement::initConstants() --> void");
		cl.def("initWaves", (void (Pythia8::HelicityMatrixElement::*)(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &)) &Pythia8::HelicityMatrixElement::initWaves, "C++: Pythia8::HelicityMatrixElement::initWaves(class std::vector<class Pythia8::HelicityParticle, class std::allocator<class Pythia8::HelicityParticle> > &) --> void", pybind11::arg(""));
		cl.def("assign", (class Pythia8::HelicityMatrixElement & (Pythia8::HelicityMatrixElement::*)(const class Pythia8::HelicityMatrixElement &)) &Pythia8::HelicityMatrixElement::operator=, "C++: Pythia8::HelicityMatrixElement::operator=(const class Pythia8::HelicityMatrixElement &) --> class Pythia8::HelicityMatrixElement &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
