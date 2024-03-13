#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/BoseEinstein.h>
#include <Pythia8/ColourTracing.h>
#include <Pythia8/DeuteronProduction.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/MiniStringFragmentation.h>
#include <Pythia8/NucleonExcitations.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/Ropewalk.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/Weights.h>
#include <cwchar>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <set>
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

// Pythia8::BeamShape file:Pythia8/BeamShape.h line:21
struct PyCallBack_Pythia8_BeamShape : public Pythia8::BeamShape {
	using Pythia8::BeamShape::BeamShape;

	void init(class Pythia8::Settings & a0, class Pythia8::Rndm * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamShape *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BeamShape::init(a0, a1);
	}
	void pick() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamShape *>(this), "pick");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BeamShape::pick();
	}
};

void bind_Pythia8_BeamShape(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::BeamShape file:Pythia8/BeamShape.h line:21
		pybind11::class_<Pythia8::BeamShape, std::shared_ptr<Pythia8::BeamShape>, PyCallBack_Pythia8_BeamShape> cl(M("Pythia8"), "BeamShape", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::BeamShape(); }, [](){ return new PyCallBack_Pythia8_BeamShape(); } ) );
		cl.def_readwrite("deltaPxA", &Pythia8::BeamShape::deltaPxA);
		cl.def_readwrite("deltaPyA", &Pythia8::BeamShape::deltaPyA);
		cl.def_readwrite("deltaPzA", &Pythia8::BeamShape::deltaPzA);
		cl.def_readwrite("deltaPxB", &Pythia8::BeamShape::deltaPxB);
		cl.def_readwrite("deltaPyB", &Pythia8::BeamShape::deltaPyB);
		cl.def_readwrite("deltaPzB", &Pythia8::BeamShape::deltaPzB);
		cl.def_readwrite("vertexX", &Pythia8::BeamShape::vertexX);
		cl.def_readwrite("vertexY", &Pythia8::BeamShape::vertexY);
		cl.def_readwrite("vertexZ", &Pythia8::BeamShape::vertexZ);
		cl.def_readwrite("vertexT", &Pythia8::BeamShape::vertexT);
		cl.def_readwrite("allowMomentumSpread", &Pythia8::BeamShape::allowMomentumSpread);
		cl.def_readwrite("allowVertexSpread", &Pythia8::BeamShape::allowVertexSpread);
		cl.def_readwrite("sigmaPxA", &Pythia8::BeamShape::sigmaPxA);
		cl.def_readwrite("sigmaPyA", &Pythia8::BeamShape::sigmaPyA);
		cl.def_readwrite("sigmaPzA", &Pythia8::BeamShape::sigmaPzA);
		cl.def_readwrite("maxDevA", &Pythia8::BeamShape::maxDevA);
		cl.def_readwrite("sigmaPxB", &Pythia8::BeamShape::sigmaPxB);
		cl.def_readwrite("sigmaPyB", &Pythia8::BeamShape::sigmaPyB);
		cl.def_readwrite("sigmaPzB", &Pythia8::BeamShape::sigmaPzB);
		cl.def_readwrite("maxDevB", &Pythia8::BeamShape::maxDevB);
		cl.def_readwrite("sigmaVertexX", &Pythia8::BeamShape::sigmaVertexX);
		cl.def_readwrite("sigmaVertexY", &Pythia8::BeamShape::sigmaVertexY);
		cl.def_readwrite("sigmaVertexZ", &Pythia8::BeamShape::sigmaVertexZ);
		cl.def_readwrite("maxDevVertex", &Pythia8::BeamShape::maxDevVertex);
		cl.def_readwrite("sigmaTime", &Pythia8::BeamShape::sigmaTime);
		cl.def_readwrite("maxDevTime", &Pythia8::BeamShape::maxDevTime);
		cl.def_readwrite("offsetX", &Pythia8::BeamShape::offsetX);
		cl.def_readwrite("offsetY", &Pythia8::BeamShape::offsetY);
		cl.def_readwrite("offsetZ", &Pythia8::BeamShape::offsetZ);
		cl.def_readwrite("offsetT", &Pythia8::BeamShape::offsetT);
		cl.def("init", (void (Pythia8::BeamShape::*)(class Pythia8::Settings &, class Pythia8::Rndm *)) &Pythia8::BeamShape::init, "C++: Pythia8::BeamShape::init(class Pythia8::Settings &, class Pythia8::Rndm *) --> void", pybind11::arg("settings"), pybind11::arg("rndmPtrIn"));
		cl.def("pick", (void (Pythia8::BeamShape::*)()) &Pythia8::BeamShape::pick, "C++: Pythia8::BeamShape::pick() --> void");
		cl.def("deltaPA", (class Pythia8::Vec4 (Pythia8::BeamShape::*)() const) &Pythia8::BeamShape::deltaPA, "C++: Pythia8::BeamShape::deltaPA() const --> class Pythia8::Vec4");
		cl.def("deltaPB", (class Pythia8::Vec4 (Pythia8::BeamShape::*)() const) &Pythia8::BeamShape::deltaPB, "C++: Pythia8::BeamShape::deltaPB() const --> class Pythia8::Vec4");
		cl.def("vertex", (class Pythia8::Vec4 (Pythia8::BeamShape::*)() const) &Pythia8::BeamShape::vertex, "C++: Pythia8::BeamShape::vertex() const --> class Pythia8::Vec4");
		cl.def("assign", (class Pythia8::BeamShape & (Pythia8::BeamShape::*)(const class Pythia8::BeamShape &)) &Pythia8::BeamShape::operator=, "C++: Pythia8::BeamShape::operator=(const class Pythia8::BeamShape &) --> class Pythia8::BeamShape &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
