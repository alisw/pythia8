#include <Pythia8/Basics.h>
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
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
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

// Pythia8::BoseEinstein file:Pythia8/BoseEinstein.h line:48
struct PyCallBack_Pythia8_BoseEinstein : public Pythia8::BoseEinstein {
	using Pythia8::BoseEinstein::BoseEinstein;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BoseEinstein *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BoseEinstein *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BoseEinstein *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BoseEinstein *>(this), "onStat");
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

// Pythia8::DeuteronProduction file:Pythia8/DeuteronProduction.h line:22
struct PyCallBack_Pythia8_DeuteronProduction : public Pythia8::DeuteronProduction {
	using Pythia8::DeuteronProduction::DeuteronProduction;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DeuteronProduction *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DeuteronProduction *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DeuteronProduction *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::DeuteronProduction *>(this), "onStat");
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
	{ // Pythia8::BoseEinsteinHadron file:Pythia8/BoseEinstein.h line:26
		pybind11::class_<Pythia8::BoseEinsteinHadron, std::shared_ptr<Pythia8::BoseEinsteinHadron>> cl(M("Pythia8"), "BoseEinsteinHadron", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::BoseEinsteinHadron(); } ) );
		cl.def( pybind11::init<int, int, class Pythia8::Vec4, double>(), pybind11::arg("idIn"), pybind11::arg("iPosIn"), pybind11::arg("pIn"), pybind11::arg("mIn") );

		cl.def_readwrite("id", &Pythia8::BoseEinsteinHadron::id);
		cl.def_readwrite("iPos", &Pythia8::BoseEinsteinHadron::iPos);
		cl.def_readwrite("p", &Pythia8::BoseEinsteinHadron::p);
		cl.def_readwrite("pShift", &Pythia8::BoseEinsteinHadron::pShift);
		cl.def_readwrite("pComp", &Pythia8::BoseEinsteinHadron::pComp);
		cl.def_readwrite("m2", &Pythia8::BoseEinsteinHadron::m2);
	}
	{ // Pythia8::BoseEinstein file:Pythia8/BoseEinstein.h line:48
		pybind11::class_<Pythia8::BoseEinstein, std::shared_ptr<Pythia8::BoseEinstein>, PyCallBack_Pythia8_BoseEinstein, Pythia8::PhysicsBase> cl(M("Pythia8"), "BoseEinstein", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::BoseEinstein(); }, [](){ return new PyCallBack_Pythia8_BoseEinstein(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_BoseEinstein const &o){ return new PyCallBack_Pythia8_BoseEinstein(o); } ) );
		cl.def( pybind11::init( [](Pythia8::BoseEinstein const &o){ return new Pythia8::BoseEinstein(o); } ) );
		cl.def("init", (bool (Pythia8::BoseEinstein::*)()) &Pythia8::BoseEinstein::init, "C++: Pythia8::BoseEinstein::init() --> bool");
		cl.def("shiftEvent", (bool (Pythia8::BoseEinstein::*)(class Pythia8::Event &)) &Pythia8::BoseEinstein::shiftEvent, "C++: Pythia8::BoseEinstein::shiftEvent(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("assign", (class Pythia8::BoseEinstein & (Pythia8::BoseEinstein::*)(const class Pythia8::BoseEinstein &)) &Pythia8::BoseEinstein::operator=, "C++: Pythia8::BoseEinstein::operator=(const class Pythia8::BoseEinstein &) --> class Pythia8::BoseEinstein &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ColourTracing file:Pythia8/ColourTracing.h line:22
		pybind11::class_<Pythia8::ColourTracing, std::shared_ptr<Pythia8::ColourTracing>> cl(M("Pythia8"), "ColourTracing", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ColourTracing(); } ) );
		cl.def( pybind11::init( [](Pythia8::ColourTracing const &o){ return new Pythia8::ColourTracing(o); } ) );
		cl.def("setupColList", (bool (Pythia8::ColourTracing::*)(class Pythia8::Event &)) &Pythia8::ColourTracing::setupColList, "C++: Pythia8::ColourTracing::setupColList(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("traceFromAcol", (bool (Pythia8::ColourTracing::*)(int, class Pythia8::Event &, int, int, class std::vector<int, class std::allocator<int> > &)) &Pythia8::ColourTracing::traceFromAcol, "C++: Pythia8::ColourTracing::traceFromAcol(int, class Pythia8::Event &, int, int, class std::vector<int, class std::allocator<int> > &) --> bool", pybind11::arg("indxCol"), pybind11::arg("event"), pybind11::arg("iJun"), pybind11::arg("iCol"), pybind11::arg("iParton"));
		cl.def("traceFromCol", (bool (Pythia8::ColourTracing::*)(int, class Pythia8::Event &, int, int, class std::vector<int, class std::allocator<int> > &)) &Pythia8::ColourTracing::traceFromCol, "C++: Pythia8::ColourTracing::traceFromCol(int, class Pythia8::Event &, int, int, class std::vector<int, class std::allocator<int> > &) --> bool", pybind11::arg("indxCol"), pybind11::arg("event"), pybind11::arg("iJun"), pybind11::arg("iCol"), pybind11::arg("iParton"));
		cl.def("traceInLoop", (bool (Pythia8::ColourTracing::*)(class Pythia8::Event &, class std::vector<int, class std::allocator<int> > &)) &Pythia8::ColourTracing::traceInLoop, "C++: Pythia8::ColourTracing::traceInLoop(class Pythia8::Event &, class std::vector<int, class std::allocator<int> > &) --> bool", pybind11::arg("event"), pybind11::arg("iParton"));
		cl.def("finished", (bool (Pythia8::ColourTracing::*)()) &Pythia8::ColourTracing::finished, "C++: Pythia8::ColourTracing::finished() --> bool");
		cl.def("colFinished", (bool (Pythia8::ColourTracing::*)()) &Pythia8::ColourTracing::colFinished, "C++: Pythia8::ColourTracing::colFinished() --> bool");
		cl.def("getJunChains", (class std::vector<class std::vector<int, class std::allocator<int> >, class std::allocator<class std::vector<int, class std::allocator<int> > > > (Pythia8::ColourTracing::*)(class Pythia8::Event &)) &Pythia8::ColourTracing::getJunChains, "C++: Pythia8::ColourTracing::getJunChains(class Pythia8::Event &) --> class std::vector<class std::vector<int, class std::allocator<int> >, class std::allocator<class std::vector<int, class std::allocator<int> > > >", pybind11::arg("event"));
		cl.def("assign", (class Pythia8::ColourTracing & (Pythia8::ColourTracing::*)(const class Pythia8::ColourTracing &)) &Pythia8::ColourTracing::operator=, "C++: Pythia8::ColourTracing::operator=(const class Pythia8::ColourTracing &) --> class Pythia8::ColourTracing &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::DeuteronProduction file:Pythia8/DeuteronProduction.h line:22
		pybind11::class_<Pythia8::DeuteronProduction, std::shared_ptr<Pythia8::DeuteronProduction>, PyCallBack_Pythia8_DeuteronProduction, Pythia8::PhysicsBase> cl(M("Pythia8"), "DeuteronProduction", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::DeuteronProduction(); }, [](){ return new PyCallBack_Pythia8_DeuteronProduction(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_DeuteronProduction const &o){ return new PyCallBack_Pythia8_DeuteronProduction(o); } ) );
		cl.def( pybind11::init( [](Pythia8::DeuteronProduction const &o){ return new Pythia8::DeuteronProduction(o); } ) );
		cl.def("init", (bool (Pythia8::DeuteronProduction::*)()) &Pythia8::DeuteronProduction::init, "C++: Pythia8::DeuteronProduction::init() --> bool");
		cl.def("combine", (bool (Pythia8::DeuteronProduction::*)(class Pythia8::Event &)) &Pythia8::DeuteronProduction::combine, "C++: Pythia8::DeuteronProduction::combine(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("assign", (class Pythia8::DeuteronProduction & (Pythia8::DeuteronProduction::*)(const class Pythia8::DeuteronProduction &)) &Pythia8::DeuteronProduction::operator=, "C++: Pythia8::DeuteronProduction::operator=(const class Pythia8::DeuteronProduction &) --> class Pythia8::DeuteronProduction &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ColSinglet file:Pythia8/FragmentationSystems.h line:31
		pybind11::class_<Pythia8::ColSinglet, std::shared_ptr<Pythia8::ColSinglet>> cl(M("Pythia8"), "ColSinglet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ColSinglet(); } ) );
		cl.def( pybind11::init( [](class std::vector<int, class std::allocator<int> > & a0, class Pythia8::Vec4 const & a1, double const & a2, double const & a3){ return new Pythia8::ColSinglet(a0, a1, a2, a3); } ), "doc" , pybind11::arg("iPartonIn"), pybind11::arg("pSumIn"), pybind11::arg("massIn"), pybind11::arg("massExcessIn"));
		cl.def( pybind11::init( [](class std::vector<int, class std::allocator<int> > & a0, class Pythia8::Vec4 const & a1, double const & a2, double const & a3, bool const & a4){ return new Pythia8::ColSinglet(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("iPartonIn"), pybind11::arg("pSumIn"), pybind11::arg("massIn"), pybind11::arg("massExcessIn"), pybind11::arg("hasJunctionIn"));
		cl.def( pybind11::init( [](class std::vector<int, class std::allocator<int> > & a0, class Pythia8::Vec4 const & a1, double const & a2, double const & a3, bool const & a4, bool const & a5){ return new Pythia8::ColSinglet(a0, a1, a2, a3, a4, a5); } ), "doc" , pybind11::arg("iPartonIn"), pybind11::arg("pSumIn"), pybind11::arg("massIn"), pybind11::arg("massExcessIn"), pybind11::arg("hasJunctionIn"), pybind11::arg("isClosedIn"));
		cl.def( pybind11::init<class std::vector<int, class std::allocator<int> > &, class Pythia8::Vec4, double, double, bool, bool, bool>(), pybind11::arg("iPartonIn"), pybind11::arg("pSumIn"), pybind11::arg("massIn"), pybind11::arg("massExcessIn"), pybind11::arg("hasJunctionIn"), pybind11::arg("isClosedIn"), pybind11::arg("isCollectedIn") );

		cl.def( pybind11::init( [](Pythia8::ColSinglet const &o){ return new Pythia8::ColSinglet(o); } ) );
		cl.def_readwrite("iParton", &Pythia8::ColSinglet::iParton);
		cl.def_readwrite("pSum", &Pythia8::ColSinglet::pSum);
		cl.def_readwrite("mass", &Pythia8::ColSinglet::mass);
		cl.def_readwrite("massExcess", &Pythia8::ColSinglet::massExcess);
		cl.def_readwrite("hasJunction", &Pythia8::ColSinglet::hasJunction);
		cl.def_readwrite("isClosed", &Pythia8::ColSinglet::isClosed);
		cl.def_readwrite("isCollected", &Pythia8::ColSinglet::isCollected);
		cl.def("size", (int (Pythia8::ColSinglet::*)() const) &Pythia8::ColSinglet::size, "C++: Pythia8::ColSinglet::size() const --> int");
		cl.def("assign", (class Pythia8::ColSinglet & (Pythia8::ColSinglet::*)(const class Pythia8::ColSinglet &)) &Pythia8::ColSinglet::operator=, "C++: Pythia8::ColSinglet::operator=(const class Pythia8::ColSinglet &) --> class Pythia8::ColSinglet &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ColConfig file:Pythia8/FragmentationSystems.h line:60
		pybind11::class_<Pythia8::ColConfig, std::shared_ptr<Pythia8::ColConfig>> cl(M("Pythia8"), "ColConfig", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ColConfig(); } ) );
		cl.def( pybind11::init( [](Pythia8::ColConfig const &o){ return new Pythia8::ColConfig(o); } ) );
		cl.def_readwrite("rapPairs", &Pythia8::ColConfig::rapPairs);
		cl.def("init", [](Pythia8::ColConfig &o, class Pythia8::Info * a0, class Pythia8::StringFlav * a1) -> void { return o.init(a0, a1); }, "", pybind11::arg("infoPtrIn"), pybind11::arg("flavSelPtrIn"));
		cl.def("init", (void (Pythia8::ColConfig::*)(class Pythia8::Info *, class Pythia8::StringFlav *, double)) &Pythia8::ColConfig::init, "C++: Pythia8::ColConfig::init(class Pythia8::Info *, class Pythia8::StringFlav *, double) --> void", pybind11::arg("infoPtrIn"), pybind11::arg("flavSelPtrIn"), pybind11::arg("mVecRatio"));
		cl.def("size", (int (Pythia8::ColConfig::*)() const) &Pythia8::ColConfig::size, "C++: Pythia8::ColConfig::size() const --> int");
		cl.def("__getitem__", (class Pythia8::ColSinglet & (Pythia8::ColConfig::*)(int)) &Pythia8::ColConfig::operator[], "C++: Pythia8::ColConfig::operator[](int) --> class Pythia8::ColSinglet &", pybind11::return_value_policy::reference, pybind11::arg("iSub"));
		cl.def("clear", (void (Pythia8::ColConfig::*)()) &Pythia8::ColConfig::clear, "C++: Pythia8::ColConfig::clear() --> void");
		cl.def("insert", (bool (Pythia8::ColConfig::*)(class std::vector<int, class std::allocator<int> > &, class Pythia8::Event &)) &Pythia8::ColConfig::insert, "C++: Pythia8::ColConfig::insert(class std::vector<int, class std::allocator<int> > &, class Pythia8::Event &) --> bool", pybind11::arg("iPartonIn"), pybind11::arg("event"));
		cl.def("simpleInsert", [](Pythia8::ColConfig &o, class std::vector<int, class std::allocator<int> > & a0, class Pythia8::Event & a1) -> bool { return o.simpleInsert(a0, a1); }, "", pybind11::arg("iPartonIn"), pybind11::arg("event"));
		cl.def("simpleInsert", (bool (Pythia8::ColConfig::*)(class std::vector<int, class std::allocator<int> > &, class Pythia8::Event &, bool)) &Pythia8::ColConfig::simpleInsert, "C++: Pythia8::ColConfig::simpleInsert(class std::vector<int, class std::allocator<int> > &, class Pythia8::Event &, bool) --> bool", pybind11::arg("iPartonIn"), pybind11::arg("event"), pybind11::arg("fixOrder"));
		cl.def("erase", (void (Pythia8::ColConfig::*)(int)) &Pythia8::ColConfig::erase, "C++: Pythia8::ColConfig::erase(int) --> void", pybind11::arg("iSub"));
		cl.def("collect", [](Pythia8::ColConfig &o, int const & a0, class Pythia8::Event & a1) -> void { return o.collect(a0, a1); }, "", pybind11::arg("iSub"), pybind11::arg("event"));
		cl.def("collect", (void (Pythia8::ColConfig::*)(int, class Pythia8::Event &, bool)) &Pythia8::ColConfig::collect, "C++: Pythia8::ColConfig::collect(int, class Pythia8::Event &, bool) --> void", pybind11::arg("iSub"), pybind11::arg("event"), pybind11::arg("skipTrivial"));
		cl.def("findSinglet", (int (Pythia8::ColConfig::*)(int)) &Pythia8::ColConfig::findSinglet, "C++: Pythia8::ColConfig::findSinglet(int) --> int", pybind11::arg("i"));
		cl.def("list", (void (Pythia8::ColConfig::*)() const) &Pythia8::ColConfig::list, "C++: Pythia8::ColConfig::list() const --> void");
		cl.def("assign", (class Pythia8::ColConfig & (Pythia8::ColConfig::*)(const class Pythia8::ColConfig &)) &Pythia8::ColConfig::operator=, "C++: Pythia8::ColConfig::operator=(const class Pythia8::ColConfig &) --> class Pythia8::ColConfig &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
