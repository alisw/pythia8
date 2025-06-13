#include <Pythia8/Basics.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/Info.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/Merging.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonLevel.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/RHadrons.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Ropewalk.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SpaceShower.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/TimeShower.h>
#include <istream>
#include <iterator>
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

// Pythia8::Merging file:Pythia8/Merging.h line:33
struct PyCallBack_Pythia8_Merging : public Pythia8::Merging {
	using Pythia8::Merging::Merging;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Merging *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Merging::init();
	}
	void statistics() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Merging *>(this), "statistics");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Merging::statistics();
	}
	int mergeProcess(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Merging *>(this), "mergeProcess");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return Merging::mergeProcess(a0);
	}
	double generateSingleSudakov(double a0, double a1, double a2, int a3, int a4, double a5, double a6) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Merging *>(this), "generateSingleSudakov");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Merging::generateSingleSudakov(a0, a1, a2, a3, a4, a5, a6);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Merging *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Merging *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Merging *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Merging *>(this), "onStat");
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

// Pythia8::Ropewalk file:Pythia8/Ropewalk.h line:211
struct PyCallBack_Pythia8_Ropewalk : public Pythia8::Ropewalk {
	using Pythia8::Ropewalk::Ropewalk;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Ropewalk *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Ropewalk::init();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Ropewalk *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Ropewalk *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Ropewalk *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::Ropewalk *>(this), "onStat");
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

void bind_Pythia8_Merging(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::Merging file:Pythia8/Merging.h line:33
		pybind11::class_<Pythia8::Merging, std::shared_ptr<Pythia8::Merging>, PyCallBack_Pythia8_Merging, Pythia8::PhysicsBase> cl(M("Pythia8"), "Merging", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Merging(); }, [](){ return new PyCallBack_Pythia8_Merging(); } ) );
		cl.def_readwrite("lhaPtr", &Pythia8::Merging::lhaPtr);
		cl.def_readwrite("mergingHooksPtr", &Pythia8::Merging::mergingHooksPtr);
		cl.def_readwrite("tmsNowMin", &Pythia8::Merging::tmsNowMin);
		cl.def_readwrite("stoppingScalesSave", &Pythia8::Merging::stoppingScalesSave);
		cl.def_readwrite("mDipSave", &Pythia8::Merging::mDipSave);
		cl.def_readwrite("radSave", &Pythia8::Merging::radSave);
		cl.def_readwrite("emtSave", &Pythia8::Merging::emtSave);
		cl.def_readwrite("recSave", &Pythia8::Merging::recSave);
		cl.def_readwrite("isInDeadzone", &Pythia8::Merging::isInDeadzone);
		cl.def("initPtrs", (void (Pythia8::Merging::*)(class std::shared_ptr<class Pythia8::MergingHooks>, class Pythia8::PartonLevel *)) &Pythia8::Merging::initPtrs, "C++: Pythia8::Merging::initPtrs(class std::shared_ptr<class Pythia8::MergingHooks>, class Pythia8::PartonLevel *) --> void", pybind11::arg("mergingHooksPtrIn"), pybind11::arg("trialPartonLevelPtrIn"));
		cl.def("init", (void (Pythia8::Merging::*)()) &Pythia8::Merging::init, "C++: Pythia8::Merging::init() --> void");
		cl.def("statistics", (void (Pythia8::Merging::*)()) &Pythia8::Merging::statistics, "C++: Pythia8::Merging::statistics() --> void");
		cl.def("mergeProcess", (int (Pythia8::Merging::*)(class Pythia8::Event &)) &Pythia8::Merging::mergeProcess, "C++: Pythia8::Merging::mergeProcess(class Pythia8::Event &) --> int", pybind11::arg("process"));
		cl.def("generateSingleSudakov", [](Pythia8::Merging &o, double const & a0, double const & a1, double const & a2, int const & a3, int const & a4) -> double { return o.generateSingleSudakov(a0, a1, a2, a3, a4); }, "", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("idA"), pybind11::arg("type"));
		cl.def("generateSingleSudakov", [](Pythia8::Merging &o, double const & a0, double const & a1, double const & a2, int const & a3, int const & a4, double const & a5) -> double { return o.generateSingleSudakov(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("idA"), pybind11::arg("type"), pybind11::arg("s"));
		cl.def("generateSingleSudakov", (double (Pythia8::Merging::*)(double, double, double, int, int, double, double)) &Pythia8::Merging::generateSingleSudakov, "C++: Pythia8::Merging::generateSingleSudakov(double, double, double, int, int, double, double) --> double", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("idA"), pybind11::arg("type"), pybind11::arg("s"), pybind11::arg("x"));
		cl.def("setLHAPtr", (void (Pythia8::Merging::*)(class std::shared_ptr<class Pythia8::LHEF3FromPythia8>)) &Pythia8::Merging::setLHAPtr, "C++: Pythia8::Merging::setLHAPtr(class std::shared_ptr<class Pythia8::LHEF3FromPythia8>) --> void", pybind11::arg("lhaUpIn"));
		cl.def("mergeProcessCKKWL", (int (Pythia8::Merging::*)(class Pythia8::Event &)) &Pythia8::Merging::mergeProcessCKKWL, "C++: Pythia8::Merging::mergeProcessCKKWL(class Pythia8::Event &) --> int", pybind11::arg("process"));
		cl.def("mergeProcessUMEPS", (int (Pythia8::Merging::*)(class Pythia8::Event &)) &Pythia8::Merging::mergeProcessUMEPS, "C++: Pythia8::Merging::mergeProcessUMEPS(class Pythia8::Event &) --> int", pybind11::arg("process"));
		cl.def("mergeProcessNL3", (int (Pythia8::Merging::*)(class Pythia8::Event &)) &Pythia8::Merging::mergeProcessNL3, "C++: Pythia8::Merging::mergeProcessNL3(class Pythia8::Event &) --> int", pybind11::arg("process"));
		cl.def("mergeProcessUNLOPS", (int (Pythia8::Merging::*)(class Pythia8::Event &)) &Pythia8::Merging::mergeProcessUNLOPS, "C++: Pythia8::Merging::mergeProcessUNLOPS(class Pythia8::Event &) --> int", pybind11::arg("process"));
		cl.def("cutOnProcess", (bool (Pythia8::Merging::*)(class Pythia8::Event &)) &Pythia8::Merging::cutOnProcess, "C++: Pythia8::Merging::cutOnProcess(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("clearInfos", (void (Pythia8::Merging::*)()) &Pythia8::Merging::clearInfos, "C++: Pythia8::Merging::clearInfos() --> void");
		cl.def("clusterAndStore", (int (Pythia8::Merging::*)(class Pythia8::Event &)) &Pythia8::Merging::clusterAndStore, "C++: Pythia8::Merging::clusterAndStore(class Pythia8::Event &) --> int", pybind11::arg("process"));
		cl.def("getDipoles", (void (Pythia8::Merging::*)(int, int, int, const class Pythia8::Event &, class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > &)) &Pythia8::Merging::getDipoles, "C++: Pythia8::Merging::getDipoles(int, int, int, const class Pythia8::Event &, class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > &) --> void", pybind11::arg("iRad"), pybind11::arg("colTag"), pybind11::arg("colSign"), pybind11::arg("event"), pybind11::arg("dipEnds"));
		cl.def("assign", (class Pythia8::Merging & (Pythia8::Merging::*)(const class Pythia8::Merging &)) &Pythia8::Merging::operator=, "C++: Pythia8::Merging::operator=(const class Pythia8::Merging &) --> class Pythia8::Merging &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::RopeDipoleEnd file:Pythia8/Ropewalk.h line:35
		pybind11::class_<Pythia8::RopeDipoleEnd, std::shared_ptr<Pythia8::RopeDipoleEnd>> cl(M("Pythia8"), "RopeDipoleEnd", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::RopeDipoleEnd(); } ) );
		cl.def( pybind11::init<class Pythia8::Event *, int>(), pybind11::arg("eIn"), pybind11::arg("neIn") );

		cl.def( pybind11::init( [](Pythia8::RopeDipoleEnd const &o){ return new Pythia8::RopeDipoleEnd(o); } ) );
		cl.def("getParticlePtr", (class Pythia8::Particle * (Pythia8::RopeDipoleEnd::*)()) &Pythia8::RopeDipoleEnd::getParticlePtr, "C++: Pythia8::RopeDipoleEnd::getParticlePtr() --> class Pythia8::Particle *", pybind11::return_value_policy::automatic);
		cl.def("getNe", (int (Pythia8::RopeDipoleEnd::*)()) &Pythia8::RopeDipoleEnd::getNe, "C++: Pythia8::RopeDipoleEnd::getNe() --> int");
		cl.def("labrap", (double (Pythia8::RopeDipoleEnd::*)()) &Pythia8::RopeDipoleEnd::labrap, "C++: Pythia8::RopeDipoleEnd::labrap() --> double");
		cl.def("rap", (double (Pythia8::RopeDipoleEnd::*)(double)) &Pythia8::RopeDipoleEnd::rap, "C++: Pythia8::RopeDipoleEnd::rap(double) --> double", pybind11::arg("m0"));
		cl.def("rap", (double (Pythia8::RopeDipoleEnd::*)(double, class Pythia8::RotBstMatrix &)) &Pythia8::RopeDipoleEnd::rap, "C++: Pythia8::RopeDipoleEnd::rap(double, class Pythia8::RotBstMatrix &) --> double", pybind11::arg("m0"), pybind11::arg("r"));
		cl.def("assign", (class Pythia8::RopeDipoleEnd & (Pythia8::RopeDipoleEnd::*)(const class Pythia8::RopeDipoleEnd &)) &Pythia8::RopeDipoleEnd::operator=, "C++: Pythia8::RopeDipoleEnd::operator=(const class Pythia8::RopeDipoleEnd &) --> class Pythia8::RopeDipoleEnd &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::OverlappingRopeDipole file:Pythia8/Ropewalk.h line:70
		pybind11::class_<Pythia8::OverlappingRopeDipole, std::shared_ptr<Pythia8::OverlappingRopeDipole>> cl(M("Pythia8"), "OverlappingRopeDipole", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::RopeDipole *, double, class Pythia8::RotBstMatrix &>(), pybind11::arg("d"), pybind11::arg("m0"), pybind11::arg("r") );

		cl.def( pybind11::init( [](Pythia8::OverlappingRopeDipole const &o){ return new Pythia8::OverlappingRopeDipole(o); } ) );
		cl.def_readwrite("dir", &Pythia8::OverlappingRopeDipole::dir);
		cl.def_readwrite("y1", &Pythia8::OverlappingRopeDipole::y1);
		cl.def_readwrite("y2", &Pythia8::OverlappingRopeDipole::y2);
		cl.def_readwrite("b1", &Pythia8::OverlappingRopeDipole::b1);
		cl.def_readwrite("b2", &Pythia8::OverlappingRopeDipole::b2);
		cl.def("overlap", (bool (Pythia8::OverlappingRopeDipole::*)(double, class Pythia8::Vec4, double)) &Pythia8::OverlappingRopeDipole::overlap, "C++: Pythia8::OverlappingRopeDipole::overlap(double, class Pythia8::Vec4, double) --> bool", pybind11::arg("y"), pybind11::arg("ba"), pybind11::arg("r0"));
		cl.def("hadronized", (bool (Pythia8::OverlappingRopeDipole::*)()) &Pythia8::OverlappingRopeDipole::hadronized, "C++: Pythia8::OverlappingRopeDipole::hadronized() --> bool");
		cl.def("assign", (class Pythia8::OverlappingRopeDipole & (Pythia8::OverlappingRopeDipole::*)(const class Pythia8::OverlappingRopeDipole &)) &Pythia8::OverlappingRopeDipole::operator=, "C++: Pythia8::OverlappingRopeDipole::operator=(const class Pythia8::OverlappingRopeDipole &) --> class Pythia8::OverlappingRopeDipole &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::RopeDipole file:Pythia8/Ropewalk.h line:103
		pybind11::class_<Pythia8::RopeDipole, std::shared_ptr<Pythia8::RopeDipole>> cl(M("Pythia8"), "RopeDipole", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](Pythia8::RopeDipole const &o){ return new Pythia8::RopeDipole(o); } ) );
		cl.def("addExcitation", (void (Pythia8::RopeDipole::*)(double, class Pythia8::Particle *)) &Pythia8::RopeDipole::addExcitation, "C++: Pythia8::RopeDipole::addExcitation(double, class Pythia8::Particle *) --> void", pybind11::arg("ylab"), pybind11::arg("ex"));
		cl.def("d1Ptr", (class Pythia8::RopeDipoleEnd * (Pythia8::RopeDipole::*)()) &Pythia8::RopeDipole::d1Ptr, "C++: Pythia8::RopeDipole::d1Ptr() --> class Pythia8::RopeDipoleEnd *", pybind11::return_value_policy::automatic);
		cl.def("d2Ptr", (class Pythia8::RopeDipoleEnd * (Pythia8::RopeDipole::*)()) &Pythia8::RopeDipole::d2Ptr, "C++: Pythia8::RopeDipole::d2Ptr() --> class Pythia8::RopeDipoleEnd *", pybind11::return_value_policy::automatic);
		cl.def("getDipoleRestFrame", (class Pythia8::RotBstMatrix (Pythia8::RopeDipole::*)()) &Pythia8::RopeDipole::getDipoleRestFrame, "C++: Pythia8::RopeDipole::getDipoleRestFrame() --> class Pythia8::RotBstMatrix");
		cl.def("getDipoleLabFrame", (class Pythia8::RotBstMatrix (Pythia8::RopeDipole::*)()) &Pythia8::RopeDipole::getDipoleLabFrame, "C++: Pythia8::RopeDipole::getDipoleLabFrame() --> class Pythia8::RotBstMatrix");
		cl.def("dipoleMomentum", (class Pythia8::Vec4 (Pythia8::RopeDipole::*)()) &Pythia8::RopeDipole::dipoleMomentum, "C++: Pythia8::RopeDipole::dipoleMomentum() --> class Pythia8::Vec4");
		cl.def("bInterpolateDip", (class Pythia8::Vec4 (Pythia8::RopeDipole::*)(double, double)) &Pythia8::RopeDipole::bInterpolateDip, "C++: Pythia8::RopeDipole::bInterpolateDip(double, double) --> class Pythia8::Vec4", pybind11::arg("y"), pybind11::arg("m0"));
		cl.def("bInterpolateLab", (class Pythia8::Vec4 (Pythia8::RopeDipole::*)(double, double)) &Pythia8::RopeDipole::bInterpolateLab, "C++: Pythia8::RopeDipole::bInterpolateLab(double, double) --> class Pythia8::Vec4", pybind11::arg("y"), pybind11::arg("m0"));
		cl.def("bInterpolate", (class Pythia8::Vec4 (Pythia8::RopeDipole::*)(double, class Pythia8::RotBstMatrix, double)) &Pythia8::RopeDipole::bInterpolate, "C++: Pythia8::RopeDipole::bInterpolate(double, class Pythia8::RotBstMatrix, double) --> class Pythia8::Vec4", pybind11::arg("y"), pybind11::arg("rb"), pybind11::arg("m0"));
		cl.def("getOverlaps", (struct std::pair<int, int> (Pythia8::RopeDipole::*)(double, double, double)) &Pythia8::RopeDipole::getOverlaps, "C++: Pythia8::RopeDipole::getOverlaps(double, double, double) --> struct std::pair<int, int>", pybind11::arg("yfrac"), pybind11::arg("m0"), pybind11::arg("r0"));
		cl.def("addOverlappingDipole", (void (Pythia8::RopeDipole::*)(class Pythia8::OverlappingRopeDipole &)) &Pythia8::RopeDipole::addOverlappingDipole, "C++: Pythia8::RopeDipole::addOverlappingDipole(class Pythia8::OverlappingRopeDipole &) --> void", pybind11::arg("d"));
		cl.def("maxRapidity", (double (Pythia8::RopeDipole::*)(double)) &Pythia8::RopeDipole::maxRapidity, "C++: Pythia8::RopeDipole::maxRapidity(double) --> double", pybind11::arg("m0"));
		cl.def("minRapidity", (double (Pythia8::RopeDipole::*)(double)) &Pythia8::RopeDipole::minRapidity, "C++: Pythia8::RopeDipole::minRapidity(double) --> double", pybind11::arg("m0"));
		cl.def("maxRapidity", (double (Pythia8::RopeDipole::*)(double, class Pythia8::RotBstMatrix &)) &Pythia8::RopeDipole::maxRapidity, "C++: Pythia8::RopeDipole::maxRapidity(double, class Pythia8::RotBstMatrix &) --> double", pybind11::arg("m0"), pybind11::arg("r"));
		cl.def("minRapidity", (double (Pythia8::RopeDipole::*)(double, class Pythia8::RotBstMatrix &)) &Pythia8::RopeDipole::minRapidity, "C++: Pythia8::RopeDipole::minRapidity(double, class Pythia8::RotBstMatrix &) --> double", pybind11::arg("m0"), pybind11::arg("r"));
		cl.def("propagateInit", (void (Pythia8::RopeDipole::*)(double)) &Pythia8::RopeDipole::propagateInit, "C++: Pythia8::RopeDipole::propagateInit(double) --> void", pybind11::arg("deltat"));
		cl.def("propagate", (void (Pythia8::RopeDipole::*)(double, double)) &Pythia8::RopeDipole::propagate, "C++: Pythia8::RopeDipole::propagate(double, double) --> void", pybind11::arg("deltat"), pybind11::arg("m0"));
		cl.def("splitMomentum", [](Pythia8::RopeDipole &o, class Pythia8::Vec4 const & a0, class Pythia8::Particle * a1, class Pythia8::Particle * a2) -> void { return o.splitMomentum(a0, a1, a2); }, "", pybind11::arg("mom"), pybind11::arg("p1"), pybind11::arg("p2"));
		cl.def("splitMomentum", (void (Pythia8::RopeDipole::*)(class Pythia8::Vec4, class Pythia8::Particle *, class Pythia8::Particle *, double)) &Pythia8::RopeDipole::splitMomentum, "C++: Pythia8::RopeDipole::splitMomentum(class Pythia8::Vec4, class Pythia8::Particle *, class Pythia8::Particle *, double) --> void", pybind11::arg("mom"), pybind11::arg("p1"), pybind11::arg("p2"), pybind11::arg("frac"));
		cl.def("excitationsToString", (void (Pythia8::RopeDipole::*)(double, class Pythia8::Event &)) &Pythia8::RopeDipole::excitationsToString, "C++: Pythia8::RopeDipole::excitationsToString(double, class Pythia8::Event &) --> void", pybind11::arg("m0"), pybind11::arg("event"));
		cl.def("hadronized", (bool (Pythia8::RopeDipole::*)()) &Pythia8::RopeDipole::hadronized, "C++: Pythia8::RopeDipole::hadronized() --> bool");
		cl.def("index", (int (Pythia8::RopeDipole::*)()) &Pythia8::RopeDipole::index, "C++: Pythia8::RopeDipole::index() --> int");
		cl.def("recoil", [](Pythia8::RopeDipole &o, class Pythia8::Vec4 & a0) -> bool { return o.recoil(a0); }, "", pybind11::arg("pg"));
		cl.def("recoil", (bool (Pythia8::RopeDipole::*)(class Pythia8::Vec4 &, bool)) &Pythia8::RopeDipole::recoil, "C++: Pythia8::RopeDipole::recoil(class Pythia8::Vec4 &, bool) --> bool", pybind11::arg("pg"), pybind11::arg("dummy"));
		cl.def("hadronized", (void (Pythia8::RopeDipole::*)(bool)) &Pythia8::RopeDipole::hadronized, "C++: Pythia8::RopeDipole::hadronized(bool) --> void", pybind11::arg("h"));
		cl.def("nExcitations", (int (Pythia8::RopeDipole::*)()) &Pythia8::RopeDipole::nExcitations, "C++: Pythia8::RopeDipole::nExcitations() --> int");
	}
	{ // Pythia8::Ropewalk file:Pythia8/Ropewalk.h line:211
		pybind11::class_<Pythia8::Ropewalk, std::shared_ptr<Pythia8::Ropewalk>, PyCallBack_Pythia8_Ropewalk, Pythia8::StringInteractions> cl(M("Pythia8"), "Ropewalk", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Ropewalk(); }, [](){ return new PyCallBack_Pythia8_Ropewalk(); } ) );
		cl.def("init", (bool (Pythia8::Ropewalk::*)()) &Pythia8::Ropewalk::init, "C++: Pythia8::Ropewalk::init() --> bool");
		cl.def("extractDipoles", (bool (Pythia8::Ropewalk::*)(class Pythia8::Event &, class Pythia8::ColConfig &)) &Pythia8::Ropewalk::extractDipoles, "C++: Pythia8::Ropewalk::extractDipoles(class Pythia8::Event &, class Pythia8::ColConfig &) --> bool", pybind11::arg("event"), pybind11::arg("colConfig"));
		cl.def("calculateOverlaps", (bool (Pythia8::Ropewalk::*)()) &Pythia8::Ropewalk::calculateOverlaps, "C++: Pythia8::Ropewalk::calculateOverlaps() --> bool");
		cl.def("getKappaHere", (double (Pythia8::Ropewalk::*)(int, int, double)) &Pythia8::Ropewalk::getKappaHere, "C++: Pythia8::Ropewalk::getKappaHere(int, int, double) --> double", pybind11::arg("e1"), pybind11::arg("e2"), pybind11::arg("yfrac"));
		cl.def("multiplicity", (double (Pythia8::Ropewalk::*)(double, double)) &Pythia8::Ropewalk::multiplicity, "C++: Pythia8::Ropewalk::multiplicity(double, double) --> double", pybind11::arg("p"), pybind11::arg("q"));
		cl.def("averageKappa", (double (Pythia8::Ropewalk::*)()) &Pythia8::Ropewalk::averageKappa, "C++: Pythia8::Ropewalk::averageKappa() --> double");
		cl.def("select", (struct std::pair<int, int> (Pythia8::Ropewalk::*)(int, int, class Pythia8::Rndm *)) &Pythia8::Ropewalk::select, "C++: Pythia8::Ropewalk::select(int, int, class Pythia8::Rndm *) --> struct std::pair<int, int>", pybind11::arg("m"), pybind11::arg("n"), pybind11::arg("rndm"));
		cl.def("shoveTheDipoles", (void (Pythia8::Ropewalk::*)(class Pythia8::Event &)) &Pythia8::Ropewalk::shoveTheDipoles, "C++: Pythia8::Ropewalk::shoveTheDipoles(class Pythia8::Event &) --> void", pybind11::arg("event"));
	}
}
