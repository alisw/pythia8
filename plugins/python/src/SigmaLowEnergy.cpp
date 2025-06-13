#include <Pythia8/Basics.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/Info.h>
#include <Pythia8/LowEnergyProcess.h>
#include <Pythia8/MiniStringFragmentation.h>
#include <Pythia8/NucleonExcitations.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/StringFragmentation.h>
#include <Pythia8/StringInteractions.h>
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

// Pythia8::SigmaCombined file:Pythia8/SigmaLowEnergy.h line:135
struct PyCallBack_Pythia8_SigmaCombined : public Pythia8::SigmaCombined {
	using Pythia8::SigmaCombined::SigmaCombined;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaCombined *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaCombined *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaCombined *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SigmaCombined *>(this), "onStat");
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

// Pythia8::LowEnergyProcess file:Pythia8/LowEnergyProcess.h line:28
struct PyCallBack_Pythia8_LowEnergyProcess : public Pythia8::LowEnergyProcess {
	using Pythia8::LowEnergyProcess::LowEnergyProcess;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LowEnergyProcess *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LowEnergyProcess *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LowEnergyProcess *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LowEnergyProcess *>(this), "onStat");
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

// Pythia8::PartonVertex file:Pythia8/PartonVertex.h line:24
struct PyCallBack_Pythia8_PartonVertex : public Pythia8::PartonVertex {
	using Pythia8::PartonVertex::PartonVertex;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartonVertex::init();
	}
	void vertexBeam(int a0, class std::vector<int, class std::allocator<int> > & a1, class std::vector<int, class std::allocator<int> > & a2, class Pythia8::Event & a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "vertexBeam");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartonVertex::vertexBeam(a0, a1, a2, a3);
	}
	void vertexMPI(int a0, int a1, double a2, class Pythia8::Event & a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "vertexMPI");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartonVertex::vertexMPI(a0, a1, a2, a3);
	}
	void vertexFSR(int a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "vertexFSR");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartonVertex::vertexFSR(a0, a1);
	}
	void vertexISR(int a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "vertexISR");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartonVertex::vertexISR(a0, a1);
	}
	void vertexHadrons(int a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "vertexHadrons");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartonVertex::vertexHadrons(a0, a1);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PartonVertex *>(this), "onStat");
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

void bind_Pythia8_SigmaLowEnergy(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::SigmaCombined file:Pythia8/SigmaLowEnergy.h line:135
		pybind11::class_<Pythia8::SigmaCombined, std::shared_ptr<Pythia8::SigmaCombined>, PyCallBack_Pythia8_SigmaCombined, Pythia8::PhysicsBase> cl(M("Pythia8"), "SigmaCombined", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SigmaCombined(); }, [](){ return new PyCallBack_Pythia8_SigmaCombined(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_SigmaCombined const &o){ return new PyCallBack_Pythia8_SigmaCombined(o); } ) );
		cl.def( pybind11::init( [](Pythia8::SigmaCombined const &o){ return new Pythia8::SigmaCombined(o); } ) );
		cl.def("init", (void (Pythia8::SigmaCombined::*)(class Pythia8::SigmaLowEnergy *)) &Pythia8::SigmaCombined::init, "C++: Pythia8::SigmaCombined::init(class Pythia8::SigmaLowEnergy *) --> void", pybind11::arg("sigmaLowPtrIn"));
		cl.def("sigmaTotal", [](Pythia8::SigmaCombined &o, int const & a0, int const & a1, double const & a2, double const & a3, double const & a4) -> double { return o.sigmaTotal(a0, a1, a2, a3, a4); }, "", pybind11::arg("id1"), pybind11::arg("id2"), pybind11::arg("eCM12"), pybind11::arg("m1"), pybind11::arg("m2"));
		cl.def("sigmaTotal", (double (Pythia8::SigmaCombined::*)(int, int, double, double, double, int)) &Pythia8::SigmaCombined::sigmaTotal, "C++: Pythia8::SigmaCombined::sigmaTotal(int, int, double, double, double, int) --> double", pybind11::arg("id1"), pybind11::arg("id2"), pybind11::arg("eCM12"), pybind11::arg("m1"), pybind11::arg("m2"), pybind11::arg("mixLoHi"));
		cl.def("sigmaPartial", [](Pythia8::SigmaCombined &o, int const & a0, int const & a1, double const & a2, double const & a3, double const & a4, int const & a5) -> double { return o.sigmaPartial(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("id1"), pybind11::arg("id2"), pybind11::arg("eCM12"), pybind11::arg("m1"), pybind11::arg("m2"), pybind11::arg("type"));
		cl.def("sigmaPartial", (double (Pythia8::SigmaCombined::*)(int, int, double, double, double, int, int)) &Pythia8::SigmaCombined::sigmaPartial, "C++: Pythia8::SigmaCombined::sigmaPartial(int, int, double, double, double, int, int) --> double", pybind11::arg("id1"), pybind11::arg("id2"), pybind11::arg("eCM12"), pybind11::arg("m1"), pybind11::arg("m2"), pybind11::arg("type"), pybind11::arg("mixLoHi"));
		cl.def("assign", (class Pythia8::SigmaCombined & (Pythia8::SigmaCombined::*)(const class Pythia8::SigmaCombined &)) &Pythia8::SigmaCombined::operator=, "C++: Pythia8::SigmaCombined::operator=(const class Pythia8::SigmaCombined &) --> class Pythia8::SigmaCombined &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LowEnergyProcess file:Pythia8/LowEnergyProcess.h line:28
		pybind11::class_<Pythia8::LowEnergyProcess, std::shared_ptr<Pythia8::LowEnergyProcess>, PyCallBack_Pythia8_LowEnergyProcess, Pythia8::PhysicsBase> cl(M("Pythia8"), "LowEnergyProcess", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LowEnergyProcess(); }, [](){ return new PyCallBack_Pythia8_LowEnergyProcess(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_LowEnergyProcess const &o){ return new PyCallBack_Pythia8_LowEnergyProcess(o); } ) );
		cl.def( pybind11::init( [](Pythia8::LowEnergyProcess const &o){ return new Pythia8::LowEnergyProcess(o); } ) );
		cl.def_readwrite("leEvent", &Pythia8::LowEnergyProcess::leEvent);
		cl.def("init", (void (Pythia8::LowEnergyProcess::*)(class Pythia8::StringFlav *, class Pythia8::StringFragmentation *, class Pythia8::MiniStringFragmentation *, class Pythia8::SigmaLowEnergy *, class Pythia8::NucleonExcitations *)) &Pythia8::LowEnergyProcess::init, "C++: Pythia8::LowEnergyProcess::init(class Pythia8::StringFlav *, class Pythia8::StringFragmentation *, class Pythia8::MiniStringFragmentation *, class Pythia8::SigmaLowEnergy *, class Pythia8::NucleonExcitations *) --> void", pybind11::arg("flavSelPtrIn"), pybind11::arg("stringFragPtrIn"), pybind11::arg("ministringFragPtrIn"), pybind11::arg("sigmaLowEnergyPtrIn"), pybind11::arg("nucleonExcitationsPtrIn"));
		cl.def("collide", [](Pythia8::LowEnergyProcess &o, int const & a0, int const & a1, int const & a2, class Pythia8::Event & a3) -> bool { return o.collide(a0, a1, a2, a3); }, "", pybind11::arg("i1"), pybind11::arg("i2"), pybind11::arg("typeIn"), pybind11::arg("event"));
		cl.def("collide", [](Pythia8::LowEnergyProcess &o, int const & a0, int const & a1, int const & a2, class Pythia8::Event & a3, class Pythia8::Vec4 const & a4) -> bool { return o.collide(a0, a1, a2, a3, a4); }, "", pybind11::arg("i1"), pybind11::arg("i2"), pybind11::arg("typeIn"), pybind11::arg("event"), pybind11::arg("vtx"));
		cl.def("collide", [](Pythia8::LowEnergyProcess &o, int const & a0, int const & a1, int const & a2, class Pythia8::Event & a3, class Pythia8::Vec4 const & a4, class Pythia8::Vec4 const & a5) -> bool { return o.collide(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("i1"), pybind11::arg("i2"), pybind11::arg("typeIn"), pybind11::arg("event"), pybind11::arg("vtx"), pybind11::arg("vtx1"));
		cl.def("collide", (bool (Pythia8::LowEnergyProcess::*)(int, int, int, class Pythia8::Event &, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4)) &Pythia8::LowEnergyProcess::collide, "C++: Pythia8::LowEnergyProcess::collide(int, int, int, class Pythia8::Event &, class Pythia8::Vec4, class Pythia8::Vec4, class Pythia8::Vec4) --> bool", pybind11::arg("i1"), pybind11::arg("i2"), pybind11::arg("typeIn"), pybind11::arg("event"), pybind11::arg("vtx"), pybind11::arg("vtx1"), pybind11::arg("vtx2"));
		cl.def("bSlope", [](Pythia8::LowEnergyProcess &o, int const & a0, int const & a1, double const & a2, double const & a3, double const & a4) -> double { return o.bSlope(a0, a1, a2, a3, a4); }, "", pybind11::arg("id1In"), pybind11::arg("id2In"), pybind11::arg("eCMIn"), pybind11::arg("mAIn"), pybind11::arg("mBIn"));
		cl.def("bSlope", (double (Pythia8::LowEnergyProcess::*)(int, int, double, double, double, int)) &Pythia8::LowEnergyProcess::bSlope, "C++: Pythia8::LowEnergyProcess::bSlope(int, int, double, double, double, int) --> double", pybind11::arg("id1In"), pybind11::arg("id2In"), pybind11::arg("eCMIn"), pybind11::arg("mAIn"), pybind11::arg("mBIn"), pybind11::arg("typeIn"));
		cl.def("assign", (class Pythia8::LowEnergyProcess & (Pythia8::LowEnergyProcess::*)(const class Pythia8::LowEnergyProcess &)) &Pythia8::LowEnergyProcess::operator=, "C++: Pythia8::LowEnergyProcess::operator=(const class Pythia8::LowEnergyProcess &) --> class Pythia8::LowEnergyProcess &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PartonSystem file:Pythia8/PartonSystems.h line:22
		pybind11::class_<Pythia8::PartonSystem, std::shared_ptr<Pythia8::PartonSystem>> cl(M("Pythia8"), "PartonSystem", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PartonSystem(); } ) );
		cl.def( pybind11::init( [](Pythia8::PartonSystem const &o){ return new Pythia8::PartonSystem(o); } ) );
		cl.def_readwrite("hard", &Pythia8::PartonSystem::hard);
		cl.def_readwrite("iInA", &Pythia8::PartonSystem::iInA);
		cl.def_readwrite("iInB", &Pythia8::PartonSystem::iInB);
		cl.def_readwrite("iInRes", &Pythia8::PartonSystem::iInRes);
		cl.def_readwrite("iOut", &Pythia8::PartonSystem::iOut);
		cl.def_readwrite("sHat", &Pythia8::PartonSystem::sHat);
		cl.def_readwrite("pTHat", &Pythia8::PartonSystem::pTHat);
		cl.def("assign", (class Pythia8::PartonSystem & (Pythia8::PartonSystem::*)(const class Pythia8::PartonSystem &)) &Pythia8::PartonSystem::operator=, "C++: Pythia8::PartonSystem::operator=(const class Pythia8::PartonSystem &) --> class Pythia8::PartonSystem &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::PartonSystems file:Pythia8/PartonSystems.h line:42
		pybind11::class_<Pythia8::PartonSystems, std::shared_ptr<Pythia8::PartonSystems>> cl(M("Pythia8"), "PartonSystems", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PartonSystems(); } ) );
		cl.def( pybind11::init( [](Pythia8::PartonSystems const &o){ return new Pythia8::PartonSystems(o); } ) );
		cl.def("clear", (void (Pythia8::PartonSystems::*)()) &Pythia8::PartonSystems::clear, "C++: Pythia8::PartonSystems::clear() --> void");
		cl.def("addSys", (int (Pythia8::PartonSystems::*)()) &Pythia8::PartonSystems::addSys, "C++: Pythia8::PartonSystems::addSys() --> int");
		cl.def("sizeSys", (int (Pythia8::PartonSystems::*)() const) &Pythia8::PartonSystems::sizeSys, "C++: Pythia8::PartonSystems::sizeSys() const --> int");
		cl.def("setHard", (void (Pythia8::PartonSystems::*)(int, bool)) &Pythia8::PartonSystems::setHard, "C++: Pythia8::PartonSystems::setHard(int, bool) --> void", pybind11::arg("iSys"), pybind11::arg("hard"));
		cl.def("setInA", (void (Pythia8::PartonSystems::*)(int, int)) &Pythia8::PartonSystems::setInA, "C++: Pythia8::PartonSystems::setInA(int, int) --> void", pybind11::arg("iSys"), pybind11::arg("iPos"));
		cl.def("setInB", (void (Pythia8::PartonSystems::*)(int, int)) &Pythia8::PartonSystems::setInB, "C++: Pythia8::PartonSystems::setInB(int, int) --> void", pybind11::arg("iSys"), pybind11::arg("iPos"));
		cl.def("setInRes", (void (Pythia8::PartonSystems::*)(int, int)) &Pythia8::PartonSystems::setInRes, "C++: Pythia8::PartonSystems::setInRes(int, int) --> void", pybind11::arg("iSys"), pybind11::arg("iPos"));
		cl.def("addOut", (void (Pythia8::PartonSystems::*)(int, int)) &Pythia8::PartonSystems::addOut, "C++: Pythia8::PartonSystems::addOut(int, int) --> void", pybind11::arg("iSys"), pybind11::arg("iPos"));
		cl.def("popBackOut", (void (Pythia8::PartonSystems::*)(int)) &Pythia8::PartonSystems::popBackOut, "C++: Pythia8::PartonSystems::popBackOut(int) --> void", pybind11::arg("iSys"));
		cl.def("setOut", (void (Pythia8::PartonSystems::*)(int, int, int)) &Pythia8::PartonSystems::setOut, "C++: Pythia8::PartonSystems::setOut(int, int, int) --> void", pybind11::arg("iSys"), pybind11::arg("iMem"), pybind11::arg("iPos"));
		cl.def("replace", (void (Pythia8::PartonSystems::*)(int, int, int)) &Pythia8::PartonSystems::replace, "C++: Pythia8::PartonSystems::replace(int, int, int) --> void", pybind11::arg("iSys"), pybind11::arg("iPosOld"), pybind11::arg("iPosNew"));
		cl.def("setSHat", (void (Pythia8::PartonSystems::*)(int, double)) &Pythia8::PartonSystems::setSHat, "C++: Pythia8::PartonSystems::setSHat(int, double) --> void", pybind11::arg("iSys"), pybind11::arg("sHatIn"));
		cl.def("setPTHat", (void (Pythia8::PartonSystems::*)(int, double)) &Pythia8::PartonSystems::setPTHat, "C++: Pythia8::PartonSystems::setPTHat(int, double) --> void", pybind11::arg("iSys"), pybind11::arg("pTHatIn"));
		cl.def("setSizeSys", (void (Pythia8::PartonSystems::*)(int)) &Pythia8::PartonSystems::setSizeSys, "C++: Pythia8::PartonSystems::setSizeSys(int) --> void", pybind11::arg("iSize"));
		cl.def("hasInAB", (bool (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::hasInAB, "C++: Pythia8::PartonSystems::hasInAB(int) const --> bool", pybind11::arg("iSys"));
		cl.def("hasInRes", (bool (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::hasInRes, "C++: Pythia8::PartonSystems::hasInRes(int) const --> bool", pybind11::arg("iSys"));
		cl.def("getHard", (bool (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::getHard, "C++: Pythia8::PartonSystems::getHard(int) const --> bool", pybind11::arg("iSys"));
		cl.def("getInA", (int (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::getInA, "C++: Pythia8::PartonSystems::getInA(int) const --> int", pybind11::arg("iSys"));
		cl.def("getInB", (int (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::getInB, "C++: Pythia8::PartonSystems::getInB(int) const --> int", pybind11::arg("iSys"));
		cl.def("getInRes", (int (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::getInRes, "C++: Pythia8::PartonSystems::getInRes(int) const --> int", pybind11::arg("iSys"));
		cl.def("sizeOut", (int (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::sizeOut, "C++: Pythia8::PartonSystems::sizeOut(int) const --> int", pybind11::arg("iSys"));
		cl.def("getOut", (int (Pythia8::PartonSystems::*)(int, int) const) &Pythia8::PartonSystems::getOut, "C++: Pythia8::PartonSystems::getOut(int, int) const --> int", pybind11::arg("iSys"), pybind11::arg("iMem"));
		cl.def("sizeAll", (int (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::sizeAll, "C++: Pythia8::PartonSystems::sizeAll(int) const --> int", pybind11::arg("iSys"));
		cl.def("getAll", (int (Pythia8::PartonSystems::*)(int, int) const) &Pythia8::PartonSystems::getAll, "C++: Pythia8::PartonSystems::getAll(int, int) const --> int", pybind11::arg("iSys"), pybind11::arg("iMem"));
		cl.def("getSHat", (double (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::getSHat, "C++: Pythia8::PartonSystems::getSHat(int) const --> double", pybind11::arg("iSys"));
		cl.def("getPTHat", (double (Pythia8::PartonSystems::*)(int) const) &Pythia8::PartonSystems::getPTHat, "C++: Pythia8::PartonSystems::getPTHat(int) const --> double", pybind11::arg("iSys"));
		cl.def("getSystemOf", [](Pythia8::PartonSystems const &o, int const & a0) -> int { return o.getSystemOf(a0); }, "", pybind11::arg("iPos"));
		cl.def("getSystemOf", (int (Pythia8::PartonSystems::*)(int, bool) const) &Pythia8::PartonSystems::getSystemOf, "C++: Pythia8::PartonSystems::getSystemOf(int, bool) const --> int", pybind11::arg("iPos"), pybind11::arg("alsoIn"));
		cl.def("getIndexOfOut", (int (Pythia8::PartonSystems::*)(int, int) const) &Pythia8::PartonSystems::getIndexOfOut, "C++: Pythia8::PartonSystems::getIndexOfOut(int, int) const --> int", pybind11::arg("iSys"), pybind11::arg("iPos"));
		cl.def("list", (void (Pythia8::PartonSystems::*)() const) &Pythia8::PartonSystems::list, "C++: Pythia8::PartonSystems::list() const --> void");
		cl.def("popBack", (void (Pythia8::PartonSystems::*)()) &Pythia8::PartonSystems::popBack, "C++: Pythia8::PartonSystems::popBack() --> void");
	}
	{ // Pythia8::PartonVertex file:Pythia8/PartonVertex.h line:24
		pybind11::class_<Pythia8::PartonVertex, std::shared_ptr<Pythia8::PartonVertex>, PyCallBack_Pythia8_PartonVertex, Pythia8::PhysicsBase> cl(M("Pythia8"), "PartonVertex", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PartonVertex(); }, [](){ return new PyCallBack_Pythia8_PartonVertex(); } ) );
		cl.def("init", (void (Pythia8::PartonVertex::*)()) &Pythia8::PartonVertex::init, "C++: Pythia8::PartonVertex::init() --> void");
		cl.def("vertexBeam", (void (Pythia8::PartonVertex::*)(int, class std::vector<int, class std::allocator<int> > &, class std::vector<int, class std::allocator<int> > &, class Pythia8::Event &)) &Pythia8::PartonVertex::vertexBeam, "C++: Pythia8::PartonVertex::vertexBeam(int, class std::vector<int, class std::allocator<int> > &, class std::vector<int, class std::allocator<int> > &, class Pythia8::Event &) --> void", pybind11::arg("iBeam"), pybind11::arg("iRemn"), pybind11::arg("iInit"), pybind11::arg("event"));
		cl.def("vertexMPI", (void (Pythia8::PartonVertex::*)(int, int, double, class Pythia8::Event &)) &Pythia8::PartonVertex::vertexMPI, "C++: Pythia8::PartonVertex::vertexMPI(int, int, double, class Pythia8::Event &) --> void", pybind11::arg("iBeg"), pybind11::arg("nAdd"), pybind11::arg("bNowIn"), pybind11::arg("event"));
		cl.def("vertexFSR", (void (Pythia8::PartonVertex::*)(int, class Pythia8::Event &)) &Pythia8::PartonVertex::vertexFSR, "C++: Pythia8::PartonVertex::vertexFSR(int, class Pythia8::Event &) --> void", pybind11::arg("iNow"), pybind11::arg("event"));
		cl.def("vertexISR", (void (Pythia8::PartonVertex::*)(int, class Pythia8::Event &)) &Pythia8::PartonVertex::vertexISR, "C++: Pythia8::PartonVertex::vertexISR(int, class Pythia8::Event &) --> void", pybind11::arg("iNow"), pybind11::arg("event"));
		cl.def("vertexHadrons", (void (Pythia8::PartonVertex::*)(int, class Pythia8::Event &)) &Pythia8::PartonVertex::vertexHadrons, "C++: Pythia8::PartonVertex::vertexHadrons(int, class Pythia8::Event &) --> void", pybind11::arg("nBefFrag"), pybind11::arg("event"));
		cl.def("assign", (class Pythia8::PartonVertex & (Pythia8::PartonVertex::*)(const class Pythia8::PartonVertex &)) &Pythia8::PartonVertex::operator=, "C++: Pythia8::PartonVertex::operator=(const class Pythia8::PartonVertex &) --> class Pythia8::PartonVertex &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
