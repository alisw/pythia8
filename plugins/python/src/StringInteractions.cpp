#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamRemnants.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/ColourReconnection.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/Weights.h>
#include <functional>
#include <istream>
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

// Pythia8::StringRepulsionBase file:Pythia8/StringInteractions.h line:156
struct PyCallBack_Pythia8_StringRepulsionBase : public Pythia8::StringRepulsionBase {
	using Pythia8::StringRepulsionBase::StringRepulsionBase;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringRepulsionBase *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringRepulsionBase::init();
	}
	bool stringRepulsion(class Pythia8::Event & a0, class Pythia8::ColConfig & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringRepulsionBase *>(this), "stringRepulsion");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"StringRepulsionBase::stringRepulsion\"");
	}
	bool hadronRepulsion(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringRepulsionBase *>(this), "hadronRepulsion");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringRepulsionBase::hadronRepulsion(a0);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringRepulsionBase *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringRepulsionBase *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringRepulsionBase *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringRepulsionBase *>(this), "onStat");
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

// Pythia8::FragmentationModifierBase file:Pythia8/StringInteractions.h line:182
struct PyCallBack_Pythia8_FragmentationModifierBase : public Pythia8::FragmentationModifierBase {
	using Pythia8::FragmentationModifierBase::FragmentationModifierBase;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FragmentationModifierBase *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return FragmentationModifierBase::init();
	}
	bool initEvent(class Pythia8::Event & a0, class Pythia8::ColConfig & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FragmentationModifierBase *>(this), "initEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"FragmentationModifierBase::initEvent\"");
	}
	bool doChangeFragPar(class Pythia8::StringFlav * a0, class Pythia8::StringZ * a1, class Pythia8::StringPT * a2, double a3, class std::vector<int, class std::allocator<int> > a4, int a5) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FragmentationModifierBase *>(this), "doChangeFragPar");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"FragmentationModifierBase::doChangeFragPar\"");
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FragmentationModifierBase *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FragmentationModifierBase *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FragmentationModifierBase *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::FragmentationModifierBase *>(this), "onStat");
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

// Pythia8::ColourParticle file:Pythia8/ColourReconnection.h line:142
struct PyCallBack_Pythia8_ColourParticle : public Pythia8::ColourParticle {
	using Pythia8::ColourParticle::ColourParticle;

	int index() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourParticle *>(this), "index");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return Particle::index();
	}
};

// Pythia8::ColourReconnection file:Pythia8/ColourReconnection.h line:167
struct PyCallBack_Pythia8_ColourReconnection : public Pythia8::ColourReconnection {
	using Pythia8::ColourReconnection::ColourReconnection;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnection *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ColourReconnection::init();
	}
	void reassignBeamPtrs(class Pythia8::BeamParticle * a0, class Pythia8::BeamParticle * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnection *>(this), "reassignBeamPtrs");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColourReconnection::reassignBeamPtrs(a0, a1);
	}
	bool next(class Pythia8::Event & a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnection *>(this), "next");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ColourReconnection::next(a0, a1);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnection *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnection *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnection *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ColourReconnection *>(this), "onStat");
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

// Pythia8::BeamRemnants file:Pythia8/BeamRemnants.h line:35
struct PyCallBack_Pythia8_BeamRemnants : public Pythia8::BeamRemnants {
	using Pythia8::BeamRemnants::BeamRemnants;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamRemnants *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BeamRemnants::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamRemnants *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamRemnants *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BeamRemnants *>(this), "onStat");
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

void bind_Pythia8_StringInteractions(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::StringRepulsionBase file:Pythia8/StringInteractions.h line:156
		pybind11::class_<Pythia8::StringRepulsionBase, std::shared_ptr<Pythia8::StringRepulsionBase>, PyCallBack_Pythia8_StringRepulsionBase, Pythia8::PhysicsBase> cl(M("Pythia8"), "StringRepulsionBase", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_StringRepulsionBase(); } ) );
		cl.def(pybind11::init<PyCallBack_Pythia8_StringRepulsionBase const &>());
		cl.def("init", (bool (Pythia8::StringRepulsionBase::*)()) &Pythia8::StringRepulsionBase::init, "C++: Pythia8::StringRepulsionBase::init() --> bool");
		cl.def("stringRepulsion", (bool (Pythia8::StringRepulsionBase::*)(class Pythia8::Event &, class Pythia8::ColConfig &)) &Pythia8::StringRepulsionBase::stringRepulsion, "C++: Pythia8::StringRepulsionBase::stringRepulsion(class Pythia8::Event &, class Pythia8::ColConfig &) --> bool", pybind11::arg("event"), pybind11::arg("colConfig"));
		cl.def("hadronRepulsion", (bool (Pythia8::StringRepulsionBase::*)(class Pythia8::Event &)) &Pythia8::StringRepulsionBase::hadronRepulsion, "C++: Pythia8::StringRepulsionBase::hadronRepulsion(class Pythia8::Event &) --> bool", pybind11::arg(""));
		cl.def("assign", (class Pythia8::StringRepulsionBase & (Pythia8::StringRepulsionBase::*)(const class Pythia8::StringRepulsionBase &)) &Pythia8::StringRepulsionBase::operator=, "C++: Pythia8::StringRepulsionBase::operator=(const class Pythia8::StringRepulsionBase &) --> class Pythia8::StringRepulsionBase &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::FragmentationModifierBase file:Pythia8/StringInteractions.h line:182
		pybind11::class_<Pythia8::FragmentationModifierBase, std::shared_ptr<Pythia8::FragmentationModifierBase>, PyCallBack_Pythia8_FragmentationModifierBase, Pythia8::PhysicsBase> cl(M("Pythia8"), "FragmentationModifierBase", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_FragmentationModifierBase(); } ) );
		cl.def(pybind11::init<PyCallBack_Pythia8_FragmentationModifierBase const &>());
		cl.def("init", (bool (Pythia8::FragmentationModifierBase::*)()) &Pythia8::FragmentationModifierBase::init, "C++: Pythia8::FragmentationModifierBase::init() --> bool");
		cl.def("initEvent", (bool (Pythia8::FragmentationModifierBase::*)(class Pythia8::Event &, class Pythia8::ColConfig &)) &Pythia8::FragmentationModifierBase::initEvent, "C++: Pythia8::FragmentationModifierBase::initEvent(class Pythia8::Event &, class Pythia8::ColConfig &) --> bool", pybind11::arg("event"), pybind11::arg("colConfig"));
		cl.def("doChangeFragPar", (bool (Pythia8::FragmentationModifierBase::*)(class Pythia8::StringFlav *, class Pythia8::StringZ *, class Pythia8::StringPT *, double, class std::vector<int, class std::allocator<int> >, int)) &Pythia8::FragmentationModifierBase::doChangeFragPar, "C++: Pythia8::FragmentationModifierBase::doChangeFragPar(class Pythia8::StringFlav *, class Pythia8::StringZ *, class Pythia8::StringPT *, double, class std::vector<int, class std::allocator<int> >, int) --> bool", pybind11::arg("flavPtr"), pybind11::arg("zPtr"), pybind11::arg("pTPtr"), pybind11::arg("m2Had"), pybind11::arg("iParton"), pybind11::arg("endId"));
		cl.def("assign", (class Pythia8::FragmentationModifierBase & (Pythia8::FragmentationModifierBase::*)(const class Pythia8::FragmentationModifierBase &)) &Pythia8::FragmentationModifierBase::operator=, "C++: Pythia8::FragmentationModifierBase::operator=(const class Pythia8::FragmentationModifierBase &) --> class Pythia8::FragmentationModifierBase &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ColourDipole file:Pythia8/ColourReconnection.h line:36
		pybind11::class_<Pythia8::ColourDipole, std::shared_ptr<Pythia8::ColourDipole>> cl(M("Pythia8"), "ColourDipole", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ColourDipole(); } ), "doc" );
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::ColourDipole(a0); } ), "doc" , pybind11::arg("colIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Pythia8::ColourDipole(a0, a1); } ), "doc" , pybind11::arg("colIn"), pybind11::arg("iColIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2){ return new Pythia8::ColourDipole(a0, a1, a2); } ), "doc" , pybind11::arg("colIn"), pybind11::arg("iColIn"), pybind11::arg("iAcolIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3){ return new Pythia8::ColourDipole(a0, a1, a2, a3); } ), "doc" , pybind11::arg("colIn"), pybind11::arg("iColIn"), pybind11::arg("iAcolIn"), pybind11::arg("colReconnectionIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, bool const & a4){ return new Pythia8::ColourDipole(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("colIn"), pybind11::arg("iColIn"), pybind11::arg("iAcolIn"), pybind11::arg("colReconnectionIn"), pybind11::arg("isJunIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, bool const & a4, bool const & a5){ return new Pythia8::ColourDipole(a0, a1, a2, a3, a4, a5); } ), "doc" , pybind11::arg("colIn"), pybind11::arg("iColIn"), pybind11::arg("iAcolIn"), pybind11::arg("colReconnectionIn"), pybind11::arg("isJunIn"), pybind11::arg("isAntiJunIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3, bool const & a4, bool const & a5, bool const & a6){ return new Pythia8::ColourDipole(a0, a1, a2, a3, a4, a5, a6); } ), "doc" , pybind11::arg("colIn"), pybind11::arg("iColIn"), pybind11::arg("iAcolIn"), pybind11::arg("colReconnectionIn"), pybind11::arg("isJunIn"), pybind11::arg("isAntiJunIn"), pybind11::arg("isActiveIn"));
		cl.def( pybind11::init<int, int, int, int, bool, bool, bool, bool>(), pybind11::arg("colIn"), pybind11::arg("iColIn"), pybind11::arg("iAcolIn"), pybind11::arg("colReconnectionIn"), pybind11::arg("isJunIn"), pybind11::arg("isAntiJunIn"), pybind11::arg("isActiveIn"), pybind11::arg("isRealIn") );

		cl.def( pybind11::init( [](Pythia8::ColourDipole const &o){ return new Pythia8::ColourDipole(o); } ) );
		cl.def_readwrite("col", &Pythia8::ColourDipole::col);
		cl.def_readwrite("iCol", &Pythia8::ColourDipole::iCol);
		cl.def_readwrite("iAcol", &Pythia8::ColourDipole::iAcol);
		cl.def_readwrite("iColLeg", &Pythia8::ColourDipole::iColLeg);
		cl.def_readwrite("iAcolLeg", &Pythia8::ColourDipole::iAcolLeg);
		cl.def_readwrite("colReconnection", &Pythia8::ColourDipole::colReconnection);
		cl.def_readwrite("isJun", &Pythia8::ColourDipole::isJun);
		cl.def_readwrite("isAntiJun", &Pythia8::ColourDipole::isAntiJun);
		cl.def_readwrite("isActive", &Pythia8::ColourDipole::isActive);
		cl.def_readwrite("isReal", &Pythia8::ColourDipole::isReal);
		cl.def_readwrite("printed", &Pythia8::ColourDipole::printed);
		cl.def_readwrite("leftDip", &Pythia8::ColourDipole::leftDip);
		cl.def_readwrite("rightDip", &Pythia8::ColourDipole::rightDip);
		cl.def_readwrite("colDips", &Pythia8::ColourDipole::colDips);
		cl.def_readwrite("acolDips", &Pythia8::ColourDipole::acolDips);
		cl.def_readwrite("p1p2", &Pythia8::ColourDipole::p1p2);
		cl.def_readwrite("dipoleMomentum", &Pythia8::ColourDipole::dipoleMomentum);
		cl.def_readwrite("ciCol", &Pythia8::ColourDipole::ciCol);
		cl.def_readwrite("ciAcol", &Pythia8::ColourDipole::ciAcol);
		cl.def_readwrite("pCalculated", &Pythia8::ColourDipole::pCalculated);
		cl.def_readwrite("index", &Pythia8::ColourDipole::index);
		cl.def("mDip", (double (Pythia8::ColourDipole::*)(class Pythia8::Event &)) &Pythia8::ColourDipole::mDip, "C++: Pythia8::ColourDipole::mDip(class Pythia8::Event &) --> double", pybind11::arg("event"));
		cl.def("list", (void (Pythia8::ColourDipole::*)() const) &Pythia8::ColourDipole::list, "C++: Pythia8::ColourDipole::list() const --> void");
	}
	{ // Pythia8::ColourJunction file:Pythia8/ColourReconnection.h line:80
		pybind11::class_<Pythia8::ColourJunction, std::shared_ptr<Pythia8::ColourJunction>, Pythia8::Junction> cl(M("Pythia8"), "ColourJunction", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<const class Pythia8::Junction &>(), pybind11::arg("ju") );

		cl.def( pybind11::init( [](Pythia8::ColourJunction const &o){ return new Pythia8::ColourJunction(o); } ) );
		cl.def("assign", (class Pythia8::ColourJunction & (Pythia8::ColourJunction::*)(const class Pythia8::ColourJunction &)) &Pythia8::ColourJunction::operator=, "C++: Pythia8::ColourJunction::operator=(const class Pythia8::ColourJunction &) --> class Pythia8::ColourJunction &", pybind11::return_value_policy::reference, pybind11::arg("ju"));
		cl.def("getColDip", (class std::shared_ptr<class Pythia8::ColourDipole> (Pythia8::ColourJunction::*)(int)) &Pythia8::ColourJunction::getColDip, "C++: Pythia8::ColourJunction::getColDip(int) --> class std::shared_ptr<class Pythia8::ColourDipole>", pybind11::arg("i"));
		cl.def("setColDip", (void (Pythia8::ColourJunction::*)(int, class std::shared_ptr<class Pythia8::ColourDipole>)) &Pythia8::ColourJunction::setColDip, "C++: Pythia8::ColourJunction::setColDip(int, class std::shared_ptr<class Pythia8::ColourDipole>) --> void", pybind11::arg("i"), pybind11::arg("dip"));
		cl.def("list", (void (Pythia8::ColourJunction::*)() const) &Pythia8::ColourJunction::list, "C++: Pythia8::ColourJunction::list() const --> void");
	}
	{ // Pythia8::TrialReconnection file:Pythia8/ColourReconnection.h line:112
		pybind11::class_<Pythia8::TrialReconnection, std::shared_ptr<Pythia8::TrialReconnection>> cl(M("Pythia8"), "TrialReconnection", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::TrialReconnection(); } ), "doc" );
		cl.def( pybind11::init( [](class std::shared_ptr<class Pythia8::ColourDipole> const & a0){ return new Pythia8::TrialReconnection(a0); } ), "doc" , pybind11::arg("dip1In"));
		cl.def( pybind11::init( [](class std::shared_ptr<class Pythia8::ColourDipole> const & a0, class std::shared_ptr<class Pythia8::ColourDipole> const & a1){ return new Pythia8::TrialReconnection(a0, a1); } ), "doc" , pybind11::arg("dip1In"), pybind11::arg("dip2In"));
		cl.def( pybind11::init( [](class std::shared_ptr<class Pythia8::ColourDipole> const & a0, class std::shared_ptr<class Pythia8::ColourDipole> const & a1, class std::shared_ptr<class Pythia8::ColourDipole> const & a2){ return new Pythia8::TrialReconnection(a0, a1, a2); } ), "doc" , pybind11::arg("dip1In"), pybind11::arg("dip2In"), pybind11::arg("dip3In"));
		cl.def( pybind11::init( [](class std::shared_ptr<class Pythia8::ColourDipole> const & a0, class std::shared_ptr<class Pythia8::ColourDipole> const & a1, class std::shared_ptr<class Pythia8::ColourDipole> const & a2, class std::shared_ptr<class Pythia8::ColourDipole> const & a3){ return new Pythia8::TrialReconnection(a0, a1, a2, a3); } ), "doc" , pybind11::arg("dip1In"), pybind11::arg("dip2In"), pybind11::arg("dip3In"), pybind11::arg("dip4In"));
		cl.def( pybind11::init( [](class std::shared_ptr<class Pythia8::ColourDipole> const & a0, class std::shared_ptr<class Pythia8::ColourDipole> const & a1, class std::shared_ptr<class Pythia8::ColourDipole> const & a2, class std::shared_ptr<class Pythia8::ColourDipole> const & a3, int const & a4){ return new Pythia8::TrialReconnection(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("dip1In"), pybind11::arg("dip2In"), pybind11::arg("dip3In"), pybind11::arg("dip4In"), pybind11::arg("modeIn"));
		cl.def( pybind11::init<class std::shared_ptr<class Pythia8::ColourDipole>, class std::shared_ptr<class Pythia8::ColourDipole>, class std::shared_ptr<class Pythia8::ColourDipole>, class std::shared_ptr<class Pythia8::ColourDipole>, int, double>(), pybind11::arg("dip1In"), pybind11::arg("dip2In"), pybind11::arg("dip3In"), pybind11::arg("dip4In"), pybind11::arg("modeIn"), pybind11::arg("lambdaDiffIn") );

		cl.def_readwrite("dips", &Pythia8::TrialReconnection::dips);
		cl.def_readwrite("mode", &Pythia8::TrialReconnection::mode);
		cl.def_readwrite("lambdaDiff", &Pythia8::TrialReconnection::lambdaDiff);
		cl.def("list", (void (Pythia8::TrialReconnection::*)()) &Pythia8::TrialReconnection::list, "C++: Pythia8::TrialReconnection::list() --> void");
	}
	{ // Pythia8::ColourParticle file:Pythia8/ColourReconnection.h line:142
		pybind11::class_<Pythia8::ColourParticle, std::shared_ptr<Pythia8::ColourParticle>, PyCallBack_Pythia8_ColourParticle, Pythia8::Particle> cl(M("Pythia8"), "ColourParticle", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<const class Pythia8::Particle &>(), pybind11::arg("ju") );

		cl.def_readwrite("dips", &Pythia8::ColourParticle::dips);
		cl.def_readwrite("colEndIncluded", &Pythia8::ColourParticle::colEndIncluded);
		cl.def_readwrite("acolEndIncluded", &Pythia8::ColourParticle::acolEndIncluded);
		cl.def_readwrite("activeDips", &Pythia8::ColourParticle::activeDips);
		cl.def_readwrite("isJun", &Pythia8::ColourParticle::isJun);
		cl.def_readwrite("junKind", &Pythia8::ColourParticle::junKind);
		cl.def("listParticle", (void (Pythia8::ColourParticle::*)()) &Pythia8::ColourParticle::listParticle, "C++: Pythia8::ColourParticle::listParticle() --> void");
		cl.def("listActiveDips", (void (Pythia8::ColourParticle::*)()) &Pythia8::ColourParticle::listActiveDips, "C++: Pythia8::ColourParticle::listActiveDips() --> void");
		cl.def("listDips", (void (Pythia8::ColourParticle::*)()) &Pythia8::ColourParticle::listDips, "C++: Pythia8::ColourParticle::listDips() --> void");
		cl.def("assign", (class Pythia8::ColourParticle & (Pythia8::ColourParticle::*)(const class Pythia8::ColourParticle &)) &Pythia8::ColourParticle::operator=, "C++: Pythia8::ColourParticle::operator=(const class Pythia8::ColourParticle &) --> class Pythia8::ColourParticle &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ColourReconnection file:Pythia8/ColourReconnection.h line:167
		pybind11::class_<Pythia8::ColourReconnection, std::shared_ptr<Pythia8::ColourReconnection>, PyCallBack_Pythia8_ColourReconnection, Pythia8::ColourReconnectionBase> cl(M("Pythia8"), "ColourReconnection", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ColourReconnection(); }, [](){ return new PyCallBack_Pythia8_ColourReconnection(); } ) );
		cl.def("init", (bool (Pythia8::ColourReconnection::*)()) &Pythia8::ColourReconnection::init, "C++: Pythia8::ColourReconnection::init() --> bool");
		cl.def("reassignBeamPtrs", (void (Pythia8::ColourReconnection::*)(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *)) &Pythia8::ColourReconnection::reassignBeamPtrs, "C++: Pythia8::ColourReconnection::reassignBeamPtrs(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *) --> void", pybind11::arg("beamAPtrIn"), pybind11::arg("beamBPtrIn"));
		cl.def("next", (bool (Pythia8::ColourReconnection::*)(class Pythia8::Event &, int)) &Pythia8::ColourReconnection::next, "C++: Pythia8::ColourReconnection::next(class Pythia8::Event &, int) --> bool", pybind11::arg("event"), pybind11::arg("oldSize"));
		cl.def("assign", (class Pythia8::ColourReconnection & (Pythia8::ColourReconnection::*)(const class Pythia8::ColourReconnection &)) &Pythia8::ColourReconnection::operator=, "C++: Pythia8::ColourReconnection::operator=(const class Pythia8::ColourReconnection &) --> class Pythia8::ColourReconnection &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::BeamRemnants file:Pythia8/BeamRemnants.h line:35
		pybind11::class_<Pythia8::BeamRemnants, std::shared_ptr<Pythia8::BeamRemnants>, PyCallBack_Pythia8_BeamRemnants, Pythia8::PhysicsBase> cl(M("Pythia8"), "BeamRemnants", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::BeamRemnants(); }, [](){ return new PyCallBack_Pythia8_BeamRemnants(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_BeamRemnants const &o){ return new PyCallBack_Pythia8_BeamRemnants(o); } ) );
		cl.def( pybind11::init( [](Pythia8::BeamRemnants const &o){ return new Pythia8::BeamRemnants(o); } ) );
		cl.def("init", (bool (Pythia8::BeamRemnants::*)(class std::shared_ptr<class Pythia8::PartonVertex>, class std::shared_ptr<class Pythia8::ColourReconnectionBase>)) &Pythia8::BeamRemnants::init, "C++: Pythia8::BeamRemnants::init(class std::shared_ptr<class Pythia8::PartonVertex>, class std::shared_ptr<class Pythia8::ColourReconnectionBase>) --> bool", pybind11::arg("partonVertexPtrIn"), pybind11::arg("colourReconnectionPtrIn"));
		cl.def("reassignBeamPtrs", (void (Pythia8::BeamRemnants::*)(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *, int)) &Pythia8::BeamRemnants::reassignBeamPtrs, "C++: Pythia8::BeamRemnants::reassignBeamPtrs(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *, int) --> void", pybind11::arg("beamAPtrIn"), pybind11::arg("beamBPtrIn"), pybind11::arg("iDSin"));
		cl.def("add", [](Pythia8::BeamRemnants &o, class Pythia8::Event & a0) -> bool { return o.add(a0); }, "", pybind11::arg("event"));
		cl.def("add", [](Pythia8::BeamRemnants &o, class Pythia8::Event & a0, int const & a1) -> bool { return o.add(a0, a1); }, "", pybind11::arg("event"), pybind11::arg("iFirst"));
		cl.def("add", (bool (Pythia8::BeamRemnants::*)(class Pythia8::Event &, int, bool)) &Pythia8::BeamRemnants::add, "C++: Pythia8::BeamRemnants::add(class Pythia8::Event &, int, bool) --> bool", pybind11::arg("event"), pybind11::arg("iFirst"), pybind11::arg("doDiffCR"));
		cl.def("onInitInfoPtr", (void (Pythia8::BeamRemnants::*)()) &Pythia8::BeamRemnants::onInitInfoPtr, "C++: Pythia8::BeamRemnants::onInitInfoPtr() --> void");
		cl.def("assign", (class Pythia8::BeamRemnants & (Pythia8::BeamRemnants::*)(const class Pythia8::BeamRemnants &)) &Pythia8::BeamRemnants::operator=, "C++: Pythia8::BeamRemnants::operator=(const class Pythia8::BeamRemnants &) --> class Pythia8::BeamRemnants &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
