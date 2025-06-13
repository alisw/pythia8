#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/Info.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/SimpleSpaceShower.h>
#include <Pythia8/SimpleTimeShower.h>
#include <Pythia8/SpaceShower.h>
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

// Pythia8::SimpleSpaceShower file:Pythia8/SimpleSpaceShower.h line:74
struct PyCallBack_Pythia8_SimpleSpaceShower : public Pythia8::SimpleSpaceShower {
	using Pythia8::SimpleSpaceShower::SimpleSpaceShower;

	void init(class Pythia8::BeamParticle * a0, class Pythia8::BeamParticle * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleSpaceShower::init(a0, a1);
	}
	bool limitPTmax(class Pythia8::Event & a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "limitPTmax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleSpaceShower::limitPTmax(a0, a1, a2);
	}
	void prepare(int a0, class Pythia8::Event & a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "prepare");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleSpaceShower::prepare(a0, a1, a2);
	}
	void update(int a0, class Pythia8::Event & a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleSpaceShower::update(a0, a1, a2);
	}
	double pTnext(class Pythia8::Event & a0, double a1, double a2, int a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "pTnext");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SimpleSpaceShower::pTnext(a0, a1, a2, a3, a4);
	}
	bool branch(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "branch");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleSpaceShower::branch(a0);
	}
	void list() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "list");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleSpaceShower::list();
	}
	bool initUncertainties() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "initUncertainties");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleSpaceShower::initUncertainties();
	}
	bool initEnhancements() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "initEnhancements");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleSpaceShower::initEnhancements();
	}
	bool doRestart() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "doRestart");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleSpaceShower::doRestart();
	}
	bool wasGamma2qqbar() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "wasGamma2qqbar");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleSpaceShower::wasGamma2qqbar();
	}
	bool getHasWeaklyRadiated() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "getHasWeaklyRadiated");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleSpaceShower::getHasWeaklyRadiated();
	}
	int system() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "system");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return SimpleSpaceShower::system();
	}
	double enhancePTmax() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "enhancePTmax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SimpleSpaceShower::enhancePTmax();
	}
	double noEmissionProbability(double a0, double a1, double a2, int a3, int a4, double a5, double a6) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "noEmissionProbability");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SimpleSpaceShower::noEmissionProbability(a0, a1, a2, a3, a4, a5, a6);
	}
	class Pythia8::Event clustered(const class Pythia8::Event & a0, int a1, int a2, int a3, class std::basic_string<char> a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "clustered");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Event>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Event> caster;
				return pybind11::detail::cast_ref<class Pythia8::Event>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Event>(std::move(o));
		}
		return SpaceShower::clustered(a0, a1, a2, a3, a4);
	}
	using _binder_ret_0 = class std::map<class std::basic_string<char>, double, struct std::less<class std::basic_string<char> >, class std::allocator<struct std::pair<const class std::basic_string<char>, double> > >;
	_binder_ret_0 getStateVariables(const class Pythia8::Event & a0, int a1, int a2, int a3, class std::basic_string<char> a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "getStateVariables");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return SpaceShower::getStateVariables(a0, a1, a2, a3, a4);
	}
	bool isSpacelike(const class Pythia8::Event & a0, int a1, int a2, int a3, class std::basic_string<char> a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "isSpacelike");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SpaceShower::isSpacelike(a0, a1, a2, a3, a4);
	}
	using _binder_ret_1 = class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > >;
	_binder_ret_1 getSplittingName(const class Pythia8::Event & a0, int a1, int a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "getSplittingName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		return SpaceShower::getSplittingName(a0, a1, a2, a3);
	}
	double getSplittingProb(const class Pythia8::Event & a0, int a1, int a2, int a3, class std::basic_string<char> a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "getSplittingProb");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SpaceShower::getSplittingProb(a0, a1, a2, a3, a4);
	}
	bool allowedSplitting(const class Pythia8::Event & a0, int a1, int a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "allowedSplitting");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SpaceShower::allowedSplitting(a0, a1, a2);
	}
	using _binder_ret_2 = class std::vector<int, class std::allocator<int> >;
	_binder_ret_2 getRecoilers(const class Pythia8::Event & a0, int a1, int a2, class std::basic_string<char> a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "getRecoilers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		return SpaceShower::getRecoilers(a0, a1, a2, a3);
	}
	double enhanceFactor(const class std::basic_string<char> & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "enhanceFactor");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SpaceShower::enhanceFactor(a0);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleSpaceShower *>(this), "onStat");
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

void bind_Pythia8_SimpleSpaceShower(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::SimpleSpaceShower file:Pythia8/SimpleSpaceShower.h line:74
		pybind11::class_<Pythia8::SimpleSpaceShower, std::shared_ptr<Pythia8::SimpleSpaceShower>, PyCallBack_Pythia8_SimpleSpaceShower, Pythia8::SpaceShower> cl(M("Pythia8"), "SimpleSpaceShower", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SimpleSpaceShower(); }, [](){ return new PyCallBack_Pythia8_SimpleSpaceShower(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_SimpleSpaceShower const &o){ return new PyCallBack_Pythia8_SimpleSpaceShower(o); } ) );
		cl.def( pybind11::init( [](Pythia8::SimpleSpaceShower const &o){ return new Pythia8::SimpleSpaceShower(o); } ) );
		cl.def_readwrite("pdfMode", &Pythia8::SimpleSpaceShower::pdfMode);
		cl.def("init", (void (Pythia8::SimpleSpaceShower::*)(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *)) &Pythia8::SimpleSpaceShower::init, "C++: Pythia8::SimpleSpaceShower::init(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *) --> void", pybind11::arg("beamAPtrIn"), pybind11::arg("beamBPtrIn"));
		cl.def("limitPTmax", [](Pythia8::SimpleSpaceShower &o, class Pythia8::Event & a0) -> bool { return o.limitPTmax(a0); }, "", pybind11::arg("event"));
		cl.def("limitPTmax", [](Pythia8::SimpleSpaceShower &o, class Pythia8::Event & a0, double const & a1) -> bool { return o.limitPTmax(a0, a1); }, "", pybind11::arg("event"), pybind11::arg("Q2Fac"));
		cl.def("limitPTmax", (bool (Pythia8::SimpleSpaceShower::*)(class Pythia8::Event &, double, double)) &Pythia8::SimpleSpaceShower::limitPTmax, "C++: Pythia8::SimpleSpaceShower::limitPTmax(class Pythia8::Event &, double, double) --> bool", pybind11::arg("event"), pybind11::arg("Q2Fac"), pybind11::arg("Q2Ren"));
		cl.def("prepare", [](Pythia8::SimpleSpaceShower &o, int const & a0, class Pythia8::Event & a1) -> void { return o.prepare(a0, a1); }, "", pybind11::arg("iSys"), pybind11::arg("event"));
		cl.def("prepare", (void (Pythia8::SimpleSpaceShower::*)(int, class Pythia8::Event &, bool)) &Pythia8::SimpleSpaceShower::prepare, "C++: Pythia8::SimpleSpaceShower::prepare(int, class Pythia8::Event &, bool) --> void", pybind11::arg("iSys"), pybind11::arg("event"), pybind11::arg("limitPTmaxIn"));
		cl.def("update", [](Pythia8::SimpleSpaceShower &o, int const & a0, class Pythia8::Event & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("iSys"), pybind11::arg("event"));
		cl.def("update", (void (Pythia8::SimpleSpaceShower::*)(int, class Pythia8::Event &, bool)) &Pythia8::SimpleSpaceShower::update, "C++: Pythia8::SimpleSpaceShower::update(int, class Pythia8::Event &, bool) --> void", pybind11::arg("iSys"), pybind11::arg("event"), pybind11::arg("hasWeakRad"));
		cl.def("pTnext", [](Pythia8::SimpleSpaceShower &o, class Pythia8::Event & a0, double const & a1, double const & a2) -> double { return o.pTnext(a0, a1, a2); }, "", pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"));
		cl.def("pTnext", [](Pythia8::SimpleSpaceShower &o, class Pythia8::Event & a0, double const & a1, double const & a2, int const & a3) -> double { return o.pTnext(a0, a1, a2, a3); }, "", pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("nRadIn"));
		cl.def("pTnext", (double (Pythia8::SimpleSpaceShower::*)(class Pythia8::Event &, double, double, int, bool)) &Pythia8::SimpleSpaceShower::pTnext, "C++: Pythia8::SimpleSpaceShower::pTnext(class Pythia8::Event &, double, double, int, bool) --> double", pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("nRadIn"), pybind11::arg("doTrialIn"));
		cl.def("branch", (bool (Pythia8::SimpleSpaceShower::*)(class Pythia8::Event &)) &Pythia8::SimpleSpaceShower::branch, "C++: Pythia8::SimpleSpaceShower::branch(class Pythia8::Event &) --> bool", pybind11::arg("event"));
		cl.def("list", (void (Pythia8::SimpleSpaceShower::*)() const) &Pythia8::SimpleSpaceShower::list, "C++: Pythia8::SimpleSpaceShower::list() const --> void");
		cl.def("initUncertainties", (bool (Pythia8::SimpleSpaceShower::*)()) &Pythia8::SimpleSpaceShower::initUncertainties, "C++: Pythia8::SimpleSpaceShower::initUncertainties() --> bool");
		cl.def("initEnhancements", (bool (Pythia8::SimpleSpaceShower::*)()) &Pythia8::SimpleSpaceShower::initEnhancements, "C++: Pythia8::SimpleSpaceShower::initEnhancements() --> bool");
		cl.def("doRestart", (bool (Pythia8::SimpleSpaceShower::*)() const) &Pythia8::SimpleSpaceShower::doRestart, "C++: Pythia8::SimpleSpaceShower::doRestart() const --> bool");
		cl.def("wasGamma2qqbar", (bool (Pythia8::SimpleSpaceShower::*)()) &Pythia8::SimpleSpaceShower::wasGamma2qqbar, "C++: Pythia8::SimpleSpaceShower::wasGamma2qqbar() --> bool");
		cl.def("getHasWeaklyRadiated", (bool (Pythia8::SimpleSpaceShower::*)()) &Pythia8::SimpleSpaceShower::getHasWeaklyRadiated, "C++: Pythia8::SimpleSpaceShower::getHasWeaklyRadiated() --> bool");
		cl.def("system", (int (Pythia8::SimpleSpaceShower::*)() const) &Pythia8::SimpleSpaceShower::system, "C++: Pythia8::SimpleSpaceShower::system() const --> int");
		cl.def("enhancePTmax", (double (Pythia8::SimpleSpaceShower::*)() const) &Pythia8::SimpleSpaceShower::enhancePTmax, "C++: Pythia8::SimpleSpaceShower::enhancePTmax() const --> double");
		cl.def("noEmissionProbability", [](Pythia8::SimpleSpaceShower &o, double const & a0, double const & a1, double const & a2, int const & a3, int const & a4) -> double { return o.noEmissionProbability(a0, a1, a2, a3, a4); }, "", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"));
		cl.def("noEmissionProbability", [](Pythia8::SimpleSpaceShower &o, double const & a0, double const & a1, double const & a2, int const & a3, int const & a4, double const & a5) -> double { return o.noEmissionProbability(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"), pybind11::arg("s"));
		cl.def("noEmissionProbability", (double (Pythia8::SimpleSpaceShower::*)(double, double, double, int, int, double, double)) &Pythia8::SimpleSpaceShower::noEmissionProbability, "C++: Pythia8::SimpleSpaceShower::noEmissionProbability(double, double, double, int, int, double, double) --> double", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"), pybind11::arg("s"), pybind11::arg("x"));
		cl.def("pTnext", [](Pythia8::SimpleSpaceShower &o, class std::vector<class Pythia8::SpaceDipoleEnd, class std::allocator<class Pythia8::SpaceDipoleEnd> > const & a0, class Pythia8::Event const & a1, double const & a2, double const & a3, double const & a4, int const & a5, int const & a6) -> double { return o.pTnext(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("dipEnds"), pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"));
		cl.def("pTnext", [](Pythia8::SimpleSpaceShower &o, class std::vector<class Pythia8::SpaceDipoleEnd, class std::allocator<class Pythia8::SpaceDipoleEnd> > const & a0, class Pythia8::Event const & a1, double const & a2, double const & a3, double const & a4, int const & a5, int const & a6, double const & a7) -> double { return o.pTnext(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("dipEnds"), pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"), pybind11::arg("s"));
		cl.def("pTnext", (double (Pythia8::SimpleSpaceShower::*)(class std::vector<class Pythia8::SpaceDipoleEnd, class std::allocator<class Pythia8::SpaceDipoleEnd> >, class Pythia8::Event, double, double, double, int, int, double, double)) &Pythia8::SimpleSpaceShower::pTnext, "C++: Pythia8::SimpleSpaceShower::pTnext(class std::vector<class Pythia8::SpaceDipoleEnd, class std::allocator<class Pythia8::SpaceDipoleEnd> >, class Pythia8::Event, double, double, double, int, int, double, double) --> double", pybind11::arg("dipEnds"), pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"), pybind11::arg("s"), pybind11::arg("x"));
		cl.def("assign", (class Pythia8::SimpleSpaceShower & (Pythia8::SimpleSpaceShower::*)(const class Pythia8::SimpleSpaceShower &)) &Pythia8::SimpleSpaceShower::operator=, "C++: Pythia8::SimpleSpaceShower::operator=(const class Pythia8::SimpleSpaceShower &) --> class Pythia8::SimpleSpaceShower &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::TimeDipoleEnd file:Pythia8/SimpleTimeShower.h line:25
		pybind11::class_<Pythia8::TimeDipoleEnd, std::shared_ptr<Pythia8::TimeDipoleEnd>> cl(M("Pythia8"), "TimeDipoleEnd", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::TimeDipoleEnd(); } ) );
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Pythia8::TimeDipoleEnd(a0, a1); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2){ return new Pythia8::TimeDipoleEnd(a0, a1, a2); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11, int const & a12){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11, int const & a12, bool const & a13){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"), pybind11::arg("isHiddenValleyIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11, int const & a12, bool const & a13, int const & a14){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"), pybind11::arg("isHiddenValleyIn"), pybind11::arg("colvTypeIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11, int const & a12, bool const & a13, int const & a14, double const & a15){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"), pybind11::arg("isHiddenValleyIn"), pybind11::arg("colvTypeIn"), pybind11::arg("MEmixIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11, int const & a12, bool const & a13, int const & a14, double const & a15, bool const & a16){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"), pybind11::arg("isHiddenValleyIn"), pybind11::arg("colvTypeIn"), pybind11::arg("MEmixIn"), pybind11::arg("MEorderIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11, int const & a12, bool const & a13, int const & a14, double const & a15, bool const & a16, bool const & a17){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"), pybind11::arg("isHiddenValleyIn"), pybind11::arg("colvTypeIn"), pybind11::arg("MEmixIn"), pybind11::arg("MEorderIn"), pybind11::arg("MEsplitIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11, int const & a12, bool const & a13, int const & a14, double const & a15, bool const & a16, bool const & a17, bool const & a18){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"), pybind11::arg("isHiddenValleyIn"), pybind11::arg("colvTypeIn"), pybind11::arg("MEmixIn"), pybind11::arg("MEorderIn"), pybind11::arg("MEsplitIn"), pybind11::arg("MEgluinoRecIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, double const & a2, int const & a3, int const & a4, int const & a5, int const & a6, int const & a7, int const & a8, int const & a9, int const & a10, int const & a11, int const & a12, bool const & a13, int const & a14, double const & a15, bool const & a16, bool const & a17, bool const & a18, bool const & a19){ return new Pythia8::TimeDipoleEnd(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19); } ), "doc" , pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"), pybind11::arg("isHiddenValleyIn"), pybind11::arg("colvTypeIn"), pybind11::arg("MEmixIn"), pybind11::arg("MEorderIn"), pybind11::arg("MEsplitIn"), pybind11::arg("MEgluinoRecIn"), pybind11::arg("isFlexibleIn"));
		cl.def( pybind11::init<int, int, double, int, int, int, int, int, int, int, int, int, int, bool, int, double, bool, bool, bool, bool, bool>(), pybind11::arg("iRadiatorIn"), pybind11::arg("iRecoilerIn"), pybind11::arg("pTmaxIn"), pybind11::arg("colIn"), pybind11::arg("chgIn"), pybind11::arg("gamIn"), pybind11::arg("weakTypeIn"), pybind11::arg("isrIn"), pybind11::arg("systemIn"), pybind11::arg("MEtypeIn"), pybind11::arg("iMEpartnerIn"), pybind11::arg("weakPolIn"), pybind11::arg("oniumTypeIn"), pybind11::arg("isHiddenValleyIn"), pybind11::arg("colvTypeIn"), pybind11::arg("MEmixIn"), pybind11::arg("MEorderIn"), pybind11::arg("MEsplitIn"), pybind11::arg("MEgluinoRecIn"), pybind11::arg("isFlexibleIn"), pybind11::arg("hasJunctionIn") );

		cl.def_readwrite("iRadiator", &Pythia8::TimeDipoleEnd::iRadiator);
		cl.def_readwrite("iRecoiler", &Pythia8::TimeDipoleEnd::iRecoiler);
		cl.def_readwrite("pTmax", &Pythia8::TimeDipoleEnd::pTmax);
		cl.def_readwrite("colType", &Pythia8::TimeDipoleEnd::colType);
		cl.def_readwrite("chgType", &Pythia8::TimeDipoleEnd::chgType);
		cl.def_readwrite("gamType", &Pythia8::TimeDipoleEnd::gamType);
		cl.def_readwrite("weakType", &Pythia8::TimeDipoleEnd::weakType);
		cl.def_readwrite("isrType", &Pythia8::TimeDipoleEnd::isrType);
		cl.def_readwrite("system", &Pythia8::TimeDipoleEnd::system);
		cl.def_readwrite("systemRec", &Pythia8::TimeDipoleEnd::systemRec);
		cl.def_readwrite("MEtype", &Pythia8::TimeDipoleEnd::MEtype);
		cl.def_readwrite("iMEpartner", &Pythia8::TimeDipoleEnd::iMEpartner);
		cl.def_readwrite("weakPol", &Pythia8::TimeDipoleEnd::weakPol);
		cl.def_readwrite("oniumType", &Pythia8::TimeDipoleEnd::oniumType);
		cl.def_readwrite("isHiddenValley", &Pythia8::TimeDipoleEnd::isHiddenValley);
		cl.def_readwrite("colvType", &Pythia8::TimeDipoleEnd::colvType);
		cl.def_readwrite("MEmix", &Pythia8::TimeDipoleEnd::MEmix);
		cl.def_readwrite("MEorder", &Pythia8::TimeDipoleEnd::MEorder);
		cl.def_readwrite("MEsplit", &Pythia8::TimeDipoleEnd::MEsplit);
		cl.def_readwrite("MEgluinoRec", &Pythia8::TimeDipoleEnd::MEgluinoRec);
		cl.def_readwrite("isFlexible", &Pythia8::TimeDipoleEnd::isFlexible);
		cl.def_readwrite("hasJunction", &Pythia8::TimeDipoleEnd::hasJunction);
		cl.def_readwrite("flavour", &Pythia8::TimeDipoleEnd::flavour);
		cl.def_readwrite("iAunt", &Pythia8::TimeDipoleEnd::iAunt);
		cl.def_readwrite("mRad", &Pythia8::TimeDipoleEnd::mRad);
		cl.def_readwrite("m2Rad", &Pythia8::TimeDipoleEnd::m2Rad);
		cl.def_readwrite("mRec", &Pythia8::TimeDipoleEnd::mRec);
		cl.def_readwrite("m2Rec", &Pythia8::TimeDipoleEnd::m2Rec);
		cl.def_readwrite("mDip", &Pythia8::TimeDipoleEnd::mDip);
		cl.def_readwrite("m2Dip", &Pythia8::TimeDipoleEnd::m2Dip);
		cl.def_readwrite("m2DipCorr", &Pythia8::TimeDipoleEnd::m2DipCorr);
		cl.def_readwrite("pT2", &Pythia8::TimeDipoleEnd::pT2);
		cl.def_readwrite("m2", &Pythia8::TimeDipoleEnd::m2);
		cl.def_readwrite("z", &Pythia8::TimeDipoleEnd::z);
		cl.def_readwrite("mFlavour", &Pythia8::TimeDipoleEnd::mFlavour);
		cl.def_readwrite("asymPol", &Pythia8::TimeDipoleEnd::asymPol);
		cl.def_readwrite("flexFactor", &Pythia8::TimeDipoleEnd::flexFactor);
		cl.def_readwrite("pAccept", &Pythia8::TimeDipoleEnd::pAccept);
		cl.def_readwrite("m2A", &Pythia8::TimeDipoleEnd::m2A);
		cl.def_readwrite("m2B", &Pythia8::TimeDipoleEnd::m2B);
		cl.def_readwrite("m2C", &Pythia8::TimeDipoleEnd::m2C);
		cl.def_readwrite("m2gg", &Pythia8::TimeDipoleEnd::m2gg);
		cl.def_readwrite("emissionPtr", &Pythia8::TimeDipoleEnd::emissionPtr);
	}
}
