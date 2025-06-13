#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/FragmentationSystems.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/MiniStringFragmentation.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/StringFragmentation.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/Weights.h>
#include <functional>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
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

// Pythia8::LundFragmentation file:Pythia8/FragmentationModel.h line:68
struct PyCallBack_Pythia8_LundFragmentation : public Pythia8::LundFragmentation {
	using Pythia8::LundFragmentation::LundFragmentation;

	bool init(class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2, class std::shared_ptr<class Pythia8::FragmentationModifierBase> a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LundFragmentation *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LundFragmentation::init(a0, a1, a2, a3);
	}
	bool fragment(int a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LundFragmentation *>(this), "fragment");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LundFragmentation::fragment(a0, a1, a2, a3, a4);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LundFragmentation *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LundFragmentation *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LundFragmentation *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::LundFragmentation *>(this), "onStat");
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

// Pythia8::MiniStringFragmentation file:Pythia8/MiniStringFragmentation.h line:22
struct PyCallBack_Pythia8_MiniStringFragmentation : public Pythia8::MiniStringFragmentation {
	using Pythia8::MiniStringFragmentation::MiniStringFragmentation;

	bool init(class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2, class std::shared_ptr<class Pythia8::FragmentationModifierBase> a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MiniStringFragmentation *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return MiniStringFragmentation::init(a0, a1, a2, a3);
	}
	bool fragment(int a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MiniStringFragmentation *>(this), "fragment");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return MiniStringFragmentation::fragment(a0, a1, a2, a3, a4);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MiniStringFragmentation *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MiniStringFragmentation *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MiniStringFragmentation *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::MiniStringFragmentation *>(this), "onStat");
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

// Pythia8::StringFragmentation file:Pythia8/StringFragmentation.h line:106
struct PyCallBack_Pythia8_StringFragmentation : public Pythia8::StringFragmentation {
	using Pythia8::StringFragmentation::StringFragmentation;

	bool init(class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2, class std::shared_ptr<class Pythia8::FragmentationModifierBase> a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFragmentation *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringFragmentation::init(a0, a1, a2, a3);
	}
	bool fragment(int a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFragmentation *>(this), "fragment");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StringFragmentation::fragment(a0, a1, a2, a3, a4);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFragmentation *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFragmentation *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFragmentation *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::StringFragmentation *>(this), "onStat");
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

void bind_Pythia8_FragmentationModel(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::LundFragmentation file:Pythia8/FragmentationModel.h line:68
		pybind11::class_<Pythia8::LundFragmentation, std::shared_ptr<Pythia8::LundFragmentation>, PyCallBack_Pythia8_LundFragmentation, Pythia8::FragmentationModel> cl(M("Pythia8"), "LundFragmentation", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LundFragmentation(); }, [](){ return new PyCallBack_Pythia8_LundFragmentation(); } ) );
		cl.def("init", [](Pythia8::LundFragmentation &o) -> bool { return o.init(); }, "");
		cl.def("init", [](Pythia8::LundFragmentation &o, class Pythia8::StringFlav * a0) -> bool { return o.init(a0); }, "", pybind11::arg("flavSelPtrIn"));
		cl.def("init", [](Pythia8::LundFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1) -> bool { return o.init(a0, a1); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"));
		cl.def("init", [](Pythia8::LundFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2) -> bool { return o.init(a0, a1, a2); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"));
		cl.def("init", (bool (Pythia8::LundFragmentation::*)(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>)) &Pythia8::LundFragmentation::init, "C++: Pythia8::LundFragmentation::init(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>) --> bool", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"), pybind11::arg("fragModPtrIn"));
		cl.def("fragment", [](Pythia8::LundFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2) -> bool { return o.fragment(a0, a1, a2); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"));
		cl.def("fragment", [](Pythia8::LundFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool const & a3) -> bool { return o.fragment(a0, a1, a2, a3); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"));
		cl.def("fragment", (bool (Pythia8::LundFragmentation::*)(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool)) &Pythia8::LundFragmentation::fragment, "C++: Pythia8::LundFragmentation::fragment(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool) --> bool", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"), pybind11::arg("systemRecoil"));
		cl.def("assign", (class Pythia8::LundFragmentation & (Pythia8::LundFragmentation::*)(const class Pythia8::LundFragmentation &)) &Pythia8::LundFragmentation::operator=, "C++: Pythia8::LundFragmentation::operator=(const class Pythia8::LundFragmentation &) --> class Pythia8::LundFragmentation &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::MiniStringFragmentation file:Pythia8/MiniStringFragmentation.h line:22
		pybind11::class_<Pythia8::MiniStringFragmentation, std::shared_ptr<Pythia8::MiniStringFragmentation>, PyCallBack_Pythia8_MiniStringFragmentation, Pythia8::FragmentationModel> cl(M("Pythia8"), "MiniStringFragmentation", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::MiniStringFragmentation(); }, [](){ return new PyCallBack_Pythia8_MiniStringFragmentation(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_MiniStringFragmentation const &o){ return new PyCallBack_Pythia8_MiniStringFragmentation(o); } ) );
		cl.def( pybind11::init( [](Pythia8::MiniStringFragmentation const &o){ return new Pythia8::MiniStringFragmentation(o); } ) );
		cl.def("init", [](Pythia8::MiniStringFragmentation &o) -> bool { return o.init(); }, "");
		cl.def("init", [](Pythia8::MiniStringFragmentation &o, class Pythia8::StringFlav * a0) -> bool { return o.init(a0); }, "", pybind11::arg("flavSelPtrIn"));
		cl.def("init", [](Pythia8::MiniStringFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1) -> bool { return o.init(a0, a1); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"));
		cl.def("init", [](Pythia8::MiniStringFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2) -> bool { return o.init(a0, a1, a2); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"));
		cl.def("init", (bool (Pythia8::MiniStringFragmentation::*)(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>)) &Pythia8::MiniStringFragmentation::init, "C++: Pythia8::MiniStringFragmentation::init(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>) --> bool", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"), pybind11::arg("fragModPtrIn"));
		cl.def("fragment", [](Pythia8::MiniStringFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2) -> bool { return o.fragment(a0, a1, a2); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"));
		cl.def("fragment", [](Pythia8::MiniStringFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool const & a3) -> bool { return o.fragment(a0, a1, a2, a3); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"));
		cl.def("fragment", (bool (Pythia8::MiniStringFragmentation::*)(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool)) &Pythia8::MiniStringFragmentation::fragment, "C++: Pythia8::MiniStringFragmentation::fragment(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool) --> bool", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"), pybind11::arg("systemRecoil"));
		cl.def("setMVecRatio", (void (Pythia8::MiniStringFragmentation::*)(double)) &Pythia8::MiniStringFragmentation::setMVecRatio, "C++: Pythia8::MiniStringFragmentation::setMVecRatio(double) --> void", pybind11::arg("mVecRatioIn"));
		cl.def("assign", (class Pythia8::MiniStringFragmentation & (Pythia8::MiniStringFragmentation::*)(const class Pythia8::MiniStringFragmentation &)) &Pythia8::MiniStringFragmentation::operator=, "C++: Pythia8::MiniStringFragmentation::operator=(const class Pythia8::MiniStringFragmentation &) --> class Pythia8::MiniStringFragmentation &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::StringEnd file:Pythia8/StringFragmentation.h line:23
		pybind11::class_<Pythia8::StringEnd, std::shared_ptr<Pythia8::StringEnd>> cl(M("Pythia8"), "StringEnd", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::StringEnd(); } ) );
		cl.def( pybind11::init( [](Pythia8::StringEnd const &o){ return new Pythia8::StringEnd(o); } ) );
		cl.def_readwrite("flavSelNow", &Pythia8::StringEnd::flavSelNow);
		cl.def_readwrite("fromPos", &Pythia8::StringEnd::fromPos);
		cl.def_readwrite("thermalModel", &Pythia8::StringEnd::thermalModel);
		cl.def_readwrite("mT2suppression", &Pythia8::StringEnd::mT2suppression);
		cl.def_readwrite("closePacking", &Pythia8::StringEnd::closePacking);
		cl.def_readwrite("iEnd", &Pythia8::StringEnd::iEnd);
		cl.def_readwrite("iMax", &Pythia8::StringEnd::iMax);
		cl.def_readwrite("idHad", &Pythia8::StringEnd::idHad);
		cl.def_readwrite("iPosOld", &Pythia8::StringEnd::iPosOld);
		cl.def_readwrite("iNegOld", &Pythia8::StringEnd::iNegOld);
		cl.def_readwrite("iPosNew", &Pythia8::StringEnd::iPosNew);
		cl.def_readwrite("iNegNew", &Pythia8::StringEnd::iNegNew);
		cl.def_readwrite("hadSoFar", &Pythia8::StringEnd::hadSoFar);
		cl.def_readwrite("colOld", &Pythia8::StringEnd::colOld);
		cl.def_readwrite("colNew", &Pythia8::StringEnd::colNew);
		cl.def_readwrite("pxOld", &Pythia8::StringEnd::pxOld);
		cl.def_readwrite("pyOld", &Pythia8::StringEnd::pyOld);
		cl.def_readwrite("pxNew", &Pythia8::StringEnd::pxNew);
		cl.def_readwrite("pyNew", &Pythia8::StringEnd::pyNew);
		cl.def_readwrite("pxHad", &Pythia8::StringEnd::pxHad);
		cl.def_readwrite("pyHad", &Pythia8::StringEnd::pyHad);
		cl.def_readwrite("mHad", &Pythia8::StringEnd::mHad);
		cl.def_readwrite("mT2Had", &Pythia8::StringEnd::mT2Had);
		cl.def_readwrite("zHad", &Pythia8::StringEnd::zHad);
		cl.def_readwrite("GammaOld", &Pythia8::StringEnd::GammaOld);
		cl.def_readwrite("GammaNew", &Pythia8::StringEnd::GammaNew);
		cl.def_readwrite("xPosOld", &Pythia8::StringEnd::xPosOld);
		cl.def_readwrite("xPosNew", &Pythia8::StringEnd::xPosNew);
		cl.def_readwrite("xPosHad", &Pythia8::StringEnd::xPosHad);
		cl.def_readwrite("xNegOld", &Pythia8::StringEnd::xNegOld);
		cl.def_readwrite("xNegNew", &Pythia8::StringEnd::xNegNew);
		cl.def_readwrite("xNegHad", &Pythia8::StringEnd::xNegHad);
		cl.def_readwrite("aLund", &Pythia8::StringEnd::aLund);
		cl.def_readwrite("bLund", &Pythia8::StringEnd::bLund);
		cl.def_readwrite("iPosOldPrev", &Pythia8::StringEnd::iPosOldPrev);
		cl.def_readwrite("iNegOldPrev", &Pythia8::StringEnd::iNegOldPrev);
		cl.def_readwrite("colOldPrev", &Pythia8::StringEnd::colOldPrev);
		cl.def_readwrite("pxOldPrev", &Pythia8::StringEnd::pxOldPrev);
		cl.def_readwrite("pyOldPrev", &Pythia8::StringEnd::pyOldPrev);
		cl.def_readwrite("GammaOldPrev", &Pythia8::StringEnd::GammaOldPrev);
		cl.def_readwrite("xPosOldPrev", &Pythia8::StringEnd::xPosOldPrev);
		cl.def_readwrite("xNegOldPrev", &Pythia8::StringEnd::xNegOldPrev);
		cl.def_readwrite("mVecRatio", &Pythia8::StringEnd::mVecRatio);
		cl.def_readwrite("tinyEq", &Pythia8::StringEnd::tinyEq);
		cl.def_readwrite("pT2tiny", &Pythia8::StringEnd::pT2tiny);
		cl.def_readwrite("flavOld", &Pythia8::StringEnd::flavOld);
		cl.def_readwrite("flavNew", &Pythia8::StringEnd::flavNew);
		cl.def_readwrite("flavOldPrev", &Pythia8::StringEnd::flavOldPrev);
		cl.def_readwrite("pHad", &Pythia8::StringEnd::pHad);
		cl.def_readwrite("pSoFar", &Pythia8::StringEnd::pSoFar);
		cl.def("init", (void (Pythia8::StringEnd::*)(class Pythia8::ParticleData *, class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class Pythia8::Settings &)) &Pythia8::StringEnd::init, "C++: Pythia8::StringEnd::init(class Pythia8::ParticleData *, class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class Pythia8::Settings &) --> void", pybind11::arg("particleDataPtrIn"), pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"), pybind11::arg("settings"));
		cl.def("setUp", (void (Pythia8::StringEnd::*)(bool, int, int, int, double, double, double, double, double, int, double)) &Pythia8::StringEnd::setUp, "C++: Pythia8::StringEnd::setUp(bool, int, int, int, double, double, double, double, double, int, double) --> void", pybind11::arg("fromPosIn"), pybind11::arg("iEndIn"), pybind11::arg("idOldIn"), pybind11::arg("iMaxIn"), pybind11::arg("pxIn"), pybind11::arg("pyIn"), pybind11::arg("GammaIn"), pybind11::arg("xPosIn"), pybind11::arg("xNegIn"), pybind11::arg("colIn"), pybind11::arg("mVecRatioIn"));
		cl.def("newHadron", [](Pythia8::StringEnd &o, double const & a0) -> void { return o.newHadron(a0); }, "", pybind11::arg("kappaModifier"));
		cl.def("newHadron", [](Pythia8::StringEnd &o, double const & a0, bool const & a1) -> void { return o.newHadron(a0, a1); }, "", pybind11::arg("kappaModifier"), pybind11::arg("forbidPopcornNow"));
		cl.def("newHadron", [](Pythia8::StringEnd &o, double const & a0, bool const & a1, double const & a2) -> void { return o.newHadron(a0, a1, a2); }, "", pybind11::arg("kappaModifier"), pybind11::arg("forbidPopcornNow"), pybind11::arg("strangeJunc"));
		cl.def("newHadron", (void (Pythia8::StringEnd::*)(double, bool, double, double)) &Pythia8::StringEnd::newHadron, "C++: Pythia8::StringEnd::newHadron(double, bool, double, double) --> void", pybind11::arg("kappaModifier"), pybind11::arg("forbidPopcornNow"), pybind11::arg("strangeJunc"), pybind11::arg("probQQmod"));
		cl.def("kinematicsHadron", (class Pythia8::Vec4 (Pythia8::StringEnd::*)(class Pythia8::StringSystem &, class Pythia8::StringVertex &, double)) &Pythia8::StringEnd::kinematicsHadron, "C++: Pythia8::StringEnd::kinematicsHadron(class Pythia8::StringSystem &, class Pythia8::StringVertex &, double) --> class Pythia8::Vec4", pybind11::arg("system"), pybind11::arg("newVertex"), pybind11::arg("zHadIn"));
		cl.def("kinematicsHadronTmp", (class Pythia8::Vec4 (Pythia8::StringEnd::*)(class Pythia8::StringSystem, class Pythia8::Vec4, double, double)) &Pythia8::StringEnd::kinematicsHadronTmp, "C++: Pythia8::StringEnd::kinematicsHadronTmp(class Pythia8::StringSystem, class Pythia8::Vec4, double, double) --> class Pythia8::Vec4", pybind11::arg("system"), pybind11::arg("pRem"), pybind11::arg("phi"), pybind11::arg("mult"));
		cl.def("update", (void (Pythia8::StringEnd::*)()) &Pythia8::StringEnd::update, "C++: Pythia8::StringEnd::update() --> void");
		cl.def("storePrev", (void (Pythia8::StringEnd::*)()) &Pythia8::StringEnd::storePrev, "C++: Pythia8::StringEnd::storePrev() --> void");
		cl.def("updateToPrev", (void (Pythia8::StringEnd::*)()) &Pythia8::StringEnd::updateToPrev, "C++: Pythia8::StringEnd::updateToPrev() --> void");
		cl.def("assign", (class Pythia8::StringEnd & (Pythia8::StringEnd::*)(const class Pythia8::StringEnd &)) &Pythia8::StringEnd::operator=, "C++: Pythia8::StringEnd::operator=(const class Pythia8::StringEnd &) --> class Pythia8::StringEnd &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::StringFragmentation file:Pythia8/StringFragmentation.h line:106
		pybind11::class_<Pythia8::StringFragmentation, std::shared_ptr<Pythia8::StringFragmentation>, PyCallBack_Pythia8_StringFragmentation, Pythia8::FragmentationModel> cl(M("Pythia8"), "StringFragmentation", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::StringFragmentation(); }, [](){ return new PyCallBack_Pythia8_StringFragmentation(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_StringFragmentation const &o){ return new PyCallBack_Pythia8_StringFragmentation(o); } ) );
		cl.def( pybind11::init( [](Pythia8::StringFragmentation const &o){ return new Pythia8::StringFragmentation(o); } ) );
		cl.def_readwrite("flavSelNow", &Pythia8::StringFragmentation::flavSelNow);
		cl.def("init", [](Pythia8::StringFragmentation &o) -> bool { return o.init(); }, "");
		cl.def("init", [](Pythia8::StringFragmentation &o, class Pythia8::StringFlav * a0) -> bool { return o.init(a0); }, "", pybind11::arg("flavSelPtrIn"));
		cl.def("init", [](Pythia8::StringFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1) -> bool { return o.init(a0, a1); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"));
		cl.def("init", [](Pythia8::StringFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2) -> bool { return o.init(a0, a1, a2); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"));
		cl.def("init", (bool (Pythia8::StringFragmentation::*)(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>)) &Pythia8::StringFragmentation::init, "C++: Pythia8::StringFragmentation::init(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>) --> bool", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"), pybind11::arg("fragModPtrIn"));
		cl.def("fragment", [](Pythia8::StringFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2) -> bool { return o.fragment(a0, a1, a2); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"));
		cl.def("fragment", [](Pythia8::StringFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool const & a3) -> bool { return o.fragment(a0, a1, a2, a3); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"));
		cl.def("fragment", (bool (Pythia8::StringFragmentation::*)(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool)) &Pythia8::StringFragmentation::fragment, "C++: Pythia8::StringFragmentation::fragment(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool) --> bool", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"), pybind11::arg("systemRecoil"));
		cl.def("junctionRestFrame", [](Pythia8::StringFragmentation const &o, const class Pythia8::Vec4 & a0, const class Pythia8::Vec4 & a1, const class Pythia8::Vec4 & a2) -> Pythia8::Vec4 { return o.junctionRestFrame(a0, a1, a2); }, "", pybind11::arg("p0"), pybind11::arg("p1"), pybind11::arg("p2"));
		cl.def("junctionRestFrame", (class Pythia8::Vec4 (Pythia8::StringFragmentation::*)(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const bool) const) &Pythia8::StringFragmentation::junctionRestFrame, "C++: Pythia8::StringFragmentation::junctionRestFrame(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const bool) const --> class Pythia8::Vec4", pybind11::arg("p0"), pybind11::arg("p1"), pybind11::arg("p2"), pybind11::arg("angleCheck"));
		cl.def("setMVecRatio", (void (Pythia8::StringFragmentation::*)(double)) &Pythia8::StringFragmentation::setMVecRatio, "C++: Pythia8::StringFragmentation::setMVecRatio(double) --> void", pybind11::arg("mVecRatioIn"));
		cl.def("assign", (class Pythia8::StringFragmentation & (Pythia8::StringFragmentation::*)(const class Pythia8::StringFragmentation &)) &Pythia8::StringFragmentation::operator=, "C++: Pythia8::StringFragmentation::operator=(const class Pythia8::StringFragmentation &) --> class Pythia8::StringFragmentation &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
