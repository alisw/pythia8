#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/GammaKinematics.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/ProcessContainer.h>
#include <Pythia8/ProcessLevel.h>
#include <Pythia8/ResonanceDecays.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/SLHAinterface.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaProcess.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
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

// Pythia8::ProcessContainer file:Pythia8/ProcessContainer.h line:39
struct PyCallBack_Pythia8_ProcessContainer : public Pythia8::ProcessContainer {
	using Pythia8::ProcessContainer::ProcessContainer;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ProcessContainer *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ProcessContainer *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ProcessContainer *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ProcessContainer *>(this), "onStat");
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

// Pythia8::ProcessLevel file:Pythia8/ProcessLevel.h line:36
struct PyCallBack_Pythia8_ProcessLevel : public Pythia8::ProcessLevel {
	using Pythia8::ProcessLevel::ProcessLevel;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ProcessLevel *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ProcessLevel::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ProcessLevel *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ProcessLevel *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ProcessLevel *>(this), "onStat");
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

void bind_Pythia8_ProcessContainer(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::ProcessContainer file:Pythia8/ProcessContainer.h line:39
		pybind11::class_<Pythia8::ProcessContainer, std::shared_ptr<Pythia8::ProcessContainer>, PyCallBack_Pythia8_ProcessContainer, Pythia8::PhysicsBase> cl(M("Pythia8"), "ProcessContainer", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ProcessContainer(); }, [](){ return new PyCallBack_Pythia8_ProcessContainer(); } ), "doc");
		cl.def( pybind11::init( [](class std::shared_ptr<class Pythia8::SigmaProcess> const & a0){ return new Pythia8::ProcessContainer(a0); }, [](class std::shared_ptr<class Pythia8::SigmaProcess> const & a0){ return new PyCallBack_Pythia8_ProcessContainer(a0); } ), "doc");
		cl.def( pybind11::init<class std::shared_ptr<class Pythia8::SigmaProcess>, class std::shared_ptr<class Pythia8::PhaseSpace>>(), pybind11::arg("sigmaProcessPtrIn"), pybind11::arg("phaseSpacePtrIn") );

		cl.def( pybind11::init( [](PyCallBack_Pythia8_ProcessContainer const &o){ return new PyCallBack_Pythia8_ProcessContainer(o); } ) );
		cl.def( pybind11::init( [](Pythia8::ProcessContainer const &o){ return new Pythia8::ProcessContainer(o); } ) );
		cl.def("init", (bool (Pythia8::ProcessContainer::*)(bool, class Pythia8::ResonanceDecays *, class Pythia8::SLHAinterface *, class Pythia8::GammaKinematics *)) &Pythia8::ProcessContainer::init, "C++: Pythia8::ProcessContainer::init(bool, class Pythia8::ResonanceDecays *, class Pythia8::SLHAinterface *, class Pythia8::GammaKinematics *) --> bool", pybind11::arg("isFirst"), pybind11::arg("resDecaysPtrIn"), pybind11::arg("slhaInterfacePtr"), pybind11::arg("gammaKinPtrIn"));
		cl.def("setLHAPtr", [](Pythia8::ProcessContainer &o, class std::shared_ptr<class Pythia8::LHAup> const & a0) -> void { return o.setLHAPtr(a0); }, "", pybind11::arg("lhaUpPtrIn"));
		cl.def("setLHAPtr", [](Pythia8::ProcessContainer &o, class std::shared_ptr<class Pythia8::LHAup> const & a0, class Pythia8::ParticleData * a1) -> void { return o.setLHAPtr(a0, a1); }, "", pybind11::arg("lhaUpPtrIn"), pybind11::arg("particleDataPtrIn"));
		cl.def("setLHAPtr", [](Pythia8::ProcessContainer &o, class std::shared_ptr<class Pythia8::LHAup> const & a0, class Pythia8::ParticleData * a1, class Pythia8::Settings * a2) -> void { return o.setLHAPtr(a0, a1, a2); }, "", pybind11::arg("lhaUpPtrIn"), pybind11::arg("particleDataPtrIn"), pybind11::arg("settingsPtrIn"));
		cl.def("setLHAPtr", (void (Pythia8::ProcessContainer::*)(class std::shared_ptr<class Pythia8::LHAup>, class Pythia8::ParticleData *, class Pythia8::Settings *, class Pythia8::Rndm *)) &Pythia8::ProcessContainer::setLHAPtr, "C++: Pythia8::ProcessContainer::setLHAPtr(class std::shared_ptr<class Pythia8::LHAup>, class Pythia8::ParticleData *, class Pythia8::Settings *, class Pythia8::Rndm *) --> void", pybind11::arg("lhaUpPtrIn"), pybind11::arg("particleDataPtrIn"), pybind11::arg("settingsPtrIn"), pybind11::arg("rndmPtrIn"));
		cl.def("updateBeamIDs", (void (Pythia8::ProcessContainer::*)()) &Pythia8::ProcessContainer::updateBeamIDs, "C++: Pythia8::ProcessContainer::updateBeamIDs() --> void");
		cl.def("newECM", (void (Pythia8::ProcessContainer::*)(double)) &Pythia8::ProcessContainer::newECM, "C++: Pythia8::ProcessContainer::newECM(double) --> void", pybind11::arg("eCM"));
		cl.def("trialProcess", (bool (Pythia8::ProcessContainer::*)()) &Pythia8::ProcessContainer::trialProcess, "C++: Pythia8::ProcessContainer::trialProcess() --> bool");
		cl.def("constructState", (bool (Pythia8::ProcessContainer::*)()) &Pythia8::ProcessContainer::constructState, "C++: Pythia8::ProcessContainer::constructState() --> bool");
		cl.def("constructProcess", [](Pythia8::ProcessContainer &o, class Pythia8::Event & a0) -> bool { return o.constructProcess(a0); }, "", pybind11::arg("process"));
		cl.def("constructProcess", (bool (Pythia8::ProcessContainer::*)(class Pythia8::Event &, bool)) &Pythia8::ProcessContainer::constructProcess, "C++: Pythia8::ProcessContainer::constructProcess(class Pythia8::Event &, bool) --> bool", pybind11::arg("process"), pybind11::arg("isHardest"));
		cl.def("constructDecays", (bool (Pythia8::ProcessContainer::*)(class Pythia8::Event &)) &Pythia8::ProcessContainer::constructDecays, "C++: Pythia8::ProcessContainer::constructDecays(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("decayResonances", (bool (Pythia8::ProcessContainer::*)(class Pythia8::Event &)) &Pythia8::ProcessContainer::decayResonances, "C++: Pythia8::ProcessContainer::decayResonances(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("accumulate", (void (Pythia8::ProcessContainer::*)()) &Pythia8::ProcessContainer::accumulate, "C++: Pythia8::ProcessContainer::accumulate() --> void");
		cl.def("reset", (void (Pythia8::ProcessContainer::*)()) &Pythia8::ProcessContainer::reset, "C++: Pythia8::ProcessContainer::reset() --> void");
		cl.def("setBeamModes", [](Pythia8::ProcessContainer &o) -> void { return o.setBeamModes(); }, "");
		cl.def("setBeamModes", [](Pythia8::ProcessContainer &o, bool const & a0) -> void { return o.setBeamModes(a0); }, "", pybind11::arg("setVMD"));
		cl.def("setBeamModes", (void (Pythia8::ProcessContainer::*)(bool, bool)) &Pythia8::ProcessContainer::setBeamModes, "C++: Pythia8::ProcessContainer::setBeamModes(bool, bool) --> void", pybind11::arg("setVMD"), pybind11::arg("isSampled"));
		cl.def("name", (std::string (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::name, "C++: Pythia8::ProcessContainer::name() const --> std::string");
		cl.def("code", (int (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::code, "C++: Pythia8::ProcessContainer::code() const --> int");
		cl.def("nFinal", (int (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::nFinal, "C++: Pythia8::ProcessContainer::nFinal() const --> int");
		cl.def("isSUSY", (bool (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::isSUSY, "C++: Pythia8::ProcessContainer::isSUSY() const --> bool");
		cl.def("isNonDiffractive", (bool (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::isNonDiffractive, "C++: Pythia8::ProcessContainer::isNonDiffractive() const --> bool");
		cl.def("isSoftQCD", (bool (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::isSoftQCD, "C++: Pythia8::ProcessContainer::isSoftQCD() const --> bool");
		cl.def("isElastic", (bool (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::isElastic, "C++: Pythia8::ProcessContainer::isElastic() const --> bool");
		cl.def("newSigmaMax", (bool (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::newSigmaMax, "C++: Pythia8::ProcessContainer::newSigmaMax() const --> bool");
		cl.def("sigmaMax", (double (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::sigmaMax, "C++: Pythia8::ProcessContainer::sigmaMax() const --> double");
		cl.def("nTried", (long (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::nTried, "C++: Pythia8::ProcessContainer::nTried() const --> long");
		cl.def("nSelected", (long (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::nSelected, "C++: Pythia8::ProcessContainer::nSelected() const --> long");
		cl.def("nAccepted", (long (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::nAccepted, "C++: Pythia8::ProcessContainer::nAccepted() const --> long");
		cl.def("weightSum", (double (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::weightSum, "C++: Pythia8::ProcessContainer::weightSum() const --> double");
		cl.def("sigmaSelMC", [](Pythia8::ProcessContainer &o) -> double { return o.sigmaSelMC(); }, "");
		cl.def("sigmaSelMC", (double (Pythia8::ProcessContainer::*)(bool)) &Pythia8::ProcessContainer::sigmaSelMC, "C++: Pythia8::ProcessContainer::sigmaSelMC(bool) --> double", pybind11::arg("doAccumulate"));
		cl.def("sigmaMC", [](Pythia8::ProcessContainer &o) -> double { return o.sigmaMC(); }, "");
		cl.def("sigmaMC", (double (Pythia8::ProcessContainer::*)(bool)) &Pythia8::ProcessContainer::sigmaMC, "C++: Pythia8::ProcessContainer::sigmaMC(bool) --> double", pybind11::arg("doAccumulate"));
		cl.def("deltaMC", [](Pythia8::ProcessContainer &o) -> double { return o.deltaMC(); }, "");
		cl.def("deltaMC", (double (Pythia8::ProcessContainer::*)(bool)) &Pythia8::ProcessContainer::deltaMC, "C++: Pythia8::ProcessContainer::deltaMC(bool) --> double", pybind11::arg("doAccumulate"));
		cl.def("sigmaMaxSwitch", (double (Pythia8::ProcessContainer::*)()) &Pythia8::ProcessContainer::sigmaMaxSwitch, "C++: Pythia8::ProcessContainer::sigmaMaxSwitch() --> double");
		cl.def("id1", (int (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::id1, "C++: Pythia8::ProcessContainer::id1() const --> int");
		cl.def("id2", (int (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::id2, "C++: Pythia8::ProcessContainer::id2() const --> int");
		cl.def("x1", (double (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::x1, "C++: Pythia8::ProcessContainer::x1() const --> double");
		cl.def("x2", (double (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::x2, "C++: Pythia8::ProcessContainer::x2() const --> double");
		cl.def("Q2Fac", (double (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::Q2Fac, "C++: Pythia8::ProcessContainer::Q2Fac() const --> double");
		cl.def("mHat", (double (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::mHat, "C++: Pythia8::ProcessContainer::mHat() const --> double");
		cl.def("pTHat", (double (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::pTHat, "C++: Pythia8::ProcessContainer::pTHat() const --> double");
		cl.def("isLHAContainer", (bool (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::isLHAContainer, "C++: Pythia8::ProcessContainer::isLHAContainer() const --> bool");
		cl.def("lhaStrategy", (int (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::lhaStrategy, "C++: Pythia8::ProcessContainer::lhaStrategy() const --> int");
		cl.def("codeLHASize", (int (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::codeLHASize, "C++: Pythia8::ProcessContainer::codeLHASize() const --> int");
		cl.def("subCodeLHA", (int (Pythia8::ProcessContainer::*)(int) const) &Pythia8::ProcessContainer::subCodeLHA, "C++: Pythia8::ProcessContainer::subCodeLHA(int) const --> int", pybind11::arg("i"));
		cl.def("nTriedLHA", (long (Pythia8::ProcessContainer::*)(int) const) &Pythia8::ProcessContainer::nTriedLHA, "C++: Pythia8::ProcessContainer::nTriedLHA(int) const --> long", pybind11::arg("i"));
		cl.def("nSelectedLHA", (long (Pythia8::ProcessContainer::*)(int) const) &Pythia8::ProcessContainer::nSelectedLHA, "C++: Pythia8::ProcessContainer::nSelectedLHA(int) const --> long", pybind11::arg("i"));
		cl.def("nAcceptedLHA", (long (Pythia8::ProcessContainer::*)(int) const) &Pythia8::ProcessContainer::nAcceptedLHA, "C++: Pythia8::ProcessContainer::nAcceptedLHA(int) const --> long", pybind11::arg("i"));
		cl.def("isSame", (void (Pythia8::ProcessContainer::*)(bool)) &Pythia8::ProcessContainer::isSame, "C++: Pythia8::ProcessContainer::isSame(bool) --> void", pybind11::arg("isSameIn"));
		cl.def("isSame", (bool (Pythia8::ProcessContainer::*)() const) &Pythia8::ProcessContainer::isSame, "C++: Pythia8::ProcessContainer::isSame() const --> bool");
		cl.def("assign", (class Pythia8::ProcessContainer & (Pythia8::ProcessContainer::*)(const class Pythia8::ProcessContainer &)) &Pythia8::ProcessContainer::operator=, "C++: Pythia8::ProcessContainer::operator=(const class Pythia8::ProcessContainer &) --> class Pythia8::ProcessContainer &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SetupContainers file:Pythia8/ProcessContainer.h line:233
		pybind11::class_<Pythia8::SetupContainers, std::shared_ptr<Pythia8::SetupContainers>> cl(M("Pythia8"), "SetupContainers", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SetupContainers(); } ) );
		cl.def("init", (bool (Pythia8::SetupContainers::*)(class std::vector<class Pythia8::ProcessContainer *, class std::allocator<class Pythia8::ProcessContainer *> > &, class Pythia8::Info *)) &Pythia8::SetupContainers::init, "C++: Pythia8::SetupContainers::init(class std::vector<class Pythia8::ProcessContainer *, class std::allocator<class Pythia8::ProcessContainer *> > &, class Pythia8::Info *) --> bool", pybind11::arg("containerPtrs"), pybind11::arg("infoPtr"));
		cl.def("init2", (bool (Pythia8::SetupContainers::*)(class std::vector<class Pythia8::ProcessContainer *, class std::allocator<class Pythia8::ProcessContainer *> > &, class Pythia8::Info *)) &Pythia8::SetupContainers::init2, "C++: Pythia8::SetupContainers::init2(class std::vector<class Pythia8::ProcessContainer *, class std::allocator<class Pythia8::ProcessContainer *> > &, class Pythia8::Info *) --> bool", pybind11::arg("container2Ptrs"), pybind11::arg("infoPtr"));
	}
	{ // Pythia8::ProcessLevel file:Pythia8/ProcessLevel.h line:36
		pybind11::class_<Pythia8::ProcessLevel, std::shared_ptr<Pythia8::ProcessLevel>, PyCallBack_Pythia8_ProcessLevel, Pythia8::PhysicsBase> cl(M("Pythia8"), "ProcessLevel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ProcessLevel(); }, [](){ return new PyCallBack_Pythia8_ProcessLevel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_ProcessLevel const &o){ return new PyCallBack_Pythia8_ProcessLevel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::ProcessLevel const &o){ return new Pythia8::ProcessLevel(o); } ) );
		cl.def("init", (bool (Pythia8::ProcessLevel::*)(bool, class Pythia8::SLHAinterface *, class std::vector<class std::shared_ptr<class Pythia8::SigmaProcess>, class std::allocator<class std::shared_ptr<class Pythia8::SigmaProcess> > > &, class std::vector<class std::shared_ptr<class Pythia8::PhaseSpace>, class std::allocator<class std::shared_ptr<class Pythia8::PhaseSpace> > > &)) &Pythia8::ProcessLevel::init, "C++: Pythia8::ProcessLevel::init(bool, class Pythia8::SLHAinterface *, class std::vector<class std::shared_ptr<class Pythia8::SigmaProcess>, class std::allocator<class std::shared_ptr<class Pythia8::SigmaProcess> > > &, class std::vector<class std::shared_ptr<class Pythia8::PhaseSpace>, class std::allocator<class std::shared_ptr<class Pythia8::PhaseSpace> > > &) --> bool", pybind11::arg("doLHAin"), pybind11::arg("slhaInterfacePtrIn"), pybind11::arg("sigmaPtrs"), pybind11::arg("phaseSpacePtrs"));
		cl.def("setLHAPtr", (void (Pythia8::ProcessLevel::*)(class std::shared_ptr<class Pythia8::LHAup>)) &Pythia8::ProcessLevel::setLHAPtr, "C++: Pythia8::ProcessLevel::setLHAPtr(class std::shared_ptr<class Pythia8::LHAup>) --> void", pybind11::arg("lhaUpPtrIn"));
		cl.def("updateBeamIDs", (void (Pythia8::ProcessLevel::*)()) &Pythia8::ProcessLevel::updateBeamIDs, "C++: Pythia8::ProcessLevel::updateBeamIDs() --> void");
		cl.def("next", [](Pythia8::ProcessLevel &o, class Pythia8::Event & a0) -> bool { return o.next(a0); }, "", pybind11::arg("process"));
		cl.def("next", (bool (Pythia8::ProcessLevel::*)(class Pythia8::Event &, int)) &Pythia8::ProcessLevel::next, "C++: Pythia8::ProcessLevel::next(class Pythia8::Event &, int) --> bool", pybind11::arg("process"), pybind11::arg("procTypeIn"));
		cl.def("nextLHAdec", (bool (Pythia8::ProcessLevel::*)(class Pythia8::Event &)) &Pythia8::ProcessLevel::nextLHAdec, "C++: Pythia8::ProcessLevel::nextLHAdec(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("accumulate", [](Pythia8::ProcessLevel &o) -> void { return o.accumulate(); }, "");
		cl.def("accumulate", (void (Pythia8::ProcessLevel::*)(bool)) &Pythia8::ProcessLevel::accumulate, "C++: Pythia8::ProcessLevel::accumulate(bool) --> void", pybind11::arg("doAccumulate"));
		cl.def("statistics", [](Pythia8::ProcessLevel &o) -> void { return o.statistics(); }, "");
		cl.def("statistics", (void (Pythia8::ProcessLevel::*)(bool)) &Pythia8::ProcessLevel::statistics, "C++: Pythia8::ProcessLevel::statistics(bool) --> void", pybind11::arg("reset"));
		cl.def("resetStatistics", (void (Pythia8::ProcessLevel::*)()) &Pythia8::ProcessLevel::resetStatistics, "C++: Pythia8::ProcessLevel::resetStatistics() --> void");
		cl.def("findJunctions", (void (Pythia8::ProcessLevel::*)(class Pythia8::Event &)) &Pythia8::ProcessLevel::findJunctions, "C++: Pythia8::ProcessLevel::findJunctions(class Pythia8::Event &) --> void", pybind11::arg("junEvent"));
		cl.def("initDecays", (void (Pythia8::ProcessLevel::*)(class std::shared_ptr<class Pythia8::LHAup>)) &Pythia8::ProcessLevel::initDecays, "C++: Pythia8::ProcessLevel::initDecays(class std::shared_ptr<class Pythia8::LHAup>) --> void", pybind11::arg("lhaUpPtrIn"));
		cl.def("nextDecays", (bool (Pythia8::ProcessLevel::*)(class Pythia8::Event &)) &Pythia8::ProcessLevel::nextDecays, "C++: Pythia8::ProcessLevel::nextDecays(class Pythia8::Event &) --> bool", pybind11::arg("process"));
		cl.def("onInitInfoPtr", (void (Pythia8::ProcessLevel::*)()) &Pythia8::ProcessLevel::onInitInfoPtr, "C++: Pythia8::ProcessLevel::onInitInfoPtr() --> void");
		cl.def("assign", (class Pythia8::ProcessLevel & (Pythia8::ProcessLevel::*)(const class Pythia8::ProcessLevel &)) &Pythia8::ProcessLevel::operator=, "C++: Pythia8::ProcessLevel::operator=(const class Pythia8::ProcessLevel &) --> class Pythia8::ProcessLevel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
