#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/MathTools.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
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

// Pythia8::PhysicsBase file:Pythia8/PhysicsBase.h line:23
struct PyCallBack_Pythia8_PhysicsBase : public Pythia8::PhysicsBase {
	using Pythia8::PhysicsBase::PhysicsBase;

	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhysicsBase *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhysicsBase *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhysicsBase *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::PhysicsBase *>(this), "onStat");
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

void bind_Pythia8_MathTools(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Pythia8::factorial(const int) file:Pythia8/MathTools.h line:44
	M("Pythia8").def("factorial", (double (*)(const int)) &Pythia8::factorial, "C++: Pythia8::factorial(const int) --> double", pybind11::arg(""));

	// Pythia8::binomial(const int, int) file:Pythia8/MathTools.h line:47
	M("Pythia8").def("binomial", (int (*)(const int, int)) &Pythia8::binomial, "C++: Pythia8::binomial(const int, int) --> int", pybind11::arg(""), pybind11::arg(""));

	// Pythia8::lambertW(const double) file:Pythia8/MathTools.h line:50
	M("Pythia8").def("lambertW", (double (*)(const double)) &Pythia8::lambertW, "C++: Pythia8::lambertW(const double) --> double", pybind11::arg("x"));

	// Pythia8::kallenFunction(const double, const double, const double) file:Pythia8/MathTools.h line:53
	M("Pythia8").def("kallenFunction", (double (*)(const double, const double, const double)) &Pythia8::kallenFunction, "C++: Pythia8::kallenFunction(const double, const double, const double) --> double", pybind11::arg("x"), pybind11::arg("y"), pybind11::arg("z"));

	// Pythia8::linSpace(int, double, double) file:Pythia8/MathTools.h line:56
	M("Pythia8").def("linSpace", (class std::vector<double, class std::allocator<double> > (*)(int, double, double)) &Pythia8::linSpace, "C++: Pythia8::linSpace(int, double, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("nPts"), pybind11::arg("xMin"), pybind11::arg("xMax"));

	// Pythia8::logSpace(int, double, double) file:Pythia8/MathTools.h line:57
	M("Pythia8").def("logSpace", (class std::vector<double, class std::allocator<double> > (*)(int, double, double)) &Pythia8::logSpace, "C++: Pythia8::logSpace(int, double, double) --> class std::vector<double, class std::allocator<double> >", pybind11::arg("nPts"), pybind11::arg("xMin"), pybind11::arg("xMax"));

	{ // Pythia8::LinearInterpolator file:Pythia8/MathTools.h line:65
		pybind11::class_<Pythia8::LinearInterpolator, std::shared_ptr<Pythia8::LinearInterpolator>> cl(M("Pythia8"), "LinearInterpolator", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LinearInterpolator(); } ) );
		cl.def( pybind11::init<double, double, class std::vector<double, class std::allocator<double> >>(), pybind11::arg("leftIn"), pybind11::arg("rightIn"), pybind11::arg("ysIn") );

		cl.def( pybind11::init( [](Pythia8::LinearInterpolator const &o){ return new Pythia8::LinearInterpolator(o); } ) );
		cl.def("data", (const class std::vector<double, class std::allocator<double> > & (Pythia8::LinearInterpolator::*)() const) &Pythia8::LinearInterpolator::data, "C++: Pythia8::LinearInterpolator::data() const --> const class std::vector<double, class std::allocator<double> > &", pybind11::return_value_policy::reference);
		cl.def("left", (double (Pythia8::LinearInterpolator::*)() const) &Pythia8::LinearInterpolator::left, "C++: Pythia8::LinearInterpolator::left() const --> double");
		cl.def("right", (double (Pythia8::LinearInterpolator::*)() const) &Pythia8::LinearInterpolator::right, "C++: Pythia8::LinearInterpolator::right() const --> double");
		cl.def("dx", (double (Pythia8::LinearInterpolator::*)() const) &Pythia8::LinearInterpolator::dx, "C++: Pythia8::LinearInterpolator::dx() const --> double");
		cl.def("xi", (double (Pythia8::LinearInterpolator::*)(int) const) &Pythia8::LinearInterpolator::xi, "C++: Pythia8::LinearInterpolator::xi(int) const --> double", pybind11::arg("i"));
		cl.def("at", (double (Pythia8::LinearInterpolator::*)(double) const) &Pythia8::LinearInterpolator::at, "C++: Pythia8::LinearInterpolator::at(double) const --> double", pybind11::arg("x"));
		cl.def("__call__", (double (Pythia8::LinearInterpolator::*)(double) const) &Pythia8::LinearInterpolator::operator(), "C++: Pythia8::LinearInterpolator::operator()(double) const --> double", pybind11::arg("x"));
		cl.def("plot", (class Pythia8::Hist (Pythia8::LinearInterpolator::*)(std::string) const) &Pythia8::LinearInterpolator::plot, "C++: Pythia8::LinearInterpolator::plot(std::string) const --> class Pythia8::Hist", pybind11::arg("title"));
		cl.def("plot", (class Pythia8::Hist (Pythia8::LinearInterpolator::*)(std::string, double, double) const) &Pythia8::LinearInterpolator::plot, "C++: Pythia8::LinearInterpolator::plot(std::string, double, double) const --> class Pythia8::Hist", pybind11::arg("title"), pybind11::arg("xMin"), pybind11::arg("xMax"));
		cl.def("sample", (double (Pythia8::LinearInterpolator::*)(class Pythia8::Rndm &) const) &Pythia8::LinearInterpolator::sample, "C++: Pythia8::LinearInterpolator::sample(class Pythia8::Rndm &) const --> double", pybind11::arg("rndm"));
		cl.def("assign", (class Pythia8::LinearInterpolator & (Pythia8::LinearInterpolator::*)(const class Pythia8::LinearInterpolator &)) &Pythia8::LinearInterpolator::operator=, "C++: Pythia8::LinearInterpolator::operator=(const class Pythia8::LinearInterpolator &) --> class Pythia8::LinearInterpolator &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LogInterpolator file:Pythia8/MathTools.h line:109
		pybind11::class_<Pythia8::LogInterpolator, std::shared_ptr<Pythia8::LogInterpolator>> cl(M("Pythia8"), "LogInterpolator", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LogInterpolator(); } ) );
		cl.def( pybind11::init<double, double, class std::vector<double, class std::allocator<double> >>(), pybind11::arg("leftIn"), pybind11::arg("rightIn"), pybind11::arg("ysIn") );

		cl.def("data", (const class std::vector<double, class std::allocator<double> > & (Pythia8::LogInterpolator::*)() const) &Pythia8::LogInterpolator::data, "C++: Pythia8::LogInterpolator::data() const --> const class std::vector<double, class std::allocator<double> > &", pybind11::return_value_policy::reference);
		cl.def("left", (double (Pythia8::LogInterpolator::*)() const) &Pythia8::LogInterpolator::left, "C++: Pythia8::LogInterpolator::left() const --> double");
		cl.def("right", (double (Pythia8::LogInterpolator::*)() const) &Pythia8::LogInterpolator::right, "C++: Pythia8::LogInterpolator::right() const --> double");
		cl.def("rx", (double (Pythia8::LogInterpolator::*)() const) &Pythia8::LogInterpolator::rx, "C++: Pythia8::LogInterpolator::rx() const --> double");
		cl.def("at", (double (Pythia8::LogInterpolator::*)(double) const) &Pythia8::LogInterpolator::at, "C++: Pythia8::LogInterpolator::at(double) const --> double", pybind11::arg("x"));
		cl.def("__call__", (double (Pythia8::LogInterpolator::*)(double) const) &Pythia8::LogInterpolator::operator(), "C++: Pythia8::LogInterpolator::operator()(double) const --> double", pybind11::arg("x"));
		cl.def("plot", (class Pythia8::Hist (Pythia8::LogInterpolator::*)(std::string, int, double, double) const) &Pythia8::LogInterpolator::plot, "C++: Pythia8::LogInterpolator::plot(std::string, int, double, double) const --> class Pythia8::Hist", pybind11::arg("title"), pybind11::arg("nBins"), pybind11::arg("xMin"), pybind11::arg("xMax"));
	}
	{ // Pythia8::HungarianAlgorithm file:Pythia8/MathTools.h line:187
		pybind11::class_<Pythia8::HungarianAlgorithm, std::shared_ptr<Pythia8::HungarianAlgorithm>> cl(M("Pythia8"), "HungarianAlgorithm", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HungarianAlgorithm(); } ) );
		cl.def("solve", (double (Pythia8::HungarianAlgorithm::*)(class std::vector<class std::vector<double, class std::allocator<double> >, class std::allocator<class std::vector<double, class std::allocator<double> > > > &, class std::vector<int, class std::allocator<int> > &)) &Pythia8::HungarianAlgorithm::solve, "C++: Pythia8::HungarianAlgorithm::solve(class std::vector<class std::vector<double, class std::allocator<double> >, class std::allocator<class std::vector<double, class std::allocator<double> > > > &, class std::vector<int, class std::allocator<int> > &) --> double", pybind11::arg("distMatrix"), pybind11::arg("assignment"));
	}
	{ // Pythia8::PhysicsBase file:Pythia8/PhysicsBase.h line:23
		pybind11::class_<Pythia8::PhysicsBase, std::shared_ptr<Pythia8::PhysicsBase>, PyCallBack_Pythia8_PhysicsBase> cl(M("Pythia8"), "PhysicsBase", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PhysicsBase(); }, [](){ return new PyCallBack_Pythia8_PhysicsBase(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_PhysicsBase const &o){ return new PyCallBack_Pythia8_PhysicsBase(o); } ) );
		cl.def( pybind11::init( [](Pythia8::PhysicsBase const &o){ return new Pythia8::PhysicsBase(o); } ) );

		pybind11::enum_<Pythia8::PhysicsBase::Status>(cl, "Status", pybind11::arithmetic(), "")
			.value("INCOMPLETE", Pythia8::PhysicsBase::Status::INCOMPLETE)
			.value("COMPLETE", Pythia8::PhysicsBase::Status::COMPLETE)
			.value("CONSTRUCTOR_FAILED", Pythia8::PhysicsBase::Status::CONSTRUCTOR_FAILED)
			.value("INIT_FAILED", Pythia8::PhysicsBase::Status::INIT_FAILED)
			.value("LHEF_END", Pythia8::PhysicsBase::Status::LHEF_END)
			.value("LOWENERGY_FAILED", Pythia8::PhysicsBase::Status::LOWENERGY_FAILED)
			.value("PROCESSLEVEL_FAILED", Pythia8::PhysicsBase::Status::PROCESSLEVEL_FAILED)
			.value("PROCESSLEVEL_USERVETO", Pythia8::PhysicsBase::Status::PROCESSLEVEL_USERVETO)
			.value("MERGING_FAILED", Pythia8::PhysicsBase::Status::MERGING_FAILED)
			.value("PARTONLEVEL_FAILED", Pythia8::PhysicsBase::Status::PARTONLEVEL_FAILED)
			.value("PARTONLEVEL_USERVETO", Pythia8::PhysicsBase::Status::PARTONLEVEL_USERVETO)
			.value("HADRONLEVEL_FAILED", Pythia8::PhysicsBase::Status::HADRONLEVEL_FAILED)
			.value("CHECK_FAILED", Pythia8::PhysicsBase::Status::CHECK_FAILED)
			.value("OTHER_UNPHYSICAL", Pythia8::PhysicsBase::Status::OTHER_UNPHYSICAL)
			.value("HEAVYION_FAILED", Pythia8::PhysicsBase::Status::HEAVYION_FAILED)
			.value("HADRONLEVEL_USERVETO", Pythia8::PhysicsBase::Status::HADRONLEVEL_USERVETO)
			.export_values();

		cl.def_readwrite("subObjects", &Pythia8::PhysicsBase::subObjects);
		cl.def_readwrite("userHooksPtr", &Pythia8::PhysicsBase::userHooksPtr);
		cl.def("initInfoPtr", (void (Pythia8::PhysicsBase::*)(class Pythia8::Info &)) &Pythia8::PhysicsBase::initInfoPtr, "C++: Pythia8::PhysicsBase::initInfoPtr(class Pythia8::Info &) --> void", pybind11::arg("infoPtrIn"));
		cl.def("flag", (bool (Pythia8::PhysicsBase::*)(std::string) const) &Pythia8::PhysicsBase::flag, "C++: Pythia8::PhysicsBase::flag(std::string) const --> bool", pybind11::arg("key"));
		cl.def("mode", (int (Pythia8::PhysicsBase::*)(std::string) const) &Pythia8::PhysicsBase::mode, "C++: Pythia8::PhysicsBase::mode(std::string) const --> int", pybind11::arg("key"));
		cl.def("parm", (double (Pythia8::PhysicsBase::*)(std::string) const) &Pythia8::PhysicsBase::parm, "C++: Pythia8::PhysicsBase::parm(std::string) const --> double", pybind11::arg("key"));
		cl.def("word", (std::string (Pythia8::PhysicsBase::*)(std::string) const) &Pythia8::PhysicsBase::word, "C++: Pythia8::PhysicsBase::word(std::string) const --> std::string", pybind11::arg("key"));
		cl.def("fvec", (class std::vector<bool, class std::allocator<bool> > (Pythia8::PhysicsBase::*)(std::string) const) &Pythia8::PhysicsBase::fvec, "C++: Pythia8::PhysicsBase::fvec(std::string) const --> class std::vector<bool, class std::allocator<bool> >", pybind11::arg("key"));
		cl.def("mvec", (class std::vector<int, class std::allocator<int> > (Pythia8::PhysicsBase::*)(std::string) const) &Pythia8::PhysicsBase::mvec, "C++: Pythia8::PhysicsBase::mvec(std::string) const --> class std::vector<int, class std::allocator<int> >", pybind11::arg("key"));
		cl.def("pvec", (class std::vector<double, class std::allocator<double> > (Pythia8::PhysicsBase::*)(std::string) const) &Pythia8::PhysicsBase::pvec, "C++: Pythia8::PhysicsBase::pvec(std::string) const --> class std::vector<double, class std::allocator<double> >", pybind11::arg("key"));
		cl.def("wvec", (class std::vector<std::string, class std::allocator<std::string > > (Pythia8::PhysicsBase::*)(std::string) const) &Pythia8::PhysicsBase::wvec, "C++: Pythia8::PhysicsBase::wvec(std::string) const --> class std::vector<std::string, class std::allocator<std::string > >", pybind11::arg("key"));
		cl.def("onInitInfoPtr", (void (Pythia8::PhysicsBase::*)()) &Pythia8::PhysicsBase::onInitInfoPtr, "C++: Pythia8::PhysicsBase::onInitInfoPtr() --> void");
		cl.def("onBeginEvent", (void (Pythia8::PhysicsBase::*)()) &Pythia8::PhysicsBase::onBeginEvent, "C++: Pythia8::PhysicsBase::onBeginEvent() --> void");
		cl.def("onEndEvent", (void (Pythia8::PhysicsBase::*)(enum Pythia8::PhysicsBase::Status)) &Pythia8::PhysicsBase::onEndEvent, "C++: Pythia8::PhysicsBase::onEndEvent(enum Pythia8::PhysicsBase::Status) --> void", pybind11::arg(""));
		cl.def("onStat", (void (Pythia8::PhysicsBase::*)()) &Pythia8::PhysicsBase::onStat, "C++: Pythia8::PhysicsBase::onStat() --> void");
		cl.def("registerSubObject", (void (Pythia8::PhysicsBase::*)(class Pythia8::PhysicsBase &)) &Pythia8::PhysicsBase::registerSubObject, "C++: Pythia8::PhysicsBase::registerSubObject(class Pythia8::PhysicsBase &) --> void", pybind11::arg("pb"));
		cl.def("assign", (class Pythia8::PhysicsBase & (Pythia8::PhysicsBase::*)(const class Pythia8::PhysicsBase &)) &Pythia8::PhysicsBase::operator=, "C++: Pythia8::PhysicsBase::operator=(const class Pythia8::PhysicsBase &) --> class Pythia8::PhysicsBase &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::FlavContainer file:Pythia8/FragmentationFlavZpT.h line:33
		pybind11::class_<Pythia8::FlavContainer, std::shared_ptr<Pythia8::FlavContainer>> cl(M("Pythia8"), "FlavContainer", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::FlavContainer(); } ), "doc" );
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::FlavContainer(a0); } ), "doc" , pybind11::arg("idIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Pythia8::FlavContainer(a0, a1); } ), "doc" , pybind11::arg("idIn"), pybind11::arg("rankIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2){ return new Pythia8::FlavContainer(a0, a1, a2); } ), "doc" , pybind11::arg("idIn"), pybind11::arg("rankIn"), pybind11::arg("nPopIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1, int const & a2, int const & a3){ return new Pythia8::FlavContainer(a0, a1, a2, a3); } ), "doc" , pybind11::arg("idIn"), pybind11::arg("rankIn"), pybind11::arg("nPopIn"), pybind11::arg("idPopIn"));
		cl.def( pybind11::init<int, int, int, int, int>(), pybind11::arg("idIn"), pybind11::arg("rankIn"), pybind11::arg("nPopIn"), pybind11::arg("idPopIn"), pybind11::arg("idVtxIn") );

		cl.def( pybind11::init( [](Pythia8::FlavContainer const &o){ return new Pythia8::FlavContainer(o); } ) );
		cl.def_readwrite("id", &Pythia8::FlavContainer::id);
		cl.def_readwrite("rank", &Pythia8::FlavContainer::rank);
		cl.def_readwrite("nPop", &Pythia8::FlavContainer::nPop);
		cl.def_readwrite("idPop", &Pythia8::FlavContainer::idPop);
		cl.def_readwrite("idVtx", &Pythia8::FlavContainer::idVtx);
		cl.def("assign", (class Pythia8::FlavContainer & (Pythia8::FlavContainer::*)(const class Pythia8::FlavContainer &)) &Pythia8::FlavContainer::operator=, "C++: Pythia8::FlavContainer::operator=(const class Pythia8::FlavContainer &) --> class Pythia8::FlavContainer &", pybind11::return_value_policy::reference, pybind11::arg("flav"));
		cl.def("anti", (class Pythia8::FlavContainer & (Pythia8::FlavContainer::*)()) &Pythia8::FlavContainer::anti, "C++: Pythia8::FlavContainer::anti() --> class Pythia8::FlavContainer &", pybind11::return_value_policy::reference);
		cl.def("copy", (class Pythia8::FlavContainer & (Pythia8::FlavContainer::*)(const class Pythia8::FlavContainer &)) &Pythia8::FlavContainer::copy, "C++: Pythia8::FlavContainer::copy(const class Pythia8::FlavContainer &) --> class Pythia8::FlavContainer &", pybind11::return_value_policy::reference, pybind11::arg("flav"));
		cl.def("anti", (class Pythia8::FlavContainer & (Pythia8::FlavContainer::*)(const class Pythia8::FlavContainer &)) &Pythia8::FlavContainer::anti, "C++: Pythia8::FlavContainer::anti(const class Pythia8::FlavContainer &) --> class Pythia8::FlavContainer &", pybind11::return_value_policy::reference, pybind11::arg("flav"));
		cl.def("isDiquark", (bool (Pythia8::FlavContainer::*)()) &Pythia8::FlavContainer::isDiquark, "C++: Pythia8::FlavContainer::isDiquark() --> bool");
	}
}
