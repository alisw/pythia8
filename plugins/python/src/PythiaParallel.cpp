#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/HIBasics.h>
#include <Pythia8/HINucleusModel.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PythiaParallel.h>
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

// Pythia8::NucleusModel file:Pythia8/HINucleusModel.h line:187
struct PyCallBack_Pythia8_NucleusModel : public Pythia8::NucleusModel {
	using Pythia8::NucleusModel::NucleusModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleusModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::init();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleusModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::initGeometry();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleusModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleusModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleusModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NucleusModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"NucleusModel::generate\"");
	}
};

void bind_Pythia8_PythiaParallel(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::PythiaParallel file:Pythia8/PythiaParallel.h line:18
		pybind11::class_<Pythia8::PythiaParallel, std::shared_ptr<Pythia8::PythiaParallel>> cl(M("Pythia8"), "PythiaParallel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::PythiaParallel(); } ), "doc" );
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0){ return new Pythia8::PythiaParallel(a0); } ), "doc" , pybind11::arg("xmlDir"));
		cl.def( pybind11::init<std::string, bool>(), pybind11::arg("xmlDir"), pybind11::arg("printBanner") );

		cl.def("readString", [](Pythia8::PythiaParallel &o, class std::basic_string<char> const & a0) -> bool { return o.readString(a0); }, "", pybind11::arg("setting"));
		cl.def("readString", (bool (Pythia8::PythiaParallel::*)(std::string, bool)) &Pythia8::PythiaParallel::readString, "C++: Pythia8::PythiaParallel::readString(std::string, bool) --> bool", pybind11::arg("setting"), pybind11::arg("warn"));
		cl.def("readFile", [](Pythia8::PythiaParallel &o, class std::basic_string<char> const & a0) -> bool { return o.readFile(a0); }, "", pybind11::arg("fileName"));
		cl.def("readFile", [](Pythia8::PythiaParallel &o, class std::basic_string<char> const & a0, bool const & a1) -> bool { return o.readFile(a0, a1); }, "", pybind11::arg("fileName"), pybind11::arg("warn"));
		cl.def("readFile", (bool (Pythia8::PythiaParallel::*)(std::string, bool, int)) &Pythia8::PythiaParallel::readFile, "C++: Pythia8::PythiaParallel::readFile(std::string, bool, int) --> bool", pybind11::arg("fileName"), pybind11::arg("warn"), pybind11::arg("subrun"));
		cl.def("readFile", (bool (Pythia8::PythiaParallel::*)(std::string, int)) &Pythia8::PythiaParallel::readFile, "C++: Pythia8::PythiaParallel::readFile(std::string, int) --> bool", pybind11::arg("fileName"), pybind11::arg("subrun"));
		cl.def("readFile", [](Pythia8::PythiaParallel &o) -> bool { return o.readFile(); }, "");
		cl.def("readFile", [](Pythia8::PythiaParallel &o, class std::basic_istream<char> & a0) -> bool { return o.readFile(a0); }, "", pybind11::arg("is"));
		cl.def("readFile", [](Pythia8::PythiaParallel &o, class std::basic_istream<char> & a0, bool const & a1) -> bool { return o.readFile(a0, a1); }, "", pybind11::arg("is"), pybind11::arg("warn"));
		cl.def("readFile", (bool (Pythia8::PythiaParallel::*)(class std::basic_istream<char> &, bool, int)) &Pythia8::PythiaParallel::readFile, "C++: Pythia8::PythiaParallel::readFile(class std::basic_istream<char> &, bool, int) --> bool", pybind11::arg("is"), pybind11::arg("warn"), pybind11::arg("subrun"));
		cl.def("readFile", (bool (Pythia8::PythiaParallel::*)(class std::basic_istream<char> &, int)) &Pythia8::PythiaParallel::readFile, "C++: Pythia8::PythiaParallel::readFile(class std::basic_istream<char> &, int) --> bool", pybind11::arg("is"), pybind11::arg("subrun"));
		cl.def("init", (bool (Pythia8::PythiaParallel::*)()) &Pythia8::PythiaParallel::init, pybind11::call_guard<pybind11::gil_scoped_release>(), "C++: Pythia8::PythiaParallel::init() --> bool");
		cl.def("init", (bool (Pythia8::PythiaParallel::*)(class std::function<bool (class Pythia8::Pythia *)>)) &Pythia8::PythiaParallel::init, pybind11::call_guard<pybind11::gil_scoped_release>(), "C++: Pythia8::PythiaParallel::init(class std::function<bool (class Pythia8::Pythia *)>) --> bool", pybind11::arg("additionalSetup"));
		cl.def("foreach", (void (Pythia8::PythiaParallel::*)(class std::function<void (class Pythia8::Pythia *)>)) &Pythia8::PythiaParallel::foreach, "C++: Pythia8::PythiaParallel::foreach(class std::function<void (class Pythia8::Pythia *)>) --> void", pybind11::arg("action"));
		cl.def("foreachAsync", (void (Pythia8::PythiaParallel::*)(class std::function<void (class Pythia8::Pythia *)>)) &Pythia8::PythiaParallel::foreachAsync, "C++: Pythia8::PythiaParallel::foreachAsync(class std::function<void (class Pythia8::Pythia *)>) --> void", pybind11::arg("action"));
		cl.def("stat", (void (Pythia8::PythiaParallel::*)()) &Pythia8::PythiaParallel::stat, "C++: Pythia8::PythiaParallel::stat() --> void");
		cl.def("run", (class std::vector<long, class std::allocator<long> > (Pythia8::PythiaParallel::*)(long, class std::function<void (class Pythia8::Pythia *)>)) &Pythia8::PythiaParallel::run, pybind11::call_guard<pybind11::gil_scoped_release>(), "C++: Pythia8::PythiaParallel::run(long, class std::function<void (class Pythia8::Pythia *)>) --> class std::vector<long, class std::allocator<long> >", pybind11::arg("nEvents"), pybind11::arg("callback"));
		cl.def("run", (class std::vector<long, class std::allocator<long> > (Pythia8::PythiaParallel::*)(class std::function<void (class Pythia8::Pythia *)>)) &Pythia8::PythiaParallel::run, pybind11::call_guard<pybind11::gil_scoped_release>(), "C++: Pythia8::PythiaParallel::run(class std::function<void (class Pythia8::Pythia *)>) --> class std::vector<long, class std::allocator<long> >", pybind11::arg("callback"));
		cl.def("sigmaGen", (double (Pythia8::PythiaParallel::*)() const) &Pythia8::PythiaParallel::sigmaGen, "C++: Pythia8::PythiaParallel::sigmaGen() const --> double");
		cl.def("weightSum", (double (Pythia8::PythiaParallel::*)() const) &Pythia8::PythiaParallel::weightSum, "C++: Pythia8::PythiaParallel::weightSum() const --> double");
	}
	{ // Pythia8::EventInfo file:Pythia8/HIBasics.h line:25
		pybind11::class_<Pythia8::EventInfo, std::shared_ptr<Pythia8::EventInfo>> cl(M("Pythia8"), "EventInfo", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::EventInfo(); } ) );
		cl.def( pybind11::init( [](Pythia8::EventInfo const &o){ return new Pythia8::EventInfo(o); } ) );
		cl.def_readwrite("event", &Pythia8::EventInfo::event);
		cl.def_readwrite("info", &Pythia8::EventInfo::info);
		cl.def_readwrite("code", &Pythia8::EventInfo::code);
		cl.def_readwrite("ordering", &Pythia8::EventInfo::ordering);
		cl.def_readwrite("ok", &Pythia8::EventInfo::ok);
		cl.def_readwrite("projs", &Pythia8::EventInfo::projs);
		cl.def_readwrite("targs", &Pythia8::EventInfo::targs);
	}
	{ // Pythia8::Nucleon file:Pythia8/HINucleusModel.h line:28
		pybind11::class_<Pythia8::Nucleon, std::shared_ptr<Pythia8::Nucleon>> cl(M("Pythia8"), "Nucleon", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Nucleon(); } ), "doc" );
		cl.def( pybind11::init( [](int const & a0){ return new Pythia8::Nucleon(a0); } ), "doc" , pybind11::arg("idIn"));
		cl.def( pybind11::init( [](int const & a0, int const & a1){ return new Pythia8::Nucleon(a0, a1); } ), "doc" , pybind11::arg("idIn"), pybind11::arg("indexIn"));
		cl.def( pybind11::init<int, int, const class Pythia8::Vec4 &>(), pybind11::arg("idIn"), pybind11::arg("indexIn"), pybind11::arg("pos") );

		cl.def( pybind11::init( [](Pythia8::Nucleon const &o){ return new Pythia8::Nucleon(o); } ) );

		pybind11::enum_<Pythia8::Nucleon::Status>(cl, "Status", pybind11::arithmetic(), "")
			.value("UNWOUNDED", Pythia8::Nucleon::Status::UNWOUNDED)
			.value("ELASTIC", Pythia8::Nucleon::Status::ELASTIC)
			.value("DIFF", Pythia8::Nucleon::Status::DIFF)
			.value("ABS", Pythia8::Nucleon::Status::ABS)
			.export_values();

		cl.def("id", (int (Pythia8::Nucleon::*)() const) &Pythia8::Nucleon::id, "C++: Pythia8::Nucleon::id() const --> int");
		cl.def("index", (int (Pythia8::Nucleon::*)() const) &Pythia8::Nucleon::index, "C++: Pythia8::Nucleon::index() const --> int");
		cl.def("nPos", (const class Pythia8::Vec4 & (Pythia8::Nucleon::*)() const) &Pythia8::Nucleon::nPos, "C++: Pythia8::Nucleon::nPos() const --> const class Pythia8::Vec4 &", pybind11::return_value_policy::reference);
		cl.def("bPos", (const class Pythia8::Vec4 & (Pythia8::Nucleon::*)() const) &Pythia8::Nucleon::bPos, "C++: Pythia8::Nucleon::bPos() const --> const class Pythia8::Vec4 &", pybind11::return_value_policy::reference);
		cl.def("bShift", (void (Pythia8::Nucleon::*)(const class Pythia8::Vec4 &)) &Pythia8::Nucleon::bShift, "C++: Pythia8::Nucleon::bShift(const class Pythia8::Vec4 &) --> void", pybind11::arg("bvec"));
		cl.def("status", (enum Pythia8::Nucleon::Status (Pythia8::Nucleon::*)() const) &Pythia8::Nucleon::status, "C++: Pythia8::Nucleon::status() const --> enum Pythia8::Nucleon::Status");
		cl.def("done", (bool (Pythia8::Nucleon::*)() const) &Pythia8::Nucleon::done, "C++: Pythia8::Nucleon::done() const --> bool");
		cl.def("event", (class Pythia8::EventInfo * (Pythia8::Nucleon::*)() const) &Pythia8::Nucleon::event, "C++: Pythia8::Nucleon::event() const --> class Pythia8::EventInfo *", pybind11::return_value_policy::automatic);
		cl.def("state", (const class std::vector<double, class std::allocator<double> > & (Pythia8::Nucleon::*)() const) &Pythia8::Nucleon::state, "C++: Pythia8::Nucleon::state() const --> const class std::vector<double, class std::allocator<double> > &", pybind11::return_value_policy::reference);
		cl.def("altState", [](Pythia8::Nucleon &o) -> const std::vector<double, class std::allocator<double> > & { return o.altState(); }, "", pybind11::return_value_policy::reference);
		cl.def("altState", (const class std::vector<double, class std::allocator<double> > & (Pythia8::Nucleon::*)(int)) &Pythia8::Nucleon::altState, "C++: Pythia8::Nucleon::altState(int) --> const class std::vector<double, class std::allocator<double> > &", pybind11::return_value_policy::reference, pybind11::arg("i"));
		cl.def("status", (void (Pythia8::Nucleon::*)(enum Pythia8::Nucleon::Status)) &Pythia8::Nucleon::status, "C++: Pythia8::Nucleon::status(enum Pythia8::Nucleon::Status) --> void", pybind11::arg("s"));
		cl.def("state", (void (Pythia8::Nucleon::*)(class std::vector<double, class std::allocator<double> >)) &Pythia8::Nucleon::state, "C++: Pythia8::Nucleon::state(class std::vector<double, class std::allocator<double> >) --> void", pybind11::arg("s"));
		cl.def("addAltState", (void (Pythia8::Nucleon::*)(class std::vector<double, class std::allocator<double> >)) &Pythia8::Nucleon::addAltState, "C++: Pythia8::Nucleon::addAltState(class std::vector<double, class std::allocator<double> >) --> void", pybind11::arg("s"));
		cl.def("select", (void (Pythia8::Nucleon::*)(class Pythia8::EventInfo &, enum Pythia8::Nucleon::Status)) &Pythia8::Nucleon::select, "C++: Pythia8::Nucleon::select(class Pythia8::EventInfo &, enum Pythia8::Nucleon::Status) --> void", pybind11::arg("evp"), pybind11::arg("s"));
		cl.def("select", (void (Pythia8::Nucleon::*)()) &Pythia8::Nucleon::select, "C++: Pythia8::Nucleon::select() --> void");
		cl.def("debug", (void (Pythia8::Nucleon::*)()) &Pythia8::Nucleon::debug, "C++: Pythia8::Nucleon::debug() --> void");
		cl.def("reset", (void (Pythia8::Nucleon::*)()) &Pythia8::Nucleon::reset, "C++: Pythia8::Nucleon::reset() --> void");
		cl.def("assign", (class Pythia8::Nucleon & (Pythia8::Nucleon::*)(const class Pythia8::Nucleon &)) &Pythia8::Nucleon::operator=, "C++: Pythia8::Nucleon::operator=(const class Pythia8::Nucleon &) --> class Pythia8::Nucleon &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::Nucleus file:Pythia8/HINucleusModel.h line:152
		pybind11::class_<Pythia8::Nucleus, std::shared_ptr<Pythia8::Nucleus>> cl(M("Pythia8"), "Nucleus", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Nucleus(); } ) );
		cl.def( pybind11::init<class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >, class Pythia8::Vec4>(), pybind11::arg("nucleons"), pybind11::arg("bPos") );

		cl.def("assign", (class Pythia8::Nucleus & (Pythia8::Nucleus::*)(const class Pythia8::Nucleus &)) &Pythia8::Nucleus::operator=, "C++: Pythia8::Nucleus::operator=(const class Pythia8::Nucleus &) --> class Pythia8::Nucleus &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::NucleusModel file:Pythia8/HINucleusModel.h line:187
		pybind11::class_<Pythia8::NucleusModel, std::shared_ptr<Pythia8::NucleusModel>, PyCallBack_Pythia8_NucleusModel> cl(M("Pythia8"), "NucleusModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_NucleusModel(); } ) );
		cl.def(pybind11::init<PyCallBack_Pythia8_NucleusModel const &>());
		cl.def_readwrite("isProj", &Pythia8::NucleusModel::isProj);
		cl.def_readwrite("idSave", &Pythia8::NucleusModel::idSave);
		cl.def_readwrite("ISave", &Pythia8::NucleusModel::ISave);
		cl.def_readwrite("ASave", &Pythia8::NucleusModel::ASave);
		cl.def_readwrite("ZSave", &Pythia8::NucleusModel::ZSave);
		cl.def_readwrite("LSave", &Pythia8::NucleusModel::LSave);
		cl.def_readwrite("RSave", &Pythia8::NucleusModel::RSave);
		cl.def_readwrite("mSave", &Pythia8::NucleusModel::mSave);
		cl.def_readwrite("pNSave", &Pythia8::NucleusModel::pNSave);
		cl.def_readwrite("mNSave", &Pythia8::NucleusModel::mNSave);
		cl.def_readwrite("idNSave", &Pythia8::NucleusModel::idNSave);
		cl.def_static("create", (class std::shared_ptr<class Pythia8::NucleusModel> (*)(int)) &Pythia8::NucleusModel::create, "C++: Pythia8::NucleusModel::create(int) --> class std::shared_ptr<class Pythia8::NucleusModel>", pybind11::arg("model"));
		cl.def("initPtr", (void (Pythia8::NucleusModel::*)(int, bool, class Pythia8::Info &)) &Pythia8::NucleusModel::initPtr, "C++: Pythia8::NucleusModel::initPtr(int, bool, class Pythia8::Info &) --> void", pybind11::arg("idIn"), pybind11::arg("isProjIn"), pybind11::arg("infoIn"));
		cl.def("init", (bool (Pythia8::NucleusModel::*)()) &Pythia8::NucleusModel::init, "C++: Pythia8::NucleusModel::init() --> bool");
		cl.def("initGeometry", (bool (Pythia8::NucleusModel::*)()) &Pythia8::NucleusModel::initGeometry, "C++: Pythia8::NucleusModel::initGeometry() --> bool");
		cl.def("setParticle", (void (Pythia8::NucleusModel::*)(int)) &Pythia8::NucleusModel::setParticle, "C++: Pythia8::NucleusModel::setParticle(int) --> void", pybind11::arg("idIn"));
		cl.def("setPN", (void (Pythia8::NucleusModel::*)(const class Pythia8::Vec4 &)) &Pythia8::NucleusModel::setPN, "C++: Pythia8::NucleusModel::setPN(const class Pythia8::Vec4 &) --> void", pybind11::arg("pNIn"));
		cl.def("setMN", (void (Pythia8::NucleusModel::*)(double)) &Pythia8::NucleusModel::setMN, "C++: Pythia8::NucleusModel::setMN(double) --> void", pybind11::arg("mNIn"));
		cl.def("produceIon", (class Pythia8::Particle (Pythia8::NucleusModel::*)()) &Pythia8::NucleusModel::produceIon, "C++: Pythia8::NucleusModel::produceIon() --> class Pythia8::Particle");
		cl.def("generate", (class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> > (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::generate, "C++: Pythia8::NucleusModel::generate() const --> class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >");
		cl.def("id", (int (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::id, "C++: Pythia8::NucleusModel::id() const --> int");
		cl.def("I", (int (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::I, "C++: Pythia8::NucleusModel::I() const --> int");
		cl.def("A", (int (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::A, "C++: Pythia8::NucleusModel::A() const --> int");
		cl.def("Z", (int (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::Z, "C++: Pythia8::NucleusModel::Z() const --> int");
		cl.def("L", (int (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::L, "C++: Pythia8::NucleusModel::L() const --> int");
		cl.def("R", (double (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::R, "C++: Pythia8::NucleusModel::R() const --> double");
		cl.def("idN", (int (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::idN, "C++: Pythia8::NucleusModel::idN() const --> int");
		cl.def("pN", (const class Pythia8::Vec4 & (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::pN, "C++: Pythia8::NucleusModel::pN() const --> const class Pythia8::Vec4 &", pybind11::return_value_policy::reference);
		cl.def("mN", (double (Pythia8::NucleusModel::*)() const) &Pythia8::NucleusModel::mN, "C++: Pythia8::NucleusModel::mN() const --> double");
		cl.def("assign", (class Pythia8::NucleusModel & (Pythia8::NucleusModel::*)(const class Pythia8::NucleusModel &)) &Pythia8::NucleusModel::operator=, "C++: Pythia8::NucleusModel::operator=(const class Pythia8::NucleusModel &) --> class Pythia8::NucleusModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
