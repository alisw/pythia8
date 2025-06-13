#include <Pythia8/SusyLesHouches.h>
#include <cwchar>
#include <ios>
#include <istream>
#include <iterator>
#include <memory>
#include <sstream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
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

void bind_Pythia8_SusyLesHouches(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::LHgenericBlock file:Pythia8/SusyLesHouches.h line:111
		pybind11::class_<Pythia8::LHgenericBlock, std::shared_ptr<Pythia8::LHgenericBlock>, Pythia8::LHblock<std::string>> cl(M("Pythia8"), "LHgenericBlock", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHgenericBlock(); } ) );
		cl.def( pybind11::init( [](Pythia8::LHgenericBlock const &o){ return new Pythia8::LHgenericBlock(o); } ) );
		cl.def("set", (int (Pythia8::LHgenericBlock::*)(std::string)) &Pythia8::LHgenericBlock::set, "C++: Pythia8::LHgenericBlock::set(std::string) --> int", pybind11::arg("lineIn"));
		cl.def("assign", (class Pythia8::LHgenericBlock & (Pythia8::LHgenericBlock::*)(const class Pythia8::LHgenericBlock &)) &Pythia8::LHgenericBlock::operator=, "C++: Pythia8::LHgenericBlock::operator=(const class Pythia8::LHgenericBlock &) --> class Pythia8::LHgenericBlock &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LHdecayChannel file:Pythia8/SusyLesHouches.h line:295
		pybind11::class_<Pythia8::LHdecayChannel, std::shared_ptr<Pythia8::LHdecayChannel>> cl(M("Pythia8"), "LHdecayChannel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHdecayChannel(); } ) );
		cl.def( pybind11::init( [](double const & a0, int const & a1, class std::vector<int, class std::allocator<int> > const & a2){ return new Pythia8::LHdecayChannel(a0, a1, a2); } ), "doc" , pybind11::arg("bratIn"), pybind11::arg("nDaIn"), pybind11::arg("idDaIn"));
		cl.def( pybind11::init<double, int, class std::vector<int, class std::allocator<int> >, std::string>(), pybind11::arg("bratIn"), pybind11::arg("nDaIn"), pybind11::arg("idDaIn"), pybind11::arg("cIn") );

		cl.def( pybind11::init( [](Pythia8::LHdecayChannel const &o){ return new Pythia8::LHdecayChannel(o); } ) );
		cl.def("setChannel", [](Pythia8::LHdecayChannel &o, double const & a0, int const & a1, class std::vector<int, class std::allocator<int> > const & a2) -> void { return o.setChannel(a0, a1, a2); }, "", pybind11::arg("bratIn"), pybind11::arg("nDaIn"), pybind11::arg("idDaIn"));
		cl.def("setChannel", (void (Pythia8::LHdecayChannel::*)(double, int, class std::vector<int, class std::allocator<int> >, std::string)) &Pythia8::LHdecayChannel::setChannel, "C++: Pythia8::LHdecayChannel::setChannel(double, int, class std::vector<int, class std::allocator<int> >, std::string) --> void", pybind11::arg("bratIn"), pybind11::arg("nDaIn"), pybind11::arg("idDaIn"), pybind11::arg("cIn"));
		cl.def("setBrat", (void (Pythia8::LHdecayChannel::*)(double)) &Pythia8::LHdecayChannel::setBrat, "C++: Pythia8::LHdecayChannel::setBrat(double) --> void", pybind11::arg("bratIn"));
		cl.def("setIdDa", (void (Pythia8::LHdecayChannel::*)(class std::vector<int, class std::allocator<int> >)) &Pythia8::LHdecayChannel::setIdDa, "C++: Pythia8::LHdecayChannel::setIdDa(class std::vector<int, class std::allocator<int> >) --> void", pybind11::arg("idDaIn"));
		cl.def("getBrat", (double (Pythia8::LHdecayChannel::*)()) &Pythia8::LHdecayChannel::getBrat, "C++: Pythia8::LHdecayChannel::getBrat() --> double");
		cl.def("getNDa", (int (Pythia8::LHdecayChannel::*)()) &Pythia8::LHdecayChannel::getNDa, "C++: Pythia8::LHdecayChannel::getNDa() --> int");
		cl.def("getIdDa", (class std::vector<int, class std::allocator<int> > (Pythia8::LHdecayChannel::*)()) &Pythia8::LHdecayChannel::getIdDa, "C++: Pythia8::LHdecayChannel::getIdDa() --> class std::vector<int, class std::allocator<int> >");
		cl.def("getComment", (std::string (Pythia8::LHdecayChannel::*)()) &Pythia8::LHdecayChannel::getComment, "C++: Pythia8::LHdecayChannel::getComment() --> std::string");
		cl.def("assign", (class Pythia8::LHdecayChannel & (Pythia8::LHdecayChannel::*)(const class Pythia8::LHdecayChannel &)) &Pythia8::LHdecayChannel::operator=, "C++: Pythia8::LHdecayChannel::operator=(const class Pythia8::LHdecayChannel &) --> class Pythia8::LHdecayChannel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LHdecayTable file:Pythia8/SusyLesHouches.h line:328
		pybind11::class_<Pythia8::LHdecayTable, std::shared_ptr<Pythia8::LHdecayTable>> cl(M("Pythia8"), "LHdecayTable", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHdecayTable(); } ) );
		cl.def( pybind11::init<int>(), pybind11::arg("idIn") );

		cl.def( pybind11::init<int, double>(), pybind11::arg("idIn"), pybind11::arg("widthIn") );

		cl.def("getId", (int (Pythia8::LHdecayTable::*)()) &Pythia8::LHdecayTable::getId, "C++: Pythia8::LHdecayTable::getId() --> int");
		cl.def("getWidth", (double (Pythia8::LHdecayTable::*)()) &Pythia8::LHdecayTable::getWidth, "C++: Pythia8::LHdecayTable::getWidth() --> double");
		cl.def("setId", (void (Pythia8::LHdecayTable::*)(int)) &Pythia8::LHdecayTable::setId, "C++: Pythia8::LHdecayTable::setId(int) --> void", pybind11::arg("idIn"));
		cl.def("setWidth", (void (Pythia8::LHdecayTable::*)(double)) &Pythia8::LHdecayTable::setWidth, "C++: Pythia8::LHdecayTable::setWidth(double) --> void", pybind11::arg("widthIn"));
		cl.def("reset", [](Pythia8::LHdecayTable &o) -> void { return o.reset(); }, "");
		cl.def("reset", (void (Pythia8::LHdecayTable::*)(double)) &Pythia8::LHdecayTable::reset, "C++: Pythia8::LHdecayTable::reset(double) --> void", pybind11::arg("widthIn"));
		cl.def("addChannel", (void (Pythia8::LHdecayTable::*)(class Pythia8::LHdecayChannel)) &Pythia8::LHdecayTable::addChannel, "C++: Pythia8::LHdecayTable::addChannel(class Pythia8::LHdecayChannel) --> void", pybind11::arg("channelIn"));
		cl.def("addChannel", [](Pythia8::LHdecayTable &o, double const & a0, int const & a1, class std::vector<int, class std::allocator<int> > const & a2) -> void { return o.addChannel(a0, a1, a2); }, "", pybind11::arg("bratIn"), pybind11::arg("nDaIn"), pybind11::arg("idDaIn"));
		cl.def("addChannel", (void (Pythia8::LHdecayTable::*)(double, int, class std::vector<int, class std::allocator<int> >, std::string)) &Pythia8::LHdecayTable::addChannel, "C++: Pythia8::LHdecayTable::addChannel(double, int, class std::vector<int, class std::allocator<int> >, std::string) --> void", pybind11::arg("bratIn"), pybind11::arg("nDaIn"), pybind11::arg("idDaIn"), pybind11::arg("cIn"));
		cl.def("size", (int (Pythia8::LHdecayTable::*)()) &Pythia8::LHdecayTable::size, "C++: Pythia8::LHdecayTable::size() --> int");
		cl.def("getBrat", (double (Pythia8::LHdecayTable::*)(int)) &Pythia8::LHdecayTable::getBrat, "C++: Pythia8::LHdecayTable::getBrat(int) --> double", pybind11::arg("iChannel"));
		cl.def("getIdDa", (class std::vector<int, class std::allocator<int> > (Pythia8::LHdecayTable::*)(int)) &Pythia8::LHdecayTable::getIdDa, "C++: Pythia8::LHdecayTable::getIdDa(int) --> class std::vector<int, class std::allocator<int> >", pybind11::arg("iChannel"));
		cl.def("getChannel", (class Pythia8::LHdecayChannel (Pythia8::LHdecayTable::*)(int)) &Pythia8::LHdecayTable::getChannel, "C++: Pythia8::LHdecayTable::getChannel(int) --> class Pythia8::LHdecayChannel", pybind11::arg("iChannel"));
	}
	{ // Pythia8::SusyLesHouches file:Pythia8/SusyLesHouches.h line:393
		pybind11::class_<Pythia8::SusyLesHouches, std::shared_ptr<Pythia8::SusyLesHouches>> cl(M("Pythia8"), "SusyLesHouches", "");
		pybind11::handle cl_type = cl;

		{ // Pythia8::SusyLesHouches::Entry file:Pythia8/SusyLesHouches.h line:425
			auto & enclosing_class = cl;
			pybind11::class_<Pythia8::SusyLesHouches::Entry, std::shared_ptr<Pythia8::SusyLesHouches::Entry>> cl(enclosing_class, "Entry", "");
			pybind11::handle cl_type = cl;

			cl.def( pybind11::init( [](){ return new Pythia8::SusyLesHouches::Entry(); } ) );
			cl.def("isInt", (bool (Pythia8::SusyLesHouches::Entry::*)()) &Pythia8::SusyLesHouches::Entry::isInt, "C++: Pythia8::SusyLesHouches::Entry::isInt() --> bool");
			cl.def("isDouble", (bool (Pythia8::SusyLesHouches::Entry::*)()) &Pythia8::SusyLesHouches::Entry::isDouble, "C++: Pythia8::SusyLesHouches::Entry::isDouble() --> bool");
			cl.def("isString", (bool (Pythia8::SusyLesHouches::Entry::*)()) &Pythia8::SusyLesHouches::Entry::isString, "C++: Pythia8::SusyLesHouches::Entry::isString() --> bool");
			cl.def("assign", (class Pythia8::SusyLesHouches::Entry & (Pythia8::SusyLesHouches::Entry::*)(double &)) &Pythia8::SusyLesHouches::Entry::operator=, "C++: Pythia8::SusyLesHouches::Entry::operator=(double &) --> class Pythia8::SusyLesHouches::Entry &", pybind11::return_value_policy::reference, pybind11::arg("val"));
			cl.def("assign", (class Pythia8::SusyLesHouches::Entry & (Pythia8::SusyLesHouches::Entry::*)(int &)) &Pythia8::SusyLesHouches::Entry::operator=, "C++: Pythia8::SusyLesHouches::Entry::operator=(int &) --> class Pythia8::SusyLesHouches::Entry &", pybind11::return_value_policy::reference, pybind11::arg("val"));
			cl.def("assign", (class Pythia8::SusyLesHouches::Entry & (Pythia8::SusyLesHouches::Entry::*)(std::string &)) &Pythia8::SusyLesHouches::Entry::operator=, "C++: Pythia8::SusyLesHouches::Entry::operator=(std::string &) --> class Pythia8::SusyLesHouches::Entry &", pybind11::return_value_policy::reference, pybind11::arg("val"));
			cl.def("setComment", (void (Pythia8::SusyLesHouches::Entry::*)(std::string)) &Pythia8::SusyLesHouches::Entry::setComment, "C++: Pythia8::SusyLesHouches::Entry::setComment(std::string) --> void", pybind11::arg("comment"));
			cl.def("getComment", (void (Pythia8::SusyLesHouches::Entry::*)(std::string &)) &Pythia8::SusyLesHouches::Entry::getComment, "C++: Pythia8::SusyLesHouches::Entry::getComment(std::string &) --> void", pybind11::arg("comment"));
			cl.def("get", (bool (Pythia8::SusyLesHouches::Entry::*)(int &)) &Pythia8::SusyLesHouches::Entry::get, "C++: Pythia8::SusyLesHouches::Entry::get(int &) --> bool", pybind11::arg("val"));
			cl.def("get", (bool (Pythia8::SusyLesHouches::Entry::*)(double &)) &Pythia8::SusyLesHouches::Entry::get, "C++: Pythia8::SusyLesHouches::Entry::get(double &) --> bool", pybind11::arg("val"));
			cl.def("get", (bool (Pythia8::SusyLesHouches::Entry::*)(std::string &)) &Pythia8::SusyLesHouches::Entry::get, "C++: Pythia8::SusyLesHouches::Entry::get(std::string &) --> bool", pybind11::arg("val"));
		}

		cl.def( pybind11::init( [](){ return new Pythia8::SusyLesHouches(); } ), "doc" );
		cl.def( pybind11::init<int>(), pybind11::arg("verboseIn") );

		cl.def( pybind11::init( [](class std::basic_string<char> const & a0){ return new Pythia8::SusyLesHouches(a0); } ), "doc" , pybind11::arg("filename"));
		cl.def( pybind11::init<std::string, int>(), pybind11::arg("filename"), pybind11::arg("verboseIn") );

		cl.def( pybind11::init( [](Pythia8::SusyLesHouches const &o){ return new Pythia8::SusyLesHouches(o); } ) );
		cl.def_readwrite("slhaFile", &Pythia8::SusyLesHouches::slhaFile);
		cl.def_readwrite("modsel", &Pythia8::SusyLesHouches::modsel);
		cl.def_readwrite("modsel21", &Pythia8::SusyLesHouches::modsel21);
		cl.def_readwrite("modsel12", &Pythia8::SusyLesHouches::modsel12);
		cl.def_readwrite("minpar", &Pythia8::SusyLesHouches::minpar);
		cl.def_readwrite("extpar", &Pythia8::SusyLesHouches::extpar);
		cl.def_readwrite("sminputs", &Pythia8::SusyLesHouches::sminputs);
		cl.def_readwrite("spinfo", &Pythia8::SusyLesHouches::spinfo);
		cl.def_readwrite("spinfo3", &Pythia8::SusyLesHouches::spinfo3);
		cl.def_readwrite("spinfo4", &Pythia8::SusyLesHouches::spinfo4);
		cl.def_readwrite("dcinfo", &Pythia8::SusyLesHouches::dcinfo);
		cl.def_readwrite("dcinfo3", &Pythia8::SusyLesHouches::dcinfo3);
		cl.def_readwrite("dcinfo4", &Pythia8::SusyLesHouches::dcinfo4);
		cl.def_readwrite("mass", &Pythia8::SusyLesHouches::mass);
		cl.def_readwrite("nmix", &Pythia8::SusyLesHouches::nmix);
		cl.def_readwrite("umix", &Pythia8::SusyLesHouches::umix);
		cl.def_readwrite("vmix", &Pythia8::SusyLesHouches::vmix);
		cl.def_readwrite("stopmix", &Pythia8::SusyLesHouches::stopmix);
		cl.def_readwrite("sbotmix", &Pythia8::SusyLesHouches::sbotmix);
		cl.def_readwrite("staumix", &Pythia8::SusyLesHouches::staumix);
		cl.def_readwrite("alpha", &Pythia8::SusyLesHouches::alpha);
		cl.def_readwrite("hmix", &Pythia8::SusyLesHouches::hmix);
		cl.def_readwrite("gauge", &Pythia8::SusyLesHouches::gauge);
		cl.def_readwrite("msoft", &Pythia8::SusyLesHouches::msoft);
		cl.def_readwrite("au", &Pythia8::SusyLesHouches::au);
		cl.def_readwrite("ad", &Pythia8::SusyLesHouches::ad);
		cl.def_readwrite("ae", &Pythia8::SusyLesHouches::ae);
		cl.def_readwrite("yu", &Pythia8::SusyLesHouches::yu);
		cl.def_readwrite("yd", &Pythia8::SusyLesHouches::yd);
		cl.def_readwrite("ye", &Pythia8::SusyLesHouches::ye);
		cl.def_readwrite("decays", &Pythia8::SusyLesHouches::decays);
		cl.def_readwrite("decayIndices", &Pythia8::SusyLesHouches::decayIndices);
		cl.def_readwrite("qnumbers", &Pythia8::SusyLesHouches::qnumbers);
		cl.def_readwrite("qnumbersName", &Pythia8::SusyLesHouches::qnumbersName);
		cl.def_readwrite("qnumbersAntiName", &Pythia8::SusyLesHouches::qnumbersAntiName);
		cl.def_readwrite("qextpar", &Pythia8::SusyLesHouches::qextpar);
		cl.def_readwrite("vckmin", &Pythia8::SusyLesHouches::vckmin);
		cl.def_readwrite("upmnsin", &Pythia8::SusyLesHouches::upmnsin);
		cl.def_readwrite("msq2in", &Pythia8::SusyLesHouches::msq2in);
		cl.def_readwrite("msu2in", &Pythia8::SusyLesHouches::msu2in);
		cl.def_readwrite("msd2in", &Pythia8::SusyLesHouches::msd2in);
		cl.def_readwrite("msl2in", &Pythia8::SusyLesHouches::msl2in);
		cl.def_readwrite("mse2in", &Pythia8::SusyLesHouches::mse2in);
		cl.def_readwrite("tuin", &Pythia8::SusyLesHouches::tuin);
		cl.def_readwrite("tdin", &Pythia8::SusyLesHouches::tdin);
		cl.def_readwrite("tein", &Pythia8::SusyLesHouches::tein);
		cl.def_readwrite("vckm", &Pythia8::SusyLesHouches::vckm);
		cl.def_readwrite("upmns", &Pythia8::SusyLesHouches::upmns);
		cl.def_readwrite("msq2", &Pythia8::SusyLesHouches::msq2);
		cl.def_readwrite("msu2", &Pythia8::SusyLesHouches::msu2);
		cl.def_readwrite("msd2", &Pythia8::SusyLesHouches::msd2);
		cl.def_readwrite("msl2", &Pythia8::SusyLesHouches::msl2);
		cl.def_readwrite("mse2", &Pythia8::SusyLesHouches::mse2);
		cl.def_readwrite("tu", &Pythia8::SusyLesHouches::tu);
		cl.def_readwrite("td", &Pythia8::SusyLesHouches::td);
		cl.def_readwrite("te", &Pythia8::SusyLesHouches::te);
		cl.def_readwrite("usqmix", &Pythia8::SusyLesHouches::usqmix);
		cl.def_readwrite("dsqmix", &Pythia8::SusyLesHouches::dsqmix);
		cl.def_readwrite("selmix", &Pythia8::SusyLesHouches::selmix);
		cl.def_readwrite("snumix", &Pythia8::SusyLesHouches::snumix);
		cl.def_readwrite("snsmix", &Pythia8::SusyLesHouches::snsmix);
		cl.def_readwrite("snamix", &Pythia8::SusyLesHouches::snamix);
		cl.def_readwrite("rvlamllein", &Pythia8::SusyLesHouches::rvlamllein);
		cl.def_readwrite("rvlamlqdin", &Pythia8::SusyLesHouches::rvlamlqdin);
		cl.def_readwrite("rvlamuddin", &Pythia8::SusyLesHouches::rvlamuddin);
		cl.def_readwrite("rvtllein", &Pythia8::SusyLesHouches::rvtllein);
		cl.def_readwrite("rvtlqdin", &Pythia8::SusyLesHouches::rvtlqdin);
		cl.def_readwrite("rvtuddin", &Pythia8::SusyLesHouches::rvtuddin);
		cl.def_readwrite("rvkappain", &Pythia8::SusyLesHouches::rvkappain);
		cl.def_readwrite("rvdin", &Pythia8::SusyLesHouches::rvdin);
		cl.def_readwrite("rvm2lh1in", &Pythia8::SusyLesHouches::rvm2lh1in);
		cl.def_readwrite("rvsnvevin", &Pythia8::SusyLesHouches::rvsnvevin);
		cl.def_readwrite("rvlamlle", &Pythia8::SusyLesHouches::rvlamlle);
		cl.def_readwrite("rvlamlqd", &Pythia8::SusyLesHouches::rvlamlqd);
		cl.def_readwrite("rvlamudd", &Pythia8::SusyLesHouches::rvlamudd);
		cl.def_readwrite("rvtlle", &Pythia8::SusyLesHouches::rvtlle);
		cl.def_readwrite("rvtlqd", &Pythia8::SusyLesHouches::rvtlqd);
		cl.def_readwrite("rvtudd", &Pythia8::SusyLesHouches::rvtudd);
		cl.def_readwrite("rvkappa", &Pythia8::SusyLesHouches::rvkappa);
		cl.def_readwrite("rvd", &Pythia8::SusyLesHouches::rvd);
		cl.def_readwrite("rvm2lh1", &Pythia8::SusyLesHouches::rvm2lh1);
		cl.def_readwrite("rvsnvev", &Pythia8::SusyLesHouches::rvsnvev);
		cl.def_readwrite("rvnmix", &Pythia8::SusyLesHouches::rvnmix);
		cl.def_readwrite("rvumix", &Pythia8::SusyLesHouches::rvumix);
		cl.def_readwrite("rvvmix", &Pythia8::SusyLesHouches::rvvmix);
		cl.def_readwrite("rvhmix", &Pythia8::SusyLesHouches::rvhmix);
		cl.def_readwrite("rvamix", &Pythia8::SusyLesHouches::rvamix);
		cl.def_readwrite("rvlmix", &Pythia8::SusyLesHouches::rvlmix);
		cl.def_readwrite("imminpar", &Pythia8::SusyLesHouches::imminpar);
		cl.def_readwrite("imextpar", &Pythia8::SusyLesHouches::imextpar);
		cl.def_readwrite("cvhmix", &Pythia8::SusyLesHouches::cvhmix);
		cl.def_readwrite("imcvhmix", &Pythia8::SusyLesHouches::imcvhmix);
		cl.def_readwrite("imau", &Pythia8::SusyLesHouches::imau);
		cl.def_readwrite("imad", &Pythia8::SusyLesHouches::imad);
		cl.def_readwrite("imae", &Pythia8::SusyLesHouches::imae);
		cl.def_readwrite("imhmix", &Pythia8::SusyLesHouches::imhmix);
		cl.def_readwrite("immsoft", &Pythia8::SusyLesHouches::immsoft);
		cl.def_readwrite("immsq2in", &Pythia8::SusyLesHouches::immsq2in);
		cl.def_readwrite("immsu2in", &Pythia8::SusyLesHouches::immsu2in);
		cl.def_readwrite("immsd2in", &Pythia8::SusyLesHouches::immsd2in);
		cl.def_readwrite("immsl2in", &Pythia8::SusyLesHouches::immsl2in);
		cl.def_readwrite("immse2in", &Pythia8::SusyLesHouches::immse2in);
		cl.def_readwrite("imtuin", &Pythia8::SusyLesHouches::imtuin);
		cl.def_readwrite("imtdin", &Pythia8::SusyLesHouches::imtdin);
		cl.def_readwrite("imtein", &Pythia8::SusyLesHouches::imtein);
		cl.def_readwrite("imvckm", &Pythia8::SusyLesHouches::imvckm);
		cl.def_readwrite("imupmns", &Pythia8::SusyLesHouches::imupmns);
		cl.def_readwrite("immsq2", &Pythia8::SusyLesHouches::immsq2);
		cl.def_readwrite("immsu2", &Pythia8::SusyLesHouches::immsu2);
		cl.def_readwrite("immsd2", &Pythia8::SusyLesHouches::immsd2);
		cl.def_readwrite("immsl2", &Pythia8::SusyLesHouches::immsl2);
		cl.def_readwrite("immse2", &Pythia8::SusyLesHouches::immse2);
		cl.def_readwrite("imtu", &Pythia8::SusyLesHouches::imtu);
		cl.def_readwrite("imtd", &Pythia8::SusyLesHouches::imtd);
		cl.def_readwrite("imte", &Pythia8::SusyLesHouches::imte);
		cl.def_readwrite("imusqmix", &Pythia8::SusyLesHouches::imusqmix);
		cl.def_readwrite("imdsqmix", &Pythia8::SusyLesHouches::imdsqmix);
		cl.def_readwrite("imselmix", &Pythia8::SusyLesHouches::imselmix);
		cl.def_readwrite("imsnumix", &Pythia8::SusyLesHouches::imsnumix);
		cl.def_readwrite("imnmix", &Pythia8::SusyLesHouches::imnmix);
		cl.def_readwrite("imumix", &Pythia8::SusyLesHouches::imumix);
		cl.def_readwrite("imvmix", &Pythia8::SusyLesHouches::imvmix);
		cl.def_readwrite("nmssmrun", &Pythia8::SusyLesHouches::nmssmrun);
		cl.def_readwrite("nmhmix", &Pythia8::SusyLesHouches::nmhmix);
		cl.def_readwrite("nmamix", &Pythia8::SusyLesHouches::nmamix);
		cl.def_readwrite("nmnmix", &Pythia8::SusyLesHouches::nmnmix);
		cl.def_readwrite("imnmnmix", &Pythia8::SusyLesHouches::imnmnmix);
		cl.def_readwrite("genericBlocks", &Pythia8::SusyLesHouches::genericBlocks);
		cl.def("readFile", [](Pythia8::SusyLesHouches &o) -> int { return o.readFile(); }, "");
		cl.def("readFile", [](Pythia8::SusyLesHouches &o, class std::basic_string<char> const & a0) -> int { return o.readFile(a0); }, "", pybind11::arg("slhaFileIn"));
		cl.def("readFile", [](Pythia8::SusyLesHouches &o, class std::basic_string<char> const & a0, int const & a1) -> int { return o.readFile(a0, a1); }, "", pybind11::arg("slhaFileIn"), pybind11::arg("verboseIn"));
		cl.def("readFile", (int (Pythia8::SusyLesHouches::*)(std::string, int, bool)) &Pythia8::SusyLesHouches::readFile, "C++: Pythia8::SusyLesHouches::readFile(std::string, int, bool) --> int", pybind11::arg("slhaFileIn"), pybind11::arg("verboseIn"), pybind11::arg("useDecayIn"));
		cl.def("readFile", [](Pythia8::SusyLesHouches &o, class std::basic_istream<char> & a0) -> int { return o.readFile(a0); }, "", pybind11::arg(""));
		cl.def("readFile", [](Pythia8::SusyLesHouches &o, class std::basic_istream<char> & a0, int const & a1) -> int { return o.readFile(a0, a1); }, "", pybind11::arg(""), pybind11::arg("verboseIn"));
		cl.def("readFile", (int (Pythia8::SusyLesHouches::*)(class std::basic_istream<char> &, int, bool)) &Pythia8::SusyLesHouches::readFile, "C++: Pythia8::SusyLesHouches::readFile(class std::basic_istream<char> &, int, bool) --> int", pybind11::arg(""), pybind11::arg("verboseIn"), pybind11::arg("useDecayIn"));
		cl.def("listHeader", (void (Pythia8::SusyLesHouches::*)()) &Pythia8::SusyLesHouches::listHeader, "C++: Pythia8::SusyLesHouches::listHeader() --> void");
		cl.def("listFooter", (void (Pythia8::SusyLesHouches::*)()) &Pythia8::SusyLesHouches::listFooter, "C++: Pythia8::SusyLesHouches::listFooter() --> void");
		cl.def("listSpectrum", [](Pythia8::SusyLesHouches &o) -> void { return o.listSpectrum(); }, "");
		cl.def("listSpectrum", (void (Pythia8::SusyLesHouches::*)(int)) &Pythia8::SusyLesHouches::listSpectrum, "C++: Pythia8::SusyLesHouches::listSpectrum(int) --> void", pybind11::arg("ifail"));
		cl.def("checkSpectrum", (int (Pythia8::SusyLesHouches::*)()) &Pythia8::SusyLesHouches::checkSpectrum, "C++: Pythia8::SusyLesHouches::checkSpectrum() --> int");
		cl.def("verbose", (int (Pythia8::SusyLesHouches::*)()) &Pythia8::SusyLesHouches::verbose, "C++: Pythia8::SusyLesHouches::verbose() --> int");
		cl.def("verbose", (void (Pythia8::SusyLesHouches::*)(int)) &Pythia8::SusyLesHouches::verbose, "C++: Pythia8::SusyLesHouches::verbose(int) --> void", pybind11::arg("verboseIn"));
		cl.def("message", [](Pythia8::SusyLesHouches &o, int const & a0, class std::basic_string<char> const & a1, class std::basic_string<char> const & a2) -> void { return o.message(a0, a1, a2); }, "", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""));
		cl.def("message", (void (Pythia8::SusyLesHouches::*)(int, std::string, std::string, int)) &Pythia8::SusyLesHouches::message, "C++: Pythia8::SusyLesHouches::message(int, std::string, std::string, int) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg(""), pybind11::arg("line"));
	}
}
