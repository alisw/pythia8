#include <Pythia8/Basics.h>
#include <cwchar>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <locale>
#include <memory>
#include <ostream>
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

void bind_Pythia8_Basics_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::RndmState file:Pythia8/Basics.h line:373
		pybind11::class_<Pythia8::RndmState, std::shared_ptr<Pythia8::RndmState>> cl(M("Pythia8"), "RndmState", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](Pythia8::RndmState const &o){ return new Pythia8::RndmState(o); } ) );
		cl.def( pybind11::init( [](){ return new Pythia8::RndmState(); } ) );
		cl.def_readwrite("i97", &Pythia8::RndmState::i97);
		cl.def_readwrite("j97", &Pythia8::RndmState::j97);
		cl.def_readwrite("seed", &Pythia8::RndmState::seed);
		cl.def_readwrite("sequence", &Pythia8::RndmState::sequence);
		cl.def_readwrite("c", &Pythia8::RndmState::c);
		cl.def_readwrite("cd", &Pythia8::RndmState::cd);
		cl.def_readwrite("cm", &Pythia8::RndmState::cm);
		cl.def("__eq__", (bool (Pythia8::RndmState::*)(const struct Pythia8::RndmState &) const) &Pythia8::RndmState::operator==, "C++: Pythia8::RndmState::operator==(const struct Pythia8::RndmState &) const --> bool", pybind11::arg("other"));
		cl.def("assign", (struct Pythia8::RndmState & (Pythia8::RndmState::*)(const struct Pythia8::RndmState &)) &Pythia8::RndmState::operator=, "C++: Pythia8::RndmState::operator=(const struct Pythia8::RndmState &) --> struct Pythia8::RndmState &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::Rndm file:Pythia8/Basics.h line:388
		pybind11::class_<Pythia8::Rndm, std::shared_ptr<Pythia8::Rndm>> cl(M("Pythia8"), "Rndm", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Rndm(); } ) );
		cl.def( pybind11::init<int>(), pybind11::arg("seedIn") );

		cl.def( pybind11::init( [](Pythia8::Rndm const &o){ return new Pythia8::Rndm(o); } ) );
		cl.def("rndmEnginePtr", (bool (Pythia8::Rndm::*)(class std::shared_ptr<class Pythia8::RndmEngine>)) &Pythia8::Rndm::rndmEnginePtr, "C++: Pythia8::Rndm::rndmEnginePtr(class std::shared_ptr<class Pythia8::RndmEngine>) --> bool", pybind11::arg("rndmEngPtrIn"));
		cl.def("init", [](Pythia8::Rndm &o) -> void { return o.init(); }, "");
		cl.def("init", (void (Pythia8::Rndm::*)(int)) &Pythia8::Rndm::init, "C++: Pythia8::Rndm::init(int) --> void", pybind11::arg("seedIn"));
		cl.def("flat", (double (Pythia8::Rndm::*)()) &Pythia8::Rndm::flat, "C++: Pythia8::Rndm::flat() --> double");
		cl.def("exp", (double (Pythia8::Rndm::*)()) &Pythia8::Rndm::exp, "C++: Pythia8::Rndm::exp() --> double");
		cl.def("xexp", (double (Pythia8::Rndm::*)()) &Pythia8::Rndm::xexp, "C++: Pythia8::Rndm::xexp() --> double");
		cl.def("gauss", (double (Pythia8::Rndm::*)()) &Pythia8::Rndm::gauss, "C++: Pythia8::Rndm::gauss() --> double");
		cl.def("gauss2", (struct std::pair<double, double> (Pythia8::Rndm::*)()) &Pythia8::Rndm::gauss2, "C++: Pythia8::Rndm::gauss2() --> struct std::pair<double, double>");
		cl.def("gamma", (double (Pythia8::Rndm::*)(double, double)) &Pythia8::Rndm::gamma, "C++: Pythia8::Rndm::gamma(double, double) --> double", pybind11::arg("k0"), pybind11::arg("r0"));
		cl.def("phaseSpace2", (struct std::pair<class Pythia8::Vec4, class Pythia8::Vec4> (Pythia8::Rndm::*)(double, double, double)) &Pythia8::Rndm::phaseSpace2, "C++: Pythia8::Rndm::phaseSpace2(double, double, double) --> struct std::pair<class Pythia8::Vec4, class Pythia8::Vec4>", pybind11::arg("eCM"), pybind11::arg("m1"), pybind11::arg("m2"));
		cl.def("pick", (int (Pythia8::Rndm::*)(const class std::vector<double, class std::allocator<double> > &)) &Pythia8::Rndm::pick, "C++: Pythia8::Rndm::pick(const class std::vector<double, class std::allocator<double> > &) --> int", pybind11::arg("prob"));
		cl.def("peekFlat", (double (Pythia8::Rndm::*)()) &Pythia8::Rndm::peekFlat, "C++: Pythia8::Rndm::peekFlat() --> double");
		cl.def("dumpState", (bool (Pythia8::Rndm::*)(std::string)) &Pythia8::Rndm::dumpState, "C++: Pythia8::Rndm::dumpState(std::string) --> bool", pybind11::arg("fileName"));
		cl.def("readState", (bool (Pythia8::Rndm::*)(std::string)) &Pythia8::Rndm::readState, "C++: Pythia8::Rndm::readState(std::string) --> bool", pybind11::arg("fileName"));
		cl.def("getState", (struct Pythia8::RndmState (Pythia8::Rndm::*)() const) &Pythia8::Rndm::getState, "C++: Pythia8::Rndm::getState() const --> struct Pythia8::RndmState");
		cl.def("setState", (void (Pythia8::Rndm::*)(const struct Pythia8::RndmState &)) &Pythia8::Rndm::setState, "C++: Pythia8::Rndm::setState(const struct Pythia8::RndmState &) --> void", pybind11::arg("state"));
	}
	{ // Pythia8::Hist file:Pythia8/Basics.h line:485
		pybind11::class_<Pythia8::Hist, std::shared_ptr<Pythia8::Hist>> cl(M("Pythia8"), "Hist", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Hist(); } ) );
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0){ return new Pythia8::Hist(a0); } ), "doc" , pybind11::arg("titleIn"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, int const & a1){ return new Pythia8::Hist(a0, a1); } ), "doc" , pybind11::arg("titleIn"), pybind11::arg("nBinIn"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, int const & a1, double const & a2){ return new Pythia8::Hist(a0, a1, a2); } ), "doc" , pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, int const & a1, double const & a2, double const & a3){ return new Pythia8::Hist(a0, a1, a2, a3); } ), "doc" , pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"), pybind11::arg("xMaxIn"));
		cl.def( pybind11::init( [](class std::basic_string<char> const & a0, int const & a1, double const & a2, double const & a3, bool const & a4){ return new Pythia8::Hist(a0, a1, a2, a3, a4); } ), "doc" , pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"), pybind11::arg("xMaxIn"), pybind11::arg("logXIn"));
		cl.def( pybind11::init<std::string, int, double, double, bool, bool>(), pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"), pybind11::arg("xMaxIn"), pybind11::arg("logXIn"), pybind11::arg("doStatsIn") );

		cl.def( pybind11::init( [](Pythia8::Hist const &o){ return new Pythia8::Hist(o); } ) );
		cl.def( pybind11::init<std::string, const class Pythia8::Hist &>(), pybind11::arg("titleIn"), pybind11::arg("h") );

		cl.def("assign", (class Pythia8::Hist & (Pythia8::Hist::*)(const class Pythia8::Hist &)) &Pythia8::Hist::operator=, "C++: Pythia8::Hist::operator=(const class Pythia8::Hist &) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("h"));
		cl.def_static("plotFunc", [](class std::function<double (double)> const & a0, class std::basic_string<char> const & a1, int const & a2, double const & a3, double const & a4) -> Pythia8::Hist { return Pythia8::Hist::plotFunc(a0, a1, a2, a3, a4); }, "", pybind11::arg("f"), pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"), pybind11::arg("xMaxIn"));
		cl.def_static("plotFunc", (class Pythia8::Hist (*)(class std::function<double (double)>, std::string, int, double, double, bool)) &Pythia8::Hist::plotFunc, "C++: Pythia8::Hist::plotFunc(class std::function<double (double)>, std::string, int, double, double, bool) --> class Pythia8::Hist", pybind11::arg("f"), pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"), pybind11::arg("xMaxIn"), pybind11::arg("logXIn"));
		cl.def("book", [](Pythia8::Hist &o) -> void { return o.book(); }, "");
		cl.def("book", [](Pythia8::Hist &o, class std::basic_string<char> const & a0) -> void { return o.book(a0); }, "", pybind11::arg("titleIn"));
		cl.def("book", [](Pythia8::Hist &o, class std::basic_string<char> const & a0, int const & a1) -> void { return o.book(a0, a1); }, "", pybind11::arg("titleIn"), pybind11::arg("nBinIn"));
		cl.def("book", [](Pythia8::Hist &o, class std::basic_string<char> const & a0, int const & a1, double const & a2) -> void { return o.book(a0, a1, a2); }, "", pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"));
		cl.def("book", [](Pythia8::Hist &o, class std::basic_string<char> const & a0, int const & a1, double const & a2, double const & a3) -> void { return o.book(a0, a1, a2, a3); }, "", pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"), pybind11::arg("xMaxIn"));
		cl.def("book", [](Pythia8::Hist &o, class std::basic_string<char> const & a0, int const & a1, double const & a2, double const & a3, bool const & a4) -> void { return o.book(a0, a1, a2, a3, a4); }, "", pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"), pybind11::arg("xMaxIn"), pybind11::arg("logXIn"));
		cl.def("book", (void (Pythia8::Hist::*)(std::string, int, double, double, bool, bool)) &Pythia8::Hist::book, "C++: Pythia8::Hist::book(std::string, int, double, double, bool, bool) --> void", pybind11::arg("titleIn"), pybind11::arg("nBinIn"), pybind11::arg("xMinIn"), pybind11::arg("xMaxIn"), pybind11::arg("logXIn"), pybind11::arg("doStatsIn"));
		cl.def("title", [](Pythia8::Hist &o) -> void { return o.title(); }, "");
		cl.def("title", (void (Pythia8::Hist::*)(std::string)) &Pythia8::Hist::title, "C++: Pythia8::Hist::title(std::string) --> void", pybind11::arg("titleIn"));
		cl.def("null", (void (Pythia8::Hist::*)()) &Pythia8::Hist::null, "C++: Pythia8::Hist::null() --> void");
		cl.def("fill", [](Pythia8::Hist &o, double const & a0) -> void { return o.fill(a0); }, "", pybind11::arg("x"));
		cl.def("fill", (void (Pythia8::Hist::*)(double, double)) &Pythia8::Hist::fill, "C++: Pythia8::Hist::fill(double, double) --> void", pybind11::arg("x"), pybind11::arg("w"));
		cl.def("table", [](Pythia8::Hist const &o) -> void { return o.table(); }, "");
		cl.def("table", [](Pythia8::Hist const &o, class std::basic_ostream<char> & a0) -> void { return o.table(a0); }, "", pybind11::arg("os"));
		cl.def("table", [](Pythia8::Hist const &o, class std::basic_ostream<char> & a0, bool const & a1) -> void { return o.table(a0, a1); }, "", pybind11::arg("os"), pybind11::arg("printOverUnder"));
		cl.def("table", [](Pythia8::Hist const &o, class std::basic_ostream<char> & a0, bool const & a1, bool const & a2) -> void { return o.table(a0, a1, a2); }, "", pybind11::arg("os"), pybind11::arg("printOverUnder"), pybind11::arg("xMidBin"));
		cl.def("table", (void (Pythia8::Hist::*)(std::ostream &, bool, bool, bool) const) &Pythia8::Hist::table, "C++: Pythia8::Hist::table(std::ostream &, bool, bool, bool) const --> void", pybind11::arg("os"), pybind11::arg("printOverUnder"), pybind11::arg("xMidBin"), pybind11::arg("printError"));
		cl.def("table", [](Pythia8::Hist const &o, class std::basic_string<char> const & a0) -> void { return o.table(a0); }, "", pybind11::arg("fileName"));
		cl.def("table", [](Pythia8::Hist const &o, class std::basic_string<char> const & a0, bool const & a1) -> void { return o.table(a0, a1); }, "", pybind11::arg("fileName"), pybind11::arg("printOverUnder"));
		cl.def("table", [](Pythia8::Hist const &o, class std::basic_string<char> const & a0, bool const & a1, bool const & a2) -> void { return o.table(a0, a1, a2); }, "", pybind11::arg("fileName"), pybind11::arg("printOverUnder"), pybind11::arg("xMidBin"));
		cl.def("table", (void (Pythia8::Hist::*)(std::string, bool, bool, bool) const) &Pythia8::Hist::table, "C++: Pythia8::Hist::table(std::string, bool, bool, bool) const --> void", pybind11::arg("fileName"), pybind11::arg("printOverUnder"), pybind11::arg("xMidBin"), pybind11::arg("printError"));
		cl.def("rivetTable", [](Pythia8::Hist const &o) -> void { return o.rivetTable(); }, "");
		cl.def("rivetTable", [](Pythia8::Hist const &o, class std::basic_ostream<char> & a0) -> void { return o.rivetTable(a0); }, "", pybind11::arg("os"));
		cl.def("rivetTable", (void (Pythia8::Hist::*)(std::ostream &, bool) const) &Pythia8::Hist::rivetTable, "C++: Pythia8::Hist::rivetTable(std::ostream &, bool) const --> void", pybind11::arg("os"), pybind11::arg("printError"));
		cl.def("rivetTable", [](Pythia8::Hist const &o, class std::basic_string<char> const & a0) -> void { return o.rivetTable(a0); }, "", pybind11::arg("fileName"));
		cl.def("rivetTable", (void (Pythia8::Hist::*)(std::string, bool) const) &Pythia8::Hist::rivetTable, "C++: Pythia8::Hist::rivetTable(std::string, bool) const --> void", pybind11::arg("fileName"), pybind11::arg("printError"));
		cl.def("pyplotTable", [](Pythia8::Hist const &o) -> void { return o.pyplotTable(); }, "");
		cl.def("pyplotTable", [](Pythia8::Hist const &o, class std::basic_ostream<char> & a0) -> void { return o.pyplotTable(a0); }, "", pybind11::arg("os"));
		cl.def("pyplotTable", [](Pythia8::Hist const &o, class std::basic_ostream<char> & a0, bool const & a1) -> void { return o.pyplotTable(a0, a1); }, "", pybind11::arg("os"), pybind11::arg("isHist"));
		cl.def("pyplotTable", (void (Pythia8::Hist::*)(std::ostream &, bool, bool) const) &Pythia8::Hist::pyplotTable, "C++: Pythia8::Hist::pyplotTable(std::ostream &, bool, bool) const --> void", pybind11::arg("os"), pybind11::arg("isHist"), pybind11::arg("printError"));
		cl.def("pyplotTable", [](Pythia8::Hist const &o, class std::basic_string<char> const & a0) -> void { return o.pyplotTable(a0); }, "", pybind11::arg("fileName"));
		cl.def("pyplotTable", [](Pythia8::Hist const &o, class std::basic_string<char> const & a0, bool const & a1) -> void { return o.pyplotTable(a0, a1); }, "", pybind11::arg("fileName"), pybind11::arg("isHist"));
		cl.def("pyplotTable", (void (Pythia8::Hist::*)(std::string, bool, bool) const) &Pythia8::Hist::pyplotTable, "C++: Pythia8::Hist::pyplotTable(std::string, bool, bool) const --> void", pybind11::arg("fileName"), pybind11::arg("isHist"), pybind11::arg("printError"));
		cl.def("fillTable", [](Pythia8::Hist &o) -> void { return o.fillTable(); }, "");
		cl.def("fillTable", (void (Pythia8::Hist::*)(class std::basic_istream<char> &)) &Pythia8::Hist::fillTable, "C++: Pythia8::Hist::fillTable(class std::basic_istream<char> &) --> void", pybind11::arg("is"));
		cl.def("fillTable", (void (Pythia8::Hist::*)(std::string)) &Pythia8::Hist::fillTable, "C++: Pythia8::Hist::fillTable(std::string) --> void", pybind11::arg("fileName"));
		cl.def("getTitle", (std::string (Pythia8::Hist::*)() const) &Pythia8::Hist::getTitle, "C++: Pythia8::Hist::getTitle() const --> std::string");
		cl.def("getBinNumber", (int (Pythia8::Hist::*)() const) &Pythia8::Hist::getBinNumber, "C++: Pythia8::Hist::getBinNumber() const --> int");
		cl.def("getNonFinite", (int (Pythia8::Hist::*)() const) &Pythia8::Hist::getNonFinite, "C++: Pythia8::Hist::getNonFinite() const --> int");
		cl.def("getLinX", (bool (Pythia8::Hist::*)() const) &Pythia8::Hist::getLinX, "C++: Pythia8::Hist::getLinX() const --> bool");
		cl.def("getXMin", (double (Pythia8::Hist::*)() const) &Pythia8::Hist::getXMin, "C++: Pythia8::Hist::getXMin() const --> double");
		cl.def("getXMax", (double (Pythia8::Hist::*)() const) &Pythia8::Hist::getXMax, "C++: Pythia8::Hist::getXMax() const --> double");
		cl.def("getYMin", (double (Pythia8::Hist::*)() const) &Pythia8::Hist::getYMin, "C++: Pythia8::Hist::getYMin() const --> double");
		cl.def("getYMax", (double (Pythia8::Hist::*)() const) &Pythia8::Hist::getYMax, "C++: Pythia8::Hist::getYMax() const --> double");
		cl.def("getYAbsMin", (double (Pythia8::Hist::*)() const) &Pythia8::Hist::getYAbsMin, "C++: Pythia8::Hist::getYAbsMin() const --> double");
		cl.def("getXMean", [](Pythia8::Hist const &o) -> double { return o.getXMean(); }, "");
		cl.def("getXMean", (double (Pythia8::Hist::*)(bool) const) &Pythia8::Hist::getXMean, "C++: Pythia8::Hist::getXMean(bool) const --> double", pybind11::arg("unbinned"));
		cl.def("getXMeanErr", [](Pythia8::Hist const &o) -> double { return o.getXMeanErr(); }, "");
		cl.def("getXMeanErr", (double (Pythia8::Hist::*)(bool) const) &Pythia8::Hist::getXMeanErr, "C++: Pythia8::Hist::getXMeanErr(bool) const --> double", pybind11::arg("unbinned"));
		cl.def("getXMedian", [](Pythia8::Hist const &o) -> double { return o.getXMedian(); }, "");
		cl.def("getXMedian", (double (Pythia8::Hist::*)(bool) const) &Pythia8::Hist::getXMedian, "C++: Pythia8::Hist::getXMedian(bool) const --> double", pybind11::arg("includeOverUnder"));
		cl.def("getXMedianErr", [](Pythia8::Hist const &o) -> double { return o.getXMedianErr(); }, "");
		cl.def("getXMedianErr", (double (Pythia8::Hist::*)(bool) const) &Pythia8::Hist::getXMedianErr, "C++: Pythia8::Hist::getXMedianErr(bool) const --> double", pybind11::arg("unbinned"));
		cl.def("getYMean", (double (Pythia8::Hist::*)() const) &Pythia8::Hist::getYMean, "C++: Pythia8::Hist::getYMean() const --> double");
		cl.def("getXRMN", [](Pythia8::Hist const &o) -> double { return o.getXRMN(); }, "");
		cl.def("getXRMN", [](Pythia8::Hist const &o, int const & a0) -> double { return o.getXRMN(a0); }, "", pybind11::arg("n"));
		cl.def("getXRMN", (double (Pythia8::Hist::*)(int, bool) const) &Pythia8::Hist::getXRMN, "C++: Pythia8::Hist::getXRMN(int, bool) const --> double", pybind11::arg("n"), pybind11::arg("unbinned"));
		cl.def("getXRMS", [](Pythia8::Hist const &o) -> double { return o.getXRMS(); }, "");
		cl.def("getXRMS", (double (Pythia8::Hist::*)(bool) const) &Pythia8::Hist::getXRMS, "C++: Pythia8::Hist::getXRMS(bool) const --> double", pybind11::arg("unbinned"));
		cl.def("getXRMNErr", [](Pythia8::Hist const &o) -> double { return o.getXRMNErr(); }, "");
		cl.def("getXRMNErr", [](Pythia8::Hist const &o, int const & a0) -> double { return o.getXRMNErr(a0); }, "", pybind11::arg("n"));
		cl.def("getXRMNErr", (double (Pythia8::Hist::*)(int, bool) const) &Pythia8::Hist::getXRMNErr, "C++: Pythia8::Hist::getXRMNErr(int, bool) const --> double", pybind11::arg("n"), pybind11::arg("unbinned"));
		cl.def("getXRMSErr", [](Pythia8::Hist const &o) -> double { return o.getXRMSErr(); }, "");
		cl.def("getXRMSErr", (double (Pythia8::Hist::*)(bool) const) &Pythia8::Hist::getXRMSErr, "C++: Pythia8::Hist::getXRMSErr(bool) const --> double", pybind11::arg("unbinned"));
		cl.def("getBinContent", (double (Pythia8::Hist::*)(int) const) &Pythia8::Hist::getBinContent, "C++: Pythia8::Hist::getBinContent(int) const --> double", pybind11::arg("iBin"));
		cl.def("getBinEdge", (double (Pythia8::Hist::*)(int) const) &Pythia8::Hist::getBinEdge, "C++: Pythia8::Hist::getBinEdge(int) const --> double", pybind11::arg("iBin"));
		cl.def("getBinWidth", [](Pythia8::Hist const &o) -> double { return o.getBinWidth(); }, "");
		cl.def("getBinWidth", (double (Pythia8::Hist::*)(int) const) &Pythia8::Hist::getBinWidth, "C++: Pythia8::Hist::getBinWidth(int) const --> double", pybind11::arg("iBin"));
		cl.def("getBinContents", (class std::vector<double, class std::allocator<double> > (Pythia8::Hist::*)() const) &Pythia8::Hist::getBinContents, "C++: Pythia8::Hist::getBinContents() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("getBinEdges", (class std::vector<double, class std::allocator<double> > (Pythia8::Hist::*)() const) &Pythia8::Hist::getBinEdges, "C++: Pythia8::Hist::getBinEdges() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("getEntries", [](Pythia8::Hist const &o) -> int { return o.getEntries(); }, "");
		cl.def("getEntries", (int (Pythia8::Hist::*)(bool) const) &Pythia8::Hist::getEntries, "C++: Pythia8::Hist::getEntries(bool) const --> int", pybind11::arg("alsoNonFinite"));
		cl.def("getWeightSum", [](Pythia8::Hist const &o) -> double { return o.getWeightSum(); }, "");
		cl.def("getWeightSum", (double (Pythia8::Hist::*)(bool) const) &Pythia8::Hist::getWeightSum, "C++: Pythia8::Hist::getWeightSum(bool) const --> double", pybind11::arg("alsoOverUnder"));
		cl.def("getNEffective", (double (Pythia8::Hist::*)() const) &Pythia8::Hist::getNEffective, "C++: Pythia8::Hist::getNEffective() const --> double");
		cl.def("sameSize", (bool (Pythia8::Hist::*)(const class Pythia8::Hist &) const) &Pythia8::Hist::sameSize, "C++: Pythia8::Hist::sameSize(const class Pythia8::Hist &) const --> bool", pybind11::arg("h"));
		cl.def("takeFunc", (void (Pythia8::Hist::*)(class std::function<double (double)>)) &Pythia8::Hist::takeFunc, "C++: Pythia8::Hist::takeFunc(class std::function<double (double)>) --> void", pybind11::arg("func"));
		cl.def("takeLog", [](Pythia8::Hist &o) -> void { return o.takeLog(); }, "");
		cl.def("takeLog", (void (Pythia8::Hist::*)(bool)) &Pythia8::Hist::takeLog, "C++: Pythia8::Hist::takeLog(bool) --> void", pybind11::arg("tenLog"));
		cl.def("takeSqrt", (void (Pythia8::Hist::*)()) &Pythia8::Hist::takeSqrt, "C++: Pythia8::Hist::takeSqrt() --> void");
		cl.def("normalize", [](Pythia8::Hist &o) -> void { return o.normalize(); }, "");
		cl.def("normalize", [](Pythia8::Hist &o, double const & a0) -> void { return o.normalize(a0); }, "", pybind11::arg("f"));
		cl.def("normalize", (void (Pythia8::Hist::*)(double, bool)) &Pythia8::Hist::normalize, "C++: Pythia8::Hist::normalize(double, bool) --> void", pybind11::arg("f"), pybind11::arg("overflow"));
		cl.def("normalizeIntegral", [](Pythia8::Hist &o) -> void { return o.normalizeIntegral(); }, "");
		cl.def("normalizeIntegral", [](Pythia8::Hist &o, double const & a0) -> void { return o.normalizeIntegral(a0); }, "", pybind11::arg("f"));
		cl.def("normalizeIntegral", (void (Pythia8::Hist::*)(double, bool)) &Pythia8::Hist::normalizeIntegral, "C++: Pythia8::Hist::normalizeIntegral(double, bool) --> void", pybind11::arg("f"), pybind11::arg("overflow"));
		cl.def("normalizeSpectrum", (void (Pythia8::Hist::*)(double)) &Pythia8::Hist::normalizeSpectrum, "C++: Pythia8::Hist::normalizeSpectrum(double) --> void", pybind11::arg("wtSum"));
		cl.def("__iadd__", (class Pythia8::Hist & (Pythia8::Hist::*)(const class Pythia8::Hist &)) &Pythia8::Hist::operator+=, "C++: Pythia8::Hist::operator+=(const class Pythia8::Hist &) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("h"));
		cl.def("__isub__", (class Pythia8::Hist & (Pythia8::Hist::*)(const class Pythia8::Hist &)) &Pythia8::Hist::operator-=, "C++: Pythia8::Hist::operator-=(const class Pythia8::Hist &) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("h"));
		cl.def("__imul__", (class Pythia8::Hist & (Pythia8::Hist::*)(const class Pythia8::Hist &)) &Pythia8::Hist::operator*=, "C++: Pythia8::Hist::operator*=(const class Pythia8::Hist &) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("h"));
		cl.def("__idiv__", (class Pythia8::Hist & (Pythia8::Hist::*)(const class Pythia8::Hist &)) &Pythia8::Hist::operator/=, "C++: Pythia8::Hist::operator/=(const class Pythia8::Hist &) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("h"));
		cl.def("__iadd__", (class Pythia8::Hist & (Pythia8::Hist::*)(double)) &Pythia8::Hist::operator+=, "C++: Pythia8::Hist::operator+=(double) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("f"));
		cl.def("__isub__", (class Pythia8::Hist & (Pythia8::Hist::*)(double)) &Pythia8::Hist::operator-=, "C++: Pythia8::Hist::operator-=(double) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("f"));
		cl.def("__imul__", (class Pythia8::Hist & (Pythia8::Hist::*)(double)) &Pythia8::Hist::operator*=, "C++: Pythia8::Hist::operator*=(double) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("f"));
		cl.def("__idiv__", (class Pythia8::Hist & (Pythia8::Hist::*)(double)) &Pythia8::Hist::operator/=, "C++: Pythia8::Hist::operator/=(double) --> class Pythia8::Hist &", pybind11::return_value_policy::reference, pybind11::arg("f"));
		cl.def("__add__", (class Pythia8::Hist (Pythia8::Hist::*)(double) const) &Pythia8::Hist::operator+, "C++: Pythia8::Hist::operator+(double) const --> class Pythia8::Hist", pybind11::arg("f"));
		cl.def("__add__", (class Pythia8::Hist (Pythia8::Hist::*)(const class Pythia8::Hist &) const) &Pythia8::Hist::operator+, "C++: Pythia8::Hist::operator+(const class Pythia8::Hist &) const --> class Pythia8::Hist", pybind11::arg("h2"));
		cl.def("__sub__", (class Pythia8::Hist (Pythia8::Hist::*)(double) const) &Pythia8::Hist::operator-, "C++: Pythia8::Hist::operator-(double) const --> class Pythia8::Hist", pybind11::arg("f"));
		cl.def("__sub__", (class Pythia8::Hist (Pythia8::Hist::*)(const class Pythia8::Hist &) const) &Pythia8::Hist::operator-, "C++: Pythia8::Hist::operator-(const class Pythia8::Hist &) const --> class Pythia8::Hist", pybind11::arg("h2"));
		cl.def("__mul__", (class Pythia8::Hist (Pythia8::Hist::*)(double) const) &Pythia8::Hist::operator*, "C++: Pythia8::Hist::operator*(double) const --> class Pythia8::Hist", pybind11::arg("f"));
		cl.def("__mul__", (class Pythia8::Hist (Pythia8::Hist::*)(const class Pythia8::Hist &) const) &Pythia8::Hist::operator*, "C++: Pythia8::Hist::operator*(const class Pythia8::Hist &) const --> class Pythia8::Hist", pybind11::arg("h2"));
		cl.def("__div__", (class Pythia8::Hist (Pythia8::Hist::*)(double) const) &Pythia8::Hist::operator/, "C++: Pythia8::Hist::operator/(double) const --> class Pythia8::Hist", pybind11::arg("f"));
		cl.def("__div__", (class Pythia8::Hist (Pythia8::Hist::*)(const class Pythia8::Hist &) const) &Pythia8::Hist::operator/, "C++: Pythia8::Hist::operator/(const class Pythia8::Hist &) const --> class Pythia8::Hist", pybind11::arg("h2"));

		cl.def("__str__", [](Pythia8::Hist const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
}
