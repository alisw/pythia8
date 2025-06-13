#include <Pythia8/LHEF3.h>
#include <Pythia8/Streams.h>
#include <cwchar>
#include <ios>
#include <iterator>
#include <memory>
#include <ostream>
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

void bind_Pythia8_Streams(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::DummyForStreams file:Pythia8/Streams.h line:118
		pybind11::class_<Pythia8::DummyForStreams, std::shared_ptr<Pythia8::DummyForStreams>> cl(M("Pythia8"), "DummyForStreams", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::DummyForStreams(); } ) );
		cl.def("xtox", (double (Pythia8::DummyForStreams::*)(double)) &Pythia8::DummyForStreams::xtox, "C++: Pythia8::DummyForStreams::xtox(double) --> double", pybind11::arg("x"));
	}
	{ // Pythia8::XMLTag file:Pythia8/LHEF3.h line:27
		pybind11::class_<Pythia8::XMLTag, std::shared_ptr<Pythia8::XMLTag>> cl(M("Pythia8"), "XMLTag", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::XMLTag(); } ) );
		cl.def( pybind11::init( [](Pythia8::XMLTag const &o){ return new Pythia8::XMLTag(o); } ) );
		cl.def_readwrite("name", &Pythia8::XMLTag::name);
		cl.def_readwrite("attr", &Pythia8::XMLTag::attr);
		cl.def_readwrite("tags", &Pythia8::XMLTag::tags);
		cl.def_readwrite("contents", &Pythia8::XMLTag::contents);
		cl.def("getattr", (bool (Pythia8::XMLTag::*)(std::string, double &) const) &Pythia8::XMLTag::getattr, "C++: Pythia8::XMLTag::getattr(std::string, double &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (Pythia8::XMLTag::*)(std::string, bool &) const) &Pythia8::XMLTag::getattr, "C++: Pythia8::XMLTag::getattr(std::string, bool &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (Pythia8::XMLTag::*)(std::string, long &) const) &Pythia8::XMLTag::getattr, "C++: Pythia8::XMLTag::getattr(std::string, long &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (Pythia8::XMLTag::*)(std::string, int &) const) &Pythia8::XMLTag::getattr, "C++: Pythia8::XMLTag::getattr(std::string, int &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def("getattr", (bool (Pythia8::XMLTag::*)(std::string, std::string &) const) &Pythia8::XMLTag::getattr, "C++: Pythia8::XMLTag::getattr(std::string, std::string &) const --> bool", pybind11::arg("n"), pybind11::arg("v"));
		cl.def_static("findXMLTags", [](class std::basic_string<char> const & a0) -> std::vector<struct Pythia8::XMLTag *, class std::allocator<struct Pythia8::XMLTag *> > { return Pythia8::XMLTag::findXMLTags(a0); }, "", pybind11::arg("str"));
		cl.def_static("findXMLTags", (class std::vector<struct Pythia8::XMLTag *, class std::allocator<struct Pythia8::XMLTag *> > (*)(std::string, std::string *)) &Pythia8::XMLTag::findXMLTags, "C++: Pythia8::XMLTag::findXMLTags(std::string, std::string *) --> class std::vector<struct Pythia8::XMLTag *, class std::allocator<struct Pythia8::XMLTag *> >", pybind11::arg("str"), pybind11::arg("leftover"));
		cl.def("list", (void (Pythia8::XMLTag::*)(std::ostream &) const) &Pythia8::XMLTag::list, "C++: Pythia8::XMLTag::list(std::ostream &) const --> void", pybind11::arg("os"));
	}
	{ // Pythia8::LHAweights file:Pythia8/LHEF3.h line:245
		pybind11::class_<Pythia8::LHAweights, std::shared_ptr<Pythia8::LHAweights>> cl(M("Pythia8"), "LHAweights", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAweights(); } ) );
		cl.def( pybind11::init<const struct Pythia8::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](Pythia8::LHAweights const &o){ return new Pythia8::LHAweights(o); } ) );
		cl.def_readwrite("weights", &Pythia8::LHAweights::weights);
		cl.def_readwrite("attributes", &Pythia8::LHAweights::attributes);
		cl.def_readwrite("contents", &Pythia8::LHAweights::contents);
		cl.def("list", (void (Pythia8::LHAweights::*)(std::ostream &) const) &Pythia8::LHAweights::list, "C++: Pythia8::LHAweights::list(std::ostream &) const --> void", pybind11::arg("file"));
		cl.def("clear", (void (Pythia8::LHAweights::*)()) &Pythia8::LHAweights::clear, "C++: Pythia8::LHAweights::clear() --> void");
		cl.def("size", (int (Pythia8::LHAweights::*)() const) &Pythia8::LHAweights::size, "C++: Pythia8::LHAweights::size() const --> int");
	}
	{ // Pythia8::LHAscales file:Pythia8/LHEF3.h line:281
		pybind11::class_<Pythia8::LHAscales, std::shared_ptr<Pythia8::LHAscales>> cl(M("Pythia8"), "LHAscales", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAscales(); } ), "doc" );
		cl.def( pybind11::init<double>(), pybind11::arg("defscale") );

		cl.def( pybind11::init( [](const struct Pythia8::XMLTag & a0){ return new Pythia8::LHAscales(a0); } ), "doc" , pybind11::arg("tag"));
		cl.def( pybind11::init<const struct Pythia8::XMLTag &, double>(), pybind11::arg("tag"), pybind11::arg("defscale") );

		cl.def( pybind11::init( [](Pythia8::LHAscales const &o){ return new Pythia8::LHAscales(o); } ) );
		cl.def_readwrite("muf", &Pythia8::LHAscales::muf);
		cl.def_readwrite("mur", &Pythia8::LHAscales::mur);
		cl.def_readwrite("mups", &Pythia8::LHAscales::mups);
		cl.def_readwrite("attributes", &Pythia8::LHAscales::attributes);
		cl.def_readwrite("SCALUP", &Pythia8::LHAscales::SCALUP);
		cl.def_readwrite("contents", &Pythia8::LHAscales::contents);
		cl.def("list", (void (Pythia8::LHAscales::*)(std::ostream &) const) &Pythia8::LHAscales::list, "C++: Pythia8::LHAscales::list(std::ostream &) const --> void", pybind11::arg("file"));
		cl.def("clear", (void (Pythia8::LHAscales::*)()) &Pythia8::LHAscales::clear, "C++: Pythia8::LHAscales::clear() --> void");
	}
	{ // Pythia8::LHAgenerator file:Pythia8/LHEF3.h line:325
		pybind11::class_<Pythia8::LHAgenerator, std::shared_ptr<Pythia8::LHAgenerator>> cl(M("Pythia8"), "LHAgenerator", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAgenerator(); } ) );
		cl.def( pybind11::init( [](const struct Pythia8::XMLTag & a0){ return new Pythia8::LHAgenerator(a0); } ), "doc" , pybind11::arg("tag"));
		cl.def( pybind11::init<const struct Pythia8::XMLTag &, std::string>(), pybind11::arg("tag"), pybind11::arg("defname") );

		cl.def( pybind11::init( [](Pythia8::LHAgenerator const &o){ return new Pythia8::LHAgenerator(o); } ) );
		cl.def_readwrite("name", &Pythia8::LHAgenerator::name);
		cl.def_readwrite("version", &Pythia8::LHAgenerator::version);
		cl.def_readwrite("attributes", &Pythia8::LHAgenerator::attributes);
		cl.def_readwrite("contents", &Pythia8::LHAgenerator::contents);
		cl.def("list", (void (Pythia8::LHAgenerator::*)(std::ostream &) const) &Pythia8::LHAgenerator::list, "C++: Pythia8::LHAgenerator::list(std::ostream &) const --> void", pybind11::arg("file"));
		cl.def("clear", (void (Pythia8::LHAgenerator::*)()) &Pythia8::LHAgenerator::clear, "C++: Pythia8::LHAgenerator::clear() --> void");
		cl.def("assign", (struct Pythia8::LHAgenerator & (Pythia8::LHAgenerator::*)(const struct Pythia8::LHAgenerator &)) &Pythia8::LHAgenerator::operator=, "C++: Pythia8::LHAgenerator::operator=(const struct Pythia8::LHAgenerator &) --> struct Pythia8::LHAgenerator &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::LHAwgt file:Pythia8/LHEF3.h line:363
		pybind11::class_<Pythia8::LHAwgt, std::shared_ptr<Pythia8::LHAwgt>> cl(M("Pythia8"), "LHAwgt", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAwgt(); } ), "doc" );
		cl.def( pybind11::init<double>(), pybind11::arg("defwgt") );

		cl.def( pybind11::init( [](const struct Pythia8::XMLTag & a0){ return new Pythia8::LHAwgt(a0); } ), "doc" , pybind11::arg("tag"));
		cl.def( pybind11::init<const struct Pythia8::XMLTag &, double>(), pybind11::arg("tag"), pybind11::arg("defwgt") );

		cl.def( pybind11::init( [](Pythia8::LHAwgt const &o){ return new Pythia8::LHAwgt(o); } ) );
		cl.def_readwrite("id", &Pythia8::LHAwgt::id);
		cl.def_readwrite("attributes", &Pythia8::LHAwgt::attributes);
		cl.def_readwrite("contents", &Pythia8::LHAwgt::contents);
		cl.def("list", (void (Pythia8::LHAwgt::*)(std::ostream &) const) &Pythia8::LHAwgt::list, "C++: Pythia8::LHAwgt::list(std::ostream &) const --> void", pybind11::arg("file"));
		cl.def("clear", (void (Pythia8::LHAwgt::*)()) &Pythia8::LHAwgt::clear, "C++: Pythia8::LHAwgt::clear() --> void");
	}
	{ // Pythia8::LHAweight file:Pythia8/LHEF3.h line:397
		pybind11::class_<Pythia8::LHAweight, std::shared_ptr<Pythia8::LHAweight>> cl(M("Pythia8"), "LHAweight", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAweight(); } ), "doc" );
		cl.def( pybind11::init<std::string>(), pybind11::arg("defname") );

		cl.def( pybind11::init( [](const struct Pythia8::XMLTag & a0){ return new Pythia8::LHAweight(a0); } ), "doc" , pybind11::arg("tag"));
		cl.def( pybind11::init<const struct Pythia8::XMLTag &, std::string>(), pybind11::arg("tag"), pybind11::arg("defname") );

		cl.def( pybind11::init( [](Pythia8::LHAweight const &o){ return new Pythia8::LHAweight(o); } ) );
		cl.def_readwrite("id", &Pythia8::LHAweight::id);
		cl.def_readwrite("attributes", &Pythia8::LHAweight::attributes);
		cl.def_readwrite("contents", &Pythia8::LHAweight::contents);
		cl.def("list", (void (Pythia8::LHAweight::*)(std::ostream &) const) &Pythia8::LHAweight::list, "C++: Pythia8::LHAweight::list(std::ostream &) const --> void", pybind11::arg("file"));
		cl.def("clear", (void (Pythia8::LHAweight::*)()) &Pythia8::LHAweight::clear, "C++: Pythia8::LHAweight::clear() --> void");
	}
	{ // Pythia8::LHAweightgroup file:Pythia8/LHEF3.h line:431
		pybind11::class_<Pythia8::LHAweightgroup, std::shared_ptr<Pythia8::LHAweightgroup>> cl(M("Pythia8"), "LHAweightgroup", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAweightgroup(); } ) );
		cl.def( pybind11::init<const struct Pythia8::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](Pythia8::LHAweightgroup const &o){ return new Pythia8::LHAweightgroup(o); } ) );
		cl.def_readwrite("contents", &Pythia8::LHAweightgroup::contents);
		cl.def_readwrite("name", &Pythia8::LHAweightgroup::name);
		cl.def_readwrite("weights", &Pythia8::LHAweightgroup::weights);
		cl.def_readwrite("weightsKeys", &Pythia8::LHAweightgroup::weightsKeys);
		cl.def_readwrite("attributes", &Pythia8::LHAweightgroup::attributes);
		cl.def("list", (void (Pythia8::LHAweightgroup::*)(std::ostream &) const) &Pythia8::LHAweightgroup::list, "C++: Pythia8::LHAweightgroup::list(std::ostream &) const --> void", pybind11::arg("file"));
		cl.def("clear", (void (Pythia8::LHAweightgroup::*)()) &Pythia8::LHAweightgroup::clear, "C++: Pythia8::LHAweightgroup::clear() --> void");
		cl.def("size", (int (Pythia8::LHAweightgroup::*)() const) &Pythia8::LHAweightgroup::size, "C++: Pythia8::LHAweightgroup::size() const --> int");
	}
	{ // Pythia8::LHArwgt file:Pythia8/LHEF3.h line:473
		pybind11::class_<Pythia8::LHArwgt, std::shared_ptr<Pythia8::LHArwgt>> cl(M("Pythia8"), "LHArwgt", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHArwgt(); } ) );
		cl.def( pybind11::init<const struct Pythia8::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](Pythia8::LHArwgt const &o){ return new Pythia8::LHArwgt(o); } ) );
		cl.def_readwrite("contents", &Pythia8::LHArwgt::contents);
		cl.def_readwrite("wgts", &Pythia8::LHArwgt::wgts);
		cl.def_readwrite("wgtsKeys", &Pythia8::LHArwgt::wgtsKeys);
		cl.def_readwrite("attributes", &Pythia8::LHArwgt::attributes);
		cl.def("list", (void (Pythia8::LHArwgt::*)(std::ostream &) const) &Pythia8::LHArwgt::list, "C++: Pythia8::LHArwgt::list(std::ostream &) const --> void", pybind11::arg("file"));
		cl.def("clear", (void (Pythia8::LHArwgt::*)()) &Pythia8::LHArwgt::clear, "C++: Pythia8::LHArwgt::clear() --> void");
		cl.def("size", (int (Pythia8::LHArwgt::*)() const) &Pythia8::LHArwgt::size, "C++: Pythia8::LHArwgt::size() const --> int");
	}
	{ // Pythia8::LHAinitrwgt file:Pythia8/LHEF3.h line:511
		pybind11::class_<Pythia8::LHAinitrwgt, std::shared_ptr<Pythia8::LHAinitrwgt>> cl(M("Pythia8"), "LHAinitrwgt", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::LHAinitrwgt(); } ) );
		cl.def( pybind11::init<const struct Pythia8::XMLTag &>(), pybind11::arg("tag") );

		cl.def( pybind11::init( [](Pythia8::LHAinitrwgt const &o){ return new Pythia8::LHAinitrwgt(o); } ) );
		cl.def_readwrite("contents", &Pythia8::LHAinitrwgt::contents);
		cl.def_readwrite("weights", &Pythia8::LHAinitrwgt::weights);
		cl.def_readwrite("weightsKeys", &Pythia8::LHAinitrwgt::weightsKeys);
		cl.def_readwrite("weightgroups", &Pythia8::LHAinitrwgt::weightgroups);
		cl.def_readwrite("weightgroupsKeys", &Pythia8::LHAinitrwgt::weightgroupsKeys);
		cl.def_readwrite("attributes", &Pythia8::LHAinitrwgt::attributes);
		cl.def("list", (void (Pythia8::LHAinitrwgt::*)(std::ostream &) const) &Pythia8::LHAinitrwgt::list, "C++: Pythia8::LHAinitrwgt::list(std::ostream &) const --> void", pybind11::arg("file"));
		cl.def("clear", (void (Pythia8::LHAinitrwgt::*)()) &Pythia8::LHAinitrwgt::clear, "C++: Pythia8::LHAinitrwgt::clear() --> void");
		cl.def("size", (int (Pythia8::LHAinitrwgt::*)() const) &Pythia8::LHAinitrwgt::size, "C++: Pythia8::LHAinitrwgt::size() const --> int");
		cl.def("sizeWeightGroups", (int (Pythia8::LHAinitrwgt::*)() const) &Pythia8::LHAinitrwgt::sizeWeightGroups, "C++: Pythia8::LHAinitrwgt::sizeWeightGroups() const --> int");
		cl.def("assign", (struct Pythia8::LHAinitrwgt & (Pythia8::LHAinitrwgt::*)(const struct Pythia8::LHAinitrwgt &)) &Pythia8::LHAinitrwgt::operator=, "C++: Pythia8::LHAinitrwgt::operator=(const struct Pythia8::LHAinitrwgt &) --> struct Pythia8::LHAinitrwgt &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HEPRUP file:Pythia8/LHEF3.h line:562
		pybind11::class_<Pythia8::HEPRUP, std::shared_ptr<Pythia8::HEPRUP>> cl(M("Pythia8"), "HEPRUP", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HEPRUP(); } ) );
		cl.def( pybind11::init( [](Pythia8::HEPRUP const &o){ return new Pythia8::HEPRUP(o); } ) );
		cl.def_readwrite("IDBMUP", &Pythia8::HEPRUP::IDBMUP);
		cl.def_readwrite("EBMUP", &Pythia8::HEPRUP::EBMUP);
		cl.def_readwrite("PDFGUP", &Pythia8::HEPRUP::PDFGUP);
		cl.def_readwrite("PDFSUP", &Pythia8::HEPRUP::PDFSUP);
		cl.def_readwrite("IDWTUP", &Pythia8::HEPRUP::IDWTUP);
		cl.def_readwrite("NPRUP", &Pythia8::HEPRUP::NPRUP);
		cl.def_readwrite("XSECUP", &Pythia8::HEPRUP::XSECUP);
		cl.def_readwrite("XERRUP", &Pythia8::HEPRUP::XERRUP);
		cl.def_readwrite("XMAXUP", &Pythia8::HEPRUP::XMAXUP);
		cl.def_readwrite("LPRUP", &Pythia8::HEPRUP::LPRUP);
		cl.def_readwrite("initrwgt", &Pythia8::HEPRUP::initrwgt);
		cl.def_readwrite("generators", &Pythia8::HEPRUP::generators);
		cl.def_readwrite("weightgroups", &Pythia8::HEPRUP::weightgroups);
		cl.def_readwrite("weights", &Pythia8::HEPRUP::weights);
		cl.def("assign", (class Pythia8::HEPRUP & (Pythia8::HEPRUP::*)(const class Pythia8::HEPRUP &)) &Pythia8::HEPRUP::operator=, "C++: Pythia8::HEPRUP::operator=(const class Pythia8::HEPRUP &) --> class Pythia8::HEPRUP &", pybind11::return_value_policy::reference, pybind11::arg("x"));
		cl.def("resize", (void (Pythia8::HEPRUP::*)(int)) &Pythia8::HEPRUP::resize, "C++: Pythia8::HEPRUP::resize(int) --> void", pybind11::arg("nrup"));
		cl.def("resize", (void (Pythia8::HEPRUP::*)()) &Pythia8::HEPRUP::resize, "C++: Pythia8::HEPRUP::resize() --> void");
		cl.def("clear", (void (Pythia8::HEPRUP::*)()) &Pythia8::HEPRUP::clear, "C++: Pythia8::HEPRUP::clear() --> void");
	}
	{ // Pythia8::HEPEUP file:Pythia8/LHEF3.h line:671
		pybind11::class_<Pythia8::HEPEUP, std::shared_ptr<Pythia8::HEPEUP>> cl(M("Pythia8"), "HEPEUP", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HEPEUP(); } ) );
		cl.def( pybind11::init( [](Pythia8::HEPEUP const &o){ return new Pythia8::HEPEUP(o); } ) );
		cl.def_readwrite("NUP", &Pythia8::HEPEUP::NUP);
		cl.def_readwrite("IDPRUP", &Pythia8::HEPEUP::IDPRUP);
		cl.def_readwrite("XWGTUP", &Pythia8::HEPEUP::XWGTUP);
		cl.def_readwrite("XPDWUP", &Pythia8::HEPEUP::XPDWUP);
		cl.def_readwrite("SCALUP", &Pythia8::HEPEUP::SCALUP);
		cl.def_readwrite("AQEDUP", &Pythia8::HEPEUP::AQEDUP);
		cl.def_readwrite("AQCDUP", &Pythia8::HEPEUP::AQCDUP);
		cl.def_readwrite("IDUP", &Pythia8::HEPEUP::IDUP);
		cl.def_readwrite("ISTUP", &Pythia8::HEPEUP::ISTUP);
		cl.def_readwrite("MOTHUP", &Pythia8::HEPEUP::MOTHUP);
		cl.def_readwrite("ICOLUP", &Pythia8::HEPEUP::ICOLUP);
		cl.def_readwrite("PUP", &Pythia8::HEPEUP::PUP);
		cl.def_readwrite("VTIMUP", &Pythia8::HEPEUP::VTIMUP);
		cl.def_readwrite("SPINUP", &Pythia8::HEPEUP::SPINUP);
		cl.def_readwrite("weights_detailed", &Pythia8::HEPEUP::weights_detailed);
		cl.def_readwrite("weights_compressed", &Pythia8::HEPEUP::weights_compressed);
		cl.def_readwrite("scalesSave", &Pythia8::HEPEUP::scalesSave);
		cl.def_readwrite("weightsSave", &Pythia8::HEPEUP::weightsSave);
		cl.def_readwrite("rwgtSave", &Pythia8::HEPEUP::rwgtSave);
		cl.def_readwrite("attributes", &Pythia8::HEPEUP::attributes);
		cl.def("setEvent", (class Pythia8::HEPEUP & (Pythia8::HEPEUP::*)(const class Pythia8::HEPEUP &)) &Pythia8::HEPEUP::setEvent, "C++: Pythia8::HEPEUP::setEvent(const class Pythia8::HEPEUP &) --> class Pythia8::HEPEUP &", pybind11::return_value_policy::reference, pybind11::arg("x"));
		cl.def("assign", (class Pythia8::HEPEUP & (Pythia8::HEPEUP::*)(const class Pythia8::HEPEUP &)) &Pythia8::HEPEUP::operator=, "C++: Pythia8::HEPEUP::operator=(const class Pythia8::HEPEUP &) --> class Pythia8::HEPEUP &", pybind11::return_value_policy::reference, pybind11::arg("x"));
		cl.def("reset", (void (Pythia8::HEPEUP::*)()) &Pythia8::HEPEUP::reset, "C++: Pythia8::HEPEUP::reset() --> void");
		cl.def("clear", (void (Pythia8::HEPEUP::*)()) &Pythia8::HEPEUP::clear, "C++: Pythia8::HEPEUP::clear() --> void");
		cl.def("resize", (void (Pythia8::HEPEUP::*)(int)) &Pythia8::HEPEUP::resize, "C++: Pythia8::HEPEUP::resize(int) --> void", pybind11::arg("nup"));
		cl.def("weight", (double (Pythia8::HEPEUP::*)() const) &Pythia8::HEPEUP::weight, "C++: Pythia8::HEPEUP::weight() const --> double");
		cl.def("resize", (void (Pythia8::HEPEUP::*)()) &Pythia8::HEPEUP::resize, "C++: Pythia8::HEPEUP::resize() --> void");
	}
}
