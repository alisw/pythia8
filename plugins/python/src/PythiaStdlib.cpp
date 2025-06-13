#include <Pythia8/Basics.h>
#include <Pythia8/PythiaStdlib.h>
#include <cwchar>
#include <ios>
#include <iterator>
#include <locale>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>

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

void bind_Pythia8_PythiaStdlib(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Pythia8::pow2(const double &) file:Pythia8/PythiaStdlib.h line:173
	M("Pythia8").def("pow2", (double (*)(const double &)) &Pythia8::pow2, "C++: Pythia8::pow2(const double &) --> double", pybind11::arg("x"));

	// Pythia8::pow3(const double &) file:Pythia8/PythiaStdlib.h line:174
	M("Pythia8").def("pow3", (double (*)(const double &)) &Pythia8::pow3, "C++: Pythia8::pow3(const double &) --> double", pybind11::arg("x"));

	// Pythia8::pow4(const double &) file:Pythia8/PythiaStdlib.h line:175
	M("Pythia8").def("pow4", (double (*)(const double &)) &Pythia8::pow4, "C++: Pythia8::pow4(const double &) --> double", pybind11::arg("x"));

	// Pythia8::pow5(const double &) file:Pythia8/PythiaStdlib.h line:176
	M("Pythia8").def("pow5", (double (*)(const double &)) &Pythia8::pow5, "C++: Pythia8::pow5(const double &) --> double", pybind11::arg("x"));

	// Pythia8::pow6(const double &) file:Pythia8/PythiaStdlib.h line:177
	M("Pythia8").def("pow6", (double (*)(const double &)) &Pythia8::pow6, "C++: Pythia8::pow6(const double &) --> double", pybind11::arg("x"));

	// Pythia8::pow7(const double &) file:Pythia8/PythiaStdlib.h line:178
	M("Pythia8").def("pow7", (double (*)(const double &)) &Pythia8::pow7, "C++: Pythia8::pow7(const double &) --> double", pybind11::arg("x"));

	// Pythia8::pow8(const double &) file:Pythia8/PythiaStdlib.h line:179
	M("Pythia8").def("pow8", (double (*)(const double &)) &Pythia8::pow8, "C++: Pythia8::pow8(const double &) --> double", pybind11::arg("x"));

	// Pythia8::sqrtpos(const double &) file:Pythia8/PythiaStdlib.h line:182
	M("Pythia8").def("sqrtpos", (double (*)(const double &)) &Pythia8::sqrtpos, "C++: Pythia8::sqrtpos(const double &) --> double", pybind11::arg("x"));

	// Pythia8::sqrtnan(const double &) file:Pythia8/PythiaStdlib.h line:185
	M("Pythia8").def("sqrtnan", (double (*)(const double &)) &Pythia8::sqrtnan, "C++: Pythia8::sqrtnan(const double &) --> double", pybind11::arg("x"));

	// Pythia8::clamp(const double &, const double &, const double &) file:Pythia8/PythiaStdlib.h line:189
	M("Pythia8").def("clamp", (double (*)(const double &, const double &, const double &)) &Pythia8::clamp, "C++: Pythia8::clamp(const double &, const double &, const double &) --> double", pybind11::arg("x"), pybind11::arg("xmin"), pybind11::arg("xmax"));

	// Pythia8::toLower(const std::string &, bool) file:Pythia8/PythiaStdlib.h line:194
	M("Pythia8").def("toLower", [](const class std::basic_string<char> & a0) -> std::string { return Pythia8::toLower(a0); }, "", pybind11::arg("name"));
	M("Pythia8").def("toLower", (std::string (*)(const std::string &, bool)) &Pythia8::toLower, "C++: Pythia8::toLower(const std::string &, bool) --> std::string", pybind11::arg("name"), pybind11::arg("trim"));

	// Pythia8::toLowerRep(std::string &, bool) file:Pythia8/PythiaStdlib.h line:197
	M("Pythia8").def("toLowerRep", [](class std::basic_string<char> & a0) -> void { return Pythia8::toLowerRep(a0); }, "", pybind11::arg("name"));
	M("Pythia8").def("toLowerRep", (void (*)(std::string &, bool)) &Pythia8::toLowerRep, "C++: Pythia8::toLowerRep(std::string &, bool) --> void", pybind11::arg("name"), pybind11::arg("trim"));

	// Pythia8::trimString(const std::string &) file:Pythia8/PythiaStdlib.h line:201
	M("Pythia8").def("trimString", (std::string (*)(const std::string &)) &Pythia8::trimString, "C++: Pythia8::trimString(const std::string &) --> std::string", pybind11::arg("name"));

	// Pythia8::trimStringRep(std::string &) file:Pythia8/PythiaStdlib.h line:204
	M("Pythia8").def("trimStringRep", (void (*)(std::string &)) &Pythia8::trimStringRep, "C++: Pythia8::trimStringRep(std::string &) --> void", pybind11::arg("name"));

	// Pythia8::toString(bool) file:Pythia8/PythiaStdlib.h line:208
	M("Pythia8").def("toString", (std::string (*)(bool)) &Pythia8::toString, "C++: Pythia8::toString(bool) --> std::string", pybind11::arg("val"));

	// Pythia8::toString(int) file:Pythia8/PythiaStdlib.h line:211
	M("Pythia8").def("toString", (std::string (*)(int)) &Pythia8::toString, "C++: Pythia8::toString(int) --> std::string", pybind11::arg("val"));

	// Pythia8::toString(double) file:Pythia8/PythiaStdlib.h line:214
	M("Pythia8").def("toString", (std::string (*)(double)) &Pythia8::toString, "C++: Pythia8::toString(double) --> std::string", pybind11::arg("val"));

	// Pythia8::methodName(const std::string &, bool) file:Pythia8/PythiaStdlib.h line:281
	M("Pythia8").def("methodName", [](const class std::basic_string<char> & a0) -> std::string { return Pythia8::methodName(a0); }, "", pybind11::arg("prettyFunction"));
	M("Pythia8").def("methodName", (std::string (*)(const std::string &, bool)) &Pythia8::methodName, "C++: Pythia8::methodName(const std::string &, bool) --> std::string", pybind11::arg("prettyFunction"), pybind11::arg("withNamespace"));

	{ // Pythia8::Vec4 file:Pythia8/Basics.h line:32
		pybind11::class_<Pythia8::Vec4, std::shared_ptr<Pythia8::Vec4>> cl(M("Pythia8"), "Vec4", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::Vec4(); } ), "doc" );
		cl.def( pybind11::init( [](double const & a0){ return new Pythia8::Vec4(a0); } ), "doc" , pybind11::arg("xIn"));
		cl.def( pybind11::init( [](double const & a0, double const & a1){ return new Pythia8::Vec4(a0, a1); } ), "doc" , pybind11::arg("xIn"), pybind11::arg("yIn"));
		cl.def( pybind11::init( [](double const & a0, double const & a1, double const & a2){ return new Pythia8::Vec4(a0, a1, a2); } ), "doc" , pybind11::arg("xIn"), pybind11::arg("yIn"), pybind11::arg("zIn"));
		cl.def( pybind11::init<double, double, double, double>(), pybind11::arg("xIn"), pybind11::arg("yIn"), pybind11::arg("zIn"), pybind11::arg("tIn") );

		cl.def( pybind11::init( [](Pythia8::Vec4 const &o){ return new Pythia8::Vec4(o); } ) );
		cl.def("assign", (class Pythia8::Vec4 & (Pythia8::Vec4::*)(const class Pythia8::Vec4 &)) &Pythia8::Vec4::operator=, "C++: Pythia8::Vec4::operator=(const class Pythia8::Vec4 &) --> class Pythia8::Vec4 &", pybind11::return_value_policy::reference, pybind11::arg("v"));
		cl.def("assign", (class Pythia8::Vec4 & (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::operator=, "C++: Pythia8::Vec4::operator=(double) --> class Pythia8::Vec4 &", pybind11::return_value_policy::reference, pybind11::arg("value"));
		cl.def("reset", (void (Pythia8::Vec4::*)()) &Pythia8::Vec4::reset, "C++: Pythia8::Vec4::reset() --> void");
		cl.def("p", (void (Pythia8::Vec4::*)(double, double, double, double)) &Pythia8::Vec4::p, "C++: Pythia8::Vec4::p(double, double, double, double) --> void", pybind11::arg("xIn"), pybind11::arg("yIn"), pybind11::arg("zIn"), pybind11::arg("tIn"));
		cl.def("p", (void (Pythia8::Vec4::*)(class Pythia8::Vec4)) &Pythia8::Vec4::p, "C++: Pythia8::Vec4::p(class Pythia8::Vec4) --> void", pybind11::arg("pIn"));
		cl.def("px", (void (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::px, "C++: Pythia8::Vec4::px(double) --> void", pybind11::arg("xIn"));
		cl.def("py", (void (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::py, "C++: Pythia8::Vec4::py(double) --> void", pybind11::arg("yIn"));
		cl.def("pz", (void (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::pz, "C++: Pythia8::Vec4::pz(double) --> void", pybind11::arg("zIn"));
		cl.def("e", (void (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::e, "C++: Pythia8::Vec4::e(double) --> void", pybind11::arg("tIn"));
		cl.def("px", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::px, "C++: Pythia8::Vec4::px() const --> double");
		cl.def("py", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::py, "C++: Pythia8::Vec4::py() const --> double");
		cl.def("pz", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::pz, "C++: Pythia8::Vec4::pz() const --> double");
		cl.def("e", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::e, "C++: Pythia8::Vec4::e() const --> double");
		cl.def("__getitem__", (double & (Pythia8::Vec4::*)(int)) &Pythia8::Vec4::operator[], "C++: Pythia8::Vec4::operator[](int) --> double &", pybind11::return_value_policy::reference, pybind11::arg("i"));
		cl.def("mCalc", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::mCalc, "C++: Pythia8::Vec4::mCalc() const --> double");
		cl.def("m2Calc", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::m2Calc, "C++: Pythia8::Vec4::m2Calc() const --> double");
		cl.def("pT", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::pT, "C++: Pythia8::Vec4::pT() const --> double");
		cl.def("pT2", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::pT2, "C++: Pythia8::Vec4::pT2() const --> double");
		cl.def("pAbs", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::pAbs, "C++: Pythia8::Vec4::pAbs() const --> double");
		cl.def("pAbs2", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::pAbs2, "C++: Pythia8::Vec4::pAbs2() const --> double");
		cl.def("eT", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::eT, "C++: Pythia8::Vec4::eT() const --> double");
		cl.def("eT2", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::eT2, "C++: Pythia8::Vec4::eT2() const --> double");
		cl.def("theta", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::theta, "C++: Pythia8::Vec4::theta() const --> double");
		cl.def("phi", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::phi, "C++: Pythia8::Vec4::phi() const --> double");
		cl.def("thetaXZ", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::thetaXZ, "C++: Pythia8::Vec4::thetaXZ() const --> double");
		cl.def("pPos", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::pPos, "C++: Pythia8::Vec4::pPos() const --> double");
		cl.def("pNeg", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::pNeg, "C++: Pythia8::Vec4::pNeg() const --> double");
		cl.def("rap", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::rap, "C++: Pythia8::Vec4::rap() const --> double");
		cl.def("eta", (double (Pythia8::Vec4::*)() const) &Pythia8::Vec4::eta, "C++: Pythia8::Vec4::eta() const --> double");
		cl.def("rescale3", (void (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::rescale3, "C++: Pythia8::Vec4::rescale3(double) --> void", pybind11::arg("fac"));
		cl.def("rescale4", (void (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::rescale4, "C++: Pythia8::Vec4::rescale4(double) --> void", pybind11::arg("fac"));
		cl.def("flip3", (void (Pythia8::Vec4::*)()) &Pythia8::Vec4::flip3, "C++: Pythia8::Vec4::flip3() --> void");
		cl.def("flip4", (void (Pythia8::Vec4::*)()) &Pythia8::Vec4::flip4, "C++: Pythia8::Vec4::flip4() --> void");
		cl.def("rot", (void (Pythia8::Vec4::*)(double, double)) &Pythia8::Vec4::rot, "C++: Pythia8::Vec4::rot(double, double) --> void", pybind11::arg("thetaIn"), pybind11::arg("phiIn"));
		cl.def("rotaxis", (void (Pythia8::Vec4::*)(double, double, double, double)) &Pythia8::Vec4::rotaxis, "C++: Pythia8::Vec4::rotaxis(double, double, double, double) --> void", pybind11::arg("phiIn"), pybind11::arg("nx"), pybind11::arg("ny"), pybind11::arg("nz"));
		cl.def("rotaxis", (void (Pythia8::Vec4::*)(double, const class Pythia8::Vec4 &)) &Pythia8::Vec4::rotaxis, "C++: Pythia8::Vec4::rotaxis(double, const class Pythia8::Vec4 &) --> void", pybind11::arg("phiIn"), pybind11::arg("n"));
		cl.def("bst", (void (Pythia8::Vec4::*)(double, double, double)) &Pythia8::Vec4::bst, "C++: Pythia8::Vec4::bst(double, double, double) --> void", pybind11::arg("betaX"), pybind11::arg("betaY"), pybind11::arg("betaZ"));
		cl.def("bst", (void (Pythia8::Vec4::*)(double, double, double, double)) &Pythia8::Vec4::bst, "C++: Pythia8::Vec4::bst(double, double, double, double) --> void", pybind11::arg("betaX"), pybind11::arg("betaY"), pybind11::arg("betaZ"), pybind11::arg("gamma"));
		cl.def("bst", (void (Pythia8::Vec4::*)(const class Pythia8::Vec4 &)) &Pythia8::Vec4::bst, "C++: Pythia8::Vec4::bst(const class Pythia8::Vec4 &) --> void", pybind11::arg("pIn"));
		cl.def("bst", (void (Pythia8::Vec4::*)(const class Pythia8::Vec4 &, double)) &Pythia8::Vec4::bst, "C++: Pythia8::Vec4::bst(const class Pythia8::Vec4 &, double) --> void", pybind11::arg("pIn"), pybind11::arg("mIn"));
		cl.def("bstback", (void (Pythia8::Vec4::*)(const class Pythia8::Vec4 &)) &Pythia8::Vec4::bstback, "C++: Pythia8::Vec4::bstback(const class Pythia8::Vec4 &) --> void", pybind11::arg("pIn"));
		cl.def("bstback", (void (Pythia8::Vec4::*)(const class Pythia8::Vec4 &, double)) &Pythia8::Vec4::bstback, "C++: Pythia8::Vec4::bstback(const class Pythia8::Vec4 &, double) --> void", pybind11::arg("pIn"), pybind11::arg("mIn"));
		cl.def("rotbst", (void (Pythia8::Vec4::*)(const class Pythia8::RotBstMatrix &)) &Pythia8::Vec4::rotbst, "C++: Pythia8::Vec4::rotbst(const class Pythia8::RotBstMatrix &) --> void", pybind11::arg("M"));
		cl.def("eInFrame", (double (Pythia8::Vec4::*)(const class Pythia8::Vec4 &) const) &Pythia8::Vec4::eInFrame, "C++: Pythia8::Vec4::eInFrame(const class Pythia8::Vec4 &) const --> double", pybind11::arg("pIn"));
		cl.def("__sub__", (class Pythia8::Vec4 (Pythia8::Vec4::*)() const) &Pythia8::Vec4::operator-, "C++: Pythia8::Vec4::operator-() const --> class Pythia8::Vec4");
		cl.def("__iadd__", (class Pythia8::Vec4 & (Pythia8::Vec4::*)(const class Pythia8::Vec4 &)) &Pythia8::Vec4::operator+=, "C++: Pythia8::Vec4::operator+=(const class Pythia8::Vec4 &) --> class Pythia8::Vec4 &", pybind11::return_value_policy::reference, pybind11::arg("v"));
		cl.def("__isub__", (class Pythia8::Vec4 & (Pythia8::Vec4::*)(const class Pythia8::Vec4 &)) &Pythia8::Vec4::operator-=, "C++: Pythia8::Vec4::operator-=(const class Pythia8::Vec4 &) --> class Pythia8::Vec4 &", pybind11::return_value_policy::reference, pybind11::arg("v"));
		cl.def("__imul__", (class Pythia8::Vec4 & (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::operator*=, "C++: Pythia8::Vec4::operator*=(double) --> class Pythia8::Vec4 &", pybind11::return_value_policy::reference, pybind11::arg("f"));
		cl.def("__idiv__", (class Pythia8::Vec4 & (Pythia8::Vec4::*)(double)) &Pythia8::Vec4::operator/=, "C++: Pythia8::Vec4::operator/=(double) --> class Pythia8::Vec4 &", pybind11::return_value_policy::reference, pybind11::arg("f"));
		cl.def("__add__", (class Pythia8::Vec4 (Pythia8::Vec4::*)(const class Pythia8::Vec4 &) const) &Pythia8::Vec4::operator+, "C++: Pythia8::Vec4::operator+(const class Pythia8::Vec4 &) const --> class Pythia8::Vec4", pybind11::arg("v"));
		cl.def("__sub__", (class Pythia8::Vec4 (Pythia8::Vec4::*)(const class Pythia8::Vec4 &) const) &Pythia8::Vec4::operator-, "C++: Pythia8::Vec4::operator-(const class Pythia8::Vec4 &) const --> class Pythia8::Vec4", pybind11::arg("v"));
		cl.def("__mul__", (class Pythia8::Vec4 (Pythia8::Vec4::*)(double) const) &Pythia8::Vec4::operator*, "C++: Pythia8::Vec4::operator*(double) const --> class Pythia8::Vec4", pybind11::arg("f"));
		cl.def("__div__", (class Pythia8::Vec4 (Pythia8::Vec4::*)(double) const) &Pythia8::Vec4::operator/, "C++: Pythia8::Vec4::operator/(double) const --> class Pythia8::Vec4", pybind11::arg("f"));
		cl.def("__mul__", (double (Pythia8::Vec4::*)(const class Pythia8::Vec4 &) const) &Pythia8::Vec4::operator*, "C++: Pythia8::Vec4::operator*(const class Pythia8::Vec4 &) const --> double", pybind11::arg("v"));

		cl.def("__str__", [](Pythia8::Vec4 const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	// Pythia8::m(const class Pythia8::Vec4 &) file:Pythia8/Basics.h line:149
	M("Pythia8").def("m", (double (*)(const class Pythia8::Vec4 &)) &Pythia8::m, "C++: Pythia8::m(const class Pythia8::Vec4 &) --> double", pybind11::arg("v1"));

	// Pythia8::m(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) file:Pythia8/Basics.h line:150
	M("Pythia8").def("m", (double (*)(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &)) &Pythia8::m, "C++: Pythia8::m(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) --> double", pybind11::arg("v1"), pybind11::arg("v2"));

	// Pythia8::m2(const class Pythia8::Vec4 &) file:Pythia8/Basics.h line:151
	M("Pythia8").def("m2", (double (*)(const class Pythia8::Vec4 &)) &Pythia8::m2, "C++: Pythia8::m2(const class Pythia8::Vec4 &) --> double", pybind11::arg("v1"));

	// Pythia8::m2(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) file:Pythia8/Basics.h line:152
	M("Pythia8").def("m2", (double (*)(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &)) &Pythia8::m2, "C++: Pythia8::m2(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) --> double", pybind11::arg("v1"), pybind11::arg("v2"));

	// Pythia8::m2(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) file:Pythia8/Basics.h line:153
	M("Pythia8").def("m2", (double (*)(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &)) &Pythia8::m2, "C++: Pythia8::m2(const class Pythia8::Vec4 &, const class Pythia8::Vec4 &, const class Pythia8::Vec4 &) --> double", pybind11::arg("v1"), pybind11::arg("v2"), pybind11::arg("v3"));

}
