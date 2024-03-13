#include <cwchar>
#include <ios>
#include <sstream> // __str__

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

void bind_std_postypes(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // std::fpos file:bits/postypes.h line:112
		pybind11::class_<std::fpos<__mbstate_t>, std::shared_ptr<std::fpos<__mbstate_t>>> cl(M("std"), "fpos___mbstate_t_t", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new std::fpos<__mbstate_t>(); } ) );
		cl.def( pybind11::init<std::streamoff>(), pybind11::arg("__off") );

		cl.def( pybind11::init( [](std::fpos<__mbstate_t> const &o){ return new std::fpos<__mbstate_t>(o); } ) );
		cl.def("__iadd__", (class std::fpos<__mbstate_t> & (std::fpos<__mbstate_t>::*)(std::streamoff)) &std::fpos<__mbstate_t>::operator+=, "C++: std::fpos<__mbstate_t>::operator+=(std::streamoff) --> class std::fpos<__mbstate_t> &", pybind11::return_value_policy::reference, pybind11::arg("__off"));
		cl.def("__isub__", (class std::fpos<__mbstate_t> & (std::fpos<__mbstate_t>::*)(std::streamoff)) &std::fpos<__mbstate_t>::operator-=, "C++: std::fpos<__mbstate_t>::operator-=(std::streamoff) --> class std::fpos<__mbstate_t> &", pybind11::return_value_policy::reference, pybind11::arg("__off"));
		cl.def("__add__", (class std::fpos<__mbstate_t> (std::fpos<__mbstate_t>::*)(std::streamoff) const) &std::fpos<__mbstate_t>::operator+, "C++: std::fpos<__mbstate_t>::operator+(std::streamoff) const --> class std::fpos<__mbstate_t>", pybind11::arg("__off"));
		cl.def("__sub__", (class std::fpos<__mbstate_t> (std::fpos<__mbstate_t>::*)(std::streamoff) const) &std::fpos<__mbstate_t>::operator-, "C++: std::fpos<__mbstate_t>::operator-(std::streamoff) const --> class std::fpos<__mbstate_t>", pybind11::arg("__off"));
		cl.def("__sub__", (std::streamoff (std::fpos<__mbstate_t>::*)(const class std::fpos<__mbstate_t> &) const) &std::fpos<__mbstate_t>::operator-, "C++: std::fpos<__mbstate_t>::operator-(const class std::fpos<__mbstate_t> &) const --> std::streamoff", pybind11::arg("__other"));
	}
}
