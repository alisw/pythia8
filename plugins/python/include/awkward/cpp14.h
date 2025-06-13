// cpp11.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Helper file to provide C++14 features needed by Awkward, when
// building with C++11.

#ifndef awkward_cpp14_H
#define awkward_cpp14_H

// Only provide these features if not natively available.
#if __cplusplus < 201402L

//==========================================================================

// Internal header code from GCC 5.5.0 used to provide C++14
// functionality. These are wrapped in the internal cpp14 namespace.

namespace cpp14 {

//--------------------------------------------------------------------------

// Copied from line 203 of include/c++/5/utility.

// Stores a tuple of indices.  Used by tuple and pair, and by bind() to
// extract the elements in a tuple.
template<size_t... _Indexes> struct _Index_tuple { };

// Concatenates two _Index_tuples.
template<typename _Itup1, typename _Itup2> struct _Itup_cat;

template<size_t... _Ind1, size_t... _Ind2>
  struct _Itup_cat<_Index_tuple<_Ind1...>, _Index_tuple<_Ind2...>>
  {
    using __type = _Index_tuple<_Ind1..., (_Ind2 + sizeof...(_Ind1))...>;
  };

// Builds an _Index_tuple<0, 1, 2, ..., _Num-1>.
template<size_t _Num>
  struct _Build_index_tuple
  : _Itup_cat<typename _Build_index_tuple<_Num / 2>::__type,
      	typename _Build_index_tuple<_Num - _Num / 2>::__type>
  { };

template<>
  struct _Build_index_tuple<1>
  {
    typedef _Index_tuple<0> __type;
  };

template<>
  struct _Build_index_tuple<0>
  {
    typedef _Index_tuple<> __type;
  };

//--------------------------------------------------------------------------

// Copied from line 239 of include/c++/5/utility.

/// Class template integer_sequence
template<typename _Tp, _Tp... _Idx>
  struct integer_sequence
  {
    typedef _Tp value_type;
    static constexpr size_t size() { return sizeof...(_Idx); }
  };

template<typename _Tp, _Tp _Num,
         typename _ISeq = typename _Build_index_tuple<_Num>::__type>
  struct _Make_integer_sequence;

template<typename _Tp, _Tp _Num,  size_t... _Idx>
  struct _Make_integer_sequence<_Tp, _Num, _Index_tuple<_Idx...>>
  {
    static_assert( _Num >= 0,
      	     "Cannot make integer sequence of negative length" );

    typedef integer_sequence<_Tp, static_cast<_Tp>(_Idx)...> __type;
  };

}

//==========================================================================

// Expose the the relevant methods in the std namespace.

namespace std {

//--------------------------------------------------------------------------

// Copied from line 770 of include/c++/5/tuple.

// Duplicate of C++14's tuple_element_t for internal use in C++11 mode
template<std::size_t __i, typename _Tp>
  using tuple_element_t = typename tuple_element<__i, _Tp>::type;

//--------------------------------------------------------------------------

// Copied from line 264 of include/c++/5/utility.

/// Alias template make_integer_sequence
template<typename _Tp, _Tp _Num>
  using make_integer_sequence
    = typename cpp14::_Make_integer_sequence<_Tp, _Num>::__type;

/// Alias template index_sequence
template<size_t... _Idx>
  using index_sequence = cpp14::integer_sequence<size_t, _Idx...>;

/// Alias template make_index_sequence
template<size_t _Num>
  using make_index_sequence = make_integer_sequence<size_t, _Num>;

/// Alias template index_sequence_for
template<typename... _Types>
  using index_sequence_for = make_index_sequence<sizeof...(_Types)>;

//--------------------------------------------------------------------------

// Constant iterators.

template<typename T> typename T::const_iterator cbegin(T data)
  noexcept { return begin(data); }
template<typename T> typename T::const_iterator cend(T data)
  noexcept { return end(data); }

}

#endif // __cplusplus < 201402L

#endif // awkward_cpp14_H
