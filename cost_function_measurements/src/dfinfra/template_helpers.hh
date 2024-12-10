/*
 * Helper constructs for templates.
 */

/*
 * Templated "type chooser" based on boolean flag
 *
 * XXX obsolete with C++11's std::conditional<bool B, class T, class F> XXX
 *
 * Can be used to avoid code duplication when implementing const iterators
 * and non-const iterators like so:
 * Inside the iterator, we usually have typedefs for reference and pointer types:
 *   non-const               |  const
 *   using reference_t = T&  |  using reference_t = const T&
 *   using pointer_t   = T*  |  using pointer_t   = const T*
 * 
 * With the template helper, this can be "simplified" like this:
 * non-const
 *   using reference_t = typename choose_type<false, T&, const T&>::type
 *   using pointer_t   = typename choose_type<false, T*, const T*>::type
 * const
 *   using reference_t = typename choose_type<true, T&, const T&>::type
 *   using pointer_t   = typename choose_type<true, T*, const T*>::type
 */
template <bool flag, typename TypeIfTrue, typename TypeIfFalse>
struct choose_type;

template <typename TypeIfTrue, typename TypeIfFalse>
struct choose_type<true, TypeIfTrue, TypeIfFalse> {
  using type = TypeIfTrue;
};

template <typename TypeIfTrue, typename TypeIfFalse>
struct choose_type<false, TypeIfTrue, TypeIfFalse> {
  using type = TypeIfFalse;
};
