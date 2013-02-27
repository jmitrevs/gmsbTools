#ifndef GMSBTOOLS_SORTHELPERS_H  
#define GMSBTOOLS_SORTHELPERS_H 

#include <list>
#include <utility>
namespace SortHelpers {

  typedef std::list<std::pair<std::size_t, float> > sl_t;

  template <typename T>
  void sort(sl_t& sl, const T& o) {
    sl.clear();
    for (std::size_t idx = 0; idx < static_cast<std::size_t>(o.n()); idx++) {
      sl_t::iterator it = sl.begin();
      for (; it != sl.end(); ++it) {
	if (o.pt(idx) > it->second) break;
      }
      sl.insert(it, std::pair<std::size_t, float>(idx, o.pt(idx)));
    }
  }
}

#endif
