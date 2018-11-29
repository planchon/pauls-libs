#ifndef PL_MISC
#define PL_MISC

namespace plibs {
	template<typename T, int N>
	int size(T (&arr)[N]){
		return sizeof(arr) / sizeof(arr[0]);
	}	
}

#endif
