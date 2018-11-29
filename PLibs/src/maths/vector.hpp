#ifndef PL_VECTOR_H
#define PL_VECTOR_H

namespace plibs {	
	template<typename T, int _Cols>
		struct Vector {
			int data[_Cols];
		};
	
	template<typename T, int _Cols>
	inline void printVector(Vector<T, _Cols>& vec) {
		std::cout << _Cols << std::endl;
	}
}

#endif
