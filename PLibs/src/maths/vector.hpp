#ifndef PL_VECTOR_H
#define PL_VECTOR_H

namespace plibs {	
	template<typename T, int _Cols>
	struct Vector {
	public:
		T data[_Cols] = {0};

		Vector () {}
		
		Vector (std::vector<T> base) {
			if (_Cols != base.size()) {
				throw std::invalid_argument("PLIBS ERROR : You can't initalize Vector with vector of different size");
			}
			for(int i = 0; i < base.size(); i++)
				data[i] = base.at(i);
		}
		
		Vector (T* liste) {
			int size = (sizeof(liste) / sizeof(*liste)) + 1;
			
			if (_Cols != size) {
				std::cout << _Cols << ", " << size << std::endl;
				throw std::invalid_argument("PLIBS ERROR : You can't initalize Vector with list of different size");
			}
			for(int i = 0; i < size; i++)
				data[i] = liste[i];
		}
		
		friend std::ostream& operator << (std::ostream& stream, const Vector<T, _Cols>& vec) {
			std::cout << "[";
			for (int i = 0; i < _Cols - 1; i++) {
				std::cout << vec.data[i] << " ";
			}
			std::cout << vec.data[_Cols - 1] << "]";			
		}

		friend Vector<T, _Cols> operator + (const Vector<T, _Cols>& vec1, const Vector<T, _Cols>& vec2) {
			Vector<T, _Cols> finalVec;
			for (int i = 0; i < _Cols; i++) {
				finalVec.data[i] = vec1.data[i] + vec2.data[i];
			}

			return finalVec;
		}

		friend Vector<T, _Cols> operator += (const Vector<T, _Cols>& vec1, const Vector<T, _Cols>& vec2) {
			return vec1 + vec2;
		}
		
		friend Vector<T, _Cols> operator - (const Vector<T, _Cols>& vec1, const Vector<T, _Cols>& vec2) {
			Vector<T, _Cols> finalVec;
			for (int i = 0; i < _Cols; i++) {
				finalVec.data[i] = vec1.data[i] - vec2.data[i];
			}

			return finalVec;
		}

		friend Vector<T, _Cols> operator -= (const Vector<T, _Cols>& vec1, const Vector<T, _Cols>& vec2) {
			return vec1 - vec2;
		}
		

		friend Vector<T, _Cols> operator * (const Vector<T, _Cols>& vec1, const T lambda) {
			Vector<T, _Cols> finalVec;
			for (int i = 0; i < _Cols; i++) {
				finalVec.data[i] = vec1.data[i] * lambda;
			}

			return finalVec;
		}

		friend Vector<T, _Cols> operator *= (const Vector<T, _Cols>& vec1, const T lambda) {
			return vec1 * lambda;
		}

		friend Vector<T, _Cols> operator / (const Vector<T, _Cols>& vec1, const T lambda) {
			if (lambda != 0) {
				Vector<T, _Cols> finalVec;
				for (int i = 0; i < _Cols; i++) {
					finalVec.data[i] = vec1.data[i] / lambda;
				}

				return finalVec;				
			} else {
				throw std::invalid_argument("PLIBS ERROR : Division by 0");				
			}
		}

		friend Vector<T, _Cols> operator /= (const Vector<T, _Cols>& vec1, const T lambda) {
			return vec1 / lambda;
		}

		
		// dot product here
		friend T operator * (const Vector<T, _Cols>& vec1, const Vector<T, _Cols>& vec2) {
			T dot = 0;
			for (int i = 0; i < _Cols; i++) {
				dot += vec1.data[i] * vec2.data[i];
			}

			return dot;
		}

		// cross product
		friend Vector<T, 3> operator ^ (const Vector<T, _Cols>& vec1, const Vector<T, _Cols>& vec2) {
			if (_Cols != 3) {
				throw std::invalid_argument("PLIBS ERROR : Cross product can only be applied between 3D vector");
			}

			Vector<T, 3> crossed;
			crossed.data[0] = vec1.data[1] * vec2.data[2] - vec1.data[2] * vec2.data[1];
			crossed.data[1] = vec1.data[2] * vec2.data[0] - vec1.data[0] * vec2.data[2];
			crossed.data[2] = vec1.data[0] * vec2.data[1] - vec1.data[1] * vec2.data[0];
			
			return crossed;
		}

		// normalisation
		inline Vector<float, _Cols> normalize () {
			Vector<float, _Cols> normalized;
			float length = sqrt(norme());
			if (length == 0) {throw std::invalid_argument("PLIBS ERROR : Division by 0");}  
			for(int i = 0; i < _Cols; i++)
				normalized.data[i] = (float) data[i] / length;

			return normalized;
		}
		
		inline float norme () {
			float sum = 0;
			for(int i = 0; i < _Cols; i++)
				sum += data[i] * data[i];
			return sum;
		}

		inline float length () {
			return sqrt(norme());
		}
	};
}

#endif
