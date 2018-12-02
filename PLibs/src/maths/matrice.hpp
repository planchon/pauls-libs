#ifndef PL_MATRIX_H
#define PL_MATRIX_H

namespace plibs {
	template<typename T, int cols, int rows>
	struct Matrix {
		T data[cols][rows];

		Matrix () {}

		Matrix (std::vector<std::vector<T>> base) {
			if (base.size() != cols || base.at(0).size() != rows) {
				throw std::invalid_argument("PLIBS ERROR : You can't initalize Matrice with vector of different dimension");
			}

			for(int i = 0; i < cols; i++)
				for(int j = 0; j < rows; j++) 
					data[i][j] = base.at(i).at(j);
		}

		Matrix (T liste[cols][rows]) {
			for(int i = 0; i < cols; i++)
				for(int j = 0; j < rows; j++)
					data[i][j] = liste[i][j];
		}

		friend std::ostream& operator << (std::ostream& stream, const Matrix& mat) {
			// on cherche le nombre le plus grand (en taille) nombre
			int max = -1000000;
			for (int i = 0; i < cols * rows; i++) {
				if (abs(mat.data[i / cols][i % rows]) > max) {
					max = mat.data[i % cols][i / rows];
				}
			}

			int maxIntSize = 0;
			(max < 0) ? maxIntSize = floor(log10(-max)) + 1: maxIntSize = floor(log10(max));
			
			stream << "[";
			for (int i = 0; i < cols; i++) {
				(i == 0) ? stream << "[" : stream << " [";
				for (int j = 0; j < rows; j++) {
					stream << mat.data[i][j];
					for (int k = floor(log10(mat.data[i][j])); k < maxIntSize; k++)
						stream << " ";
					(j == rows - 1) ? stream << "" : stream << " ";
				}
				(i == cols - 1) ? stream << "]" : stream << "]\n";
			}
			stream << "]";
		}

		friend Matrix operator + (Matrix& mat1, Matrix& mat2) {
			Matrix mat;
			for (int i = 0; i < cols; i++)
				for (int j = 0; j < rows; j++)
					mat.data[i][j] = mat1.data[i][j] + mat2.data[i][j];
			return mat;
		}

		friend Matrix operator += (Matrix& mat1, Matrix& mat2) {
			return mat1 + mat2;
		}

		friend Matrix operator - (Matrix& mat1, Matrix& mat2) {
			Matrix mat;
			for (int i = 0; i < cols; i++)
				for (int j = 0; j < rows; j++)
					mat.data[i][j] = mat1.data[i][j] - mat2.data[i][j];
			return mat;
		}

		friend Matrix operator -= (Matrix& mat1, Matrix& mat2) {
			return mat1 - mat2;
		}

		friend Matrix operator * (Matrix& mat1, T lambda) {
			Matrix mat;
			for (int i = 0; i < cols; i++)
				for (int j = 0; j < rows; j++)
					mat.data[i][j] = mat1.data[i][j] * lambda;
			return mat;
		}

		friend Matrix operator *= (Matrix& mat1, T lambda) {
			return mat1 * lambda;
		}

		friend Matrix operator / (Matrix& mat1, T lambda) {
			if (lambda == 0) {
				throw std::invalid_argument("PLIBS ERROR : Division by 0");
			}
			Matrix<T, cols, rows> mat;
			for (int i = 0; i < cols; i++)
				for (int j = 0; j < rows; j++)
					mat.data[i][j] = mat1.data[i][j] * lambda;
			return mat;
		}

		friend Matrix operator /= (Matrix& mat1, T lambda) {
			return mat1 / lambda;
		}

		template<int cols2, int rows2>
		inline Matrix<T, cols, rows2> operator * (Matrix<T, cols2, rows2>& mat2) {
			if (rows != cols2) {
				throw std::invalid_argument("PLIBS ERROR : Matrix multiplication dimension error.");
			}

			const int Al = cols;
			const int Ac = rows;
			const int Bl = cols2;
			const int Bc = rows2;
			
			double*  A = new double[Al*Ac];
			double*  B = new double[Bl*Bc];
			double*  C = new double[Al*Bc];

			// A(m,n) m cols, n rows
    
			//Fill A and B with random numbers
			for(uint i = 0; i < Al; i++){
				for(uint j = 0; j< Ac; j++){
					A[i * Ac + j] = (double) data[i][j];
				}
			}

			for(uint i = 0; i < Bl; i++){
				for(uint j= 0; j< Bc; j++){
					B[i * Bc + j] = (double) mat2.data[i][j];
				}
			}
    
			// //Calculate A*B=C
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Al, Bc, Ac, 1.0, A, Ac, B, Bc, 0.0, C, Al);
			Matrix<T, Al, Bc> matFinal;

			for(uint i = 0; i < Al; i++){
				for(uint j = 0; j< Bc; j++){
					matFinal.data[i][j] = (T) A[i * Bc + j];
				}
			}

			return matFinal;
		}
	};
}	
#endif
