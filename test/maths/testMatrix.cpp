#include <PLibs/Maths>

int main() {
	plibs::Matrix<int, 2, 3> matrice ({{1,2,3},{1,2,3}});
	plibs::Matrix<int, 2, 3> matrice3({{1,2,3},{1,2,3}});
	plibs::Matrix<int, 3, 2> matrice2({{1,2},{3,4},{5,6}});
	std::cout << matrice  << std::endl;
	std::cout << matrice2 << std::endl;
	
	std::cout << matrice * matrice2 << std::endl;
	return 0;
}
