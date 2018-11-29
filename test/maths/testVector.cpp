#include <PLibs/Core>

int main() {
	plibs::Vector<int, 3> vec({1,2,3});
	plibs::Vector<int, 3> vec2({2,-1,3});

	std::cout << (vec ^ vec2) << std::endl;
	std::cout << vec.normalize() << std::endl;
	return 0;
}
