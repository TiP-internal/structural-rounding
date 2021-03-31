
#ifndef UTIL_H
#define UTIL_H

#include <tuple>

namespace sr_apx {
namespace util {

class GrayCode {
private:
	int index = 0;
	int digits = 4;

public:
	GrayCode() = default;
	GrayCode(int);
	GrayCode(int, int);

	GrayCode& operator++();
	GrayCode operator++(int);

	int operator*();

	bool operator==(const GrayCode&);
	bool operator!=(const GrayCode&);

	int get_code();
	int get_change();
	int get_direction();
};

int log2(int);
int pow(int, int);
std::tuple<GrayCode, GrayCode> graycode(int);

}}

#endif
