
#include "sr_apx/util/util.hpp"
#include <iostream>

namespace sr_apx {
namespace util {

int log2(int n) {
	int k = 0;
	while ((1 << k) < n) {
		++k;
	}
	return k;
}

int pow(int b, int e) {
	int x = 1;
	for (int i = 0; i < e; i++) {
		x *= b;
	}

	return x;
}

GrayCode::GrayCode(int d): digits(d) {}

GrayCode::GrayCode(int d, int i): digits(d), index(i) {}

GrayCode& GrayCode::operator++() {
	++index;
	return *this;
}

GrayCode GrayCode::operator++(int) {
	GrayCode temp = *this;
	++*this;
	return temp;
}

int GrayCode::operator*() {
	return index ^ (index >> 1);
}

bool GrayCode::operator==(const GrayCode& gc) {
	return index == gc.index && digits == gc.digits;
}

bool GrayCode::operator!=(const GrayCode& gc) {
	return !(*this == gc);
}

int GrayCode::get_code() {
	return index ^ (index >> 1);
}

int GrayCode::get_change() {
	if (index == (1 << digits)) {
		return digits - 1;
	}

	int curr = index ^ (index >> 1);
	int prev = (index - 1) ^ ((index - 1) >> 1);

	return log2(curr ^ prev);
}

int GrayCode::get_direction() {
	if (index == (1 << digits)) {
		return 0;
	}

	int change = get_change();
	return (get_code() & (1 << change)) >> change;
}

std::tuple<GrayCode, GrayCode> graycode(int digits) {
	return std::make_tuple(GrayCode(digits), GrayCode(digits, 1 << digits));
}

}}
