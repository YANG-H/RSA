
#include <iostream>
#include <sstream>
#include <random>
#include <chrono>
#include "y/math/big_num.hpp"

using namespace std;
using namespace y::math;

int	main (void)
{

	std::default_random_engine generator;
	generator.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	
	/*
	while(true){
		fixed_unsigned_int<10> a, b;
		a = gen_prime_with_certainty<10>(generator);
		b = gen_prime_with_certainty<10>(generator);
		cout << "======================================" << endl;
		cout << "a: " << a << endl;
		cout << "b: " << b << endl;
		cout << "gcd(a, b): " << a.gcd(b) << endl;
		cout << "gcd(b, a): " << b.gcd(a) << endl;
		cout << "gcd(a, b): " << a.binary_gcd(b) << endl;
		cout << "gcd(b, a): " << b.binary_gcd(a) << endl;
		cout << "a^-1 mod b: " << a.mod_inverse(b) << endl;
		cout << "b^-1 mod a: " << b.mod_inverse(a) << endl;
		cout << "gcd(a^-1 mod b, b): " << a.mod_inverse(b).gcd(b) << endl;
		cout << "gcd(b^-1 mod a, a): " << b.mod_inverse(a).gcd(a) << endl;
		cout << "gcd(a^-1 mod b, b): " << a.mod_inverse(b).binary_gcd(b) << endl;
		cout << "gcd(b^-1 mod a, a): " << b.mod_inverse(a).binary_gcd(a) << endl;

		cout << "======================================" << endl;
		for(int i = 1; i <= 10; i++){
			cout << "a^-" << i << " mod b: " << a.mod_inverse(b).mod_pow(i, b) << endl;
			cout << "a^+" << i << " mod b: " << a.mod_pow(i, b) << endl;
			cout << "a^-" << i << " * a^" << i << " mod b: " << a.mod_inverse(b).mod_pow(i, b).mod_mult(a.mod_pow(i, b), b) << endl;
			cout << "-------------------------------------" << endl;
			cout << "b^-" << i << " mod a: " << b.mod_inverse(a).mod_pow(i, a) << endl;
			cout << "b^+" << i << " mod a: " << b.mod_pow(i, a) << endl;
			cout << "b^-" << i << " * b^" << i << " mod a: " << b.mod_inverse(a).mod_pow(i, a).mod_mult(b.mod_pow(i, a), a) << endl;
		}

		std::system("pause");
	}

	while(false){
		fixed_unsigned_int<3> aa("111111112222222233333333", 16);
		fixed_unsigned_int<1> bb("aaaaaaaa", 16);
		cout << hex << "aa: " << aa << endl;
		cout << hex << "bb: " << bb << endl;
		fixed_unsigned_int<5> tmp;
		__private__::divide_and_remain_range(aa.begin(), aa.end(), bb.begin(), bb.end(), tmp.begin(), tmp.end());
		assert(aa < fixed_unsigned_int<3>(bb));
		system("pause");
	}

	while(true){
		fixed_unsigned_int<3> a, b;
		//cin >> a >> b;
		a.randomize(generator);
		b.randomize(generator);

		cout << "a: " << a << endl;
		cout << "b: " << b << endl;
		cout << "======================================" << endl;
		cout << "a + b: " << (a+b) << endl;
		cout << "-a: " << (-a) << endl;
		cout << "-b: " << (-b) << endl;
		cout << "a - b: " << (a-b) << endl;
		cout << "b - a: " << (b-a) << endl;
		cout << "++a: " << (++a) << endl;
		cout << "--a: " << (--a) << endl;
		cout << "++b: " << (++b) << endl;
		cout << "--b: " << (--b) << endl;
		cout << "a * b: " << (a*b) << endl;
		cout << "a * b [unwrapped]: " << (mult_unlimit(a, b)) << endl;
		cout << "a / b: " << (a/b) << endl;
		cout << "a mod b: " << (a.mod(b)) << endl;
		cout << "b / a: " << (b/a) << endl;
		cout << "b mod a: " << (b.mod(a)) << endl;
		cout << "a / b * b + (a mod b) - a: " << (a/b*b + a.mod(b) - a) << endl;
		cout << "b / a * a + (b mod a) - b: " << (b/a*a + b.mod(a) - b) << endl;
		cout << "gcd(a, b): " << a.gcd(b) << endl;
		cout << "gcd(b, a): " << b.gcd(a) << endl;
		cout << "gcd(gcd(a, b), a): " << a.gcd(b).gcd(a) << endl;
		cout << "gcd(gcd(b, a), b): " << a.gcd(b).gcd(a) << endl;
		cout << "a^-1 mod b: " << a.mod_inverse(b) << endl;
		cout << "b^-1 mod a: " << b.mod_inverse(a) << endl;
		cout << "gcd(a^-1 mod b, b): " << a.mod_inverse(b).gcd(b) << endl;
		cout << "gcd(b^-1 mod a, a): " << b.mod_inverse(a).gcd(a) << endl;
		cout << "======================================" << endl;

		fixed_unsigned_int<3> temp = 1;
		for(int i = 1; i <= 5; i++){
			cout << " a^" << i << " mod b " << ": " << a.mod_pow(i, b) << endl;
			temp = temp.mod_mult(a, b);
			cout << "[a^" << i << " mod b]" << ": " << temp << endl;
		}
		temp = 1;
		for(int i = 1; i <= 5; i++){
			cout << "b^" << i << " mod a" << ": " << b.mod_pow(i, a) << endl;
			temp = temp.mod_mult(b, a);
			cout << "[b^" << i << " mod a]" << ": " << temp << endl;
		}

		cout << "======================================" << endl;
		cout << "most significant non-zero unit of a: " << a.most_significant_nonzero_unit() << endl;
		cout << "most significant non-zero bit of a: " << a.most_significant_nonzero_bit() << endl;
		cout << "least significant non-zero unit of a: " << a.least_significant_nonzero_unit() << endl;
		cout << "least significant non-zero bit of a: " << a.least_significant_nonzero_bit() << endl;
		cout << "most significant non-zero unit of b: " << b.most_significant_nonzero_unit() << endl;
		cout << "most significant non-zero bit of b: " << b.most_significant_nonzero_bit() << endl;
		cout << "least significant non-zero unit of b: " << b.least_significant_nonzero_unit() << endl;
		cout << "least significant non-zero bit of b: " << b.least_significant_nonzero_bit() << endl;
		cout << "number of bits of a: " << a.number_of_bits() << endl;
		cout << "number of bits of b: " << b.number_of_bits() << endl;
		cout << "======================================" << endl;

		

		std::system("pause");
	}
	*/
	return 0;
}
