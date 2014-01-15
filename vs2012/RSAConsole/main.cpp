
#include <iostream>
#include <sstream>
//#include "rsa.h"
#include "y/math/big_num.hpp"

using namespace std;
using namespace y::math;

int	main (void)
{	
	unsigned d1 [] = {0U, 0xffffffffU, 0U, 0U, 1U};
	std::reverse(d1, d1+5);
	unsigned d2 [] = {0U, 0xffffffffU};
	std::reverse(d2, d2+2);
	__private__::subtract_range(d1, d1+5, d2, d2+2);
	cout << fixed_unsigned_int<5>(d1, d1+5) << endl;

	while(true){
		fixed_unsigned_int<3> a, b;
		cin >> a >> b;

		cout << "inputs a: " << a << " b: " << b << endl;
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
		cout << "======================================" << endl;

		/*for(int i = 1; i <= 10; i++){
		cout << "a << " << i*4 << ": " << (a.left_shifted(i*4)) << endl;			
		}
		for(int i = 1; i <= 10; i++){
		cout << "a >> " << i*4 << ": " << (a.right_shifted(i*4)) << endl;
		}*/

		cout << "most significant non-zero unit of a: " << a.most_significant_nonzero_unit() << endl;
		cout << "most significant non-zero bit of a: " << a.most_significant_nonzero_bit() << endl;
		

		/*
		big_int a = 1;
		big_int b = std::numeric_limits<unsigned long long>::max();
		big_int aa = "1";
		char bf[1024];
		big_int bb(itoa(std::numeric_limits<int>::max(), bf, 16), 16);

		cout << a.compare_to(aa) << endl;

		big_integer p = big_integer::probable_prime(32);
		big_integer q = big_integer::probable_prime(32);

		big_integer n = p * q;
		big_integer p_1 = p - 1;
		big_integer q_1 = q - 1;
		big_integer p_1_q_1 = p_1 * q_1;

		int initpubkey = 65535;
		while(true) {
		big_integer gcd = p_1_q_1.gcd(initpubkey);
		if(gcd == 1)
		break;
		initpubkey ++;
		}

		big_integer pubkey = initpubkey;
		big_integer prvkey = pubkey.mod_inverse(p_1_q_1);

		std::string msg = "23123123231sdaddfsfdsfs";
		big_integer msgval(msg, 32);

		big_integer cipherval = msgval.mod_power(pubkey, n);
		cout << cipherval << endl;

		big_integer plainval = cipherval.mod_power(prvkey, n);
		cout << plainval.to_string(32) << endl;*/

		std::system("pause");
	}

	return 0;
}
