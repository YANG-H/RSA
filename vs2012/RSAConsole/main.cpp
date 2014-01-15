
#include <iostream>
#include <sstream>

#include "y/math/big_num.hpp"
#include "y/math/rsa.hpp"


int	main (int argc, char ** argv)
{
	static const int SIZE = 4;

	y::math::cipher<SIZE> cipher;
	y::math::decipher<SIZE> decpher;
	y::math::generate_keys<SIZE>(cipher, decpher);

	while(true){
		std::string msg;
		std::cout << "input your message: ";
		std::cin >> msg;
		std::cout << msg << std::endl;

		auto encoded = cipher(msg);
		
		std::string decoded = decpher(encoded);

		std::cout << decoded << std::endl;
	}

	return 0;
}
