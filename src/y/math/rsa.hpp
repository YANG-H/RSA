#ifndef _RSA_HPP_
#define _RSA_HPP_

#include <vector>
#include "big_num.hpp"

namespace y
{

	namespace math
	{

		struct message
		{
			fixed_unsigned_int<15> msg;
			fixed_unsigned_int<15> reciever_pub_key;
		};

	}

}


#endif