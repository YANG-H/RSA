#ifndef _RSA_HPP_
#define _RSA_HPP_

#include <vector>
#include <random>
#include <chrono>
#include <iostream>

#include "big_num.hpp"

namespace y
{

	namespace math
	{

		template <int len> class cipher;
		template <int len> class decipher;


		template <int len>
		class cipher
		{
		public:
			typedef std::vector<y::math::fixed_unsigned_int<len>> encoded_sequence_t;
			encoded_sequence_t operator ()(const std::string & msg) const
			{
				if(msg.empty())
					return encoded_sequence_t();

				std::vector<char> msgv(msg.begin(), msg.end());
				int msglen = msgv.size();
				encoded_sequence_t nseq((msglen * sizeof(char) / sizeof(encoded_sequence_t::value_type) + 5));

				int* intbuffer = ((int*)&(nseq[0]));
				intbuffer[0] = msglen;
				std::memcpy(intbuffer+1, &(msgv[0]), msglen * sizeof(char));

				// encode
				for(auto & i : nseq){
					i = i.mod_pow(e, n);
				}

				return nseq;
			}

		private:
			y::math::fixed_unsigned_int<len> n, e;
			template <int l> friend void generate_keys(cipher<l>& pbk, decipher<l>& pvk);
		};


		template <int len>
		class decipher
		{
		public:
			typedef std::vector<y::math::fixed_unsigned_int<len>> encoded_sequence_t;
			std::string operator () (const encoded_sequence_t& nseq_origin) const
			{
				if(nseq_origin.empty())
					return std::string();

				// decode
				encoded_sequence_t nseq(nseq_origin);
				for(auto & i : nseq){
					i = i.mod_pow(d, n);
				}

				const int* intbuffer = ((const int*)&(nseq[0]));
				int msglen = *intbuffer;

				std::vector<char> msgv(msglen, 'N');
				std::memcpy(&msgv[0], intbuffer+1, msglen * sizeof(char));
				return std::string(msgv.begin(), msgv.end());
			}

		private:
			y::math::fixed_unsigned_int<len> n, e, d, phi, p, q;
			template <int l> friend void generate_keys(cipher<l>& pbk, decipher<l>& pvk);
		};



		template <int len>
		void generate_keys(cipher<len>& pbk, decipher<len>& pvk)
		{
			typedef y::math::fixed_unsigned_int<len/2> half_num_t;
			typedef y::math::fixed_unsigned_int<len> num_t;

			std::default_random_engine rnd;
			rnd.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

			half_num_t p, q;
			num_t n, phi, e, d;

			while(true){
				p.set_to_probable_prime(rnd);
				q.set_to_probable_prime(rnd);

				std::cout << "p: \t" << p << std::endl;
				std::cout << "q: \t" << q << std::endl;

				half_num_t p_1 = p, q_1 = q;
				-- p_1; -- q_1;

				n = y::math::mult_unlimit(p, q);
				phi = y::math::mult_unlimit(p_1, q_1);

				std::cout << "n: \t" << n << std::endl;
				std::cout << "phi: \t" << phi << std::endl;

				std::cout << "= locating e and d ... =" << std::endl;
				for(int i = 0; i < 100; i++){
					e.set_to_probable_prime(rnd);
					if(!e.is_coprime_with(phi))
						continue;
					
					std::cout << "found!" << std::endl;

					d = e.mod_inverse(phi);
					std::cout << "e: \t" << e << std::endl;
					std::cout << "d: \t" << d << std::endl;
					std::cout << "ed mod phi: \t" << e.mod_mult(d, phi) << std::endl;

					pbk.n = n; pbk.e = e;

					pvk.n = n; pvk.e = e; pvk.d = d; 
					pvk.phi = phi; 
					pvk.p = num_t(p); 
					pvk.q = num_t(q);

					
					return;
				}
			}
		}

	}

}


#endif