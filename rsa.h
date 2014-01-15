#ifndef RSA_H
#define RSA_H

#include <random>

namespace y
{

	namespace util
	{

		class big_int;
		class bit_sieve;

		class bit_sieve
		{
			std::vector<long long> bits_;
			int length_;

		public:
			typedef std::default_random_engine random_engine;

			bit_sieve(const big_int & base, int searchlen);
			// TODO // BitSieve

		public:
			big_int retrieve(const big_int & initvalue, int certainty, random_engine & random);
		};

		class big_int
		{
			int signum_;
			std::vector<int> mag_;

		public:
			typedef std::default_random_engine random_engine;

			big_int();
			big_int(int value);
			big_int(unsigned int value);
			big_int(long long value);
			big_int(unsigned long long value);

			explicit big_int(const std::vector<signed char>& val);
			explicit big_int(const std::vector<int>& val);
			explicit big_int(int signum, const std::vector<signed char>& magnitude);
			explicit big_int(int signum, const std::vector<int>& magnitude);

			big_int(const std::string & s, int radix = 10);

			explicit big_int(int numbits, random_engine & rand);

			big_int(const big_int & i);
			big_int & operator = (const big_int & i);

			static const big_int ZERO; // TODO
			static const big_int ONE; // TODO

		protected:
			explicit big_int(const std::vector<int>& magnitude, int signum);

		public:	
			int compare_magnitude(const big_int & b) const;
			int compare_to(const big_int & b) const;
			bool operator == (const big_int & i) const {return compare_to(i) == 0;}
			bool operator != (const big_int & i) const {return compare_to(i) != 0;}
			bool operator < (const big_int & i) const {return compare_to(i) == -1;}
			bool operator > (const big_int & i) const {return compare_to(i) == 1;}
			bool is_zero() const;
			bool is_null() const;

			std::vector<signed char> to_byte_array() const;
			std::string to_string(int radix = 10) const;
			long long to_longlong() const;
			int bit_length() const;
			big_int set_bit(int n) const;

		public:
			void negate();
			big_int negated() const;
			big_int abs() const;
			
			int divide_one_word(int divisor, big_int & quotient) const;
			big_int divide(const big_int & b, big_int & quotient) const;


			big_int remainder(const big_int & i) const;
			big_int gcd(const big_int & i) const;
			big_int mod_inverse(const big_int & n) const;
			big_int mod_power(const big_int & i, const big_int & n) const;
			bool prime_to_certainty(int certainty, random_engine & rand) const;

		public:
			static big_int probable_prime(int bitlen, random_engine & rand);

		private:
			static big_int small_prime(int bitlen, int certainty, random_engine & rand);
			static big_int large_prime(int biglen, int certainty, random_engine & rand);
			
			static std::vector<int> make_positive(const std::vector<signed char>& a);
			static std::vector<int> make_positive(const std::vector<int>& a);
			static std::vector<int> trusted_strip_leading_zero_ints(const std::vector<int> & a);
			static std::vector<int> strip_leading_zero_bytes(const std::vector<signed char> & a);
			static std::vector<int> strip_leading_zero_ints(const std::vector<int> & a);
			static void destructive_mul_add(std::vector<int>& x, int y, int z);
		};



		big_int operator + (const big_int& a, const big_int& b);
		big_int operator - (const big_int& a, const big_int& b);
		big_int operator * (const big_int& a, const big_int& b);

	}


}

#endif // !RSA_H
