#include "rsa.h"

#include <vector>
#include <string>

namespace y
{

	namespace util
	{
		static const int MIN_RADIX = 2;
		static const int MAX_RADIX = 36;

		static const long long LONG_MASK = 0xffffffffL;

		const big_int big_int::ZERO = big_int();
		const big_int big_int::ONE = big_int(1L);

		big_int::big_int()
			: signum_(0)
		{}

		big_int::big_int(int val)
			: signum_(0)
		{
			if(val == 0)
				return;
			if (val < 0) {
				val = - val;
				signum_ = -1;
			}else{
				signum_ = 1;
			}
			mag_.resize(1);
			mag_[0] = val;
		}

		big_int::big_int(unsigned int val)
			: signum_(0)
		{
			if (val == 0)
				return;
			signum_ = 1;
			mag_.resize(1);
			mag_[0] = int(val);
		}

		big_int::big_int(long long val)
			: signum_(0)
		{
			if (val == 0)
				return;
			if (val < 0) {
				val = - val;
				signum_ = -1;
			}else{
				signum_ = 1;
			}

			int highword = (int)(unsigned long long(val) >> 32);
			if(highword == 0) {
				mag_.resize(1);
				mag_[0] = int(val);
			}else{
				mag_.resize(2);
				mag_[0] = highword;
				mag_[1] = int(val);
			}
		}

		big_int::big_int(unsigned long long val)
			: signum_(0)
		{
			if (val == 0)
				return;
			signum_ = 1;
			int highword = (int)(unsigned long long(val) >> 32);
			if(highword == 0) {
				mag_.resize(1);
				mag_[0] = int(val);
			}else{
				mag_.resize(2);
				mag_[0] = highword;
				mag_[1] = int(val);
			}
		}


		big_int::big_int(const std::vector<signed char>& val)
		{
			if (val.size() == 0)
				throw std::runtime_error("Zero length BigInteger");

			if (val[0] < 0) {
				mag_ = make_positive(val);
				signum_ = -1;
			} else {
				mag_ = strip_leading_zero_bytes(val);
				signum_ = (mag_.size() == 0 ? 0 : 1);
			}
		}

		big_int::big_int(const std::vector<int>& val )
		{
			if (val.size() == 0)
				throw std::runtime_error("Zero length BigInteger");

			if (val[0] < 0) {
				mag_ = make_positive(val);
				signum_ = -1;
			} else {
				mag_ = trusted_strip_leading_zero_ints(val);
				signum_ = (mag_.size() == 0 ? 0 : 1);
			}
		}

		big_int::big_int(int signum, const std::vector<signed char>& magnitude)
		{
			mag_ = strip_leading_zero_bytes(magnitude);

			if (signum < -1 || signum > 1)
				throw(std::runtime_error("Invalid signum value"));

			if (mag_.size()==0) {
				signum_ = 0;
			} else {
				if (signum == 0)
					throw(std::runtime_error("signum-magnitude mismatch"));
				signum_ = signum;
			}
		}

		big_int::big_int(int signum, const std::vector<int>& magnitude)
		{
			mag_ = strip_leading_zero_ints(magnitude);

			if (signum < -1 || signum > 1)
				throw(std::runtime_error("Invalid signum value"));

			if (mag_.size()==0) {
				signum_ = 0;
			} else {
				if (signum == 0)
					throw(std::runtime_error("signum-magnitude mismatch"));
				signum_ = signum;
			}
		}



		static const long long bits_per_digit[] = {
			0, 0,
			1024, 1624, 2048, 2378, 2648, 2875, 3072, 3247, 3402, 3543, 3672,
			3790, 3899, 4001, 4096, 4186, 4271, 4350, 4426, 4498, 4567, 4633,
			4696, 4756, 4814, 4870, 4923, 4975, 5025, 5074, 5120, 5166, 5210,
			5253, 5295
		};

		static const int digits_per_int[] = {
			0, 0, 30, 19, 15, 13, 11,
			11, 10, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6,
			6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5
		};

		static const int int_radix[] = {
			0x40000000, 0x4546b3db, 0x40000000, 0x48c27395, 0x159fd800,
			0x75db9c97, 0x40000000, 0x17179149, 0x3b9aca00, 0xcc6db61,
			0x19a10000, 0x309f1021, 0x57f6c100, 0xa2f1b6f,  0x10000000,
			0x18754571, 0x247dbc80, 0x3547667b, 0x4c4b4000, 0x6b5a6e1d,
			0x6c20a40,  0x8d2d931,  0xb640000,  0xe8d4a51,  0x1269ae40,
			0x17179149, 0x1cb91000, 0x23744899, 0x2b73a840, 0x34e63b41,
			0x40000000, 0x4cfa3cc1, 0x5c13d840, 0x6d91b519, 0x39aa400
		};

		big_int::big_int(const std::string & val, int radix)
		{
			int cursor = 0, numDigits;
			const int len = val.size();

			if (radix < MIN_RADIX || radix > MAX_RADIX)
				throw std::runtime_error("Radix out of range");
			if (len == 0)
				throw std::runtime_error("Zero length BigInteger");

			// Check for at most one leading sign
			int sign = 1;
			int index1 = val.find_last_of('-');
			int index2 = val.find_last_of('+');
			//if ((index1 + index2) <= -1) {
			if(index1 == std::string::npos || index2 == std::string::npos) {
				// No leading sign character or at most one leading sign character
				if (index1 == 0 || index2 == 0) {
					cursor = 1;
					if (len == 1)
						throw std::runtime_error("Zero length BigInteger");
				}
				if (index1 == 0)
					sign = -1;
			} else
				throw std::runtime_error("Illegal embedded sign character");

			// Skip leading zeros and compute number of digits in magnitude
			while (cursor < len &&
				//Character.digit(val.charAt(cursor), radix) == 0///
					(val[cursor] - '0') == 0
				)
				cursor++;
			if (cursor == len) {
				signum_ = 0;
				mag_.clear();
				return;
			}

			numDigits = len - cursor;
			signum_ = sign;

			// Pre-allocate array of expected size. May be too large but can
			// never be too small. Typically exact.
			int numBits = (int)((unsigned long long(numDigits * bits_per_digit[radix]) >> 10) + 1);
			int numWords = unsigned(numBits + 31) >> 5;
			std::vector<int> magnitude(numWords);

			// Process first (potentially short) digit group
			int firstGroupLen = numDigits % digits_per_int[radix];
			if (firstGroupLen == 0)
				firstGroupLen = digits_per_int[radix];
			//String group = val.substring(cursor, cursor += firstGroupLen);
			std::string group = val.substr(cursor, firstGroupLen); 
			cursor += firstGroupLen;
			magnitude[numWords - 1] = std::stoi(group, nullptr, radix);
			if (magnitude[numWords - 1] < 0)
				throw std::runtime_error("Illegal digit");

			// Process remaining digit groups
			int superRadix = int_radix[radix];
			int groupVal = 0;
			while (cursor < len) {
				//group = val.substring(cursor, cursor += digitsPerInt[radix]);
				group = val.substr(cursor, digits_per_int[radix]);
				cursor += digits_per_int[radix];
				//groupVal = Integer.parseInt(group, radix);
				groupVal = std::stoi(group, nullptr, radix);
				if (groupVal < 0)
					throw std::runtime_error("Illegal digit");
				destructive_mul_add(magnitude, superRadix, groupVal);
			}
			// Required for cases where the array was overallocated.
			mag_ = trusted_strip_leading_zero_ints(magnitude);
		}



		big_int::big_int(int numbits, random_engine & rand)
		{
			// TODO
			throw std::runtime_error("not implemented");
		}

		big_int::big_int(const big_int & i)
			: signum_(i.signum_), mag_(i.mag_)
		{}

		big_int & big_int::operator = (const big_int & i)
		{
			signum_ = i.signum_;
			mag_ = i.mag_;
			return *this;
		}




		big_int::big_int(const std::vector<int>& magnitude, int signum)
			: signum_(magnitude.size() == 0 ? 0 : signum), mag_(magnitude)
		{}






		int big_int::compare_magnitude( const big_int & val ) const
		{
			auto & m1 = mag_;
			int len1 = m1.size();
			auto & m2 = val.mag_;
			int len2 = m2.size();
			if (len1 < len2)
				return -1;
			if (len1 > len2)
				return 1;
			for (int i = 0; i < len1; i++) {
				int a = m1[i];
				int b = m2[i];
				if (a != b)
					return ((a & LONG_MASK) < (b & LONG_MASK)) ? -1 : 1;
			}
			return 0;
		}

		int big_int::compare_to( const big_int & val ) const
		{
			if (signum_ == val.signum_) {
				switch (signum_) {
				case 1:
					return compare_magnitude(val);
				case -1:
					return val.compare_magnitude(*this);
				default:
					return 0;
				}
			}
			return signum_ > val.signum_ ? 1 : -1;
		}


		bool big_int::is_zero() const
		{
			if((signum_ == 0) != (mag_.size() == 0))
				throw (std::runtime_error("signum-magnitude mismatch"));
			return signum_ == 0;
		}


		bool big_int::is_null() const
		{
			// TODO
			throw std::runtime_error("not implemented");
		}



		std::vector<signed char> big_int::to_byte_array() const
		{
			throw std::runtime_error("not implemented");
		}


		static std::vector<std::string> init_zero_strings()
		{
			std::vector<std::string> zeros(64);
			zeros[63] =
				"000000000000000000000000000000000000000000000000000000000000000";
			for (int i=0; i<63; i++)
				zeros[i] = zeros[63].substr(0, i);
			return zeros;
		}
		static const std::vector<std::string> zeros = init_zero_strings();

		static const int digits_per_longlong[] = {0, 0,
			62, 39, 31, 27, 24, 22, 20, 19, 18, 18, 17, 17, 16, 16, 15, 15, 15, 14,
			14, 14, 14, 13, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12};
		static const big_int longlong_radix[] = {0L, 0L,
			0x4000000000000000L, 0x383d9170b85ff80bL,
			0x4000000000000000L, 0x6765c793fa10079dL,
			0x41c21cb8e1000000L, 0x3642798750226111L,
			0x1000000000000000L, 0x12bf307ae81ffd59L,
			0xde0b6b3a7640000L, 0x4d28cb56c33fa539L,
			0x1eca170c00000000L, 0x780c7372621bd74dL,
			0x1e39a5057d810000L, 0x5b27ac993df97701L,
			0x1000000000000000L, 0x27b95e997e21d9f1L,
			0x5da0e1e53c5c8000L,  0xb16a458ef403f19L,
			0x16bcc41e90000000L, 0x2d04b7fdd9c0ef49L,
			0x5658597bcaa24000L,  0x6feb266931a75b7L,
			0xc29e98000000000L, 0x14adf4b7320334b9L,
			0x226ed36478bfa000L, 0x383d9170b85ff80bL,
			0x5a3c23e39c000000L,  0x4e900abb53e6b71L,
			0x7600ec618141000L,  0xaee5720ee830681L,
			0x1000000000000000L, 0x172588ad4f5f0981L,
			0x211e44f7d02c1000L, 0x2ee56725f06e5c71L,
			0x41c21cb8e1000000L
		};

		std::string big_int::to_string(int radix) const
		{
			if (signum_ == 0)
				return "0";
			if (radix < MIN_RADIX || radix > MAX_RADIX)
				radix = 10;

			// Compute upper bound on number of digit groups and allocate space
			int maxNumDigitGroups = (4*mag_.size() + 6)/7;
			std::vector<std::string> digitGroup(maxNumDigitGroups);

			// Translate number to string, a digit group at a time
			big_int tmp = abs();
			int numGroups = 0;

			char buffer[1024];
			while (tmp.signum_ != 0) {
				big_int d = longlong_radix[radix];

				//MutableBigInteger q = new MutableBigInteger(),
				//	a = new MutableBigInteger(tmp.mag),
				//	b = new MutableBigInteger(d.mag);
				//MutableBigInteger r = a.divide(b, q);
				//BigInteger q2 = q.toBigInteger(tmp.signum * d.signum);
				//BigInteger r2 = r.toBigInteger(tmp.signum * d.signum);
				big_int q;
				big_int r = tmp.divide(d, q);
				q.signum_ = tmp.signum_ * d.signum_;
				r.signum_ = tmp.signum_ * d.signum_;

				digitGroup[numGroups++] = itoa(r.to_longlong(), buffer, radix);
				tmp = q;
			}

			// Put sign (if any) and first digit group into result buffer
			//StringBuilder buf = new StringBuilder(numGroups*digitsPerLong[radix]+1);
			std::string buf;
			if (signum_<0)
				buf.append("-");
			buf.append(digitGroup[numGroups-1]);

			// Append remaining digit groups padded with leading zeros
			for (int i=numGroups-2; i>=0; i--) {
				// Prepend (any) leading zeros for this digit group
				int numLeadingZeros = digits_per_longlong[radix]-digitGroup[i].length();
				if (numLeadingZeros != 0)
					buf.append(zeros[numLeadingZeros]);
				buf.append(digitGroup[i]);
			}
			return buf;
		}


		long long big_int::to_longlong() const
		{
			throw std::runtime_error("not implemented");
		}

		int big_int::bit_length() const
		{
			throw std::runtime_error("not implemented");
		}

		y::util::big_int big_int::set_bit( int n ) const
		{
			throw std::runtime_error("not implemented");
		}







		void big_int::negate()
		{
			signum_ = -signum_;
		}

		big_int big_int::negated() const
		{
			return big_int(mag_, -signum_);
		}

		big_int big_int::abs() const
		{
			return (signum_ >= 0 ? *this : negated());
		}

		static int number_of_leading_zeros(int i)
		{
			if (i == 0)
				return 32;
			int n = 1;
			if (unsigned(i) >> 16 == 0) { n += 16; i <<= 16; }
			if (unsigned(i) >> 24 == 0) { n +=  8; i <<=  8; }
			if (unsigned(i) >> 28 == 0) { n +=  4; i <<=  4; }
			if (unsigned(i) >> 30 == 0) { n +=  2; i <<=  2; }
			n -= unsigned(i) >> 31;
			return n;
		}

		int big_int::divide_one_word( int divisor, big_int & quotient ) const
		{
			//long long divisorLong = divisor & LONG_MASK;
			//// Special case of one word dividend
			//if (mag_.size() == 1) {
			//	long dividendValue = mag_[0] & LONG_MASK;
			//	int q = (int) (dividendValue / divisorLong);
			//	int r = (int) (dividendValue - q * divisorLong);
			//	//quotient.mag_[0] = q;
			//	//quotient.intLen = (q == 0) ? 0 : 1;
			//	//quotient.offset = 0;
			//	if(q == 0){
			//		quotient.mag_.clear();
			//		quotient.signum_ = 0;
			//	}else{
			//		quotient.mag_.resize(1);
			//		quotient.mag_[0] = q;
			//	}
			//	return r;
			//}

			//// Normalize the divisor
			//int shift = number_of_leading_zeros(divisor);

			//int rem = mag_[0];
			//long remLong = rem & LONG_MASK;
			//if (remLong < divisorLong) {
			//	quotient.mag_[0] = 0;
			//} else {
			//	quotient.mag_[0] = (int)(remLong / divisorLong);
			//	rem = (int) (remLong - (quotient.mag_[0] * divisorLong));
			//	remLong = rem & LONG_MASK;
			//}

			//int xlen = intLen;
			//int[] qWord = new int[2];
			//while (--xlen > 0) {
			//	long dividendEstimate = (remLong<<32) |
			//		(value[offset + intLen - xlen] & LONG_MASK);
			//	if (dividendEstimate >= 0) {
			//		qWord[0] = (int) (dividendEstimate / divisorLong);
			//		qWord[1] = (int) (dividendEstimate - qWord[0] * divisorLong);
			//	} else {
			//		divWord(qWord, dividendEstimate, divisor);
			//	}
			//	quotient.value[intLen - xlen] = qWord[0];
			//	rem = qWord[1];
			//	remLong = rem & LONG_MASK;
			//}

			//quotient.normalize();
			//// Unnormalize
			//if (shift > 0)
			//	return rem % divisor;
			//else
			//	return rem;


			throw std::runtime_error("not implemented");
		}


		big_int big_int::divide( const big_int & b, big_int & quotient ) const
		{
			if(b.is_zero())
				throw std::runtime_error("BigInteger divide by zero");

			if(is_zero()) {
				quotient = ZERO;
				return *this;
			}

			int cmp = compare_magnitude(b);
			if(cmp < 0) {
				quotient = ZERO;
				return *this;
			}

			if(cmp == 0) {
				quotient = ONE;
				return ZERO;
			}

			// cmp > 0
			quotient.mag_.clear();
			if(b.mag_.size() == 1) {
				
			}


			throw std::runtime_error("not implemented");
		}

		y::util::big_int big_int::remainder( const big_int & i ) const
		{
			throw std::runtime_error("not implemented");
		}

		y::util::big_int big_int::gcd( const big_int & i ) const
		{
			throw std::runtime_error("not implemented");
		}

		y::util::big_int big_int::mod_inverse( const big_int & n ) const
		{
			throw std::runtime_error("not implemented");
		}

		y::util::big_int big_int::mod_power( const big_int & i, const big_int & n ) const
		{
			throw std::runtime_error("not implemented");
		}

		bool big_int::prime_to_certainty( int certainty, random_engine & rand ) const
		{
			throw std::runtime_error("not implemented");
		}
		






		static const int SMALL_PRIME_THRESHOLD = 95;
		static const int DEFAULT_PRIME_CERTAINTY = 100;

		big_int big_int::probable_prime( int bitlen, random_engine & rand )
		{
			if(bitlen < 2)
				throw std::runtime_error("biglen < 2");
			return (bitlen < SMALL_PRIME_THRESHOLD ? small_prime(bitlen, DEFAULT_PRIME_CERTAINTY, rand) :
				large_prime(bitlen, DEFAULT_PRIME_CERTAINTY, rand));
		}
		
		static const big_int SMALL_PRIME_PRODUCT((long long)(3)*5*7*11*13*17*19*23*29*31*37*41);

		big_int big_int::small_prime( int bitLength, int certainty, random_engine & rnd )
		{
			int magLen = unsigned(bitLength + 31) >> 5;
			std::vector<int> temp(magLen);
			int highBit = 1 << ((bitLength+31) & 0x1f);  // High bit of high int
			int highMask = (highBit << 1) - 1;  // Bits to keep in high int

			while(true) {
				// Construct a candidate
				for (int i=0; i<magLen; i++)
					temp[i] = rnd();
				temp[0] = (temp[0] & highMask) | highBit;  // Ensure exact length
				if (bitLength > 2)
					temp[magLen-1] |= 1;  // Make odd if bitlen > 2

				big_int p(temp, 1);

				// Do cheap "pre-test" if applicable
				if (bitLength > 6) {
					long long r = p.remainder(SMALL_PRIME_PRODUCT).to_longlong();
					if ((r%3==0)  || (r%5==0)  || (r%7==0)  || (r%11==0) ||
						(r%13==0) || (r%17==0) || (r%19==0) || (r%23==0) ||
						(r%29==0) || (r%31==0) || (r%37==0) || (r%41==0))
						continue; // Candidate is composite; try another
				}

				// All candidates of bitLength 2 and 3 are prime by this point
				if (bitLength < 4)
					return p;

				// Do expensive test if we survive pre-test (or it's inapplicable)
				if (p.prime_to_certainty(certainty, rnd))
					return p;
			}
		}

		big_int big_int::large_prime( int bitLength, int certainty, random_engine & rand )
		{
			big_int p(bitLength, rand);
			p = p.set_bit(bitLength-1);
			p.mag_[p.mag_.size()-1] &= 0xfffffffe;

			// Use a sieve length likely to contain the next prime number
			int searchLen = (bitLength / 20) * 64;
			bit_sieve searchSieve(p, searchLen);
			big_int candidate = searchSieve.retrieve(p, certainty, rand);

			while ((candidate.is_null()) || (candidate.bit_length() != bitLength)) {
				p = p + big_int(2*searchLen);
				if (p.bit_length() != bitLength){
					p = big_int(bitLength, rand);
					p = p.set_bit(bitLength-1);
				}
				p.mag_[p.mag_.size()-1] &= 0xfffffffe;
				searchSieve = bit_sieve(p, searchLen);
				candidate = searchSieve.retrieve(p, certainty, rand);
			}
			return candidate;
		}

		std::vector<int> big_int::make_positive( const std::vector<signed char>& a )
		{
			int keep, k;
			int byteLength = a.size();

			// Find first non-sign (0xff) byte of input
			for (keep=0; keep<byteLength && a[keep]==-1; keep++)
				;


			/* Allocate output array.  If all non-sign bytes are 0x00, we must
			* allocate space for one extra output byte. */
			for (k=keep; k<byteLength && a[k]==0; k++)
				;

			int extraByte = (k==byteLength) ? 1 : 0;
			int intLength = ((byteLength - keep + extraByte) + 3)/4;
			std::vector<int> result(intLength);

			/* Copy one's complement of input into output, leaving extra
			* byte (if it exists) == 0x00 */
			int b = byteLength - 1;
			for (int i = intLength-1; i >= 0; i--) {
				result[i] = a[b--] & 0xff;
				int numBytesToTransfer = std::min(3, b-keep+1);
				if (numBytesToTransfer < 0)
					numBytesToTransfer = 0;
				for (int j=8; j <= 8*numBytesToTransfer; j += 8)
					result[i] |= ((a[b--] & 0xff) << j);

				// Mask indicates which bits must be complemented
				int mask = unsigned(-1) >> (8*(3-numBytesToTransfer));
				result[i] = ~result[i] & mask;
			}

			// Add one to one's complement to generate two's complement
			for (int i=result.size()-1; i>=0; i--) {
				result[i] = (int)((result[i] & LONG_MASK) + 1);
				if (result[i] != 0)
					break;
			}

			return result;
		}

		std::vector<int> big_int::make_positive( const std::vector<int>& a )
		{
			int keep, j;

			// Find first non-sign (0xffffffff) int of input
			for (keep=0; keep<a.size() && a[keep]==-1; keep++)
				;

			/* Allocate output array.  If all non-sign ints are 0x00, we must
			* allocate space for one extra output int. */
			for (j=keep; j<a.size() && a[j]==0; j++)
				;
			int extraInt = (j==a.size()? 1 : 0);
			std::vector<int> result(a.size() - keep + extraInt);

			/* Copy one's complement of input into output, leaving extra
			* int (if it exists) == 0x00 */
			for (int i = keep; i<a.size(); i++)
				result[i - keep + extraInt] = ~a[i];

			// Add one to one's complement to generate two's complement
			for (int i=result.size()-1; ++result[i]==0; i--)
				;

			return result;
		}

		std::vector<int> big_int::trusted_strip_leading_zero_ints( const std::vector<int> & val )
		{
			int vlen = val.size();
			int keep;

			// Find first nonzero byte
			for (keep = 0; keep < vlen && val[keep] == 0; keep++)
				;
			if(keep == 0)
				return val;
			return std::vector<int>(val.begin() + keep, val.end());
			//return keep == 0 ? val : java.util.Arrays.copyOfRange(val, keep, vlen);
		}

		std::vector<int> big_int::strip_leading_zero_bytes( const std::vector<signed char> & a )
		{
			int byteLength = a.size();
			int keep;

			// Find first nonzero byte
			for (keep = 0; keep < byteLength && a[keep]==0; keep++)
				;

			// Allocate new array and copy relevant part of input array
			int intLength = unsigned((byteLength - keep) + 3) >> 2;
			std::vector<int> result(intLength);
			int b = byteLength - 1;
			for (int i = intLength-1; i >= 0; i--) {
				result[i] = a[b--] & 0xff;
				int bytesRemaining = b - keep + 1;
				int bytesToTransfer = std::min(3, bytesRemaining);
				for (int j=8; j <= (bytesToTransfer << 3); j += 8)
					result[i] |= ((a[b--] & 0xff) << j);
			}
			return result;
		}

		std::vector<int> big_int::strip_leading_zero_ints( const std::vector<int> & val )
		{
			int vlen = val.size();
			int keep;

			// Find first nonzero byte
			for (keep = 0; keep < vlen && val[keep] == 0; keep++)
				;
			//return java.util.Arrays.copyOfRange(val, keep, vlen);
			return std::vector<int>(val.begin() + keep, val.end());
		}

		void big_int::destructive_mul_add( std::vector<int>& x, int y, int z )
		{
			// Perform the multiplication word by word
			long long ylong = y & LONG_MASK;
			long long zlong = z & LONG_MASK;
			int len = x.size();

			long long product = 0;
			long long carry = 0;
			for (int i = len-1; i >= 0; i--) {
				product = ylong * (x[i] & LONG_MASK) + carry;
				x[i] = (int)product;
				carry = unsigned long long(product) >> 32;
			}

			// Perform the addition
			long long sum = (x[len-1] & LONG_MASK) + zlong;
			x[len-1] = (int)sum;
			carry = unsigned long long(sum) >> 32;
			for (int i = len-2; i >= 0; i--) {
				sum = (x[i] & LONG_MASK) + carry;
				x[i] = (int)sum;
				carry = unsigned long long(sum) >> 32;
			}
		}




		big_int operator + (const big_int& a, const big_int& b)
		{
			throw std::runtime_error("not implemented");
		}

		big_int operator - (const big_int& a, const big_int& b)
		{
			throw std::runtime_error("not implemented");
		}

		big_int operator * (const big_int& a, const big_int& b)
		{
			throw std::runtime_error("not implemented");
		}



























		bit_sieve::bit_sieve( const big_int & base, int searchlen )
		{
			throw std::runtime_error("not implemented");
		}

		y::util::big_int bit_sieve::retrieve( const big_int & initvalue, int certainty, random_engine & random )
		{
			throw std::runtime_error("not implemented");
		}

	}


}