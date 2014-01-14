#ifndef _BIG_NUM_H_
#define _BIG_NUM_H_

#include <array>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdint>
#include <bitset>

namespace y
{

	namespace math
	{
		using std::operator<<;
		using std::operator>>;

		static const unsigned long long LONG_MASK = 0xffffffffL;

		namespace __private__
		{
			template <class _ritera, class _riterb>
			int compare_range(_ritera arbegin, _ritera arend, _riterb brbegin)
			{
				for (auto i = arbegin, auto j = brbegin; i != arend; ++i, ++j) {
					if (*i != *j)
						return (*i < *j) ? -1 : 1;
				}
				return 0;
			}

			// ~a
			template <class _itera>
			void twos_complement_range(_itera begin, _itera end)
			{
				for(auto i = begin; i != end; ++i)
					*i = ~(*i);
			}

			// -a
			template <class _itera>
			void negate_range(_itera begin, _itera end)
			{
				twos_complement_range(begin, end);
				plus_one_range(begin, end);
			}


			// a++
			template <class _itera>
			void plus_one_range(_itera begin, _itera end)
			{
				bool carry = true;
				for(auto i = begin; i != end; ++i){
					unsigned long long sumhere = (LONG_MASK & *i) + (carry ? 1 : 0);
					*i = unsigned(sumhere);
					carry = (sumhere >> 32) != 0;
					if(!carry)
						break;
				}
			}

			// a += b
			template <class _itera, class _iterb>
			void plus_range(_itera abegin, _itera aend, _iterb bbegin)
			{
				bool carry = false;
				_itera ai;
				_iterb bi;
				for(ai = abegin, bi = bbegin; ai != aend; ++ai, ++bi){
					unsigned long long sumhere = 
						(LONG_MASK & *ai) + (LONG_MASK & *bi) + (carry ? 1 : 0);
					*ai = unsigned(sumhere);
					carry = (sumhere >> 32) != 0;
				}
			}

			// a -= b
			template <class _itera, class _iterb>
			void subtract_range(_itera abegin, _itera aend, _iterb bbegin)
			{
				std::vector<std::iterator_traits<_itera>::value_type> temp(bbegin, bbegin + (aend - abegin));
				negate_range(temp.begin(), temp.end());
				plus_range(abegin, aend, temp.begin());
			}


			// ab = a * b mod
			template <class _itera, class _iterb, class _iterab>
			void product_range(_itera abegin, _itera aend, _iterb bbegin, _iterab abbegin)
			{
				auto bend = bbegin + (aend - abegin);
				auto abend = abbegin + (aend - abegin);
				for(auto bi = bbegin; bi != bend; ++bi){
					unsigned long long carry = 0;
					_iterab abi = abbegin + (bi - bbegin);
					for(_itera ai = abegin; 
						abi != abend; ++ ai, ++abi){
						unsigned long long product = (*ai & LONG_MASK) * (*bi & LONG_MASK) + 
							(*abi & LONG_MASK) + carry;
						*abi = (unsigned)product;
						carry = product >> 32;
					}
				}
			}

		}

		template <int num>
		class fixed_unsigned_int
		{
			static_assert(num > 0, "num should be positive");			
		public:
			typedef uint32_t unit_t;
			typedef uint64_t duint_t;
			static const int unit_num = num;
			
			fixed_unsigned_int(){std::fill(mag_.begin(), mag_.end(), 0);}
			fixed_unsigned_int(unsigned long long val) 
			{
				std::fill(mag_.begin(), mag_.end(), 0);
				unsigned highword = (unsigned)(unsigned long long(val) >> 32);
				if(highword == 0) {
					mag_[0] = unsigned(val);
				}else{
					mag_[1] = highword;
					mag_[0] = unsigned(val);
				}
			}
			explicit fixed_unsigned_int(const std::array<unsigned, unit_num> & val)
				: mag_(val)
			{}
			fixed_unsigned_int(const fixed_unsigned_int& v)
				: mag_(v.mag_)
			{}
			template<int num2> 
			explicit fixed_unsigned_int(const fixed_unsigned_int<num2> & v)
			{
				std::fill(mag_.begin(), mag_.end(), 0);
				std::copy(v.mag_.begin(), v.mag_.end(), mag_.begin());
			}
			template<class _it>
			explicit fixed_unsigned_int(_it begin, _it end)
			{
				std::fill(mag_.begin(), mag_.end(), 0);
				std::copy(begin, end, mag_.begin());
			}

			inline const std::array<unsigned, unit_num> & magnitude() const {return mag_;}
			inline unsigned mag_at(int i) const {return mag_[i];}
			inline unsigned& mag_at(int i) {return mag_[i];}
			inline unsigned operator [] (int i) const {return mag_[i];}
			inline unsigned& operator [] (int i) {return mag_[i];}

			operator unsigned () const {return mag_[0];}
			operator unsigned long long () const {return (mag_[0] & LONG_MASK) + ((mag_[1] & LONG_MASK) << 32);}

		public:
			// comparison
			int compare_to(const fixed_unsigned_int & v) const;
			inline bool operator == (const fixed_unsigned_int & v) const {return compare_to(v) == 0;}
			inline bool operator != (const fixed_unsigned_int & v) const {return compare_to(v) != 0;}
			inline bool operator < (const fixed_unsigned_int & v) const {return compare_to(v) == -1;}
			inline bool operator > (const fixed_unsigned_int & v) const {return compare_to(v) == 1;}
			inline bool operator <= (const fixed_unsigned_int & v) const {return compare_to(v) <= 0;}
			inline bool operator >= (const fixed_unsigned_int & v) const {return compare_to(v) >= 0;}
			bool is_zero() const;
			int most_significant_nonzero_unit() const;
			int most_significant_nonzero_bit() const;
			void left_shift(unsigned n);
			void right_shift(unsigned n);
			fixed_unsigned_int left_shifted(unsigned n) const {fixed_unsigned_int s(*this); s.left_shift(n); return s;}
			fixed_unsigned_int right_shifted(unsigned n) const {fixed_unsigned_int s(*this); s.right_shift(n); return s;}
			

			void twos_complement();
			void negate();
			void plus_one();
			void subtract_one();

			fixed_unsigned_int operator ++ (int) {fixed_unsigned_int o(*this); plus_one(); return o;}
			fixed_unsigned_int & operator ++ () {plus_one(); return *this;}
			fixed_unsigned_int operator -- (int) {fixed_unsigned_int o(*this); subtract_one(); return o;}
			fixed_unsigned_int & operator -- () {subtract_one(); return *this;}
			
			fixed_unsigned_int & operator += (const fixed_unsigned_int & v);
			fixed_unsigned_int & operator -= (const fixed_unsigned_int & v);
			fixed_unsigned_int & operator *= (const fixed_unsigned_int & v);

		private:
			std::array<unsigned, unit_num> mag_;
		};

		template <int num>
		int fixed_unsigned_int<num>::compare_to( const fixed_unsigned_int & v ) const
		{
			//for (int i = unit_num-1; i >= 0; i--) {
			//	if (mag_[i] != v.mag_[i])
			//		return (mag_[i] < v.mag_[i]) ? -1 : 1;
			//}
			//return 0;
			return __private__::compare_range(mag_.rbegin(), mag_.rend(), v.mag_.rbegin());
		}


		template <int num>
		bool fixed_unsigned_int<num>::is_zero() const
		{
			for (int i = unit_num-1; i >= 0; i--) {
				if (mag_[i] != 0)
					return false;
			}
			return true;
		}

		template <int num>
		int fixed_unsigned_int<num>::most_significant_nonzero_unit() const
		{
			for (int i = unit_num-1; i >= 0; i--) {
				if (mag_[i] != 0)
					return i;
			}
			return -1;
		}

		template <int num>
		int fixed_unsigned_int<num>::most_significant_nonzero_bit() const
		{
			int msn = most_significant_nonzero_unit();
			if(msn == -1)
				return -1;
			unsigned b = mag_[msn];
			int a = -1;
			for(; b != 0; a++)
				b >>= 1;
			assert(a >= 0);
			return msn * 32 + a;
		}




		template <int num>
		void fixed_unsigned_int<num>::left_shift( unsigned n )
		{
			unsigned carry = 0;
			for(int i = 0; i < num; i++) {
				unsigned long long s = ((mag_[i] & LONG_MASK) << n) + carry;
				mag_[i] = unsigned(s);
				carry = (s >> 32);
			}
		}


		template <int num>
		void fixed_unsigned_int<num>::right_shift( unsigned n )
		{
			unsigned carry = 0;
			for(int i = num-1; i >= 0; i--) {
				unsigned long long s = ((mag_[i] & LONG_MASK) + ((carry & LONG_MASK) << 32));
				mag_[i] = unsigned(s >> n);
				carry = unsigned(s - ((mag_[i] & LONG_MASK) << n));
			}
		}

		/*
		template <int num>
		inline fixed_unsigned_int<num> operator << (const fixed_unsigned_int<num>& v, unsigned n)
		{
			fixed_unsigned_int<num> s(v); 
			s.shift_left(n); 
			return s;
		}

		template <int num>
		inline fixed_unsigned_int<num> operator >> (const fixed_unsigned_int<num>& v, unsigned n) 
		{
			fixed_unsigned_int<num> s(v); 
			s.shift_right(n); 
			return s;
		}
		*/


		template <int num>
		void fixed_unsigned_int<num>::twos_complement()
		{
			/*
			for(int i = 0; i < unit_num; i++)
				mag_[i] = ~mag_[i];
				*/
			__private__::twos_complement_range(mag_.begin(), mag_.end());
		}

		template <int num>
		void fixed_unsigned_int<num>::negate()
		{
			//twos_complement();
			//plus_one();
			__private__::negate_range(mag_.begin(), mag_.end());
		}

		template <int num>
		void fixed_unsigned_int<num>::plus_one()
		{
			/*
			bool carry = true;
			for(int i = 0; i < unit_num; i++){
				unsigned long long sumhere = (LONG_MASK & mag_[i]) + (carry ? 1 : 0);
				mag_[i] = unsigned(sumhere);
				carry = (sumhere >> 32) != 0;
				if(!carry)
					break;
			}*/
			__private__::plus_one_range(mag_.begin(), mag_.end());
		}

		template <int num>
		void fixed_unsigned_int<num>::subtract_one()
		{
			bool borrow = true;
			for(int i = 0; i < unit_num; i++){
				long long subhere = (LONG_MASK & mag_[i]) - (borrow ? 1 : 0);
				mag_[i] = unsigned(subhere);
				borrow = subhere < 0;
				if(!borrow)
					break;
			}
		}


		template <int num>
		fixed_unsigned_int<num> & fixed_unsigned_int<num>::operator+=( const fixed_unsigned_int & v )
		{
			/*
			bool carry = false;
			for(int i = 0; i < unit_num; i++){
				unsigned long long sumhere = (LONG_MASK & mag_[i]) + (LONG_MASK & v.mag_[i]) + (carry ? 1 : 0);
				mag_[i] = unsigned(sumhere);
				carry = (sumhere >> 32) != 0;
			}
			return *this;*/
			__private__::plus_range(mag_.begin(), mag_.end(), v.mag_.begin());
			return *this;
		}

		template <int num>
		fixed_unsigned_int<num> & fixed_unsigned_int<num>::operator-=( const fixed_unsigned_int & v )
		{
			fixed_unsigned_int s(v);
			s.negate();
			return *this += s;
		}

		template <int num>
		fixed_unsigned_int<num> & fixed_unsigned_int<num>::operator*=( const fixed_unsigned_int & v )
		{
			*this = (*this * v);
			return *this;
		}




		template <int num>
		fixed_unsigned_int<num> operator ~ (const fixed_unsigned_int<num> & v)
		{
			fixed_unsigned_int<num> s(v);
			s.twos_complement();
			return s;
		}

		template <int num>
		fixed_unsigned_int<num> operator - (const fixed_unsigned_int<num> & v)
		{
			fixed_unsigned_int<num> s(v);
			s.negate();
			return s;
		}

		template <int num>
		fixed_unsigned_int<num> operator + (const fixed_unsigned_int<num> & a, const fixed_unsigned_int<num> & b)
		{
			fixed_unsigned_int<num> s(a);
			s += b;
			return s;
		}

		template <int num>
		fixed_unsigned_int<num> operator - (const fixed_unsigned_int<num> & a, const fixed_unsigned_int<num> & b)
		{
			fixed_unsigned_int<num> s(a);
			s -= b;
			return s;
		}

		// unwrapped
		template <int num1, int num2>
		fixed_unsigned_int<num1 + num2> mult_unlimit(const fixed_unsigned_int<num1>& xx, const fixed_unsigned_int<num2>& yy)
		{
			fixed_unsigned_int<num1 + num2> p;
			for(int j = 0; j < num2; j++){
				unsigned long long carry = 0;
				for(int i = 0; i < num1; i++){
					int k = i + j;
					unsigned long long product = (xx[i] & LONG_MASK) * (yy[j] & LONG_MASK) + (p[k] & LONG_MASK) + carry;
					p[k] = (unsigned)product;
					carry = product >> 32;
				}
				p[j + num1] = carry;
			}
			return p;
		}

		// wrapped
		template <int num>
		fixed_unsigned_int<num> operator * (const fixed_unsigned_int<num>& xx, const fixed_unsigned_int<num>& yy)
		{
			/*
			fixed_unsigned_int<num> p;
			for(int j = 0; j < num; j++){
				unsigned long long carry = 0;
				for(int i = 0; i + j < num; i++){
					int k = i + j;
					unsigned long long product = (xx[i] & LONG_MASK) * (yy[j] & LONG_MASK) + (p[k] & LONG_MASK) + carry;
					p[k] = (unsigned)product;
					carry = product >> 32;
				}
			}
			return p;*/
			std::array<unsigned, num> p;
			std::fill(p.begin(), p.end(), 0);
			__private__::product_range(xx.magnitude().begin(), xx.magnitude().end(), yy.magnitude().begin(), p.begin());
			return fixed_unsigned_int<num>(p);
		}



		// the long division
		template <int num>
		fixed_unsigned_int<num> divide( const fixed_unsigned_int<num> & n, 
			const fixed_unsigned_int<num> & d, 
			fixed_unsigned_int<num> & r )
		{
			int msn = d.most_significant_nonzero_unit();
			if(msn == -1){
				throw std::runtime_error("divided by zero");
			}

			//if(msn == 0){
			//	throw std::runtime_error("not implemented");
			//}

			fixed_unsigned_int quotient, divisor(d);
			//rem = n;

			
			fixed_unsigned_int<num+1> rem(n);

			int dmsu = divisor[msn]; // most significant unit of the divisor
			assert(dmsu);							
			int shift = 0;
			for(; dmsu!=0; shift++)
				dmsu >= 1;
			shift = 32 - shift;
			
			// move most significant bit to the left end
			divisor.left_shift(shift);
			rem.left_shift(shift);
			
			
			// [msn ... 0]
			for(int i = num - msn; i >= 0; --i) {
				int rem_from = i + msn;
				int rem_to = i;
				// rem[rem_from ... rem_to] / divisor[msn ... 0]
				bool remhereislarger = true;
				while(remhereislarger){
					remhereislarger = true;
					for(int k = 0; k <= msn; k++){
						unsigned remunit = rem[rem_from - k];
						unsigned divuint = divisor[msn - k];
						if(remunit != divuint){
							remhereislarger = remunit > divuint;
							break;
						}
					}
					
					if(remhereislarger){ // rem[rem_from ... rem_to] -= divisor[msn ... 0]
						// += ~ divisor[msn ... 0]
						bool carry = false;
						for(int remid = rem_to, divid = 0; remid <= rem_from; remid++, divid++){
							unsigned long long sumhere = 
								(LONG_MASK & rem[remid]) + (LONG_MASK & ~divisor[divid]) + (carry ? 1 : 0);
							rem[remid] = unsigned(sumhere);
							carry = (sumhere >> 32) != 0;
						}
						// += 1
						carry = true;
						for(int remid = rem_to; remid <= rem_from; remid++){
							unsigned long long sumhere = (LONG_MASK & rem[remid]) + (carry ? 1 : 0);
							rem[remid] = unsigned(sumhere);
							carry = (sumhere >> 32) != 0;
							if(!carry)
								break;
						}
						quotient[i] ++;
					}
				}
			}

		}





		// prime
		static const int DEFAULT_PRIME_CERTAINTY = 100;
		template <int num, class randgen_t>
		fixed_unsigned_int<num> gen_prime_with_certainty(int certainty, randgen_t rnd)
		{
			fixed_unsigned_int<num> p;
			throw std::runtime_error("not implemented");
		}



		template <int num>
		fixed_unsigned_int<num> gcd(const fixed_unsigned_int<num>& a, const fixed_unsigned_int<num>& b)
		{
			const fixed_unsigned_int<num> * ap = &a;
			const fixed_unsigned_int<num> * bp = &b;
			throw std::runtime_error("not implemented");
		}





		template <int num>
		std::ostream & operator << (std::ostream & s, const fixed_unsigned_int<num> & v)
		{
			s << std::hex << std::setfill('0');
			s << '[';
			for(int i = num-1; i > 0; i--)
				s << std::setw(2 * sizeof(unsigned)) << v.magnitude()[i] << '|';
			s << std::setw(2 * sizeof(unsigned)) << v.magnitude()[0] << ']';
			s << std::dec;
			s << std::setfill(' ');
			return s;
		}

		template <int num>
		std::istream & operator >> (std::istream & s, fixed_unsigned_int<num> & v)
		{
			std::cout << "input from lower bits to higher bits" << std::endl;
			std::array<unsigned, num> magnitude;
			std::fill(magnitude.begin(), magnitude.end(), 0);
			for(int i = 0; i < num; i++)
				s >> std::hex >> magnitude[i];
			v = fixed_unsigned_int<num>(magnitude);
			return s;
		}

	}


}


#endif // !_BIG_NUM_H_
