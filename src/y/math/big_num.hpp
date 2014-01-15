// copyright yanghao 2013

#ifndef _BIG_NUM_H_
#define _BIG_NUM_H_

#include <array>
#include <vector>
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

			template <class _iter>
			struct is_unsigned_int_iterator 
				: std::is_same<typename std::iterator_traits<_iter>::value_type, unsigned>
			{};

			template <class _itera>
			_itera find_least_significant_nonzero(_itera abegin, _itera aend)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				while(abegin != aend){
					if(*abegin != 0)
						return abegin;
					abegin ++;
				}
				return aend;
			}

			template <class _itera>
			_itera find_most_significant_nonzero(_itera abegin, _itera aend)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				std::reverse_iterator<_itera> arbegin(aend), arend(abegin);
				while(arbegin != arend){
					if(*arbegin != 0)
						return --arbegin.base();
					arbegin ++;
				}
				return aend;
			}

			template <class _itera, class _iterb>
			int compare_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				static_assert(is_unsigned_int_iterator<_iterb>::value, "invalid iterator type");
				int result = 0;
				auto i = abegin;
				auto j = bbegin;
				for (; i != aend && j != bend; ++i, ++j) {
					if (*i != *j)
						result = (*i < *j) ? -1 : 1;
				}
				if(result != 0)
					return result;
				if(i != aend && j == bend){
					for(; i != aend; ++i){
						if(*i != 0)
							return 1;
					}
				}
				if(i == aend && j != bend){
					for(; j != bend; ++j){
						if(*j != 0)
							return -1;
					}
				}
				return 0;
			}

			// ~a
			template <class _itera>
			void twos_complement_range(_itera begin, _itera end)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				for(auto i = begin; i != end; ++i)
					*i = ~(*i);
			}

			// -a
			template <class _itera>
			void negate_range(_itera begin, _itera end)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				twos_complement_range(begin, end);
				plus_one_range(begin, end);
			}


			// ++a
			template <class _itera>
			void plus_one_range(_itera begin, _itera end)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				bool carry = true;
				for(auto i = begin; i != end; ++i){
					unsigned long long sumhere = (LONG_MASK & *i) + (carry ? 1 : 0);
					*i = unsigned(sumhere);
					carry = (sumhere >> 32) != 0;
					if(!carry)
						break;
				}
			}

			// --a
			template <class _itera>
			void subtract_one_range(_itera begin, _itera end)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				bool borrow = true;
				for(auto i = begin; i != end; ++i){
					long long subhere = (LONG_MASK & *i) - (borrow ? 1 : 0);
					*i = unsigned(subhere);
					borrow = subhere < 0;
					if(!borrow)
						break;
				}
			}

			// a += b
			template <class _itera, class _iterb>
			void plus_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				static_assert(is_unsigned_int_iterator<_iterb>::value, "invalid iterator type");
				bool carry = false;
				_itera ai;
				_iterb bi;
				for(ai = abegin, bi = bbegin; ai != aend && bi != bend; ++ai, ++bi){
					unsigned long long sumhere = 
						(LONG_MASK & *ai) + (LONG_MASK & *bi) + (carry ? 1 : 0);
					*ai = unsigned(sumhere);
					carry = (sumhere >> 32) != 0;
				}
				if(ai != aend && carry)
					*ai += 1;
				assert(!(ai == aend && bi != bend)); // b should be shorter
			}

			// a -= b
			/*
			template <class _itera, class _iterb>
			void subtract_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				static_assert(is_unsigned_int_iterator<_iterb>::value, "invalid iterator type");
				std::vector<std::iterator_traits<_itera>::value_type> temp(bbegin, bend);
				negate_range(temp.begin(), temp.end());
				plus_range(abegin, aend, temp.begin(), temp.end());
				if(aend - abegin > bend - bbegin)
					subtract_one_range(abegin + (bend - bbegin), aend);
				// BUGGY
			}*/
			
			template <class _itera, class _iterb>
			void subtract_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				static_assert(is_unsigned_int_iterator<_iterb>::value, "invalid iterator type");
				bool borrow = false;
				_itera ai;
				_iterb bi;
				for(ai = abegin, bi = bbegin; ai != aend && bi != bend; ++ai, ++bi){
					unsigned long long along = LONG_MASK & *ai;
					unsigned long long blong = LONG_MASK & *bi;
					if(along >= blong + (borrow ? 1 : 0)){
						*ai = unsigned(along - blong - (borrow ? 1 : 0));
						borrow = false;						
					}else{
						*ai = unsigned(along + (1UL << 32U) - blong - (borrow ? 1 : 0));
						borrow = true;						
					}
				}
				//assert(!(ai == aend && bi != bend)); // b should be shorter
				if(borrow){ // still borrow
					while(ai != aend && *ai == 0){
						*ai = *ai - 1; // 0 -> 0xffffffff
					}
					if(ai != aend)
						*ai = *ai - 1;
				}
			}


			// ab = a * b
			template <class _itera, class _iterb, class _iterab>
			void mult_range(_itera abegin, _itera aend, _iterb bbegin, _iterb bend, _iterab abbegin, _iterab abend)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				static_assert(is_unsigned_int_iterator<_iterb>::value, "invalid iterator type");
				static_assert(is_unsigned_int_iterator<_iterab>::value, "invalid iterator type");
				//auto bend = bbegin + (aend - abegin);
				_iterab abihere = abbegin;
				for(auto bi = bbegin; bi != bend; ++bi, ++abihere){
					unsigned long long carry = 0;
					_itera ai = abegin;
					_iterab abi = abihere;
					for(; abi != abend && ai != aend; ++ai, ++abi){
						unsigned long long product = (*ai & LONG_MASK) * (*bi & LONG_MASK) + 
							(*abi & LONG_MASK) + carry;
						*abi = (unsigned)product;
						carry = product >> 32;
					}
					if(abi != abend)
						*abi = unsigned(carry);
				}
			}

			template <class _itera>
			void left_shift_range(_itera abegin, _itera aend, unsigned n)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				assert(n <= 32);
				unsigned carry = 0;
				for(auto i = abegin; i != aend; ++i) {
					unsigned long long s = ((*i & LONG_MASK) << n) + carry;
					*i = unsigned(s);
					carry = (s >> 32);
				}
			}

			template <class _itera>
			void right_shift_range(_itera abegin, _itera aend, unsigned n)
			{
				static_assert(is_unsigned_int_iterator<_itera>::value, "invalid iterator type");
				assert(n <= 32);
				std::reverse_iterator<_itera> arbegin(aend), arend(abegin);
				unsigned carry = 0;
				for(auto i = arbegin; i != arend; ++i) {
					unsigned long long s = ((*i & LONG_MASK) + ((carry & LONG_MASK) << 32));
					*i = unsigned(s >> n);
					carry = unsigned(s - ((*i & LONG_MASK) << n));
				}
			}


			// quotient = dividend / divisor, didivend -= quotient * divisor (remainder)
			template <class _iterDividend, class _iterDivisor, class _iterQuotient>
			void divide_range(_iterDividend ddbegin, _iterDividend ddend,
				_iterDivisor dsbegin, _iterDivisor dsend, 
				_iterQuotient qbegin, _iterQuotient qend)
			{
				static_assert(is_unsigned_int_iterator<_iterDividend>::value, "invalid iterator type");
				static_assert(is_unsigned_int_iterator<_iterDivisor>::value, "invalid iterator type");
				static_assert(is_unsigned_int_iterator<_iterQuotient>::value, "invalid iterator type");

				int cmp = compare_range(ddbegin, ddend, dsbegin, dsend);
				if(cmp < 0) { // dd < ds
					std::fill(qbegin, qend, 0u);
					return;
				}

				_iterDivisor dsbig = find_most_significant_nonzero(dsbegin, dsend);
				_iterDividend ddbig = find_most_significant_nonzero(ddbegin, ddend);
				if(dsbig == dsend) { // ds = 0
					throw std::runtime_error("divided by zero");
				}
				if(ddbig == ddend) { // dd = 0
					std::fill(qbegin, qend, 0u);
					return;
				}
				unsigned long long dshead = *dsbig & LONG_MASK;
				
				// 0 ... [0? ddtoclear ... ddcurrenton] ddcurrenton-1 ddcurrenton-2 ... ddbegin
				_iterDividend ddcurrenton = ddbig - (dsbig - dsbegin);
				typedef std::reverse_iterator<_iterDividend> _riterDividend;
				_riterDividend ddcurrenton_rev(ddcurrenton+1);
				_iterDividend ddtoclear = ddbig; // will be cleared
				_riterDividend ddtoclear_rev(ddtoclear+1);
				_riterDividend ddrevend(ddbegin);

				std::vector<unsigned> dsmultgess(dsend - dsbegin + 1, 0u); // cache for multiplication
				while (ddcurrenton_rev != ddrevend) {
					ddcurrenton = ddcurrenton_rev.base()-1;
					ddtoclear = ddtoclear_rev.base()-1;

					int cmp = compare_range(ddcurrenton, ddend, dsbegin, dsend);
					if(cmp >= 0){
						assert(ddtoclear - ddcurrenton == dsbig - dsbegin);
						unsigned long long ddtoclearhead = *ddtoclear & LONG_MASK;
						if(ddtoclear + 1 != ddend){
							ddtoclearhead += (*(ddtoclear+1) << 32);
						}
						assert(ddtoclearhead >= dshead);
						unsigned long long qgess = (dshead + 1ul == 0) ? 1ul : ddtoclearhead / (dshead + 1ul);
						qgess = std::max<unsigned long long>(1ul, qgess);
						assert(qgess >= 1ul && unsigned(qgess >> 32) == 0);
						if(qgess > 1) {
							unsigned qgess_[] = {unsigned(qgess)};
							std::fill(dsmultgess.begin(), dsmultgess.end(), 0u);
							mult_range(dsbegin, dsend, qgess_, qgess_+1, dsmultgess.begin(), dsmultgess.end());
							subtract_range(ddcurrenton, ddend, dsmultgess.begin(), dsmultgess.end());
							*(qbegin + (ddcurrenton - ddbegin)) += qgess;
						}else{
							subtract_range(ddcurrenton, ddend, dsbegin, dsbig+1);
							*(qbegin + (ddcurrenton - ddbegin)) += 1;
						}
					}else{
						ddcurrenton_rev ++;
						ddtoclear_rev ++;
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
			fixed_unsigned_int & operator /= (const fixed_unsigned_int & v);

			// this = remainder, returns the quotient
			fixed_unsigned_int divide(const fixed_unsigned_int & ds);
			fixed_unsigned_int mod(const fixed_unsigned_int & ds) const;

		private:
			std::array<unsigned, unit_num> mag_;
		};

		template <int num>
		int fixed_unsigned_int<num>::compare_to( const fixed_unsigned_int & v ) const
		{
			return __private__::compare_range(mag_.begin(), mag_.end(), v.mag_.begin(), v.mag_.end());
		}


		template <int num>
		bool fixed_unsigned_int<num>::is_zero() const
		{
			return std::none_of(mag_.begin(), mag_.end(), [](unsigned i)->bool{return i != 0;});
		}

		template <int num>
		int fixed_unsigned_int<num>::most_significant_nonzero_unit() const
		{
			auto i = __private__::find_most_significant_nonzero(mag_.begin(), mag_.end());
			return i == mag_.end() ? -1 : (i - mag_.begin());
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
			while(n != 0){
				unsigned nn = std::min(32u, n);
				__private__::left_shift_range(mag_.begin(), mag_.end(), nn);
				n -= nn;
			}
		}


		template <int num>
		void fixed_unsigned_int<num>::right_shift( unsigned n )
		{
			while(n != 0){
				unsigned nn = std::min(32u, n);
				__private__::right_shift_range(mag_.begin(), mag_.end(), nn);
				n -= nn;
			}
		}


		template <int num>
		void fixed_unsigned_int<num>::twos_complement()
		{
			__private__::twos_complement_range(mag_.begin(), mag_.end());
		}

		template <int num>
		void fixed_unsigned_int<num>::negate()
		{
			twos_complement();
			plus_one();
		}

		template <int num>
		void fixed_unsigned_int<num>::plus_one()
		{
			__private__::plus_one_range(mag_.begin(), mag_.end());
		}

		template <int num>
		void fixed_unsigned_int<num>::subtract_one()
		{
			__private__::subtract_one_range(mag_.begin(), mag_.end());
		}


		template <int num>
		fixed_unsigned_int<num> & fixed_unsigned_int<num>::operator+=( const fixed_unsigned_int & v )
		{
			__private__::plus_range(mag_.begin(), mag_.end(), v.mag_.begin(), v.mag_.end());
			return *this;
		}

		template <int num>
		fixed_unsigned_int<num> & fixed_unsigned_int<num>::operator-=( const fixed_unsigned_int & v )
		{
			__private__::subtract_range(mag_.begin(), mag_.end(), v.mag_.begin(), v.mag_.end());
			return *this;
		}

		template <int num>
		fixed_unsigned_int<num> & fixed_unsigned_int<num>::operator*=( const fixed_unsigned_int & v )
		{
			*this = (*this * v);
			return *this;
		}


		template <int num>
		fixed_unsigned_int<num> & fixed_unsigned_int<num>::operator/=( const fixed_unsigned_int & v )
		{
			*this = (*this / v);
			return *this;
		}

		template <int num>
		fixed_unsigned_int<num> fixed_unsigned_int<num>::divide( const fixed_unsigned_int & ds )
		{
			fixed_unsigned_int<num> q;
			__private__::divide_range(mag_.begin(), mag_.end(), ds.mag_.begin(), ds.mag_.end(), q.mag_.begin(), q.mag_.end());
			return q;
		}

		template <int num>
		fixed_unsigned_int<num> fixed_unsigned_int<num>::mod( const fixed_unsigned_int & ds ) const
		{
			fixed_unsigned_int<num> n(*this), q;
			__private__::divide_range(n.mag_.begin(), n.mag_.end(), ds.mag_.begin(), ds.mag_.end(), q.mag_.begin(), q.mag_.end());
			return n;
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
			std::array<unsigned, num1 + num2> p;
			std::fill(p.begin(), p.end(), 0);
			__private__::mult_range(xx.magnitude().begin(), xx.magnitude().end(), 
				yy.magnitude().begin(), yy.magnitude().end(), p.begin(), p.end());
			return fixed_unsigned_int<num1 + num2>(p);
		}

		// wrapped
		template <int num>
		fixed_unsigned_int<num> operator * (const fixed_unsigned_int<num>& xx, const fixed_unsigned_int<num>& yy)
		{
			std::array<unsigned, num> p;
			std::fill(p.begin(), p.end(), 0);
			__private__::mult_range(xx.magnitude().begin(), xx.magnitude().end(), 
				yy.magnitude().begin(), yy.magnitude().end(), p.begin(), p.end());
			return fixed_unsigned_int<num>(p);
		}

		template <int num>
		fixed_unsigned_int<num> operator / (const fixed_unsigned_int<num>& xx, const fixed_unsigned_int<num>& yy)
		{
			fixed_unsigned_int<num> xxx(xx);
			return xxx.divide(yy);
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
