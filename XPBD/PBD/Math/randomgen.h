#ifndef __RANDOM_H__
#define __RANDOM_H__

static const unsigned int SEED_X = 521288629u;
static const unsigned int SEED_Y = 362436069u;

typedef unsigned int uint;

//============================================================================
//	Classes
//
//----------------------------------------------------------------------------
//	RandomGen class definition
//!	Simple, fast, good random number generator.
//!	Perhaps the fastest of any generator that passes the Diehard tests.
//!
//!	@implementation
//!		Concatenation of following two 16-bit multiply-with-carry generators
//!		x(n)=a*x(n-1)+carry mod 2^16 and y(n)=b*y(n-1)+carry mod 2^16,
//!		number and carry packed within the same 32 bit integer.
//!		Algorithm recommended by Marsaglia<br/><br/>
//! @see
//!		http://paul.rutgers.edu/~rhoads/Code/code.html
class RandomGen {
public:
	inline RandomGen(int seed = 0);

	inline void  setSeed(int seed);

	inline int getInt() const;
	inline float getFloat() const;
	inline double getDouble() const;

private:
	mutable unsigned int seedX_m;
	mutable unsigned int seedY_m;
};


//============================================================================
//	Implementation

inline RandomGen::RandomGen(const int seed) {
	RandomGen::setSeed(seed);
}

inline void RandomGen::setSeed(const int seed) {
	seedX_m = (0 != seed) ? uint(seed) : SEED_X;
	seedY_m = (0 != seed) ? uint(seed) : SEED_Y;
}

inline int RandomGen::getInt() const {
	// Use any pair of non-equal numbers from this list for the two constants:
	// 18000 18030 18273 18513 18879 19074 19098 19164 19215 19584
	// 19599 19950 20088 20508 20544 20664 20814 20970 21153 21243
	// 21423 21723 21954 22125 22188 22293 22860 22938 22965 22974
	// 23109 23124 23163 23208 23508 23520 23553 23658 23865 24114
	// 24219 24660 24699 24864 24948 25023 25308 25443 26004 26088
	// 26154 26550 26679 26838 27183 27258 27753 27795 27810 27834
	// 27960 28320 28380 28689 28710 28794 28854 28959 28980 29013
	// 29379 29889 30135 30345 30459 30714 30903 30963 31059 31083

	seedX_m = 18000u * (seedX_m & 0xFFFFu) + (seedX_m >> 16);
	seedY_m = 30903u * (seedY_m & 0xFFFFu) + (seedY_m >> 16);

	return int((seedX_m << 16) + (seedY_m & 0xFFFFu));
}

inline float RandomGen::getFloat() const {
	return float((unsigned int)(getInt())) / 4294967296.0f;
}

inline double RandomGen::getDouble() const {
	return double((unsigned int)(getInt())) / 4294967296.0;
}

#endif//__RANDOM_H__
