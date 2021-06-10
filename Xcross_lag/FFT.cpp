// fft.cpp : Fast Fourier Transform, FFT
////////////////////////////////////////////////////////////////////////////////

#include <windows.h>    // APIENTRY, HANDLE, ...
#include <stdio.h>      // fprintf, stderr
#include <math.h>       // sin, cos, log10, fmod
#include <assert.h>     // assert

#include "fft.h"

// -----------------------------------------------------------------------------
// Magic Numbers
// -----------------------------------------------------------------------------

namespace {

static const unsigned  __FALSE  = FFT_FALSE;
static const unsigned  __TRUE   = FFT_TRUE;
static const unsigned  __DEFECT = FFT_DEFECT;

static const double    __PI     = FFT_PI;
static const double    __2PI    = 2.0 * __PI;
static const double    __R2D    = FFT_R2D;
static const double    __D2R    = FFT_D2R;

static const unsigned  __E__NOERR              = FFT__E__NOERR              ;
static const unsigned  __E__INVALID_ARGUMENT   = FFT__E__INVALID_ARGUMENT   ;
static const unsigned  __E__UNAVAILABLE_MEMORY = FFT__E__UNAVAILABLE_MEMORY ;
static const unsigned  __E__NOT_POWER_OF_TWO   = FFT__E__NOT_POWER_OF_TWO   ;
static const unsigned  __E__DEFECTIVE_HARMONIC = FFT__E__DEFECTIVE_HARMONIC ;

static const char* const  __gTableErrMsg[] =
{
    "0) No error",
    "1) Invalid argument",
    "2) Unavailable memory",
    "3) This number of samples is not power of two",
    "4) There is one or more defective hamonics",
};

} // namespace

// -----------------------------------------------------------------------------
// Macros
// -----------------------------------------------------------------------------

#define  __INVALID_ARGUMENT(x)  if ( x ) { return __E__INVALID_ARGUMENT; }
#define  __CHECK_POINTER(p)     __INVALID_ARGUMENT( 0 == p )

// -----------------------------------------------------------------------------
// fft_errmsg
// -----------------------------------------------------------------------------
/*
**  This returns a static constant error string corresponded its argument
*/
extern "C"
const char* __stdcall
fft_errmsg(
      unsigned       dwIndex        // [IN ] error index, result of the following functions
    )                               // [RET] error string corresponded its argument
{
    static const unsigned  __n = sizeof(__gTableErrMsg)/sizeof(__gTableErrMsg[0]);

    if ( __n <= dwIndex ) {

        return "Unknown error, this index is out-of-range!";
    }

    return __gTableErrMsg[dwIndex];
}
// -----------------------------------------------------------------------------
// is_power_of_two
// -----------------------------------------------------------------------------
/*
**  The following function tests whether its argument is a power of two for any
**  non-negative exponent k: x==2^k
*/
extern "C"
unsigned __stdcall
is_power_of_two(
      unsigned       x              // [IN ] number for checking, x
    )                               // [RET] 1) yes, 0) no
{
    if ( x < 2 ) {
        return __FALSE;
    }
    if ( x & (x-1) ) {              // thanks to 'byang' for this cute trick!
        return __FALSE;
    }
    return __TRUE;
}
// -----------------------------------------------------------------------------
// number_of_bits_needed
// -----------------------------------------------------------------------------
/*
**  The following function return the number of bits needed of its argument
*/
extern "C"
unsigned __stdcall
number_of_bits_needed(
      unsigned       dwPowerOfTwo   // [IN ] number, must be power of two
    )                               // [RET] number of bits needed of its argument
{
    unsigned  i = 0;

    if ( 0 == dwPowerOfTwo ) { goto _QUIT_; }

    for ( ; ; ++i ) {

        if ( dwPowerOfTwo & (1 << i) ) { break; }
    }

_QUIT_:
    return i;
}
// -----------------------------------------------------------------------------
// reverse_bits
// -----------------------------------------------------------------------------
/*
**  The following function reverses the bit order of index number
*/
/*
   Bit-reverse 位元反向
   假設資料共有8筆，第一遞迴的Butterfly operation組合是這樣的：a0跟a4、a2跟a6、a1跟a5、a3跟a7。
   0,4,2,6,1,5,3,7 這串數字您有沒有看到這隱藏著一個關係，若有學過密碼理論應該可以一眼就看出來。
   我們將它們轉成Binary code來表示：000, 100, 010, 110, 001, 101, 011, 111。因為資料共有8筆，
   任一個index以二進位表示是絕不會超過3的bit，若我們把每個Binary code反過來看會變成：
   000, 001, 010, 011, 100, 101, 110, 111。竟然會是0, 1, 2, 3, 4, 5, 6, 7的順序！
   所以，我們若將原來的0~N的資料順序，經過Bit-reverse的重新排列後，第一遞迴就可以滿足上述的規律，
   完成整個Iterative FFT了。
*/
extern "C"
unsigned __stdcall
reverse_bits(
      unsigned       dwIndex        // [IN ] index number
    , unsigned       dwNumBits      // [IN ] number of bits of index
    )                               // [RET] reversed bit order of the index number
{
    unsigned  i, rev;

    for ( i=rev=0; i < dwNumBits; ++i ) {

        rev = (rev << 1) | (dwIndex & 1);
        dwIndex >>= 1;
    }

    return rev;
}
// -----------------------------------------------------------------------------
// index_to_frequency
// -----------------------------------------------------------------------------
/*
**  The following function returns an "abstract frequency" of a given index into a
**  buffer with a given number of frequency samples.
**  Multiply return value by sampling rate to get frequency expressed in Hz.
*/
extern "C"
double __stdcall
index_to_frequency(
      unsigned       dwSamples      // [IN ] number of samples
    , unsigned       dwIndex        // [IN ] index
    )                               // [RET] abstract frequency
{
    if ( dwIndex >= dwSamples ) {
        return 0.0;
    }
    else if ( dwIndex <= dwSamples/2 ) {
        return (double)dwIndex / (double)dwSamples;
    }

    return -(double)(dwSamples-dwIndex) / (double)dwSamples;
}
// -----------------------------------------------------------------------------
// fft
// -----------------------------------------------------------------------------
// simple radix-2

template < typename T >
static
unsigned __stdcall
fft(
      unsigned  IsInverse   // [IN ] 0) forward, x) inverse transform
    , unsigned  dwSamples   // [IN ] number of samples, must be power of 2
    , const T*  lpSrcReal   // [IN ] source samples, real part
    , const T*  lpSrcImag   // [IN ] source samples, image part, could be 0
    , T*        lpTgtReal   // [OUT] target samples, real part
    , T*        lpTgtImag   // [OUT] target samples, image part
    )
{
    __INVALID_ARGUMENT(
       0 == lpSrcReal
    || 0 == lpTgtReal
    || 0 == lpTgtImag
    || lpSrcReal == lpSrcImag
    || lpTgtReal == lpTgtImag
    || lpSrcReal == lpTgtReal
    || lpSrcImag == lpTgtReal
    || lpSrcReal == lpTgtImag
    || lpSrcImag == lpTgtImag
    );

    if ( ! is_power_of_two(dwSamples) ) { return __E__NOT_POWER_OF_TWO; }

    double  angle_numerator = __2PI;    // 2.0 * __PI;

    if ( ! IsInverse ) {
        angle_numerator = -angle_numerator;
    }

    // number of bits needed to store indices
    //
    const unsigned  dwNumBits = number_of_bits_needed(dwSamples);

    /*
    **   Do simultaneous data copy and bit-reversal ordering into outputs...
    */
    unsigned  i, j;

    for ( i=0; i < dwSamples; ++i ) {

        j = reverse_bits( i, dwNumBits );
        lpTgtReal[j] = lpSrcReal[i];
        lpTgtImag[j] = (0 == lpSrcImag)? 0.0 : lpSrcImag[i];
    }

    /*
    **   Do the FFT itself...
    */
    double  tr;     // temp real
    double  ti;     // temp imaginary

    unsigned  k, n;

    unsigned  dwBlockEnd = 1;
    unsigned  dwBlockSize;
    double    ar[3], ai[3];

    for (dwBlockSize = 2; dwBlockSize <= dwSamples; dwBlockSize <<= 1) {

        const double  delta_angle = angle_numerator / (double)dwBlockSize;
        const double  sm2 = sin( -2 * delta_angle );
        const double  sm1 = sin(     -delta_angle );
        const double  cm2 = cos( -2 * delta_angle );
        const double  cm1 = cos(     -delta_angle );
        const double  w = 2 * cm1;

        for ( i=0; i < dwSamples; i += dwBlockSize ) {

            ar[2] = cm2;
            ar[1] = cm1;

            ai[2] = sm2;
            ai[1] = sm1;

            for ( j=i, n=0; n < dwBlockEnd; ++j, ++n ) {

                ar[0] = w*ar[1] - ar[2];
                ar[2] = ar[1];
                ar[1] = ar[0];

                ai[0] = w*ai[1] - ai[2];
                ai[2] = ai[1];
                ai[1] = ai[0];

                k = j + dwBlockEnd;
                tr = ar[0]*lpTgtReal[k] - ai[0]*lpTgtImag[k];
                ti = ar[0]*lpTgtImag[k] + ai[0]*lpTgtReal[k];

                lpTgtReal[k] = lpTgtReal[j] - tr;
                lpTgtImag[k] = lpTgtImag[j] - ti;

                lpTgtReal[j] += tr;
                lpTgtImag[j] += ti;
            }
        }

        dwBlockEnd = dwBlockSize;
    }

    /*
    **   Need to normalize if inverse transform...
    */
    if ( IsInverse ) {

        const double  denom = (double)dwSamples;

        for ( i=0; i < dwSamples; ++i ) {

            lpTgtReal[i] /= denom;
            lpTgtImag[i] /= denom;
        }
    }

    return __E__NOERR;
}
// -----------------------------------------------------------------------------
// fft_spectrum
// -----------------------------------------------------------------------------

template < typename T >
static
unsigned __stdcall
fft_spectrum(
      unsigned  dwSamples   // [IN ] number of samples, must be power of 2
    , const T*  lpSrcReal   // [IN ] source samples, real part
    , const T*  lpSrcImag   // [IN ] source samples, image part, could be NULL
    , T*        lpPowerSp   // [OUT] power outputs, reals only, dB
    , T*        lpPhaseSp   // [OUT] phase outputs, reals only, radians
    )
{
    unsigned  result =
    fft<T>(
      FFT_FORWARD           // [IN ] 0) forward, x) inverse transform
    , dwSamples             // [IN ] number of samples, must be power of 2
    , lpSrcReal             // [IN ] source samples, real part
    , lpSrcImag             // [IN ] source samples, image part, could be NULL
    , lpPowerSp             // [OUT] target samples, real part
    , lpPhaseSp             // [OUT] target samples, image part
    );

    if ( __E__NOERR == result ) {

        T  __real, __imag;

        const T* const  __lpPowerEnd = lpPowerSp + dwSamples;

        for ( ; lpPowerSp < __lpPowerEnd; ++lpPowerSp, ++lpPhaseSp ) {

            __real = *lpPowerSp;
            __imag = *lpPhaseSp;

            *lpPowerSp = __real*__real + __imag*__imag;
            *lpPhaseSp = atan2( __real, __imag );
        }
    }

//_QUIT_:
    return result;
}
// -----------------------------------------------------------------------------
// fft_float
// -----------------------------------------------------------------------------

extern "C"
unsigned __stdcall
fft_float(
      unsigned       IsInverse      // [IN ] 0) forward, x) inverse transform
    , unsigned       dwSamples      // [IN ] number of samples, must be power of 2
    , const float*   lpSrcReal      // [IN ] source samples, reals
    , const float*   lpSrcImag      // [IN ] source samples, imaginaries, could be NULL
    , float*         lpTgtReal      // [OUT] target outputs, reals
    , float*         lpTgtImag      // [OUT] target outputs, imaginaries
    )                               // [RET] 0) no error, x) error index
{
    return
    fft<float>(
      IsInverse
    , dwSamples
    , lpSrcReal
    , lpSrcImag
    , lpTgtReal
    , lpTgtImag
    );
}
// -----------------------------------------------------------------------------
// fft_double
// -----------------------------------------------------------------------------

extern "C"
unsigned __stdcall
fft_double(
      unsigned       IsInverse      // [IN ] 0) forward, x) inverse transform
    , unsigned       dwSamples      // [IN ] number of samples, must be power of 2
    , const double*  lpSrcReal      // [IN ] source samples, reals
    , const double*  lpSrcImag      // [IN ] source samples, imaginaries, could be NULL
    , double*        lpTgtReal      // [OUT] target outputs, reals
    , double*        lpTgtImag      // [OUT] target outputs, imaginaries
    )                               // [RET] 0) no error, x) error index
{
    return
    fft<double>(
      IsInverse
    , dwSamples
    , lpSrcReal
    , lpSrcImag
    , lpTgtReal
    , lpTgtImag
    );
}
// -----------------------------------------------------------------------------
// fft_spectrum_float
// -----------------------------------------------------------------------------

extern "C"
unsigned __stdcall
fft_spectrum_float(
      unsigned       dwSamples      // [IN ] number of samples, must be power of 2
    , const float*   lpSrcReal      // [IN ] source samples, real part
    , const float*   lpSrcImag      // [IN ] source samples, image part, could be NULL
    , float*         lpPowerSp      // [OUT] power outputs, reals only
    , float*         lpPhaseSp      // [OUT] phase outputs, reals only, radians
    )                               // [RET] 0) no error, x) error index
{
    return
    fft_spectrum<float>(
      dwSamples
    , lpSrcReal
    , lpSrcImag
    , lpPowerSp
    , lpPhaseSp
    );
}
// -----------------------------------------------------------------------------
// fft_spectrum_double
// -----------------------------------------------------------------------------

extern "C"
unsigned __stdcall
fft_spectrum_double(
      unsigned       dwSamples      // [IN ] number of samples, must be power of 2
    , const double*  lpSrcReal      // [IN ] source samples, real part
    , const double*  lpSrcImag      // [IN ] source samples, image part, could be NULL
    , double*        lpPowerSp      // [OUT] power outputs, reals only
    , double*        lpPhaseSp      // [OUT] phase outputs, reals only, radians
    )                               // [RET] 0) no error, x) error index
{
    return
    fft_spectrum<double>(
      dwSamples
    , lpSrcReal
    , lpSrcImag
    , lpPowerSp
    , lpPhaseSp
    );
}
// -----------------------------------------------------------------------------
// __round
// -----------------------------------------------------------------------------
/*
static
unsigned __fastcall
__round(
      unsigned  x
    )
{
    return unsigned( 0.5f + double(x) );
}
*/
// -----------------------------------------------------------------------------
// __max_1st
// -----------------------------------------------------------------------------
// to get the first maximum value between [p1, p2)
/*
template < typename T >
static
T __fastcall
__max_1st(
      const T*        p1
    , const T* const  p2
    )
{
    T  __max = *p1;

    if ( p2 < p1 ) { goto _QUIT_; }

    for ( ; p1<p2; ++p1 ) {

        if ( __max < *p1 ) {

            __max = *p1;
        }
    }

_QUIT_:
    return __max;
}*/
// -----------------------------------------------------------------------------
// __span_max
// -----------------------------------------------------------------------------
//
// this will find the index number which has maximum value of its fft power spectrum
// within limited range specified by its center and span numbers.
//
// [NOTE] round by half number of samples
//
static
unsigned __fastcall
__span_max(
      unsigned       dwHalf         // [IN ] half number of samples, must be power of 2
    , const double*  lpPowerSp      // [IN ] fft power spectrum
    , unsigned       dwCenter       // [IN ] center number, index
    , unsigned       dwSpan         // [IN ] span on each side
    )                               // [RET] index number has maximum value within range
{
    // round by half number of samples
    //
    dwCenter %= dwHalf;

    // center part
    //
    double  max = lpPowerSp[dwCenter];

    // span part
    //
    unsigned  r = dwCenter;
    unsigned  i = 1;
    unsigned  v;
    for ( ; i<=dwSpan; ++i ) {

        // positive way, right side
        //
        v = (dwCenter +i) %dwHalf;
        if ( max < lpPowerSp[v] ) {
            max = lpPowerSp[v];
            r = v;
        }
        // negative way, left side
        //
        v = (dwCenter +dwHalf -i) %dwHalf;
        if ( max < lpPowerSp[v] ) {
            max = lpPowerSp[v];
            r = v;
        }
    }

    assert( r < dwHalf );

    return r;
}
// -----------------------------------------------------------------------------
// __span_sum
// -----------------------------------------------------------------------------
//
// this will sum the values of its fft power spectrum within limited range
// specified by its center and span numbers.
//
// [NOTE] round by half number of samples
//
static
double __fastcall
__span_sum(
      unsigned       dwHalf         // [IN ] half number of samples, must be power of 2
    , const double*  lpSpectrum     // [IN ] fft power spectrum
    , unsigned       dwCenter       // [IN ] center number, index
    , unsigned       dwSpan         // [IN ] span on each side
    )                               // [RET] value summed within range
{
    dwCenter %= dwHalf;

    // center part
    //
    double  sum = lpSpectrum[dwCenter];

    // span part
    //
    unsigned  j, v;
    for ( j=1; j<=dwSpan; ++j ) {

        // positive way, right side
        //
        v = (dwCenter +j) %dwHalf;
        sum += lpSpectrum[v];

        // negative way, left side
        //
        v = (dwCenter +dwHalf -j) %dwHalf;
        sum += lpSpectrum[v];
    }

    return sum;
}
// -----------------------------------------------------------------------------
// __is_banded
// -----------------------------------------------------------------------------
//
// this will return nonzero value when its target index is in the range specified by
// its center and span numbers.
//
// [P.S.] range = [center-span, center+span]
//
// [NOTE] round by half number of samples
//
static
unsigned __fastcall
__is_banded(
      unsigned      dwHalf          // [IN ] half number of samples, must be power of 2
    , unsigned      dwTarget        // [IN ] target index
    , unsigned      dwCenter        // [IN ] center index
    , unsigned      dwSpan          // [IN ] span on each side
    )                               // [RET] in this range? 0) no x) yes
{
    // to ensure those rounded by half number of samples
    //
    dwTarget %= dwHalf;
    dwCenter %= dwHalf;

    // to check at center part
    //
    if ( dwTarget == dwCenter ) { return __TRUE; }

    // to check at span part
    //
    unsigned  j, v;
    for ( j=1; j<=dwSpan; ++j ) {

        // positive way, right side
        //
        v = (dwCenter +j) %dwHalf;
        if ( dwTarget == v ) { return __TRUE; }

        // negative way, left side
        //
        v = (dwCenter +dwHalf -j) %dwHalf;
        if ( dwTarget == v ) { return __TRUE; }
    }

    return __FALSE;
}
// -----------------------------------------------------------------------------
// __is_dc
// -----------------------------------------------------------------------------
//
// this will return nonzero value when the range specified by its center and span
// numbers contain some dc.
//
// [P.S.] range = [center-span, center+span]
//
// [NOTE] round by half number of samples
//
static
unsigned __fastcall
__is_dc(
      unsigned      dwHalf          // [IN ] half number of samples, must be power of 2
    , unsigned      dwDc            // [IN ] dc index
    , unsigned      dwCenter        // [IN ] center index
    , unsigned      dwSpan          // [IN ] span on each side
    )                               // [RET] in this range? 0) no x) yes
{
    unsigned  j;
    for ( j=0; j<=dwDc; ++j ) {

        if ((
            __is_banded(
                  dwHalf            // [IN ] half number of samples, must be power of 2
                , j                 // [IN ] target index
                , dwCenter          // [IN ] center index
                , dwSpan            // [IN ] span on each side
                )
           ))
        {
            return __TRUE;
        }
    }

    return __FALSE;
}
// -----------------------------------------------------------------------------
// __max_spurious
// -----------------------------------------------------------------------------
/*
**  This will search the maximum spurious bin except the signal with span on each side
*/
static
unsigned __fastcall
__max_spurious(
      unsigned       dwHalf         // [IN ] half number of samples, must be power of 2
    , const double*  lpSpectrum     // [IN ] fft power spectrum
    , unsigned       dwIndexDc      // [IN ] parameter used to avoid the dc stuff
    , unsigned       dwSignal       // [IN ] signal index
    , unsigned       dwSpan         // [IN ] span on each side
    )                               // [RET] maximum spurious bin
{
    // to search max spur between dc and main bin, left part
    //
    // %%% Range1 = notdc:1:harpos(1) - span;
    // %%% MaxSpur1 = max(spectP(Range1));
    // %%% MaxSpur1_pos = find(spectP(Range1)==MaxSpur1) +notdc -1;
    // %%% MaxSpur1_pos = MaxSpur1_pos(1);   % first one
    //
    unsigned  bin = __DEFECT;
    double    max = 0.0;
    unsigned  i;
    if ( dwSignal > dwSpan ) {
        const unsigned  __left = dwSignal -dwSpan;
        for ( i=dwIndexDc; i<__left; ++i ) {
            if ( max < lpSpectrum[i] ) {
                max = lpSpectrum[i];
                bin = i;
            }
        }
    }

    // to search max spru between main bin and the rest bin, right part
    //
    // %%% Range2 = harpos(1)+span:1:numpt/2;
    // %%% MaxSpur2 = max(spectP(Range2));
    // %%% MaxSpur2_pos = find(spectP(Range2) == MaxSpur2) +harpos(1)+span -1 ;
    //
    for ( i=dwSignal +dwSpan; i<dwHalf; ++i ) {
        if ( max < lpSpectrum[i] ) {
            max = lpSpectrum[i];
            bin = i;
        }
    }

    // %%% MaxSpur = max([MaxSpur1, MaxSpur2]);
    // %%% MaxSpur_pos = find(spectP(1:numpt/2) == MaxSpur);

    assert( __DEFECT != bin );

    return bin;
}
// -----------------------------------------------------------------------------
// __harmonics
// -----------------------------------------------------------------------------
/*
**  This finds out the harmonic frequencies (indexes) in the fft power spectrum.
**
**  [NOTE] Defective Harmonic <-- FFT_DEFECT (0xFFFFFFFF);
**         For this procedure to work, ensure the folded back high order harmonics
**         do not overlap with dc or signal or lower order harmonics.
**
*/
static
unsigned __fastcall
__harmonic(
      unsigned       dwHalf         // [IN ] half number of samples, must be power of 2
    , const double*  lpSpectrum     // [IN ] fft power spectrum
    , unsigned       dwIndexDc      // [IN ] parameter used to avoid the dc stuff
    , unsigned       dwSearch       // [IN ] approximate search span for harmonics on each side
    , unsigned       dwMainSpan     // [IN ] span of the main signal frequency on each side
    , unsigned       dwHarmSpan     // [IN ] span of the harmonic frequency on each side
    , unsigned       dwCurrent      // [IN ] current harmonic frequencies for checking
    , unsigned       dwCount        // [IN ] number of harmonic frequencies
    , unsigned*      lpHarmonics    // [IN ] harmonic frequencies
    )                               // [RET] harmonic frequency, searched
{
    unsigned  __index = __DEFECT;

    assert( dwSearch >= dwHarmSpan );

    for ( ; dwSearch>=dwHarmSpan; dwSearch>>=1 ) {

        __index =
        __span_max(
              dwHalf                // [IN ] number of samples, must be power of 2
            , lpSpectrum            // [IN ] fft power spectrum
            , dwCurrent             // [IN ] center number, index
            , dwSearch              // [IN ] approximate search span for harmonics on each side
            );

        // overlap with dc
        //
        if ((
            __is_dc(
                  dwHalf            // [IN ] half number of samples, must be power of 2
                , dwIndexDc         // [IN ] dc index
                , __index           // [IN ] center index
                , dwHarmSpan        // [IN ] span on each side
                )                   // [RET] in this range? 0) no x) yes
           ))
        {
            __index = __DEFECT;     // defective
            continue;               // search again
        }

        // overlap with signal
        //
        if ((
            __is_banded(
                  dwHalf            // [IN ] half number of samples, must be power of 2
                , __index           // [IN ] target index
                , lpHarmonics[0]    // [IN ] center index
                , dwMainSpan        // [IN ] span on each side
                )                   // [RET] in this range? 0) no x) yes
           ))
        {
            __index = __DEFECT;     // defective
            continue;               // search again
        }

        // overlap with lower order harmonics
        //
        unsigned  i;
        for ( i=1; i<dwCount; ++i ) {
            if ((
                __DEFECT != lpHarmonics[i] &&
                __is_banded(
                      dwHalf            // [IN ] half number of samples, must be power of 2
                    , __index           // [IN ] target index
                    , lpHarmonics[i]    // [IN ] center index
                    , dwHarmSpan        // [IN ] span on each side
                    )                   // [RET] in this range? 0) no x) yes
               ))
            {
                __index = __DEFECT;     // defective
                break;                  // search again, break loop i
            }
        }

        if ( __DEFECT != __index ) {

            break;                      // indefective, break loop search
        }
    }

    return __index;
}
static
unsigned __fastcall
__harmonics(
      unsigned       dwHalf         // [IN ] half number of samples, must be power of 2
    , const double*  lpSpectrum     // [IN ] fft power spectrum
    , unsigned       dwIndexDc      // [IN ] parameter used to avoid the dc stuff
    , unsigned       dwSearch       // [IN ] approximate search span for harmonics on each side
    , unsigned       dwMainSpan     // [IN ] span of the main signal frequency on each side
    , unsigned       dwHarmSpan     // [IN ] span of the harmonic frequency on each side
    , unsigned       dwCount        // [IN ] number of harmonic frequencies, required
    , unsigned*      lpHarmonics    // [OUT] harmonic frequency outputs
    )                               // [RET] number of harmonic frequencies, actual
{
    // find the signal bin number, DC = dwIndexDc
    //
    // %%% maxdB = max(Dout_dB(notdc:numpt/2));
    // %%% fin = find(Dout_dB(1:numpt/2) == maxdB);
    //
    unsigned  fin = dwIndexDc;
    double    max = lpSpectrum[dwIndexDc];
    unsigned  i;
    for ( i = dwIndexDc; i<dwHalf; ++i ) {
        if ( lpSpectrum[i] > max ) {
            fin = i;
            max = lpSpectrum[i];
        }
    }

    // The 1st element in the array represents the signal, the next element represents
    // the 2nd harmonic, etc.
    //
    lpHarmonics[0] = fin;

    // %%% span_center(har_num) = round(rem(har_num*(fin-1)+1, numpt));
    // %%% span_low = span_center(har_num)-spanh;
    // %%% span_top = span_center(har_num)+spanh;
    // %%% har_peak = max(spectP(span_low:span_top));
    // %%% har_bin  = find(spectP(span_low:span_top)==har_peak,1,'first') + span_low -1;
    //
    unsigned  cnt = 1;
    for ( i=1; i<dwCount; ++i ) {

        fin = (fin + lpHarmonics[0] +1) % dwHalf;

        // for this procedure to work, ensure the folded back high order harmonics do not
        // overlap with dc or signal or lower order harmonics
        //
        lpHarmonics[i] =
        __harmonic(
              dwHalf        // [IN ] half number of samples, must be power of 2
            , lpSpectrum    // [IN ] fft power spectrum
            , dwIndexDc     // [IN ] parameter used to avoid the dc stuff
            , dwSearch      // [IN ] approximate search span for harmonics on each side
            , dwMainSpan    // [IN ] span of the main signal frequency on each side
            , dwHarmSpan    // [IN ] span of the harmonic frequency on each side
            , fin           // [IN ] current harmonic frequencies for checking
            , i             // [IN ] number of harmonic frequencies
            , lpHarmonics   // [IN ] harmonic frequencies
            );              // [RET] harmonic frequency, searched

        // indefective
        //
        if ( __DEFECT != lpHarmonics[i] ) {

            ++cnt;
        }
    }

    return cnt;
}
// -----------------------------------------------------------------------------
// hanning
// -----------------------------------------------------------------------------
/*
**  This will do one array multiplication with the input array and the N-point
**  symmetric Hanning window in a column array.
*/
extern "C"
unsigned __stdcall
hanning(
      unsigned       dwSamples      // [IN ] number of samples, must be power of 2
    , const double*  lpSources      // [IN ] source samples, reals, could be 0
    , double*        lpTargets      // [OUT] target results, reals
    )                               // [RET] actual processing count
{
    __CHECK_POINTER( lpTargets );

    // %%% w = .5*(1 - cos(2*pi*(1:n)'/(n+1)));
    //
    const bool  __is_src_existed = lpSources;
    unsigned    i;
    for ( i=1; i<=dwSamples; ++i, ++lpTargets ) {

        const double  __w = 0.5*( 1.0 - ::cos( (__2PI*i)/(1+dwSamples) ) );

        if ( __is_src_existed ) {
            *lpTargets = __w * (*lpSources);  ++lpSources;
        }
        else {
            *lpTargets = __w;
        }
    }

    return i-1;
}
// -----------------------------------------------------------------------------
// hamming
// -----------------------------------------------------------------------------
/*
**  This will do a array multiplication with the input vector and the N-point
**  symmetric Hamming window in a column vector.
*/
extern "C"
unsigned __stdcall
hamming(
      unsigned       dwSamples      // [IN ] number of samples
    , const double*  lpSources      // [IN ] source samples, reals, could be 0
    , double*        lpTargets      // [OUT] target results, reals
    )                               // [RET] actual processing count
{
    __CHECK_POINTER( lpTargets );

    // %%% w = .54 - .46*cos(2*pi*(0:n-1)'/(n-1));
    //
    const bool  __is_src_existed = lpSources;
    unsigned    i;
    for ( i=0; i<dwSamples; ++i, ++lpTargets ) {

        const double  __w = 0.54 - 0.46*(::cos((__2PI*i)/(dwSamples-1)));

        if ( __is_src_existed ) {
            *lpTargets = __w * (*lpSources);  ++lpSources;
        }
        else {
            *lpTargets = __w;
        }
    }

    return i;
}
// -----------------------------------------------------------------------------
// dynamic_performance
// -----------------------------------------------------------------------------
/*
**  This is used to calculate SNR, SINAD, THD, ENOB and SFDR values.
**
**  SINAD  : Signal-to-Noise and Distortion Ratio
**  SNR    : Signal-to-Noise Ratio
**  THD    : Total Harmonic Distortion
**  ENOB   : Effective Number of Bits
**  SFDR   : Spurious-Free Dynamic Range
**  MSB    : Maximum Spurious Bin
**
**  There are some suggested values as follow:
**      dwMainSpan <-- max(5, dwSamples/200)
**      dwHarmSpan <-- 1
**      dwSearch   <-- max(dwHarmSpan, dwHalf >> 4)
**      range(center, span) = [center-span, center+span]
**
**  [NOTE] Defective Harmonic <-- FFT_DEFECT (0xFFFFFFFF);
**         For this procedure to work, ensure the folded back high order harmonics
**         do not overlap with dc or signal or lower order harmonics.
**
**  <<< CHECKING >>>
**  unsigned  __error_index = dynamic_performance(...);
**  if ( __error_index ) { puts(fft_errmsg(__error_index)); return false; }
*/
extern "C"
unsigned __stdcall
dynamic_performance(
      unsigned       dwHalf         // [IN ] half number of samples, must be power of 2
    , const double*  lpSpectrum     // [IN ] array, fft power spectrum
    , unsigned       dwIndexDc      // [IN ] parameter used to avoid the dc stuff
    , unsigned       dwSearch       // [IN ] approximate search span for harmonics on each side
    , unsigned       dwMainSpan     // [IN ] span of the main signal frequency on each side
    , unsigned       dwHarmSpan     // [IN ] span of the harmonic frequency on each side
    , unsigned       dwHarmCount    // [IN ] harmonic count
    , unsigned*      lpHarmIndexes  // [OUT] array, harmonic frequencies (indexes)
    , double*        lpHarmPowers   // [OUT] array, harmonic powers
    , double*        lpSINAD        // [OUT] pointer, signal-to-noise and distortion ratio
    , double*        lpSNR          // [OUT] pointer, signal-to-noise ratio
    , double*        lpTHD          // [OUT] pointer, total harmonic distortion
    , double*        lpENOB         // [OUT] pointer, effective number of bits
    , double*        lpSFDR         // [OUT] pointer, spurious-free dynamic range
    , unsigned*      lpMsbIndex     // [OUT] pointer, maximum spurious frequency (index)
    , double*        lpMsbPower     // [OUT] pointer, maximum spurious power
    )                               // [RET] 0) no error, x) error index
{
    __INVALID_ARGUMENT(
       0 == lpSpectrum
    || 0 == lpHarmIndexes
    || 0 == lpHarmPowers
    || 0 == dwHarmCount
    || 0 == dwHalf
    || dwHalf <  dwHarmCount
    || dwHalf <= dwIndexDc
    || lpSpectrum == lpHarmPowers
    );

    if ( ! is_power_of_two(dwHalf) ) { return __E__NOT_POWER_OF_TWO; }

    if ( dwSearch < dwHarmSpan ) { dwSearch = dwHarmSpan; }

    // start here
    //
    unsigned  result = __E__NOERR;

    // to search the harmonic bins
    //
    if ((
        dwHarmCount !=
        __harmonics(
              dwHalf         // [IN ] number of samples, must be power of 2
            , lpSpectrum     // [IN ] fft power spectrum
            , dwIndexDc      // [IN ] parameter used to avoid the dc stuff
            , dwSearch       // [IN ] approximate search span for harmonics on each side
            , dwMainSpan     // [IN ] span of the main signal frequency on each side
            , dwHarmSpan     // [IN ] span of the harmonic frequency on each side
            , dwHarmCount    // [IN ] number of harmonic frequencies, required
            , lpHarmIndexes  // [OUT] harmonic frequency outputs
            )                // [RET] number of harmonic frequencies, actual
       ))
    {
        // still running when this error happen
        //
        result = __E__DEFECTIVE_HARMONIC;
    }

    // to calculate power of signal and harmonics
    //
    // %%% Ph = [Ph sum(spectP(har_bin-1:har_bin+1))];
    //
    unsigned  i;
    for ( i=0; i<dwHarmCount; ++i ) {

        if ( __DEFECT == lpHarmIndexes[i] ) {

            lpHarmPowers[i] = 0;
        }
        else {

            lpHarmPowers[i] =
                __span_sum(
                  dwHalf            // [IN ] half number of samples, must be power of 2
                , lpSpectrum        // [IN ] fft power spectrum
                , lpHarmIndexes[i]  // [IN ] center number, index
                , dwHarmSpan        // [IN ] span on each side
                );

            if ( 0 == lpHarmPowers[i] ) {

                lpHarmIndexes[i] = __DEFECT;
            }
        }
        assert( 0 <= lpHarmPowers[i] );
    }

    // to extract overall signal power
    //
    // %%% Ps=sum(spectP(fin-span:fin+span));
    //
    double  Ps;
    Ps = __span_sum(
              dwHalf            // [IN ] half number of samples, must be power of 2
            , lpSpectrum        // [IN ] fft power spectrum
            , lpHarmIndexes[0]  // [IN ] center number, index
            , dwMainSpan        // [IN ] span on each side
            );
    assert( 0 <= Ps );

    // to find dc offset power
    //
    // %%% Pdc = sum(spectP(1:span));
    // %%% Pdc = sum(spectP(1:fin-span)) ???
    //
    double  Pdc;
    Pdc = 0;
    for ( i=0; i<dwMainSpan; ++i ) {
        Pdc += lpSpectrum[i];
    }
    assert( 0 <= Pdc );

    // to determine the total distortion power
    //
    // %%% Pd = sum(Ph(2:harmonic_number));
    //
    double  Pd;
    Pd = 0;
    for ( i=1; i<dwHarmCount; ++i ) {
        Pd += lpHarmPowers[i];
    }
    assert( 0 <= Pd );

    // to determine the noise power
    //
    // %%% Pn = sum(spectP(1:numpt/2))-Pdc-Ps-Pd;
    //
    double  Pn;
    Pn = -Pdc -Ps -Pd;
    for ( i=0; i<dwHalf; ++i ) {
        Pn += lpSpectrum[i];
    }

    // %%% MaxSpur_Power = max(Ph(2:harmonic_number));               <-- X
    // %%% MaxSpur_Power = sum(spectP(MaxSpur_pos-1:MaxSpur_pos+1)); <-- O
    //
    double  Pms;
    /*
    const double  Pms = __max_1st(lpHarmPowers+1, lpHarmPowers+dwHarmCount);
    */
    i = __max_spurious(
              dwHalf            // [IN ] half number of samples, must be power of 2
            , lpSpectrum        // [IN ] fft power spectrum
            , dwIndexDc         // [IN ] parameter used to avoid the dc stuff
            , lpHarmIndexes[0]  // [IN ] signal index
            , dwMainSpan        // [IN ] span on each side
            );                  // [RET] maximum spurious bin
    Pms = __span_sum(
              dwHalf            // [IN ] half number of samples, must be power of 2
            , lpSpectrum        // [IN ] fft power spectrum
            , i                 // [IN ] center number, index
            , dwHarmSpan        // [IN ] span on each side
            );
    assert( 0 <= Pms );

    if ( lpMsbIndex ) {
        *lpMsbIndex = i;
    }

    if ( lpMsbPower ) {
        *lpMsbPower = Pms;
    }


    // to calculate SNR, SINAD, THD, ENOB and SFDR values
    //
    // %%% SINAD = 10*log10(Ps/(Pn+Pd));
    // %%% SNR   = 10*log10(Ps/Pn);
    // %%% THD   = 10*log10(Pd/Ph(1));
    // %%% ENOB  = (SINAD-1.76)/6.02;
    // %%% SFDR  = 10*log10(Ph(1)/MaxSpur_Power);
    // %%% SFDR  = 10*log10(Ph(1)/max(Ph(2:10)))  ???
    //
    if ( 0 == Pn +Pd ) { lpSINAD = 0; }
    if ( 0 == Pn     ) { lpSNR   = 0; }
    if ( 0 == Pms    ) { lpSFDR  = 0; }
    if ( lpSINAD ) {
        *lpSINAD = 10.0 * ::log10(Ps/(Pn+Pd));
        if ( lpENOB ) {
            *lpENOB = (*lpSINAD-1.76)/6.02;
        }
    }
    if ( lpSNR ) {
        *lpSNR  = 10.0 * ::log10(Ps/Pn);
    }
    if ( lpTHD ) {
        *lpTHD  = 10.0 * ::log10(Pd/lpHarmPowers[0]);
    }
    if ( lpSFDR ) {
        *lpSFDR = 10.0 * ::log10(lpHarmPowers[0]/Pms);
    }

    return result;
}
// -----------------------------------------------------------------------------
// iba2cba
// -----------------------------------------------------------------------------
/*
**  iba : index-based array   <-- [ c1d1, c2d1, c3d1, ..., c1d2, c2d2, c3d2, ... ]
**  cba : channel-based array <-- [ c1d1, c1d2, ..., c2d1, c2d2, ..., c3d1, c3d2, ... ]
**
**  [NOTE] length(array) = dwChannels * dwSamples;
*/
extern "C"
void __stdcall
iba2cba(
      unsigned       dwChannels     // [IN ] number of channels
    , unsigned       dwSamples      // [IN ] number of samples for each channels
    , const double*  lpIba          // [IN ] source index-based array
    , double*        lpCba          // [OUT] target channel-based array
    )
{
    for ( unsigned c=0, cc=0; c<dwChannels; ++c, cc+=dwSamples ) {

        for ( unsigned r=0; r<dwSamples; ++r ) {

            lpCba[ cc + r ] = lpIba[ r*dwChannels + c ];
        }
    }
}
// -----------------------------------------------------------------------------
// cba2iba
// -----------------------------------------------------------------------------

extern "C"
void __stdcall
cba2iba(
      unsigned       dwChannels     // [IN ] number of channels
    , unsigned       dwSamples      // [IN ] number of samples for each channels
    , const double*  lpCba          // [IN ] source channel-based array
    , double*        lpIba          // [OUT] target index-based array
    )
{
    for ( unsigned c=0, cc=0; c<dwChannels; ++c, cc+=dwSamples ) {

        for ( unsigned r=0; r<dwSamples; ++r ) {

            lpIba[ r*dwChannels + c ] = lpCba[ cc + r ];
        }
    }
}
// -----------------------------------------------------------------------------
// log10_x
// -----------------------------------------------------------------------------
/*
**  <<< CHECKING >>>
**  double*  __ending = log10_xxx(dwLen, ..., lpTgt, ...);
**  if ( __ending != lpTgt +dwLen ) { return false; }
*/
extern "C"
double __stdcall
log10_single(
      double  x                     // [IN ] source number <-- x
    )                               // [RET] log10 value   <-- log10(x)
{
    return ::log10(x);
}

extern "C"
double* __stdcall
log10_range(
      unsigned long   dwLen         // [IN ] array length      <-- n
    , const double*   lpSrc         // [IN ] source array      <-- xx
    ,       double*   lpTgt         // [OUT] target array, log <-- log10(xx)
    )                               // [RET] ending pointer    <-- lpTgt +n
{
    const double* const  __lpEnd = lpSrc +dwLen;

    if ((
           0 == lpSrc
        || 0 == lpTgt
       ))
    {
        goto _QUIT_;
    }

    for ( ; lpSrc < __lpEnd; ++lpSrc, ++lpTgt ) {

        *lpTgt = ::log10(*lpSrc);
    }

_QUIT_:
    return lpTgt;
}

extern "C"
double* __stdcall
log10_complex(
      unsigned long   dwLen         // [IN ] array length              <-- n
    , const double*   lpSrcReal     // [IN ] source array, real part   <-- a
    , const double*   lpSrcImag     // [IN ] source array, image part  <-- b
    ,       double*   lpTgt         // [OUT] target array, log complex <-- 0.5*log10(a^2+b^2)
    )                               // [RET] ending pointer            <-- lpTgt +n
{
    return
    ::log10_power(
          dwLen         // [IN ] array length              <-- n
        , lpSrcReal     // [IN ] source reals              <-- a
        , lpSrcImag     // [IN ] source images             <-- b
        , 0             // [OUT] target powers, could be 0 <-- a^2+b^2
        , lpTgt         // [OUT] target powers, db value   <-- d+c*log10(a^2+b^2)
        , 0.5           // [OPT] constant multiplier       <-- c
        , 0.0           // [OPT] constant bias             <-- d
        );              // [RET] ending pointer            <-- lpTgtPd +n
}

extern "C"
double* __stdcall
log10_db(
      unsigned long   dwLen         // [IN ] array length           <-- n
    , const double*   lpSrc         // [IN ] source array           <-- x
    ,       double*   lpTgt         // [OUT] target array, db value <-- d+c*log10(x)
    , const double    cdMult        // [OPT] constant multiplier    <-- c
    , const double    cdBias        // [OPT] constant bias          <-- d
    )                               // [RET] ending pointer         <-- lpTgt +n
{
    const double* const  __lpEnd = lpTgt +dwLen;

    if ((
           0 == lpSrc
        || 0 == lpTgt
       ))
    {
        goto _QUIT_;
    }

    for ( ; lpTgt < __lpEnd; ++lpTgt, ++lpSrc ) {

        *lpTgt = cdBias + ( cdMult * ::log10(*lpSrc) );
    }

_QUIT_:
    return lpTgt;
}

extern "C"
double* __stdcall
log10_power(
      unsigned long   dwLen         // [IN ] array length              <-- n
    , const double*   lpSrcReal     // [IN ] source reals              <-- a
    , const double*   lpSrcImag     // [IN ] source images             <-- b
    ,       double*   lpTgtPw       // [OUT] target powers, could be 0 <-- a^2+b^2
    ,       double*   lpTgtPd       // [OUT] target powers, db value   <-- d+c*log10(a^2+b^2)
    , const double    cdMult        // [OPT] constant multiplier       <-- c
    , const double    cdBias        // [OPT] constant bias             <-- d
    )                               // [RET] ending pointer            <-- lpTgtPd +n
{
    const double* const  __lpEnd = lpTgtPd +dwLen;

    if ((
           0 == lpSrcReal
        || 0 == lpSrcImag
        || 0 == lpTgtPd
       ))
    {
        goto _QUIT_;
    }

    if ( 0 == lpTgtPw ) { lpTgtPw = lpTgtPd; }

    for ( ; lpTgtPd < __lpEnd; ++lpTgtPd, ++lpTgtPw, ++lpSrcReal, ++lpSrcImag ) {

        *lpTgtPw = (*lpSrcReal)*(*lpSrcReal) + (*lpSrcImag)*(*lpSrcImag);
        *lpTgtPd = cdBias + ( cdMult * ::log10(*lpTgtPw) );
    }

_QUIT_:
    return lpTgtPd;
}

extern "C"
double* __stdcall
log10_normalize(
      unsigned long   dwLen         // [IN ] array length                 <-- n
    , const double*   lpSrcReal     // [IN ] source reals                 <-- a
    , const double*   lpSrcImag     // [IN ] source images                <-- b
    ,       double*   lpTgtPw       // [OUT] target powers, could be 0    <-- a^2+b^2
    ,       double*   lpTgtPd       // [OUT] target powers, normalize db  <-- d-m+c*log10(a^2+b^2)
    , const double    cdMult        // [OPT] constant multiplier          <-- c
    , const double    cdBias        // [OPT] constant bias                <-- d
    ,       double*   lpMaxPd       // [OPT] maximum db value, could be 0 <-- m
    )                               // [RET] ending pointer               <-- lpTgtPd +n
{
    double* const  __lpBegin = lpTgtPd;
    double* const  __lpEnd   = __lpBegin + dwLen;

    double  __MaxPd = 0.0;

    if ((
           0 == lpSrcReal
        || 0 == lpSrcImag
        || 0 == lpTgtPd
       ))
    {
        goto _QUIT_;
    }

    if ( 0 == lpTgtPw ) { lpTgtPw = lpTgtPd; }


    for ( ; lpTgtPd < __lpEnd; ++lpTgtPd, ++lpTgtPw, ++lpSrcReal, ++lpSrcImag ) {

        *lpTgtPw = (*lpSrcReal)*(*lpSrcReal) + (*lpSrcImag)*(*lpSrcImag);
        *lpTgtPd = ( cdMult * ::log10(*lpTgtPw) );

        if ( __MaxPd < *lpTgtPd ) { __MaxPd = *lpTgtPd; }
    }
    if ( lpMaxPd ) { *lpMaxPd = __MaxPd; }


    __MaxPd = cdBias -__MaxPd;

    lpTgtPd = __lpBegin;
    for ( ; lpTgtPd < __lpEnd; ++lpTgtPd ) {

        *lpTgtPd += __MaxPd;
    }

_QUIT_:
    return lpTgtPd;
}

// -----------------------------------------------------------------------------
// DllMain
// -----------------------------------------------------------------------------
/*
BOOL APIENTRY
DllMain(
    HANDLE  hModule,        // [IN ] process handle to the dll
    DWORD   dwReason,       // [IN ] reason for calling the function
    LPVOID  lpReserved      // [IN ] reserved, not used
    )                       // [RET] TRUE) success, FALSE) failure
{
    return TRUE;
}
*/

// EOF <fft.cpp>