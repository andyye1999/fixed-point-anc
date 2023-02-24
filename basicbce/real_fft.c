

#include "real_fft.h"



// The C version FFT functions (i.e. RealForwardFFTC and
// RealInverseFFTC) are real-valued FFT wrappers for complex-valued
// FFT implementation in SPL.

int RealForwardFFTC(int order,
	const Word16* real_data_in,
	Word16* complex_data_out) {
	int i = 0;
	int j = 0;
	int result = 0;
	int n = 1 << order;
	// The complex-value FFT implementation needs a buffer to hold 2^order
	// 16-bit COMPLEX numbers, for both time and frequency data.
	Word16 complex_buffer[2 << 10];

	// Insert zeros to the imaginary parts for complex forward FFT input.
	for (i = 0, j = 0; i < n; i += 1, j += 2) {
		complex_buffer[j] = real_data_in[i];
		complex_buffer[j + 1] = 0;
	};

	ComplexBitReverse(complex_buffer, order);
	result = ComplexFFT(complex_buffer, order);

	// For real FFT output, use only the first N + 2 elements from
	// complex forward FFT.
	for (i = 0; i < n + 2; i++)
	{
		complex_data_out[i] = complex_buffer[i];
	}
	//memcpy(complex_data_out, complex_buffer, sizeof(Word16) * (n + 2));

	return result;
}



void ComplexBitReverse(Word16* complex_data, int stages) {
	/* For any specific value of stages, we know exactly the indexes that are
	 * bit reversed. Currently (Feb. 2012) in WebRTC the only possible values of
	 * stages are 7 and 8, so we use tables to save unnecessary iterations and
	 * calculations for these two cases.
	 */

	int m, length;
	const Word16* index;
	if (stages == 7)
	{
		length = 112;
		index = index_7;
	}
	else if (stages == 8)
	{
		length = 240;
		index = index_8;
	}
	/* Decimation in time. Swap the elements with bit-reversed indexes. */
	for (m = 0; m < length; m += 2) {
		/* We declare a int32_t* type pointer, to load both the 16-bit real
		 * and imaginary elements from complex_data in one instruction, reducing
		 * complexity.
		 */
		int* complex_data_ptr = (int*)complex_data;
		int temp = 0;

		temp = complex_data_ptr[index[m]];  /* Real and imaginary */
		complex_data_ptr[index[m]] = complex_data_ptr[index[m + 1]];
		complex_data_ptr[index[m + 1]] = temp;
	}
}

int ComplexFFT(short frfi[], int stages)
{
	int i, j, l, k, istep, n, m;
	short wr, wi;
	int tr32, ti32, qr32, qi32;

	/* The 1024-value is a constant given from the size of kSinTable1024[],
	 * and should not be changed depending on the input parameter 'stages'
	 */
	n = 1 << stages;
	if (n > 1024)
		return -1;

	l = 1;
	k = 10 - 1; /* Constant for given kSinTable1024[]. Do not change
		 depending on the input parameter 'stages' */


		 // mode==1: High-complexity and High-accuracy mode
	while (l < n)
	{
		istep = l << 1;

		for (m = 0; m < l; ++m)
		{
			j = m << k;

			/* The 256-value is a constant given as 1/4 of the size of
			 * kSinTable1024[], and should not be changed depending on the input
			 * parameter 'stages'. It will result in 0 <= j < N_SINE_WAVE/2
			 */
			wr = kSinTable1024[j + 256];
			wi = -kSinTable1024[j];

			for (i = m; i < n; i += istep)
			{
				j = i + l;


				tr32 = MUL_16_16(wr, frfi[2 * j])
					- MUL_16_16(wi, frfi[2 * j + 1]) + CFFTRND;

				ti32 = MUL_16_16(wr, frfi[2 * j + 1])
					+ MUL_16_16(wi, frfi[2 * j]) + CFFTRND;


				tr32 = RSHIFT_W32(tr32, 15 - CFFTSFT);
				ti32 = RSHIFT_W32(ti32, 15 - CFFTSFT);

				qr32 = ((int32_t)frfi[2 * i]) << CFFTSFT;
				qi32 = ((int32_t)frfi[2 * i + 1]) << CFFTSFT;

				frfi[2 * j] = (int16_t)RSHIFT_W32(
					(qr32 - tr32 + CFFTRND2), 1 + CFFTSFT);
				frfi[2 * j + 1] = (int16_t)RSHIFT_W32(
					(qi32 - ti32 + CFFTRND2), 1 + CFFTSFT);
				frfi[2 * i] = (int16_t)RSHIFT_W32(
					(qr32 + tr32 + CFFTRND2), 1 + CFFTSFT);
				frfi[2 * i + 1] = (int16_t)RSHIFT_W32(
					(qi32 + ti32 + CFFTRND2), 1 + CFFTSFT);
			}
		}

		--k;
		l = istep;
	}

	return 0;
}

int RealInverseFFTC(short order,
	const int16_t* complex_data_in,
	int16_t* real_data_out) {
	int i = 0;
	int j = 0;
	int result = 0;
	int n = 1 << order;
	// Create the buffer specific to complex-valued FFT implementation.
	int16_t complex_buffer[2 << 10];

	// For n-point FFT, first copy the first n + 2 elements into complex
	// FFT, then construct the remaining n - 2 elements by real FFT's
	// conjugate-symmetric properties.
	for (i = 0; i < n + 2; i++)
	{
		complex_buffer[i] = complex_data_in[i];
	}
	//memcpy(complex_buffer, complex_data_in, sizeof(int16_t) * (n + 2));
	for (i = n + 2; i < 2 * n; i += 2) {
		complex_buffer[i] = complex_data_in[2 * n - i];
		complex_buffer[i + 1] = -complex_data_in[2 * n - i + 1];
	}

	ComplexBitReverse(complex_buffer, order);
	result = ComplexIFFT(complex_buffer, order);

	// Strip out the imaginary parts of the complex inverse FFT output.
	for (i = 0, j = 0; i < n; i += 1, j += 2) {
		real_data_out[i] = complex_buffer[j];
	}

	return result;
}

int ComplexIFFT(int16_t frfi[], int stages)
{
	int i, j, l, k, istep, n, m, scale, shift;
	int16_t wr, wi;
	int32_t tr32, ti32, qr32, qi32;
	int32_t tmp32, round2;

	/* The 1024-value is a constant given from the size of kSinTable1024[],
	 * and should not be changed depending on the input parameter 'stages'
	 */
	n = 1 << stages;
	if (n > 1024)
		return -1;

	scale = 0;

	l = 1;
	k = 10 - 1; /* Constant for given kSinTable1024[]. Do not change
		 depending on the input parameter 'stages' */

	while (l < n)
	{
		// variable scaling, depending upon data
		shift = 0;
		round2 = 8192;

		tmp32 = (int32_t)MaxAbsValueW16C(frfi, 2 * n);
		if (tmp32 > 13573)
		{
			shift++;
			scale++;
			round2 <<= 1;
		}
		if (tmp32 > 27146)
		{
			shift++;
			scale++;
			round2 <<= 1;
		}

		istep = l << 1;




		// mode==1: High-complexity and High-accuracy mode

		for (m = 0; m < l; ++m)
		{
			j = m << k;

			/* The 256-value is a constant given as 1/4 of the size of
			 * kSinTable1024[], and should not be changed depending on the input
			 * parameter 'stages'. It will result in 0 <= j < N_SINE_WAVE/2
			 */
			wr = kSinTable1024[j + 256];
			wi = kSinTable1024[j];



			for (i = m; i < n; i += istep)
			{
				j = i + l;



				tr32 = MUL_16_16(wr, frfi[2 * j])
					- MUL_16_16(wi, frfi[2 * j + 1]) + CIFFTRND;

				ti32 = MUL_16_16(wr, frfi[2 * j + 1])
					+ MUL_16_16(wi, frfi[2 * j]) + CIFFTRND;

				tr32 = RSHIFT_W32(tr32, 15 - CIFFTSFT);
				ti32 = RSHIFT_W32(ti32, 15 - CIFFTSFT);

				qr32 = ((int32_t)frfi[2 * i]) << CIFFTSFT;
				qi32 = ((int32_t)frfi[2 * i + 1]) << CIFFTSFT;

				frfi[2 * j] = (int16_t)RSHIFT_W32((qr32 - tr32 + round2),
					shift + CIFFTSFT);
				frfi[2 * j + 1] = (int16_t)RSHIFT_W32(
					(qi32 - ti32 + round2), shift + CIFFTSFT);
				frfi[2 * i] = (int16_t)RSHIFT_W32((qr32 + tr32 + round2),
					shift + CIFFTSFT);
				frfi[2 * i + 1] = (int16_t)RSHIFT_W32(
					(qi32 + ti32 + round2), shift + CIFFTSFT);
			}
		}


		--k;
		l = istep;
	}
	return scale;
}
// Maximum absolute value of word16 vector. C version for generic platforms.
int16_t MaxAbsValueW16C(const int16_t* vector, int length) {
	int i = 0, absolute = 0, maximum = 0;

	if (vector == NULL || length <= 0) {
		return -1;
	}

	for (i = 0; i < length; i++) {
		absolute = abs((int)vector[i]);

		if (absolute > maximum) {
			maximum = absolute;
		}
	}

	// Guard the case for abs(-32768).
	if (maximum > 32767) {
		maximum = 32767;
	}

	return (int16_t)maximum;
}