#include "anc.h"
#include "real_fft.h"
#include "inl.h"


#define	FAST_STEP_SIZE	0.05f		 // Fast update step of active noise control.
#define	SLOW_STEP_SIZE	0.005f		 // Slow update step of active noise control.




extern Flag Overflow;
static Word16 windowW16[DELAY + FRM_LEN] = { 17, 157, 434, 844, 1380, 2032, 2790, 3640, 4568, 5558, 6593, 7656, 8727, 9790, 10825, 11815, 12743, 13593, 14351, 15003, 15539, 15949, 16226, 16366, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16366, 16226, 15949, 15539, 15003, 14351, 13593, 12743, 11815, 10825, 9790, 8727, 7656, 6593, 5558, 4568, 3640, 2790, 2032, 1380, 844, 434, 157, 17 };

void ANC_init(ANC_STRUCT *st)
{
	Word16 i;
	// Initialize
	for (st->send_pre_emp_mem = 0, i = 0; i < DELAY; i++) st->send_in_overlap[i] = 0;	
	for (st->send_de_emp_mem = 0, i = 0; i < FFT_LEN - FRM_LEN; i++) st->send_out_overlap[i] = 0;
	for (st->ref_pre_emp_mem = 0, i = 0; i < DELAY; i++) st->ref_in_overlap[i] = 0;	
	for (i = 0; i < FFT_LEN; i++) st->rab[i] = 0;
	for (i = 0; i < FFT_LEN / 2; i++) st->raa[i] = 0;
	st->normData_send = 0;
	st->normData_ref = 0;
	st->Qrab = 0;
	st->Qraa = 0;
	st->anc_flag = 0;
}

void ANC_run(ANC_STRUCT* st, short* send_in, short* ref_in, short* send_out)
{
	Word16 i, j, send_dat_buf[FFT_LEN], ref_dat_buf[FFT_LEN], vv, smax, maxWinData, send_real_fft[FFT_LEN + 2], ref_real_fft[FFT_LEN + 2], normData_send, normData_ref;
	int send_outIFFT;
	Word16 send_update_flag, index_cnt, sendshift, refshift, rabshift, alpha, err[FFT_LEN];
	// Pre-emphasize
	send_dat_buf[DELAY] = send_in[0] - (Word16)(((Word32)st->send_pre_emp_mem * EMP_FAC_Q11 + 1024L) >> 11);
	for (i = 1; i < FRM_LEN; i++) send_dat_buf[DELAY + i] = send_in[i] -  (Word16)(((Word32)send_in[i - 1] * EMP_FAC_Q11 + 1024L) >> 11);
	st->send_pre_emp_mem = send_in[FRM_LEN - 1];
	for (i = 0; i < DELAY; i++) send_dat_buf[i] = st->send_in_overlap[i];
	for (i = 0; i < DELAY; i++) st->send_in_overlap[i] = send_dat_buf[FRM_LEN + i];
	for (i = 0; i < FRM_LEN + DELAY; i++) send_dat_buf[i] = (Word16)(((Word32)send_dat_buf[i] * windowW16[i] + 8192L) >> 14);
	for (i = DELAY + FRM_LEN; i < FFT_LEN; i++) send_dat_buf[i] = 0;
	for (smax = i = 0; i < FFT_LEN; i++)
	{
		smax = max(smax, (Word16)labs(send_dat_buf[i]));
	}
	maxWinData = (short)min(smax, 32767);
	normData_send = _NormW16(maxWinData);



	// Pre-emphasize
	ref_dat_buf[DELAY] = ref_in[0] - (Word16)(((Word32)st->ref_pre_emp_mem * EMP_FAC_Q11 + 1024L) >> 11);
	for (i = 1; i < FRM_LEN; i++) ref_dat_buf[DELAY + i] = ref_in[i] - (Word16)(((Word32)ref_in[i - 1] * EMP_FAC_Q11 + 1024L) >> 11);
	st->ref_pre_emp_mem = ref_in[FRM_LEN - 1];
	for (i = 0; i < DELAY; i++) ref_dat_buf[i] = st->ref_in_overlap[i];
	for (i = 0; i < DELAY; i++) st->ref_in_overlap[i] = ref_dat_buf[FRM_LEN + i];
	for (i = 0; i < FRM_LEN + DELAY; i++) ref_dat_buf[i] = (Word16)(((Word32)ref_dat_buf[i] * windowW16[i] + 8192L) >> 14);
	for (i = DELAY + FRM_LEN; i < FFT_LEN; i++) ref_dat_buf[i] = 0;
	for (smax = i = 0; i < FFT_LEN; i++)
	{
		smax = max(smax, (Word16)labs(ref_dat_buf[i]));
	}
	maxWinData = (short)min(smax, 32767);
	normData_ref = _NormW16(maxWinData);

	for (i = 0; i < FFT_LEN; i++)
	{
		send_dat_buf[i] = send_dat_buf[i] << normData_send;
		ref_dat_buf[i] = ref_dat_buf[i] << normData_ref;
	}
	
	RealForwardFFTC(8, send_dat_buf, send_real_fft);  // Q(normsend-8)  Ô­À´µÄFFTQ(-7)
	RealForwardFFTC(8, ref_dat_buf, ref_real_fft);   // Q(normref-8)

	send_dat_buf[0] = send_real_fft[0];
	send_dat_buf[1] = send_real_fft[FFT_LEN];
	for (i = 2; i < FFT_LEN; i++)
	{
		send_dat_buf[i] = send_real_fft[i];
	}
	ref_dat_buf[0] = ref_real_fft[0];
	ref_dat_buf[1] = ref_real_fft[FFT_LEN];
	for (i = 2; i < FFT_LEN; i++)
	{
		ref_dat_buf[i] = ref_real_fft[i];
	}
	alpha = 163;
	Word16 onesubalpha = 32604;
	
	
	for (smax = 0, i = 120; i < FFT_LEN; i++)
	{
		smax= max(smax, (Word16)labs(ref_dat_buf[i]));
	}
	Word16 normfftref;
	normfftref = _NormW16(smax);
	for (smax =0, i = 120; i < FFT_LEN; i++)
	{
		ref_dat_buf[i] = shl(ref_dat_buf[i], normfftref);
	}
	for (smax = 0,i = 120; i < FFT_LEN; i++)
	{
		smax = max(smax, (Word16)labs(send_dat_buf[i]));
	}
	Word16 normfftsend;
	normfftsend = _NormW16(smax);
	for (smax = 0, i = 120; i < FFT_LEN; i++)
	{
		send_dat_buf[i] = shl(send_dat_buf[i], normfftsend);
	}
	Word32 t0,t1,t2,t3,t4,raa_div;
	Word16 b0_h, b0_l, b1_h, b1_l, b2_h, b2_l, x0, x1, x2, Qtmp1, b3_h, b3_l, Qtmp2, Qtmp3, Qtmp4, x3;
	Word16 Qrab = normData_ref + normData_send - 16 + normfftref + normfftsend;
	rabshift = st->Qrab - Qrab;
	Word16 Qraa = normData_ref + normData_ref - 16 + normfftref + normfftref;
	Word16 raashift = st->Qraa - Qraa;
	for (j = 60; j < FFT_LEN / 2; j++)
	{
		Overflow = 0;
		t0 = L_mult(ref_dat_buf[2 * j], send_dat_buf[2 * j]);
		t0 = L_mac(t0, ref_dat_buf[2 * j + 1], send_dat_buf[2 * j + 1]);
		t0 = L_shr(t0, 1);
		L_Extract(t0, &b0_h, &b0_l);
		t0 = Mpy_32_16(b0_h, b0_l, alpha);//normsend+normref-16+normfftsend+normfftref
		t1 = L_shr(st->rab[2 * j], rabshift);
		L_Extract(t1, &b1_h, &b1_l);
		t1 = Mpy_32_16(b1_h, b1_l, onesubalpha);//normsend+normref-16+normfftsend+normfftref
		st->rab[2 * j] = L_add(t0, t1);

		t0 = L_mult(ref_dat_buf[2 * j], send_dat_buf[2 * j + 1]);
		t0 = L_msu(t0, ref_dat_buf[2 * j + 1], send_dat_buf[2 * j]);
		t0 = L_shr(t0, 1);
		L_Extract(t0, &b0_h, &b0_l);
		t0 = Mpy_32_16(b0_h, b0_l, alpha);//normsend+normref-16+normfftsend+normfftref
		t1 = L_shr(st->rab[2 * j + 1], rabshift);
		L_Extract(t1, &b1_h, &b1_l);
		t1 = Mpy_32_16(b1_h, b1_l, onesubalpha);
		st->rab[2 * j + 1] = L_add(t0, t1);//normsend+normref-16+normfftsend+normfftref

		t0 = L_mult(ref_dat_buf[2 * j], ref_dat_buf[2 * j]);
		t0 = L_mac(t0, ref_dat_buf[2 * j + 1], ref_dat_buf[2 * j + 1]);
		t0 = L_shr(t0, 1);
		L_Extract(t0, &b0_h, &b0_l);
		t0 = Mpy_32_16(b0_h, b0_l, alpha);//normref+normref-16+normfftref+normfftref
		t1 = L_shr(st->raa[j], raashift);
		L_Extract(t1, &b1_h, &b1_l);
		t1 = Mpy_32_16(b1_h, b1_l, onesubalpha);
		st->raa[j] = L_add(t0, t1);//normref+normref-16+normfftref+normfftref

		t1 = L_abs(st->rab[2 * j]);
		t2 = L_abs(st->rab[2 * j + 1]);
		x1 = norm_l(t1);
		x2 = norm_l(t2);
		t1 = L_shl(st->rab[2 * j], x1);
		t2 = L_shl(st->rab[2 * j + 1], x2);
		L_Extract(t1, &b1_h, &b1_l);
		t1 = Mpy_32_16(b1_h, b1_l, ref_dat_buf[2 * j]);
		L_Extract(t2, &b2_h, &b2_l);
		t2 = Mpy_32_16(b2_h, b2_l, ref_dat_buf[2 * j + 1]);
		if (x1 > x2)
		{
			t1 = L_shr(t1, x1 - x2);
			x1 = x2;
		}
		else if (x1 < x2)
		{
			t2 = L_shr(t2, x2 - x1);
			x2 = x1;
		}
		t0 = L_sub(t1, t2);
		Qtmp1 = Qrab + x1 + normData_ref + normfftref - 8 - 15;

		t1 = L_abs(t0);
		x1 = norm_l(t1);
		Qtmp2 = Qtmp1 + x1;
		t1 = L_shl(t0, x1);
		L_Extract(t1, &b1_h, &b1_l);

		t1 = L_abs(st->raa[j]);
		x1 = norm_l(t1);
		Qtmp3 = x1;
		t1 = L_shl(t1, x1);
		L_Extract(t1, &b2_h, &b2_l);
		t4 = L_shl(1L, Qtmp3);
		t2 = Div_32(t4, b2_h, b2_l);
		raa_div = t2;
		L_Extract(t2, &b2_h, &b2_l);
		t3 = Mpy_32(b1_h, b1_l, b2_h, b2_l);
		Qtmp3 = Qtmp2 - Qraa - (normData_send + normfftsend - 8);
		t3 = L_shr(t3, Qtmp3);
		x3 = extract_l(t3);
		L_Extract(t3, &b3_h, &b3_l);
		err[2 * j] = sub(send_dat_buf[2 * j], x3);

		t1 = L_abs(st->rab[2 * j]);
		t2 = L_abs(st->rab[2 * j + 1]);
		x1 = norm_l(t1);
		x2 = norm_l(t2);
		t1 = L_shl(st->rab[2 * j], x1);
		t2 = L_shl(st->rab[2 * j + 1], x2);
		L_Extract(t1, &b1_h, &b1_l);
		t1 = Mpy_32_16(b1_h, b1_l, ref_dat_buf[2 * j + 1]);
		L_Extract(t2, &b2_h, &b2_l);
		t2 = Mpy_32_16(b2_h, b2_l, ref_dat_buf[2 * j]);
		if (x1 > x2)
		{
			t1 = L_shr(t1, x1 - x2);
			x1 = x2;
		}
		else if (x1 < x2)
		{
			t2 = L_shr(t2, x2 - x1);
			x2 = x1;
		}
		t0 = L_add(t1, t2);
		Qtmp1 = Qrab + x1 + normData_ref + normfftref - 8 - 15;

		t1 = L_abs(t0);
		x1 = norm_l(t1);
		Qtmp2 = Qtmp1 + x1;
		t1 = L_shl(t0, x1);
		L_Extract(t1, &b1_h, &b1_l);
		L_Extract(raa_div, &b2_h, &b2_l);
		t3 = Mpy_32(b1_h, b1_l, b2_h, b2_l);
		Qtmp3 = Qtmp2 - Qraa - (normData_send + normfftsend - 8);
		t3 = L_shr(t3, Qtmp3);
		x3 = extract_l(t3);
		L_Extract(t3, &b3_h, &b3_l);
		err[2 * j + 1] = sub(send_dat_buf[2 * j + 1], x3);
		
	}
	if (st->anc_flag == 1)
	{
		for (j = 60; j < FFT_LEN / 2; j++)
		{
			
			t0 = L_mult(err[2 * j], err[2 * j]);
			t0 = L_mac(t0, err[2 * j + 1], err[2 * j + 1]);
			t0 = L_shr(t0, 1);
			t1 = L_mult(send_dat_buf[2 * j], send_dat_buf[2 * j]);
			t1 = L_mac(t1, send_dat_buf[2 * j + 1], send_dat_buf[2 * j + 1]);
			t1 = L_shr(t1, 1);
			if (t0 < t1)
			{
				send_dat_buf[2 * j] = (Word16)err[2 * j];
				send_dat_buf[2 * j + 1] = (Word16)err[2 * j + 1];
			}
		}
	}
	for (i = 120; i < FFT_LEN; i++)
	{
		send_dat_buf[i] = shr(send_dat_buf[i], normfftsend);
	}

	send_real_fft[0] = send_dat_buf[0];
	send_real_fft[FFT_LEN] = send_dat_buf[1];
	for (i = 2; i < FFT_LEN; i++)
	{
		send_real_fft[i] = send_dat_buf[i];
	}

	send_outIFFT = RealInverseFFTC(8, send_real_fft, send_dat_buf);

	for (i = 0; i < FFT_LEN; i++)
	{
		send_dat_buf[i] = send_dat_buf[i] << send_outIFFT >> normData_send;
	}
	st->normData_send = normData_send;
	st->normData_ref = normData_ref;
	st->Qrab = Qrab;
	st->Qraa = Qraa;

	// Overlap add 
	for (i = 0; i < FFT_LEN - FRM_LEN; i++) send_dat_buf[i] += st->send_out_overlap[i];
	for (i = 0; i < FFT_LEN - FRM_LEN; i++) st->send_out_overlap[i] = send_dat_buf[i + FRM_LEN];
	
	for (vv = st->send_de_emp_mem, i = 0; i < FRM_LEN; i++)
	{
		vv = send_dat_buf[i] + (Word16)((vv * EMP_FAC_Q11) >> 11);
		send_out[i] = max(-32768.0f, min(32767.0f, vv));
	}
	st->send_de_emp_mem = vv;
}
