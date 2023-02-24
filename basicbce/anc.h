#ifndef ANC_H_
#define ANC_H_

#include "enhance.h"

typedef struct
{
	Word16 send_pre_emp_mem;
	Word16 send_in_overlap[DELAY];
	Word16 send_de_emp_mem;
	Word16 send_out_overlap[FFT_LEN - FRM_LEN];

	Word16 ref_pre_emp_mem;
	Word16 ref_in_overlap[DELAY];


	Word32 rab[FFT_LEN];
	Word32 raa[FFT_LEN / 2];
	Word16 anc_flag;
	Word16 normData_send;
	Word16 normData_ref;
	Word16 Qrab;
	Word16 Qraa;
}ANC_STRUCT;

void ANC_init(ANC_STRUCT *st);
void ANC_run(ANC_STRUCT *st, short *send_in, short*ref_in, short*send_out);

#endif
