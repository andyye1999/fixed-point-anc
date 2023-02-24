#include "anc.h"
#include <time.h>


void main(void)
{
	Word32 frame, i;
	Word16 in_send[FRM_LEN], in_air[FRM_LEN], in_receive[FRM_LEN], out_send[FRM_LEN], out_receive[FRM_LEN];
	FILE* fq, * fc, * fs;
	clock_t start, finish;
	double Total_time;

	ANC_STRUCT anc_st;
	

	// Initialize
	ANC_init(&anc_st);
	// Function flag
	anc_st.anc_flag = 1;
	fc = fopen("ref.pcm", "rb");		// Ref
	fq = fopen("input.pcm", "rb");	// input
	fs = fopen("out.pcm", "wb");				

	if ((fq == NULL) || (fc == NULL)) { printf("Can't open file!!!\n");  exit(0); }

	start = clock();
	for (frame = 0; ; frame++)
	{
		if (fread(in_send, sizeof(short), FRM_LEN, fq) != FRM_LEN) break;
		if (fread(in_air, sizeof(short), FRM_LEN, fc) != FRM_LEN) break;
		if ((frame & 4095L) == 0)
		{
			printf("Frame=%ld\n", frame);
		}
		ANC_run(&anc_st, in_send, in_air, out_send);
		fwrite(out_send, sizeof(short), FRM_LEN, fs);
	}
	finish = clock();
	Total_time = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("\n%f seconds\n", Total_time);
	fclose(fq);
	fclose(fc);
	fclose(fs);
}
