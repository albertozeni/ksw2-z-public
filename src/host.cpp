#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>
#include <bitset>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h> 
#include <sstream>
#include <set>
#include <memory>
#include <typeinfo>
#include <pthread.h>
#include <vector>
#include <functional>
#include <iterator>
#include <string.h>
#include "xcl2.hpp"
#include <cstring>
#include <time.h>
#include <assert.h>
#include "ap_int.h"
#include "datatypes.h"

#define NUM_KERNEL 10

//HBM Banks requirements
#define MAX_HBM_BANKCOUNT 32
#define BANK_NAME(n) n | XCL_MEM_TOPOLOGY
const int bank[MAX_HBM_BANKCOUNT] = {
        BANK_NAME(0),  BANK_NAME(1),  BANK_NAME(2),  BANK_NAME(3),  BANK_NAME(4),
        BANK_NAME(5),  BANK_NAME(6),  BANK_NAME(7),  BANK_NAME(8),  BANK_NAME(9),
        BANK_NAME(10), BANK_NAME(11), BANK_NAME(12), BANK_NAME(13), BANK_NAME(14),
        BANK_NAME(15), BANK_NAME(16), BANK_NAME(17), BANK_NAME(18), BANK_NAME(19),
        BANK_NAME(20), BANK_NAME(21), BANK_NAME(22), BANK_NAME(23), BANK_NAME(24),
        BANK_NAME(25), BANK_NAME(26), BANK_NAME(27), BANK_NAME(28), BANK_NAME(29),
        BANK_NAME(30), BANK_NAME(31)};

#define NOW std::chrono::high_resolution_clock::now()

using namespace std;


using namespace chrono;



//=======================================================================
//
//  define 
//
//=======================================================================

#define KSW_NEG_INF -0x40000000

#define KSW_EZ_SCORE_ONLY  0x01 // don't record alignment path/cigar
#define KSW_EZ_RIGHT       0x02 // right-align gaps
#define KSW_EZ_GENERIC_SC  0x04 // without this flag: match/mismatch only; last symbol is a wildcard
#define KSW_EZ_APPROX_MAX  0x08 // approximate max; this is faster with sse
#define KSW_EZ_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define KSW_EZ_EXTZ_ONLY   0x40 // only perform extension
#define KSW_EZ_REV_CIGAR   0x80 // reverse CIGAR in the output
#define KSW_EZ_SPLICE_FOR  0x100
#define KSW_EZ_SPLICE_REV  0x200
#define KSW_EZ_SPLICE_FLANK 0x400
#define N_PAIRS 	   10000 // PAIRS per core (tot pairs in total is N_PAIRS*10 = 100000

void ksw_extz2_hw(int qlen, input_datatype *query_l, int tlen, input_datatype *target_l, datatype m,
				 datatype q, datatype e, int w, int zdrop, int end_bonus, int flag, std::vector<int, aligned_allocator<int>> *v_s, int n_pairs, std::string binary_file);

typedef struct { //
	uint32_t max:31, zdropped:1;
	int max_q, max_t;      // max extension coordinate
	int mqe, mqe_t;        // max score when reaching the end of query
	int mte, mte_q;        // max score when reaching the end of target
	int score;             // max score reaching both ends; may be KSW_NEG_INF
	int m_cigar, n_cigar;
	int reach_end;
	uint32_t *cigar;
} ksw_extz_t;

void copy_ez(ksw_hw_extz_t *ez_out, ksw_extz_t *ez){
//#pragma HLS INLINE
	ez_out->max_q=ez->max_q;
	ez_out->max_t=ez->max_t;
	ez_out->mqe_t=ez->mqe_t;
	ez_out->mte_q=ez->mte_q;
	ez_out->max=ez->max;
	ez_out->score=ez->score;
	ez_out->mqe=ez->mqe;
	ez_out->mte=ez->mte;
	ez_out->m_cigar=ez->m_cigar;
	ez_out->n_cigar=ez->n_cigar;
	ez_out->zdropped=ez->zdropped;
	ez_out->reach_end=ez->reach_end;
}


static inline void ksw_reset_extz(ksw_extz_t *ez)
{
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
	ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
}

static inline int ksw_apply_zdrop(ksw_extz_t *ez, int is_rot, int32_t H, int a, int b, int zdrop, int8_t e)
{
	int r, t;
	if (is_rot) r = a, t = b;
	else r = a + b, t = a;
	if (H > (int32_t)ez->max) {
		ez->max = H, ez->max_t = t, ez->max_q = r - t;
	} else if (t >= ez->max_t && r - t >= ez->max_q) {
		int tl = t - ez->max_t, ql = (r - t) - ez->max_q, l;
		l = tl > ql? tl - ql : ql - tl;
		if (zdrop >= 0 && ez->max - H > zdrop + l * e) {
			ez->zdropped = 1;
			return 1;
		}
	}
	return 0;
}

static inline uint32_t *ksw_push_cigar_new(int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = (uint32_t*)realloc(cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

static inline void ksw_backtrack_new(int is_rot, int is_rev, int min_intron_len, const uint8_t *p, const int *off, const int *off_end, int n_col, int i0, int j0,
								 int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{ // p[] - lower 3 bits: which type gets the max; bit
	int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
	uint32_t *cigar = *cigar_, tmp;
	while (i >= 0 && j >= 0) { // at the beginning of the loop, _state_ tells us which state to check
		int force_state = -1;
		if (is_rot) {
			r = i + j;
			if (i < off[r]) force_state = 2;
			if (off_end && i > off_end[r]) force_state = 1;
			tmp = force_state < 0? p[(size_t)r * n_col + i - off[r]] : 0;
		} else {
			if (j < off[i]) force_state = 2;
			if (off_end && j > off_end[i]) force_state = 1;
			tmp = force_state < 0? p[(size_t)i * n_col + j - off[i]] : 0;
		}
		if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
		else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
		if (state == 0) state = tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
		if (force_state >= 0) state = force_state;
		if (state == 0) cigar = ksw_push_cigar_new(&n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
		else if (state == 1 || (state == 3 && min_intron_len <= 0)) cigar = ksw_push_cigar_new(&n_cigar, &m_cigar, cigar, 2, 1), --i; // deletion
		else if (state == 3 && min_intron_len > 0) cigar = ksw_push_cigar_new(&n_cigar, &m_cigar, cigar, 3, 1), --i; // intron
		else cigar = ksw_push_cigar_new(&n_cigar, &m_cigar, cigar, 1, 1), --j; // insertion
	}
	if (i >= 0) cigar = ksw_push_cigar_new(&n_cigar, &m_cigar, cigar, min_intron_len > 0 && i >= min_intron_len? 3 : 2, i + 1); // first deletion
	if (j >= 0) cigar = ksw_push_cigar_new(&n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
	if (!is_rev)
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
	*m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}

//=======================================================================
//
// Software function
//
//=======================================================================


////funzioni SSE 4.1 riscritte


void mm_set1_epi8 ( char *ris, int8_t v) {  //mette tutti i valori del vettore di char ris uguali a v
// ris = vettore in cui va il risultato
//numero da mettere nel vettore

    int i;
    for(i=0; i<16; i++ ){
	    ris[i] = v;
	}
}

void mm_blendv_epi8 ( char *ris,  char *a,  char *b,  char* mask ){ //sceglie se tenere il valore di b se mask=0 e il valore di a se mask=1 a pacchetti di 8 bit
//ris = risultato
//a e b sono i due vettori da mischiare
//mask  il vettore che fa da maschera (0 e 1)

    int i;
    for(i=0; i<16; i++){
        if (mask[i])    ris[i]=b[i];
        else            ris[i]=a[i];
    }
}

void mm_cmpeq_epi8( char* ris,  char* a,  char* b){  //compara a e b e mette ris=255 se vero e ris=1 se falso

    int i;
    for(i=0; i<16; i++){
        if(a[i]==b[i])
            ris[i]=255;
        else
            ris[i]=0;
    }
}

void mm_cmpgt_epi8(char* ris, char* a, char* b){  //compara a e b e mette ris=255 se vero e ris=1 se falso

    int i;
    for(i=0; i<16; i++){
        if(a[i]>b[i])
            ris[i]=0xff;
        else
            ris[i]=0;
    }
}

void mm_loadu_si128_8( char* ris,  char* a){    //copia a in ris
    int i;
    for(i=0; i<16; i++){
		ris[i]=a[i];
	}
}

void mm_or_si128_8(  char* ris,   char* a,  char* b){  //ora tra a e b

    int i;
    for(i=0; i<16; i++){
		ris[i]=a[i]|b[i];
	}
}

void mm_and_si128_8( char* ris,  char* a,  char* b){  //ora tra a e b

    int i;
    for(i=0; i<16; i++){
		ris[i]=a[i]&b[i];
	}
}

void mm_max_epi8( char* ris,  char*a,  char* b){       //restituisce il massimo tra a e b

    int i;
    for(i=0; i<16; i++){
        if(a[i]>=b[i])      ris[i]=a[i];
        else                ris[i]=b[i];
    }

}

void mm_min_epi8( char* ris,  char*a,  char* b){       //restituisce il minimo tra a e b

    int i;
    for(i=0; i<16; i++){
        if(a[i]<=b[i])      ris[i]=a[i];
        else                ris[i]=b[i];
    }

}

void mm_add_epi8( char* ris,  char* a,  char* b) {
    int i;
    for(i=0; i<16; i++) {
		ris[i]=a[i]+b[i];
	}
}

void mm_sub_epi8( char* ris,  char* a,  char* b) {
    int i;
    for(i=0; i<16; i++) {
		ris[i]=a[i]-b[i];
	}
}

void mm_srli_si128_8( char* ris,  char* a, int8_t num){

    int i,j;
    char tmp[16];

    for(i=0; i<16; i++){
		tmp[i]=a[i];
	}

    if (num>16)     num=16;
    else{
        for(j=0; j<num; j++){
            for(i=0; i<15; i++){
                tmp[i]=tmp[i+1];
            }
        tmp[15]=0;
        }
    }

    for(i=0; i<16; i++) {
		ris[i]=tmp[i];
	}
}

void mm_slli_si128_8( char* ris,  char* a, int8_t num){

    int i,j;
    char tmp[16];

    for(i=0; i<16; i++)  {
		tmp[i]=a[i];
	}

    if (num>16)     num=16;
    else{
        for(j=0; j<num; j++){
            for(i=15; i>0; i--){
                tmp[i]=tmp[i-1];
            }
        tmp[0]=0;
        }
    }

    for(i=0; i<16; i++){
		ris[i]=tmp[i];
	}
}


void mm_set1_epi32(int32_t *ris, int32_t a){ //
    int i;
    for(i=0; i<4; i++) {
		ris[i]=a;
	}
}

void mm_setr_epi32(int32_t *ris, int32_t *a){ //inserisce a in ris, ma in ordine inverso

    ris[3]=a[0];
    ris[2]=a[1];
    ris[1]=a[2];
    ris[0]=a[3];
}

void mm_cvtsi32_si128(int32_t * ris, int32_t a){  //copia 32 bit di a nei primi 32 bit di ris, poi tutto a 0
    ris[0]=a;
    ris[1]=0;
    ris[2]=0;
    ris[3]=0;
}

void mm_loadu_si128_32(int32_t* ris, int32_t* a){    //copia a in ris

    int i;
    for(i=0; i<4; i++)  {
		ris[i]=a[i];
	}
}

void mm_add_epi32(int32_t* ris, int32_t* a, int32_t* b) {
    int i;
    for(i=0; i<4; i++) {
		ris[i]=a[i]+b[i];
	}
}

void mm_sub_epi32(int32_t* ris, int32_t* a, int32_t* b) {
    int i;
    for(i=0; i<4; i++)  {
		ris[i]=a[i]-b[i];
	}
}

void mm_cmpgt_epi32(int32_t* ris, int32_t* a, int32_t* b){

    int i;
    for(i=0; i<4; i++){
        if(a[i]>b[i])
            ris[i]= 0xffffffff; //4294967295;
        else
            ris[i]=0;
    }
}

void uguale8(char* ris, char* a){
	int i;
	for(i=0; i<16; i++)	{
		ris[i]=a[i];
	}
}

void to_8(char* ris, int32_t* a){

    memcpy(ris, (char*)a, 4);//always copy 4 bytes
}

void to_32(int32_t* ris, char* a){

    memcpy((char*)ris, a, 4);//always copy 4 bytes
}

void mm_blendv_epi32(int32_t *ris, int32_t *a, int32_t *b, int32_t* mask ){

    int i;
    for(i=0; i<4; i++){
        if (mask[i])    ris[i]=b[i];
        else            ris[i]=a[i];
    }
}


///////

void ksw2_extz2_sw(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop,
			int end_bonus, int flag, ksw_extz_t *ez)
{

	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;

	char q_8[16], qe2_8[16], zero_8[16], flag1_8[16], flag2_8[16], flag8_8[16], flag16_8[16], sc_mch_8[16], sc_mis_8[16], sc_N_8[16], m1_8[16], max_sc_8[16];
	__int128_t *u, *v, *x, *y, *s, *p = 0; //

	ksw_reset_extz(ez);
	if (m <= 0 || qlen <= 0 || tlen <= 0) return;

	mm_set1_epi8(zero_8,0);
	mm_set1_epi8(q_8,q);
	mm_set1_epi8(qe2_8,(q+e)*2);
	mm_set1_epi8(flag1_8,1);
	mm_set1_epi8(flag2_8,2);
	mm_set1_epi8(flag8_8,0x08);
	mm_set1_epi8(flag16_8,0x10);
	mm_set1_epi8(sc_mch_8,mat[0]);
	mm_set1_epi8(sc_mis_8,mat[1]);
	mat[m*m-1]==0? mm_set1_epi8(sc_N_8,-e) : mm_set1_epi8(sc_N_8,mat[m*m-1]);
	mm_set1_epi8(m1_8,m-1);
	mm_set1_epi8(max_sc_8,(mat[0] + (q+e)*2 ));

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	mem = (uint8_t*)calloc(tlen_ * 6 + qlen_ + 1, 16);
	u = (__int128_t*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, s = y + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;


	if (!approx_max) {
		H = (int32_t*)malloc(tlen_ * 16 * 4);
		for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)malloc(((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
		p = (__int128_t*)(((size_t)mem2 + 15) >> 4 << 4);
		off = (int*)malloc((qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);


	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, v1;
		uint8_t *qrr = qr + (qlen - 1 - r), *u8 = (uint8_t*)u, *v8 = (uint8_t*)v;

		 int32_t x1_32[4], v1_32[4];
		 char x1_8[16], v1_8[16];
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en)
				x1 = ((char*)x)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			else x1 = v1 = 0; // not calculated; set to zeros
		} else x1 = 0, v1 = r? q : 0;
		if (en >= r) ((uint8_t*)y)[r] = 0, u8[r] = r? q : 0;
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			for (t = st0; t <= en0; t += 16) {

				 char sq8[16], st8[16], tmp8[16], mask8[16];
				 char extra1[16], extra2[16];

				mm_loadu_si128_8(sq8, (char*)&sf[t]);
				mm_loadu_si128_8(st8, (char*)&qrr[t]);
				mm_cmpeq_epi8(extra1, sq8, m1_8); mm_cmpeq_epi8(extra2, st8, m1_8);
				mm_or_si128_8(mask8, extra1, extra2);

				mm_cmpeq_epi8(tmp8, sq8, st8);

				mm_blendv_epi8(tmp8, sc_mis_8, sc_mch_8, tmp8);
				mm_blendv_epi8(tmp8, tmp8, sc_N_8, mask8);

				mm_loadu_si128_8((char*)s + t , tmp8);


			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}

		// core loop
		mm_cvtsi32_si128(x1_32, x1);
		mm_cvtsi32_si128(v1_32, v1);

		st_ = st / 16, en_ = en / 16;

		to_8(x1_8, x1_32);
		to_8(v1_8, v1_32);



		assert(en_ - st_ + 1 <= n_col_);
		if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__int128_t *pr = p + (size_t)r * n_col_ - st_;  //
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {

				char d8[16], d8p[16], z8[16], z8p[16], xt18[16], xt18p[16], tmp8[16], tmp8p[16], vt18[16], vt18p[16], a8[16], ut8[16], b8[16], b8p[16];

				//__dp_code_block1;

				////////////////////// __dp_block1

				mm_loadu_si128_8(z8p, (char*)&s[t]); //
				mm_add_epi8(z8, z8p, qe2_8);

				mm_loadu_si128_8(xt18, (char*)&x[t]);

				mm_srli_si128_8(tmp8, xt18, 15);

				mm_slli_si128_8(xt18p, xt18, 1);
				mm_or_si128_8(xt18, xt18p, x1_8);

				uguale8(x1_8, tmp8);

				mm_loadu_si128_8(vt18, (char*)&v[t]);

				mm_srli_si128_8(tmp8, vt18, 15);

				mm_slli_si128_8(vt18p, vt18, 1);
				mm_or_si128_8(vt18, vt18p, v1_8);

				uguale8(v1_8, tmp8);

				mm_add_epi8(a8, xt18, vt18);

				mm_loadu_si128_8(ut8, (char*)&u[t]);

				mm_loadu_si128_8(b8p,(char*)&y[t]);
				mm_add_epi8(b8, b8p, ut8);


				/////////////////////

				mm_cmpgt_epi8(d8p, a8, z8);
				mm_and_si128_8(d8, d8p, flag1_8);                          // z = z > a? z : a (signed)

				mm_max_epi8(z8, z8, a8);

				mm_cmpgt_epi8(tmp8, b8, z8);           // d = b > z? 2 : d

				mm_blendv_epi8(d8, d8, flag2_8, tmp8);

				//__dp_code_block2;

				/////////////////////////__dp_code_block2


				mm_max_epi8(z8, z8, b8);

				mm_min_epi8(z8, z8, max_sc_8);

				mm_sub_epi8(z8p, z8, vt18);
				mm_loadu_si128_8((char*)&u[t], z8p);

				mm_sub_epi8(z8p, z8, ut8);
				mm_loadu_si128_8((char*)&v[t], z8p);

				mm_sub_epi8(z8, z8, q_8);

				mm_sub_epi8(a8, a8, z8);

				mm_sub_epi8(b8, b8, z8);


				///////////////////////////////////////

				mm_cmpgt_epi8(tmp8, a8, zero_8);

				mm_and_si128_8(tmp8p, tmp8, a8);
				mm_loadu_si128_8((char*)&x[t], tmp8p);

				mm_and_si128_8(d8p, tmp8, flag8_8);
				mm_or_si128_8(d8, d8, d8p);

				mm_cmpgt_epi8(tmp8, b8, zero_8);

				mm_and_si128_8(tmp8p, tmp8, b8);
				mm_loadu_si128_8((char*)&y[t], tmp8p);

				mm_and_si128_8(tmp8p, tmp8, flag16_8);
				mm_or_si128_8(d8, d8, tmp8p);

				mm_loadu_si128_8((char*)&pr[t], d8);
			}
		}


		int32_t max_H, max_t;
		// compute H[], max_H and max_t
		if (r > 0) {
			int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
			int32_t max_H_32[4], max_t_32[4], qe_32[4];

			max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] - qe : H[en0] + v8[en0] - qe; // special casing the last element
			max_t = en0;

			mm_set1_epi32(max_H_32, max_H);

			mm_set1_epi32(max_t_32, max_t);

			mm_set1_epi32(qe_32, q+e);


			for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;


				int32_t H132[4], t_32[4], tmp32[4];

				mm_loadu_si128_32(H132, (int32_t*)&H[t]);

				t_32[0]=v8[t]; t_32[1]=v8[t+1]; t_32[2]=v8[t+2]; t_32[3]=v8[t+3];

				mm_add_epi32(H132, H132, t_32);

				mm_sub_epi32(H132, H132, qe_32);

				mm_loadu_si128_32((int32_t*)&H[t], H132);

				mm_set1_epi32(t_32, t);

				mm_cmpgt_epi32(tmp32, H132, max_H_32);

				mm_blendv_epi8((char*)max_H_32, (char*) max_H_32, (char*) H132, (char*) tmp32);

				mm_blendv_epi8((char*)max_t_32,(char*) max_t_32, (char*)t_32, (char*)tmp32);

			}

			mm_loadu_si128_8((char*)HH, (char*)max_H_32);

			mm_loadu_si128_8((char*)tt, (char*)max_t_32);


			for (i = 0; i < 4; ++i)
				if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
			for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
				H[t] += (int32_t)v8[t] - qe;
				if (H[t] > max_H)
					max_H = H[t], max_t = t;
			}
		} else H[0] = v8[0] - qe - qe, max_H = H[0], max_t = 0; // special casing r==0
		// update ez
		if (en0 == tlen - 1 && H[en0] > ez->mte)
			ez->mte = H[en0], ez->mte_q = r - en;
		if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
			ez->mqe = H[st0], ez->mqe_t = st0;
		if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e)) break;
		if (r == qlen + tlen - 2 && en0 == tlen - 1)
			ez->score = H[tlen - 1];

		last_st = st, last_en = en;
		//for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
	}


	free(mem);
	if (!approx_max) free(H);
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
			ksw_backtrack_new(1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
			ez->reach_end = 1;
			ksw_backtrack_new(1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack_new(1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		free(mem2); free(off);
	}
}

//=======================================================================
//
// Common functions
//
//=======================================================================


unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}


static void global_aln_sw(const char *qseq_, const char *tseq_, int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
					   int w, int zdrop, int flag, ksw_extz_t *ez)
{
	int i, qlen, tlen;
	uint8_t *qseq, *tseq;
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
	ez->n_cigar = 0;
	qlen = MAX_SEQ_LEN;//strlen(qseq_);
	//printf("\ntarget: %s\n",qseq_); //
	tlen = MAX_SEQ_LEN;//strlen(tseq_);
	//printf("query: %s\n",tseq_); //
	qseq = (uint8_t*)calloc(qlen + 33, 1); // 32 for gaba
	tseq = (uint8_t*)calloc(tlen + 33, 1);
	for (i = 0; i < qlen; ++i)
		qseq[i] = seq_nt4_table[(uint8_t)qseq_[i]];
	for (i = 0; i < tlen; ++i)
		tseq[i] = seq_nt4_table[(uint8_t)tseq_[i]];
	
	/*//stampa delle sequenze
			printf("\n");
			for(i=0; i<qlen; i++)	printf("%d", qseq[i]);
			printf("\n");
			for(i=0; i<tlen; i++)	printf("%d", tseq[i]);
			printf("\n"); */
	flag &= KSW_EZ_SCORE_ONLY; //we are only computing the score
	// flag |= KSW_EZ_GENERIC_SC; //we are using the scoring matrix
	ksw2_extz2_sw(qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, zdrop, 0, flag, ez);


	free(qseq); free(tseq);
}

//////hw
static void global_aln_hw(const char *qseq_, const char *tseq_, int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
					   int w, int zdrop, int flag, ksw_hw_extz_t *ez, std::vector<int, aligned_allocator<int>> *v_s, std::string binary_file)
{
	int i, qlen, tlen;
	// uint8_t *qseq, *tseq;
//	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
//	ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
//	ez->m_cigar = 0, ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
	qlen = (N_PAIRS)*MAX_SEQ_LEN;//strlen(qseq_);
	tlen = (N_PAIRS)*MAX_SEQ_LEN;//strlen(tseq_);
	std::vector<uint8_t, aligned_allocator<uint8_t>> qseq(qlen+MAX_SEQ_LEN);
	std::vector<uint8_t, aligned_allocator<uint8_t>> tseq(tlen+MAX_SEQ_LEN);
	// qseq = (uint8_t*)malloc(qlen+MAX_SEQ_LEN);// + 33, 1); // 32 for gaba
	// tseq = (uint8_t*)malloc(tlen+MAX_SEQ_LEN);// + 33, 1);
	//copy the scoring matrix into qseq, for both we skip a pair of data for the hw core
	for (i = 0; i < m*m; i++){
		qseq[i]=(datatype)mat[i];
	}
	for (i = 0; i < N_PAIRS; ++i){
		for(int j = 0; j < MAX_SEQ_LEN;j++){
			qseq[(1+i)*MAX_SEQ_LEN+j] = seq_nt4_table[(uint8_t)qseq_[(i+1)*MAX_SEQ_LEN-1-j]];//(N_PAIRS*MAX_SEQ_LEN)-1-i]];
		}
	}
	for (i = 0; i < tlen; ++i){
		tseq[i+MAX_SEQ_LEN] = seq_nt4_table[(uint8_t)tseq_[i]];//(N_PAIRS*MAX_SEQ_LEN)-1-i]];
	}
	flag &= KSW_EZ_SCORE_ONLY; //we are only computing the score
	// flag |= KSW_EZ_GENERIC_SC; //we are using the scoring matrix
	ksw_extz2_hw(MAX_SEQ_LEN, (input_datatype*)qseq.data(), MAX_SEQ_LEN, (input_datatype*)tseq.data(), m, q, e, w, zdrop, 0, flag, v_s, N_PAIRS, binary_file);

	// free(qseq); free(tseq);
}

static void print_aln(const char *tname, const char *qname, ksw_extz_t *ez)
{
	printf("%s\t%s\t%d", tname, qname, ez->score);
	printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
	if (ez->n_cigar > 0) {
		printf("\t%d", ez->n_cigar);//{
		int i;
		putchar('\t');
		for (i = 0; i < ez->n_cigar; ++i)
			printf("%d%c", ez->cigar[i]>>4, "MID"[ez->cigar[i]&0xf]);
	}
	putchar('\n');
	
}

/////////hw
static void print_aln_hw(const char *tname, const char *qname, ksw_hw_extz_t *ez)
{
	printf("%s\t%s\t%d", tname, qname, ez->score);
	printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
	if (ez->n_cigar > 0){
		printf("\t%d", ez->n_cigar);
		int i;
		putchar('\t');
//		for (i = 0; i < ez->n_cigar; ++i)
//			printf("%d%c", ez->cigar[i]>>4, "MID"[ez->cigar[i]&0xf]);
	}

	printf("\n");
}

static int controllo(int* comp_s, int* comp_t, int core_id){

	int flag=0;
	for(int i=0; i<N_PAIRS; i++){

		if(comp_s[i]!=comp_t[i]){
			printf("HW: %d, SW: %d, %d\n", comp_s[i], comp_t[i], i);
			flag=1;
		}

	}

	if(flag==1)	printf("\nERROR ON CORE %d\n\n", core_id);

	return flag;
}

static void random_generator(char* query, char* target, int qlen, int tlen)
{
	char alphabet[4] = {'A', 'C', 'G', 'T'};
	int i=0;

	for(i=0; i<qlen; i++)	query[i] = alphabet[rand()%4];
	for(i=0; i<tlen; i++)	target[i] = alphabet[rand()%4];

//	for(i=0; i<qlen; i++)	printf("%c",query[i]);
//	putchar('\n');
//	for(i=0; i<tlen; i++)	printf("%c", target[i]);
//	putchar('\n');

}

static void random_generator2(char* query, char* target) //con qlen e tlen più lunghe di MAX_SEQ_LEN
{
	char alphabet[4] = {'A', 'C', 'G', 'T'};
	int i=0, j=0;

	for(i=0; i<(MAX_SEQ_LEN*N_PAIRS); i++)	query[i] = alphabet[rand()%4];

	for(i=0; i<(MAX_SEQ_LEN*N_PAIRS); i++)	target[i] = alphabet[rand()%4];

//	for(i=0; i<(MAX_SEQ_LEN*N_PAIRS); i++)	printf("%c",query[i]);
//	printf("\n");
//	for(i=0; i<(MAX_SEQ_LEN*N_PAIRS); i++)	printf("%c", target[i]);
//	printf("\n\n");

}


#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#define EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((a), (m)); \
	} while (0)

typedef struct {
	char *name;
	char *seq;
} named_seq_t;

//=======================================================================
//
// FPGA function calls 
//
//=======================================================================

void ksw_extz2_hw(int qlen, input_datatype *query_l, int tlen, input_datatype *target_l, datatype m,
				 datatype q, datatype e, int w, int zdrop, int end_bonus, int flag, std::vector<int, aligned_allocator<int>> *v_s, int n_pairs, std::string binary_file)
{

	cl_int err;
    cl::CommandQueue queue;
    std::string krnl_name = "ksw_extz2_hw";
	std::vector<cl::Kernel> krnls(NUM_KERNEL); 
	cl::Context context;
	auto devices = xcl::get_xil_devices();

    // read_binary_file() command will find the OpenCL binary file created using the
    // V++ compiler load into OpenCL Binary and return pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binary_file);

    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    int valid_device = 0;
    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));
        OCL_CHECK(err,
                queue = cl::CommandQueue(context,
                                     device,
                                     CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE |
                                     CL_QUEUE_PROFILING_ENABLE,
                                     &err));

        std::cout << "Trying to program device[" << i
                         << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, NULL, &err);
        if (err != CL_SUCCESS) {
             std::cout << "Failed to program device[" << i
                     << "] with xclbin file!\n";
        } else {
             std::cout << "Device[" << i << "]: program successful!\n";
           //  Creating Kernel object using Compute unit names

            for (int i = 0; i < NUM_KERNEL; i++) {
                std::string cu_id = std::to_string(i + 1);
                std::string krnl_name_full =
                        krnl_name + ":{" + "ksw_extz2_hw_" + cu_id + "}";

             printf("Creating a kernel [%s] for CU(%d)\n",
                              krnl_name_full.c_str(),
                              i + 1);

                //Here Kernel object is created by specifying kernel name along with compute unit.
                //For such case, this kernel object can only access the specific Compute unit

                OCL_CHECK(err,
                krnls[i] = cl::Kernel(
                program, krnl_name_full.c_str(), &err));
            }
            valid_device++;
            break; // we break because we found a valid device
        }
    }
    if (valid_device == 0) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

	//SETUP DATA TRANSFER
	std::vector<cl_mem_ext_ptr_t> inBufExt1(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> inBufExt2(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> outBufExt(NUM_KERNEL);

    std::vector<cl::Buffer> buffer_input1(NUM_KERNEL);
    std::vector<cl::Buffer> buffer_input2(NUM_KERNEL);
    std::vector<cl::Buffer> buffer_output(NUM_KERNEL);

    for (int i = 0; i < NUM_KERNEL; i++) {

	        inBufExt1[i].obj = query_l; 
	        inBufExt1[i].param = 0;
	        inBufExt1[i].flags = bank[i * 3];

	        inBufExt2[i].obj = target_l;
	        inBufExt2[i].param = 0;
	        inBufExt2[i].flags = bank[(i * 3) + 1];

	       	outBufExt[i].obj = v_s[i].data();        //score
	        outBufExt[i].param = 0;
	        outBufExt[i].flags = bank[(i * 3) + 2];

    }

	// These commands will allocate memory on the FPGA. The cl::Buffer objects can
    // be used to reference the memory locations on the device.
    //Creating Buffers
    for (int i = 0; i < NUM_KERNEL; i++) {
	        OCL_CHECK(err,
	                    buffer_input1[i] =            
	                    cl::Buffer(context,
	                    CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
	                    CL_MEM_USE_HOST_PTR,
	                    (N_PAIRS+1)*MAX_SEQ_LEN,
	                    &inBufExt1[i],
	                    &err));
	      
		        OCL_CHECK(err,
		                    buffer_input2[i] =
		                    cl::Buffer(context,
		                    CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
		                    CL_MEM_USE_HOST_PTR,
		                    (N_PAIRS+1)*MAX_SEQ_LEN, 
		                    &inBufExt2[i],
		                    &err));
		        OCL_CHECK(err,
		                    buffer_output[i] =
		                    cl::Buffer(context,
		                    CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX |
		                    CL_MEM_USE_HOST_PTR,
		                    N_PAIRS*sizeof(int), 
		                    &outBufExt[i],
		                    &err));
    }

    // Copy input data to Device Global Memory
    for (int i = 0; i < NUM_KERNEL; i++) {
        OCL_CHECK(err,
        err = queue.enqueueMigrateMemObjects(
        {buffer_input1[i], buffer_input2[i]},
        0 /* 0 means from host*/));
    }
    queue.finish();
    double kernel_time_in_sec = 0, result = 0;
    std::chrono::duration<double> kernel_time(0);
	auto kernel_start = NOW;

    for (int i = 0; i < NUM_KERNEL; i++) {
    	//Setting kernel arguments 
	
		OCL_CHECK(err, err = krnls[i].setArg(0, qlen));   
		OCL_CHECK(err, err = krnls[i].setArg(1, buffer_input1[i])); 
		OCL_CHECK(err, err = krnls[i].setArg(2, tlen));    
		OCL_CHECK(err, err = krnls[i].setArg(3, buffer_input2[i]));
		OCL_CHECK(err, err = krnls[i].setArg(4, m));
		OCL_CHECK(err, err = krnls[i].setArg(5, q));  
		OCL_CHECK(err, err = krnls[i].setArg(6, e));
		OCL_CHECK(err, err = krnls[i].setArg(7, w)); 
		OCL_CHECK(err, err = krnls[i].setArg(8, zdrop));
		OCL_CHECK(err, err = krnls[i].setArg(9, end_bonus));
		OCL_CHECK(err, err = krnls[i].setArg(10, flag));
		OCL_CHECK(err, err = krnls[i].setArg(11, buffer_output[i]));
		OCL_CHECK(err, err = krnls[i].setArg(12, n_pairs));
		
		//void ksw_extz2_hw(int qlen, input_datatype *query_l, int tlen, input_datatype *target_l, datatype m,
		// datatype q, datatype e, int w, int zdrop, int end_bonus, int flag, int *v_s, int n_pairs)


		//Invoking the kernel
        OCL_CHECK(err, err = queue.enqueueTask(krnls[i]));	
	}

	queue.finish();
    
    auto kernel_end = NOW;	 
	kernel_time = std::chrono::duration<double>(kernel_end - kernel_start);
    kernel_time_in_sec = kernel_time.count();   
    std::cout << "\nFPGA TIME: " << kernel_time_in_sec << "s" << std::endl;      

    for (int i = 0; i < NUM_KERNEL; i++) {
        OCL_CHECK(err,
            err = queue.enqueueMigrateMemObjects(
            {buffer_output[i]},
            CL_MIGRATE_MEM_OBJECT_HOST));
    }
    queue.finish();
	
	/////////////////////////////////////////////////////////////////////////
	// printf("%d\n", v_s[0]);
	// //host mem reset
	// free(scoreLeft);
	// free(scoreRight);

}


//=======================================================================
//
// Function call main  
//
//=======================================================================


int main(int argc, char *argv[])
{

	std::string binary_file = argv[1];

	int8_t a = 2, b = 4, q = 4, e = 2, q2 = 13, e2 = 1;
	int c, i, pair = 1, w = -1, flag = 0, rep = 1, zdrop = -1, no_kalloc = 0;
	char *s;
	int8_t mat[25];
	ksw_extz_t ez;
	ksw_hw_extz_t ez_hw;

	//comparazione
	int comp_t[4];
	int comp_s[4];
	
	//generazione random
	// int seed = 1644417473;
	int seed = time(NULL);
	printf("Seed: %d\n", seed);
	srand(seed);

	int random_qlen, random_tlen;
    random_qlen = rand()%1280; 
	random_tlen = rand()%1280; 
	
	char random_query[MAX_SEQ_LEN*N_PAIRS];
	char random_target[MAX_SEQ_LEN*N_PAIRS];

	char query_test[MAX_SEQ_LEN];
	char target_test[MAX_SEQ_LEN];

	random_generator2(random_query, random_target);

	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(&ez_hw, 0, sizeof(ksw_hw_extz_t));
	ksw_gen_simple_mat(5, mat, a, -b);

	char *qpoint, *tpoint;

    qpoint = &random_query[0];
    tpoint = &random_target[0];

    // ksw_extz_t ris_test[N_PAIRS];
    // ksw_hw_extz_t ris_src[N_PAIRS];

    // int v_s[N_PAIRS];
	std::vector<int, aligned_allocator<int>> v_s[NUM_KERNEL];
	for(int i=0; i < NUM_KERNEL; i++)
		v_s[i].resize(N_PAIRS);
	std::vector<int, aligned_allocator<int>> v_t(N_PAIRS);
    // int v_t[N_PAIRS];

    copy_ez(&ez_hw,&ez);
    ez_hw.score=100;
	

    global_aln_hw(random_query, random_target, 5, mat, q, e, q2, e2, w, zdrop, flag, &ez_hw, v_s, binary_file);

	double sw_time_in_sec = 0, result = 0;
   	std::chrono::duration<double> sw_time(0);

	auto sw_start=NOW;
	for(int j=0; j<N_PAIRS; j++){

        memcpy(query_test, qpoint, sizeof(char)*MAX_SEQ_LEN);
        memcpy(target_test, tpoint, sizeof(char)*MAX_SEQ_LEN);

	/*	//stampa delle sequenze
		for(i=0; i<(MAX_SEQ_LEN); i++)	printf("%c", query_test[i]);
		printf("\n");
		for(i=0; i<(MAX_SEQ_LEN); i++)	printf("%c", target_test[i]);
		printf("\n\n");
*/

		global_aln_sw(query_test, target_test, 5, mat, q, e, q2, e2, w, zdrop, flag, &ez);

		// ris_test[j].max_q = ez.max_q;
		// ris_test[j].max_t = ez.max_t;
		// ris_test[j].mqe_t = ez.mqe_t;
		// ris_test[j].mte_q = ez.mte_q;
		// ris_test[j].max = ez.max;
		// ris_test[j].score = ez.score;
		// ris_test[j].mqe = ez.mqe;
		// ris_test[j].mte = ez.mte;
		// ris_test[j].m_cigar = ez.m_cigar;
		// ris_test[j].n_cigar = ez.n_cigar;
		// ris_test[j].zdropped = ez.zdropped;
		// ris_test[j].reach_end = ez.reach_end;

		//print_aln("test", "", &ez);

		v_t[j]=ez.score;

//		printf("score test: %d\n", v_t[j]);

		//incremento per il prossimo set di vettori
		qpoint = qpoint + (sizeof(char)*MAX_SEQ_LEN);
		tpoint = tpoint + (sizeof(char)*MAX_SEQ_LEN);
	
	}
	auto sw_end = NOW;
        sw_time = std::chrono::duration<double>(sw_end - sw_start);
    	sw_time_in_sec = sw_time.count();
	printf("SW TIME (single thread): %.5lfs\n", sw_time_in_sec*10);//the FPGA is doing 10 times more alignments than the CPU, as right now it is setup to run the same N_PAIRS alignments on a single core
	printf("SW TIME (multi thread): %.5lfs\n", sw_time_in_sec*10/16);
	int error = 0;
	for(int i=0; i < NUM_KERNEL; i++)
		error|=controllo(v_s[i].data(),v_t.data(),i);
	if(!error)
		printf("\nALL RESULTS OK\n");
	return 0;
}

