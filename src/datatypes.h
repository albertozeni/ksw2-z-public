#ifndef DATATYPES_H
#define DATATYPES_H

#define MAX_SEQ_LEN 256
#define DEPTH_STREAM 200
#define LOG_2_PORT_WIDTH 8
#define LOG_2_BITS_PER_CHAR 3
#define CHAR_PER_PACKAGE ((1<<LOG_2_PORT_WIDTH)>>3)
#define PACKAGES_PER_MAX_SEQ_LEN MAX_SEQ_LEN/CHAR_PER_PACKAGE
//#define MAX_SEQ_LEN	16384
#define M 8 //dimension of the scoring matrix is M*M

typedef uint8_t u_datatype;
typedef int8_t datatype;
typedef ap_uint<1<<LOG_2_PORT_WIDTH> input_datatype;

typedef struct { //
	uint32_t max;
	uint32_t zdropped;
	int max_q, max_t;      // max extension coordinate
	int mqe, mqe_t;        // max score when reaching the end of query
	int mte, mte_q;        // max score when reaching the end of target
	int score;             // max score reaching both ends; may be HW_KSW_NEG_INF
	int m_cigar, n_cigar;
	int reach_end;
	int pat[4];
//	uint32_t cigar[MAX_TLEN+MAX_QLEN];
} ksw_hw_extz_t;


#endif
