#include "nest_loop_join.h"
#include "x86intrin.h"
#include <immintrin.h>

#include <sched.h>              /* CPU_ZERO, CPU_SET */
#include <pthread.h>            /* pthread_* */
#include <stdlib.h>             /* malloc, posix_memalign */
#include <sys/time.h>           /* gettimeofday */
#include <stdio.h>              /* printf */

#define L1 (1 << 15)    /* Working set size for L1 cache 32KB */
#define L2 (1 << 18)    /* Working set size for L2 cache 256KB */
#define L3 ((1 << 20) * 5 / 2)    /* Working set size for L3 cache 2.5MB */
#define LLC ((1 << 20) * 55)    /* Working set size for LLC cache 55MB */
#define MAXELEMS 6000
#define random(x) (rand()%(x))

int64_t
NLJ(relation_t *relR, relation_t *relS, int nthreads)
{
    int i, j;
    int64_t result = 0; 
    for (i = 0; i < relR->num_tuples; i++) {
        for (j = 0; j < relS->num_tuples; j++) {
            if (relR->tuples[i].key == relS->tuples[j].key) result++;
        }
    }
    return result;
}

int64_t 
SNLJ(relation_t *relR, relation_t *relS, int nthreads)
{
   int64_t result = 0;

   return result;
}

int64_t 
SPNLJ(relation_t *relR, relation_t *relS, int nthreads)
{
   int64_t result = 0;
	__m256i yidOuter;
	__m256i yidInner;
	__m256i yidResult = _mm256_set1_epi32(0);
	__m256i yidTmp;
        uint32_t vec_len = relS->num_tuples;
        uint32_t n = relR->num_tuples;

	int i, j, k;
	int q[8];
	int partition = L1/ 2 / 4;
	void *p;
	for (k = 0; k < vec_len; k += partition) {
		int end = k + partition < vec_len ? k + partition : vec_len;
		int bound = (end - k) / 8 * 8 + k;
		for (i = 0; i < n; ++i) {
			yidOuter = _mm256_set1_epi32(relR->tuples[i].key);
                        int z = 0; 
                        intkey_t *tmpK = malloc(sizeof(intkey_t) * 8);
                        for (z = 0; z < 8; z++)
                        {
                            tmpK[z] = relR->tuples[z + k].key; 
                        }
                        p = tmpK;
			for (j = k; j < bound; j+= 8) {
				yidInner = _mm256_loadu_si256(p);
				yidTmp = _mm256_cmpeq_epi32(yidOuter, yidInner);
				yidResult = _mm256_sub_epi32(yidResult, yidTmp);
				p += 32;
			}
			for (; j < end; ++j) {
				if (relR->tuples[i].key == relS->tuples[j].key) {
					++result;
				}
			}
		}
	}
	_mm256_storeu_si256((void *)q, yidResult);
	result += q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
	// printf("[%d]", sink);
   return result;
}
