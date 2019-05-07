/**
 * @file    no_partitioning_join.c
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Sun Feb  5 20:16:58 2012
 * @version $Id: no_partitioning_join.c 3017 2012-12-07 10:56:20Z bcagri $
 * 
 * @brief  The implementation of NPO, No Partitioning Optimized join algortihm.
 * 
 * (c) 2012, ETH Zurich, Systems Group
 * 
 */


#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sched.h>              /* CPU_ZERO, CPU_SET */
#include <pthread.h>            /* pthread_* */
#include <string.h>             /* memset */
#include <stdio.h>              /* printf */
#include <stdlib.h>             /* memalign */
#include <sys/time.h>           /* gettimeofday */
#include <immintrin.h>

#include "nest_loop_join.h"
#include "npj_params.h"         /* constant parameters */
#include "npj_types.h"          /* bucket_t, hashtable_t, bucket_buffer_t */
#include "rdtsc.h"              /* startTimer, stopTimer */
#include "lock.h"               /* lock, unlock */
#include "cpu_mapping.h"        /* get_cpu_id */
#ifdef PERF_COUNTERS
#include "perf_counters.h"      /* PCM_x */
#endif

#include "barrier.h"            /* pthread_barrier_* */
#include "affinity.h"           /* pthread_attr_setaffinity_np */
#include "generator.h"          /* numa_localize() */

#ifndef BARRIER_ARRIVE
/** barrier wait macro */
#define BARRIER_ARRIVE(B,RV)                            \
    RV = pthread_barrier_wait(B);                       \
    if(RV !=0 && RV != PTHREAD_BARRIER_SERIAL_THREAD){  \
        printf("Couldn't wait on barrier\n");           \
        exit(EXIT_FAILURE);                             \
    }
#endif

#define L1 (1 << 15)    /* Working set size for L1 cache 32KB */
#define L2 (1 << 18)    /* Working set size for L2 cache 256KB */
#define L3 ((1 << 20) * 5 / 2)    /* Working set size for L3 cache 2.5MB */
#define LLC ((1 << 20) * 55)    /* Working set size for LLC cache 55MB */
#define MAXELEMS 6000
#define random(x) (rand()%(x))

#define RAND_RANGE(N) ((double)rand() / ((double)RAND_MAX + 1) * (N))
/** Debug msg logging method */
#ifdef DEBUG
#define DEBUGMSG(COND, MSG, ...)                                    \
    if(COND) { fprintf(stdout, "[DEBUG] "MSG, ## __VA_ARGS__); }
#else
#define DEBUGMSG(COND, MSG, ...) 
#endif
typedef struct arg_nlj arg_nlj;

struct arg_nlj {
    int32_t             tid;
    relation_t          relR;
    relation_t          relS;
    pthread_barrier_t * barrier;
    int64_t             num_results;
#ifndef NO_TIMING
    /* stats about the thread */
    uint64_t timer1;
    struct timeval start, end;
#endif
} ;

/** print out the execution time statistics of the join */
static void 
print_timing(uint64_t total, uint64_t numtuples, int64_t result,
            struct timeval * start, struct timeval * end)
{
    double diff_usec = (((*end).tv_sec*1000000L + (*end).tv_usec)
                        - ((*start).tv_sec*1000000L+(*start).tv_usec));
    double cyclestuple = total;
    cyclestuple /= numtuples;
    fprintf(stdout, "TOTAL-TUPLES  ,RESULT-TUPLES ,RUNTIME TOTAL ,TOTAL-TIME-USECS,  CYCLES-PER-TUPLE: \n");
    fprintf(stderr, "%-15llu%-15llu%-15llu%11.4lf    %11.4lf", 
            numtuples, result, total, diff_usec, cyclestuple);
    fflush(stdout);
    fflush(stderr);
    fprintf(stdout, "\n");
}

void *
snlj_thread(void *param)
{
    int rv;
    arg_nlj * args = (arg_nlj*) param;

    /* wait at a barrier until each thread starts and start timer */
    BARRIER_ARRIVE(args->barrier, rv);

#ifndef NO_TIMING
    /* the first thread checkpoints the start time */
    if(args->tid == 0){
        gettimeofday(&args->start, NULL);
        startTimer(&args->timer1);
    }
#endif

    int64_t result = 0;
    __m256i yidOuter;
    __m256i yidInner;
    __m256i yidResult = _mm256_set1_epi32(0);
    __m256i yidTmp;
    uint32_t vec_len = args->relS.num_tuples * 2;
    uint32_t n = args->relR.num_tuples;

    int i, j, k;
    int q[8];
    int partition = L1 / 2 / 4;
    void *p;
    for (k = 0; k < vec_len; k += partition) {
	int end = k + partition < vec_len ? k + partition : vec_len;
	int bound = (end - k) / 8 * 8 + k;
	for (i = 0; i < n; ++i) {
	    yidOuter = _mm256_set1_epi32(args->relR.tuples[i].key);
	    p = args->relS.tuples + k;
	    for (j = k; j < bound; j += 8) {
		yidInner = _mm256_loadu_si256(p);
	        yidTmp = _mm256_cmpeq_epi32(yidOuter, yidInner);
		yidResult = _mm256_sub_epi32(yidResult, yidTmp);
		p += 32;
	    }
	    for (; j < end; ++j) {
		if (args->relR.tuples[i].key == args->relS.tuples[j].key) {
	            ++result;
		}
	    }
	}
    }
    _mm256_storeu_si256((void *)q, yidResult);
    //result += q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
    result += q[0] + q[2] + q[4] + q[6];
    args->num_results = result;
	// printf("[%d]", sink);

#ifndef NO_TIMING
    /* for a reliable timing we have to wait until all finishes */
    BARRIER_ARRIVE(args->barrier, rv);

    if(args->tid == 0){
      stopTimer(&args->timer1); 
      gettimeofday(&args->end, NULL);
    }
#endif

    return 0;
}

int64_t 
SNLJ(relation_t *relR, relation_t *relS, int nthreads)
{
    int64_t result = 0;
    int32_t numR, numS, numRthr, numSthr; /* total and per thread num */
    int i, rv;
    cpu_set_t set;
    arg_nlj args[nthreads];
    pthread_t tid[nthreads];
    pthread_attr_t attr;
    pthread_barrier_t barrier;

    numR = relR->num_tuples;
    numS = relS->num_tuples;
    numRthr = numR / nthreads;
    
    rv = pthread_barrier_init(&barrier, NULL, nthreads);
    if(rv != 0) {
        printf("Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);
    for(i = 0; i < nthreads; i++) {
        int cpu_idx = get_cpu_id(i);

        DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].tid = i;
        args[i].barrier = &barrier;

        /* assing part of the relR for next thread */
        args[i].relR.num_tuples = (i == (nthreads-1)) ? numR : numRthr;
        args[i].relR.tuples = relR->tuples + numRthr * i;
        numR -= numRthr;

        /* assing part of the relS for next thread */
        args[i].relS.num_tuples = numS;
        args[i].relS.tuples = relS->tuples;

        rv = pthread_create(&tid[i], &attr, snlj_thread, (void*)&args[i]);
        if (rv){
            printf("ERROR; return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }

    for(i = 0; i < nthreads; i++) {
        pthread_join(tid[i], NULL);
        /* sum up results */
        result += args[i].num_results;
    }


#ifndef NO_TIMING
    /* now print the timing results: */
    print_timing(args[0].timer1, relS->num_tuples, result, &args[0].start, &args[0].end);
#endif

    return result;
}

void *
nlj_thread(void *param)
{
    int rv;
    arg_nlj * args = (arg_nlj*) param;

    /* wait at a barrier until each thread starts and start timer */
    BARRIER_ARRIVE(args->barrier, rv);

#ifndef NO_TIMING
    /* the first thread checkpoints the start time */
    if(args->tid == 0){
        gettimeofday(&args->start, NULL);
        startTimer(&args->timer1);
    }
#endif

    int i, j;
    int64_t result = 0; 
    for (i = 0; i < args->relR.num_tuples; i++) {
        for (j = 0; j < args->relS.num_tuples; j++) {
            if (args->relR.tuples[i].key == args->relS.tuples[j].key) result++;
        }
    }
    args->num_results = result;

#ifndef NO_TIMING
    /* for a reliable timing we have to wait until all finishes */
    BARRIER_ARRIVE(args->barrier, rv);

    if(args->tid == 0){
      stopTimer(&args->timer1); 
      gettimeofday(&args->end, NULL);
    }
#endif

    return 0;
}

int64_t
NLJ(relation_t *relR, relation_t *relS, int nthreads)
{
    int64_t result = 0;
    int32_t numR, numS, numRthr, numSthr; /* total and per thread num */
    int i, rv;
    cpu_set_t set;
    arg_nlj args[nthreads];
    pthread_t tid[nthreads];
    pthread_attr_t attr;
    pthread_barrier_t barrier;

    numR = relR->num_tuples;
    numS = relS->num_tuples;
    numRthr = numR / nthreads;
    
    rv = pthread_barrier_init(&barrier, NULL, nthreads);
    if(rv != 0) {
        printf("Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);
    for(i = 0; i < nthreads; i++) {
        int cpu_idx = get_cpu_id(i);

        DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].tid = i;
        args[i].barrier = &barrier;

        /* assing part of the relR for next thread */
        args[i].relR.num_tuples = (i == (nthreads-1)) ? numR : numRthr;
        args[i].relR.tuples = relR->tuples + numRthr * i;
        numR -= numRthr;

        /* assing part of the relS for next thread */
        args[i].relS.num_tuples = numS;
        args[i].relS.tuples = relS->tuples;

        rv = pthread_create(&tid[i], &attr, nlj_thread, (void*)&args[i]);
        if (rv){
            printf("ERROR; return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }

    for(i = 0; i < nthreads; i++) {
        pthread_join(tid[i], NULL);
        /* sum up results */
        result += args[i].num_results;
    }


#ifndef NO_TIMING
    /* now print the timing results: */
    print_timing(args[0].timer1, relS->num_tuples, result, &args[0].start, &args[0].end);
#endif

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
    uint32_t vec_len = relS->num_tuples * 2;
    uint32_t n = relR->num_tuples;

    int i, j, k;
    int q[8];
    int partition = L1 / 2 / 4;
    void *p;
    for (k = 0; k < vec_len; k += partition) {
	int end = k + partition < vec_len ? k + partition : vec_len;
	int bound = (end - k) / 8 * 8 + k;
	for (i = 0; i < n; ++i) {
	    yidOuter = _mm256_set1_epi32(relR->tuples[i].key);
	    p = relS->tuples + k;
	    for (j = k; j < bound; j += 8) {
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
    //result += q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
    result += q[0] + q[2] + q[4] + q[6];
	// printf("[%d]", sink);

   return result;
}
