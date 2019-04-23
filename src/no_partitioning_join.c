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

#include "no_partitioning_join.h"
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

#ifndef NEXT_POW_2
/** 
 *  compute the next number, greater than or equal to 32-bit unsigned v.
 *  taken from "bit twiddling hacks":
 *  http://graphics.stanford.edu/~seander/bithacks.html
 */
#define NEXT_POW_2(V)                           \
    do {                                        \
        V--;                                    \
        V |= V >> 1;                            \
        V |= V >> 2;                            \
        V |= V >> 4;                            \
        V |= V >> 8;                            \
        V |= V >> 16;                           \
        V++;                                    \
    } while(0)
#endif

#ifndef HASH
#define HASH(X, MASK, SKIP) (((X) & MASK) >> SKIP)
#endif
#define RAND_RANGE(N) ((double)rand() / ((double)RAND_MAX + 1) * (N))
/** Debug msg logging method */
#ifdef DEBUG
#define DEBUGMSG(COND, MSG, ...)                                    \
    if(COND) { fprintf(stdout, "[DEBUG] "MSG, ## __VA_ARGS__); }
#else
#define DEBUGMSG(COND, MSG, ...) 
#endif

/** An experimental feature to allocate input relations numa-local */
extern int numalocalize;  /* defined in generator.c */
extern int nthreads;      /* defined in generator.c */

/**
 * \ingroup NPO arguments to the threads
 */
typedef struct arg_t arg_t;
typedef struct arg_bm arg_bm;
typedef struct arg_vec arg_vec;
typedef struct arg_sj arg_sj;

struct arg_t {
    int32_t             tid;
    hashtable_t *       ht;
    relation_t          relR;
    relation_t          relS;
    pthread_barrier_t * barrier;
    int64_t             num_results;
#ifndef NO_TIMING
    /* stats about the thread */
    uint64_t timer1, timer2, timer3;
    struct timeval start, end;
#endif
} ;
//bitmap for NPO_bm by zys
struct arg_bm {
    int32_t             tid;
    hashtable_t *       ht;
    int8_t *            bitmap;
    relation_t          relR;
    relation_t          relS;
    pthread_barrier_t * barrier;
    int64_t             num_results;
#ifndef NO_TIMING
    /* stats about the thread */
    uint64_t timer1, timer2, timer3;
    struct timeval start, end;
#endif
} ;
//vector struct by zys
struct arg_vec {
    int32_t             tid;
    intvector_t *       vec;
    double               update_ratio;
    relation_t          relR;
    relation_t          relS;
    pthread_barrier_t * barrier;
    int64_t             num_results;
#ifndef NO_TIMING
    /* stats about the thread */
    uint64_t timer1, timer2, timer3;
    struct timeval start, end;
#endif
} ;
//starjoin struct by zys
struct arg_sj {
    int32_t             tid;
    int32_t             fkid;
    column_t *          fks;
    vector_t *          pks;
    vectorkey_t *          MInx;
    int32_t             FKStartIndex;
    int32_t             MIStartIndex;
    pthread_barrier_t * barrier;
    int64_t             num_results;
#ifndef NO_TIMING
    /* stats about the thread */
    uint64_t timer1, timer2, timer3;
    struct timeval start, end;
#endif
} ;

/** 
 * @defgroup OverflowBuckets Buffer management for overflowing buckets.
 * Simple buffer management for overflow-buckets organized as a 
 * linked-list of bucket_buffer_t.
 * @{
 */

/** 
 * Initializes a new bucket_buffer_t for later use in allocating 
 * buckets when overflow occurs.
 * 
 * @param ppbuf [in,out] bucket buffer to be initialized
 */
void 
init_bucket_buffer(bucket_buffer_t ** ppbuf)
{
    bucket_buffer_t * overflowbuf;
    overflowbuf = (bucket_buffer_t*) malloc(sizeof(bucket_buffer_t));
    overflowbuf->count = 0;
    overflowbuf->next  = NULL;

    *ppbuf = overflowbuf;
}

/** 
 * Returns a new bucket_t from the given bucket_buffer_t.
 * If the bucket_buffer_t does not have enough space, then allocates
 * a new bucket_buffer_t and adds to the list.
 *
 * @param result [out] the new bucket
 * @param buf [in,out] the pointer to the bucket_buffer_t pointer
 */
inline void 
get_new_bucket(bucket_t ** result, bucket_buffer_t ** buf)
{
    if((*buf)->count < OVERFLOW_BUF_SIZE) {
        *result = (*buf)->buf + (*buf)->count;
        (*buf)->count ++;
    }
    else {
        /* need to allocate new buffer */
        bucket_buffer_t * new_buf = (bucket_buffer_t*) 
                                    malloc(sizeof(bucket_buffer_t));
        new_buf->count = 1;
        new_buf->next  = *buf;
        *buf    = new_buf;
        *result = new_buf->buf;
    }
}

/** De-allocates all the bucket_buffer_t */
void
free_bucket_buffer(bucket_buffer_t * buf)
{
    do {
        bucket_buffer_t * tmp = buf->next;
        free(buf);
        buf = tmp;
    } while(buf);
}

/** @} */


/** 
 * @defgroup NPO The No Partitioning Optimized Join Implementation
 * @{
 */

/** 
 * Allocates a hashtable of NUM_BUCKETS and inits everything to 0. 
 * 
 * @param ht pointer to a hashtable_t pointer
 */
void 
allocate_hashtable(hashtable_t ** ppht, uint32_t nbuckets)
{
    hashtable_t * ht;

    ht              = (hashtable_t*)malloc(sizeof(hashtable_t));
    ht->num_buckets = nbuckets;
    NEXT_POW_2((ht->num_buckets));

    /* allocate hashtable buckets cache line aligned */
    if (posix_memalign((void**)&ht->buckets, CACHE_LINE_SIZE,
                       ht->num_buckets * sizeof(bucket_t))){
        perror("Aligned allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    /** Not an elegant way of passing whether we will numa-localize, but this
        feature is experimental anyway. */
    if(numalocalize) {
        tuple_t * mem = (tuple_t *) ht->buckets;
        uint32_t ntuples = (ht->num_buckets*sizeof(bucket_t))/sizeof(tuple_t);
        numa_localize(mem, ntuples, nthreads);
    }

    memset(ht->buckets, 0, ht->num_buckets * sizeof(bucket_t));
    ht->skip_bits = 0; /* the default for modulo hash */
    ht->hash_mask = (ht->num_buckets - 1) << ht->skip_bits;
    *ppht = ht;
}

/** 
 * Releases memory allocated for the hashtable.
 * 
 * @param ht pointer to hashtable
 */
void 
destroy_hashtable(hashtable_t * ht)
{
    free(ht->buckets);
    free(ht);
}

/** 
 * Single-thread hashtable build method, ht is pre-allocated.
 * 
 * @param ht hastable to be built
 * @param rel the build relation
 */
void 
build_hashtable_st(hashtable_t *ht, relation_t *rel)
{
    uint32_t i;
    const uint32_t hashmask = ht->hash_mask;
    const uint32_t skipbits = ht->skip_bits;

    for(i=0; i < rel->num_tuples; i++){
        tuple_t * dest;
        bucket_t * curr, * nxt;
        int32_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);

        /* copy the tuple to appropriate hash bucket */
        /* if full, follow nxt pointer to find correct place */
        curr = ht->buckets + idx;
        nxt  = curr->next;

        if(curr->count == BUCKET_SIZE) {
            if(!nxt || nxt->count == BUCKET_SIZE) {
                bucket_t * b;
                b = (bucket_t*) calloc(1, sizeof(bucket_t));
                curr->next = b;
                b->next = nxt;
                b->count = 1;
                dest = b->tuples;
            }
            else {
                dest = nxt->tuples + nxt->count;
                nxt->count ++;
            }
        }
        else {
            dest = curr->tuples + curr->count;
            curr->count ++;
        }
        *dest = rel->tuples[i];
    }
}

/** 
 * Probes the hashtable for the given outer relation, returns num results. 
 * This probing method is used for both single and multi-threaded version.
 * 
 * @param ht hashtable to be probed
 * @param rel the probing outer relation
 * 
 * @return number of matching tuples
 */
int64_t 
probe_hashtable(hashtable_t *ht, relation_t *rel)
{
    uint32_t i, j;
    int64_t matches;

    const uint32_t hashmask = ht->hash_mask;
    const uint32_t skipbits = ht->skip_bits;
#ifdef PREFETCH_NPJ    
    size_t prefetch_index = PREFETCH_DISTANCE;
#endif
    
    matches = 0;

    for (i = 0; i < rel->num_tuples; i++)
    {
#ifdef PREFETCH_NPJ        
        if (prefetch_index < rel->num_tuples) {
			intkey_t idx_prefetch = HASH(rel->tuples[prefetch_index++].key,
                                         hashmask, skipbits);
			__builtin_prefetch(ht->buckets + idx_prefetch, 0, 1);
        }
#endif
        
        intkey_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
        bucket_t * b = ht->buckets+idx;

        do {
            for(j = 0; j < b->count; j++) {
                if(rel->tuples[i].key == b->tuples[j].key){
                    matches ++;
                    /* TODO: we don't materialize the results. */
                }
            }

            b = b->next;/* follow overflow pointer */
        } while(b);
    }

    return matches;
}
int64_t 
probe_bmhashtable(hashtable_t *ht, int8_t * bm, relation_t *rel)
{
    uint32_t i, j;
    int64_t matches,counter;

    const uint32_t hashmask = ht->hash_mask;
    const uint32_t skipbits = ht->skip_bits;
#ifdef PREFETCH_NPJ    
    size_t prefetch_index = PREFETCH_DISTANCE;
#endif
    
    matches = 0;counter=0;
    //printf("\nprob_bmhashtable test:%d %d",bm[0],bm[1]);

    for (i = 0; i < rel->num_tuples; i++)
    {
      if (bm[i]==1) {  //--adding bitmap filtering by zys
#ifdef PREFETCH_NPJ        
        if (prefetch_index < rel->num_tuples) {
			intkey_t idx_prefetch = HASH(rel->tuples[prefetch_index++].key,
                                         hashmask, skipbits);
			__builtin_prefetch(ht->buckets + idx_prefetch, 0, 1);
        }
#endif
        
        intkey_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
        bucket_t * b = ht->buckets+idx;

        do {
            for(j = 0; j < b->count; j++) {
                if(rel->tuples[i].key == b->tuples[j].key){
                    matches ++;
                    /* TODO: we don't materialize the results. */
                }
            }

            b = b->next;/* follow overflow pointer */
        } while(b);
      }//--end for bitmap filtering by zys
    }
    //printf("\nbitmap counter:%d",counter);
    return matches;
}

int64_t 
probe_vector(intvector_t *vec, relation_t *rel)
{
    uint32_t i;
    int64_t matches;

    //const uint32_t hashmask = ht->hash_mask;
    //const uint32_t skipbits = ht->skip_bits;
    
    matches = 0;

    for (i = 0; i < rel->num_tuples; i++)
    {
        intkey_t idx = rel->tuples[i].key;
        if(vec[idx-1]== 1){//if(rel->tuples[i].key == vec[idx-1]){  //--get predicate vector value
                    matches ++; 
         }
        //else {printf("%d \t %d \t ", idx, vec[idx-1]);}
        
     }

    return matches;
}

int64_t 
probe_uvector(intvector_t *vec, relation_t *rel,double uratio)
{
    uint32_t i;
    int64_t matches;

    //const uint32_t hashmask = ht->hash_mask;
    //const uint32_t skipbits = ht->skip_bits;
    
    matches = 0;

    for (i = 0; i < rel->num_tuples; i++)
    {
        intkey_t idx = rel->tuples[i].key;
        if(vec[idx-1]!=0 ){ rel->tuples[i].key=vec[idx-1];//--update FK value with updated PK value
                    matches ++; 
         }
        //else {printf("%d \t %d \t ", idx, vec[idx-1]);}
        
     }

    return matches;
}


/** print out the execution time statistics of the join */
static void 
print_timing(uint64_t total, uint64_t build, uint64_t part,
            uint64_t numtuples, int64_t result,
            struct timeval * start, struct timeval * end)
{
    double diff_usec = (((*end).tv_sec*1000000L + (*end).tv_usec)
                        - ((*start).tv_sec*1000000L+(*start).tv_usec));
    double cyclestuple = total;
    cyclestuple /= numtuples;
    fprintf(stdout, "TOTAL-TUPLES  ,RESULT-TUPLES ,RUNTIME TOTAL ,BUILD        ,PART            ,PROBE        ,TOTAL-TIME-USECS,  CYCLES-PER-TUPLE: \n");
    fprintf(stderr, "%-15llu%-15llu%-15llu%-15llu%-15llu%-15llu%11.4lf    %11.4lf", 
            numtuples,result, total, build, part, (total-build-part), diff_usec, cyclestuple);
    fflush(stdout);
    fflush(stderr);
    fprintf(stdout, "\n");

}

/** \copydoc NPO_st */
int64_t 
NPO_st(relation_t *relR, relation_t *relS, int nthreads)
{
    hashtable_t * ht;
    int64_t result = 0;
#ifndef NO_TIMING
    struct timeval start, end;
    uint64_t timer1, timer2, timer3;
#endif
    uint32_t nbuckets = (relR->num_tuples / BUCKET_SIZE);
    allocate_hashtable(&ht, nbuckets);

#ifndef NO_TIMING
    gettimeofday(&start, NULL);
    startTimer(&timer1);
    startTimer(&timer2); 
    timer3 = 0; /* no partitioning */
#endif

    build_hashtable_st(ht, relR);

#ifndef NO_TIMING
    stopTimer(&timer2); /* for build */
#endif

    result = probe_hashtable(ht, relS);

#ifndef NO_TIMING
    stopTimer(&timer1); /* over all */
    gettimeofday(&end, NULL);
    /* now print the timing results: */
    print_timing(timer1, timer2, timer3, relS->num_tuples, result, &start, &end);
#endif

    destroy_hashtable(ht);

    return result;
}

/** 
 * Multi-thread hashtable build method, ht is pre-allocated.
 * Writes to buckets are synchronized via latches.
 *
 * @param ht hastable to be built
 * @param rel the build relation
 * @param overflowbuf pre-allocated chunk of buckets for overflow use.
 */
void 
build_hashtable_mt(hashtable_t *ht, relation_t *rel, 
                   bucket_buffer_t ** overflowbuf)
{
    uint32_t i;
    const uint32_t hashmask = ht->hash_mask;
    const uint32_t skipbits = ht->skip_bits;

#ifdef PREFETCH_NPJ
    size_t prefetch_index = PREFETCH_DISTANCE;
#endif
    
    for(i=0; i < rel->num_tuples; i++){
        tuple_t * dest;
        bucket_t * curr, * nxt;

#ifdef PREFETCH_NPJ
        if (prefetch_index < rel->num_tuples) {
            intkey_t idx_prefetch = HASH(rel->tuples[prefetch_index++].key,
                                         hashmask, skipbits);
			__builtin_prefetch(ht->buckets + idx_prefetch, 1, 1);
        }
#endif
        
        int32_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
        /* copy the tuple to appropriate hash bucket */
        /* if full, follow nxt pointer to find correct place */
        curr = ht->buckets+idx;
        lock(&curr->latch);
        nxt = curr->next;

        if(curr->count == BUCKET_SIZE) {
            if(!nxt || nxt->count == BUCKET_SIZE) {
                bucket_t * b;
                /* b = (bucket_t*) calloc(1, sizeof(bucket_t)); */
                /* instead of calloc() everytime, we pre-allocate */
                get_new_bucket(&b, overflowbuf);
                curr->next = b;
                b->next    = nxt;
                b->count   = 1;
                dest       = b->tuples;
            }
            else {
                dest = nxt->tuples + nxt->count;
                nxt->count ++;
            }
        }
        else {
            dest = curr->tuples + curr->count;
            curr->count ++;
        }

        *dest = rel->tuples[i];
        unlock(&curr->latch);
    }

}

void 
build_vector_mt(intvector_t *vec, relation_t *rel)
{
    uint32_t i;
    for(i=0;i<rel->num_tuples;i++){
        vec[rel->tuples[i].key-1]=(rel->tuples[i].payload==rel->tuples[i].key); 
        //--simulating predicate processing and assigning predicate result for vector by zys
        //printf("%d ", vec[i]);
    }
}

void  //--simulating updates on PK table and cascading updates on FK table by zys
build_uvector_mt(intvector_t *vec, relation_t *rel, double uratio)
{
    uint32_t i;
    for(i=0;i<rel->num_tuples;i++){
        if (RAND_RANGE(rel->num_tuples)<uratio*rel->num_tuples) {
               vec[rel->tuples[i].key-1]=RAND_RANGE(i);
               }  
             //--PK column values under update ratio are set to random positions representing moving other tuple to current positions
        else  {
                vec[rel->tuples[i].key-1]=0;
              } //--PK column values beyond update ration are set 0 to represent no updating
    }
}

/** 
 * Just a wrapper to call the build and probe for each thread.
 * 
 * @param param the parameters of the thread, i.e. tid, ht, reln, ...
 * 
 * @return 
 */
void * 
npo_thread(void * param)
{
    int rv;
    arg_t * args = (arg_t*) param;

    /* allocate overflow buffer for each thread */
    bucket_buffer_t * overflowbuf;  //--no bucket by zys
    init_bucket_buffer(&overflowbuf);  //--no bucket by zys

#ifdef PERF_COUNTERS
    if(args->tid == 0){
        PCM_initPerformanceMonitor(NULL, NULL);
        PCM_start();
    }
#endif
    
    /* wait at a barrier until each thread starts and start timer */
    BARRIER_ARRIVE(args->barrier, rv);

#ifndef NO_TIMING
    /* the first thread checkpoints the start time */
    if(args->tid == 0){
        gettimeofday(&args->start, NULL);
        startTimer(&args->timer1);
        startTimer(&args->timer2); 
        args->timer3 = 0; /* no partitionig phase */
    }
#endif

    /* insert tuples from the assigned part of relR to the ht */
    build_hashtable_mt(args->ht, &args->relR, &overflowbuf);
     /* wait at a barrier until each thread completes build phase */
    BARRIER_ARRIVE(args->barrier, rv);

#ifdef PERF_COUNTERS
    if(args->tid == 0){
      PCM_stop();
      PCM_log("========== Build phase profiling results ==========\n");
      PCM_printResults();
      PCM_start();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif


#ifndef NO_TIMING
    /* build phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
        stopTimer(&args->timer2); 
    }
#endif

    /* probe for matching tuples from the assigned part of relS */
    args->num_results = probe_hashtable(args->ht, &args->relS);

#ifndef NO_TIMING
    /* for a reliable timing we have to wait until all finishes */
    BARRIER_ARRIVE(args->barrier, rv);

    /* probe phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
      stopTimer(&args->timer1); 
      gettimeofday(&args->end, NULL);
    }
#endif

#ifdef PERF_COUNTERS
    if(args->tid == 0) {
        PCM_stop();
        PCM_log("========== Probe phase profiling results ==========\n");
        PCM_printResults();
        PCM_log("===================================================\n");
        PCM_cleanup();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif

    /* clean-up the overflow buffers */
    free_bucket_buffer(overflowbuf);

    return 0;
}
void * 
npo_bm_thread(void * param)
{
    int rv;
    arg_bm * args = (arg_bm*) param;

    /* allocate overflow buffer for each thread */
    bucket_buffer_t * overflowbuf;  //--no bucket by zys
    init_bucket_buffer(&overflowbuf);  //--no bucket by zys

#ifdef PERF_COUNTERS
    if(args->tid == 0){
        PCM_initPerformanceMonitor(NULL, NULL);
        PCM_start();
    }
#endif
    
    /* wait at a barrier until each thread starts and start timer */
    BARRIER_ARRIVE(args->barrier, rv);

#ifndef NO_TIMING
    /* the first thread checkpoints the start time */
    if(args->tid == 0){
        gettimeofday(&args->start, NULL);
        startTimer(&args->timer1);
        startTimer(&args->timer2); 
        args->timer3 = 0; /* no partitionig phase */
    }
#endif

    /* insert tuples from the assigned part of relR to the ht */
    build_hashtable_mt(args->ht, &args->relR, &overflowbuf);
     /* wait at a barrier until each thread completes build phase */
    BARRIER_ARRIVE(args->barrier, rv);

#ifdef PERF_COUNTERS
    if(args->tid == 0){
      PCM_stop();
      PCM_log("========== Build phase profiling results ==========\n");
      PCM_printResults();
      PCM_start();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif


#ifndef NO_TIMING
    /* build phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
        stopTimer(&args->timer2); 
    }
#endif

    /* probe for matching tuples from the assigned part of relS */
    args->num_results = probe_bmhashtable(args->ht, args->bitmap,&args->relS);

#ifndef NO_TIMING
    /* for a reliable timing we have to wait until all finishes */
    BARRIER_ARRIVE(args->barrier, rv);

    /* probe phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
      stopTimer(&args->timer1); 
      gettimeofday(&args->end, NULL);
    }
#endif

#ifdef PERF_COUNTERS
    if(args->tid == 0) {
        PCM_stop();
        PCM_log("========== Probe phase profiling results ==========\n");
        PCM_printResults();
        PCM_log("===================================================\n");
        PCM_cleanup();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif

    /* clean-up the overflow buffers */
    free_bucket_buffer(overflowbuf);

    return 0;
}
void * 
AIR_thread(void * param)
{
    int rv;
    arg_vec * args = (arg_vec*) param;

    /* allocate overflow buffer for each thread */
    //bucket_buffer_t * overflowbuf;  //--no bucket by zys
    //init_bucket_buffer(&overflowbuf);  //--no bucket by zys

#ifdef PERF_COUNTERS
    if(args->tid == 0){
        PCM_initPerformanceMonitor(NULL, NULL);
        PCM_start();
    }
#endif
    
    /* wait at a barrier until each thread starts and start timer */
    BARRIER_ARRIVE(args->barrier, rv);

#ifndef NO_TIMING
    /* the first thread checkpoints the start time */
    if(args->tid == 0){
        gettimeofday(&args->start, NULL);
        startTimer(&args->timer1);
        startTimer(&args->timer2); 
        args->timer3 = 0; /* no partitionig phase */
    }
#endif

    /* insert tuples from the assigned part of relR to the ht */
    //build_hashtable_mt(args->ht, &args->relR, &overflowbuf);
    build_vector_mt(args->vec, &args->relR);  //--build vector by zys
    /* wait at a barrier until each thread completes build phase */
    BARRIER_ARRIVE(args->barrier, rv);

#ifdef PERF_COUNTERS
    if(args->tid == 0){
      PCM_stop();
      PCM_log("========== Build phase profiling results ==========\n");
      PCM_printResults();
      PCM_start();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif


#ifndef NO_TIMING
    /* build phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
        stopTimer(&args->timer2); 
    }
#endif

    /* probe for matching tuples from the assigned part of relS */
    args->num_results = probe_vector(args->vec, &args->relS);

#ifndef NO_TIMING
    /* for a reliable timing we have to wait until all finishes */
    BARRIER_ARRIVE(args->barrier, rv);

    /* probe phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
      stopTimer(&args->timer1); 
      gettimeofday(&args->end, NULL);
    }
#endif

#ifdef PERF_COUNTERS
    if(args->tid == 0) {
        PCM_stop();
        PCM_log("========== Probe phase profiling results ==========\n");
        PCM_printResults();
        PCM_log("===================================================\n");
        PCM_cleanup();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif

    /* clean-up the overflow buffers */
    //free_bucket_buffer(overflowbuf);

    return 0;
}

void * 
STARJOIN_thread(void * param)
{
    int rv;
    int fk_id;
    arg_sj * args = (arg_sj*) param;
    fk_id=args->fkid; //printf("\nfk_id:%d\n",fk_id);


#ifdef PERF_COUNTERS
    if(args->tid == 0){
        PCM_initPerformanceMonitor(NULL, NULL);
        PCM_start();
    }
#endif
    
    /* wait at a barrier until each thread starts and start timer */
    BARRIER_ARRIVE(args->barrier, rv);

#ifndef NO_TIMING
    /* the first thread checkpoints the start time */
    if(args->tid == 0){
        gettimeofday(&args->start, NULL);
        startTimer(&args->timer1);
        startTimer(&args->timer2); 
        args->timer3 = 0; /* no partitionig phase */
    }
#endif

    BARRIER_ARRIVE(args->barrier, rv);

#ifdef PERF_COUNTERS
    if(args->tid == 0){
      PCM_stop();
      PCM_log("========== Build phase profiling results ==========\n");
      PCM_printResults();
      PCM_start();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif


#ifndef NO_TIMING
    /* build phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
        stopTimer(&args->timer2); 
    }
#endif

    uint32_t i;
    int64_t matches;

    
    matches = 0;int idx;
    if(fk_id==0){  //--first fk column filtering and filling measure index by zys
      for (i = 0; i < args->fks->num_tuples; i++)
          {
           idx = args->pks->column[args->fks->column[args->FKStartIndex+i]];
		   args->MInx[args->MIStartIndex+i]=idx;  //--can be replaced with multidimensional index calculating by zys
           if (idx!=0) matches++;
          }
     }
    else{ 
         for (i = 0; i < args->fks->num_tuples; i++)
          {
          if(args->MInx[args->MIStartIndex+i]!=0){
             idx = args->pks->column[args->fks->column[args->FKStartIndex+i]];
		     args->MInx[args->MIStartIndex+i]=idx;  //--can be replaced with multidimensional index calculating by zys
             if (idx!=0) matches++;
           }  // end-if for MInx filtering
          }  //end-for for fk looping
	    }  // end-else for non-1 filter
    args->num_results = matches;

#ifndef NO_TIMING
    /* for a reliable timing we have to wait until all finishes */
    BARRIER_ARRIVE(args->barrier, rv);

    /* probe phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
      stopTimer(&args->timer1); 
      gettimeofday(&args->end, NULL);
    }
#endif

#ifdef PERF_COUNTERS
    if(args->tid == 0) {
        PCM_stop();
        PCM_log("========== Probe phase profiling results ==========\n");
        PCM_printResults();
        PCM_log("===================================================\n");
        PCM_cleanup();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif

    /* clean-up the overflow buffers */
    //free_bucket_buffer(overflowbuf);

    return 0;
}

void * 
AIRU_thread(void * param)
{
    int rv;
    arg_vec * args = (arg_vec*) param;

    /* allocate overflow buffer for each thread */
    //bucket_buffer_t * overflowbuf;  //--no bucket by zys
    //init_bucket_buffer(&overflowbuf);  //--no bucket by zys

#ifdef PERF_COUNTERS
    if(args->tid == 0){
        PCM_initPerformanceMonitor(NULL, NULL);
        PCM_start();
    }
#endif
    
    /* wait at a barrier until each thread starts and start timer */
    BARRIER_ARRIVE(args->barrier, rv);

#ifndef NO_TIMING
    /* the first thread checkpoints the start time */
    if(args->tid == 0){
        gettimeofday(&args->start, NULL);
        startTimer(&args->timer1);
        startTimer(&args->timer2); 
        args->timer3 = 0; /* no partitionig phase */
    }
#endif

    /* insert tuples from the assigned part of relR to the ht */
    //build_hashtable_mt(args->ht, &args->relR, &overflowbuf);
    build_uvector_mt(args->vec, &args->relR,args->update_ratio);  //--build update vector by zys
    /* wait at a barrier until each thread completes build phase */
    BARRIER_ARRIVE(args->barrier, rv);

#ifdef PERF_COUNTERS
    if(args->tid == 0){
      PCM_stop();
      PCM_log("========== Build phase profiling results ==========\n");
      PCM_printResults();
      PCM_start();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif


#ifndef NO_TIMING
    /* build phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
        stopTimer(&args->timer2); 
    }
#endif

    /* probe for matching tuples from the assigned part of relS */
    args->num_results = probe_uvector(args->vec, &args->relS,args->update_ratio);

#ifndef NO_TIMING
    /* for a reliable timing we have to wait until all finishes */
    BARRIER_ARRIVE(args->barrier, rv);

    /* probe phase finished, thread-0 checkpoints the time */
    if(args->tid == 0){
      stopTimer(&args->timer1); 
      gettimeofday(&args->end, NULL);
    }
#endif

#ifdef PERF_COUNTERS
    if(args->tid == 0) {
        PCM_stop();
        PCM_log("========== Probe phase profiling results ==========\n");
        PCM_printResults();
        PCM_log("===================================================\n");
        PCM_cleanup();
    }
    /* Just to make sure we get consistent performance numbers */
    BARRIER_ARRIVE(args->barrier, rv);
#endif

    /* clean-up the overflow buffers */
    //free_bucket_buffer(overflowbuf);

    return 0;
}

/** \copydoc NPO */
int64_t 
NPO(relation_t *relR, relation_t *relS, int nthreads)
{
    hashtable_t * ht;
    int64_t result = 0;
    int32_t numR, numS, numRthr, numSthr; /* total and per thread num */
    int i, rv;
    cpu_set_t set;
    arg_t args[nthreads];
    pthread_t tid[nthreads];
    pthread_attr_t attr;
    pthread_barrier_t barrier;

    uint32_t nbuckets = (relR->num_tuples / BUCKET_SIZE);
    allocate_hashtable(&ht, nbuckets);

    numR = relR->num_tuples;
    numS = relS->num_tuples;
    numRthr = numR / nthreads;
    numSthr = numS / nthreads;
    
    rv = pthread_barrier_init(&barrier, NULL, nthreads);
    if(rv != 0){
        printf("Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);
    for(i = 0; i < nthreads; i++){
        int cpu_idx = get_cpu_id(i);

        DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].tid = i;
        args[i].ht = ht;
        args[i].barrier = &barrier;

        /* assing part of the relR for next thread */
        args[i].relR.num_tuples = (i == (nthreads-1)) ? numR : numRthr;
        args[i].relR.tuples = relR->tuples + numRthr * i;
        numR -= numRthr;

        /* assing part of the relS for next thread */
        args[i].relS.num_tuples = (i == (nthreads-1)) ? numS : numSthr;
        args[i].relS.tuples = relS->tuples + numSthr * i;
        numS -= numSthr;

        rv = pthread_create(&tid[i], &attr, npo_thread, (void*)&args[i]);
        if (rv){
            printf("ERROR; return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
        /* sum up results */
        result += args[i].num_results;
    }


#ifndef NO_TIMING
    /* now print the timing results: */
    print_timing(args[0].timer1, args[0].timer2, args[0].timer3,
                relS->num_tuples, result,
                &args[0].start, &args[0].end);
#endif

    destroy_hashtable(ht);

    return result;
}

/** \copydoc NPO */
int64_t 
NPO_bm(relation_t *relR, relation_t *relS, int nthreads,int8_t * bm)
{
    hashtable_t * ht;
    int64_t result = 0;
    int32_t numR, numS, numRthr, numSthr; /* total and per thread num */
    int i, rv;
    cpu_set_t set;
    arg_bm args[nthreads];
    pthread_t tid[nthreads];
    pthread_attr_t attr;
    pthread_barrier_t barrier;

    uint32_t nbuckets = (relR->num_tuples / BUCKET_SIZE);
    allocate_hashtable(&ht, nbuckets);

    numR = relR->num_tuples;
    numS = relS->num_tuples;
    numRthr = numR / nthreads;
    numSthr = numS / nthreads;
    
    rv = pthread_barrier_init(&barrier, NULL, nthreads);
    if(rv != 0){
        printf("Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);
    for(i = 0; i < nthreads; i++){
        int cpu_idx = get_cpu_id(i);

        DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].tid = i;
        args[i].ht = ht;
        args[i].barrier = &barrier;

        /* assing part of the relR for next thread */
        args[i].relR.num_tuples = (i == (nthreads-1)) ? numR : numRthr;
        args[i].relR.tuples = relR->tuples + numRthr * i;
        numR -= numRthr;

        /* assing part of the relS for next thread */
        args[i].relS.num_tuples = (i == (nthreads-1)) ? numS : numSthr;
        args[i].relS.tuples = relS->tuples + numSthr * i;
        args[i].bitmap = bm + numRthr * i;  //--assign bitmap range by zys
        numS -= numSthr;

        rv = pthread_create(&tid[i], &attr, npo_bm_thread, (void*)&args[i]);
        if (rv){
            printf("ERROR; return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
        /* sum up results */
        result += args[i].num_results;
    }


#ifndef NO_TIMING
    /* now print the timing results: */
    print_timing(args[0].timer1, args[0].timer2, args[0].timer3,
                relS->num_tuples, result,
                &args[0].start, &args[0].end);
#endif

    destroy_hashtable(ht);

    return result;
}

/** \copydoc AIR by ZYS*/
int64_t 
AIR(relation_t *relR, relation_t *relS, int nthreads)
{
    intvector_t * vec;  //--define vector instead of hash table by ZYS
    int64_t result = 0;
    int32_t numR, numS, numRthr, numSthr; /* total and per thread num */
    int i, rv;
    cpu_set_t set;
    arg_vec args[nthreads];
    pthread_t tid[nthreads];
    pthread_attr_t attr;
    pthread_barrier_t barrier;

    uint32_t ncells = (relR->num_tuples); //--calculate cells of vector by ZYS
    vec=(intvector_t*)malloc(sizeof(intvector_t)*ncells);  //--allocate vector by ZYS

    numR = relR->num_tuples;
    numS = relS->num_tuples;
    numRthr = numR / nthreads;
    numSthr = numS / nthreads;
    
    rv = pthread_barrier_init(&barrier, NULL, nthreads);
    if(rv != 0){
        printf("Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);
    for(i = 0; i < nthreads; i++){
        int cpu_idx = get_cpu_id(i);

        DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].tid = i;
        args[i].vec = vec;
        args[i].barrier = &barrier;

        /* assing part of the relR for next thread */
        args[i].relR.num_tuples = (i == (nthreads-1)) ? numR : numRthr;
        args[i].relR.tuples = relR->tuples + numRthr * i;
        numR -= numRthr;

        /* assing part of the relS for next thread */
        args[i].relS.num_tuples = (i == (nthreads-1)) ? numS : numSthr;
        args[i].relS.tuples = relS->tuples + numSthr * i;
        numS -= numSthr;

        rv = pthread_create(&tid[i], &attr, AIR_thread, (void*)&args[i]);
        if (rv){
            printf("ERROR; return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }

    //for(i=0;i<relR->num_tuples;i++){
       // printf("%d ",vec[i]);}

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
        /* sum up results */
        result += args[i].num_results;
    }


#ifndef NO_TIMING
    /* now print the timing results: */
    print_timing(args[0].timer1, args[0].timer2, args[0].timer3,
                relS->num_tuples, result,
                &args[0].start, &args[0].end);
#endif

    free(vec);

    return result;
}
/** \copydoc AIR by ZYS, as updating test for foreign key column */
int64_t 
AIRU(relation_t *relR, relation_t *relS, int nthreads, double updateratio)
{
    intvector_t * vec;  //--define vector instead of hash table by ZYS
    int64_t result = 0;
    int32_t numR, numS, numRthr, numSthr; /* total and per thread num */
    int i, rv;
    cpu_set_t set;
    arg_vec args[nthreads];
    pthread_t tid[nthreads];
    pthread_attr_t attr;
    pthread_barrier_t barrier;

    uint32_t ncells = (relR->num_tuples); //--calculate cells of vector by ZYS
    vec=(intvector_t*)malloc(sizeof(intvector_t)*ncells);  //--allocate vector by ZYS

    numR = relR->num_tuples;
    numS = relS->num_tuples;
    numRthr = numR / nthreads;
    numSthr = numS / nthreads;
    
    rv = pthread_barrier_init(&barrier, NULL, nthreads);
    if(rv != 0){
        printf("Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);
    for(i = 0; i < nthreads; i++){
        int cpu_idx = get_cpu_id(i);

        DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].tid = i;
        args[i].vec = vec;
        args[i].barrier = &barrier;
        args[i].update_ratio = updateratio; 

        /* assing part of the relR for next thread */
        args[i].relR.num_tuples = (i == (nthreads-1)) ? numR : numRthr;
        args[i].relR.tuples = relR->tuples + numRthr * i;
        numR -= numRthr;

        /* assing part of the relS for next thread */
        args[i].relS.num_tuples = (i == (nthreads-1)) ? numS : numSthr;
        args[i].relS.tuples = relS->tuples + numSthr * i;
        numS -= numSthr;

        rv = pthread_create(&tid[i], &attr, AIRU_thread, (void*)&args[i]);
        if (rv){
            printf("ERROR; return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }

    //for(i=0;i<relR->num_tuples;i++){
       // printf("%d ",vec[i]);}

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
        /* sum up results */
        result += args[i].num_results;
    }


#ifndef NO_TIMING
    /* now print the timing results: */
    print_timing(args[0].timer1, args[0].timer2, args[0].timer3,
                relS->num_tuples, result,
                &args[0].start, &args[0].end);
#endif

    free(vec);

    return result;
}

int64_t 
STARJOIN(column_t *factT, vector_t *DimVec, vector_t *MIndex,int nthreads, vector_para *parame,int *filterflag)
{

    int64_t result = 0;
    int32_t numS, numSthr,FactTuples; //numR,numRthr --total and per thread num 
    int i, rv;
    cpu_set_t set;
    arg_sj args[nthreads];
    pthread_t tid[nthreads];
    pthread_attr_t attr;
    pthread_barrier_t barrier;
    numS = factT->num_tuples;
	FactTuples= factT->num_tuples;
    numSthr = numS / nthreads;


    rv = pthread_barrier_init(&barrier, NULL, nthreads);
    if(rv != 0){
        printf("Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }
    pthread_attr_init(&attr);
    for(i = 0; i < nthreads; i++){
        int cpu_idx = get_cpu_id(i);

        DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].tid = i;
        args[i].fkid = *filterflag;
        args[i].fks = factT;
        args[i].pks = DimVec;
        args[i].MInx = MIndex->column;
        args[i].FKStartIndex = 0;
        args[i].MIStartIndex = 0;
        args[i].barrier = &barrier;

        args[i].fks->num_tuples = (i == (nthreads-1)) ? numS : numSthr;
	    args[i].FKStartIndex+=numSthr*i;        
	    args[i].MIStartIndex+=numSthr*i;
        numS -= numSthr;

        rv = pthread_create(&tid[i], &attr, STARJOIN_thread, (void*)&args[i]);
        if (rv){
            printf("ERROR; return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
        //-- sum up results
        result += args[i].num_results;
    }

#ifndef NO_TIMING
    //-- now print the timing results: 
    print_timing(args[0].timer1, args[0].timer2, args[0].timer3,
                FactTuples, result,
                &args[0].start, &args[0].end);
#endif
 // }  //--end of fks loop by zys

    return result;

}
/** @}*/
