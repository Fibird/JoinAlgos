/* @version $Id: generator.c 3031 2012-12-07 14:37:54Z bcagri $ */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>              /* perror */
#include <stdlib.h>             /* posix_memalign */
#include <math.h>               /* fmod, pow */
#include <time.h>               /* time() */
#include <unistd.h>             /* getpagesize() */
#include <string.h>             /* memcpy() */

#include "cpu_mapping.h"        /* get_cpu_id() */
#include "affinity.h"           /* pthread_attr_setaffinity_np */
#include "genzipf.h"            /* gen_zipf() */
#include "generator.h"          /* create_relation_*() */
#include "prj_params.h"         /* RELATION_PADDING for Parallel Radix */

/* return a random number in range [0,N] */
#define RAND_RANGE(N) ((double)rand() / ((double)RAND_MAX + 1) * (N))
#define RAND_RANGE48(N,STATE) ((double)nrand48(STATE)/((double)RAND_MAX+1)*(N))
#define MALLOC(SZ) alloc_aligned(SZ+RELATION_PADDING) /*malloc(SZ+RELATION_PADDING)*/ 
#define FREE(X,SZ) free(X)

/* Uncomment the following to persist input relations to disk. */
/* #define PERSIST_RELATIONS 1 */

/** An experimental feature to allocate input relations numa-local */
int numalocalize;
int nthreads;

static int seeded = 0;
static unsigned int seedValue;

void *
alloc_aligned(size_t size)
{
    void * ret;
    int rv;
    rv = posix_memalign((void**)&ret, CACHE_LINE_SIZE, size);

    if (rv) { 
        perror("alloc_aligned() failed: out of memory");
        return 0; 
    }
    
    /** Not an elegant way of passing whether we will numa-localize, but this
        feature is experimental anyway. */
    if(numalocalize) {
        tuple_t * mem = (tuple_t *) ret;
        uint32_t ntuples = size / sizeof(tuple_t);
        numa_localize(mem, ntuples, nthreads);
    }

    return ret;
}

void 
seed_generator(unsigned int seed) 
{
    srand(seed);
    seedValue = seed;
    seeded = 1;
}

/** Check wheter seeded, if not seed the generator with current time */
static void
check_seed()
{
    if(!seeded) {
        seedValue = time(NULL);
        srand(seedValue);
        seeded = 1;
    }
}


/** 
 * Shuffle tuples of the relation using Knuth shuffle.
 * 
 * @param relation 
 */
void 
knuth_shufflefk(column_t * colu)
{
    int i;
    for (i = colu->num_tuples - 1; i > 0; i--) {
        int32_t  j              = RAND_RANGE(i);
        intkey_t tmp            = colu->column[i];
        colu->column[i] = colu->column[j];
        colu->column[j] = tmp;
    }
}
void 
knuth_shuffle(relation_t * relation)
{
    int i;
    for (i = relation->num_tuples - 1; i > 0; i--) {
        int32_t  j              = RAND_RANGE(i);
        intkey_t tmp            = relation->tuples[i].key;
        relation->tuples[i].key = relation->tuples[j].key;
        relation->tuples[j].key = tmp;
    }
}

void 
knuth_shuffle48(relation_t * relation, unsigned short * state)
{
    int i;
    for (i = relation->num_tuples - 1; i > 0; i--) {
        int32_t  j              = RAND_RANGE48(i, state);
        intkey_t tmp            = relation->tuples[i].key;
        relation->tuples[i].key = relation->tuples[j].key;
        relation->tuples[j].key = tmp;
    }
}

/**
 * Generate unique tuple IDs with Knuth shuffling
 * relation must have been allocated
 */
void
random_unique_gen(relation_t *rel) 
{
    uint32_t i;

    for (i = 0; i < rel->num_tuples; i++) {
        rel->tuples[i].key = (i+1);
        rel->tuples[i].payload = (i+1);  //--assing key and payload with the same key by zys
        //printf("%d \t %d \t ", rel->tuples[i].key, rel->tuples[i].payload);
    }

    /* randomly shuffle elements */
    knuth_shuffle(rel); //--simulating surrogate key of dimension table by zys
}
void
random_unique_gen_pk(relation_t *rel) 
{
    uint32_t i;

    for (i = 0; i < rel->num_tuples; i++) {
        rel->tuples[i].key = (i+1);
        rel->tuples[i].payload = (i+1);  //--assing key and payload with the same key by zys
        //printf("%d \t %d \t ", rel->tuples[i].key, rel->tuples[i].payload);
    }

    /* randomly shuffle elements */
    //knuth_shuffle(rel); //--simulating surrogate key of dimension table by zys
}
void
random_uniquefk_gen(column_t *col) 
{
    uint32_t i;

    for (i = 0; i < col->num_tuples; i++) {
        col->column[i]= (i+1);
        //printf("%d \t %d \t ", rel->tuples[i].key, rel->tuples[i].payload);
    }

    //printf("\ncall random_uniquefk_gen!\n");
    /* randomly shuffle elements */
    knuth_shufflefk(col); //--simulating surrogate key of dimension table by zys
}

struct create_arg_t {
    relation_t rel;
    uint32_t firstkey;
};

typedef struct create_arg_t create_arg_t;

/**
 * Create random unique keys starting from firstkey 
 */
void *
random_unique_gen_thread(void * args) 
{
    create_arg_t * arg      = (create_arg_t *) args;
    relation_t *   rel      = & arg->rel;
    uint32_t       firstkey = arg->firstkey;
    uint32_t i;

    /* for randomly seeding nrand48() */
    unsigned short state[3] = {0, 0, 0};
    unsigned int seed       = time(NULL) + * (unsigned int *) pthread_self();
    memcpy(state, &seed, sizeof(seed));
    
    for (i = 0; i < rel->num_tuples; i++) {
        rel->tuples[i].key = firstkey ++;
    }

    /* randomly shuffle elements */
    knuth_shuffle48(rel, state);

    return 0;
}

/** 
 * Just initialize mem. to 0 for making sure it will be allocated numa-local 
 */
void *
numa_localize_thread(void * args) 
{
    create_arg_t * arg = (create_arg_t *) args;
    relation_t *   rel = & arg->rel;
    uint32_t i;
    
    for (i = 0; i < rel->num_tuples; i++) {
        rel->tuples[i].key = 0;
    }

    return 0;
}

/**
 * Write relation to a file.
 */
void
write_relation(relation_t * rel, char * filename)
{
    FILE * fp = fopen(filename, "w");
    uint32_t i;

    fprintf(fp, "#KEY, VAL\n");

    for (i = 0; i < rel->num_tuples; i++) {
        fprintf(fp, "%d %d\n", rel->tuples[i].key, rel->tuples[i].payload);
    }

    fclose(fp);
}

/**
 * Generate tuple IDs -> random distribution
 * relation must have been allocated
 */
void 
random_gen(relation_t *rel, const int32_t maxid) 
{
    uint32_t i;

    for (i = 0; i < rel->num_tuples; i++) {
        rel->tuples[i].key = RAND_RANGE(maxid);
    }
}

int 
create_relation_pk(relation_t *relation, int32_t num_tuples) 
{
    check_seed();

    relation->num_tuples = num_tuples;
    relation->tuples = (tuple_t*)MALLOC(relation->num_tuples * sizeof(tuple_t));
    
    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }
  
    random_unique_gen_pk(relation);

#ifdef PERSIST_RELATIONS
    write_relation(relation, "R.tbl");
#endif

    return 0;
}

//--create dimension vector by zys
int
create_vectors_pk(vector_t * DimVec,vector_para *VecParams)
{
    check_seed();
    printf("\ncreating dimension vectors...\n");
    int i,j,counter=0,veccounter=0;
    vectorkey_t * tmpvec;
    for(i=0;i<4;i++){
          if(VecParams[i].selectivity!=0){
            DimVec[i].column=(vectorkey_t*)MALLOC(VecParams[i].num_tuples * sizeof(vectorkey_t));
            veccounter++;
            tmpvec=DimVec[i].column;
            if(VecParams[i].bitmapflag==1)
             for(j=0;j<VecParams[i].num_tuples;j++){
              if((double)rand()<VecParams[i].selectivity*(double)RAND_MAX){
                 counter++;
                 tmpvec[j]=1;
                 }
              else tmpvec[j]=0;
             }//--filling dimension vectors with [0,1] for bitmap by zys
            else if(VecParams[i].selectivity==1)
                    for(j=0;j<VecParams[i].num_tuples;j++){
                        tmpvec[j]=(int32_t)RAND_RANGE(VecParams[i].num_groups)+1;counter++;
                        //printf("%d ",tmpvec[j]);
                        }
                 else for(j=0;j<VecParams[i].num_tuples;j++)
                          if((double)rand()<VecParams[i].selectivity*(double)RAND_MAX){
                               counter++;
                               tmpvec[j]=(int32_t)RAND_RANGE(VecParams[i].num_groups)+1; 
                               //printf("%d ",tmpvec[j]);
                  //--make group number starts with 1 instead of 0 by zys
                               }
                          else tmpvec[j]=0;
          //printf("%d ",DimVec[i].column[0]);
//for(j=0;j<VecParams[i].num_tuples;j++) printf("%d ",DimVec[i].column[j]);
          printf("\nvector tuples:%d,selectivity:%f,bitmap flag:%d,non-0 positions:%d\n",VecParams[i].num_tuples,VecParams[i].selectivity,VecParams[i].bitmapflag,counter);
          counter=0;
      }
      DimVec[i].column=tmpvec;
    }
    printf("\ncreate vectors for PK table, %d vectors are created.",veccounter);
    //free(tmpvec);
    //relation->num_tuples = num_tuples;
    //relation->tuples = (tuple_t*)MALLOC(relation->num_tuples * sizeof(tuple_t));
    /*
    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }
    */
    //random_unique_gen(relation);

#ifdef PERSIST_RELATIONS
    //write_relation(relation, "R.tbl");
#endif

    return 0;
}

//--create fact fk columns and measure index by zys
int
create_fact_fk(column_t * FactColumns,vector_para *VecParams,int factrows)
{
    check_seed();
    printf("\ncreating fact fk columns...\n");
    int i=0,j,tablecounter=0;

    column_t tmpcol;intkey_t * tmp;
    for(i=0;i<4;i++){
//i=2;printf("\nselectivity:%f\n",VecParams[i].selectivity);
          if(VecParams[i].selectivity!=0){
          FactColumns[i].column=(column_t*)MALLOC(factrows * sizeof(column_t));
          if (!FactColumns[i].column) { 
             perror("out of memory when creating fact fk column.");
          return -1; 
          }
          FactColumns[i].num_tuples=factrows;
          tmp=FactColumns[i].column;
          tablecounter++;

//--revised from create_relation_fk function by zys
    int32_t iters, remainder;
   //--alternative generation method 
    iters = FactColumns[i].num_tuples / VecParams[i].num_tuples;//printf("\nitres:%d\n",iters);
    for(j = 0; j < iters; j++){
        tmpcol.num_tuples = VecParams[i].num_tuples;
        tmpcol.column = tmp + VecParams[i].num_tuples * j;
        random_uniquefk_gen(&tmpcol);
   }

    //-- if num_tuples is not an exact multiple of maxid 
    remainder = FactColumns[i].num_tuples % VecParams[i].num_tuples;
    if(remainder > 0) {
        tmpcol.num_tuples = remainder;
        tmpcol.column = tmp + VecParams[i].num_tuples * iters;
        random_uniquefk_gen(&tmpcol);
    }
/**/
    FactColumns[i].column=tmp;
    printf("\ncreate column[%d] for fK table, %d tuples are created.\n",i+1,factrows);

     }//--end for if(VecParams[i].selectivity>0)    by zys 
    }//--end for table loop by zys
   // free(tmpvec);

    return 0;
}

int 
parallel_create_relation_pk(relation_t *relation, int32_t num_tuples, 
                            uint32_t nthreads) 
{
    uint32_t i, rv;
    uint32_t offset = 0;

    check_seed();

    relation->num_tuples = num_tuples;

    /* we need aligned allocation of items */
    relation->tuples = (tuple_t*) MALLOC(num_tuples * sizeof(tuple_t));

    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }

    create_arg_t args[nthreads];
    pthread_t tid[nthreads];
    cpu_set_t set;
    pthread_attr_t attr;

    unsigned int pagesize;
    unsigned int npages;
    unsigned int npages_perthr;
    unsigned int ntuples_perthr;
    unsigned int ntuples_lastthr;

    pagesize        = getpagesize();
    npages          = (num_tuples * sizeof(tuple_t)) / pagesize + 1;
    npages_perthr   = npages / nthreads;
    ntuples_perthr  = npages_perthr * (pagesize/sizeof(tuple_t));
    ntuples_lastthr = num_tuples - ntuples_perthr * (nthreads-1);

    pthread_attr_init(&attr);

    for( i = 0; i < nthreads; i++ ) {
        int cpu_idx = get_cpu_id(i);
        
        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].firstkey       = offset + 1;
        args[i].rel.tuples     = relation->tuples + offset;
        args[i].rel.num_tuples = (i == nthreads-1) ? ntuples_lastthr 
                                 : ntuples_perthr;
        offset += ntuples_perthr;

        rv = pthread_create(&tid[i], &attr, random_unique_gen_thread, 
                            (void*)&args[i]);
        if (rv){
            fprintf(stderr, "[ERROR] pthread_create() return code is %d\n", rv);
            exit(-1);
        }
    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
    }

    /* randomly shuffle elements */
    knuth_shuffle(relation);

    return 0;
}

int 
numa_localize(tuple_t * relation, int32_t num_tuples, uint32_t nthreads) 
{
    uint32_t i, rv;
    uint32_t offset = 0;

    /* we need aligned allocation of items */
    create_arg_t args[nthreads];
    pthread_t tid[nthreads];
    cpu_set_t set;
    pthread_attr_t attr;

    unsigned int pagesize;
    unsigned int npages;
    unsigned int npages_perthr;
    unsigned int ntuples_perthr;
    unsigned int ntuples_lastthr;

    pagesize        = getpagesize();
    npages          = (num_tuples * sizeof(tuple_t)) / pagesize + 1;
    npages_perthr   = npages / nthreads;
    ntuples_perthr  = npages_perthr * (pagesize/sizeof(tuple_t));
    ntuples_lastthr = num_tuples - ntuples_perthr * (nthreads-1);

    pthread_attr_init(&attr);

    for( i = 0; i < nthreads; i++ ) {
        int cpu_idx = get_cpu_id(i);
        
        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].firstkey       = offset + 1;
        args[i].rel.tuples     = relation + offset;
        args[i].rel.num_tuples = (i == nthreads-1) ? ntuples_lastthr 
                                 : ntuples_perthr;
        offset += ntuples_perthr;

        rv = pthread_create(&tid[i], &attr, numa_localize_thread, 
                            (void*)&args[i]);
        if (rv){
            fprintf(stderr, "[ERROR] pthread_create() return code is %d\n", rv);
            exit(-1);
        }
    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
    }

    return 0;
}


int 
create_relation_fk(relation_t *relation, int32_t num_tuples, const int32_t maxid)
{
    int32_t i, iters, remainder;
    relation_t tmp;

    check_seed();

    relation->num_tuples = num_tuples;
    relation->tuples = (tuple_t*)MALLOC(relation->num_tuples * sizeof(tuple_t));
      
    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }
  
    /* alternative generation method */
    iters = num_tuples / maxid;
    for(i = 0; i < iters; i++){
        tmp.num_tuples = maxid;
        tmp.tuples = relation->tuples + maxid * i;
        random_unique_gen(&tmp);
    }

    /* if num_tuples is not an exact multiple of maxid */
    remainder = num_tuples % maxid;
    if(remainder > 0) {
        tmp.num_tuples = remainder;
        tmp.tuples = relation->tuples + maxid * iters;
        random_unique_gen(&tmp);
    }

#ifdef PERSIST_RELATIONS
    write_relation(relation, "S.tbl");
#endif

    return 0;
}

int   //--creating fk table with bitmap filtering by zys
create_relation_bmfk(relation_t *relation, int8_t *bitmap, int32_t num_tuples, const int32_t maxid, double sele)
{
    int32_t i, iters, remainder;
    relation_t tmp;

    check_seed();

    relation->num_tuples = num_tuples;
    relation->tuples = (tuple_t*)MALLOC(relation->num_tuples * sizeof(tuple_t));
    //bitmap=(int8_t*)MALLOC(relation->num_tuples * sizeof(int8_t));  
    //--create bitmap for relation S by zys

    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }
     int32_t counter=0;
     for(i=0;i<relation->num_tuples;i++){
          if((double)rand()<sele*(double)RAND_MAX) {
              bitmap[i]=1;counter++;}
          else bitmap[i]=0;
     }
    printf("\nselectivity:%f, 1 position of bitmap is: %d[%d]\n",sele,counter,relation->num_tuples);
    //for(i=0;i<relation->num_tuples;i++) printf("%d ",bitmap[i]);


    /* alternative generation method */
    iters = num_tuples / maxid;
    for(i = 0; i < iters; i++){
        tmp.num_tuples = maxid;
        tmp.tuples = relation->tuples + maxid * i;
        random_unique_gen(&tmp);
    }

    /* if num_tuples is not an exact multiple of maxid */
    remainder = num_tuples % maxid;
    if(remainder > 0) {
        tmp.num_tuples = remainder;
        tmp.tuples = relation->tuples + maxid * iters;
        random_unique_gen(&tmp);
    }

#ifdef PERSIST_RELATIONS
    write_relation(relation, "S.tbl");
#endif

    return 0;
}

/** 
 * Create a foreign-key relation using the given primary-key relation and
 * foreign-key relation size. Keys in pkrel is randomly distributed in the full
 * integer range.
 * 
 * @param fkrel [output] foreign-key relation
 * @param pkrel [input] primary-key relation
 * @param num_tuples 
 * 
 * @return 
 */
int 
create_relation_fk_from_pk(relation_t *fkrel, relation_t *pkrel,
                           int32_t num_tuples) 
{
    int rv, i, iters, remainder;

    rv = posix_memalign((void**)&fkrel->tuples, CACHE_LINE_SIZE, 
                        num_tuples * sizeof(tuple_t) + RELATION_PADDING);

    if (rv) { 
        perror("aligned alloc failed: out of memory");
        return 0; 
    }

    fkrel->num_tuples = num_tuples;

    /* alternative generation method */
    iters = num_tuples / pkrel->num_tuples;
    for(i = 0; i < iters; i++){
        memcpy(fkrel->tuples + i * pkrel->num_tuples, pkrel->tuples, 
               pkrel->num_tuples * sizeof(tuple_t));
    }

    /* if num_tuples is not an exact multiple of pkrel->num_tuples */
    remainder = num_tuples % pkrel->num_tuples;
    if(remainder > 0) {
        memcpy(fkrel->tuples + i * pkrel->num_tuples, pkrel->tuples, 
               remainder * sizeof(tuple_t));
    }

    knuth_shuffle(fkrel);

    return 0;
}

int create_relation_nonunique(relation_t *relation, int32_t num_tuples,
                              const int32_t maxid) 
{
    check_seed();

    relation->num_tuples = num_tuples;
    relation->tuples = (tuple_t*)MALLOC(relation->num_tuples * sizeof(tuple_t));
    
    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }

    random_gen(relation, maxid);

    return 0;
}

double 
zipf_ggl(double * seed) 
{
    double t, d2=0.2147483647e10;
    t = *seed;
    t = fmod(0.16807e5*t, d2);
    *seed = t;
    return (t-1.0e0)/(d2-1.0e0);
}

int 
create_relation_zipf(relation_t * relation, int32_t num_tuples,
                     const int32_t maxid, const double zipf_param) 
{
    check_seed();

    relation->num_tuples = num_tuples;
    relation->tuples = (tuple_t*) MALLOC(relation->num_tuples * sizeof(tuple_t));

    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }

    gen_zipf(num_tuples, maxid, zipf_param, &relation->tuples);

    return 0;
}

void 
delete_relation(relation_t * rel) 
{
    /* clean up */
    FREE(rel->tuples, rel->num_tuples * sizeof(tuple_t));
}
