#pragma once

#include "display.h"
#include "view.h"
#include "render.h"
#include "pthread_pool.h"

#include <stdio.h>
#include <assert.h>
#include <sys/syscall.h>
#include <errno.h>
#include <stdatomic.h>
#include <string.h>

#define threads_mutex_t          pthread_mutex_t
#define threads_mutex_lock(m)    pthread_mutex_lock(m)
#define threads_mutex_unlock(m)  pthread_mutex_unlock(m)
#define threads_mutex_destroy(m) pthread_mutex_destroy(m)
#define threads_mutex_init(m, p) pthread_mutex_init(m, p)


typedef struct threads_t
{
  pthread_pool_t pool;
  pthread_pool_worker_t *worker;
  pthread_pool_task_t *task;

  uint32_t *cpuid;

  // the job queue is really just a counter, the worker function view_* will know how to interpret that.
  uint64_t counter, end;
  uint64_t px_counter, px_end;

  // sync init:
  atomic_int_fast32_t init;
}
threads_t;

// per-worker initialisation
static void *threads_tls_init_one(void *arg)
{
  // remember our thread id:
  rt_tls.tid = *(uint64_t *)arg;

  // pin ourselves to a cpu:
  common_setaffinity(rt.threads->cpuid[rt_tls.tid]);

  // now init submodules:
  rt_tls.render = render_tls_init();
  // TODO: more module initialisation here, too.

  // pid_t tid = syscall(SYS_gettid);
  // fprintf(stdout, "[worker %03u] pinned to cpu %03u by thread %03d\n", rt_tls.tid, rt.threads->cpuid[rt_tls.tid], tid);
  // make sure we only init one thread here, and don't pick another task
  rt.threads->init++;
  while(rt.threads->init < rt.num_threads) sched_yield();
  return 0;
}

static void *threads_tls_cleanup_one(void *arg)
{
  render_tls_cleanup(rt_tls.render);
  rt.threads->init++;
  while(rt.threads->init < rt.num_threads) sched_yield();
  return 0;
}

static inline threads_t *threads_init()
{
  threads_t *t = (threads_t *)malloc(sizeof(threads_t));
  pthread_pool_create(&t->pool, NULL);
  t->worker = (pthread_pool_worker_t *)malloc(sizeof(pthread_pool_worker_t)*rt.num_threads);
  for(int k=0;k<rt.num_threads;k++)
    pthread_pool_worker_init(t->worker + k, &t->pool, NULL);
  t->counter = t->end = 0;
  t->init = 0;

  t->task = (pthread_pool_task_t *)malloc(sizeof(pthread_pool_task_t)*rt.num_threads);

  const char *def_file = "affinity";
  const char *filename = def_file;
  for(int k=0;k<rt.argc;k++) if(!strcmp(rt.argv[k], "--affinity") && k < rt.argc-1)
      filename = rt.argv[++k];
  // load cpu affinity mask, if any.
  // create one of those by doing something like:
  //  for i in $(seq 0 1 11); do cat /sys/devices/system/cpu/cpu$i/topology/thread_siblings_list; done
  //  | sort -g | uniq
  // or numactl --hardware
  // and then edit it to your needs (this list will start with one thread per core, no hyperthreading used)
  t->cpuid = (uint32_t *)malloc(sizeof(uint32_t)*rt.num_threads);
  for(int k=0;k<rt.num_threads;k++)
    t->cpuid[k] = k; // default init
  FILE *f = fopen(filename, "rb");
  if(f)
  {
    int k = 0;
    while(!feof(f))
    {
      if(k >= rt.num_threads) break;
      uint32_t cpu1;
      int read = fscanf(f, "%d", &cpu1);
      if(read != 1)
      {
        fprintf(stderr, "[threads] there is a problem with your affinity file `%s' in line %d!\n", filename, k);
        break;
      }
      read = fscanf(f, "%*[^\n]");
      fgetc(f); // read newline
      t->cpuid[k++] = cpu1;
    }
    if(k < rt.num_threads)
      fprintf(stderr, "[threads] not enough entries in your affinity file `%s'! (%d/%d threads)\n", filename, k, rt.num_threads);
    fclose(f);
  }

  return t;
}

static inline void threads_tls_init(threads_t *t)
{
  // init thread local storage:
  t->init = 0;
  uint64_t tid[rt.num_threads];
  for(uint64_t k=0;k<rt.num_threads;k++)
  {
    tid[k] = k;
    pthread_pool_task_init(t->task + k, &t->pool, threads_tls_init_one, tid+k);
  }
  pthread_pool_wait(&t->pool);
}

static inline void threads_tls_cleanup(threads_t *t)
{
  t->init = 0;
  // cleanup thread local storage:
  uint64_t tid[rt.num_threads];
  for(uint64_t k=0;k<rt.num_threads;k++)
  {
    tid[k] = k;
    pthread_pool_task_init(t->task + k, &t->pool, threads_tls_cleanup_one, tid+k);
  }
  pthread_pool_wait(&t->pool);
}


static inline void threads_cleanup(threads_t *t)
{
  threads_tls_cleanup(t);
  pthread_pool_destroy(&t->pool);
  free(t->cpuid);
  free(t->task);
  free(t->worker);
  free(t);
}

