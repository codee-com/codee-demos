#ifndef __ONEMALLOC__
#define __ONEMALLOC__

#include "current.h"
#include "emf.h"
#include <stdint.h>

typedef struct {
  int size;
  int offset;

  t_fld *Ex_buf;
  t_fld *Ey_buf;
  t_fld *Ez_buf;
  t_fld *Bx_buf;
  t_fld *By_buf;
  t_fld *Bz_buf;

  t_fld *Ex;
  t_fld *Ey;
  t_fld *Ez;
  t_fld *Bx;
  t_fld *By;
  t_fld *Bz;
} t_emf_data;

typedef struct {
  int size;
  int offset;

  t_fld *Jx_buf;
  t_fld *Jy_buf;
  t_fld *Jz_buf;

  t_fld *Jx;
  t_fld *Jy;
  t_fld *Jz;
} t_current_data;

typedef struct {
  t_current_data current_data;
  t_emf_data emf_data;
} t_all_data;

void init_all_data(t_all_data *data, t_current *current, t_emf *emf);
void delete_all_data(t_all_data *data);

void copy_in_all_data(t_all_data *data, t_current *current, t_emf *emf);
void copy_out_all_data(t_all_data *data, t_current *current, t_emf *emf);

#endif //__ONEMALLOC__
