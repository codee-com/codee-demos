#include "one-malloc.h"
#include <assert.h>
#include <stdlib.h>

void init_all_data(t_all_data *data, t_current *current, t_emf *emf) {

  // t_current buffers
  data->current_data.size = current->size;
  data->current_data.offset = current->offset;

  data->current_data.Jx_buf = malloc(current->size * sizeof(t_fld));
  data->current_data.Jy_buf = malloc(current->size * sizeof(t_fld));
  data->current_data.Jz_buf = malloc(current->size * sizeof(t_fld));

  assert(data->current_data.Jx_buf);
  assert(data->current_data.Jy_buf);
  assert(data->current_data.Jz_buf);

  // Move pointers to [0][0]
  data->current_data.Jx = data->current_data.Jx_buf + data->current_data.offset;
  data->current_data.Jy = data->current_data.Jy_buf + data->current_data.offset;
  data->current_data.Jz = data->current_data.Jz_buf + data->current_data.offset;

  // t_emf buffers
  data->emf_data.size = emf->size;
  data->emf_data.offset = emf->offset;

  data->emf_data.Ex_buf = malloc(emf->size * sizeof(t_fld));
  data->emf_data.Ey_buf = malloc(emf->size * sizeof(t_fld));
  data->emf_data.Ez_buf = malloc(emf->size * sizeof(t_fld));
  data->emf_data.Bx_buf = malloc(emf->size * sizeof(t_fld));
  data->emf_data.By_buf = malloc(emf->size * sizeof(t_fld));
  data->emf_data.Bz_buf = malloc(emf->size * sizeof(t_fld));

  assert(data->emf_data.Ex_buf);
  assert(data->emf_data.Ey_buf);
  assert(data->emf_data.Ez_buf);
  assert(data->emf_data.Bx_buf);
  assert(data->emf_data.By_buf);
  assert(data->emf_data.Bz_buf);

  // Move pointers to [0][0]
  data->emf_data.Ex = data->emf_data.Ex_buf + data->emf_data.offset;
  data->emf_data.Ey = data->emf_data.Ey_buf + data->emf_data.offset;
  data->emf_data.Ez = data->emf_data.Ez_buf + data->emf_data.offset;

  data->emf_data.Bx = data->emf_data.Bx_buf + data->emf_data.offset;
  data->emf_data.By = data->emf_data.By_buf + data->emf_data.offset;
  data->emf_data.Bz = data->emf_data.Bz_buf + data->emf_data.offset;
}

void delete_all_data(t_all_data *data) {
  free(data->current_data.Jx_buf);
  free(data->current_data.Jy_buf);
  free(data->current_data.Jz_buf);

  data->current_data.Jx_buf = NULL;
  data->current_data.Jy_buf = NULL;
  data->current_data.Jz_buf = NULL;

  free(data->emf_data.Ex_buf);
  free(data->emf_data.Ey_buf);
  free(data->emf_data.Ez_buf);
  free(data->emf_data.Bx_buf);
  free(data->emf_data.By_buf);
  free(data->emf_data.Bz_buf);

  data->emf_data.Ex_buf = NULL;
  data->emf_data.Ey_buf = NULL;
  data->emf_data.Ez_buf = NULL;
  data->emf_data.Bx_buf = NULL;
  data->emf_data.By_buf = NULL;
  data->emf_data.Bz_buf = NULL;
}

void copy_in_all_data(t_all_data *data, t_current *current, t_emf *emf) {

  for (int i = 0; i < data->emf_data.size; ++i) {
    data->emf_data.Ex_buf[i] = emf->E_buf[i].x;
    data->emf_data.Ey_buf[i] = emf->E_buf[i].y;
    data->emf_data.Ez_buf[i] = emf->E_buf[i].z;

    data->emf_data.Bx_buf[i] = emf->B_buf[i].x;
    data->emf_data.By_buf[i] = emf->B_buf[i].y;
    data->emf_data.Bz_buf[i] = emf->B_buf[i].z;
  }

  for (int i = 0; i < data->current_data.size; ++i) {
    data->current_data.Jx_buf[i] = current->J_buf[i].x;
    data->current_data.Jy_buf[i] = current->J_buf[i].y;
    data->current_data.Jz_buf[i] = current->J_buf[i].z;
  }
}

void copy_out_all_data(t_all_data *data, t_current *current, t_emf *emf) {

  for (int i = 0; i < data->emf_data.size; ++i) {
    emf->E_buf[i].x = data->emf_data.Ex_buf[i];
    emf->E_buf[i].y = data->emf_data.Ey_buf[i];
    emf->E_buf[i].z = data->emf_data.Ez_buf[i];

    emf->B_buf[i].x = data->emf_data.Bx_buf[i];
    emf->B_buf[i].y = data->emf_data.By_buf[i];
    emf->B_buf[i].z = data->emf_data.Bz_buf[i];
  }

  for (int i = 0; i < data->current_data.size; ++i) {
    current->J_buf[i].x = data->current_data.Jx_buf[i];
    current->J_buf[i].y = data->current_data.Jy_buf[i];
    current->J_buf[i].z = data->current_data.Jz_buf[i];
  }
}
