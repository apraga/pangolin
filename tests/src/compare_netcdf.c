#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

/* Compare 2 netcdf files for Pangolin to ensure parallel I/O is the same as
 * sequential I/O */
typedef enum {
  id_nb_cells=0, id_lat=1, id_lon=2, id_q=3
} id_type ;

void handle_error(int e) {
  if (!e) return;
  printf("Error: %s\n", nc_strerror(e)); 
  exit(EXIT_FAILURE);
}

#define check_alloc(q) { if (q == NULL) { printf("Malloc failed\n"); exit(EXIT_FAILURE); }}

int main(int argc, char* argv[]) {

  if (argc != 3) {
    printf("Need 2 arguments\n");
    exit(EXIT_FAILURE);
  }
  char* fname1 = argv[1];
  char* fname2 = argv[2];
  int ncid[2];
  size_t nb_cells[2];
  /* Array of ids for every thing */
  int ids[4][2];
  int status;
  int i, k, offset;
  double *q, *lat, *lon;
  double error;
  int nb_errors[3] = {0,0,0};
  id_type id;

  status = nc_open(fname1, NC_NOWRITE, &ncid[0]);
  handle_error(status);
  status = nc_open(fname2, NC_NOWRITE, &ncid[1]);
  printf("File 1: %s\n", fname1);
  printf("File 2: %s\n", fname2);
  handle_error(status);

  /* Get number of cells */
  id = id_nb_cells;
  for (i = 0; i < 2; ++i) {
    status = nc_inq_dimid(ncid[i], "nb_cells", &ids[id][i]);
    handle_error(status);

    status = nc_inq_dimlen(ncid[i], ids[id][i], &nb_cells[i]);
    handle_error(status);
  }

  printf("nb cells %d, %d\n", nb_cells[0], nb_cells[1]);
  if (nb_cells[0] != nb_cells[1]) {
    printf("Different nb cells %d, %d\n", nb_cells[0], nb_cells[1]);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < 2; ++i) {
    status = nc_inq_varid(ncid[i], "center_lat", &ids[id_lat][i]) ;
    handle_error(status);
    status = nc_inq_varid(ncid[i], "center_lon", &ids[id_lon][i]) ;
    handle_error(status);
    status = nc_inq_varid(ncid[i], "q", &ids[id_q][i]) ;
    handle_error(status);
  }

  /* Each variable contains both array */
  q = (double*) malloc(2*sizeof(double)*nb_cells[0]);
  check_alloc(q);
  lat = (double*) malloc(2*sizeof(double)*nb_cells[0]);
  check_alloc(lat);
  lon = (double*) malloc(2*sizeof(double)*nb_cells[0]);
  check_alloc(lon);

  /* Read each variable in the proper half of the array */
  for (i = 0; i < 2; ++i) {
    offset = i*nb_cells[0];
    status = nc_get_var_double(ncid[i], ids[id_lat][i], lat + offset ) ;
    handle_error(status);
    status = nc_get_var_double(ncid[i], ids[id_lon][i], lon + offset) ;
    handle_error(status);
    status = nc_get_var_double(ncid[i], ids[id_q][i], q + offset) ;
    handle_error(status);
  }

  /* Compare values */
  offset = nb_cells[0];
  error = 1e-13;
  for (i = 0; i < offset; ++i) {
    if (fabs(lat[i] - lat[offset+i]) > error) {
      printf("lat: %.15f %.15f \n", lat[i], lat[offset+i]);
      nb_errors[0]++;
    }
    if (fabs(lon[i] - lon[offset+i]) > error) {
      printf("lon: %.15f %.15f \n", lon[i], lon[offset+i]);
      nb_errors[1]++;
    }
    if (fabs(q[i] - q[offset+i]) > error){
      printf("lat: %.15f %.15f\n", q[i], q[offset+i]);
      nb_errors[2]++;
    }
  }
  if (nb_errors[0] || nb_errors[1] || nb_errors[2]) {
    if (nb_errors[0] > 0) 
      printf("Difference in nb_lat : %d \n", nb_errors[0] );

    if (nb_errors[1] > 0) 
      printf("Difference in nb_lon : %d \n", nb_errors[1] );

    if (nb_errors[2] > 0) 
      printf("Difference in q : %d \n", nb_errors[2] );

    exit(EXIT_FAILURE);
  }
  else {
    printf("OK\n");
    exit(EXIT_SUCCESS);
  }
  

  free(q);
  free(lat);
  free(lon);

}
