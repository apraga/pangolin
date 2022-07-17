#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Compare 2 hdf5 files for Pangolin to ensure parallel I/O is the same as
 * sequential I/O 
 * Format must be lat, lon, data (in different datasets) */

typedef struct list_strings {
  char* list[3];
  int pos;
} list_strings;

/* Operation function, to store the names of all 3 datasets */
herr_t get_names (hid_t loc_id, const char *name, const H5L_info_t *info,
    void *op_data)
{
  H5G_stat_t statbuf;
  list_strings *names = (list_strings *) op_data;
  int cur = names->pos;
  int n ;

  /* Only read 3 datasets */
  H5Gget_objinfo(loc_id, name, 0, &statbuf);
  if ( statbuf.type == H5G_DATASET) {
    n = strlen(name);
    if (cur > 2) {
      printf("Error, more that 3 datasets in file\n");
      exit(1);
    }
    /* Copies the name into the list */
    names->list[cur] = (char*) malloc(n*sizeof(char));
    strcpy(names->list[cur], name);
    names->pos++;
  }
  else {
    printf(" Object with name %s is not a dataset, skipping \n", name);
  }
  return 0;
}

/* Check the datasests have the same names */
void check_datasets_names(list_strings names[2]) {
  int i, j, nb_errors;

  nb_errors = 0;
  for (j = 0; j < 3; ++j) {
    if (strcmp(names[0].list[j], names[1].list[j]) != 0) {
      nb_errors++;
    }
  }
  if (nb_errors > 0) {
    printf("Datasets do not match.\n");
    for (i = 0; i < 2; ++i) {
      printf("File %d has ", i);
      for (j = 0; j < 3; ++j) {
        printf("%s ", names[i].list[j]);
      }
      printf("\n");
    }
    exit(1);
  }
  else {
    printf("Datasets names match: ");
    for (j = 0; j < 3; ++j) {
      printf("%s ", names[i].list[j]);
    }
    printf("\n");
  }
}

void check_dimensions(int rank[2], hsize_t dims_out[2][1]) {
  if (rank[0] != rank[1]) {
    printf("Different ranks %d %d\n", rank[0], rank[1]);
    exit(EXIT_FAILURE);
  }
  if (dims_out[0][0] != dims_out[0][1] ) {
    printf("Different dims %d %d\n", dims_out[0][0], dims_out[0][1]);
    exit(EXIT_FAILURE);
  }
  printf("dims %d %d\n", dims_out[0][0], dims_out[0][1]);
}

void print_results(int nb_errors[3], list_strings names[2]) {
  int j;

  if (nb_errors[0] || nb_errors[1] || nb_errors[2]) {
    for (j = 0; j < 3; j++) {
      if (nb_errors[j] > 0) 
        printf("Difference in %s : %d \n", names[0].list[j], nb_errors[j] );
      else
        printf("%s : ok \n", names[0].list[j]);
    }

    exit(EXIT_FAILURE);
  }
  else {
    printf("OK\n");
    exit(EXIT_SUCCESS);
  }
}

/*  Open each file with its datasets  */
void open_datasets(hid_t file[2], char* argv[], list_strings names[2], 
    hid_t dataspace[2], hid_t dataset[2][3], int rank[2], hsize_t dims_out[2][1]) {
  int i, j, status_n;

  for (i = 0; i < 2; ++i) {
    file[i] = H5Fopen(argv[i+1], H5F_ACC_RDONLY, H5P_DEFAULT);
    for (j = 0; j < 3; ++j) {
      dataset[i][j] = H5Dopen(file[i], names[i].list[j], H5P_DEFAULT);
    }

    dataspace[i] = H5Dget_space(dataset[i][0]);    /* dataspace handle */
    rank[i] = H5Sget_simple_extent_ndims(dataspace[i]);
    status_n  = H5Sget_simple_extent_dims(dataspace[i], dims_out[i], NULL);
  }
}

/* Filenames are given as arguments*/
int main(int argc, char* argv[]) {
  hid_t file[2], dataset[2][3];         /* handles */
  hid_t dataspace[2];   
  hsize_t dims_out[2][1];           /* dataset dimensions */      
  herr_t status;                             
  double *data, *lat, *lon; 
  double precision = 2e-12;
  int i, j, rank[2], status_n;
  int nb_cells, nb_errors[3];
  struct list_strings names[2];

  if (argc != 3) {
    printf("Need 2 arguments\n");
    exit(EXIT_FAILURE);
  }

  /* First we search the names of the 3 datasets for each file*/
  for (i = 0; i < 2; ++i) {
    names[i].pos = 0;
    printf("Reading %s...\n", argv[i+1]);
    file[i] = H5Fopen(argv[i+1], H5F_ACC_RDONLY, H5P_DEFAULT);

    H5Literate(file[i], H5_INDEX_NAME, H5_ITER_INC, NULL, get_names, (void *)&names[i]);
  }

  check_datasets_names(names);
  open_datasets(file, argv, names, dataspace, dataset, rank, dims_out);
  check_dimensions(rank, dims_out);
  
  nb_cells = dims_out[0][0];

  /* Init data */
  data = (double*) malloc(2*nb_cells*sizeof(double*));
  lat = (double*) malloc(2*nb_cells*sizeof(double*));
  lon = (double*) malloc(2*nb_cells*sizeof(double*));

  for (j = 0; j < 2*nb_cells; j++) {
    data[j] = 0;
    lat[j] = 0;
    lon[j] = 0;
  }

  /* Read data */
  for (i = 0; i < 2; ++i) {
    status = H5Dread(dataset[i][0], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &lat[i*nb_cells]);
    status = H5Dread(dataset[i][1], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &lon[i*nb_cells]);
    status = H5Dread(dataset[i][2], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &data[i*nb_cells]);
  }

  for (i = 0; i < 3; ++i) 
    nb_errors[i] = 0;

  /* Compare the arrays */
  for (j = 0; j < nb_cells; j++) {
    if (fabs(lat[j] - lat[j+nb_cells]) > precision) {
      nb_errors[0]++;
      printf("Error for %s line %d: %.10f %.10f\n", names[0].list[0], j, lat[j], lat[j+nb_cells]);
    }
  }

  for (j = 0; j < nb_cells; j++) {
    if (fabs(lon[j] - lon[j+nb_cells]) > precision) {
      nb_errors[1]++;
      printf("Error for %s line %d: %.10f %.10f\n", names[0].list[1], j, lon[j], lon[j+nb_cells]);
    }
  }

  for (j = 0; j < nb_cells; j++) {
    if (fabs(data[j] - data[j+nb_cells]) > precision) {
      nb_errors[2]++;
      printf("Error for %s line %d: %.10f %.10f", names[0].list[2], j, data[j], data[j+nb_cells]);
      printf(" of %.12f\n", fabs(data[j] - data[j+nb_cells]));
    }
  }

  print_results(nb_errors, names);
 
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 3; ++j) {
      H5Dclose(dataset[i][j]);
    }
    H5Sclose(dataspace[i]);
    H5Fclose(file[i]);
  }

  free(data);
  free(lat);
  free(lon);

  return 0;
}     


