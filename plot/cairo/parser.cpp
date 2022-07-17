#include "parser.h"

/* Returns the type corresponding to a data with concentration / id + cell
 * coordinates 
 * File is already, we do not close it and rewind to first readable element */
int check_file_type(ifstream &file) {
  int type;
  string dummy, column;
  string first;
  int nb_col;

  getline(file, dummy);
  /* If first line is a comment, skip it */
  if (dummy.find("#") != string::npos) {
    getline(file, dummy);
    cout << "File starts with comment." << endl;
  }

  int first_elt = file.tellg();

  /* Check the number of columns */
  stringstream tmp(dummy);
  nb_col = 0;
  while (tmp) {
    tmp >> column;
    if (nb_col == 0) first = column;
    nb_col++;
  }
  /* Remove newline character */
  nb_col--;

  if (nb_col != 9)
    throw "Wrong format : int or float with 8 cell coordinates";

  /* Check if the first is integer or float */
  if (first.find(".") != string::npos) {
    type = CELL_RATIO;
    cout << "File type : cell ratio " << endl;
  }
  else {
    type = PART_ID;
    cout << "File type : partition id " << endl;
  }

  /* Rewind */
  file.seekg(first_elt);
  return type;

}

/* Find max and min for latitude, longitude. Assumes 8 columns */ 
void get_extrema(ifstream &file, float& lat_max, float& lat_min, 
    float& lon_max, float& lon_min)
{
  float cur[8];
  int i;

  lat_max = -9999;
  lat_min = 9999;
  lon_max = -9999;
  lon_min = 9999;
  while (file.good()) {
    /* Skip first value */
    file >> cur[0];
    for (i = 0; i < 8; ++i){
      file >> cur[i];
    }
    for (i = 0; i < 8; i = i+2){
      lat_max = max(lat_max, cur[i]);
      lat_min = min(lat_min, cur[i]);
    }
    for (i = 1; i < 8; i = i+2){
      lon_max = max(lon_max, cur[i]);
      lon_min = min(lon_min, cur[i]);
    }
  }
  //rewind(input);
}

