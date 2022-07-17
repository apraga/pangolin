#ifndef __PARSER_CAIRO_
#define __PARSER_CAIRO_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#define CELL_RATIO 1
#define PART_ID 2

using namespace std;

int check_file_type(ifstream &file);
void get_extrema(ifstream &file, float& lat_max, float& lat_min, 
    float& lon_max, float& lon_min);

#endif
