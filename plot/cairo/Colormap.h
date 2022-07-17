#ifndef __COLORMAP_CAIRO
#define __COLORMAP_CAIRO

#include <iostream>
#include <stdlib.h>
#include "parser.h"

using namespace std;

class Color {
  public:
    float r;
    float g;
    float b;
    Color(): r(0), g(0), b(0){}; 
    Color(double r1, double g1, double b1): r(r1), g(g1), b(b1){}; 
};

class Colormap {
  public:
  float q_max;
  float q_min;
  int nb_colors;
  int type;
  vector<float> colors;

  void colormapPartitions();
  void colormapJet();
  void getColor(Color &c, float val);
  //void getColor(Color &c, int id);

  Colormap();
  void set(int type);
  ~Colormap();
};

#endif


