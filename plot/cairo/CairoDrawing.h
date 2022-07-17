#ifndef __DRAWING_CAIRO
#define __DRAWING_CAIRO

#include <cairo.h>
#include <cairo-pdf.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include "Colormap.h"
#include "parser.h"

using namespace std;

class CairoDrawing {
  /* Pointer for cairo */
  cairo_surface_t *cs;
  cairo_t *cr;
  /* Offet and scale for drawing */
  float offset[2];
  float scale[2];
  /* Height and width of the cairo data */
  int height;
  int width;
  /* Global margin */
  int margin;
  /* Margin for the colormap */
  int wmargin;
  Colormap cmap;
  public:

  CairoDrawing();
  void init(string input, string output);
  void close();
  void setOffsetScale(ifstream &file);
  void setFromFileformat(string input);
  ~CairoDrawing();

  void drawColormap(int x0, int y0, int dx, int dy);
  void drawRectangle(vector<float> &pos, Color c);
  void adjust_position(vector<float> &pos);

  template <typename T> void drawCell(vector<float>& pos, T val);

  //void drawCell(vector<float>& pos, float val);
  void drawColormap();
};

#endif
