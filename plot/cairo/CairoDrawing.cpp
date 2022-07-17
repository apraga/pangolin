#include "CairoDrawing.h"

CairoDrawing::CairoDrawing() { }

/* Create output file : width, height 
 * Check filetype
 * Create colormap */
void CairoDrawing::init(string input, string output) {

  height = 600;
  width = 800;
  margin = 7;
  wmargin = 100;

  /* Cairo initialization */
  cs = cairo_pdf_surface_create(output.c_str(), width + wmargin, height);
  cr = cairo_create (cs);
  cairo_set_line_width(cr,0.1);
  /* Some margin */
  cairo_translate(cr, margin, margin);
  height -= 2*margin;
  width -= 2*margin;

  setFromFileformat(input);
}

void CairoDrawing::setFromFileformat(string input)  {
  ifstream file;
  int type;

  /* Find file format and define offset, scale and colormap */
  file.open(input.c_str());
  if (!file.is_open()) {throw string("Unable to open input file ");}

  type = check_file_type(file);
  cout << "type " << type << endl;
  if (type != CELL_RATIO && type != PART_ID) 
    throw  "Wrong file type : should be concentration/id + 8 cell coordinates";

  cmap.set(type);
  setOffsetScale(file);
  file.close();
}

/* Find offset and scale on the file (already open and rewound) */
void CairoDrawing::setOffsetScale(ifstream &file) {
  float lat_max, lat_min;
  float lon_max, lon_min;
  float tmp;

  get_extrema(file, lat_max, lat_min, lon_max, lon_min);
  cout << "lat max min " << lat_max << " " << lat_min << endl;
  cout << "lon max min " << lon_max << " " << lon_min << endl;

  /* (x - offset)*scale must be between 0 and height */
  offset[0] = lat_min;
  scale[0] = ((float) height) / (lat_max - lat_min);

  offset[1] = lon_min;
  scale[1] = ((float) width) / (lon_max - lon_min);
  tmp = ((float) width) / (lon_max - lat_min);
  //cout << "width lon " << width << " " << lon_max << " " << lon_min << endl;
  //cout << "offset " << offset[0] << " " << offset[1] << endl;
  //cout << "scale " << scale[0] << " " << scale[1] << endl;
}

/* Draw to output file and close it */
void CairoDrawing::close() {
  /* Needed for PDF output */
  cairo_show_page(cr);
  cairo_destroy(cr);

  cairo_surface_flush(cs);
  cairo_surface_destroy(cs);
}

CairoDrawing::~CairoDrawing() {}

/* Draw an arrow */
void draw_arrow(float x0, float y0, float x1, float y1, cairo_t* cr)
{
  /* Head length and angle */
  float length = 1.;
  float theta = M_PI/4.;
  float beta, k;
  float x2, y2;
  float diff_x, diff_y;
  /* Draw the line */
  y0 = 90. - y0;
  y1 = 90. - y1;
  diff_x = x1 - x0;
  diff_y = y1 - y0;
  x1 = x0 + 0.9*diff_x;
  y1 = y0 + 0.9*diff_y;
  x0 = x0 + 0.1*diff_x;
  y0 = y0 + 0.1*diff_y;

  cairo_move_to(cr, x1, y1);
  cairo_line_to(cr, x0, y0);
  cairo_set_line_width(cr, 0.3);
  cairo_stroke(cr);

  /* Draw the head */
  beta = M_PI + atan2(y1 - y0, x1 - x0) + theta;
  k = length/cos(theta);
  x2 = x1 + k*cos(beta);
  y2 = y1 + k*sin(beta);
  /* Cleaner arrows : close the head*/
  diff_x = 0.1*(x2 - x1);
  diff_y = 0.1*(y2 - y1);
  cairo_move_to(cr, x1 - diff_x, y1 - diff_y);
  cairo_line_to(cr, x2, y2);
  cairo_stroke(cr);

  beta -= 2*theta;
  x2 = x1 + k*cos(beta);
  y2 = y1 + k*sin(beta);
  cairo_move_to(cr, x1, y1);
  cairo_line_to(cr, x2, y2);
  cairo_stroke(cr);

}


/* Create a rectangle from its position with a color */
/* Data is stored like this :
 * (y,x) ---- (y,x)
 *   |          |
 * (y,x) ---- (y,x) */
void CairoDrawing::drawRectangle(vector<float> &pos, Color c)
{
  /* Lat-lon format, so input is (y,x) */
  float rwidth = fabs(pos[3] - pos[1]);
  float rheight = fabs(pos[4] - pos[2]);
  float x0,y0;

  /* Start from upper left corner */
  x0 = pos[1];
  /* y is already reversed*/
  y0 = pos[0];
  cairo_rectangle (cr, x0, y0,  rwidth, rheight);

  //cout << "rectangle at " << x0 << " " << y0 << endl;
  //cout << "width height " << width << " " << height << endl;
  cairo_set_source_rgb(cr, c.r, c.g, c.b);
  //cout << "rgb " << c.r << " " << c.g << " " << c.b << endl;

  /* No cell borders */
  //cairo_fill(cr);

  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_stroke(cr);
}

void CairoDrawing::drawColormap(int x0, int y0, int dx, int dy)
{
  int i, k;
  int step;
  char test[5];

  step = ((float) dy) / cmap.nb_colors;
  k = 0;
  for (i = cmap.nb_colors-1; i > -1; i--) {
    cairo_set_source_rgb(cr, cmap.colors[3*i], cmap.colors[3*i+1], cmap.colors[3*i+2]);
    cairo_rectangle(cr, x0, y0 + k*step, dx, step);
    cairo_fill(cr);
    k++;
  }
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_move_to(cr, x0+dx, y0+10);
  sprintf(test, "%2.3f",cmap.q_max);
  cairo_show_text (cr, test);

  cairo_move_to(cr, x0+dx, y0+step*cmap.nb_colors);
  sprintf(test, "%2.3f",cmap.q_min);
  cairo_show_text (cr, test);
  //printf("%d %d %d %d \n", x0, y0, dx, dy);
}

/* Adjust position to fit on the plot. Also revert latitude coordinates as the
 * plot is downwards. */ 
void CairoDrawing::adjust_position(vector<float> &pos)
{
  //cout << "offset scale" << offset[0] << " " << scale[0] << endl;
  //cout << "offset scale" << offset[1] << " " << scale[1] << endl;
  int i;
  /* Data is on the format (latitude, longitude) */
  for (i = 0; i < 8; i = i+2)
  {
    pos[i] = (pos[i] - offset[0])*scale[0];
    /* Revert y */
    pos[i] = ((float) height) - pos[i];
    pos[i + 1] = (pos[i + 1] - offset[1])*scale[1];
  }
}


/* Draw a cell from its corner and value at the center */
template <typename T> 
void CairoDrawing::drawCell(vector<float> &pos, T val) {
  Color c;

  cmap.getColor(c, val);
  adjust_position(pos);
  drawRectangle(pos, c);
}

template void CairoDrawing::drawCell<int>(vector<float>& pos, int val);
template void CairoDrawing::drawCell<float>(vector<float>& pos, float val);

/* Draw the colormap with label on a corner (of width = wmargin)*/
void CairoDrawing::drawColormap()
{
  int i, k;
  int step;
  char test[5];

  int dx = wmargin/3;
  int x0 = width + margin;
  int y0 = 0;
  int dy = height/2;;

  step = ((float) dy) / cmap.nb_colors;
  k = 0;
  for (i = cmap.nb_colors-1; i > -1; i--) {
    cairo_set_source_rgb(cr, cmap.colors[3*i], cmap.colors[3*i+1], cmap.colors[3*i+2]);
    cairo_rectangle(cr, x0, y0 + k*step, dx, step);
    cairo_fill(cr);
    k++;
  }
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_move_to(cr, x0+dx, y0+10);
  sprintf(test, "%2.3f",cmap.q_max);
  cairo_show_text (cr, test);

  cairo_move_to(cr, x0+dx, y0+step*cmap.nb_colors);
  sprintf(test, "%2.3f",cmap.q_min);
  cairo_show_text (cr, test);
  //printf("%d %d %d %d \n", x0, y0, dx, dy);
}
