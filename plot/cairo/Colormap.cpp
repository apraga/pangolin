#include "Colormap.h"

Colormap::Colormap() { }

void Colormap::set(int type_i) {
  type = type_i;
  if (type == PART_ID) 
    colormapPartitions();
  else
    colormapJet();
}

Colormap::~Colormap() {
}

/* One color for each partition */
void Colormap::colormapPartitions() {
  nb_colors = 16;
  colors.resize(3*nb_colors);
  cout << "colormap :partitions" << endl;
  /* Red*/
  colors[0] = 1.;
  colors[1] = 0.;
  colors[2] = 0.;
  /* Green */
  colors[3] = 0.;
  colors[4] = 1.;
  colors[5] = 0.;
  /* Yellow       */
  colors[6] = 1.;
  colors[7] = 1.;
  colors[8] = 0.;
  /* Blue         */
  colors[9] = 0.;
  colors[10] = 0.;
  colors[11] = 1.;
  /* Magenta      */
  colors[12] = 1.;
  colors[13] = 0.;
  colors[14] = 1.;
  /* Cyan         */
  colors[15] = 0.;
  colors[16] = 1.;
  colors[17] = 1.;
  /* Orange       */
  colors[18] = 1.;
  colors[19] = 0.5;
  colors[20] = 0.2;
  /* Olive        */
  colors[21] = 0.3;
  colors[22] = 0.55;
  colors[23] = 0.;
  /* Dark pink    */
  colors[24] = 0.72;
  colors[25] = 0.47;
  colors[26] = 0.47;
  /* Sea blue     */
  colors[27] = 0.33;
  colors[28] = 0.33;
  colors[29] = 0.81;
  /* Pink         */
  colors[30] = 1.;
  colors[31] = 0.63;
  colors[32] = 0.63;
  /* Violet       */
  colors[33] = 0.62;
  colors[34] = 0.44;
  colors[35] = 0.65;
  /* Pale green   */
  colors[36] = 0.6;
  colors[37] = 0.8;
  colors[38] = 0.7;
  /* Brown        */
  colors[39] = 0.47;
  colors[40] = 0.2;
  colors[41] = 0.;
  /*Turquoise    */
  colors[42] = 0.;
  colors[43] = 0.68;
  colors[44] = 0.68;
  /* Purple     */
  colors[45] = 0.81;
  colors[46] = 0.;
  colors[47] = 0.4;

}

/* Jet color map with 256 colors */
void Colormap::colormapJet() {
  int step, i;
  float alpha;

  /* Set by user here */
  q_max = 1.;
  q_min = 0.;
  
  nb_colors = 124;
  colors.resize(3*nb_colors);
  cout << "colormap : jet" << endl;

  step = nb_colors/3;

  /* Dark blue to blue lagoon */
  for (i = 0; i < step; i ++) {
    alpha = ((float) i)/step;
    colors[3*i] = 0;
    colors[3*i + 1] = alpha*239/255.;
    colors[3*i + 2] = (1-alpha)*143/255. + alpha;
  }

  /* Blue lagoon to yellow */
  for (i = step; i < 2*step; i ++) {
    alpha = ((float) (i-step))/step;
    colors[3*i] = alpha;
    colors[3*i + 1] = 239/255.;
    colors[3*i + 2] = (1-alpha);
  }

  /* Yellow to dark red */
  for (i = 2*step; i < nb_colors; i ++) {
    alpha = ((float) (i-2*step))/step;
    colors[3*i] = (1-alpha) + alpha*127/255.;
    colors[3*i + 1] = (1-alpha)*239/255.;
    colors[3*i + 2] = 0;
  }
}

/* val is either a value (float) or the partition id (int) */
void Colormap::getColor(Color &c, float val){
  int col;

  if (type == PART_ID) {
    col = ((int) val) % nb_colors;
  }
  else {
    float alpha;
    alpha = (val - q_min)/(q_max-q_min);
    col = alpha*(nb_colors-1);
  }

  if (col < 0 || 3*col+2 > colors.size()-1){
    c = Color(1., 1., 1.);
    return;
  }

  c.r = colors[3*col];
  c.g = colors[3*col+1];
  c.b = colors[3*col+2];
}
