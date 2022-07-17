#include "CairoDrawing.h"

/* Plot a grid with the format : id with the 8 corners coordinates */
void plot_grid_with_id(CairoDrawing drawing, string input) {
  ifstream file;
  string dummy;
  int i;
  float val;
  vector<float> pos(8);

  file.open(input.c_str());
  if (!file.is_open()) {throw string("Unable to open input file ");}

  /* If first line is not a comment, go back */
  getline(file, dummy);
  if (dummy.find("#") == string::npos) {
    file.seekg (0, file.beg);
  }


  while (file.good()) {
    /* Skip first value */
    file >> val;
    for (i = 0; i < 8; ++i){
      file >> pos[i];
    }

    drawing.drawCell(pos, val);
  }
}


int main (int argc, char* argv[]){
  CairoDrawing drawing;
  string input, output;

  try {

    if ( argc == 2) {
      cout << "Setting output to grid.pdf\n" << endl;
      output = "grid.pdf";
    }
    else if ( argc == 3) 
      output = argv[2];
    else 
      throw "Usage: " + string(argv[0]) + " input output";

    input = argv[1];
    drawing.init(input, output);
    plot_grid_with_id(drawing, input);
    drawing.drawColormap();

    drawing.close();
  }
  catch (char const * str){ cout << "Error: " << str << endl; }
  catch(const string str) { cout << "Error: " << str << endl; }


}
