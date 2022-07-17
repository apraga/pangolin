CFLAGS=-Wall -pedantic -g -I/usr/X11R6/include `pkg-config --cflags cairo`
LDFLAGS=`pkg-config --libs cairo` -L/usr/X11R6/lib -lX11 -lm
CPPFLAGS=-Wall -g `pkg-config --cflags cairo`
CXX = g++

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

all: plot_grid

exec: 
	@echo "###############################################################################"
	@echo "# Reading ../output/grid.dat and output plot in grid.pdf"
	@echo "###############################################################################"
	./plot_grid ../../output/grid_init.dat grid.pdf

plot_grid: $(OBJ)
	$(CXX) -o $@ ${LDFLAGS} $^

#plot_grid.o: parser.o
#plot_grid.o: parser.h
#plot_grid.o: drawing.h
#drawing.o: colormap.h

clean:
	rm *.o
