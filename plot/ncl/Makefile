NCL=/usr/local/ncl/bin/ncl

# Generate all contours for gaussian hills
gaussian:
	make create_contour TIME=0 TIME2=0 TITLE=0
	make create_contour TIME=0 TIME2=T2 TITLE="T\/2"
	make create_contour TIME=T2 TIME2=T TITLE=T

# Important, must use double quote for splitting string !
# Substitute variables
create_contour:
	perl -pi -e '/ratio_160lat_${TIME}/ && s/lat_${TIME}/lat_${TIME2}/' plot_contour.ncl
	perl -pi -e "/MainString/ && s/= \"Gaussian.*$$/= \"Gaussian hills, ${TITLE} \
	  (0.56x0.37:F34:0)\"/" plot_contour.ncl
	${NCL} plot_contour.ncl
	perl rotate_pdf.pl ratio_160lat_${TIME2}.pdf ratio_160lat_${TIME2}.pdf

install:
	perl -MCPAN -e 'install PDF::API2'    
