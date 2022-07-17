# Html generated with doxygen (default)
all: html

html:
	perl -pi -e 's/(?<=GENERATE_LATEX)\s*=\s*\w+/ = NO/' doxygen.conf
	doxygen doxygen.conf > /dev/null

pdf:
	perl -pi -e 's/(?<=GENERATE_LATEX)\s*=\s*\w+/ = YES/' doxygen.conf
	doxygen doxygen.conf > /dev/null
	cd doxygen_latex && pdflatex refman.tex
	mv doxygen_latex/refman.pdf pangolin_ref.pdf
