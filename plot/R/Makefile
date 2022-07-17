all: cv_rates filaments

# Generate two figures for convergence rate
cv_rates: cv_rate1 cv_rate2

cv_rate1: SRC=cv_rate.ps
cv_rate1: FIRST=TRUE

cv_rate2: SRC=cv_rate2.ps
cv_rate2: FIRST=FALSE

# Plot either one of the 2 possible plots
cv_rate1 cv_rate2: 
	@echo "Plotting ${SRC}"; \
	perl -pi -e 's/first_plot <-.*$$/first_plot <- ${FIRST}/' plot_cv_rate.r ; \
	perl -pi -e 's/fname <-.*$$/fname <- "${SRC}"/' plot_cv_rate.r ; \
	R < plot_cv_rate.r --no-save > /dev/null; \
	make -s SRC=${SRC} pdf

# Generate two figures for filaments diagnostics
filaments: filament1 filament2

filament1: SRC=filament.ps
filament1: FIRST=TRUE

filament2: SRC=filament2.ps
filament2: FIRST=FALSE

# Plot either one of the 2 possible plots
filament1 filament2: 
	@echo "Plotting ${SRC}"; \
	perl -pi -e 's/first_plot <-.*$$/first_plot <- ${FIRST}/' plot_filaments.r ; \
	perl -pi -e 's/fname <-.*$$/fname <- "${SRC}"/' plot_filaments.r ; \
	R < plot_filaments.r --no-save > /dev/null; \
	make -s SRC=${SRC} pdf

# Convert a postscript to pdf (with correct boundaries)
BASE=$(SRC:.ps=)
DEST=$(SRC:.ps=.pdf)
pdf:
	@echo "Converting ${SRC} to pdf"; \
	test  -f ${SRC} || { echo "File does not exist"; exit 1; }; \
	ps2pdf -dEPSCrop ${SRC} ; \
	pdfcrop ${DEST} ; \
	echo "Renaming..."; \
	mv ${BASE}"-crop.pdf" ${DEST}; \
