EXTRA_DIST = analytical_case.f90 pangolingrid.py main.py

all:
	gfortran analytical_case.f90 -shared -fPIC -o libanalytical.so

doc:
	pydoc PangolinGrid
exec:
	python main.py

clean:
	rm -f *.pyc
