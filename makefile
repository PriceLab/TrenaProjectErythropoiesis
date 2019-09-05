all:  docs install

docs:
	R -e "devtools::document()"
build:
	(R CMD build --no-build-vignettes TrenaProjectErythropoiesis)

install:
	(cd ..; R CMD INSTALL --no-test-load TrenaProjectErythropoiesis)

check:
	(cd ..; R CMD check `ls -t TrenaProjectErythropoiesis) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

