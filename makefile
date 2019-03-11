all:  docs install

docs:
	R -e "devtools::document()"
build:
	(R CMD build --no-build-vignettes trenaProjectErythropoiesis)

install:
	(cd ..; R CMD INSTALL --no-test-load trenaProjectErythropoiesis)

check:
	(cd ..; R CMD check `ls -t trenaProjectErythropoiesis) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

