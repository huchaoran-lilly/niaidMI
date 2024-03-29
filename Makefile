objects := $(wildcard R/*.R) DESCRIPTION
dir := $(shell pwd)
version := $(shell grep "Version" DESCRIPTION | sed "s/Version: //")
pkg := $(shell grep "Package" DESCRIPTION | sed "s/Package: //")
tar := $(pkg)_$(version).tar.gz
tests := $(wildcard tests/testthat/*.R)
checkLog := $(pkg).Rcheck/00check.log

yr := $(shell date +"%Y")
dt := $(shell date +"%Y-%m-%d")


.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)


$(tar): $(objects)
	Rscript -e "library(methods); devtools::document();"
	R CMD build $(dir)

$(checkLog): $(tar) $(tests)
	R CMD check --as-cran $(tar)

.PHONY: newVersion
newVersion:
	@read -p "new version number: " NEWVER;\
	sed -i 's/^Version: [0-9]\.[0-9]\.[0-9]\.*[0-9]*[0-9]*[0-9]*[0-9]*/Version: '$$NEWVER'/' DESCRIPTION;\
	sed -i 's/Version: [0-9]\.[0-9]\.[0-9]\.*[0-9]*[0-9]*[0-9]*[0-9]*/Version: '$$NEWVER'/' README.md

	@echo "NEWS.md and cran-comment should be modified manually."

.PHONY: clean
clean:
	rm -rf *~ */*~ *.Rhistroy *.tar.gz *.Rcheck/ .\#* src/*.so src/*.o vignettes/*.html
	rm -rf config.log config.status src/Makevars autom4te.cache/
	rm -rf tests/testthat/_snaps CRAN-RELEASE