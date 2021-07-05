.PHONY: all docs build check

all: for-commit

for-commit: docs build check

docs:
	rm -rf NAMESPACE man docs
	R -e 'library(devtools); document()'
	R -e 'pkgdown::build_site()'

build:
	rm -f ../chameleon_*.tar.gz
	mv Makefile ../Makefile.chameleon
	test ! -d docs || mv docs ../docs.chameleon
	cd .. && R CMD build --resave-data chameleon
	mv ../Makefile.chameleon Makefile
	test ! -d ../docs.chameleon || mv ../docs.chameleon docs

check:
	R CMD check --as-cran ../chameleon_*.tar.gz
