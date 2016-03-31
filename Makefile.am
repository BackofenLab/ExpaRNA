# this is example-file: Makefile.am

# the subdirectories of the project to go into

SUBDIRS = src
EXTRA_DIST = Examples/*.fa Examples/*.align

run: $(EXECS)
	./src/ExpaRNA -v -e -a -i -o -L -s 2 Examples/16SrRNA_Ecoli_Ddiscoid.fa # 16SrRNA_Ecoli_Ddiscoid.mcs # HCVirus_IRES_RNAs.mcs # 16SrRNA_Ecoli_Ddiscoid.mcs

