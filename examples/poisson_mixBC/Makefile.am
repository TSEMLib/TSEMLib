exampledir = $(datadir)/TSEMLib/examples/poisson_mixBC

dist_example_DATA = cube.mesh Makefile.sample  poisson.cpp  run.sh.sample  test.sh

install-data-hook:
	chmod 666 $(DESTDIR)$(exampledir)/Makefile.sample
	sed -i '/CPPFLAGS =/ c CPPFLAGS=-I$(prefix)/include -I$(AFEPackPath)/include -I$(dealIIPath)/include -I$(dealIIPath)/include/deal.II' $(DESTDIR)$(exampledir)/Makefile.sample
	sed -i '/LDFLAGS =/ c LDFLAGS=-L$(prefix)/lib -L$(AFEPackPath)/lib -L$(AFEPackPath)/lib' $(DESTDIR)$(exampledir)/Makefile.sample
	chmod 777 $(DESTDIR)$(exampledir)/run.sh.sample
	chmod 777 $(DESTDIR)$(exampledir)/test.sh
	sed -i '/export LD_LIBRARY_PATH/ c export LD_LIBRARY_PATH=$(prefix)/lib:$(AFEPackPath)/lib' $(DESTDIR)$(exampledir)/run.sh.sample
	sed -i '/export AFEPACK_PATH/ c export AFEPACK_PATH=$(AFEPackPath)/include/AFEPack' $(DESTDIR)$(exampledir)/run.sh.sample

