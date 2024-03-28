headerincludedir = $(includedir)

headerinclude_HEADERS = Correspondence.h\
			Multiindex.h\
			TetrahedralSEM.h\
			TSEMSolver.h\
			TSEMVisualize.h

install-data-hook:
	@sed -i 's|\.\./\.\./include|$(prefix)/include|g' $(prefix)/include/TetrahedralSEM.h
	cp -r quad_info $(prefix)/include/
