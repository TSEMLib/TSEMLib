#include "../include/TetrahedralSEM.h"
// This file contains functions for initialization
// Functions:
// init()


#define TEMPLATE_TSEM template<typename valuetype>
#define THIS_TSEM TSEM<valuetype>


TEMPLATE_TSEM
void THIS_TSEM::init(int polynomial_order,
		     std::vector<int> &n_quadrature_point, std::vector<std::string> &quad_filename,
		     const std::string &interp_filename,
		     RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    setup();
    read_parameter(polynomial_order);
    read_quad_info(n_quadrature_point, quad_filename);
    read_interp_info(interp_filename);
    build_geometry_info(mesh, fem_space);
    build_basis_value();
    build_local_transform();
    build_global_transform(mesh, fem_space);
    build_stiff_matrix_init();
    build_stiff_matrix(fem_space);
    build_mass_matrix_init();
    build_mass_matrix(fem_space);
};

TEMPLATE_TSEM
void THIS_TSEM::init(int polynomial_order, Mesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    setup();
    read_parameter(polynomial_order);
    read_quad_info(n_quadPoint, path_quadInfo);
    read_interp_info(path_interpInfo);
    build_geometry_info(mesh, fem_space);
    build_basis_value();
    build_local_transform();
    build_global_transform(mesh, fem_space);
    build_stiff_matrix_init();
    build_stiff_matrix(fem_space);
    build_mass_matrix_init();
    build_mass_matrix(fem_space);
};

TEMPLATE_TSEM
void THIS_TSEM::init(int polynomial_order, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    setup();
    read_parameter(polynomial_order);
    read_quad_info(n_quadPoint, path_quadInfo);
    read_interp_info(path_interpInfo);
    build_geometry_info(mesh, fem_space);
    build_basis_value();
    build_local_transform();
    build_global_transform(mesh, fem_space);
    build_stiff_matrix_init();
    build_stiff_matrix(fem_space);
    build_mass_matrix_init();
    build_mass_matrix(fem_space);
;}

TEMPLATE_TSEM
void THIS_TSEM::init(int polynomial_order, int qinfo_order, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    setup();
    read_parameter(polynomial_order);
    build_quad_info(qinfo_order);
    read_interp_info(path_interpInfo);
    build_geometry_info(mesh, fem_space);
    build_basis_value();
    build_local_transform();
    build_global_transform(mesh, fem_space);
    build_stiff_matrix_init();
    build_stiff_matrix(fem_space);
    build_mass_matrix_init();
    build_mass_matrix(fem_space);

};

TEMPLATE_TSEM
void THIS_TSEM::init_lazy(int polynomial_order, int qinfo_order,
			  RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    setup();
    read_parameter(polynomial_order);
    build_quad_info(qinfo_order);
    read_interp_info(path_interpInfo);
    build_geometry_info(mesh, fem_space);
    build_basis_value();
    build_local_transform();
    build_global_transform(mesh, fem_space);
};

TEMPLATE_TSEM
void THIS_TSEM::init_lazy(int polynomial_order,
			  std::vector<int> &n_quadrature_point, std::vector<std::string> &quad_filename,
			  const std::string &interp_filename,
			  RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    setup();
    read_parameter(polynomial_order);
    read_quad_info(n_quadrature_point, quad_filename);
    read_interp_info(interp_filename);
    build_geometry_info(mesh, fem_space);
    build_basis_value();
    build_local_transform();
    build_global_transform(mesh, fem_space);
};

TEMPLATE_TSEM
void THIS_TSEM::setup()
{
    // setup for multiindex
    for (unsigned int i = 0; i < 3; ++i){
	for (unsigned int j = 0; j < 3; ++j) Unitary_Multiindex[i].index[j] = 0;
	Unitary_Multiindex[i].index[i] = 1;
	Zero_Multiindex.index[i] = 0;
    }
}

TEMPLATE_TSEM
void THIS_TSEM::read_parameter(int polynomial_order)
{
    M = polynomial_order;
    correspondence.init(polynomial_order);
    n_index = correspondence.n_index();
    std::cerr << "read coefficient, n_index = " << n_index << '\n';
    std::cerr << "set polynomial order for TSEM, M = " << M << '\n';
}

#undef TEMPLATE_TSEM
#undef THIS_TSEM
