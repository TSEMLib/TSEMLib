/*
 * realization of tetrahedral spectral element method
 *
 * used correspondence, which is declared in Correspondence.h
 * require Unitary_Multiindex and Zero_Multiindex defined in Multiindex.h
 *     and been initialized before using TSEM
 */
#ifndef __TETRAHEDRALSEM_
#define __TETRAHEDRALSEM_

#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "AFEPack/FEMSpace.h"
#include "AFEPack/Geometry.h"
#include "lac/sparsity_pattern.h"
#include "lac/sparse_matrix.h"
#include "lac/vector.h"
#include "Multiindex.h"
#include "Correspondence.h"

#include <iomanip>
//#include "Option.h"

// #define VALUETYPE_ITERATOR_TSEM double

template <typename valuetype> class TSEM
{
protected:
    const valuetype tol_zero = 1.0e-12;
    Correspondence<3> correspondence;
    Multiindex<3> Unitary_Multiindex[3];
    Multiindex<3> Zero_Multiindex;
public:
    int M, n_index; // polynomial order, number of polynomials

public:
    // int n_quad_accuracy[3];
    std::vector<int> n_q_point; // [i]: number of quadrature point in i dimensional rule
    std::vector<AFEPack::Point<3> > QPoint; // 3-d quadrature point
    std::vector<std::vector<std::vector<valuetype> > > QPoint_Barycentric; // 1 & 2-d quadrature point in barycenter coordinate
    std::vector<std::vector<valuetype> > Weight; // 1, 2 & 3-d quadrature weight

protected:
    int n_dof_edge, n_dof_face;
    // number of nonzero contribution on edges
    std::vector<std::vector<std::vector<int> > > n_expr_edge;
    // index of nonzero contribution on edges
    std::vector<std::vector<std::vector<std::vector<int> > > > ind_expr_edge;
    // value of nonzero contribution on edges
    std::vector<std::vector<std::vector<std::vector<valuetype> > > > val_expr_edge;
    // number of nonzero contribution on face
    std::vector<std::vector<int> > n_expr_face;
    // index of nonzero contribution on face
    std::vector<std::vector<std::vector<int> > > ind_expr_face;
    // value of nonzero contribution on face
    std::vector<std::vector<std::vector<valuetype> > > val_expr_face;

    std::vector<std::vector<AFEPack::Point<3> > > point_ref_mesh; // [i = 0:3], barycenter of i dimensional geometry
public:
    int n_dof_geometry[4];
    int n_geometry[4];
    int n_element, n_geometry_total, n_dof_total;
    // number_node[ind+1=1:2][order+1=1:ind+1][]: the node number of ind+1-dimensional geometry
    std::vector<std::vector<std::vector<int> > > number_node;
    std::vector<valuetype> val_volume;
protected:
public:
    std::vector<std::vector<int> > number_edge;
    // flag that whether order of start/end point on global face is the same as that on global edge
    std::vector<std::vector<bool> > flag_sameorder_edgeonface;
    // type of projection for edge and face on each element
    std::vector<std::vector<std::vector<int> > > type_projection;
    std::vector<std::vector<int> > type_projection_inv; // inverse projection type, for face on each element
    std::vector<std::vector<std::vector<int> > > index_geometry_onelement;

    // int n_dof;
    /* {weight_local, geometry_dimension, geometry_order}:
     *    weight_local:       a number for sorting, read from fem_space, determine the order of all dimensional geomtry
     *    geometry_dimension: the dimension of geometry corresponds to this set
     *    geometry_order:     the order of geometry in the same dimensional ones corresponds to this set
     */
    std::vector<valuetype> weight_location; // the weight_location of 0: 3 dimensional geometry in turns
    std::vector<int> geometry_dimension;
    std::vector<int> geometry_order; // [i = 0:n_geometry_total-1] the order of i-th entry in weight_location, whose dimension is geometry_dimension[i]
    // correspondence between fem dof/element and geometry
    // transform_femdof2geometry: index from fem dof to 0 & 1 dimensional geometry, for 1 dimensional geometry, its index plus mesh.n_geometry(0)
    std::vector<int> transform_femdof2geometry;
    // location_geometry: location of all geometry (0:3 dimensional) according to increasing order of weight_location
    std::vector<int> location_geometry;
    // location_actualdof: start index of geometry in actual discretized matrix
    std::vector<int> location_actualdof;
public:
    // expression of local coefficient by linear summation of global ones
    std::vector<std::vector<int> >                     transform_n_global2local;
    std::vector<std::vector<std::vector<int> > >       transform_ind_global2local;
    std::vector<std::vector<std::vector<valuetype> > > transform_val_global2local;
protected:
    // expression of global coefficient on each-dimensional geometry by local ones, equivalent to column compress of transformation matrix
    std::vector<std::vector<int> >                     transform_n_local2global;
    std::vector<std::vector<std::vector<int> > >       transform_ind_local2global;
    std::vector<std::vector<std::vector<valuetype> > > transform_val_local2global;
    
    std::vector<std::vector<std::vector<valuetype> > > basis_value; // [ind = 0:2][p = 0:n_q_point[ind]-1][]: basis_value for ind+1 dimensional quadrature region
    std::vector<std::vector<std::vector<valuetype> > > basis_value_interp; // basis value for interpolation
    std::vector<std::vector<std::vector<valuetype> > > basis_value_addition; // additional basis value for interpolation on face
    std::vector<std::vector<std::vector<valuetype> > > basis_gradient; // [p = 0:n_q_point[2]-1][ind_index = 0:n_index-1][3]: gradient of basis function on quadrature points

    std::vector<int> n_transform_local;
    std::vector<std::vector<int> > transform_local;
    std::vector<std::vector<valuetype> > weight_transform_local;
public:
    std::vector<std::vector<valuetype> > basis_value_actual; // [p = 0:n_q_point[ind]-1][]: basis_value for 3 dimensional quadrature region
    std::vector<std::vector<std::vector<valuetype> > > basis_gradient_actual; // actual gradient of basis on quadrature points

protected:
    std::vector<int> n_index_variation; // number of the coefficient for the derivative of generalized jacobi polynomial
    std::vector<std::vector<Multiindex<3> > > index_variation; // variation of multiindex corresponding to the coefficients
    std::vector<std::vector<std::vector<valuetype> > > coefficient_derivative;
    std::vector<SparsityPattern> sp_stiff_matrix_element_actual;
    std::vector<SparseMatrix<valuetype> > stiff_matrix_element_actual;
public:
    std::vector<unsigned int> n_nonzero_per_row_stiff;
    SparsityPattern sp_stiff_matrix;

    SparsityPattern sp_mass_matrix_element_actual;
    SparseMatrix<valuetype> mass_matrix_element_actual;

    std::vector<unsigned int> n_nonzero_per_row_mass;
    SparsityPattern sp_mass_matrix;
    SparsityPattern sp_mass_matrix_essential;

public:
    std::vector<std::vector<unsigned int> > flag_bm; // boundary flag for all dimensional geometry
    std::vector<unsigned int> flag_bm_dof;
    valuetype Rmin, Rmax;
    std::vector<std::vector<unsigned int> > flag_mask; // 0: inner; 1: transition; 2: outer

protected:
    std::vector<std::vector<unsigned int> > conversion;

    std::vector<int> n_quadPoint{20, 130, 923};
    std::vector<std::string> path_quadInfo{"../../include/quad_info/quad_info_1d_p38.dat",
					   "../../include/quad_info/quad_2d_deg26_n130.dat",
					   "../../include/quad_info/quad_3d_deg26_n923.dat"};
    std::string path_interpInfo = "../../include/quad_info/Info_Interpolate_More";
    std::string path_1DQuadInfo = "../../include/quad_info/qinfo_1d/";

protected:
    valuetype calc_generalized_jacobi_polynomial(int alpha, int beta, int k, valuetype x);
    valuetype calc_coefficient_a(int alpha, int beta, int ind, int k);
    valuetype calc_coefficient_b(int alpha, int beta, int ind, int k);
    valuetype calc_coefficient_c(int alpha, int beta, int ind, int k);
    valuetype calc_coefficient_d(int alpha, int beta, int k);
    valuetype calc_coefficient_e(int alpha, int beta, int ind, int k);
    valuetype calc_coefficient_g(int alpha, int beta, int ind, int k);
    valuetype calc_coefficient_rho(Multiindex<3> index);
    valuetype calc_coefficient_kappa(Multiindex<3> index);
    valuetype calc_coefficient_theta(Multiindex<3> index);
    valuetype calc_coefficient_D(int ind_derivative, int ind_variation, Multiindex<3> index);
    int calc_delta(int i, int j);
    int calc_binomial_coefficient(int n, int m);
    AFEPack::Point<3> calc_cross_product(AFEPack::Point<3> p1, AFEPack::Point<3> p2);
    double calc_inner_product(AFEPack::Point<3> p1, AFEPack::Point<3> p2);
    
    bool is_on_boundary(AFEPack::Point<3> &p, valuetype bnd_left, valuetype bnd_right);
    bool is_on_boundary(AFEPack::Point<3> &p, valuetype radius);
    valuetype calc_length_line(AFEPack::Point<3> &p0, AFEPack::Point<3> &p1);
    valuetype calc_area_triangle(AFEPack::Point<3> &p0, AFEPack::Point<3> &p1, AFEPack::Point<3> &p2);
public:
    valuetype calc_volume_tetrahedron(AFEPack::Point<3> &p0, AFEPack::Point<3> &p1, AFEPack::Point<3> &p2, AFEPack::Point<3> &p3);
public:
    TSEM(){};
    TSEM(int polynomial_order,
	 std::vector<int> &n_quadrature_point, std::vector<std::string> &quad_filename,
	 const std::string &interp_filename,
	 RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space){
	init(polynomial_order,
	     n_quadrature_point, quad_filename,
	     interp_filename,
	     mesh, fem_space);	       
    };
    TSEM(int polynomial_order,
	 std::vector<int> &n_quadrature_point, std::vector<std::string> &quad_filename,
	 RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space){
	init(polynomial_order,
	     n_quadrature_point, quad_filename,
	     path_interpInfo,
	     mesh, fem_space);	       
    };
    TSEM(int polynomial_order, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space){
	init(polynomial_order, mesh, fem_space);
    };
    TSEM(int polynomial_order, Mesh<3> &mesh, FEMSpace<double, 3> &fem_space){
	init(polynomial_order, mesh, fem_space);
    };
    TSEM(int polynomial_order, int qinfo_order, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space){
	init(polynomial_order, qinfo_order, mesh, fem_space);
    };

    void init(int polynomial_order,
	      std::vector<int> &n_quadrature_point, std::vector<std::string> &quad_filename,
	      const std::string &interp_filename,
	      RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void init(int polynomial_order, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void init(int polynomial_order, int qinfo_order, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void init(int polynomial_order, Mesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void init_lazy(int polynomial_order,
		   std::vector<int> &n_quadrature_point, std::vector<std::string> &quad_filename,
		   const std::string &interp_filename,
		   RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void init_lazy(int polynomial_order, int qinfo_order,
		   RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space);

    SparseMatrix<valuetype> stiff_matrix;
    SparseMatrix<valuetype> mass_matrix;
    SparseMatrix<valuetype> mass_matrix_essential;

    void setup();
    void read_parameter(int polynomial_order);
    void read_quad_info(std::vector<int> &n_quadrature_point, std::vector<std::string> &quad_filename);
    void build_quad_info(int qinfo_order); // build quadrature info from 1d qinfo
    void read_interp_info(const std::string &interp_filename);
    void build_geometry_info(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void build_geometry_info(Mesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void build_basis_value();
    void build_local_transform();
    void build_global_transform(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void build_global_transform(Mesh<3> &mesh, FEMSpace<double, 3> &fem_space);
    void build_stiff_matrix_init();
    void build_stiff_matrix(FEMSpace<double, 3> &fem_space);
    void build_mass_matrix_init();
    void build_mass_matrix(FEMSpace<double, 3> &fem_space);
    void build_mass_V_matrix(std::vector<unsigned int> &n_nonzero_per_row_mass_V, SparsityPattern &sp_pattern, SparseMatrix<valuetype> &sp_matrix,
			     std::vector<std::vector<valuetype> > &value_V);
    void build_mass_V_matrix(std::vector<unsigned int> &n_nonzero_per_row_mass_V, SparsityPattern &sp_pattern, SparseMatrix<valuetype> &sp_matrix,
			     FEMSpace<double, 3> &fem_space, valuetype(*func)(double*));
    
    void calc_rhs(Vector<valuetype> &rhs, std::vector<std::vector<valuetype> > &value_f);
    void calc_rhs(Vector<valuetype> &rhs, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*));
    void calc_rhs(Vector<valuetype> &rhs, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*, double), double t);
    void build_flag_bm(valuetype bnd_left, valuetype bnd_right);
    void build_flag_bm(valuetype radius);
    void build_flag_bm(valuetype bnd_left, valuetype bnd_right, valuetype axis[]);
    void build_flag_bm(RegularMesh<3> &mesh);
    void build_flag_mask(RegularMesh<3> &mesh, valuetype Rinner, valuetype Router);
    void impose_zero_boundary_condition(SparseMatrix<valuetype> &sp_matrix, unsigned int flag_boundary = 1);
    void impose_zero_boundary_condition_complex(SparseMatrix<valuetype> &sp_matrix, unsigned int flag_boundary = 1);
    void impose_zero_boundary_condition(Vector<valuetype> &rhs, unsigned int flag_boundary = 1);
    void impose_zero_boundary_condition_complex(Vector<valuetype> &rhs, unsigned int flag_boundary = 1);
    void impose_boundary_condition_rowOnly(SparseMatrix<valuetype> &sp_matrix, unsigned int flag_boundary = 1);
    void impose_boundary_condition_rowOnly(SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &rhs, std::vector<std::vector<std::vector<valuetype> > > &value_bnd, unsigned int flag_boundary = 1);
    void impose_boundary_condition(SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &rhs, std::vector<std::vector<std::vector<valuetype> > > &value_bnd,
				   bool flag_modify_matrix, unsigned int flag_boundary = 1);
    void impose_boundary_condition(SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &rhs,
				   RegularMesh<3> &mesh, valuetype(*func)(double*),
				   bool flag_modify_matrix, unsigned int flag_boundary = 1);
    void impose_boundary_condition(SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &rhs,
				   Mesh<3>& mesh, valuetype(*func)(double*),
				   bool flag_modify_matrix, unsigned int flag_boundary = 1);
    void impose_boundary_condition(Vector<double>& rhs, RegularMesh<3>& mesh, valuetype(*func)(double*), unsigned int flag_boundary = 2); // impose nuemann boundary condition
    void impose_mask(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*), std::vector<Vector<valuetype> > &psi);

    void calc_coef_onElement(Vector<valuetype> &src, std::vector<valuetype> &dst, unsigned int ind_ele);
    void calc_val_qp_onElement(Vector<valuetype> &src, std::vector<valuetype> &dst, unsigned int ind_ele);
    void calc_val_gradient_qp_onElement(Vector<valuetype> &src, std::vector<std::vector<double> > &dst, RegularMesh<3> &mesh, unsigned int idx_ele);
    void calc_basis_gradient(Vector<valuetype> &src, std::vector<std::vector<std::vector<double> > > &dst, RegularMesh<3> &mesh, unsigned int idx_ele);

    valuetype calc_l2_error(Vector<valuetype> &sol, std::vector<std::vector<valuetype> > &val_exc);
    valuetype calc_l2_error(Vector<valuetype> &sol, FEMSpace<double, 3> &fem_space, valuetype (*func)(double *));
    valuetype calc_l2_error(Vector<valuetype> &sol, FEMSpace<double, 3> &fem_space, valuetype (*func)(double *, double), valuetype t);
    valuetype calc_l2_error_gradient(RegularMesh<3> &mesh,
				     Vector<valuetype> &sol, std::vector<std::vector<std::vector<valuetype> > > &val_g_exc);
    valuetype calc_l2_error_gradient(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space,
				     Vector<valuetype> &sol, void (*func)(double*, valuetype*));
    valuetype calc_l2_error_gradient(Mesh<3> &mesh,
				     Vector<valuetype> &sol, std::vector<std::vector<std::vector<valuetype> > > &val_g_exc);
    valuetype calc_l2_error_gradient(Mesh<3> &mesh,
				     Vector<valuetype> &sol, std::vector<std::vector<valuetype> > &val_g_exc, unsigned int idx_ele);
    valuetype calc_h1_error(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space,
			    Vector<valuetype> &sol, valuetype (*func)(double*), void (*func_g)(double*, valuetype*));
    valuetype calc_l2_error_gradient_difference(RegularMesh<3> &mesh,
						Vector<valuetype> &u, Vector<valuetype> &v);
    valuetype calc_l2_density_difference(Vector<valuetype> &u, std::vector<std::vector<valuetype> > &v);
    valuetype calc_l2_density_difference(Vector<valuetype> &u, Vector<valuetype> &v);
    valuetype calc_l2_density_difference(std::vector<Vector<valuetype> > &psi, FEMSpace<double, 3> &fem_space, double(*psi_re)(double *, double), double(*psi_im)(double *, double), double t);
    valuetype calc_l2_difference(Vector<valuetype> &u, Vector<valuetype> &v);
    valuetype calc_l2_difference(Vector<valuetype> &u, Vector<valuetype> &v, int polynomial_order, int flag);
    valuetype calc_h1_difference(RegularMesh<3> &mesh, Vector<valuetype> &u, Vector<valuetype> &v);
    valuetype calc_l2_error_component(Vector<valuetype> &u, int polynomial_order_less, int polynomial_order_more);
    valuetype calc_density_l2_difference(std::vector<Vector<valuetype> > &u, std::vector<Vector<valuetype> > &v, std::vector<valuetype> &n_occupation);

    void calc_interpolation(Vector<valuetype> &interp, std::vector<std::vector<std::vector<valuetype> > > &val_interp);
    void calc_interpolation(Vector<valuetype> &interp, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*));
    void calc_interpolation(Vector<valuetype> &interp, Mesh<3> &mesh, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*));
    void calc_interpolation(Vector<valuetype> &interp, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*, double), double t);
    void calc_interpolation(std::vector<Vector<valuetype> > &psi, std::vector<valuetype> &n_occ, RegularMesh<3> &mesh, const std::string &filename_interpMesh, std::string &filename_out);

    void calc_density(std::vector<Vector<valuetype> > &src, std::vector<valuetype> &dst, std::vector<valuetype> &n_occupation, unsigned int idx_ele);
    void calc_density(std::vector<Vector<valuetype> > &src, std::vector<std::vector<valuetype> > &dst, std::vector<valuetype> &n_occupation);

    valuetype calc_dipole(std::vector<Vector<valuetype> > &psi, FEMSpace<double, 3> &fem_space, int idx_dim);

    void read_coef(Vector<valuetype> &dst, std::string filename);
    void read_coef(std::vector<Vector<valuetype> > &dst, std::string filename);
    void write_coef(Vector<valuetype> &src, std::string filename);
    void write_coef(std::vector<Vector<valuetype> > &src, std::string filename);

    void calc_sum_dof(std::vector<unsigned int> &sum_dof);

    // void calc_val_qpoint(Vector<valuetype> &src, int ind_element, std::vector<valuetype> &dst);

    void calc_dependence(RegularMesh<3> &mesh, std::vector<AFEPack::Point<3> > &point, std::vector<int> &n_dep, std::vector<std::vector<int> > &dep_idx, std::vector<std::vector<double> > &dep_para);
    void calc_dependence(RegularMesh<3> &mesh, std::vector<AFEPack::Point<3> > &point, std::vector<int> &idx_ele, std::vector<std::vector<double> > &dep_basis, std::vector<std::vector<std::vector<double> > > &dep_gradient);
    void calc_dependence(RegularMesh<3> &mesh, int np, std::vector<std::vector<int> > &idx_eleInCube, std::vector<AFEPack::Point<3> > &point, std::vector<int> &idx_ele, std::vector<std::vector<double> > &dep_basis, std::vector<std::vector<std::vector<double> > > &dep_gradient);

    // bool is_on_element(RegularMesh<3> &mesh, int ind_dim, int ind_geo, AFEPack::Point<3> &pos);
    valuetype calc_val_inElement(RegularMesh<3> &mesh, int ind_ele, AFEPack::Point<3> &pos, Vector<valuetype> &psi);
    valuetype calc_val_inElement(RegularMesh<3> &mesh, int ind_ele, AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi);
    valuetype calc_val_inElement(RegularMesh<3> &mesh, int ind_ele, AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi, std::vector<valuetype>& n_occupation);
    valuetype calc_val_point(RegularMesh<3> &mesh, AFEPack::Point<3> &pos, Vector<valuetype> &u);
    valuetype calc_val_point(Mesh<3> &mesh, AFEPack::Point<3> &pos, Vector<valuetype> &u);
    valuetype calc_val_point(RegularMesh<3> &mesh, AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi, std::vector<valuetype> &n_occupation);
    valuetype calc_val_point(Mesh<3> &mesh, AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi, std::vector<valuetype> &n_occupation);
    valuetype calc_val_point_insideElement(RegularMesh<3> &mesh, int ind_ele,
    					   AFEPack::Point<3> &pos, Vector<valuetype> &u);
    valuetype calc_val_point_insideElement(RegularMesh<3> &mesh, int ind_ele,
    					   AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi, std::vector<valuetype> &n_occupation);
};


//#define TEMPLATE_TSEM template<typename valuetype>
//#define THIS_TSEM TSEM<valuetype>

// TEMPLATE_TSEM
// void THIS_TSEM::calc_val_qpoint(Vector<valuetype> &src, int ind_element, std::vector<valuetype> &dst)
// { // recover function value on ind_element-th element in fem_space, generate dst from whole Vector src
//     std::vector<valuetype> coef_local(n_index, 0);
//     for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
// 	for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_element][ind_index]; ++ind_nnz)
// 	    coef_local[ind_index] += src(transform_ind_global2local[ind_element][ind_index][ind_nnz])
// 		* transform_val_global2local[ind_element][ind_index][ind_nnz];

//     // dst.resize(n_q_point[2]);
//     for (unsigned int p = 0; p < n_q_point[2]; ++p){
// 	valuetype val_tmp = 0;
// 	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
// 	    val_tmp += coef_local[ind_index] * basis_value_actual[p][ind_index];
// 	dst[p] = val_tmp;
//     }
// }

// Explicit instantiation for the types you want to use
template class TSEM<double>;
template class TSEM<float>;

#endif
