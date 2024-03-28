/*
 * solve Poisson problem on \Omega = [0, 1]^3 with mixed boundary conditions
 *     -\nabla^2 u = f,
 *     (\partial u/\partial n)_{B_1} = g, u_{B_2} = u_{ext}
 *         where B_1 = \partial\Omega\cap\{(x,y,z):0<x,y<1,z=0\}. B_2 = \partial\Omega\backslash B_1
 */

#include <iostream>
#include <iomanip>
#include <math.h>

#include "AFEPack/Geometry.h"
#include "AFEPack/FEMSpace.h"

#include "lac/sparsity_pattern.h"
#include "lac/sparse_matrix.h"
#include "lac/vector.h"

#include "TetrahedralSEM.h"
#include "TSEMSolver.h"

#define DIM 3
#define PI (4.0*atan(1.0))

double u(double* p){ return cos(PI*p[0])*sin(1.5*PI*p[1])*sin(0.25*PI+2*PI*p[2]); }
double f(double* p){ return (7.25*pow(PI,2)) * u(p); }
double g(double* p){ return -2*PI*cos(PI*p[0])*sin(1.5*PI*p[1])*cos(0.25*PI+2*PI*p[2]); }
void calc_gradient(double* p, std::vector<double>& dst){
    double& x = p[0], y = p[1], z = p[2];
    dst[0] =     -PI * sin(PI*x) * sin(1.5*PI*y) * sin(0.25*PI+2*PI*z);
    dst[1] =  1.5*PI * cos(PI*x) * cos(1.5*PI*y) * sin(0.25*PI+2*PI*z);
    dst[2] =    2*PI * cos(PI*x) * sin(1.5*PI*y) * cos(0.25*PI+2*PI*z);
}

int main(int argc, char *argv[])
{
    if (argc != 4){
	std::cout << "Usage: " << argv[0] << " filename_mesh polynomial_order n_global_refinement\n";
	return 1;
    }
    

    // read mesh from .mesh file
    HGeometryTree<DIM> h_tree; h_tree.readMesh(argv[1]);
    IrregularMesh<DIM> *irregular_mesh = new IrregularMesh<DIM>;
    irregular_mesh->reinit(h_tree);
    if (atoi(argv[3]) > 0) irregular_mesh->globalRefine(atoi(argv[3]));
    irregular_mesh->semiregularize();
    irregular_mesh->regularize(false);
    RegularMesh<DIM> &mesh = irregular_mesh->regularMesh();
    
    // generate 3D template
    TemplateGeometry<DIM> template_geometry;
    CoordTransform<DIM, DIM> coord_transform;
    TemplateDOF<DIM> template_dof;
    BasisFunctionAdmin<double, DIM, DIM> basis_function;
    template_geometry.readData("tetrahedron.tmp_geo");
    coord_transform.readData("tetrahedron.crd_trs");
    template_dof.reinit(template_geometry); template_dof.readData("tetrahedron.2.tmp_dof");
    basis_function.reinit(template_dof);    basis_function.readData("tetrahedron.2.bas_fun");
    std::vector<TemplateElement<double, DIM, DIM> > template_element(1);
    template_element[0].reinit(template_geometry, template_dof, coord_transform, basis_function);
    // construct finite element space
    FEMSpace<double, DIM> fem_space(mesh, template_element);
    int n_element = mesh.n_geometry(DIM);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i) fem_space.element(i).reinit(fem_space, i, 0);
    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();

    // construct tsem space, with quadrature info order 30 in each direction
    int M = atoi(argv[2]);
    TSEM<double> tsem(M, 30, mesh, fem_space);
    tsem.build_flag_bm(mesh);
    int n_dof_total = tsem.n_dof_total;

    // build linear system
    SparseMatrix<double> sp_matrix(tsem.stiff_matrix.get_sparsity_pattern());
    sp_matrix.copy_from(tsem.stiff_matrix);
    Vector<double> rhs(n_dof_total);
    tsem.calc_rhs(rhs, fem_space, &f);
    tsem.impose_boundary_condition(rhs, mesh, &g, 2);
    tsem.impose_boundary_condition(sp_matrix, rhs, mesh, &u, true, 1);

    // solve
    TSEMSolver<double> solver;
    Vector<double> solution(n_dof_total);
    solver.solve_PCG(sp_matrix, solution, rhs, 1.0e-14, 100000);

    // calculate error
    Vector<double> interpolation(n_dof_total); tsem.calc_interpolation(interpolation, mesh, fem_space, &u);
    std::vector<std::vector<std::vector<double> > > value_gradient_exact(n_element, std::vector<std::vector<double> > (tsem.n_q_point[2], std::vector<double> (3)));
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
    	std::vector<AFEPack::Point<3> > q_point = fem_space.element(ind_ele).local_to_global(tsem.QPoint);
    	for (int p = 0; p < tsem.n_q_point[2]; ++p) calc_gradient(q_point[p], value_gradient_exact[ind_ele][p]);
    }
    double err_l2  = tsem.calc_l2_error(solution, fem_space, &u);
    double err_l2g = tsem.calc_l2_error_gradient(mesh, solution, value_gradient_exact);
    double erri_l2  = tsem.calc_l2_error(interpolation, fem_space, &u);
    double erri_l2g = tsem.calc_l2_error_gradient(mesh, interpolation, value_gradient_exact);
    std::cout.precision(16);
    std::cout << "n_element = " << n_element
	      << ", M = " << M << ", n_dof = " << n_dof_total << '\n'
	      << "\terr_l2  = " << err_l2 << ",\terr_l2g  = " << err_l2g << '\n'
	      << "\terri_l2 = " << erri_l2 << ",\terri_l2g = " << erri_l2g
	      << '\n';

    
    return 0;
}
