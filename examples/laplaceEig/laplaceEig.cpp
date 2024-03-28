/*
 * solve Laplace eigenvalue problem on \Omega = [0, 1]^3 with zero boundary
 *     -\nabla^2 u = \lambda u, u_{\partial\Omega} = 0
 * use inveres power method, find the first eigenpair, represented as rayleigh quotient of the calculated state
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
#define PI (4*atan(1.0))

double calc_rq(Vector<double>& x, SparseMatrix<double>& A, SparseMatrix<double>& B){
    Vector<double> tmpv(x.size());
    double nu = 0, de = 0;
    A.vmult(tmpv, x); for (int i = 0; i < x.size(); ++i) nu += tmpv(i) * x(i);
    B.vmult(tmpv, x); for (int i = 0; i < x.size(); ++i) de += tmpv(i) * x(i);
    x /= sqrt(de);
    return nu / de;
}

int main(int argc, char *argv[])
{
    if (argc != 3){
	std::cout << "Usage: " << argv[0] << " polynomial_order n_global_refinement\n";
	return 1;
    }

    
    // read mesh from .mesh file
    HGeometryTree<DIM> h_tree; h_tree.readMesh("./cube.mesh");
    IrregularMesh<DIM> *irregular_mesh = new IrregularMesh<DIM>;
    irregular_mesh->reinit(h_tree);
    if (atoi(argv[2]) > 0) irregular_mesh->globalRefine(atoi(argv[2]));
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

    // construct tsem space
    int M = atoi(argv[1]);
    TSEM<double> tsem(M, mesh, fem_space);
    tsem.build_flag_bm(mesh);
    int n_dof_total = tsem.n_dof_total;

    // solve the eigenvector corresponds to the least eigenvalue, using random initial guess
    SparseMatrix<double> stiff(tsem.stiff_matrix.get_sparsity_pattern());
    stiff.copy_from(tsem.stiff_matrix);
    tsem.impose_zero_boundary_condition(stiff);
    SparseMatrix<double> mass(tsem.mass_matrix.get_sparsity_pattern());
    mass.copy_from(tsem.mass_matrix);
    tsem.impose_zero_boundary_condition(mass);
    Vector<double> solution(n_dof_total);
    for (int i = 0; i < n_dof_total; ++i) solution(i) = rand() * 1. / RAND_MAX;
    tsem.impose_zero_boundary_condition(solution);
    TSEMSolver<double> solver;
    solver.solve_InversePower(solution, stiff, mass, false, 1.0e-12, 1000);

    // calculate eigenvalue error, represented as Rayleigh quotient of solution
    double numerator = 0, denominator = 0;
    Vector<double> tmpv(n_dof_total);
    tsem.stiff_matrix.vmult(tmpv, solution);
    for (int i = 0; i < n_dof_total; ++i) numerator   += tmpv(i) * solution(i);
    tsem.mass_matrix.vmult(tmpv, solution);
    for (int i = 0; i < n_dof_total; ++i) denominator += tmpv(i) * solution(i);
    double eig = numerator / denominator, eig_ref = 3*pow(PI,2);
    double err = fabs(eig - eig_ref);
    std::cout << std::setprecision(16);
    std::cout << "n_element = " << n_element << ", M = " << M << ", n_dof = " << n_dof_total
	      << ", err_eig = " << err << '\n';
    
    
    return 0;
}
