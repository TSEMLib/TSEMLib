#include "Poisson.h"

const double Para = 20;
double u(double *p){ return atan(Para*(p[0]+p[1])); }
double f(double *p){ return 4*pow(Para,3)*(p[0]+p[1])/pow(1+pow(Para*(p[0]+p[1]),2),2); }

Poisson::Poisson(const std::string& filename_mesh, int polynomial_order, int move_step_number) : path_mesh(filename_mesh), order_tsem(polynomial_order), n_move(move_step_number)
{
    init();
}
Poisson::~Poisson(){};
void Poisson::init()
{
    // init for moving mesh
    readDomain(path_mesh);
    
    // read mesh
    h_tree.readMesh(path_mesh);
    irregular_mesh = new IrregularMesh<DIM>;
    irregular_mesh->reinit(h_tree);
    irregular_mesh->semiregularize();
    irregular_mesh->regularize(false);
    RegularMesh<DIM>& mesh = irregular_mesh->regularMesh();
    std::cout << "read mesh with file name: " << path_mesh << " with element number: " << mesh.n_geometry(3) << '\n';

    // template
    template_geometry.readData("tetrahedron.tmp_geo");
    coord_transform.readData("tetrahedron.crd_trs");
    template_dof.reinit(template_geometry);     template_dof.readData("tetrahedron.2.tmp_dof");
    basis_function.reinit(template_dof);        basis_function.readData("tetrahedron.2.bas_fun");
    template_element.resize(1);
    template_element[0].reinit(template_geometry, template_dof, coord_transform, basis_function);

    // fem space
    fem_space.reinit(*this, template_element);
    n_element = this->n_geometry(DIM);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i) fem_space.element(i).reinit(fem_space, i, 0);
    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();

    // tsem space
    tsem.init(order_tsem, *this, fem_space);
    tsem.build_flag_bm(mesh);
    n_dof_total = tsem.n_dof_total;
    std::cout << "order_tsem = " << order_tsem << ", n_dof_total = " << n_dof_total << '\n';
    solution.reinit(n_dof_total);
    rhs.reinit(n_dof_total);

    tolerence() = 10000; // such that move single step per iteration
}

void Poisson::run()
{
    solve();
    n_ite = 0; outputSolution();
 
    for (n_ite = 1; n_ite <= n_move; ++n_ite){
	moveMesh();
	outputSolution();
    }
}

void Poisson::solve()
{
    tsem.build_geometry_info(*this, fem_space);
    tsem.build_stiff_matrix(fem_space);
    sp_matrix.reinit(tsem.stiff_matrix.get_sparsity_pattern());
    sp_matrix.copy_from(tsem.stiff_matrix);
    tsem.calc_rhs(rhs, fem_space, &f);
    tsem.impose_boundary_condition(sp_matrix, rhs, *this, &u, true, -1);
    TSEMSolver<double> solver;
    solver.solve_PCG(sp_matrix, solution, rhs, 1.0e-12, 100000);
}


void Poisson::outputSolution()
{
    std::ostringstream oss;
    oss << "./result/mesh" << n_ite << ".mesh";
    outputPhysicalMesh(oss.str());
    std::ostringstream oss_mesh;
    oss_mesh << "./result/FEMFunction" << n_ite << ".dx";
    FEMFunction<double, 3> fem_function(fem_space);
    fem_function.writeOpenDXData(oss_mesh.str());
    output_err(n_ite, oss.str());
}

void Poisson::output_err(int n_ite, std::string str)
{
    // read mesh from .mesh file
    HGeometryTree<DIM> h_tree_t; h_tree_t.readMesh(str);
    IrregularMesh<DIM> *irregular_mesh_t = new IrregularMesh<DIM>;
    irregular_mesh_t->reinit(h_tree_t);
    irregular_mesh_t->semiregularize();
    irregular_mesh_t->regularize(false);
    RegularMesh<DIM> &mesh = irregular_mesh_t->regularMesh();
    
    // construct finite element space
    FEMSpace<double, DIM> fem_space_t(mesh, template_element);
    fem_space_t.element().resize(n_element);
    for (int i = 0; i < n_element; ++i) fem_space_t.element(i).reinit(fem_space_t, i, 0);
    fem_space_t.buildElement();
    fem_space_t.buildDof();
    fem_space_t.buildDofBoundaryMark();

    // construct tsem space, with quadrature info order 30 in each direction
    TSEM<double> tsem_t(order_tsem, mesh, fem_space_t);
    tsem_t.build_flag_bm(mesh);
    int n_dof_total_t = tsem_t.n_dof_total;

    // build linear system
    SparseMatrix<double> sp_matrix_t(tsem_t.stiff_matrix.get_sparsity_pattern());
    sp_matrix_t.copy_from(tsem_t.stiff_matrix);
    Vector<double> rhs_t(n_dof_total_t);
    tsem_t.calc_rhs(rhs_t, fem_space_t, &f);
    tsem_t.impose_boundary_condition(sp_matrix_t, rhs_t, mesh, &u, true, -1);

    // solve
    TSEMSolver<double> solver;
    Vector<double> solution_t(n_dof_total_t);
    solver.solve_PCG(sp_matrix_t, solution_t, rhs_t, 1.0e-12, 100000);

    // output error
    double err_l2  = tsem_t.calc_l2_error(solution_t, fem_space_t, &u);
    std::cout.precision(16);
    std::cout << "n_ite = " << n_ite
    	      << ", n_element = " << n_element
    	      << ", M = " << order_tsem << ", n_dof = " << n_dof_total_t << '\n'
    	      << "\terr_l2  = " << err_l2 
    	      << '\n';
}

void Poisson::getMonitor()
{
    solve();
    
    Vector<double> res(n_dof_total);
    sp_matrix.vmult(res, solution);
    res.add(-1, rhs);
    std::vector<std::vector<double> > tmp(tsem.n_q_point[2], std::vector<double> (tsem.n_index, 0));
    for (unsigned int idx_ele = 0; idx_ele < n_element; ++idx_ele)
	monitor(idx_ele) = tsem.calc_l2_error_gradient(*this, solution, tmp, idx_ele);
    double C = pow(0.5, n_ite);
    for (int i = 0;i < n_geometry(3);i ++)
	monitor(i) = C/(C + monitor(i));
}

void Poisson::updateSolution(){}
