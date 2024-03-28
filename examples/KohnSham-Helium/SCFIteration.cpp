#include "SCFIteration.h"


double wave_init(double *p){ double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]); return exp(-2*r); }

void SCFIteration::normalize()
{
    Vector<double> tmpv(n_dof_total);
    double coef_nor = 0;
    tsem.mass_matrix.vmult(tmpv, psi[0]);
    for (unsigned int idx_dof = 0; idx_dof < n_dof_total; ++idx_dof) coef_nor += tmpv(idx_dof) * psi[0](idx_dof);
    psi[0] /= sqrt(coef_nor);
}

void SCFIteration::assemble_matrix(std::vector<Vector<double> >&psi_now, std::vector<unsigned int>&nnz, SparsityPattern&spp, SparseMatrix<double>&spm)
{
    std::vector<std::vector<double> > val_V(n_element, std::vector<double> (n_qp, 0));
    for (unsigned int idx_ele; idx_ele < n_element; ++idx_ele){
	// recover density from psi_now
	std::vector<double> den_local(n_qp);
	tsem.calc_density(psi_now, den_local, n_occupation, idx_ele);
	// recover V_Har
	std::vector<double> v_Har_local(n_qp);
	tsem.calc_val_qp_onElement(v_Har, v_Har_local, idx_ele);
	// count contribution from each component of potential
	for (unsigned int idx_qp = 0; idx_qp < n_qp; ++idx_qp){
	    val_V[idx_ele][idx_qp] += val_V_Ext[idx_ele][idx_qp] + v_Har_local[idx_qp] + calc_V_XC(den_local[idx_qp]);
	}
    }
    
    // assign sparse matrix M_V
    std::vector<unsigned int> nnz_MV;
    SparsityPattern spp_MV;
    SparseMatrix<double> mass_V;
    tsem.build_mass_V_matrix(nnz_MV, spp_MV, mass_V, val_V);

    // assemble matrix spm = 0.5S + M_V
    for (unsigned int idx_dof = 0; idx_dof < n_dof_total; ++idx_dof) nnz[idx_dof] = nnz_MV[idx_dof] + tsem.n_nonzero_per_row_stiff[idx_dof];
    spp.reinit(n_dof_total, n_dof_total, nnz);
    SparseMatrix<double>::iterator spm_ite = tsem.stiff_matrix.begin(0);
    SparseMatrix<double>::iterator spm_end = tsem.stiff_matrix.end(n_dof_total-1);
    for (; spm_ite != spm_end; ++spm_ite) spp.add(spm_ite->row(), spm_ite->column());
    spm_ite = mass_V.begin(0); spm_end = mass_V.end(n_dof_total-1);
    for (; spm_ite != spm_end; ++spm_ite) spp.add(spm_ite->row(), spm_ite->column());
    spp.compress();
    spm.reinit(spp);
    spm_ite = mass_V.begin(0); spm_end = mass_V.end(n_dof_total-1);
    for (; spm_ite != spm_end; ++spm_ite) spm.add(spm_ite->row(), spm_ite->column(), spm_ite->value());
    spm_ite = tsem.stiff_matrix.begin(0); spm_end = tsem.stiff_matrix.end(n_dof_total-1);
    for (; spm_ite != spm_end; ++spm_ite) spm.add(spm_ite->row(), spm_ite->column(), spm_ite->value() * 0.5);
}

void SCFIteration::update_density()
{
    // update kinetic energy
    double ek = 0;
    Vector<double> tmpv(n_dof_total); tsem.stiff_matrix.vmult(tmpv, psi[0]);
    for (unsigned int idx_dof = 0; idx_dof < n_dof_total; ++idx_dof) ek += tmpv(idx_dof) * psi[0](idx_dof);
    ene_kin += coef_mix * (ek - ene_kin);
    
    // update density
    std::vector<std::vector<double> > density_new(n_element, std::vector<double> (n_qp));
    tsem.calc_density(psi, density_new, n_occupation);
    for (unsigned int idx_ele = 0; idx_ele < n_element; ++idx_ele)
	for (unsigned int idx_qp = 0; idx_qp < n_qp; ++idx_qp)
	    density[idx_ele][idx_qp] += coef_mix * (density_new[idx_ele][idx_qp] - density[idx_ele][idx_qp]);
}

// functions for calculating Hartree potential
void SCFIteration::init_coef_vH_me()
{
    coef_vH_me.resize(10, std::vector<std::vector<double> > (n_element, std::vector<double> (n_qp, 0)));
    for (unsigned int idx_ele = 0; idx_ele < n_element; ++idx_ele){
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(idx_ele).local_to_global(tsem.QPoint);
	std::vector<double> &weight = tsem.Weight[2];
	double &volume = tsem.val_volume[idx_ele];
	// coef_vH_me[0], for <1>
	for (unsigned int p = 0; p < n_qp; ++p) coef_vH_me[0][idx_ele][p] = weight[p] * volume;
	// coef_vH_me[1:3], for <x1>, <x2>, <x3>
	for (unsigned int idx_c = 0; idx_c < 3; ++idx_c)
	    for (unsigned int p = 0; p < n_qp; ++p) coef_vH_me[1+idx_c][idx_ele][p] = weight[p] * q_point[p][idx_c] * volume;
	// coef_vH_me[4:9], for <xixj> with i<j
	int p_coef_me = 4;
	for (unsigned int idx_ci = 0; idx_ci < 3; ++idx_ci)
	    for(unsigned int idx_cj = idx_ci; idx_cj < 3; ++idx_cj){
		for (unsigned int p = 0; p < n_qp; ++p) coef_vH_me[p_coef_me][idx_ele][p] = weight[p] * q_point[p][idx_ci]*q_point[p][idx_cj] * volume;
		p_coef_me++;
	    }
    }
}
void SCFIteration::calc_coef_me()
{
    for (unsigned int idx_cme = 0; idx_cme < 10; ++idx_cme)
	for (unsigned int idx_ele = 0; idx_ele < n_element; ++idx_ele)
	    for (unsigned int idx_qp = 0; idx_qp < n_qp; ++idx_qp)
		coef_me[idx_cme] += coef_vH_me[idx_cme][idx_ele][idx_qp] * density[idx_ele][idx_qp];
}
double SCFIteration::calc_multipole_expansion(AFEPack::Point<DIM> &pos)
{
    double x = pos[0], y = pos[1], z = pos[2];
    double r = sqrt(x*x + y*y + z*z);
    double xr = x/r, yr = y/r, zr = z/r;
    double r_ = 1.0 / r;
    return coef_me[0] * r_
	+ coef_me[1] * xr * pow(r_, 2)
	+ coef_me[2] * yr * pow(r_, 2)
	+ coef_me[3] * zr * pow(r_, 2)
	+ 0.5 * coef_me[4] * (3 * pow(xr, 2) - 1) * pow(r_, 3)
	+ 0.5 * coef_me[7] * (3 * pow(yr, 2) - 1) * pow(r_, 3)
	+ 0.5 * coef_me[9] * (3 * pow(zr, 2) - 1) * pow(r_, 3)
	+ coef_me[5] * 3 * xr * yr * pow(r_, 3)
	+ coef_me[6] * 3 * xr * zr * pow(r_, 3)
	+ coef_me[8] * 3 * yr * zr * pow(r_, 3);
}
void SCFIteration::calc_V_Har_bnd(std::vector<std::vector<std::vector<double> > >&val_bnd)
{
    RegularMesh<DIM>& mesh = irregular_mesh->regularMesh();
    for (unsigned int idx_dim = 0; idx_dim < 3; ++idx_dim){
	val_bnd[idx_dim].resize(tsem.n_geometry[idx_dim]);
	for (unsigned int ind_geo = 0; ind_geo < tsem.n_geometry[idx_dim]; ++ind_geo){
	    if (!tsem.flag_bm[idx_dim][ind_geo]) continue;
	    if (idx_dim == 0) val_bnd[idx_dim][ind_geo].resize(1);
	    else val_bnd[idx_dim][ind_geo].resize(tsem.n_q_point[idx_dim-1]);
	    // 0-d
	    if (idx_dim == 0) val_bnd[idx_dim][ind_geo][0] = calc_multipole_expansion(mesh.point(ind_geo));
	    // 1-d
	    if (idx_dim == 1){
		int ind_point_s = tsem.number_node[0][ind_geo][0], ind_point_e = tsem.number_node[0][ind_geo][1];
		for (unsigned int p = 0; p < tsem.n_q_point[0]; ++p){
		    AFEPack::Point<DIM> p_tmp = mesh.point(ind_point_s), p_ttmp = mesh.point(ind_point_e);
		    p_tmp  *= tsem.QPoint_Barycentric[0][p][1];
		    p_ttmp *= tsem.QPoint_Barycentric[0][p][0];
		    p_tmp  += p_ttmp;
		    val_bnd[idx_dim][ind_geo][p] = calc_multipole_expansion(p_tmp);
		}
	    }
	    // 2-d
	    if (idx_dim == 2){
		for (unsigned int p = 0; p < tsem.n_q_point[1]; ++p){
		    AFEPack::Point<DIM> p_tmp;
		    for (unsigned int ind = 0; ind < DIM; ++ind) p_tmp[ind] = 0;
		    for (unsigned int ind_p = 0; ind_p <= 2; ++ind_p){ // vertex[0]=(0,0)<->[2]=1-x-y; vertex[1]=(1,0)<->[0]=x; vertex[2]=(0,1)<->[1]=y
			int ind_vertex = tsem.number_node[1][ind_geo][ind_p];
			AFEPack::Point<DIM> p_ttmp = mesh.point(ind_vertex);
			p_ttmp *= tsem.QPoint_Barycentric[1][p][(ind_p+2)%3];
			p_tmp  += p_ttmp;
		    }
		    val_bnd[idx_dim][ind_geo][p] = calc_multipole_expansion(p_tmp);
		}
	    }
	}
    }
}
void SCFIteration::calc_V_Har()
{// calculate density
    // calculate coefficients for multipole expansion
    for (unsigned int i = 0; i < 10; ++i) coef_me[i] = 0;
    calc_coef_me();
    
    // calculate right-hand-side
    Vector<double> rhs;
    tsem.calc_rhs(rhs, density);
    rhs *= 4 * PI;

    // calculate Dirichlet boundary
    std::vector<std::vector<std::vector<double> > > val_bnd(3);
    calc_V_Har_bnd(val_bnd);

    // copy stiff matrix, impose boundary condition
    SparseMatrix<double> stiff_Har(tsem.stiff_matrix.get_sparsity_pattern());
    stiff_Har.copy_from(tsem.stiff_matrix);
    tsem.impose_boundary_condition(stiff_Har, rhs, val_bnd, true);

    // solve
    TSEMSolver<double> solver;
    solver.solve_PCG(stiff_Har, v_Har, rhs, tolerance_solver, max_step_solver);
}

// functions for calculating potential and energy
double SCFIteration::calc_V_Ext(AFEPack::Point<DIM> &p, std::vector<int> &charge, std::vector<std::vector<double> > &pos)
{
    double dis, pot = 0;
    for (unsigned int ind_nu = 0; ind_nu < charge.size(); ++ind_nu){
	dis = 0;
	for (unsigned int ind_c = 0; ind_c < 3; ++ind_c)
	    dis += pow(pos[ind_nu][ind_c] - p[ind_c], 2);
	dis = sqrt(dis);
	pot -= charge[ind_nu] * 1.0 / dis;
    }
    return pot;
}
// use LDA from VWN 1980
const double coef_A = 0.0621814, coef_x0 = -0.10498, coef_b = 3.72744, coef_c = 12.9352;
const double coef_Q = sqrt(4*coef_c-pow(coef_b,2)), coef_X0 = pow(coef_x0,2)+coef_b*coef_x0+coef_c;
double SCFIteration::calc_E_C(double rho)
{
    const double A = coef_A, x0 = coef_x0, b = coef_b, c = coef_c, Q = coef_Q, X0 = coef_X0;
    if (rho < Tol_Zero) return rho * A * 0.5*PI*(2*b/Q - b*x0/X0*2*(b+2*x0)/Q) * 0.5;
    double rs = cbrt(0.75/(PI*rho)), x = sqrt(rs);
    double X = pow(x,2) + b*x + c;
    double e_C = A * (log(pow(x,2)/X) + 2*b/Q*atan(Q/(2*x+b))
		      - b*x0/X0*(log(pow(x-x0,2)/X)
				 + 2*(b+2*x0)/Q*atan(Q/(2*x+b))));
    return e_C * rho * 0.5;
}
double SCFIteration::calc_V_C(double rho)
{
    const double A = coef_A, x0 = coef_x0, b = coef_b, c = coef_c, Q = coef_Q, X0 = coef_X0;
    if (rho < Tol_Zero) return A * 0.5*PI*(2*b/Q - b*x0/X0*2*(b+2*x0)/Q) * 0.5;
    double rs = cbrt(0.75/(PI*rho)), x = sqrt(rs);
    double X = pow(x,2) + b*x + c, dX_dx = 2*x + b;
    double e_C = A * (log(pow(x,2)/X) + 2*b/Q*atan(Q/(2*x+b))
		      - b*x0/X0*(log(pow(x-x0,2)/X)
				 + 2*(b+2*x0)/Q*atan(Q/(2*x+b))));
    double dx_drs = 0.5/x;
    double drs_drho = cbrt(0.75/PI) * 1./3.;
    double tmp = X/pow(x,2) * (2*x/X - pow(x/X,2)*dX_dx)
	+ 2*b/Q * 1./(pow(Q/(2*x+b),2) + 1) * -Q/pow(2*x+b,2) * 2
	- b*x0/X0 * (X/pow(x-x0,2) * (2*(x-x0)/X
				      - pow((x-x0)/X,2) * dX_dx)
		     + 2*(b+2*x0)/Q * 1./(pow(Q/(2*x+b),2)+1) * -Q/pow(2*x+b,2) * 2);
    return (e_C + coef_A * tmp * dx_drs * drs_drho * cbrt(rho)) * 0.5; // combine rho in prodcut and rho^-2/3 in drs_drho
}
double SCFIteration::calc_V_X(double rho){ return -cbrt(3./PI * rho); }
double SCFIteration::calc_E_X(double rho){ return -0.75 * cbrt(3./PI) * pow(cbrt(rho), 4); }
double SCFIteration::calc_V_XC(double rho){ return calc_V_X(rho) + calc_V_C(rho); }
double SCFIteration::calc_E_XC(double rho){ return calc_E_X(rho) + calc_E_C(rho); }
double SCFIteration::calc_energy(std::vector<Vector<double> > &psi_now, bool flag_output)
{
    double ene_ext = 0, ene_Har = 0, ene_X = 0, ene_C = 0;
    std::vector<double> val_psi_qp(n_qp), vHar_local(n_qp);
    for (unsigned int idx_ele = 0; idx_ele < n_element; ++idx_ele){
	// recover density on quadrature points
	std::vector<double> den_local(n_qp);
	tsem.calc_density(psi_now, den_local, n_occupation, idx_ele);
	// recover v_Har on quadrature points
	tsem.calc_val_qp_onElement(v_Har, vHar_local, idx_ele);
	double cnt_ext = 0, cnt_h = 0, cnt_x = 0, cnt_c = 0;
	for (unsigned int p = 0; p < n_qp; ++p){
	    cnt_ext += tsem.Weight[2][p] * val_V_Ext[idx_ele][p] * den_local[p];
	    cnt_h   += tsem.Weight[2][p] * vHar_local[p]         * den_local[p];
	    cnt_x   += tsem.Weight[2][p] * calc_E_X(den_local[p]);
	    cnt_c   += tsem.Weight[2][p] * calc_E_C(den_local[p]);
	}
	ene_ext += cnt_ext * tsem.val_volume[idx_ele];
	ene_Har += cnt_h   * tsem.val_volume[idx_ele];
	ene_X   += cnt_x   * tsem.val_volume[idx_ele];
	ene_C   += cnt_c   * tsem.val_volume[idx_ele];
    }
    ene_Har *= 0.5;
    if (flag_output) std::cout << "ene_kin = " << ene_kin << ",\tene_ext = " << ene_ext << ",\tene_Har = " << ene_Har << ",\tene_X = " << ene_X << ",\tene_C = " << ene_C << "\n\t";
    return ene_kin + ene_ext + ene_Har + ene_X + ene_C;
}

// basic functions
SCFIteration::SCFIteration(const std::string& filename_evolve) : info_evolve(filename_evolve)
{
    init();
};
SCFIteration::~SCFIteration()
{
    delete irregular_mesh;
}
void SCFIteration::init()
{
    // mesh
    std::string filename_mesh;
    std::ifstream input(info_evolve);
    getline(input, filename_mesh);
    h_tree.readMesh(filename_mesh);
    irregular_mesh = new IrregularMesh<DIM>;
    irregular_mesh->reinit(h_tree);
    irregular_mesh->semiregularize();
    irregular_mesh->regularize(false);
    RegularMesh<DIM>& mesh = irregular_mesh->regularMesh();
    std::cout << "read mesh with file name: " << filename_mesh << " with element number: " << mesh.n_geometry(3) << '\n';

    // template
    template_geometry.readData("tetrahedron.tmp_geo");
    coord_transform.readData("tetrahedron.crd_trs");
    template_dof.reinit(template_geometry);     template_dof.readData("tetrahedron.2.tmp_dof");
    basis_function.reinit(template_dof);        basis_function.readData("tetrahedron.2.bas_fun");
    template_element.resize(1);
    template_element[0].reinit(template_geometry, template_dof, coord_transform, basis_function);

    // fem space
    fem_space.reinit(mesh, template_element);
    n_element = mesh.n_geometry(DIM);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i) fem_space.element(i).reinit(fem_space, i, 0);
    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();

    // tsem space
    input >> order_tsem;
    tsem.init(order_tsem, mesh, fem_space);
    tsem.build_flag_bm(mesh);
    n_dof_total = tsem.n_dof_total;
    n_qp = tsem.n_q_point[2];
    std::cout << "order_tsem = " << order_tsem << ", n_dof_total = " << n_dof_total << '\n';

    // evolution info
    input >> coef_mix
	  >> tolerance_scf >> max_step_scf
	  >> tolerance_lobpcg >> max_step_lobpcg
	  >> tolerance_solver >> max_step_solver;
    std::cout << "coef_mix = " << coef_mix
	      << ", tolerance_SCF = " << tolerance_scf << ", max_step_SCF = " << max_step_scf
	      << ", tolerance_LOBPCG = " << tolerance_lobpcg << ", max_step_LOBPCG = " << max_step_lobpcg
	      << ", tolerance_solver = " << tolerance_solver << ", max_step_solver = " << max_step_solver
	      << '\n';

    // output nucleus info
    std::cout << "n_orbital = " << n_orbital << ", n_occupation = [" << n_occupation[0] << "]\n";

    // io info
    input >> flag_read >> flag_write >> flag_write_geometry;
    std::string tmps; getline(input, tmps);
    if (flag_read){
	getline(input, filename_read);
	std::cout << "read psi from file: " << filename_read << '\n';
    }
    if (flag_write){
	getline(input, filename_write);
	std::cout << "write psi to file: " << filename_write << std::endl;
    }
    if (flag_write_geometry) std::cout << "write geometry info to mesh.dx and density to density.dx\n";
    input >> flag_output;
    getline(input, tmps);
    if (flag_output){
	getline(input, filename_output_energy);
	getline(input, filename_output_orthonormal);
	std::cout << "output energy to file: " << filename_output_energy
		  << ", output orthonormal properties to file: " << filename_output_orthonormal
		  << std::endl;
    }
    
    input.close();

    
    // normalized initial value
    psi.resize(n_orbital, Vector<double> (n_dof_total));
    SparseMatrix<double> &mass = tsem.mass_matrix;
    Vector<double> tmpv(n_dof_total);
    if (!flag_read){
	tsem.calc_interpolation(psi[0], mesh, fem_space, &wave_init);
	normalize();
    }
    else
	tsem.read_coef(psi, filename_read);
    // calculate density
    density.resize(n_element, std::vector<double> (n_qp));
    tsem.calc_density(psi, density, n_occupation);
    // calculate kinetic energy
    tsem.stiff_matrix.vmult(tmpv, psi[0]);
    for (unsigned int idx_dof = 0; idx_dof < n_dof_total; ++idx_dof) ene_kin += tmpv(idx_dof) * psi[0](idx_dof);


    // calculate external potential on quadrature point
    val_V_Ext.resize(n_element, std::vector<double> (n_qp));
    for (unsigned int idx_ele = 0; idx_ele < n_element; ++idx_ele){	
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(idx_ele).local_to_global(tsem.QPoint);
	for (unsigned int p = 0; p < n_qp; ++p) val_V_Ext[idx_ele][p] = calc_V_Ext(q_point[p], charge_nu, pos_nu);
    }
    // initialize v_Har, calculate coef_vH_me
    v_Har.reinit(n_dof_total);
    init_coef_vH_me();
}
void SCFIteration::run()
{
    solve();
    outputSolution();
}
void SCFIteration::solve()
{
    // calculate initial energy
    calc_V_Har();
    double ene_last = calc_energy(psi, false);
    std::cout << "initial energy: " << ene_last << '\n';

    TSEMSolver<double> solver;
    for (unsigned int n_ite = 1; n_ite <= max_step_scf; ++n_ite){
	// assign matrix
	std::vector<unsigned int> n_nonzero(n_dof_total, 0);
	SparsityPattern sp_pattern;
	SparseMatrix<double> sp_matrix;
	assemble_matrix(psi, n_nonzero, sp_pattern, sp_matrix);

	// solve linear eigenvalue problem
	solver.solve_LOBPCG(tsem, psi, sp_matrix, tsem.mass_matrix, false, tolerance_solver, max_step_solver);
	normalize();

	// calculate v_Har, energy and error
	update_density();
	// tsem.calc_density(psi, density, n_occupation);
	calc_V_Har();
	double ene = calc_energy(psi, true);
	double err = fabs(ene-ene_last);
	std::cout.precision(8);
	std::cout << "n_ite = " << n_ite << ", ene = " << ene << ", err_ene = " << err << '\n';

	// break critirion
	if (err < tolerance_scf) break;
	ene_last = ene;
    }
}
void SCFIteration::outputSolution()
{
    // TSEM solution 
    if (flag_write) tsem.write_coef(psi, filename_write);
    if (flag_write_geometry){
	// mesh
	FEMFunction<double, 3> fem_function(fem_space);
	fem_function.writeOpenDXData("result/mesh.dx");
	// RegularMesh<3>& mesh = irregular_mesh->regularMesh();
	// mesh.writeOpenDXData("mesh.dx");
	// visualization of density
	Visualize visualize;
	visualize.read(&tsem, &h_tree, &psi, &n_occupation);
	visualize.build();
	visualize.write("result/density.dx");
    }
}


Visualize::Visualize()
{
    n_globalRefine() = 2;
    n_localRefine() = 3;
    rateRefine() = 0.1;
    n_smoothIndicator() = 1;
    Tolerance() = 0;
}
Visualize::~Visualize(){}

void Visualize::get_indicator()
{ // use l2 norm of gradient as indicator
    RegularMesh<3>& mesh = irregular_mesh->regularMesh();
    val_indicator.resize(mesh.n_geometry(3), 0);
    std::vector<std::vector<int> > idx_vtx = {{0,1,2,3}, {0,1,3,4}};
    for (unsigned int idx_ele = 0; idx_ele < mesh.n_geometry(3); ++idx_ele){
	int idx = mesh.geometry(3,idx_ele).n_vertex() == 5 ? 1: 0;
	std::vector<int> idx_vertex = {mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][0]),
				       mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][1]),
				       mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][2]),
				       mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][3])};
	// std::cout << idx_ele << ':' << mesh.geometry(3,idx_ele).n_vertex() << ',' << idx << '-' << idx_vertex[0] << ',' << idx_vertex[1] << ',' << idx_vertex[2] << ',' << idx_vertex[3] << '\t';
	AFEPack::Point<3>& x0 = mesh.point(idx_vertex[0]);
	AFEPack::Point<3>& x1 = mesh.point(idx_vertex[1]);
	AFEPack::Point<3>& x2 = mesh.point(idx_vertex[2]);
	AFEPack::Point<3>& x3 = mesh.point(idx_vertex[3]);
	double det = fabs((x1[0] - x0[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
			  (x3[0] - x0[0])*(x1[1] - x0[1])*(x2[2] - x0[2]) +
			  (x2[0] - x0[0])*(x3[1] - x0[1])*(x1[2] - x0[2]) -
			  (x1[0] - x0[0])*(x3[1] - x0[1])*(x2[2] - x0[2]) -
			  (x2[0] - x0[0])*(x1[1] - x0[1])*(x3[2] - x0[2]) -
			  (x3[0] - x0[0])*(x2[1] - x0[1])*(x1[2] - x0[2]));
	std::vector<double> fi = {val_interp[idx_vertex[0]], val_interp[idx_vertex[1]], val_interp[idx_vertex[2]], val_interp[idx_vertex[3]]};
	double detx = fabs((fi[1] - fi[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
			   (fi[3] - fi[0])*(x1[1] - x0[1])*(x2[2] - x0[2]) +
			   (fi[2] - fi[0])*(x3[1] - x0[1])*(x1[2] - x0[2]) -
			   (fi[1] - fi[0])*(x3[1] - x0[1])*(x2[2] - x0[2]) -
			   (fi[2] - fi[0])*(x1[1] - x0[1])*(x3[2] - x0[2]) -
			   (fi[3] - fi[0])*(x2[1] - x0[1])*(x1[2] - x0[2]));
	double dety = fabs((x1[0] - x0[0])*(fi[2] - fi[0])*(x3[2] - x0[2]) +
			   (x3[0] - x0[0])*(fi[1] - fi[0])*(x2[2] - x0[2]) +
			   (x2[0] - x0[0])*(fi[3] - fi[0])*(x1[2] - x0[2]) -
			   (x1[0] - x0[0])*(fi[3] - fi[0])*(x2[2] - x0[2]) -
			   (x2[0] - x0[0])*(fi[1] - fi[0])*(x3[2] - x0[2]) -
			   (x3[0] - x0[0])*(fi[2] - fi[0])*(x1[2] - x0[2]));
	double detz = fabs((x1[0] - x0[0])*(x2[1] - x0[1])*(fi[3] - fi[0]) +
			   (x3[0] - x0[0])*(x1[1] - x0[1])*(fi[2] - fi[0]) +
			   (x2[0] - x0[0])*(x3[1] - x0[1])*(fi[1] - fi[0]) -
			   (x1[0] - x0[0])*(x3[1] - x0[1])*(fi[2] - fi[0]) -
			   (x2[0] - x0[0])*(x1[1] - x0[1])*(fi[3] - fi[0]) -
			   (x3[0] - x0[0])*(x2[1] - x0[1])*(fi[1] - fi[0]));
	val_indicator[idx_ele] = sqrt((pow(detx,2)+pow(dety,2)+pow(detz,2)) / pow(det,2) * det/6.); // L2 integral of gradient
    }
}
