#include "../include/TetrahedralSEM.h"

// This file contains IO functions and ...
// read_coef()
// write_coef()
//
// build_flag_mask()
// calc_area_triangle()
// calc_binomial_coefficient()
// calc_coefficient_a()
// calc_coefficient_b()
// calc_coefficient_c()
// calc_coefficient_d()
// calc_coefficient_e()
// calc_coefficient_g()
// calc_coefficient_rho()
// calc_coefficient_kappa()
// calc_coefficient_theta()
// calc_coefficient_D()
// calc_cross_product()
// calc_delta()
// calc_dipole()
// calc_generalized_jacobi_polynomial
// calc_inner_product()
// calc_length_line()
// calc_sum_dof()
// calc_volume_tetrahedron()
// impose_mask()


#define TEMPLATE_TSEM template<typename valuetype>
#define THIS_TSEM TSEM<valuetype>

TEMPLATE_TSEM
void THIS_TSEM::read_coef(Vector<valuetype> &dst, std::string filename)
{
    dst.reinit(n_dof_total);
    std::ifstream input_coef(filename);
    int M_interp, n_dof_total_interp;
    input_coef >> M_interp;
    input_coef >> n_dof_total_interp;
    // construct pointer
    std::vector<int> n_dof_geometry_interp(4); // number of dof location on each dimensional geoemtry for interpolation
    for (int ind = 0; ind <= 3; ++ind) // number of dof on ind dimensional geometry
	n_dof_geometry_interp[ind] = calc_binomial_coefficient(M_interp-1, ind);
    std::vector<int> location_actualdof_interp(n_geometry_total); // start index of geometry in actual discretized matrix for interpolation;
    location_actualdof_interp[0] = 0;
    for (int i = 1; i < n_geometry_total; ++i)
	location_actualdof_interp[i] = location_actualdof_interp[i-1] + n_dof_geometry_interp[geometry_dimension[i-1]];
    // read dof info
    Vector<valuetype> src(n_dof_total_interp);
    for (int i = 0; i < n_dof_total_interp; ++i)
	input_coef >> src(i);
    // initialize
    for (int i = 0; i < n_dof_total; ++i)
	dst(i) = 0;
    // interpolate
    for (int ind = 0; ind <= 3; ++ind)
	for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
	    int location = ind_geo + ((ind >= 1) ? n_geometry[0] : 0) + ((ind >= 2) ? n_geometry[1] : 0) + ((ind >= 3) ? n_geometry[2] : 0);
	    int position =        location_actualdof[location_geometry[location]];
	    int position_interp = location_actualdof_interp[location_geometry[location]];
	    for (int ind_index = 0; ind_index < n_dof_geometry_interp[ind]; ++ind_index)
		dst(position + ind_index) = src(position_interp + ind_index);
	}
    input_coef.close();
}

TEMPLATE_TSEM
void THIS_TSEM::read_coef(std::vector<Vector<valuetype> > &dst, std::string filename)
{ // suppose the length of dst and its element match the info in filename
    std::ifstream input_coef(filename);
    int M_interp, n_dof_total_interp;
    input_coef >> M_interp;
    input_coef >> n_dof_total_interp;
    // construct pointer
    std::vector<int> n_dof_geometry_interp(4); // number of dof location on each dimensional geoemtry for interpolation
    for (int ind = 0; ind <= 3; ++ind) // number of dof on ind dimensional geometry
	n_dof_geometry_interp[ind] = calc_binomial_coefficient(M_interp-1, ind);
    std::vector<int> location_actualdof_interp(n_geometry_total); // start index of geometry in actual discretized matrix for interpolation;
    location_actualdof_interp[0] = 0;
    for (int i = 1; i < n_geometry_total; ++i)
	location_actualdof_interp[i] = location_actualdof_interp[i-1] + n_dof_geometry_interp[geometry_dimension[i-1]];
    
    for (int ind_vec = 0; ind_vec < dst.size(); ++ind_vec){
	// read dof info
	Vector<valuetype> src(n_dof_total_interp);
	for (int i = 0; i < n_dof_total_interp; ++i)
	    input_coef >> src(i);
	// initialize
	for (int i = 0; i < n_dof_total; ++i)
	    dst[ind_vec](i) = 0;
	// interpolate
	for (int ind = 0; ind <= 3; ++ind)
	    for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
		int location = ind_geo + ((ind >= 1) ? n_geometry[0] : 0) + ((ind >= 2) ? n_geometry[1] : 0) + ((ind >= 3) ? n_geometry[2] : 0);
		int position =        location_actualdof[location_geometry[location]];
		int position_interp = location_actualdof_interp[location_geometry[location]];
		for (int ind_index = 0; ind_index < n_dof_geometry_interp[ind]; ++ind_index)
		    dst[ind_vec](position + ind_index) = src(position_interp + ind_index);
	    }
    }
    input_coef.close();
}

TEMPLATE_TSEM
void THIS_TSEM::write_coef(Vector<valuetype> &src, std::string filename)
{
    std::ofstream output_coef(filename);
    output_coef << M << '\t' << n_dof_total << '\n';
    for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
	output_coef << std::setprecision(14) << src(ind_dof) << '\n';
    output_coef.close();
}

TEMPLATE_TSEM
void THIS_TSEM::write_coef(std::vector<Vector<valuetype> > &src, std::string filename)
{
    std::ofstream output_coef(filename);
    output_coef << M << '\t' << n_dof_total << '\n';
    for (int ind_vec = 0; ind_vec < src.size(); ++ind_vec)
	for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
	    output_coef << std::setprecision(14) << src[ind_vec](ind_dof) << '\n';
    output_coef.close();
}


TEMPLATE_TSEM
void THIS_TSEM::build_flag_mask(RegularMesh<3> &mesh, valuetype Rinner, valuetype Router)
{
    Rmin = Rinner; Rmax = Router;
    flag_mask.resize(4);
    // point
    flag_mask[0].resize(n_geometry[0]);
    for (unsigned int idx_p = 0; idx_p < n_geometry[0]; ++idx_p){
	AFEPack::Point<3> &p = mesh.point(idx_p);
	valuetype dis = sqrt(pow(p[0],2)+pow(p[1],2)+pow(p[2],2));
	if (dis < Rmin){ flag_mask[0][idx_p] = 0; continue; }
	if (dis > Rmax){ flag_mask[0][idx_p] = 2; continue; }
	flag_mask[0][idx_p] = 1;
    }
    // others
    unsigned int tf;
    for (unsigned int idx = 1; idx <= 3; ++idx){
	flag_mask[idx].resize(n_geometry[idx], 0);
	for (unsigned int idx_geo = 0; idx_geo < n_geometry[idx]; ++idx_geo){
	    tf = flag_mask[idx-1][mesh.geometry(idx, idx_geo).boundary(0)];
	    for (unsigned int idx_bnd = 1; idx_bnd <= idx+1; ++idx_bnd)
		if (tf != flag_mask[idx-1][mesh.geometry(idx, idx_geo).boundary(idx_bnd)]){
		    tf = 1; break; }
	    flag_mask[idx][idx_geo] = tf;
	}
    }
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_area_triangle(AFEPack::Point<3> &p0, AFEPack::Point<3> &p1, AFEPack::Point<3> &p2)
{
    return 0.5 * sqrt(pow((p1[1]-p0[1])*(p2[2]-p0[2]) - (p1[2]-p0[2])*(p2[1]-p0[1]), 2) +
		      pow((p1[0]-p0[0])*(p2[2]-p0[2]) - (p1[2]-p0[2])*(p2[0]-p0[0]), 2) +
		      pow((p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]), 2));
}

TEMPLATE_TSEM
int THIS_TSEM::calc_binomial_coefficient(int n, int m)
{// calculate binomial coefficient $\binom{n}{m}$
    int val = 1;
    if (n > 0 && n >= m){
        if (n-m < m) m = n-m;
        for (int i = 0; i < m; ++i){val *= (n-i); val /= (i+1);}}
    return val;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_a(int alpha, int beta, int ind, int k)
{// calculate coefficient $a_{ind, k}^{alpha, beta}$
    valuetype ans;
    if (k < 0){
	ans = (valuetype) 0;
	return ans;
    }
    switch (ind){
    case 1:
	ans = ((valuetype) (2 * (k+1) * (k+alpha+beta+1))) / ((2*k+alpha+beta+1) * (2*k+alpha+beta+2));
        if (alpha == -1 && beta == -1){
            switch (k){
            case 0: ans = (valuetype) 1;   break;
            case 1: ans = (valuetype) 4;   break;
            case 2: ans = (valuetype) 0.5; break;
            }
	    break;
	}
        if (k == 0 && alpha + beta != -2){
            ans = ((valuetype) 2) / (alpha + beta + 2);
	    break;
	}
	break;
    case 2:
        if (0 <= k && k <= 2 && alpha == -1 && beta == -1){
	    ans = (valuetype) 0;
	    break;
	}
        if (k == 0 && alpha + beta != -2){
            ans = ((valuetype) (beta - alpha)) / (alpha + beta + 2);
	    break;
	}
	ans = ((valuetype) (pow(beta,2)-pow(alpha,2))) / ((2*k+alpha+beta) * (2*k+alpha+beta+2));
	break;
    case 3:
	ans = ((valuetype) (2 * (k+alpha) * (k+beta))) / ((2*k+alpha+beta) * (2*k+alpha+beta+1));
        if (alpha == -1 && beta == -1){
            if (k == 0) ans = (valuetype) 0;
            if (k == 1) ans = (valuetype) 1;
            if (k == 2) ans = (valuetype) 0;
	    break;
        }
        if (k == 0 && alpha+beta != -2){
            ans = (valuetype) 0;
	    break;
	}
	break;
    }
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_b(int alpha, int beta, int ind, int k)
{// calculate coefficient $b_{ind, k}^{alpha, beta}$
    valuetype ans;
    if (k < 0){
	ans = (valuetype) 0;
	return ans;
    }
    switch (ind){
    case 1:
        if (k == 0 && alpha >= -1 && beta >= -1){
	    ans = (valuetype) 1;
	    break;
	}
        if (k == 1 && alpha == -1 && beta == -1){
	    ans = (valuetype) 2;
	    break;
	}
	ans = ((valuetype) (k+alpha+beta+1)) / (2*k+alpha+beta+1);
	break;
    case 2:
        if (k == 0 && alpha >= -1 && beta >= -1){
	    ans = (valuetype) 0;
	    break;
	}
        if (k == 1 && alpha == -1 && beta == -1){
	    ans = (valuetype) -1;
	    break;
	}
        ans = ((valuetype) -(k+beta)) / (2*k+alpha+beta+1);
        break;
    }
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_c(int alpha, int beta, int ind, int k)
{// calculate coefficient $c_{ind, k}^{alpha, beta}$
    valuetype ans;
    switch (ind){
    case 1:
	ans = calc_coefficient_b(alpha, beta, 1, k) * calc_coefficient_b(alpha+1, beta, 1, k);
	break;
    case 2:
	ans = calc_coefficient_b(alpha, beta, 1, k) * calc_coefficient_b(alpha+1, beta, 2, k) +
	      calc_coefficient_b(alpha, beta, 2, k) * calc_coefficient_b(alpha+1, beta, 1, k-1);
	break;
    case 3:
	ans = calc_coefficient_b(alpha, beta, 2, k) * calc_coefficient_b(alpha+1, beta, 2, k-1);
        break;
    }
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_d(int alpha, int beta, int k)
{// calculate coefficient $d_k^{alpha, beta}$
    valuetype ans;
    if (k <= 0){
	ans = (valuetype) 0;
	return ans;
    }
    if (k == 1 && alpha == -1 && beta == -1){
	ans = (valuetype) 1;
	return ans;
    }
    ans = ((valuetype) (k+alpha+beta+1)) / 2;
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_e(int alpha, int beta, int ind, int k)
{// calculate coefficient $e_{ind, k}^{alpha, beta}$
    valuetype ans;
    if (k < 0){
	ans = (valuetype) 0;
	return ans;
    }
    switch (ind){
    case 1:
        if (k == 0 && alpha == -1 && beta == -1){
	    ans = (valuetype) 0.5;
	    break;
	}
        if (k == 1 && alpha == -1 && beta == -1){
	    ans = (valuetype) 0;
	    break;
	}
        ans = ((valuetype) (k+alpha+1)) / (2*k+alpha+beta+2);
	break;
    case 2:
        if (k == 0 && alpha == -1 && beta == -1){
	    ans = (valuetype) -0.5;
	    break;
	}
        if (k == 1 && alpha == -1 && beta == -1){
	    ans = (valuetype) -1;
	    break;
	}
        ans = ((valuetype) -(k+1)) / (2*k+alpha+beta+2);
        break;
    }
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_g(int alpha, int beta, int ind, int k)
{// calculate coefficient $g_{ind, k}^{alpha, beta}$
    valuetype ans;
    switch (ind){
    case 1:
	ans = calc_coefficient_e(alpha+1, beta, 2, k) * calc_coefficient_e(alpha, beta, 2, k+1);
	break;
    case 2:
	ans = calc_coefficient_e(alpha+1, beta, 1, k) * calc_coefficient_e(alpha, beta, 2, k) +
	      calc_coefficient_e(alpha+1, beta, 2, k) * calc_coefficient_e(alpha, beta, 1, k+1);
	break;
    case 3:
	ans = calc_coefficient_e(alpha+1, beta, 1, k) * calc_coefficient_e(alpha, beta, 1, k);
        break;
    }
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_rho(Multiindex<3> index)
{// calculate coefficient $\rho_{index}^{-\mathds{1}}$
    valuetype ans;
    int l1 = index.index[0], l2 = index.index[1];
    ans = 2 * calc_coefficient_d(-1, -1, l1) * calc_coefficient_e(-1, 0, 1, l1-1)
        - l1 * calc_coefficient_b(-1, -1, 2, l1);
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_kappa(Multiindex<3> index)
{// calculate coefficient $\kappa_{index}^{-\mathds{1}}$
    valuetype ans;
    int l1 = index.index[0], l2 = index.index[1];
    ans = l1 * calc_coefficient_b(-1, -1, 2, l1)
        - 2 * calc_coefficient_d(-1, -1, l1) * calc_coefficient_e(-1, 0, 1, l1-1);
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_theta(Multiindex<3> index)
{// calculate coefficient $\theta_{index}^{-\mathds{1}}$
    valuetype ans;
    int l1 = index.index[0], l2 = index.index[1];
    ans = 2 * calc_coefficient_d(2*l1-1, -1, l2) * calc_coefficient_e(-1, 2*l1, 1, l2-1)
        - l2 * calc_coefficient_b(2*l1-1, -1, 2, l2);
    return ans;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_coefficient_D(int ind_derivative, int ind_variation, Multiindex<3> index)
{// calculate coefficient for the derivative of generalized jacobi polynomial
    valuetype ans;
    int l1 = index.index[0], l2 = index.index[1], l3 = index.index[2];
    int aph2 = 2 * l1, aph3 = 2 * l1 + 2 * l2;
    switch (ind_derivative){
    case 0: // $\partial x1$,               or multiindex (2, 0, 0)
        ans = 2 * calc_coefficient_d(-1, -1, l1);
	break;
    case 1: // $\partial x2 - \partial x1$, or multiindex (1, 1, 0)
        switch (ind_variation){
        case 0: // $D^{21}_0 =  D^{2}_0$
            ans = 2 * calc_coefficient_d(aph2-1, -1, l2) * calc_coefficient_b(-1, -1, 1, l1);
	    break;
        case 1: // $D^{21}_1 = -D^{2}_1$
            ans =  -((2*calc_coefficient_d(-1,-1,l1)*calc_coefficient_e(-1,0,1,l1-1) - l1*calc_coefficient_b(-1,-1,2,l1)) * calc_coefficient_b(aph2-1,-1,1,l2)
		     + 2*calc_coefficient_d(aph2-1,-1,l2)*calc_coefficient_b(-1,-1,2,l1)*calc_coefficient_e(aph2-1,0,2,l2-1)) / calc_coefficient_b(aph2-2,0,1,l2);
	    break;
        }
	break;
    case 2: // $\partial x1 - \partial x3$, or multiindex (1, 0, 1)
        switch (ind_variation){
        case 0: // $D^{31}_{0,0} = -D^3_{0,0}$
            ans = -2 * calc_coefficient_d(aph3-1,-1,l3) * calc_coefficient_b(-1,-1,1,l1) * calc_coefficient_b(aph2-1,-1,1,l2);
	    break;
        case 1: // $D^{31}_{1,0} =  D^3_{1,0}$
            ans = 2 * calc_coefficient_d(aph3-1,-1,l3) * calc_coefficient_b(-1,-1,2,l1) * calc_coefficient_e(aph2-2,-1,2,l2);
	    break;
        case 2: // $D^{31}_{0,1} = -D^3_{0,1}$
            ans = -(calc_coefficient_b(-1,-1,1,l1)*calc_coefficient_theta(index)*calc_coefficient_b(aph3-1,-1,1,l3)
		    + 2*calc_coefficient_d(aph3-1,-1,l3)*calc_coefficient_b(-1,-1,1,l1)*calc_coefficient_b(aph2-1,-1,2,l2)
		    *calc_coefficient_e(aph3-1,0,2,l3-1)) / calc_coefficient_b(aph3-2,0,1,l3);
	    break;
        case 3:
            ans = ((calc_coefficient_b(-1,-1,2,l1)*calc_coefficient_e(aph2-1,-1,2,l2-1)*calc_coefficient_theta(index) - calc_coefficient_kappa(index))
		   * calc_coefficient_b(aph3-1,-1,1,l3) + 2*calc_coefficient_b(aph2-2,-1,1,l2)*calc_coefficient_d(aph3-1,-1,l3)
		   *calc_coefficient_b(-1,-1,2,l1)*calc_coefficient_e(aph2-2,-1,1,l2)*calc_coefficient_e(aph3-1,0,2,l3-1))
                / (calc_coefficient_b(aph2-2,-1,1,l2)*calc_coefficient_b(aph3-2,0,1,l3));
	    break;
        }
	break;
    case 3: // $\partial x2$,               or multiindex (0, 2, 0)
        switch (ind_variation){
        case 0:
	    ans = 2 * calc_coefficient_d(aph2-1, -1, l2) * calc_coefficient_b(-1, -1, 1, l1);
	    break;
        case 1:
	    ans = ((2*calc_coefficient_d(-1,-1,l1)*calc_coefficient_e(-1,0,1,l1-1) - l1*calc_coefficient_b(-1,-1,2,l1)) * calc_coefficient_b(aph2-1,-1,1,l2)
		   + 2*calc_coefficient_d(aph2-1,-1,l2)*calc_coefficient_b(-1,-1,2,l1)*calc_coefficient_e(aph2-1,0,2,l2-1)) / calc_coefficient_b(aph2-2,0,1,l2);
	    break;
        }
	break;
    case 4: // $\partial x3 - \partial x2$, or multiindex (0, 1, 1)
        switch (ind_variation){
        case 0:
	    ans = 2 * calc_coefficient_d(aph3-1,-1,l3) * calc_coefficient_b(aph2-1,-1,1,l2);
	    break;
        case 1:
            ans = ((l2*calc_coefficient_b(-1,aph2-1,2,l2) - 2*calc_coefficient_d(aph2-1,-1,l2)*calc_coefficient_e(aph2-1,0,1,l2-1))
		   * calc_coefficient_b(aph3-1,-1,1,l3) - 2*calc_coefficient_d(aph3-1,-1,l3)*calc_coefficient_b(-1,aph2-1,2,l2)
		   *calc_coefficient_e(aph3-1,0,2,l3-1)) / calc_coefficient_b(aph3-2,0,1,l3);
	    break;
        }
	break;
    case 5: // $\partial x3$,               or multiindex (0, 0, 2)
        switch (ind_variation){
        case 0:
	    ans = 2 * calc_coefficient_d(aph3-1,-1,l3) * calc_coefficient_b(-1,-1,1,l1) * calc_coefficient_b(aph2-1,-1,1,l2);
	    break;
        case 1:
	    ans = 2 * calc_coefficient_d(aph3-1,-1,l3) * calc_coefficient_b(-1,-1,2,l1) * calc_coefficient_e(aph2-2,-1,2,l2);
	    break;
        case 2:
            ans = (calc_coefficient_b(-1,-1,1,l1)*calc_coefficient_theta(index)*calc_coefficient_b(aph3-1,-1,1,l3)
		   + 2*calc_coefficient_d(aph3-1,-1,l3)*calc_coefficient_b(-1,-1,1,l1)*calc_coefficient_b(aph2-1,-1,2,l2)
		   *calc_coefficient_e(aph3-1,0,2,l3-1)) / calc_coefficient_b(aph3-2,0,1,l3);
	    break;
        case 3:
            ans = ((calc_coefficient_rho(index) + calc_coefficient_b(-1,-1,2,l1)*calc_coefficient_e(aph2-1,-1,2,l2-1)*calc_coefficient_theta(index))
		   * calc_coefficient_b(aph3-1,-1,1,l3) + 2*calc_coefficient_b(aph2-2,-1,1,l2)*calc_coefficient_d(aph3-1,-1,l3)
		   *calc_coefficient_b(-1,-1,2,l1)*calc_coefficient_e(aph2-2,-1,1,l2)
		   *calc_coefficient_e(aph3-1,0,2,l3-1)) / (calc_coefficient_b(aph2-2,-1,1,l2)*calc_coefficient_b(aph3-2,0,1,l3));
	    break;
        }
	break;
    }
    return ans;
}

TEMPLATE_TSEM
AFEPack::Point<3> THIS_TSEM::calc_cross_product(AFEPack::Point<3> p1, AFEPack::Point<3> p2)
{// calculate the cross product between 0p1 and 0p2
    AFEPack::Point<3> ans;
    ans[0] =  p1[1] * p2[2] - p1[2] * p2[1];
    ans[1] = -p1[0] * p2[2] + p1[2] * p2[0];
    ans[2] =  p1[0] * p2[1] - p1[1] * p2[0];
    return ans;
}

TEMPLATE_TSEM
int THIS_TSEM::calc_delta(int i, int j)
{
    return (i==j)? 1: 0;
};

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_dipole(std::vector<Vector<valuetype> > &psi, FEMSpace<double, 3> &fem_space, int idx_dim)
{
    valuetype dipole = 0, count;
    std::vector<valuetype> val_psi_qp(n_q_point[2]);
    for (unsigned int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	count = 0;
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(ind_ele).local_to_global(QPoint);
	// summation over real part
	calc_val_qp_onElement(psi[0], val_psi_qp, ind_ele);
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    count += Weight[2][p] * pow(val_psi_qp[p], 2) * q_point[p][idx_dim];
	// summation over imaginary part
	calc_val_qp_onElement(psi[1], val_psi_qp, ind_ele);
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    count += Weight[2][p] * pow(val_psi_qp[p], 2) * q_point[p][idx_dim];
	// add to dipole
	dipole += count * val_volume[ind_ele];
    }
    
    return dipole;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_generalized_jacobi_polynomial(int alpha, int beta, int k, valuetype x)
{// calculate k-th generalized jocobi polynomial with index (alpha, beta) at point x, where $alpha, beta >= -1$
    valuetype ans;
    if (k == 0){
	ans = (valuetype) 1;
	return ans;
    }
    if (alpha + beta == -2){
        if (k == 1) ans = x;
        else ans = ((valuetype) 0.25) * (x-1) * (x+1) * calc_generalized_jacobi_polynomial(1, 1, k-2, x);
	return ans;
    }
    if (alpha == -1){
        ans = (k+beta) * (x-1) / (k*2) * calc_generalized_jacobi_polynomial(1, beta, k-1, x);
	return ans;
    }
    if (beta == -1){
        ans = (k+alpha) * (x+1) / (k*2) * calc_generalized_jacobi_polynomial(alpha, 1, k-1, x);
	return ans;
    }
    valuetype tmp_power = (valuetype) 1;
    ans = (valuetype) 0;
    for (int j = 0; j <= k; ++j){
        valuetype factor = (valuetype) 1;
        for (int i = 0; i < k-j; ++i)
            factor *= ((valuetype) (alpha+j+1 + i)) / (i+1);
        for (int i = 0; i < j; ++i)
            factor *= ((valuetype) (k+alpha+beta+1 + i)) / (i+1);
        ans += factor * tmp_power;
        tmp_power *= (x - 1) / 2;
    }
    return ans;
}

TEMPLATE_TSEM
double THIS_TSEM::calc_inner_product(AFEPack::Point<3> p1, AFEPack::Point<3> p2)
{// calculate the inner product between two vectors, whose enties are given by points p1 and p2
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_length_line(AFEPack::Point<3> &p0, AFEPack::Point<3> &p1)
{
    return sqrt(pow(p0[0]-p1[0], 2) + pow(p0[1]-p1[1], 2) + pow(p0[2]-p1[2], 2));
}

TEMPLATE_TSEM
void THIS_TSEM::calc_sum_dof(std::vector<unsigned int> &sum_dof)
{
    // 0-d
    for (unsigned int ind_geo = 0; ind_geo < n_geometry[0]; ++ind_geo)
	sum_dof[location_actualdof[location_geometry[ind_geo]]] = 1;
    // 1-d
    for (unsigned int ind_geo = 0; ind_geo < n_geometry[1]; ++ind_geo){
	int &pos = location_actualdof[location_geometry[ind_geo + n_geometry[0]]];
	for (int l = 2; l <= M; ++l)
	    sum_dof[pos + l-2] = l;
    }
    // 2-d
    for (unsigned int ind_geo = 0; ind_geo < n_geometry[2]; ++ind_geo){
	int &pos = location_actualdof[location_geometry[ind_geo + n_geometry[0] + n_geometry[1]]];
	for (unsigned int l1 = 2; l1 <= M; ++l1)
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		sum_dof[pos + ind_index] = l1 + l2;
	    }
    }
    // 3-d
    for (unsigned int ind_geo = 0; ind_geo < n_geometry[3]; ++ind_geo){
	int &pos = location_actualdof[location_geometry[ind_geo + n_geometry[0] + n_geometry[1] + n_geometry[2]]];
	for (unsigned int l1 = 2; l1 <= M; ++l1)
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2)
		for (unsigned int l3 = 1; l3 <= M-l1-l2; ++l3){
		    Multiindex<3> index_now = Unitary_Multiindex[0]*(l1-2) + Unitary_Multiindex[1]*(l2-1) + Unitary_Multiindex[2]*(l3-1);
		    sum_dof[pos + correspondence.index2number(index_now)] = l1 + l2 + l3;
		}
    }
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_volume_tetrahedron(AFEPack::Point<3> &p0, AFEPack::Point<3> &p1, AFEPack::Point<3> &p2, AFEPack::Point<3> &p3)
{
    return fabs(((p1[0] - p0[0])*(p2[1] - p0[1])*(p3[2] - p0[2])
		 + (p1[1] - p0[1])*(p2[2] - p0[2])*(p3[0] - p0[0])
		 + (p1[2] - p0[2])*(p2[0] - p0[0])*(p3[1] - p0[1])
		 - (p1[0] - p0[0])*(p2[2] - p0[2])*(p3[1] - p0[1])
		 - (p1[1] - p0[1])*(p2[0] - p0[0])*(p3[2] - p0[2])
                 - (p1[2] - p0[2])*(p2[1] - p0[1])*(p3[0] - p0[0]))/6.);
}

TEMPLATE_TSEM
void THIS_TSEM::impose_mask(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*), std::vector<Vector<valuetype> > &psi)
{    
    Vector<double> interp_re(n_dof_total), interp_im(n_dof_total);

    // traverse 0 dimensional geometry, correspond to vertex in fem dof
    for (unsigned int ind_point = 0; ind_point < n_geometry[0]; ++ind_point){
	int &pos = location_actualdof[location_geometry[ind_point]];
	interp_re(pos) = psi[0](pos) * func(mesh.point(ind_point));
	interp_im(pos) = psi[1](pos) * func(mesh.point(ind_point));
    }

    // traverse 1 dimensional geometry
    std::vector<double> val_interp1(n_q_point[0]), val_psi_re1(n_q_point[0]), val_psi_im1(n_q_point[0]);
    for (unsigned int ind_edge = 0; ind_edge < n_geometry[1]; ++ind_edge){
	int &ind_point_s = number_node[0][ind_edge][0], &ind_point_e = number_node[0][ind_edge][1];
	int &pos_s = location_actualdof[location_geometry[ind_point_s]], &pos_e = location_actualdof[location_geometry[ind_point_e]];
	valuetype p_s_re = psi[0](pos_s), p_e_re = psi[0](pos_e);
	valuetype p_s_im = psi[1](pos_s), p_e_im = psi[1](pos_e);
        valuetype c_s_re = interp_re(pos_s), c_e_re = interp_re(pos_e);
        valuetype c_s_im = interp_im(pos_s), c_e_im = interp_im(pos_e);
	for (unsigned int p = 0; p < n_q_point[0]; ++p){
	    AFEPack::Point<3> p_tmp = mesh.point(ind_point_s), p_ttmp = mesh.point(ind_point_e);
	    p_tmp  *= QPoint_Barycentric[0][p][1];
	    p_ttmp *= QPoint_Barycentric[0][p][0];
	    p_tmp  += p_ttmp;
	    val_interp1[p] = func(p_tmp);
	}
	for (unsigned int p = 0; p < n_q_point[0]; ++p){
	    val_psi_re1[p] = (p_s_re*val_interp1[p]-c_s_re)*QPoint_Barycentric[0][p][1] + (p_e_re*val_interp1[p]-c_e_re)*QPoint_Barycentric[0][p][0];
	    val_psi_im1[p] = (p_s_im*val_interp1[p]-c_s_im)*QPoint_Barycentric[0][p][1] + (p_e_im*val_interp1[p]-c_e_im)*QPoint_Barycentric[0][p][0];
	    valuetype count_re = 0, count_im = 0;
	    for (int l = 2; l <= M; ++l){
		int &location_dof = location_actualdof[location_geometry[ind_edge+n_geometry[0]]];
		count_re += psi[0](location_dof + l-2) * basis_value[0][p][l-2];
		count_im += psi[1](location_dof + l-2) * basis_value[0][p][l-2];
	    }
	    val_psi_re1[p] += count_re * val_interp1[p];
	    val_psi_im1[p] += count_im * val_interp1[p];
	}
        for (int l = 2; l <= M; ++l){ // dof locate on this edge
            int location_dof = location_actualdof[location_geometry[ind_edge+n_geometry[0]]] + l-2; // position of this dof
            valuetype count_re = 0, count_im = 0;
            for (unsigned int p = 0; p < n_q_point[0]; ++p){
                count_re += val_psi_re1[p] * basis_value_interp[0][p][l-2] * Weight[0][p];
                count_im += val_psi_im1[p] * basis_value_interp[0][p][l-2] * Weight[0][p];
	    }
            interp_re(location_dof) = count_re * -l * (2*l - 1) / (2 * (l-1));
            interp_im(location_dof) = count_im * -l * (2*l - 1) / (2 * (l-1));
        }
    }
    
    // traverse 2 dimensional geometry
    std::vector<double> val_interp2(n_q_point[1]), val_psi_re2(n_q_point[1]), val_psi_im2(n_q_point[1]);
    for (unsigned int ind_face = 0; ind_face < n_geometry[2]; ++ind_face){
	int location_dof_edge[3]; // start location of dof on each edge
	for (unsigned int ind_edge = 0; ind_edge < 3; ++ind_edge)
	    location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_face][ind_edge] + n_geometry[0]]];
	std::vector<std::vector<valuetype> > val_edgedof_re(3, std::vector<valuetype> (n_dof_edge));
	std::vector<std::vector<valuetype> > val_edgedof_im(3, std::vector<valuetype> (n_dof_edge));
	std::vector<std::vector<valuetype> > val_interp_edgedof_re(3, std::vector<valuetype> (n_dof_edge));
	std::vector<std::vector<valuetype> > val_interp_edgedof_im(3, std::vector<valuetype> (n_dof_edge));
	for (unsigned int ind_e = 0; ind_e < 3; ++ind_e)
	    for (unsigned int l = 2; l <= M; ++l){
		val_edgedof_re[ind_e][l-2] = psi[0](location_dof_edge[ind_e] + l-2);
		val_edgedof_im[ind_e][l-2] = psi[1](location_dof_edge[ind_e] + l-2);
		val_interp_edgedof_re[ind_e][l-2] = interp_re(location_dof_edge[ind_e] + l-2);
		val_interp_edgedof_im[ind_e][l-2] = interp_im(location_dof_edge[ind_e] + l-2);
		if (!flag_sameorder_edgeonface[ind_face][ind_e] && l % 2 == 1){
		    val_edgedof_re[ind_e][l-2] *= -1;
		    val_edgedof_im[ind_e][l-2] *= -1;
		    val_interp_edgedof_re[ind_e][l-2] *= -1;
		    val_interp_edgedof_im[ind_e][l-2] *= -1;
		}
	    }
	for (unsigned int p = 0; p < n_q_point[1]; ++p){
	    val_psi_re2[p] = val_psi_im2[p] = 0;
	    valuetype &xp = QPoint_Barycentric[1][p][0], &yp = QPoint_Barycentric[1][p][1], &rp = QPoint_Barycentric[1][p][2];
	    for (int ind = 0; ind < 3; ++ind){ // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
		int &pos = location_actualdof[location_geometry[number_node[1][ind_face][ind]]];
		val_psi_re2[p] += psi[0](pos) * QPoint_Barycentric[1][p][(ind+2)%3];
		val_psi_im2[p] += psi[1](pos) * QPoint_Barycentric[1][p][(ind+2)%3];
	    }
	    for (unsigned int l = 2; l <= M; ++l){ // substract the contribution from edge
		val_psi_re2[p] -= 2 * (xp*yp * val_edgedof_re[0][l-2] * basis_value_addition[p][l-2][1] +
				       yp*rp * val_edgedof_re[1][l-2] * basis_value_addition[p][l-2][1] +
				       xp*rp * val_edgedof_re[2][l-2] * basis_value_addition[p][l-2][0]);
		val_psi_im2[p] -= 2 * (xp*yp * val_edgedof_im[0][l-2] * basis_value_addition[p][l-2][1] +
				       yp*rp * val_edgedof_im[1][l-2] * basis_value_addition[p][l-2][1] +
				       xp*rp * val_edgedof_im[2][l-2] * basis_value_addition[p][l-2][0]);
	    }
	    for (int l1 = 2; l1 <= M; ++l1)
		for (int l2 = 1; l2 <= M-l1; ++l2){
		    int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
		    int location_dof = location_actualdof[location_geometry[ind_face+n_geometry[0]+n_geometry[1]]] + ind_index;
		    val_psi_re2[p] += psi[0](location_dof) * basis_value[1][p][ind_index];
		    val_psi_im2[p] += psi[1](location_dof) * basis_value[1][p][ind_index];
		}
	}
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2){
		for (unsigned int p = 0; p < n_q_point[1]; ++p){
    		    AFEPack::Point<3> p_tmp;
    		    for (unsigned int ind = 0; ind < 3; ++ind) p_tmp[ind] = 0;
    		    for (int ind_p = 0; ind_p <= 2; ++ind_p){ // vertex[0]=(0,0)<->[2]=1-x-y; vertex[1]=(1,0)<->[0]=x; vertex[2]=(0,1)<->[1]=y
    			int &ind_vertex = number_node[1][ind_face][ind_p];
    			AFEPack::Point<3> p_ttmp = mesh.point(ind_vertex);
    			p_ttmp *= QPoint_Barycentric[1][p][(ind_p+2)%3];
    			p_tmp  += p_ttmp;
    		    }
    		    val_interp2[p] = func(p_tmp);
    		}
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                int location_dof = location_actualdof[location_geometry[ind_face+n_geometry[0]+n_geometry[1]]] + ind_index;
                valuetype cnt_re = 0, cnt_im = 0;
                for (unsigned int p = 0; p < n_q_point[1]; ++p){
                    valuetype &xp = QPoint_Barycentric[1][p][0], &yp = QPoint_Barycentric[1][p][1], &rp = QPoint_Barycentric[1][p][2];
                    valuetype count_re = val_psi_re2[p]*val_interp2[p], count_im = val_psi_im2[p]*val_interp2[p];
                    for (int ind = 0; ind < 3; ++ind){
                        count_re -= interp_re(location_actualdof[location_geometry[number_node[1][ind_face][ind]]]) * QPoint_Barycentric[1][p][(ind+2)%3];
                        count_im -= interp_im(location_actualdof[location_geometry[number_node[1][ind_face][ind]]]) * QPoint_Barycentric[1][p][(ind+2)%3];
		    }
                    for (unsigned int l = 2; l <= M; ++l){ // substract the contribution from edge
                        count_re += 2 * (xp*yp * val_interp_edgedof_re[0][l-2] * basis_value_addition[p][l-2][1] +
					 yp*rp * val_interp_edgedof_re[1][l-2] * basis_value_addition[p][l-2][1] +
					 xp*rp * val_interp_edgedof_re[2][l-2] * basis_value_addition[p][l-2][0]);
                        count_im += 2 * (xp*yp * val_interp_edgedof_im[0][l-2] * basis_value_addition[p][l-2][1] +
					 yp*rp * val_interp_edgedof_im[1][l-2] * basis_value_addition[p][l-2][1] +
					 xp*rp * val_interp_edgedof_im[2][l-2] * basis_value_addition[p][l-2][0]);
                    }
                    cnt_re += Weight[1][p] * count_re * basis_value_interp[1][p][ind_index];
                    cnt_im += Weight[1][p] * count_im * basis_value_interp[1][p][ind_index];
                }
                interp_re(location_dof) = cnt_re * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
                interp_im(location_dof) = cnt_im * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
            }
    }
    
    // traverse fem element, assign value for interior dof
    std::vector<valuetype> val_interp3(n_q_point[2]), val_psi_re3(n_q_point[2]), val_psi_im3(n_q_point[2]);
    for (unsigned int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    val_interp3[p] = func(fem_space.element(ind_ele).local_to_global(QPoint[p]));
	calc_val_qp_onElement(psi[0], val_psi_re3, ind_ele);
	calc_val_qp_onElement(psi[1], val_psi_im3, ind_ele);
        // get local coefficients for vertex, edge and face
        std::vector<valuetype> val_coef_re(n_index, 0), val_coef_im(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            Multiindex<3> index_now = correspondence.number2index(ind_index);
            if (index_now.index[0] >= 2 && index_now.index[1] >= 1 && index_now.index[2] >= 1)
                continue;
            for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz){
                val_coef_re[ind_index] += interp_re(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
                val_coef_im[ind_index] += interp_im(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
	    }
        }
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2)
                for (int l3 = 1; l3 <= M-l1-l2; ++l3){
                    int ind_index = correspondence.index2number(Unitary_Multiindex[0] * (l1-2) + Unitary_Multiindex[1] * (l2-1) + Unitary_Multiindex[2] * (l3-1));
                    int location_dof = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]] + ind_index;
                    valuetype cnt_re = 0, cnt_im = 0;
                    for (unsigned int p = 0; p < n_q_point[2]; ++p){
                        valuetype count_re = val_psi_re3[p]*val_interp3[p], count_im = val_psi_im3[p]*val_interp3[p];
                        for (int j = 0; j < n_index; ++j){ // traverse multiindex for vertex, edge and face
                            Multiindex<3> index_tmp = correspondence.number2index(j);
                            if (index_tmp.index[0] >= 2 && index_tmp.index[1] >= 1 && index_tmp.index[2] >= 1) // the interior model function are orthogonal
                                continue;
                            count_re -= val_coef_re[j] * basis_value_actual[p][j];
                            count_im -= val_coef_im[j] * basis_value_actual[p][j];
                        }
                        cnt_re += Weight[2][p] * count_re *  basis_value_interp[2][p][ind_index];
                        cnt_im += Weight[2][p] * count_im *  basis_value_interp[2][p][ind_index];
                    }
                    // for interior dof, global coefficient is exactly local one
                    interp_re(location_dof) = cnt_re * -l1 * (2*l1-1) * (2*l1+2*l2-1) * (2*l1+2*l2+2*l3-1) / (pow(2, 2*l1+l2-5) * (l1-1) * 6);
                    interp_im(location_dof) = cnt_im * -l1 * (2*l1-1) * (2*l1+2*l2-1) * (2*l1+2*l2+2*l3-1) / (pow(2, 2*l1+l2-5) * (l1-1) * 6);
                }
    }

    // double dif = 0;
    // for (int idx = 0; idx < mesh.n_geometry(1); ++idx){
    // 	// int &pos = location_actualdof[location_geometry[idx]];
    //     // for (unsigned int l = 2; l <= M; ++l){ // dof locate on this edge
    //     //     int pos = location_actualdof[location_geometry[idx+n_geometry[0]]] + l-2; // position of this dof
    // 	//     dif += fabs(psi[0](pos) - interp_re(pos));
    // 	// }
    // 	for (int l1 = 2; l1 <= M; ++l1)
    // 	    for (int l2 = 1; l2 <= M-l1; ++l2){
    //             int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
    //             int pos = location_actualdof[location_geometry[idx+n_geometry[0]+n_geometry[1]]] + ind_index;
    // 		dif += fabs(psi[0](pos) - interp_re(pos));
    // 	    }
    // }
    // Vector<valuetype> tmpd(n_dof_total);
    // tmpd = psi[0]; tmpd -= interp_re;
    // std::cout << "dif for all dof: " << tmpd.l1_norm() << '\n';

    psi[0] = interp_re; psi[1] = interp_im;
    
    std::cerr << "impose mask function\n";
}

#undef TEMPLATE_TSEM
#undef THIS_TSEM
