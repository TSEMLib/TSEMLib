#include "../include/TetrahedralSEM.h"
// This file contains functions for calculating the TSEM solution value at the specified point
// Functions:
// calc_basis_gradient()
// calc_coef_onElement()
// calc_density()
// calc_dependence()
// calc_val_gradient_qp_onElement()
// calc_val_inElement()
// calc_val_point()
// calc_val_point_insideElement()
// calc_val_qp_onElement()


#define TEMPLATE_TSEM template<typename valuetype>
#define THIS_TSEM TSEM<valuetype>


TEMPLATE_TSEM
void THIS_TSEM::calc_basis_gradient(Vector<valuetype> &src, std::vector<std::vector<std::vector<double> > > &dst, RegularMesh<3> &mesh, unsigned int idx_ele)
{ // calculate gradient of basis function on the idx_ele-th element in mesh
    // initialize
    int &n_qp = n_q_point[2];
    dst.resize(n_qp, std::vector<std::vector<double> > (n_index, std::vector<double> (3, 0)));
    
    // calculate
    // prepare coordinates of element vertices
    valuetype xj_k[3][4];
    for (unsigned int k = 0; k < 4; ++k){
	int &ind_node = index_geometry_onelement[idx_ele][0][k];
	for (unsigned int j = 0; j < 3; ++j) xj_k[j][k] = mesh.point(ind_node)[j];
    }
    valuetype d_jk[3][3];
    for (unsigned int j = 0; j < 3; ++j)
	for (unsigned int k = 0; k < 3; ++k)
	    d_jk[j][k] = xj_k[j][k+1] - xj_k[j][0];
    valuetype val_V = d_jk[0][0]*d_jk[1][1]*d_jk[2][2] + d_jk[0][1]*d_jk[1][2]*d_jk[2][0] + d_jk[0][2]*d_jk[1][0]*d_jk[2][1]
	- d_jk[0][0]*d_jk[1][2]*d_jk[2][1] - d_jk[0][1]*d_jk[1][0]*d_jk[2][2] - d_jk[0][2]*d_jk[1][1]*d_jk[2][0];
    // calculate \partial x_j/\partial \hat{x}_k
    valuetype parxj_parxk[3][3];
    parxj_parxk[0][0] = (d_jk[1][1]*d_jk[2][2] - d_jk[1][2]*d_jk[2][1]) / val_V;
    parxj_parxk[0][1] = (d_jk[0][2]*d_jk[2][1] - d_jk[0][1]*d_jk[2][2]) / val_V;
    parxj_parxk[0][2] = (d_jk[0][1]*d_jk[1][2] - d_jk[0][2]*d_jk[1][1]) / val_V;
    parxj_parxk[1][0] = (d_jk[1][2]*d_jk[2][0] - d_jk[1][0]*d_jk[2][2]) / val_V;
    parxj_parxk[1][1] = (d_jk[0][0]*d_jk[2][2] - d_jk[0][2]*d_jk[2][0]) / val_V;
    parxj_parxk[1][2] = (d_jk[0][2]*d_jk[1][0] - d_jk[0][0]*d_jk[1][2]) / val_V;
    parxj_parxk[2][0] = (d_jk[1][0]*d_jk[2][1] - d_jk[1][1]*d_jk[2][0]) / val_V;
    parxj_parxk[2][1] = (d_jk[0][1]*d_jk[2][0] - d_jk[0][0]*d_jk[2][1]) / val_V;
    parxj_parxk[2][2] = (d_jk[0][0]*d_jk[1][1] - d_jk[0][1]*d_jk[1][0]) / val_V;
    // traverse quadrature point, assign corresponding value
    for (unsigned int p = 0; p < n_qp; ++p){
	// recover the value of basis gradient
	for (unsigned int idx_index = 0; idx_index < n_index; ++idx_index)
	    for (unsigned int j = 0; j < 3; ++j)
		for (unsigned int k = 0; k < 3; ++k)
		    dst[p][idx_index][j] += basis_gradient_actual[p][idx_index][k] * parxj_parxk[k][j];
    }    
}

TEMPLATE_TSEM
void THIS_TSEM::calc_coef_onElement(Vector<valuetype> &src, std::vector<valuetype> &dst, unsigned int ind_ele)
{ // calculate coeffcieint of sem solution on ind_ele-th element from source src to vector dst, suppose dst has size n_index
    // dst.resize(n_index);
    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	dst[ind_index] = 0;
    for (int ind_index = 0; ind_index < n_index; ++ind_index)
	for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
	    dst[ind_index] += src(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
		* transform_val_global2local[ind_ele][ind_index][ind_nnz];
}

TEMPLATE_TSEM
void THIS_TSEM::calc_density(std::vector<Vector<valuetype> > &src, std::vector<valuetype> &dst, std::vector<valuetype> &n_occupation, unsigned int idx_ele)
{
    // initialize
    for (unsigned int idx_qp = 0; idx_qp < dst.size(); ++idx_qp) dst[idx_qp] = 0;
    // assign density from psi on several orbitals with occupation number n_occupation
    std::vector<valuetype> psi_local(dst.size());
    for (unsigned int idx_orb = 0; idx_orb < src.size(); ++idx_orb){
	calc_val_qp_onElement(src[idx_orb], psi_local, idx_ele);
	for (unsigned int idx_qp = 0; idx_qp < dst.size(); ++idx_qp) dst[idx_qp] += n_occupation[idx_orb] * pow(psi_local[idx_qp], 2);
    }
}

TEMPLATE_TSEM
void THIS_TSEM::calc_density(std::vector<Vector<valuetype> > &src, std::vector<std::vector<valuetype> > &dst, std::vector<valuetype> &n_occupation)
{
    for (unsigned int idx_ele = 0; idx_ele < n_element; ++idx_ele)
	calc_density(src, dst[idx_ele], n_occupation, idx_ele);
}

TEMPLATE_TSEM
void THIS_TSEM::calc_dependence(RegularMesh<3> &mesh, std::vector<AFEPack::Point<3> > &point, std::vector<int> &n_dep, std::vector<std::vector<int> > &dep_idx, std::vector<std::vector<double> > &dep_para)
{
    int np = point.size();
    for (unsigned int idx_point = 0; idx_point < np; ++idx_point){
	AFEPack::Point<3> &pos = point[idx_point];
	int &nd = n_dep[idx_point];
	std::vector<int> &di = dep_idx[idx_point];
	std::vector<double> &dp = dep_para[idx_point];
	// check 0d geoemtries
	for (unsigned int idx_p = 0; idx_p < n_geometry[0]; ++idx_p)
	    if (calc_length_line(mesh.point(idx_p), pos) < tol_zero){
		nd = 1; di.resize(nd); dp.resize(nd);
		di[0] = location_actualdof[location_geometry[idx_p]];
		dp[0] = 1;
	    }
	// check 1-d geometries
	for (unsigned int ind_e = 0; ind_e < n_geometry[1]; ++ind_e){
	    int &ind_point0 = number_node[0][ind_e][0], &ind_point1 = number_node[0][ind_e][1];
	    valuetype dis = calc_length_line(mesh.point(ind_point0), mesh.point(ind_point1));
	    valuetype dis0 = calc_length_line(mesh.point(ind_point0), pos), dis1 = calc_length_line(mesh.point(ind_point1), pos);
	    if (fabs(dis0+dis1 - dis) > tol_zero) continue;
	    
	    int &loc0 = location_actualdof[location_geometry[ind_point0]], loc1 = location_actualdof[location_geometry[ind_point1]];
	    int &loc = location_actualdof[location_geometry[n_geometry[0] + ind_e]];
	
	    valuetype x = dis0 / dis, r = 1 - x, xi = 2 * x - 1;
	    valuetype Jxi[M+1]; Jxi[0] = 1; Jxi[1] = xi;
	    for (unsigned int l1 = 1; l1 < M; ++l1)
		Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1]) / calc_coefficient_a(-1, -1, 1, l1);

	    nd = M+1; // 2 + n_dof_geometry[1]
	    di.resize(nd); dp.resize(nd);
	    di[0] = loc0; dp[0] = r; di[1] = loc1; dp[1] = x;
	    for (unsigned int l1 = 2; l1 <= M; ++l1){
		di[l1] = loc+l1-2; dp[l1] = 2*Jxi[l1]; }
	}    
	// check 2-d geometries
	for (unsigned int ind_f = 0; ind_f < n_geometry[2]; ++ind_f){
	    int &ind_point0 = number_node[1][ind_f][0], &ind_point1 = number_node[1][ind_f][1], &ind_point2 = number_node[1][ind_f][2];
	    valuetype area = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), mesh.point(ind_point2));
	    valuetype area0 = calc_area_triangle(pos, mesh.point(ind_point1), mesh.point(ind_point2));
	    valuetype area1 = calc_area_triangle(mesh.point(ind_point0), pos, mesh.point(ind_point2));
	    valuetype area2 = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), pos);
	    if (fabs(area0+area1+area2 - area) > tol_zero) continue;
	    
	    int &ind_edge0 = number_edge[ind_f][0], &ind_edge1 = number_edge[ind_f][1], &ind_edge2 = number_edge[ind_f][2];
	    int &loc = location_actualdof[location_geometry[n_geometry[0] + n_geometry[1] + ind_f]];
	    int &locp0 = location_actualdof[location_geometry[ind_point0]], &locp1 = location_actualdof[location_geometry[ind_point1]], &locp2 = location_actualdof[location_geometry[ind_point2]];
	    int &loce0 = location_actualdof[location_geometry[n_geometry[0] + ind_edge0]];
	    int &loce1 = location_actualdof[location_geometry[n_geometry[0] + ind_edge1]];
	    int &loce2 = location_actualdof[location_geometry[n_geometry[0] + ind_edge2]];
	
	    std::vector<valuetype> bas_val(n_dof_geometry[2], 0);
	    valuetype x = area1 / area, y = area2 / area, r = 1 - x - y;
	    valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
	    valuetype Jxi[M+1], Jeta[M+1]; Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
	    for (unsigned int l1 = 1; l1 < M; ++l1)
		Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1]) / calc_coefficient_a(-1, -1, 1, l1);
	    for (unsigned int l1 = 2; l1 <= M-1; ++l1){
		int aph2 = 2 * l1 - 1;
		Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
		for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		    Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1]) / calc_coefficient_a( aph2, -1, 1, l2);
		for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		    unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		    bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
		}
	    }
	    std::vector<valuetype> bas_val_add_xi( n_dof_geometry[1], 0);
	    std::vector<valuetype> bas_val_add_eta(n_dof_geometry[1], 0);
	    Jxi[0] = Jeta[0] = 1;
	    Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
	    Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
	    for (unsigned int l = 1; l < M; ++l){
		Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
		Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
	    }
	    for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
		bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
		bas_val_add_eta[l] = Jeta[l];
	    }

	    nd = 3 + n_dof_geometry[1]*3 + n_dof_geometry[2];
	    di.resize(nd); dp.resize(nd);
	    di[0] = locp0; dp[0] = r; di[1] = locp1; dp[1] = x; di[2] = locp2; dp[2] = y;
	    for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
		di[l*3+3] = loce0+l; dp[l*3+3] = -2*x*y*bas_val_add_eta[l];
		di[l*3+4] = loce1+l; dp[l*3+4] = -2*y*r*bas_val_add_eta[l];
		di[l*3+5] = loce2+l; dp[l*3+5] = -2*x*r*bas_val_add_xi[l];
	    }
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		    unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		    di[n_dof_geometry[1]+3+ind_index] = loc + ind_index;
		    dp[n_dof_geometry[1]+3+ind_index] = bas_val[ind_index];
		}
	}    
	// check 3-d geometries
	for (unsigned int ind_e = 0; ind_e < n_geometry[3]; ++ind_e){
	    int &ind_v0 = index_geometry_onelement[ind_e][0][0], &ind_v1 = index_geometry_onelement[ind_e][0][1], &ind_v2 = index_geometry_onelement[ind_e][0][2], &ind_v3 = index_geometry_onelement[ind_e][0][3];
	    valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
	    valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
	    if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) continue;

	    valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
	    valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
	    std::vector<valuetype> basis(n_index);
	    // calculate immediate variable
	    valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
	    Jxi[0] = Jeta[0] = Jzeta[0] = 1; Jxi[1]  = xi;
	    for (unsigned int l1 = 1; l1 < M; ++l1)
		Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])	/ calc_coefficient_a(-1, -1, 1, l1);
	    for (unsigned int l1 = 0; l1 <= M; ++l1){
		int aph2 = 2 * l1 - 1;
		Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
		for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		    Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1]) / calc_coefficient_a( aph2, -1, 1, l2);
		for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		    int aph3 = 2 * l1 + 2 * l2 - 1;
		    Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
		    for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
			Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1]) / calc_coefficient_a( aph3, -1, 1, l3);
		    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
			Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
			int ind_index = correspondence.index2number(index_now);
			basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
		    }
		}
	    }
	    std::vector<valuetype> basis_actual(n_index, 0);
	    for (int ind_index = 0; ind_index < n_index; ++ind_index)
		for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		    basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];

	    // specially, nd = -1 for flag_insideElement, di[0] = ind_e
	    nd = -1; di.resize(1); di[0] = ind_e;
	    dp.resize(n_index);
	    for (unsigned int idx_index = 0; idx_index < n_index; ++idx_index)
		dp[idx_index] = basis_actual[idx_index];
	}
    }
}

TEMPLATE_TSEM
void THIS_TSEM::calc_dependence(RegularMesh<3> &mesh, std::vector<AFEPack::Point<3> > &point, std::vector<int> &idx_ele, std::vector<std::vector<double> > &dep_basis, std::vector<std::vector<std::vector<double> > > &dep_gradient)
{ // suppose searched point is inside some element
    int np = point.size();
    for (unsigned int idx_point = 0; idx_point < np; ++idx_point){
	AFEPack::Point<3> &pos = point[idx_point];
	for (unsigned int ind_e = 0; ind_e < n_geometry[3]; ++ind_e){
	    int &ind_v0 = index_geometry_onelement[ind_e][0][0], &ind_v1 = index_geometry_onelement[ind_e][0][1], &ind_v2 = index_geometry_onelement[ind_e][0][2], &ind_v3 = index_geometry_onelement[ind_e][0][3];
	    valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
	    valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
	    if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) continue;

	    idx_ele[idx_point] = ind_e;

	    // calculate Jacobi polynomial value for basis expression
	    valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
	    valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
	    std::vector<valuetype> basis(n_index);
	    // calculate immediate variable
	    valuetype Jxi_[M+1], Jeta_[M+1], Jzeta_[M+1];
	    Jxi_[0] = Jeta_[0] = Jzeta_[0] = 1; Jxi_[1]  = xi;
	    for (unsigned int l1 = 1; l1 < M; ++l1)
		Jxi_[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi_[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi_[l1-1])	/ calc_coefficient_a(-1, -1, 1, l1);
	    for (unsigned int l1 = 0; l1 <= M; ++l1){
		int aph2 = 2 * l1 - 1;
		Jeta_[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
		for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		    Jeta_[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta_[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta_[l2-1]) / calc_coefficient_a( aph2, -1, 1, l2);
		for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		    int aph3 = 2 * l1 + 2 * l2 - 1;
		    Jzeta_[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
		    for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
			Jzeta_[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta_[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta_[l3-1]) / calc_coefficient_a( aph3, -1, 1, l3);
		    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
			Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
			int ind_index = correspondence.index2number(index_now);
			basis[ind_index] = pow(1-y-z, l1) * Jxi_[l1] * pow(1-z, l2) * Jeta_[l2] * Jzeta_[l3];
		    }
		}
	    }
	    std::vector<valuetype> basis_actual(n_index, 0);
	    for (int ind_index = 0; ind_index < n_index; ++ind_index)
		for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		    basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];

	    for (unsigned int idx_index = 0; idx_index < n_index; ++idx_index)
		dep_basis[idx_point][idx_index] = basis_actual[idx_index];

	    // calculate Jacobi polynomial value for gradient expression
	    valuetype coeff_D1[n_index];
	    valuetype coeff_D2[n_index][2]; // p = 0, 1
	    valuetype coeff_D3[n_index][4]; // (p, q) = (0, 0) (1, 0) (0, 1) (1, 1)
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index){
		Multiindex<3> index_now = correspondence.number2index(ind_index);
		coeff_D1[ind_index] = calc_coefficient_D(0, 0, index_now);
		coeff_D2[ind_index][0] = calc_coefficient_D(3, 0, index_now);
		coeff_D2[ind_index][1] = calc_coefficient_D(3, 1, index_now);
		coeff_D3[ind_index][0] = calc_coefficient_D(5, 0, index_now);
		coeff_D3[ind_index][1] = calc_coefficient_D(5, 1, index_now);
		coeff_D3[ind_index][2] = calc_coefficient_D(5, 2, index_now);
		coeff_D3[ind_index][3] = calc_coefficient_D(5, 3, index_now);
	    }
	    std::vector<std::vector<valuetype> > curlicue_J(3, std::vector<valuetype> (n_index));
	    // assign curlicue J by calculating generalize Jacobi polynomials
	    std::vector<std::vector<double> > basisGradient(n_index, std::vector<double> (3));
	    valuetype x1 = x, x2 = y, x3 = z;
	    valuetype t1 = 1-x2-x3, t2 = 1-x3;
	    xi = 2*x1/t1-1, eta = 2*x2/t2-1, zeta = 2*x3-1;
	    valuetype pow_t1[M+1], pow_t2[M+1];
	    pow_t1[0] = (valuetype) 1; pow_t2[0] = (valuetype) 1;
	    for (unsigned int l = 1; l <= M; ++l){
		pow_t1[l] = pow_t1[l-1] * t1; pow_t2[l] = pow_t2[l-1] * t2; }
	    valuetype Jxi  [2][M+1]; // [0][]: ^{0,0};          [1][]: ^{0,-1}
	    valuetype Jeta [3][M+1]; // [0][]: ^{2l1+1,-1};     [1][]: ^{2l1,0};     [2][]: ^{2l1,-1}
	    valuetype Jzeta[2][M+1]; // [0][]: ^{2l1+2l2+1,-1}; [1][]: ^{2l1+2l2,0}
	    Jxi[0][0] = Jxi[1][0] = (valuetype) 1;
	    Jeta[0][0] = Jeta[1][0] = Jeta[2][0] = (valuetype) 1;
	    Jzeta[0][0] = Jzeta[1][0] = (valuetype) 1;
	    Jxi[0][1] = calc_generalized_jacobi_polynomial(0,  0, 1, xi);
	    Jxi[1][1] = calc_generalized_jacobi_polynomial(0, -1, 1, xi);
	    for (unsigned int l = 1; l < M; ++l){
		Jxi[0][l+1] = ((xi - calc_coefficient_a(0, 0,2,l)) * Jxi[0][l] - calc_coefficient_a(0, 0,3,l) * Jxi[0][l-1]) / calc_coefficient_a(0, 0,1,l);
		Jxi[1][l+1] = ((xi - calc_coefficient_a(0,-1,2,l)) * Jxi[1][l] - calc_coefficient_a(0,-1,3,l) * Jxi[1][l-1]) / calc_coefficient_a(0,-1,1,l);
	    }
	    for (unsigned int l1 = 0; l1 <= M; ++l1){
		int aph2_0 = 2*l1+1, aph2_1 = 2*l1;
		Jeta[0][1] = calc_generalized_jacobi_polynomial(aph2_0, -1, 1, eta);
		Jeta[1][1] = calc_generalized_jacobi_polynomial(aph2_1,  0, 1, eta);
		Jeta[2][1] = calc_generalized_jacobi_polynomial(aph2_1, -1, 1, eta);
		for (unsigned int l = 1; l < M-l1; ++l){
		    Jeta[0][l+1] = ((eta - calc_coefficient_a(aph2_0,-1,2,l)) * Jeta[0][l] - calc_coefficient_a(aph2_0,-1,3,l) * Jeta[0][l-1]) / calc_coefficient_a(aph2_0,-1,1,l);
		    Jeta[1][l+1] = ((eta - calc_coefficient_a(aph2_1, 0,2,l)) * Jeta[1][l] - calc_coefficient_a(aph2_1, 0,3,l) * Jeta[1][l-1]) / calc_coefficient_a(aph2_1, 0,1,l);
		    Jeta[2][l+1] = ((eta - calc_coefficient_a(aph2_1,-1,2,l)) * Jeta[2][l] - calc_coefficient_a(aph2_1,-1,3,l) * Jeta[2][l-1]) / calc_coefficient_a(aph2_1,-1,1,l);
		}
		for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		    int aph3_0 = 2*l1+2*l2+1, aph3_1 = 2*l1+2*l2;
		    Jzeta[0][1] = calc_generalized_jacobi_polynomial(aph3_0, -1, 1, zeta);
		    Jzeta[1][1] = calc_generalized_jacobi_polynomial(aph3_1,  0, 1, zeta);
		    for (unsigned int l = 1; l < M-l1-l2; ++l){
			Jzeta[0][l+1] = ((zeta-calc_coefficient_a(aph3_0,-1,2,l)) * Jzeta[0][l] - calc_coefficient_a(aph3_0,-1,3,l) * Jzeta[0][l-1]) / calc_coefficient_a(aph3_0,-1,1,l);
			Jzeta[1][l+1] = ((zeta-calc_coefficient_a(aph3_1, 0,2,l)) * Jzeta[1][l] - calc_coefficient_a(aph3_1, 0,3,l) * Jzeta[1][l-1]) / calc_coefficient_a(aph3_1, 0,1,l);
		    }
		    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
			Multiindex<3> index_now = Unitary_Multiindex[0]*l1 + Unitary_Multiindex[1]*l2 + Unitary_Multiindex[2]*l3;
			int ind_index = correspondence.index2number(index_now);
			curlicue_J[0][ind_index] = pow_t1[l1]*Jxi[0][l1] * pow_t2[l2]*Jeta[0][l2] * Jzeta[0][l3];
			curlicue_J[1][ind_index] = pow_t1[l1]*Jxi[1][l1] * pow_t2[l2]*Jeta[1][l2] * Jzeta[0][l3];
			curlicue_J[2][ind_index] = pow_t1[l1]*Jxi[1][l1] * pow_t2[l2]*Jeta[2][l2] * Jzeta[1][l3];
		    }
		}
	    }
	    // assign basis_gradient[p][][]
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index){
		Multiindex<3> index_tmp = correspondence.number2index(ind_index);
		// 1st component
		index_tmp.index[0]--;
		if (0 <= index_tmp.index[0])
		    basisGradient[ind_index][0] += coeff_D1[ind_index]    * curlicue_J[0][correspondence.index2number(index_tmp)];
		index_tmp.index[0]++;
		// 2nd component
		// p = 0
		index_tmp.index[1]--;
		if (0 <= index_tmp.index[1])
		    basisGradient[ind_index][1] += coeff_D2[ind_index][0] * curlicue_J[1][correspondence.index2number(index_tmp)];
		index_tmp.index[1]++;
		// p = 1;
		index_tmp.index[0]--;
		if (0 <= index_tmp.index[0])
		    basisGradient[ind_index][1] += coeff_D2[ind_index][1] * curlicue_J[1][correspondence.index2number(index_tmp)];
		index_tmp.index[0]++;
		// 3rd component
		// p = 0, q = 0
		index_tmp.index[2]--;
		if (0 <= index_tmp.index[2])
		    basisGradient[ind_index][2] += coeff_D3[ind_index][0] * curlicue_J[2][correspondence.index2number(index_tmp)];
		index_tmp.index[2]++;
		// p = 1, q = 0
		index_tmp.index[0]--; index_tmp.index[1]++; index_tmp.index[2]--;
		if (0 <= index_tmp.index[0] && 0 <= index_tmp.index[2])
		    basisGradient[ind_index][2] += coeff_D3[ind_index][1] * curlicue_J[2][correspondence.index2number(index_tmp)];
		index_tmp.index[0]++; index_tmp.index[1]--; index_tmp.index[2]++;
		// p = 0, q = 1
		index_tmp.index[1]--;
		if (0 <= index_tmp.index[1])
		    basisGradient[ind_index][2] += coeff_D3[ind_index][2] * curlicue_J[2][correspondence.index2number(index_tmp)];
		index_tmp.index[1]++;
		// p = 1, q = 1
		index_tmp.index[0]--;
		if (0 <= index_tmp.index[0])
		    basisGradient[ind_index][2] += coeff_D3[ind_index][3] * curlicue_J[2][correspondence.index2number(index_tmp)];
	    }
	    // assign actual gradient to dep_gradient
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
		for (unsigned int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		    for (unsigned int ind = 0; ind < 3; ++ind)
			dep_gradient[idx_point][transform_local[ind_index][ind_tl]][ind] += weight_transform_local[ind_index][ind_tl] * basisGradient[ind_index][ind];
	}	
    }
}


TEMPLATE_TSEM
void THIS_TSEM::calc_dependence(RegularMesh<3> &mesh, int n_partition, std::vector<std::vector<int> > &idx_eleInCube, std::vector<AFEPack::Point<3> > &point, std::vector<int> &idx_ele, std::vector<std::vector<double> > &dep_basis, std::vector<std::vector<std::vector<double> > > &dep_gradient)
{ // suppose searched point is inside some element, only works for equal partition mesh, whose cube contains 6 elements, while domain is [-1, 1]^3
    int np = point.size();
    for (unsigned int idx_point = 0; idx_point < np; ++idx_point){
	AFEPack::Point<3> &pos = point[idx_point];
	std::vector<int> index(3);
	for (unsigned int idx = 0; idx < 3; ++idx) index[idx] = floor((pos[idx]+1)*0.5*n_partition); // (coord-(-1)) / mesh_size, mesh_size = 2/np
	int idx_cube = (index[0]*n_partition + index[1])*n_partition + index[2];
	for (unsigned int idx_eiq = 0; idx_eiq < 6; ++idx_eiq){ // eiq for element in cube
	    int &ind_e = idx_eleInCube[idx_cube][idx_eiq];
	// for (unsigned int ind_e = 0; ind_e < n_geometry[3]; ++ind_e){
	    int &ind_v0 = index_geometry_onelement[ind_e][0][0], &ind_v1 = index_geometry_onelement[ind_e][0][1], &ind_v2 = index_geometry_onelement[ind_e][0][2], &ind_v3 = index_geometry_onelement[ind_e][0][3];
	    valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
	    valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
	    valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
	    if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) continue;

	    idx_ele[idx_point] = ind_e;

	    // calculate Jacobi polynomial value for basis expression
	    valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
	    valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
	    std::vector<valuetype> basis(n_index);
	    // calculate immediate variable
	    valuetype Jxi_[M+1], Jeta_[M+1], Jzeta_[M+1];
	    Jxi_[0] = Jeta_[0] = Jzeta_[0] = 1; Jxi_[1]  = xi;
	    for (unsigned int l1 = 1; l1 < M; ++l1)
		Jxi_[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi_[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi_[l1-1])	/ calc_coefficient_a(-1, -1, 1, l1);
	    for (unsigned int l1 = 0; l1 <= M; ++l1){
		int aph2 = 2 * l1 - 1;
		Jeta_[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
		for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		    Jeta_[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta_[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta_[l2-1]) / calc_coefficient_a( aph2, -1, 1, l2);
		for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		    int aph3 = 2 * l1 + 2 * l2 - 1;
		    Jzeta_[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
		    for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
			Jzeta_[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta_[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta_[l3-1]) / calc_coefficient_a( aph3, -1, 1, l3);
		    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
			Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
			int ind_index = correspondence.index2number(index_now);
			basis[ind_index] = pow(1-y-z, l1) * Jxi_[l1] * pow(1-z, l2) * Jeta_[l2] * Jzeta_[l3];
		    }
		}
	    }
	    std::vector<valuetype> basis_actual(n_index, 0);
	    for (int ind_index = 0; ind_index < n_index; ++ind_index)
		for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		    basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];

	    for (unsigned int idx_index = 0; idx_index < n_index; ++idx_index)
		dep_basis[idx_point][idx_index] = basis_actual[idx_index];

	    // calculate Jacobi polynomial value for gradient expression
	    valuetype coeff_D1[n_index];
	    valuetype coeff_D2[n_index][2]; // p = 0, 1
	    valuetype coeff_D3[n_index][4]; // (p, q) = (0, 0) (1, 0) (0, 1) (1, 1)
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index){
		Multiindex<3> index_now = correspondence.number2index(ind_index);
		coeff_D1[ind_index] = calc_coefficient_D(0, 0, index_now);
		coeff_D2[ind_index][0] = calc_coefficient_D(3, 0, index_now);
		coeff_D2[ind_index][1] = calc_coefficient_D(3, 1, index_now);
		coeff_D3[ind_index][0] = calc_coefficient_D(5, 0, index_now);
		coeff_D3[ind_index][1] = calc_coefficient_D(5, 1, index_now);
		coeff_D3[ind_index][2] = calc_coefficient_D(5, 2, index_now);
		coeff_D3[ind_index][3] = calc_coefficient_D(5, 3, index_now);
	    }
	    std::vector<std::vector<valuetype> > curlicue_J(3, std::vector<valuetype> (n_index));
	    // assign curlicue J by calculating generalize Jacobi polynomials
	    std::vector<std::vector<double> > basisGradient(n_index, std::vector<double> (3));
	    valuetype x1 = x, x2 = y, x3 = z;
	    valuetype t1 = 1-x2-x3, t2 = 1-x3;
	    xi = 2*x1/t1-1, eta = 2*x2/t2-1, zeta = 2*x3-1;
	    valuetype pow_t1[M+1], pow_t2[M+1];
	    pow_t1[0] = (valuetype) 1; pow_t2[0] = (valuetype) 1;
	    for (unsigned int l = 1; l <= M; ++l){
		pow_t1[l] = pow_t1[l-1] * t1; pow_t2[l] = pow_t2[l-1] * t2; }
	    valuetype Jxi  [2][M+1]; // [0][]: ^{0,0};          [1][]: ^{0,-1}
	    valuetype Jeta [3][M+1]; // [0][]: ^{2l1+1,-1};     [1][]: ^{2l1,0};     [2][]: ^{2l1,-1}
	    valuetype Jzeta[2][M+1]; // [0][]: ^{2l1+2l2+1,-1}; [1][]: ^{2l1+2l2,0}
	    Jxi[0][0] = Jxi[1][0] = (valuetype) 1;
	    Jeta[0][0] = Jeta[1][0] = Jeta[2][0] = (valuetype) 1;
	    Jzeta[0][0] = Jzeta[1][0] = (valuetype) 1;
	    Jxi[0][1] = calc_generalized_jacobi_polynomial(0,  0, 1, xi);
	    Jxi[1][1] = calc_generalized_jacobi_polynomial(0, -1, 1, xi);
	    for (unsigned int l = 1; l < M; ++l){
		Jxi[0][l+1] = ((xi - calc_coefficient_a(0, 0,2,l)) * Jxi[0][l] - calc_coefficient_a(0, 0,3,l) * Jxi[0][l-1]) / calc_coefficient_a(0, 0,1,l);
		Jxi[1][l+1] = ((xi - calc_coefficient_a(0,-1,2,l)) * Jxi[1][l] - calc_coefficient_a(0,-1,3,l) * Jxi[1][l-1]) / calc_coefficient_a(0,-1,1,l);
	    }
	    for (unsigned int l1 = 0; l1 <= M; ++l1){
		int aph2_0 = 2*l1+1, aph2_1 = 2*l1;
		Jeta[0][1] = calc_generalized_jacobi_polynomial(aph2_0, -1, 1, eta);
		Jeta[1][1] = calc_generalized_jacobi_polynomial(aph2_1,  0, 1, eta);
		Jeta[2][1] = calc_generalized_jacobi_polynomial(aph2_1, -1, 1, eta);
		for (unsigned int l = 1; l < M-l1; ++l){
		    Jeta[0][l+1] = ((eta - calc_coefficient_a(aph2_0,-1,2,l)) * Jeta[0][l] - calc_coefficient_a(aph2_0,-1,3,l) * Jeta[0][l-1]) / calc_coefficient_a(aph2_0,-1,1,l);
		    Jeta[1][l+1] = ((eta - calc_coefficient_a(aph2_1, 0,2,l)) * Jeta[1][l] - calc_coefficient_a(aph2_1, 0,3,l) * Jeta[1][l-1]) / calc_coefficient_a(aph2_1, 0,1,l);
		    Jeta[2][l+1] = ((eta - calc_coefficient_a(aph2_1,-1,2,l)) * Jeta[2][l] - calc_coefficient_a(aph2_1,-1,3,l) * Jeta[2][l-1]) / calc_coefficient_a(aph2_1,-1,1,l);
		}
		for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		    int aph3_0 = 2*l1+2*l2+1, aph3_1 = 2*l1+2*l2;
		    Jzeta[0][1] = calc_generalized_jacobi_polynomial(aph3_0, -1, 1, zeta);
		    Jzeta[1][1] = calc_generalized_jacobi_polynomial(aph3_1,  0, 1, zeta);
		    for (unsigned int l = 1; l < M-l1-l2; ++l){
			Jzeta[0][l+1] = ((zeta-calc_coefficient_a(aph3_0,-1,2,l)) * Jzeta[0][l] - calc_coefficient_a(aph3_0,-1,3,l) * Jzeta[0][l-1]) / calc_coefficient_a(aph3_0,-1,1,l);
			Jzeta[1][l+1] = ((zeta-calc_coefficient_a(aph3_1, 0,2,l)) * Jzeta[1][l] - calc_coefficient_a(aph3_1, 0,3,l) * Jzeta[1][l-1]) / calc_coefficient_a(aph3_1, 0,1,l);
		    }
		    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
			Multiindex<3> index_now = Unitary_Multiindex[0]*l1 + Unitary_Multiindex[1]*l2 + Unitary_Multiindex[2]*l3;
			int ind_index = correspondence.index2number(index_now);
			curlicue_J[0][ind_index] = pow_t1[l1]*Jxi[0][l1] * pow_t2[l2]*Jeta[0][l2] * Jzeta[0][l3];
			curlicue_J[1][ind_index] = pow_t1[l1]*Jxi[1][l1] * pow_t2[l2]*Jeta[1][l2] * Jzeta[0][l3];
			curlicue_J[2][ind_index] = pow_t1[l1]*Jxi[1][l1] * pow_t2[l2]*Jeta[2][l2] * Jzeta[1][l3];
		    }
		}
	    }
	    // assign basis_gradient[p][][]
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index){
		Multiindex<3> index_tmp = correspondence.number2index(ind_index);
		// 1st component
		index_tmp.index[0]--;
		if (0 <= index_tmp.index[0])
		    basisGradient[ind_index][0] += coeff_D1[ind_index]    * curlicue_J[0][correspondence.index2number(index_tmp)];
		index_tmp.index[0]++;
		// 2nd component
		// p = 0
		index_tmp.index[1]--;
		if (0 <= index_tmp.index[1])
		    basisGradient[ind_index][1] += coeff_D2[ind_index][0] * curlicue_J[1][correspondence.index2number(index_tmp)];
		index_tmp.index[1]++;
		// p = 1;
		index_tmp.index[0]--;
		if (0 <= index_tmp.index[0])
		    basisGradient[ind_index][1] += coeff_D2[ind_index][1] * curlicue_J[1][correspondence.index2number(index_tmp)];
		index_tmp.index[0]++;
		// 3rd component
		// p = 0, q = 0
		index_tmp.index[2]--;
		if (0 <= index_tmp.index[2])
		    basisGradient[ind_index][2] += coeff_D3[ind_index][0] * curlicue_J[2][correspondence.index2number(index_tmp)];
		index_tmp.index[2]++;
		// p = 1, q = 0
		index_tmp.index[0]--; index_tmp.index[1]++; index_tmp.index[2]--;
		if (0 <= index_tmp.index[0] && 0 <= index_tmp.index[2])
		    basisGradient[ind_index][2] += coeff_D3[ind_index][1] * curlicue_J[2][correspondence.index2number(index_tmp)];
		index_tmp.index[0]++; index_tmp.index[1]--; index_tmp.index[2]++;
		// p = 0, q = 1
		index_tmp.index[1]--;
		if (0 <= index_tmp.index[1])
		    basisGradient[ind_index][2] += coeff_D3[ind_index][2] * curlicue_J[2][correspondence.index2number(index_tmp)];
		index_tmp.index[1]++;
		// p = 1, q = 1
		index_tmp.index[0]--;
		if (0 <= index_tmp.index[0])
		    basisGradient[ind_index][2] += coeff_D3[ind_index][3] * curlicue_J[2][correspondence.index2number(index_tmp)];
	    }
	    // assign actual gradient to dep_gradient
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
		for (unsigned int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		    for (unsigned int ind = 0; ind < 3; ++ind)
			dep_gradient[idx_point][transform_local[ind_index][ind_tl]][ind] += weight_transform_local[ind_index][ind_tl] * basisGradient[ind_index][ind];
	}	
    }
}

TEMPLATE_TSEM
void THIS_TSEM::calc_val_gradient_qp_onElement(Vector<valuetype> &src, std::vector<std::vector<double> > &dst, RegularMesh<3> &mesh, unsigned int idx_ele)
{
    int &n_qp = n_q_point[2];
    std::vector<valuetype> val_coef(n_index);
    calc_coef_onElement(src, val_coef, idx_ele);
    for (unsigned int idx_qp = 0; idx_qp < n_qp; ++idx_qp)
	for (unsigned int idx = 0; idx < 3; ++idx)
	    dst[idx_qp][idx] = 0;
    // prepare coordinates of element vertices
    valuetype xj_k[3][4];
    for (unsigned int k = 0; k < 4; ++k){
	int &ind_node = index_geometry_onelement[idx_ele][0][k];
	for (unsigned int j = 0; j < 3; ++j) xj_k[j][k] = mesh.point(ind_node)[j];
    }
    valuetype d_jk[3][3];
    for (unsigned int j = 0; j < 3; ++j)
	for (unsigned int k = 0; k < 3; ++k)
	    d_jk[j][k] = xj_k[j][k+1] - xj_k[j][0];
    valuetype val_V = d_jk[0][0]*d_jk[1][1]*d_jk[2][2] + d_jk[0][1]*d_jk[1][2]*d_jk[2][0] + d_jk[0][2]*d_jk[1][0]*d_jk[2][1]
	- d_jk[0][0]*d_jk[1][2]*d_jk[2][1] - d_jk[0][1]*d_jk[1][0]*d_jk[2][2] - d_jk[0][2]*d_jk[1][1]*d_jk[2][0];
    // calculate \partial x_j/\partial \hat{x}_k
    valuetype parxj_parxk[3][3];
    parxj_parxk[0][0] = (d_jk[1][1]*d_jk[2][2] - d_jk[1][2]*d_jk[2][1]) / val_V;
    parxj_parxk[0][1] = (d_jk[0][2]*d_jk[2][1] - d_jk[0][1]*d_jk[2][2]) / val_V;
    parxj_parxk[0][2] = (d_jk[0][1]*d_jk[1][2] - d_jk[0][2]*d_jk[1][1]) / val_V;
    parxj_parxk[1][0] = (d_jk[1][2]*d_jk[2][0] - d_jk[1][0]*d_jk[2][2]) / val_V;
    parxj_parxk[1][1] = (d_jk[0][0]*d_jk[2][2] - d_jk[0][2]*d_jk[2][0]) / val_V;
    parxj_parxk[1][2] = (d_jk[0][2]*d_jk[1][0] - d_jk[0][0]*d_jk[1][2]) / val_V;
    parxj_parxk[2][0] = (d_jk[1][0]*d_jk[2][1] - d_jk[1][1]*d_jk[2][0]) / val_V;
    parxj_parxk[2][1] = (d_jk[0][1]*d_jk[2][0] - d_jk[0][0]*d_jk[2][1]) / val_V;
    parxj_parxk[2][2] = (d_jk[0][0]*d_jk[1][1] - d_jk[0][1]*d_jk[1][0]) / val_V;
    // traverse quadrature point, assign corresponding value
    for (unsigned int p = 0; p < n_qp; ++p)
	// recover the value of basis gradient
	for (int ind = 0; ind < 3; ++ind)
	    for (int ind_index = 0; ind_index < n_index; ++ind_index)
		for (int indt = 0; indt < 3; ++indt)
		    dst[p][ind] += parxj_parxk[indt][ind] * val_coef[ind_index] * basis_gradient_actual[p][ind_index][indt];
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_inElement(RegularMesh<3> &mesh, int ind_ele, AFEPack::Point<3> &pos, Vector<valuetype> &psi)
{ // calculate the value of point with coordinate pos on ind_ele-th element
    int n_orbital = psi.size();
    // vertex
    for (unsigned int ind_p = 0; ind_p < 4; ++ind_p){
	int &ind_point = index_geometry_onelement[ind_ele][0][ind_p];
	if (calc_length_line(pos, mesh.point(ind_point)) < tol_zero) return psi(location_actualdof[location_geometry[ind_point]]);
    }
    // edge
    for (unsigned int ind_e = 0; ind_e < 6; ++ind_e){
	int &ind_edge = index_geometry_onelement[ind_ele][1][ind_e];
	int &ind_point0 = number_node[0][ind_edge][0];
	int &ind_point1 = number_node[0][ind_edge][1];
	valuetype dis = calc_length_line(mesh.point(ind_point0), mesh.point(ind_point1));
	valuetype dis0 = calc_length_line(pos, mesh.point(ind_point0));
	valuetype dis1 = calc_length_line(pos, mesh.point(ind_point1));
	if (fabs(dis0+dis1 - dis) > tol_zero) continue;
	int &loc0 = location_actualdof[location_geometry[ind_point0]];
	int &loc1 = location_actualdof[location_geometry[ind_point1]];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + ind_edge]];
	
	valuetype x = dis0 / dis, r = 1 - x;
	valuetype xi = 2 * x - 1;
	valuetype Jxi[M+1];
	Jxi[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);

	valuetype psi_local = psi(loc0)*r + psi(loc1)*x;
	for (unsigned int l1 = 2; l1 <= M; ++l1)
	    psi_local += 2 * Jxi[l1] * psi(loc + l1-2);
	return psi_local;
    }
    // face
    for (unsigned int ind_f = 0; ind_f < 4; ++ind_f){
	int &ind_face = index_geometry_onelement[ind_ele][2][ind_f];
	int &ind_point0 = number_node[1][ind_face][0];
	int &ind_point1 = number_node[1][ind_face][1];
	int &ind_point2 = number_node[1][ind_face][2];
	valuetype area = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area0 = calc_area_triangle(pos, mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area1 = calc_area_triangle(mesh.point(ind_point0), pos, mesh.point(ind_point2));
	valuetype area2 = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), pos);
	if (fabs(area0+area1+area2 - area) > tol_zero) continue;
	int &ind_edge0 = number_edge[ind_face][0];
	int &ind_edge1 = number_edge[ind_face][1];
	int &ind_edge2 = number_edge[ind_face][2];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + n_geometry[1] + ind_face]];
	int &locp0 = location_actualdof[location_geometry[ind_point0]];
	int &locp1 = location_actualdof[location_geometry[ind_point1]];
	int &locp2 = location_actualdof[location_geometry[ind_point2]];
	int &loce0 = location_actualdof[location_geometry[n_geometry[0] + ind_edge0]];
	int &loce1 = location_actualdof[location_geometry[n_geometry[0] + ind_edge1]];
	int &loce2 = location_actualdof[location_geometry[n_geometry[0] + ind_edge2]];
	
	std::vector<valuetype> bas_val(n_dof_geometry[2], 0);
	valuetype x = area1 / area, y = area2 / area, r = 1 - x - y;
	valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
	valuetype Jxi[M+1], Jeta[M+1];
	Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 2; l1 <= M-1; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
	    }
	}
	std::vector<valuetype> bas_val_add_xi( n_dof_geometry[1], 0);
	std::vector<valuetype> bas_val_add_eta(n_dof_geometry[1], 0);
	Jxi[0] = Jeta[0] = 1;
	Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
	Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
	for (unsigned int l = 1; l < M; ++l){
	    Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
	    Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
	}
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
	    bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
	    bas_val_add_eta[l] = Jeta[l];
	}

	valuetype psi_local = psi(locp0) * r + psi(locp1) * x + psi(locp2) * y;
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l)
	    psi_local -= 2 * (x*y * psi(loce0+l) * bas_val_add_eta[l] +
			      y*r * psi(loce1+l) * bas_val_add_eta[l] +
			      x*r * psi(loce2+l) * bas_val_add_xi[l]);
	for (unsigned int l1 = 2; l1 <= M; ++l1)
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		psi_local += bas_val[ind_index] * psi(loc + ind_index);
	    }
	return psi_local;
    }
    // interior
    int &ind_v0 = index_geometry_onelement[ind_ele][0][0];
    int &ind_v1 = index_geometry_onelement[ind_ele][0][1];
    int &ind_v2 = index_geometry_onelement[ind_ele][0][2];
    int &ind_v3 = index_geometry_onelement[ind_ele][0][3];
    valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
    valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
    if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero)
	std::cout << "error: interpolation from SEM solution to FEM solution\n";

    valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
    valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
    std::vector<valuetype> basis(n_index);
    // calculate immediate variable
    valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
    Jxi[0] = Jeta[0] = Jzeta[0] = 1;
    Jxi[1]  = xi;
    for (unsigned int l1 = 1; l1 < M; ++l1)
	Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
	    / calc_coefficient_a(-1, -1, 1, l1);
    for (unsigned int l1 = 0; l1 <= M; ++l1){
	int aph2 = 2 * l1 - 1;
	Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	for (unsigned int l2 = 1; l2 < M-l1; ++l2)
	    Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		/ calc_coefficient_a( aph2, -1, 1, l2);
	for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
	    int aph3 = 2 * l1 + 2 * l2 - 1;
	    Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
	    for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
		    / calc_coefficient_a( aph3, -1, 1, l3);
	    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		int ind_index = correspondence.index2number(index_now);
		basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
	    }
	}
    }
    std::vector<valuetype> basis_actual(n_index, 0);
    for (int ind_index = 0; ind_index < n_index; ++ind_index)
	for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
	    basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];

    std::vector<valuetype> coef_local(n_index, 0);
    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
	    coef_local[ind_index] += psi(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
		* transform_val_global2local[ind_ele][ind_index][ind_nnz];
    valuetype psi_local = 0;
    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	psi_local += coef_local[ind_index] * basis_actual[ind_index];
    return psi_local;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_inElement(RegularMesh<3> &mesh, int ind_ele, AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi)
{ // calculate the value of point with coordinate pos on ind_ele-th element for closed shell system
    int n_orbital = psi.size();
    // vertex
    for (unsigned int ind_p = 0; ind_p < 4; ++ind_p){
	int &ind_point = index_geometry_onelement[ind_ele][0][ind_p];
	if (calc_length_line(pos, mesh.point(ind_point)) < tol_zero){
	    valuetype rho_local = 0;
	    for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
		valuetype n_occupation = 2.0;
		int &loc = location_actualdof[location_geometry[ind_point]];
		rho_local += n_occupation * pow(psi[ind_orbital](loc), 2);
	    }
	    return rho_local;
	}
    }
    // edge
    for (unsigned int ind_e = 0; ind_e < 6; ++ind_e){
	int &ind_edge = index_geometry_onelement[ind_ele][1][ind_e];
	int &ind_point0 = number_node[0][ind_edge][0];
	int &ind_point1 = number_node[0][ind_edge][1];
	valuetype dis = calc_length_line(mesh.point(ind_point0), mesh.point(ind_point1));
	valuetype dis0 = calc_length_line(pos, mesh.point(ind_point0));
	valuetype dis1 = calc_length_line(pos, mesh.point(ind_point1));
	if (fabs(dis0+dis1 - dis) > tol_zero) continue;
	int &loc0 = location_actualdof[location_geometry[ind_point0]];
	int &loc1 = location_actualdof[location_geometry[ind_point1]];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + ind_edge]];
	
	valuetype x = dis0 / dis, r = 1 - x;
	valuetype xi = 2 * x - 1;
	valuetype Jxi[M+1];
	Jxi[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	
	valuetype rho_local = 0;
	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    valuetype n_occupation = 2.0;
	    valuetype psi_local = psi[ind_orbital](loc0)*r + psi[ind_orbital](loc1)*x;
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		psi_local += 2 * Jxi[l1] * psi[ind_orbital](loc + l1-2);
	    rho_local += n_occupation * pow(psi_local, 2);
	}
	return rho_local;
    }
    // face
    for (unsigned int ind_f = 0; ind_f < 4; ++ind_f){
	int &ind_face = index_geometry_onelement[ind_ele][2][ind_f];
	int &ind_point0 = number_node[1][ind_face][0];
	int &ind_point1 = number_node[1][ind_face][1];
	int &ind_point2 = number_node[1][ind_face][2];
	valuetype area = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area0 = calc_area_triangle(pos, mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area1 = calc_area_triangle(mesh.point(ind_point0), pos, mesh.point(ind_point2));
	valuetype area2 = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), pos);
	if (fabs(area0+area1+area2 - area) > tol_zero) continue;
	int &ind_edge0 = number_edge[ind_face][0];
	int &ind_edge1 = number_edge[ind_face][1];
	int &ind_edge2 = number_edge[ind_face][2];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + n_geometry[1] + ind_face]];
	int &locp0 = location_actualdof[location_geometry[ind_point0]];
	int &locp1 = location_actualdof[location_geometry[ind_point1]];
	int &locp2 = location_actualdof[location_geometry[ind_point2]];
	int &loce0 = location_actualdof[location_geometry[n_geometry[0] + ind_edge0]];
	int &loce1 = location_actualdof[location_geometry[n_geometry[0] + ind_edge1]];
	int &loce2 = location_actualdof[location_geometry[n_geometry[0] + ind_edge2]];
	
	std::vector<valuetype> bas_val(n_dof_geometry[2], 0);
	valuetype x = area1 / area, y = area2 / area, r = 1 - x - y;
	valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
	valuetype Jxi[M+1], Jeta[M+1];
	Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 2; l1 <= M-1; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
	    }
	}
	std::vector<valuetype> bas_val_add_xi( n_dof_geometry[1], 0);
	std::vector<valuetype> bas_val_add_eta(n_dof_geometry[1], 0);
	Jxi[0] = Jeta[0] = 1;
	Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
	Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
	for (unsigned int l = 1; l < M; ++l){
	    Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
	    Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
	}
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
	    bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
	    bas_val_add_eta[l] = Jeta[l];
	}
	
	valuetype rho_local = 0;
	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    valuetype n_occupation = 2.0;
	    valuetype psi_local;
	    psi_local = psi[ind_orbital](locp0) * r + psi[ind_orbital](locp1) * x + psi[ind_orbital](locp2) * y;
	    for (unsigned int l = 0; l < n_dof_geometry[1]; ++l)
		psi_local -= 2 * (x*y * psi[ind_orbital](loce0+l) * bas_val_add_eta[l] +
				  y*r * psi[ind_orbital](loce1+l) * bas_val_add_eta[l] +
				  x*r * psi[ind_orbital](loce2+l) * bas_val_add_xi[l]);
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		    unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		    psi_local += bas_val[ind_index] * psi[ind_orbital](loc + ind_index);
		}
	    rho_local += n_occupation * pow(psi_local, 2);
	}
	return rho_local;
    }
    // interior
    int &ind_v0 = index_geometry_onelement[ind_ele][0][0];
    int &ind_v1 = index_geometry_onelement[ind_ele][0][1];
    int &ind_v2 = index_geometry_onelement[ind_ele][0][2];
    int &ind_v3 = index_geometry_onelement[ind_ele][0][3];
    valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
    valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
    if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero)
	std::cout << "error: interpolation from SEM solution to FEM solution\n";

    valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
    valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
    std::vector<valuetype> basis(n_index);
    // calculate immediate variable
    valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
    Jxi[0] = Jeta[0] = Jzeta[0] = 1;
    Jxi[1]  = xi;
    for (unsigned int l1 = 1; l1 < M; ++l1)
	Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
	    / calc_coefficient_a(-1, -1, 1, l1);
    for (unsigned int l1 = 0; l1 <= M; ++l1){
	int aph2 = 2 * l1 - 1;
	Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	for (unsigned int l2 = 1; l2 < M-l1; ++l2)
	    Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		/ calc_coefficient_a( aph2, -1, 1, l2);
	for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
	    int aph3 = 2 * l1 + 2 * l2 - 1;
	    Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
	    for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
		    / calc_coefficient_a( aph3, -1, 1, l3);
	    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		int ind_index = correspondence.index2number(index_now);
		basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
	    }
	}
    }
    std::vector<valuetype> basis_actual(n_index, 0);
    for (int ind_index = 0; ind_index < n_index; ++ind_index)
	for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
	    basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];
    
    valuetype rho_local = 0;
    for (int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	valuetype n_occupation = 2.0;
	// calculate function value
	std::vector<valuetype> coef_local(n_index, 0);
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
		coef_local[ind_index] += psi[ind_orbital](transform_ind_global2local[ind_ele][ind_index][ind_nnz])
		    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	valuetype psi_local = 0;
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    psi_local += coef_local[ind_index] * basis_actual[ind_index];
	// add contribution
	rho_local += n_occupation * pow(psi_local, 2);
    }
    return rho_local;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_inElement(RegularMesh<3> &mesh, int ind_ele, AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi, std::vector<valuetype>& n_occupation)
{ // calculate the value of point with coordinate pos on ind_ele-th element
    int n_orbital = psi.size();
    // vertex
    for (unsigned int ind_p = 0; ind_p < 4; ++ind_p){
	int &ind_point = index_geometry_onelement[ind_ele][0][ind_p];
	if (calc_length_line(pos, mesh.point(ind_point)) < tol_zero){
	    valuetype rho_local = 0;
	    for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
		int &loc = location_actualdof[location_geometry[ind_point]];
		rho_local += n_occupation[ind_orbital] * pow(psi[ind_orbital](loc), 2);
	    }
	    return rho_local;
	}
    }
    // edge
    for (unsigned int ind_e = 0; ind_e < 6; ++ind_e){
	int &ind_edge = index_geometry_onelement[ind_ele][1][ind_e];
	int &ind_point0 = number_node[0][ind_edge][0];
	int &ind_point1 = number_node[0][ind_edge][1];
	valuetype dis = calc_length_line(mesh.point(ind_point0), mesh.point(ind_point1));
	valuetype dis0 = calc_length_line(pos, mesh.point(ind_point0));
	valuetype dis1 = calc_length_line(pos, mesh.point(ind_point1));
	if (fabs(dis0+dis1 - dis) > tol_zero) continue;
	int &loc0 = location_actualdof[location_geometry[ind_point0]];
	int &loc1 = location_actualdof[location_geometry[ind_point1]];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + ind_edge]];
	
	valuetype x = dis0 / dis, r = 1 - x;
	valuetype xi = 2 * x - 1;
	valuetype Jxi[M+1];
	Jxi[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	
	valuetype rho_local = 0;
	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    valuetype psi_local = psi[ind_orbital](loc0)*r + psi[ind_orbital](loc1)*x;
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		psi_local += 2 * Jxi[l1] * psi[ind_orbital](loc + l1-2);
	    rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
	}
	return rho_local;
    }
    // face
    for (unsigned int ind_f = 0; ind_f < 4; ++ind_f){
	int &ind_face = index_geometry_onelement[ind_ele][2][ind_f];
	int &ind_point0 = number_node[1][ind_face][0];
	int &ind_point1 = number_node[1][ind_face][1];
	int &ind_point2 = number_node[1][ind_face][2];
	valuetype area = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area0 = calc_area_triangle(pos, mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area1 = calc_area_triangle(mesh.point(ind_point0), pos, mesh.point(ind_point2));
	valuetype area2 = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), pos);
	if (fabs(area0+area1+area2 - area) > tol_zero) continue;
	int &ind_edge0 = number_edge[ind_face][0];
	int &ind_edge1 = number_edge[ind_face][1];
	int &ind_edge2 = number_edge[ind_face][2];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + n_geometry[1] + ind_face]];
	int &locp0 = location_actualdof[location_geometry[ind_point0]];
	int &locp1 = location_actualdof[location_geometry[ind_point1]];
	int &locp2 = location_actualdof[location_geometry[ind_point2]];
	int &loce0 = location_actualdof[location_geometry[n_geometry[0] + ind_edge0]];
	int &loce1 = location_actualdof[location_geometry[n_geometry[0] + ind_edge1]];
	int &loce2 = location_actualdof[location_geometry[n_geometry[0] + ind_edge2]];
	
	std::vector<valuetype> bas_val(n_dof_geometry[2], 0);
	valuetype x = area1 / area, y = area2 / area, r = 1 - x - y;
	valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
	valuetype Jxi[M+1], Jeta[M+1];
	Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 2; l1 <= M-1; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
	    }
	}
	std::vector<valuetype> bas_val_add_xi( n_dof_geometry[1], 0);
	std::vector<valuetype> bas_val_add_eta(n_dof_geometry[1], 0);
	Jxi[0] = Jeta[0] = 1;
	Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
	Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
	for (unsigned int l = 1; l < M; ++l){
	    Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
	    Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
	}
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
	    bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
	    bas_val_add_eta[l] = Jeta[l];
	}
	
	valuetype rho_local = 0;
	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    valuetype psi_local;
	    psi_local = psi[ind_orbital](locp0) * r + psi[ind_orbital](locp1) * x + psi[ind_orbital](locp2) * y;
	    for (unsigned int l = 0; l < n_dof_geometry[1]; ++l)
		psi_local -= 2 * (x*y * psi[ind_orbital](loce0+l) * bas_val_add_eta[l] +
				  y*r * psi[ind_orbital](loce1+l) * bas_val_add_eta[l] +
				  x*r * psi[ind_orbital](loce2+l) * bas_val_add_xi[l]);
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		    unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		    psi_local += bas_val[ind_index] * psi[ind_orbital](loc + ind_index);
		}
	    rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
	}
	return rho_local;
    }
    // interior
    int &ind_v0 = index_geometry_onelement[ind_ele][0][0];
    int &ind_v1 = index_geometry_onelement[ind_ele][0][1];
    int &ind_v2 = index_geometry_onelement[ind_ele][0][2];
    int &ind_v3 = index_geometry_onelement[ind_ele][0][3];
    valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
    valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
    if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero)
	std::cout << "error: interpolation from SEM solution to FEM solution\n";

    valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
    valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
    std::vector<valuetype> basis(n_index);
    // calculate immediate variable
    valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
    Jxi[0] = Jeta[0] = Jzeta[0] = 1;
    Jxi[1]  = xi;
    for (unsigned int l1 = 1; l1 < M; ++l1)
	Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
	    / calc_coefficient_a(-1, -1, 1, l1);
    for (unsigned int l1 = 0; l1 <= M; ++l1){
	int aph2 = 2 * l1 - 1;
	Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	for (unsigned int l2 = 1; l2 < M-l1; ++l2)
	    Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		/ calc_coefficient_a( aph2, -1, 1, l2);
	for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
	    int aph3 = 2 * l1 + 2 * l2 - 1;
	    Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
	    for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
		    / calc_coefficient_a( aph3, -1, 1, l3);
	    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		int ind_index = correspondence.index2number(index_now);
		basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
	    }
	}
    }
    std::vector<valuetype> basis_actual(n_index, 0);
    for (int ind_index = 0; ind_index < n_index; ++ind_index)
	for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
	    basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];
    
    valuetype rho_local = 0;
    for (int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	// calculate function value
	std::vector<valuetype> coef_local(n_index, 0);
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
		coef_local[ind_index] += psi[ind_orbital](transform_ind_global2local[ind_ele][ind_index][ind_nnz])
		    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	valuetype psi_local = 0;
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    psi_local += coef_local[ind_index] * basis_actual[ind_index];
	// add contribution
	rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
    }
    return rho_local;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_point(RegularMesh<3> &mesh, AFEPack::Point<3> &pos, Vector<valuetype> &u)
{
    // traverse 0-d geometries
    for (unsigned int ind_p = 0; ind_p < n_geometry[0]; ++ind_p)
	if (calc_length_line(mesh.point(ind_p), pos) < tol_zero){
	    int &loc = location_actualdof[location_geometry[ind_p]];
	    valuetype u_local = u(loc);
	    return u_local;
	}
    // std::cout << "searched point\n";

    // traverse 1-d geometries
    for (unsigned int ind_e = 0; ind_e < n_geometry[1]; ++ind_e){
	int &ind_point0 = number_node[0][ind_e][0];
	int &ind_point1 = number_node[0][ind_e][1];
	valuetype dis = calc_length_line(mesh.point(ind_point0), mesh.point(ind_point1));
	valuetype dis0 = calc_length_line(mesh.point(ind_point0), pos);
	valuetype dis1 = calc_length_line(mesh.point(ind_point1), pos);
	if (fabs(dis0+dis1 - dis) > tol_zero) continue;
	int &loc0 = location_actualdof[location_geometry[ind_point0]];
	int &loc1 = location_actualdof[location_geometry[ind_point1]];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + ind_e]];
	
	valuetype x = dis0 / dis, r = 1 - x;
	valuetype xi = 2 * x - 1;
	valuetype Jxi[M+1];
	Jxi[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	
	valuetype u_local = u(loc0)*r + u(loc1)*x;
	for (unsigned int l1 = 2; l1 <= M; ++l1) u_local += 2 * Jxi[l1] * u(loc+l1-2);
	return u_local;
    }
    // std::cout << "searched edge\n";
    
    // traverse 2-d geometries
    for (unsigned int ind_f = 0; ind_f < n_geometry[2]; ++ind_f){
	int &ind_point0 = number_node[1][ind_f][0];
	int &ind_point1 = number_node[1][ind_f][1];
	int &ind_point2 = number_node[1][ind_f][2];
	valuetype area = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area0 = calc_area_triangle(pos, mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area1 = calc_area_triangle(mesh.point(ind_point0), pos, mesh.point(ind_point2));
	valuetype area2 = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), pos);
	if (fabs(area0+area1+area2 - area) > tol_zero) continue;
	int &ind_edge0 = number_edge[ind_f][0];
	int &ind_edge1 = number_edge[ind_f][1];
	int &ind_edge2 = number_edge[ind_f][2];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + n_geometry[1] + ind_f]];
	int &locp0 = location_actualdof[location_geometry[ind_point0]];
	int &locp1 = location_actualdof[location_geometry[ind_point1]];
	int &locp2 = location_actualdof[location_geometry[ind_point2]];
	int &loce0 = location_actualdof[location_geometry[n_geometry[0] + ind_edge0]];
	int &loce1 = location_actualdof[location_geometry[n_geometry[0] + ind_edge1]];
	int &loce2 = location_actualdof[location_geometry[n_geometry[0] + ind_edge2]];
	
	std::vector<valuetype> bas_val(n_dof_geometry[2], 0);
	valuetype x = area1 / area, y = area2 / area, r = 1 - x - y;
	valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
	valuetype Jxi[M+1], Jeta[M+1];
	Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 2; l1 <= M-1; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
	    }
	}
	std::vector<valuetype> bas_val_add_xi( n_dof_geometry[1], 0);
	std::vector<valuetype> bas_val_add_eta(n_dof_geometry[1], 0);
	Jxi[0] = Jeta[0] = 1;
	Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
	Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
	for (unsigned int l = 1; l < M; ++l){
	    Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
	    Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
	}
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
	    bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
	    bas_val_add_eta[l] = Jeta[l];
	}
	
	valuetype u_local = u(locp0)*r + u(locp1)*x + u(locp2)*y;
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l)
	    u_local -= 2 * (x*y * u(loce0+l) * bas_val_add_eta[l] +
			    y*r * u(loce1+l) * bas_val_add_eta[l] +
			    x*r * u(loce2+l) * bas_val_add_xi[l]);
	for (unsigned int l1 = 2; l1 <= M; ++l1)
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		u_local += bas_val[ind_index] * u(loc + ind_index);
	    }
	return u_local;
    }
    // std::cout << "searched face\n";
    
    // traverse 3-d geometries
    for (unsigned int ind_e = 0; ind_e < n_geometry[3]; ++ind_e){
	int &ind_v0 = index_geometry_onelement[ind_e][0][0];
	int &ind_v1 = index_geometry_onelement[ind_e][0][1];
	int &ind_v2 = index_geometry_onelement[ind_e][0][2];
	int &ind_v3 = index_geometry_onelement[ind_e][0][3];
	valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
	valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
	if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) continue;
    
	// valuetype u_local = calc_val_point_insideElement(correspondence, mesh, ind_e, pos, u);
	valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
	valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
	std::vector<valuetype> basis(n_index);
	// calculate immediate variable
	valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
	Jxi[0] = Jeta[0] = Jzeta[0] = 1;
	Jxi[1]  = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 0; l1 <= M; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		int aph3 = 2 * l1 + 2 * l2 - 1;
		Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
		for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		    Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
			/ calc_coefficient_a( aph3, -1, 1, l3);
		for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		    Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		    int ind_index = correspondence.index2number(index_now);
		    basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
		}
	    }
	}
	std::vector<valuetype> basis_actual(n_index, 0);
	for (int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];
    
	valuetype u_local = 0;
	// calculate function value
	std::vector<valuetype> coef_local(n_index, 0);
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_e][ind_index]; ++ind_nnz)
		coef_local[ind_index] += u(transform_ind_global2local[ind_e][ind_index][ind_nnz])
		    * transform_val_global2local[ind_e][ind_index][ind_nnz];

	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    u_local += coef_local[ind_index] * basis_actual[ind_index];
	return u_local;
    }

    std::cout << "error, cant find position of given pos in current SEM mesh, pos = (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    return 0;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_point(Mesh<3> &mesh, AFEPack::Point<3> &pos, Vector<valuetype> &u)
{
    // traverse 0-d geometries
    for (unsigned int ind_p = 0; ind_p < n_geometry[0]; ++ind_p)
	if (calc_length_line(mesh.point(ind_p), pos) < tol_zero){
	    int &loc = location_actualdof[location_geometry[ind_p]];
	    valuetype u_local = u(loc);
	    return u_local;
	}
    // std::cout << "searched point\n";

    // traverse 1-d geometries
    for (unsigned int ind_e = 0; ind_e < n_geometry[1]; ++ind_e){
	int &ind_point0 = number_node[0][ind_e][0];
	int &ind_point1 = number_node[0][ind_e][1];
	valuetype dis = calc_length_line(mesh.point(ind_point0), mesh.point(ind_point1));
	valuetype dis0 = calc_length_line(mesh.point(ind_point0), pos);
	valuetype dis1 = calc_length_line(mesh.point(ind_point1), pos);
	if (fabs(dis0+dis1 - dis) > tol_zero) continue;
	int &loc0 = location_actualdof[location_geometry[ind_point0]];
	int &loc1 = location_actualdof[location_geometry[ind_point1]];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + ind_e]];
	
	valuetype x = dis0 / dis, r = 1 - x;
	valuetype xi = 2 * x - 1;
	valuetype Jxi[M+1];
	Jxi[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	
	valuetype u_local = u(loc0)*r + u(loc1)*x;
	for (unsigned int l1 = 2; l1 <= M; ++l1) u_local += 2 * Jxi[l1] * u(loc+l1-2);
	return u_local;
    }
    // std::cout << "searched edge\n";
    
    // traverse 2-d geometries
    for (unsigned int ind_f = 0; ind_f < n_geometry[2]; ++ind_f){
	int &ind_point0 = number_node[1][ind_f][0];
	int &ind_point1 = number_node[1][ind_f][1];
	int &ind_point2 = number_node[1][ind_f][2];
	valuetype area = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area0 = calc_area_triangle(pos, mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area1 = calc_area_triangle(mesh.point(ind_point0), pos, mesh.point(ind_point2));
	valuetype area2 = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), pos);
	if (fabs(area0+area1+area2 - area) > tol_zero) continue;
	int &ind_edge0 = number_edge[ind_f][0];
	int &ind_edge1 = number_edge[ind_f][1];
	int &ind_edge2 = number_edge[ind_f][2];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + n_geometry[1] + ind_f]];
	int &locp0 = location_actualdof[location_geometry[ind_point0]];
	int &locp1 = location_actualdof[location_geometry[ind_point1]];
	int &locp2 = location_actualdof[location_geometry[ind_point2]];
	int &loce0 = location_actualdof[location_geometry[n_geometry[0] + ind_edge0]];
	int &loce1 = location_actualdof[location_geometry[n_geometry[0] + ind_edge1]];
	int &loce2 = location_actualdof[location_geometry[n_geometry[0] + ind_edge2]];
	
	std::vector<valuetype> bas_val(n_dof_geometry[2], 0);
	valuetype x = area1 / area, y = area2 / area, r = 1 - x - y;
	valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
	valuetype Jxi[M+1], Jeta[M+1];
	Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 2; l1 <= M-1; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
	    }
	}
	std::vector<valuetype> bas_val_add_xi( n_dof_geometry[1], 0);
	std::vector<valuetype> bas_val_add_eta(n_dof_geometry[1], 0);
	Jxi[0] = Jeta[0] = 1;
	Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
	Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
	for (unsigned int l = 1; l < M; ++l){
	    Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
	    Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
	}
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
	    bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
	    bas_val_add_eta[l] = Jeta[l];
	}
	
	valuetype u_local = u(locp0)*r + u(locp1)*x + u(locp2)*y;
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l)
	    u_local -= 2 * (x*y * u(loce0+l) * bas_val_add_eta[l] +
			    y*r * u(loce1+l) * bas_val_add_eta[l] +
			    x*r * u(loce2+l) * bas_val_add_xi[l]);
	for (unsigned int l1 = 2; l1 <= M; ++l1)
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		u_local += bas_val[ind_index] * u(loc + ind_index);
	    }
	return u_local;
    }
    // std::cout << "searched face\n";
    
    // traverse 3-d geometries
    for (unsigned int ind_e = 0; ind_e < n_geometry[3]; ++ind_e){
	int &ind_v0 = index_geometry_onelement[ind_e][0][0];
	int &ind_v1 = index_geometry_onelement[ind_e][0][1];
	int &ind_v2 = index_geometry_onelement[ind_e][0][2];
	int &ind_v3 = index_geometry_onelement[ind_e][0][3];
	valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
	valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
	if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) continue;
    
	// valuetype u_local = calc_val_point_insideElement(correspondence, mesh, ind_e, pos, u);
	valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
	valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
	std::vector<valuetype> basis(n_index);
	// calculate immediate variable
	valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
	Jxi[0] = Jeta[0] = Jzeta[0] = 1;
	Jxi[1]  = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 0; l1 <= M; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		int aph3 = 2 * l1 + 2 * l2 - 1;
		Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
		for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		    Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
			/ calc_coefficient_a( aph3, -1, 1, l3);
		for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		    Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		    int ind_index = correspondence.index2number(index_now);
		    basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
		}
	    }
	}
	std::vector<valuetype> basis_actual(n_index, 0);
	for (int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];
    
	valuetype u_local = 0;
	// calculate function value
	std::vector<valuetype> coef_local(n_index, 0);
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_e][ind_index]; ++ind_nnz)
		coef_local[ind_index] += u(transform_ind_global2local[ind_e][ind_index][ind_nnz])
		    * transform_val_global2local[ind_e][ind_index][ind_nnz];

	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    u_local += coef_local[ind_index] * basis_actual[ind_index];
	return u_local;
    }

    std::cout << "error, cant find position of given pos in current SEM mesh, pos = (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    return 0;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_point(RegularMesh<3> &mesh, AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi, std::vector<valuetype> &n_occupation)
{
    int n_orbital = psi.size();
    
    // traverse 0-d geometries
    for (unsigned int ind_p = 0; ind_p < n_geometry[0]; ++ind_p){
	if (calc_length_line(mesh.point(ind_p), pos) < tol_zero){
	    int &loc = location_actualdof[location_geometry[ind_p]];
	    valuetype rho_local = 0;
	    for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
		// valuetype n_occupation = 2.0;
		rho_local += n_occupation[ind_orbital] * pow(psi[ind_orbital](loc), 2);
	    }
	    if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_p = " << ind_p << '\n';
	    return rho_local;
	}
    }
    // std::cout << "searched point\n";

    // traverse 1-d geometries
    for (unsigned int ind_e = 0; ind_e < n_geometry[1]; ++ind_e){
	int &ind_point0 = number_node[0][ind_e][0];
	int &ind_point1 = number_node[0][ind_e][1];
	valuetype dis = calc_length_line(mesh.point(ind_point0), mesh.point(ind_point1));
	valuetype dis0 = calc_length_line(mesh.point(ind_point0), pos);
	valuetype dis1 = calc_length_line(mesh.point(ind_point1), pos);
	if (fabs(dis0+dis1 - dis) > tol_zero) continue;
	int &loc0 = location_actualdof[location_geometry[ind_point0]];
	int &loc1 = location_actualdof[location_geometry[ind_point1]];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + ind_e]];
	
	valuetype x = dis0 / dis, r = 1 - x;
	valuetype xi = 2 * x - 1;
	valuetype Jxi[M+1];
	Jxi[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	
	valuetype rho_local = 0;
	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    // valuetype n_occupation = 2.0;
	    valuetype psi_local = psi[ind_orbital](loc0)*r + psi[ind_orbital](loc1)*x;
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		psi_local += 2 * Jxi[l1] * psi[ind_orbital](loc + l1-2);
	    rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
	}
	if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_edge = " << ind_e << '\n';
	return rho_local;
    }
    // std::cout << "searched edge\n";
    
    // traverse 2-d geometries
    for (unsigned int ind_f = 0; ind_f < n_geometry[2]; ++ind_f){
	int &ind_point0 = number_node[1][ind_f][0];
	int &ind_point1 = number_node[1][ind_f][1];
	int &ind_point2 = number_node[1][ind_f][2];
	valuetype area = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area0 = calc_area_triangle(pos, mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area1 = calc_area_triangle(mesh.point(ind_point0), pos, mesh.point(ind_point2));
	valuetype area2 = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), pos);
	if (fabs(area0+area1+area2 - area) > tol_zero) continue;
	int &ind_edge0 = number_edge[ind_f][0];
	int &ind_edge1 = number_edge[ind_f][1];
	int &ind_edge2 = number_edge[ind_f][2];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + n_geometry[1] + ind_f]];
	int &locp0 = location_actualdof[location_geometry[ind_point0]];
	int &locp1 = location_actualdof[location_geometry[ind_point1]];
	int &locp2 = location_actualdof[location_geometry[ind_point2]];
	int &loce0 = location_actualdof[location_geometry[n_geometry[0] + ind_edge0]];
	int &loce1 = location_actualdof[location_geometry[n_geometry[0] + ind_edge1]];
	int &loce2 = location_actualdof[location_geometry[n_geometry[0] + ind_edge2]];
	
	std::vector<valuetype> bas_val(n_dof_geometry[2], 0);
	valuetype x = area1 / area, y = area2 / area, r = 1 - x - y;
	valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
	valuetype Jxi[M+1], Jeta[M+1];
	Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 2; l1 <= M-1; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
	    }
	}
	std::vector<valuetype> bas_val_add_xi( n_dof_geometry[1], 0);
	std::vector<valuetype> bas_val_add_eta(n_dof_geometry[1], 0);
	Jxi[0] = Jeta[0] = 1;
	Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
	Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
	for (unsigned int l = 1; l < M; ++l){
	    Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
	    Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
	}
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
	    bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
	    bas_val_add_eta[l] = Jeta[l];
	}
	
	valuetype rho_local = 0;
	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    // valuetype n_occupation = 2.0;
	    valuetype psi_local;
	    psi_local = psi[ind_orbital](locp0) * r + psi[ind_orbital](locp1) * x + psi[ind_orbital](locp2) * y;
	    for (unsigned int l = 0; l < n_dof_geometry[1]; ++l)
		psi_local -= 2 * (x*y * psi[ind_orbital](loce0+l) * bas_val_add_eta[l] +
				  y*r * psi[ind_orbital](loce1+l) * bas_val_add_eta[l] +
				  x*r * psi[ind_orbital](loce2+l) * bas_val_add_xi[l]);
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		    unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		    psi_local += bas_val[ind_index] * psi[ind_orbital](loc + ind_index);
		}
	    rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
	}
	if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_f = " << ind_f << '\n';
	return rho_local;
    }
    // std::cout << "searched face\n";
    
    // traverse 3-d geometries
    for (unsigned int ind_e = 0; ind_e < n_geometry[3]; ++ind_e){
	int &ind_v0 = index_geometry_onelement[ind_e][0][0];
	int &ind_v1 = index_geometry_onelement[ind_e][0][1];
	int &ind_v2 = index_geometry_onelement[ind_e][0][2];
	int &ind_v3 = index_geometry_onelement[ind_e][0][3];
	valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
	valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
	if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) continue;
	
	valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
	valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
	std::vector<valuetype> basis(n_index);
	// calculate immediate variable
	valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
	Jxi[0] = Jeta[0] = Jzeta[0] = 1;
	Jxi[1]  = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 0; l1 <= M; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		int aph3 = 2 * l1 + 2 * l2 - 1;
		Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
		for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		    Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
			/ calc_coefficient_a( aph3, -1, 1, l3);
		for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		    Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		    int ind_index = correspondence.index2number(index_now);
		    basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
		}
	    }
	}
	std::vector<valuetype> basis_actual(n_index, 0);
	for (int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];
    
	valuetype rho_local = 0;
	for (int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    // valuetype n_occupation = 2.0;
	    // calculate function value
	    std::vector<valuetype> coef_local(n_index, 0);
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
		for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_e][ind_index]; ++ind_nnz)
		    coef_local[ind_index] += psi[ind_orbital](transform_ind_global2local[ind_e][ind_index][ind_nnz])
			* transform_val_global2local[ind_e][ind_index][ind_nnz];

	    valuetype psi_local = 0;
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
		psi_local += coef_local[ind_index] * basis_actual[ind_index];
	    // add contribution
	    rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
	}
	if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_ele = " << ind_e << '\n';
	return rho_local;
    }

    std::cout << "error, cant find position of given pos in current SEM mesh, pos = (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    return 0;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_point(Mesh<3> &mesh, AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi, std::vector<valuetype> &n_occupation)
{
    int n_orbital = psi.size();
    
    // traverse 0-d geometries
    for (unsigned int ind_p = 0; ind_p < n_geometry[0]; ++ind_p)
	if (calc_length_line(mesh.point(ind_p), pos) < tol_zero){
	    int &loc = location_actualdof[location_geometry[ind_p]];
	    valuetype rho_local = 0;
	    for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
		// valuetype n_occupation = 2.0;
		rho_local += n_occupation[ind_orbital] * pow(psi[ind_orbital](loc), 2);
	    }
	    if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_p = " << ind_p << '\n';
	    return rho_local;
	}
    // std::cout << "searched point\n";

    // traverse 1-d geometries
    for (unsigned int ind_e = 0; ind_e < n_geometry[1]; ++ind_e){
	int &ind_point0 = number_node[0][ind_e][0];
	int &ind_point1 = number_node[0][ind_e][1];
	valuetype dis = calc_length_line(mesh.point(ind_point0), mesh.point(ind_point1));
	valuetype dis0 = calc_length_line(mesh.point(ind_point0), pos);
	valuetype dis1 = calc_length_line(mesh.point(ind_point1), pos);
	if (fabs(dis0+dis1 - dis) > tol_zero) continue;
	int &loc0 = location_actualdof[location_geometry[ind_point0]];
	int &loc1 = location_actualdof[location_geometry[ind_point1]];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + ind_e]];
	
	valuetype x = dis0 / dis, r = 1 - x;
	valuetype xi = 2 * x - 1;
	valuetype Jxi[M+1];
	Jxi[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	
	valuetype rho_local = 0;
	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    // valuetype n_occupation = 2.0;
	    valuetype psi_local = psi[ind_orbital](loc0)*r + psi[ind_orbital](loc1)*x;
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		psi_local += 2 * Jxi[l1] * psi[ind_orbital](loc + l1-2);
	    rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
	}
	if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_edge = " << ind_e << '\n';
	return rho_local;
    }
    // std::cout << "searched edge\n";
    
    // traverse 2-d geometries
    for (unsigned int ind_f = 0; ind_f < n_geometry[2]; ++ind_f){
	int &ind_point0 = number_node[1][ind_f][0];
	int &ind_point1 = number_node[1][ind_f][1];
	int &ind_point2 = number_node[1][ind_f][2];
	valuetype area = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area0 = calc_area_triangle(pos, mesh.point(ind_point1), mesh.point(ind_point2));
	valuetype area1 = calc_area_triangle(mesh.point(ind_point0), pos, mesh.point(ind_point2));
	valuetype area2 = calc_area_triangle(mesh.point(ind_point0), mesh.point(ind_point1), pos);
	if (fabs(area0+area1+area2 - area) > tol_zero) continue;
	int &ind_edge0 = number_edge[ind_f][0];
	int &ind_edge1 = number_edge[ind_f][1];
	int &ind_edge2 = number_edge[ind_f][2];
	int &loc = location_actualdof[location_geometry[n_geometry[0] + n_geometry[1] + ind_f]];
	int &locp0 = location_actualdof[location_geometry[ind_point0]];
	int &locp1 = location_actualdof[location_geometry[ind_point1]];
	int &locp2 = location_actualdof[location_geometry[ind_point2]];
	int &loce0 = location_actualdof[location_geometry[n_geometry[0] + ind_edge0]];
	int &loce1 = location_actualdof[location_geometry[n_geometry[0] + ind_edge1]];
	int &loce2 = location_actualdof[location_geometry[n_geometry[0] + ind_edge2]];
	
	std::vector<valuetype> bas_val(n_dof_geometry[2], 0);
	valuetype x = area1 / area, y = area2 / area, r = 1 - x - y;
	valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
	valuetype Jxi[M+1], Jeta[M+1];
	Jxi[0] = Jeta[0] = 1; Jxi[1] = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 2; l1 <= M-1; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		bas_val[ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
	    }
	}
	std::vector<valuetype> bas_val_add_xi( n_dof_geometry[1], 0);
	std::vector<valuetype> bas_val_add_eta(n_dof_geometry[1], 0);
	Jxi[0] = Jeta[0] = 1;
	Jxi[1]  = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
	Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
	for (unsigned int l = 1; l < M; ++l){
	    Jxi[ l+1] = ((xi  - calc_coefficient_a(1, 1, 2, l)) *  Jxi[l] - calc_coefficient_a(1, 1, 3, l) *  Jxi[l-1]) / calc_coefficient_a(1, 1, 1, l);
	    Jeta[l+1] = ((eta - calc_coefficient_a(1, 1, 2, l)) * Jeta[l] - calc_coefficient_a(1, 1, 3, l) * Jeta[l-1]) / calc_coefficient_a(1, 1, 1, l);
	}
	for (unsigned int l = 0; l < n_dof_geometry[1]; ++l){
	    bas_val_add_xi[l]  = Jxi[l] * pow(1-y, l);
	    bas_val_add_eta[l] = Jeta[l];
	}
	
	valuetype rho_local = 0;
	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    // valuetype n_occupation = 2.0;
	    valuetype psi_local;
	    psi_local = psi[ind_orbital](locp0) * r + psi[ind_orbital](locp1) * x + psi[ind_orbital](locp2) * y;
	    for (unsigned int l = 0; l < n_dof_geometry[1]; ++l)
		psi_local -= 2 * (x*y * psi[ind_orbital](loce0+l) * bas_val_add_eta[l] +
				  y*r * psi[ind_orbital](loce1+l) * bas_val_add_eta[l] +
				  x*r * psi[ind_orbital](loce2+l) * bas_val_add_xi[l]);
	    for (unsigned int l1 = 2; l1 <= M; ++l1)
		for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
		    unsigned int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
		    psi_local += bas_val[ind_index] * psi[ind_orbital](loc + ind_index);
		}
	    rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
	}
	if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_f = " << ind_f << '\n';
	return rho_local;
    }
    // std::cout << "searched face\n";
    
    // traverse 3-d geometries
    for (unsigned int ind_e = 0; ind_e < n_geometry[3]; ++ind_e){
	int &ind_v0 = index_geometry_onelement[ind_e][0][0];
	int &ind_v1 = index_geometry_onelement[ind_e][0][1];
	int &ind_v2 = index_geometry_onelement[ind_e][0][2];
	int &ind_v3 = index_geometry_onelement[ind_e][0][3];
	valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
	valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
	valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
	if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) continue;

	valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
	valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
	std::vector<valuetype> basis(n_index);
	// calculate immediate variable
	valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
	Jxi[0] = Jeta[0] = Jzeta[0] = 1;
	Jxi[1]  = xi;
	for (unsigned int l1 = 1; l1 < M; ++l1)
	    Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
		/ calc_coefficient_a(-1, -1, 1, l1);
	for (unsigned int l1 = 0; l1 <= M; ++l1){
	    int aph2 = 2 * l1 - 1;
	    Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	    for (unsigned int l2 = 1; l2 < M-l1; ++l2)
		Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		    / calc_coefficient_a( aph2, -1, 1, l2);
	    for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
		int aph3 = 2 * l1 + 2 * l2 - 1;
		Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
		for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		    Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
			/ calc_coefficient_a( aph3, -1, 1, l3);
		for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		    Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		    int ind_index = correspondence.index2number(index_now);
		    basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
		}
	    }
	}
	std::vector<valuetype> basis_actual(n_index, 0);
	for (int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];
    
	valuetype rho_local = 0;
	for (int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    // valuetype n_occupation = 2.0;
	    // calculate function value
	    std::vector<valuetype> coef_local(n_index, 0);
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
		for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_e][ind_index]; ++ind_nnz)
		    coef_local[ind_index] += psi[ind_orbital](transform_ind_global2local[ind_e][ind_index][ind_nnz])
			* transform_val_global2local[ind_e][ind_index][ind_nnz];

	    valuetype psi_local = 0;
	    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
		psi_local += coef_local[ind_index] * basis_actual[ind_index];
	    // add contribution
	    rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
	}
	if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_ele = " << ind_e << '\n';
	return rho_local;
    }

    std::cout << "error, cant find position of given pos in current SEM mesh, pos = (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    return 0;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_point_insideElement(RegularMesh<3> &mesh, int ind_ele,
						  AFEPack::Point<3> &pos, Vector<valuetype> &u)
{
    int &ind_v0 = index_geometry_onelement[ind_ele][0][0];
    int &ind_v1 = index_geometry_onelement[ind_ele][0][1];
    int &ind_v2 = index_geometry_onelement[ind_ele][0][2];
    int &ind_v3 = index_geometry_onelement[ind_ele][0][3];
    valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
    valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
    if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) 
	std::cout << "error: interpolation from SEM solution to FEM solution\n";

    valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
    valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
    std::vector<valuetype> basis(n_index);
    // calculate immediate variable
    valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
    Jxi[0] = Jeta[0] = Jzeta[0] = 1;
    Jxi[1]  = xi;
    for (unsigned int l1 = 1; l1 < M; ++l1)
	Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
	    / calc_coefficient_a(-1, -1, 1, l1);
    for (unsigned int l1 = 0; l1 <= M; ++l1){
	int aph2 = 2 * l1 - 1;
	Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	for (unsigned int l2 = 1; l2 < M-l1; ++l2)
	    Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		/ calc_coefficient_a( aph2, -1, 1, l2);
	for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
	    int aph3 = 2 * l1 + 2 * l2 - 1;
	    Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
	    for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
		    / calc_coefficient_a( aph3, -1, 1, l3);
	    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		int ind_index = correspondence.index2number(index_now);
		basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
	    }
	}
    }
    std::vector<valuetype> basis_actual(n_index, 0);
    for (int ind_index = 0; ind_index < n_index; ++ind_index)
	for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
	    basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];
    
    valuetype u_local = 0;
    // calculate function value
    std::vector<valuetype> coef_local(n_index, 0);
    calc_coef_onElement(u, coef_local, ind_ele);
    for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	u_local += coef_local[ind_index] * basis_actual[ind_index];
    return u_local;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_val_point_insideElement(RegularMesh<3> &mesh, int ind_ele,
						  AFEPack::Point<3> &pos, std::vector<Vector<valuetype> > &psi, std::vector<valuetype> &n_occupation)
{
    int n_orbital = psi.size();

    int &ind_v0 = index_geometry_onelement[ind_ele][0][0];
    int &ind_v1 = index_geometry_onelement[ind_ele][0][1];
    int &ind_v2 = index_geometry_onelement[ind_ele][0][2];
    int &ind_v3 = index_geometry_onelement[ind_ele][0][3];
    valuetype volume = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume0 = calc_volume_tetrahedron(pos, mesh.point(ind_v1), mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume1 = calc_volume_tetrahedron(mesh.point(ind_v0), pos, mesh.point(ind_v2), mesh.point(ind_v3));
    valuetype volume2 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), pos, mesh.point(ind_v3));
    valuetype volume3 = calc_volume_tetrahedron(mesh.point(ind_v0), mesh.point(ind_v1), mesh.point(ind_v2), pos);
    if (fabs(volume0+volume1+volume2+volume3 - volume) > tol_zero) 
	std::cout << "error: interpolation from SEM solution to FEM solution\n";

    valuetype x = volume1 / volume, y = volume2 / volume, z = volume3 / volume;
    valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
    std::vector<valuetype> basis(n_index);
    // calculate immediate variable
    valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
    Jxi[0] = Jeta[0] = Jzeta[0] = 1;
    Jxi[1]  = xi;
    for (unsigned int l1 = 1; l1 < M; ++l1)
	Jxi[l1+1]  = ((xi - calc_coefficient_a(-1, -1, 2, l1)) *  Jxi[l1] - calc_coefficient_a(-1, -1, 3, l1) *  Jxi[l1-1])
	    / calc_coefficient_a(-1, -1, 1, l1);
    for (unsigned int l1 = 0; l1 <= M; ++l1){
	int aph2 = 2 * l1 - 1;
	Jeta[1]  = calc_generalized_jacobi_polynomial( aph2, -1, 1, eta);
	for (unsigned int l2 = 1; l2 < M-l1; ++l2)
	    Jeta[l2+1]  = ((eta - calc_coefficient_a( aph2, -1, 2, l2)) *  Jeta[l2] - calc_coefficient_a( aph2, -1, 3, l2) *  Jeta[l2-1])
		/ calc_coefficient_a( aph2, -1, 1, l2);
	for (unsigned int l2 = 0; l2 <= M-l1; ++l2){
	    int aph3 = 2 * l1 + 2 * l2 - 1;
	    Jzeta[1]  = calc_generalized_jacobi_polynomial( aph3, -1, 1, zeta);
	    for (unsigned int l3 = 1; l3 < M-l1-l2; ++l3)
		Jzeta[l3+1]  = ((zeta - calc_coefficient_a( aph3, -1, 2, l3)) *  Jzeta[l3] - calc_coefficient_a( aph3, -1, 3, l3) *  Jzeta[l3-1])
		    / calc_coefficient_a( aph3, -1, 1, l3);
	    for (unsigned int l3 = 0; l3 <= M-l1-l2; ++l3){
		Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
		int ind_index = correspondence.index2number(index_now);
		basis[ind_index] = pow(1-y-z, l1) * Jxi[l1] * pow(1-z, l2) * Jeta[l2] * Jzeta[l3];
	    }
	}
    }
    std::vector<valuetype> basis_actual(n_index, 0);
    for (int ind_index = 0; ind_index < n_index; ++ind_index)
	for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
	    basis_actual[transform_local[ind_index][ind_tl]] += weight_transform_local[ind_index][ind_tl] * basis[ind_index];
    
    valuetype rho_local = 0;
    for (int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	// valuetype n_occupation = 2.0;
	// calculate function value
	std::vector<valuetype> coef_local(n_index, 0);
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (unsigned int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
		coef_local[ind_index] += psi[ind_orbital](transform_ind_global2local[ind_ele][ind_index][ind_nnz])
		    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	valuetype psi_local = 0;
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    psi_local += coef_local[ind_index] * basis_actual[ind_index];
	// add contribution
	rho_local += n_occupation[ind_orbital] * pow(psi_local, 2);
    }
    if (rho_local > 1.0e4) std::cout << "rho_local = " << rho_local << ", ind_ele = " << ind_ele << '\n';
    return rho_local;
}

TEMPLATE_TSEM
void THIS_TSEM::calc_val_qp_onElement(Vector<valuetype> &src, std::vector<valuetype> &dst, unsigned int ind_ele)
{ // calculate value on quadrature point for sem solution on ind_ele-th element from source src to vector dst, suppose dst has size n_q_point[2]
    std::vector<valuetype> val_coef(n_index);
    calc_coef_onElement(src, val_coef, ind_ele);
    for (unsigned int p = 0; p < n_q_point[2]; ++p)
	dst[p] = 0;
    for (unsigned int p = 0; p < n_q_point[2]; ++p)
	for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
	    dst[p] += val_coef[ind_index] * basis_value_actual[p][ind_index];
}

#undef TEMPLATE_TSEM
#undef THIS_TSEM
