#include "../include/TetrahedralSEM.h"
// This file contains functions for building TSEM space
// Functions:
// build_basis_value()
// build_local_transform()
// build_geometry_info()
// build_global_transform()
// read_interp_info()


#define TEMPLATE_TSEM template<typename valuetype>
#define THIS_TSEM TSEM<valuetype>


TEMPLATE_TSEM
void THIS_TSEM::build_basis_value()
{ // calculate basis value on quadrature points

    // calculate the value of basis function at each quadrature point
    basis_value.resize(3);
    for (int ind = 0; ind < 3; ++ind)
        basis_value[ind].resize(n_q_point[ind]);
    // basis function value at 3-d quadrature points
    for (int p = 0; p < n_q_point[2]; ++p){
        valuetype x = QPoint[p][0], y = QPoint[p][1], z = QPoint[p][2];
        valuetype xi = 2*x/(1-y-z)-1, eta = 2*y/(1-z)-1, zeta = 2*z-1;
        basis_value[2][p].resize(n_index);
        // calculate $J_k^{-1,-1}(xi)$ for k = 0: M
        valuetype Jxi[M+1], Jeta[M+1], Jzeta[M+1];
        Jxi[0] = Jeta[0] = Jzeta[0] = 1;
        Jxi[1] = xi; // as J_1^{-1,-1}(xi) = xi
        for (int k = 1; k < M; ++k)
            Jxi[k+1] = ((xi - calc_coefficient_a(-1,-1,2,k)) * Jxi[k] - calc_coefficient_a(-1,-1,3,k) * Jxi[k-1]) / calc_coefficient_a(-1,-1,1,k);
        // traverse first component of multiindex
        for (int l1 = 0; l1 <= M; ++l1){
            int aph2 = 2 * l1 - 1; // alpha for the generalized jacobi polynomial of eta
            // calculate value of generalized jacobi polynomial Jeta
            Jeta[1] = calc_generalized_jacobi_polynomial(aph2, -1, 1, eta);
            for (int k = 1; k < M-l1; ++k)
                Jeta[k+1] = ((eta - calc_coefficient_a(aph2,-1,2,k)) * Jeta[k] - calc_coefficient_a(aph2,-1,3,k) * Jeta[k-1]) / calc_coefficient_a(aph2,-1,1,k);
            // traverse second component
            for (int l2 = 0; l2 <= M-l1; ++l2){
                int aph3 = 2 * l1 + 2 * l2 - 1;
                Jzeta[1] = calc_generalized_jacobi_polynomial(aph3, -1, 1, zeta);
                for (int k = 1; k < M-l1-l2; ++k)
                    Jzeta[k+1] = ((zeta - calc_coefficient_a(aph3,-1,2,k)) * Jzeta[k] - calc_coefficient_a(aph3,-1,3,k) * Jzeta[k-1]) / calc_coefficient_a(aph3,-1,1,k);
                // traverse third component
                for (int l3 = 0; l3 <= M-l1-l2; ++l3){
                    Multiindex<3> index_now = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
                    basis_value[2][p][correspondence.index2number(index_now)] = pow(1-y-z,l1)*Jxi[l1] * pow(1-z,l2)*Jeta[l2] * Jzeta[l3];
		    // if (fabs(basis_value[2][p][correspondence.index2number(index_now)]) > 1.0e8)
		    // 	std::cout << "find large basis_value, Jxi[" << l1 << "] = " << Jxi[l1]
		    // 		  << ", Jeta[" << l2 << "] = " << Jeta[l2]
		    // 		  << ", Jzeta[" << l3 << "] = " << Jzeta[l3]
		    // 		  << ", a(" << aph2 << ", -1, 1, " << l2-1 << ") = " << calc_coefficient_a(aph2, -1, 1, l2-1)
		    // 		  << ", a(" << aph3 << ", -1, 1, " << l3-1 << ") = " << calc_coefficient_a(aph3, -1, 1, l3-1)
		    // 		  << '\n';
                }
            }
        }
    }
    // 2-d quadrature point
    for (int p = 0; p < n_q_point[1]; ++p){
        basis_value[1][p].resize(n_dof_geometry[2]);
        valuetype x = QPoint_Barycentric[1][p][0], y = QPoint_Barycentric[1][p][1];
        valuetype xi = 2*x/(1-y)-1, eta = 2*y-1;
        valuetype Jxi[M+1], Jeta[M+1];
        Jxi[0] = Jeta[0] = 1;
        Jxi[1] = xi;
        for (int l = 1; l < M; ++l)
            Jxi[l+1] = ((xi - calc_coefficient_a(-1,-1,2,l))*Jxi[l] - calc_coefficient_a(-1,-1,3,l)*Jxi[l-1]) / calc_coefficient_a(-1,-1,1,l);
        for (int l1 = 2; l1 <= M; ++l1){
            int aph = 2 * l1 - 1;
            Jeta[1] = calc_generalized_jacobi_polynomial(aph, -1, 1, eta);
            for (int l = 1; l < M; ++l)
                Jeta[l+1] = ((eta - calc_coefficient_a(aph,-1,2,l))*Jeta[l] - calc_coefficient_a(aph,-1,3,l)*Jeta[l-1]) / calc_coefficient_a(aph,-1,1,l);
            for (int l2 = 1; l2 <= M-l1; ++l2){
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2 - 1;
                basis_value[1][p][ind_index] = 2 * pow(1-y, l1) * Jxi[l1] * Jeta[l2];
            }
        }
    }
    // 1-d quadrature point
    for (int p = 0; p < n_q_point[0]; ++p){
        basis_value[0][p].resize(n_dof_geometry[1]);
        valuetype xi = QPoint_Barycentric[0][p][0] * 2 - 1;
        valuetype Jxi[M+1];
        Jxi[0] = 1;
        Jxi[1] = xi;
        for (int l = 1; l < M; ++l)
            Jxi[l+1] = ((xi - calc_coefficient_a(-1,-1,2,l))*Jxi[l] - calc_coefficient_a(-1,-1,3,l)*Jxi[l-1]) / calc_coefficient_a(-1,-1,1,l);
        for (int l = 2; l <= M; ++l)
            basis_value[0][p][l-2] = 2 * Jxi[l];
    }
    std::cerr << "calculate function value of generalized jacobi polynomial at quadrature points\n";
    // for (int ind_index = 0; ind_index < n_index; ++ind_index){
    // 	std::cout << "ind_index = " << ind_index << ", basis_value:";
    // 	for (int p = 0; p < n_q_point[2]; ++p)
    // 	    std::cout << ' ' << basis_value[2][p][ind_index];
    // 	std::cout << '\n';
    // }

    // calculate the value of generalized jacobi polynomial for the interpolation of function
    basis_value_interp.resize(3); // basis value for interpolation
    for (int ind = 0; ind < 3; ++ind){
        basis_value_interp[ind].resize(n_q_point[ind]);
        for (int p = 0; p < n_q_point[ind]; ++p)
            basis_value_interp[ind][p].resize(n_dof_geometry[ind+1]);
    }
    // for the interpolation of the dof on edge
    for (int p = 0; p < n_q_point[0]; ++p){
        valuetype xp = 2 * QPoint_Barycentric[0][p][0] - 1;
        basis_value_interp[0][p][0] = 1;
        basis_value_interp[0][p][1] = calc_generalized_jacobi_polynomial(1, 1, 1, xp);
        for (int l = 1; l < M-2; ++l)
            basis_value_interp[0][p][l+1] = ((xp - calc_coefficient_a(1,1,2,l)) * basis_value_interp[0][p][l]
                                             - calc_coefficient_a(1,1,3,l) * basis_value_interp[0][p][l-1]) / calc_coefficient_a(1,1,1,l);
    }
    // for the interpolation of the dof on face
    // additional basis value, [p][l][0]: $(1-y)^l*J_{l-2}^{1,1}(\xi_p)$, [p][l][1]: $J_{l-2}^{1,1}(\eta_p)$
    basis_value_addition.resize(n_q_point[1]);
    for (int p = 0; p < n_q_point[1]; ++p){
        basis_value_addition[p].resize(n_dof_geometry[1]);
        for (int l = 0; l < n_dof_geometry[1]; ++l)
            basis_value_addition[p][l].resize(2);
    }
    for (int p = 0; p < n_q_point[1]; ++p){
        valuetype xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1];
        valuetype xi = 2*xp/(1-yp) - 1, eta = 2*yp - 1;
        valuetype Jxi[M-1], Jeta[M];
        // evaluate basis_value_interp
        Jxi[0] = Jeta[0] = 1;
        Jxi[1] = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
        for (int l = 1; l < M-2; ++l) // evaluate Jxi by recursion relation
            Jxi[l+1] = ((xi - calc_coefficient_a(1,1,2,l)) * Jxi[l] - calc_coefficient_a(1,1,3,l) * Jxi[l-1]) / calc_coefficient_a(1,1,1,l);
        for (int l1 = 2; l1 <= M; ++l1){
            Jeta[1] = calc_generalized_jacobi_polynomial(2*l1-1, 1, 1, eta);
            for (int l2 = 1; l2 < M-l1; ++l2) // evaluate Jeta by recursion relation, Jeta[M-l1] is un-used
                Jeta[l2+1] = ((eta - calc_coefficient_a(2*l1-1,1,2,l2)) * Jeta[l2] - calc_coefficient_a(2*l1-1,1,3,l2) * Jeta[l2-1]) / calc_coefficient_a(2*l1-1,1,1,l2);
            for (int l2 = 1; l2 <= M-l1; ++l2){ // assign basis_value_interp, corresponds to (l1-2, l2-1)
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1;
                basis_value_interp[1][p][ind_index] = pow(1-eta, l1-2) * Jxi[l1-2] * Jeta[l2-1];
            }
        }
        // evaluate basis_value_addition, which correspond to $(1-yp)^{l-2}J_{l-2}^{1,1}(xi)$ and $J_{l-2}^{1,1}(eta)$ for l = 2: M
        for (int l = 2; l <= M; ++l)
            basis_value_addition[p][l-2][0] = pow(1-yp, l-2) * Jxi[l-2]; // as Jxi is the same
        Jeta[1] = calc_generalized_jacobi_polynomial(1, 1, 1, eta);
        for (int l = 1; l < M-2; ++l)
            Jeta[l+1] = ((eta - calc_coefficient_a(1,1,2,l)) * Jeta[l] - calc_coefficient_a(1,1,3,l) * Jeta[l-1]) / calc_coefficient_a(1,1,1,l);
        for (int l = 2; l <= M; ++l)
            basis_value_addition[p][l-2][1] = Jeta[l-2];
    }
    // for the interpolation of the interior dof
    for (int p = 0; p < n_q_point[2]; ++p){
        valuetype xp = QPoint[p][0], yp = QPoint[p][1], zp = QPoint[p][2];
        valuetype xi = 2*xp/(1-yp-zp)-1, eta = 2*yp/(1-zp)-1, zeta = 2*zp-1;
        valuetype Jxi[M-1], Jeta[M], Jzeta[M];
        Jxi[0] = Jeta[0] = Jzeta[0] = 1;
        Jxi[1] = calc_generalized_jacobi_polynomial(1, 1, 1, xi);
        for (int l = 1; l < M-2; ++l) // l1 = 2: M -> l1-2 = 0: M-2
            Jxi[l+1] = ((xi - calc_coefficient_a(1,1,2,l)) * Jxi[l] - calc_coefficient_a(1,1,3,l) * Jxi[l-1]) / calc_coefficient_a(1,1,1,l);
        for (int l1 = 2; l1 <= M; ++l1){
            Jeta[1] = calc_generalized_jacobi_polynomial(2*l1-1, 1, 1, eta);
            for (int l = 1; l < M - l1; ++l) // in fact, we consider l2 = 1: M-l1 -> l2-1 = 0: M-l1-1, so Jeta[M-l1] is un-used
                Jeta[l+1] = ((eta - calc_coefficient_a(2*l1-1,1,2,l)) * Jeta[l] - calc_coefficient_a(2*l1-1,1,3,l) * Jeta[l-1]) / calc_coefficient_a(2*l1-1,1,1,l);
            for (int l2 = 1; l2 <= M - l1; ++l2){
                int aph = 2 * l1 + 2 * l2 - 1;
                Jzeta[1] = calc_generalized_jacobi_polynomial(aph, 1, 1, zeta);
                for (int l = 1; l < M - l1 - l2; ++l) // similarly, Jzeta[M-l1-l2] is un-used
                    Jzeta[l+1] = ((zeta - calc_coefficient_a(aph,1,2,l)) * Jzeta[l] - calc_coefficient_a(aph,1,3,l) * Jzeta[l-1]) / calc_coefficient_a(aph,1,1,l);
                for (int l3 = 1; l3 <= M - l1 - l2; ++l3){
                    int ind_index = correspondence.index2number(Unitary_Multiindex[0] * (l1-2) + Unitary_Multiindex[1] * (l2-1) + Unitary_Multiindex[2] * (l3-1));
                    basis_value_interp[2][p][ind_index] = Jxi[l1-2] * Jeta[l2-1]*pow(1-eta,l1-2) * Jzeta[l3-1]*pow(1-zeta,l1+l2-3);
                }
            }
        }
    }
    std::cerr << "calculate the basis value for interpolation\n";

    // calculate gradient of basis function on quadrature points
    basis_gradient.resize(n_q_point[2], std::vector<std::vector<valuetype> > (n_index, std::vector<valuetype> (3, 0)));
    // prepare auxiliary varaible: coefficient and point value of generalize Jacobi polynomial
    valuetype coeff_D1[n_index];
    valuetype coeff_D2[n_index][2]; // p = 0, 1
    valuetype coeff_D3[n_index][4]; // (p, q) = (0, 0) (1, 0) (0, 1) (1, 1)
    for (int ind_index = 0; ind_index < n_index; ++ind_index){
	Multiindex<3> index_now = correspondence.number2index(ind_index);
	coeff_D1[ind_index] = calc_coefficient_D(0, 0, index_now);
	coeff_D2[ind_index][0] = calc_coefficient_D(3, 0, index_now);
	coeff_D2[ind_index][1] = calc_coefficient_D(3, 1, index_now);
	coeff_D3[ind_index][0] = calc_coefficient_D(5, 0, index_now);
	coeff_D3[ind_index][1] = calc_coefficient_D(5, 1, index_now);
	coeff_D3[ind_index][2] = calc_coefficient_D(5, 2, index_now);
	coeff_D3[ind_index][3] = calc_coefficient_D(5, 3, index_now);
    }
    // traverse quadrature points, calculate val_Jacobi, then assign basis_gradient[][p][]
    for (int p = 0; p < n_q_point[2]; ++p){
	// setup
	std::vector<std::vector<valuetype> > curlicue_J(3, std::vector<valuetype> (n_index));
	// assign curlicue J by calculating generalize Jacobi polynomials
	valuetype x1 = QPoint[p][0], x2 = QPoint[p][1], x3 = QPoint[p][2];
	valuetype t1 = 1-x2-x3, t2 = 1-x3;
	valuetype xi = 2*x1/t1-1, eta = 2*x2/t2-1, zeta = 2*x3-1;
	valuetype pow_t1[M+1], pow_t2[M+1];
	pow_t1[0] = (valuetype) 1; pow_t2[0] = (valuetype) 1;
	for (int l = 1; l <= M; ++l){
	    pow_t1[l] = pow_t1[l-1] * t1;
	    pow_t2[l] = pow_t2[l-1] * t2;
	}
	valuetype Jxi  [2][M+1]; // [0][]: ^{0,0};          [1][]: ^{0,-1}
	valuetype Jeta [3][M+1]; // [0][]: ^{2l1+1,-1};     [1][]: ^{2l1,0};     [2][]: ^{2l1,-1}
	valuetype Jzeta[2][M+1]; // [0][]: ^{2l1+2l2+1,-1}; [1][]: ^{2l1+2l2,0}
	Jxi[0][0] = Jxi[1][0] = (valuetype) 1;
	Jeta[0][0] = Jeta[1][0] = Jeta[2][0] = (valuetype) 1;
	Jzeta[0][0] = Jzeta[1][0] = (valuetype) 1;
	Jxi[0][1] = calc_generalized_jacobi_polynomial(0,  0, 1, xi);
	Jxi[1][1] = calc_generalized_jacobi_polynomial(0, -1, 1, xi);
	for (int l = 1; l < M; ++l){
	    Jxi[0][l+1] = ((xi - calc_coefficient_a(0, 0,2,l)) * Jxi[0][l] - calc_coefficient_a(0, 0,3,l) * Jxi[0][l-1]) / calc_coefficient_a(0, 0,1,l);
	    Jxi[1][l+1] = ((xi - calc_coefficient_a(0,-1,2,l)) * Jxi[1][l] - calc_coefficient_a(0,-1,3,l) * Jxi[1][l-1]) / calc_coefficient_a(0,-1,1,l);
	}
	for (int l1 = 0; l1 <= M; ++l1){
	    int aph2_0 = 2*l1+1, aph2_1 = 2*l1;
	    Jeta[0][1] = calc_generalized_jacobi_polynomial(aph2_0, -1, 1, eta);
	    Jeta[1][1] = calc_generalized_jacobi_polynomial(aph2_1,  0, 1, eta);
	    Jeta[2][1] = calc_generalized_jacobi_polynomial(aph2_1, -1, 1, eta);
	    for (int l = 1; l < M-l1; ++l){
		Jeta[0][l+1] = ((eta - calc_coefficient_a(aph2_0,-1,2,l)) * Jeta[0][l] - calc_coefficient_a(aph2_0,-1,3,l) * Jeta[0][l-1]) / calc_coefficient_a(aph2_0,-1,1,l);
		Jeta[1][l+1] = ((eta - calc_coefficient_a(aph2_1, 0,2,l)) * Jeta[1][l] - calc_coefficient_a(aph2_1, 0,3,l) * Jeta[1][l-1]) / calc_coefficient_a(aph2_1, 0,1,l);
		Jeta[2][l+1] = ((eta - calc_coefficient_a(aph2_1,-1,2,l)) * Jeta[2][l] - calc_coefficient_a(aph2_1,-1,3,l) * Jeta[2][l-1]) / calc_coefficient_a(aph2_1,-1,1,l);
	    }
	    for (int l2 = 0; l2 <= M-l1; ++l2){
		int aph3_0 = 2*l1+2*l2+1, aph3_1 = 2*l1+2*l2;
		Jzeta[0][1] = calc_generalized_jacobi_polynomial(aph3_0, -1, 1, zeta);
		Jzeta[1][1] = calc_generalized_jacobi_polynomial(aph3_1,  0, 1, zeta);
		for (int l = 1; l < M-l1-l2; ++l){
                    Jzeta[0][l+1] = ((zeta-calc_coefficient_a(aph3_0,-1,2,l)) * Jzeta[0][l] - calc_coefficient_a(aph3_0,-1,3,l) * Jzeta[0][l-1]) / calc_coefficient_a(aph3_0,-1,1,l);
                    Jzeta[1][l+1] = ((zeta-calc_coefficient_a(aph3_1, 0,2,l)) * Jzeta[1][l] - calc_coefficient_a(aph3_1, 0,3,l) * Jzeta[1][l-1]) / calc_coefficient_a(aph3_1, 0,1,l);
		}
		for (int l3 = 0; l3 <= M-l1-l2; ++l3){
		    Multiindex<3> index_now = Unitary_Multiindex[0]*l1 + Unitary_Multiindex[1]*l2 + Unitary_Multiindex[2]*l3;
		    int ind_index = correspondence.index2number(index_now);
		    curlicue_J[0][ind_index] = pow_t1[l1]*Jxi[0][l1] * pow_t2[l2]*Jeta[0][l2] * Jzeta[0][l3];
		    curlicue_J[1][ind_index] = pow_t1[l1]*Jxi[1][l1] * pow_t2[l2]*Jeta[1][l2] * Jzeta[0][l3];
		    curlicue_J[2][ind_index] = pow_t1[l1]*Jxi[1][l1] * pow_t2[l2]*Jeta[2][l2] * Jzeta[1][l3];
		}
	    }
	}
	// assign basis_gradient[p][][]
	for (int ind_index = 0; ind_index < n_index; ++ind_index){
	    Multiindex<3> index_tmp = correspondence.number2index(ind_index);
	    // 1st component
	    index_tmp.index[0]--;
	    if (0 <= index_tmp.index[0])
		basis_gradient[p][ind_index][0] += coeff_D1[ind_index]    * curlicue_J[0][correspondence.index2number(index_tmp)];
	    index_tmp.index[0]++;
	    // 2nd component
	    // p = 0
	    index_tmp.index[1]--;
	    if (0 <= index_tmp.index[1])
		basis_gradient[p][ind_index][1] += coeff_D2[ind_index][0] * curlicue_J[1][correspondence.index2number(index_tmp)];
	    index_tmp.index[1]++;
	    // p = 1;
	    index_tmp.index[0]--;
	    if (0 <= index_tmp.index[0])
		basis_gradient[p][ind_index][1] += coeff_D2[ind_index][1] * curlicue_J[1][correspondence.index2number(index_tmp)];
	    index_tmp.index[0]++;
	    // 3rd component
            // p = 0, q = 0
	    index_tmp.index[2]--;
	    if (0 <= index_tmp.index[2])
		basis_gradient[p][ind_index][2] += coeff_D3[ind_index][0] * curlicue_J[2][correspondence.index2number(index_tmp)];
	    index_tmp.index[2]++;
	    // p = 1, q = 0
	    index_tmp.index[0]--; index_tmp.index[1]++; index_tmp.index[2]--;
	    if (0 <= index_tmp.index[0] && 0 <= index_tmp.index[2])
		basis_gradient[p][ind_index][2] += coeff_D3[ind_index][1] * curlicue_J[2][correspondence.index2number(index_tmp)];
	    index_tmp.index[0]++; index_tmp.index[1]--; index_tmp.index[2]++;
	    // p = 0, q = 1
	    index_tmp.index[1]--;
	    if (0 <= index_tmp.index[1])
		basis_gradient[p][ind_index][2] += coeff_D3[ind_index][2] * curlicue_J[2][correspondence.index2number(index_tmp)];
	    index_tmp.index[1]++;
	    // p = 1, q = 1
	    index_tmp.index[0]--;
	    if (0 <= index_tmp.index[0])
		basis_gradient[p][ind_index][2] += coeff_D3[ind_index][3] * curlicue_J[2][correspondence.index2number(index_tmp)];
	}
    }
    std::cerr << "calculate the basis gradient\n";
}


TEMPLATE_TSEM
void THIS_TSEM::build_local_transform()
{
    // construct tranform_local, weight_transform_local: generalize jacobi polynomial -> basis function
    n_transform_local.resize(n_index);
    transform_local.resize(n_index);
    weight_transform_local.resize(n_index);
    // vertex modes
    n_transform_local[0] = 4; n_transform_local[1] = 2; n_transform_local[2] = 3; n_transform_local[3] = 4;
    for (int ind = 0; ind <= 3; ++ind){
        transform_local[ind].resize(n_transform_local[ind]);
        weight_transform_local[ind].resize(n_transform_local[ind]);
    }
    transform_local[0][0] = 0; weight_transform_local[0][0] = (valuetype) 0.125; // $J_{0,0,0} -> \varphi_{0,0,0}$
    transform_local[0][1] = 1; weight_transform_local[0][1] = (valuetype) 0.125; // $J_{0,0,0} -> \varphi_{1,0,0}$
    transform_local[0][2] = 2; weight_transform_local[0][2] = (valuetype) 0.25;  // $J_{0,0,0} -> \varphi_{0,1,0}$
    transform_local[0][3] = 3; weight_transform_local[0][3] = (valuetype) 0.5;   // $J_{0,0,0} -> \varphi_{0,0,1}$
    transform_local[1][0] = 0; weight_transform_local[1][0] = (valuetype) -0.5;  // $J_{1,0,0} -> \varphi_{0,0,0}$
    transform_local[1][1] = 1; weight_transform_local[1][1] = (valuetype) 0.5;   // $J_{1,0,0} -> \varphi_{1,0,0}$
    transform_local[2][0] = 0; weight_transform_local[2][0] = (valuetype) -0.25; // $J_{0,1,0} -> \varphi_{0,0,0}$
    transform_local[2][1] = 1; weight_transform_local[2][1] = (valuetype) -0.25; // $J_{0,1,0} -> \varphi_{1,0,0}$
    transform_local[2][2] = 2; weight_transform_local[2][2] = (valuetype) 0.5;   // $J_{1,1,0} -> \varphi_{0,1,0}$
    transform_local[3][0] = 0; weight_transform_local[3][0] = (valuetype) -0.125;// $J_{0,0,1} -> \varphi_{0,0,0}$
    transform_local[3][1] = 1; weight_transform_local[3][1] = (valuetype) -0.125;// $J_{0,0,1} -> \varphi_{1,0,0}$
    transform_local[3][2] = 2; weight_transform_local[3][2] = (valuetype) -0.25; // $J_{0,0,1} -> \varphi_{0,1,0}$
    transform_local[3][3] = 3; weight_transform_local[3][3] = (valuetype) 0.5;   // $J_{0,0,1} -> \varphi_{0,0,1}$
    // edge modes
    for (int l = 2; l <= M; ++l){
        // $J_{l,0,0}$ -> $\varphi_{l,0,0}$ 
        Multiindex<3> index_tmp3 = Unitary_Multiindex[0] * l;
        int ind_tmp3 = correspondence.index2number(index_tmp3);
        n_transform_local[ind_tmp3] = 1;
        transform_local[ind_tmp3].resize(n_transform_local[ind_tmp3]);
        weight_transform_local[ind_tmp3].resize(n_transform_local[ind_tmp3]);
        transform_local[ind_tmp3][0] = ind_tmp3; weight_transform_local[ind_tmp3][0] = (valuetype) 2;
        // $J_{0,l,0}$, $J_{1,l-1,0}$ -> $\varphi_{0,l,0}$, $\varphi_{1,l-1,0}$
        Multiindex<3> index_tmp1 = Unitary_Multiindex[1] * l;
        Multiindex<3> index_tmp2 = Unitary_Multiindex[0] + Unitary_Multiindex[1] * (l-1);
        int ind_tmp1 = correspondence.index2number(index_tmp1);
        int ind_tmp2 = correspondence.index2number(index_tmp2);
        n_transform_local[ind_tmp1] = 2;
        n_transform_local[ind_tmp2] = 2;
        transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
        transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
        weight_transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
        weight_transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
        transform_local[ind_tmp1][0] = ind_tmp1; weight_transform_local[ind_tmp1][0] = (valuetype) 1;         // $J_{0,l,0}$   -> $\varphi_{0,l,0}$
        transform_local[ind_tmp1][1] = ind_tmp2; weight_transform_local[ind_tmp1][1] = (valuetype) 1;         // $J_{0,l,0}$   -> $\varphi_{1,l-1,0}$
        transform_local[ind_tmp2][0] = ind_tmp1; weight_transform_local[ind_tmp2][0] = ((valuetype) (l-1))/l; // $J_{1,l-1,0}$ -> $\varphi_{0,l,0}$
        transform_local[ind_tmp2][1] = ind_tmp2; weight_transform_local[ind_tmp2][1] =((valuetype) -(l-1))/l; // $J_{1,l-1,0}$ -> $\varphi_{1,l-1,0}$
        // $J_{0,0,l}$, $J_{1,0,l-1}$, $J_{0,1,l-1}$ -> $\varphi_{0,0,l}$, $varphi_{1,0,l-1}$, $varphi_{0,1,l-1}$
        index_tmp1 = Unitary_Multiindex[2] * l;
        index_tmp2 = Unitary_Multiindex[0] + Unitary_Multiindex[2] * (l-1);
        index_tmp3 = Unitary_Multiindex[1] + Unitary_Multiindex[2] * (l-1);
        ind_tmp1 = correspondence.index2number(index_tmp1);
        ind_tmp2 = correspondence.index2number(index_tmp2);
        ind_tmp3 = correspondence.index2number(index_tmp3);
        n_transform_local[ind_tmp1] = 3;
        n_transform_local[ind_tmp2] = 2;
        n_transform_local[ind_tmp3] = 3;
        transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
        transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
        transform_local[ind_tmp3].resize(n_transform_local[ind_tmp3]);
        weight_transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
        weight_transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
        weight_transform_local[ind_tmp3].resize(n_transform_local[ind_tmp3]);
        transform_local[ind_tmp1][0] = ind_tmp1; weight_transform_local[ind_tmp1][0] = (valuetype) 0.5;           // $J_{0,0,l}$   -> $\varphi_{0,0,l}$
        transform_local[ind_tmp1][1] = ind_tmp2; weight_transform_local[ind_tmp1][1] = (valuetype) 0.5;           // $J_{0,0,l}$   -> $\varphi_{1,0,l-1}$
        transform_local[ind_tmp1][2] = ind_tmp3; weight_transform_local[ind_tmp1][2] = (valuetype) 1;             // $J_{0,0,l}$   -> $\varphi_{0,1,l-1}$
        transform_local[ind_tmp2][0] = ind_tmp1; weight_transform_local[ind_tmp2][0] = ((valuetype) (l-1))/l;     // $J_{1,0,l-1}$ -> $\varphi_{0,0,l}$
        transform_local[ind_tmp2][1] = ind_tmp2; weight_transform_local[ind_tmp2][1] =((valuetype) -(l-1))/l;     // $J_{1,0,l-1}$ -> $\varphi_{1,0,l-1}$
        transform_local[ind_tmp3][0] = ind_tmp1; weight_transform_local[ind_tmp3][0] = ((valuetype) (l-1))/(2*l); // $J_{0,1,l-1}$ -> $\varphi_{0,0,l}$
        transform_local[ind_tmp3][1] = ind_tmp2; weight_transform_local[ind_tmp3][1] = ((valuetype) (l-1))/(2*l); // $J_{0,1,l-1}$ -> $\varphi_{1,0,l-1}$
        transform_local[ind_tmp3][2] = ind_tmp3; weight_transform_local[ind_tmp3][2] =((valuetype) -(l-1))/l;     // $J_{0,1,l-1}$ -> $\varphi_{0,1,l-1}$
    }
    // face modes
    for (int l1 = 2; l1 <= M; ++l1)
        for (int l2 = 1; l2 <= M-l1; ++l2){
            // $J_{0,l1,l2}$, $J_{1,l1-1,l2}$ -> $\varphi_{0,l1,l2}$, $varphi_{1,l1-1,l2}$
            Multiindex<3> index_tmp1 = Unitary_Multiindex[1] * l1 + Unitary_Multiindex[2] * l2;
            Multiindex<3> index_tmp2 = Unitary_Multiindex[0] + Unitary_Multiindex[1] * (l1-1) + Unitary_Multiindex[2] * l2;
            int ind_tmp1 = correspondence.index2number(index_tmp1);
            int ind_tmp2 = correspondence.index2number(index_tmp2);
            n_transform_local[ind_tmp1] = (valuetype) 2;
            n_transform_local[ind_tmp2] = (valuetype) 2;
            transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
            transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
            weight_transform_local[ind_tmp1].resize(n_transform_local[ind_tmp1]);
            weight_transform_local[ind_tmp2].resize(n_transform_local[ind_tmp2]);
            transform_local[ind_tmp1][0] = ind_tmp1; weight_transform_local[ind_tmp1][0] = (valuetype) 1;           // $J_{0,l1,l2}$   -> $\varphi_{0,l1,l2}$
            transform_local[ind_tmp1][1] = ind_tmp2; weight_transform_local[ind_tmp1][1] = (valuetype) 1;           // $J_{0,l1,l2}$   -> $\varphi_{1,l1-1,l2}$
            transform_local[ind_tmp2][0] = ind_tmp1; weight_transform_local[ind_tmp2][0] = ((valuetype) (l1-1))/l1; // $J_{1,l1-1,l2}$ -> $\varphi_{0,l1,l2}$
            transform_local[ind_tmp2][1] = ind_tmp2; weight_transform_local[ind_tmp2][1] =((valuetype) -(l1-1))/l1; // $J_{1,l1-1,l2}$ -> $\varphi_{1,l1-1,l2}$
            // $J_{l1,0,l2}$ -> $\varphi_{l1,0,l2}$ and $J_{l1,l2,0}$ -> $\varphi_{l1,l2,0}$
            for (int ind = 1; ind <= 2; ++ind){
                Multiindex<3> index_tmp = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[ind] * l2;
                int ind_tmp = correspondence.index2number(index_tmp);
                n_transform_local[ind_tmp] = (valuetype) 1;
                transform_local[ind_tmp].resize(n_transform_local[ind_tmp]);
                weight_transform_local[ind_tmp].resize(n_transform_local[ind_tmp]);
                transform_local[ind_tmp][0] = ind_tmp; weight_transform_local[ind_tmp][0] = (valuetype) 2;
            }
        }
    // interior modes
    for (int l1 = 2; l1 <= M; ++l1)
        for (int l2 = 1; l2 <= M-l1; ++l2)
            for (int l3 = 1; l3 <= M-l1-l2; ++l3){
                Multiindex<3> index_tmp = Unitary_Multiindex[0] * l1 + Unitary_Multiindex[1] * l2 + Unitary_Multiindex[2] * l3;
                int ind_tmp = correspondence.index2number(index_tmp);
                n_transform_local[ind_tmp] = 1;
                transform_local[ind_tmp].resize(n_transform_local[ind_tmp]);
                weight_transform_local[ind_tmp].resize(n_transform_local[ind_tmp]);
                transform_local[ind_tmp][0] = ind_tmp; weight_transform_local[ind_tmp][0] = (valuetype) 1;
            }

    // for (int ind_index = 0; ind_index < n_index; ++ind_index){
    // 	std::cout << "ind_index = " << ind_index << ", n_transform_local = " << n_transform_local[ind_index]
    // 		  << ", weight_transform_local:";
    // 	for (int i = 0; i < n_transform_local[ind_index]; ++i)
    // 	    std::cout << ' ' << weight_transform_local[ind_index][i];
    // 	std::cout << '\n';
    // }
    // generate actual basis function value at local fem element, by transform_local and weigth_transform_local
    basis_value_actual.resize(n_q_point[2], std::vector<valuetype> (n_index, 0));
    // traverse 3-d quadrature point, use transform_local and weight_transform_local, calculate basis_value_actual
    for (int p = 0; p < n_q_point[2]; ++p)
        for (int i = 0; i < n_index; ++i)
            for (int j = 0; j < n_transform_local[i]; ++j)
                basis_value_actual[p][transform_local[i][j]] += weight_transform_local[i][j] * basis_value[2][p][i];
    std::cerr << "assigned basis_value_actual\n";

    basis_gradient_actual.resize(n_q_point[2], std::vector<std::vector<valuetype> > (n_index, std::vector<valuetype> (3, 0)));
    for (int p = 0; p < n_q_point[2]; ++p)
	for (int ind_index = 0; ind_index < n_index; ++ind_index)
	    for (int ind_tl = 0; ind_tl < n_transform_local[ind_index]; ++ind_tl)
		for (int ind = 0; ind < 3; ++ind)
		    basis_gradient_actual[p][transform_local[ind_index][ind_tl]][ind] += weight_transform_local[ind_index][ind_tl] * basis_gradient[p][ind_index][ind];
    std::cerr << "assigned basis_gradient_actual\n";
}

TEMPLATE_TSEM
void THIS_TSEM::build_geometry_info(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    // setup, basic variables
    point_ref_mesh.resize(4); // [i = 0:3], barycenter of i dimensional geometry
    for (int ind = 0; ind <= 3; ++ind)
        n_geometry[ind] = mesh.n_geometry(ind);
    n_element = n_geometry[3];
    std::cerr << "read number of geometries, n_geometry[0:3] = [" << n_geometry[0] << ", " << n_geometry[1] << ", " << n_geometry[2] << ", " << n_geometry[3] << "]\n";
    std::cerr << "n_element = " << n_element << '\n';
    n_geometry_total = 0;
    for (int ind = 0; ind <= 3; ++ind)
        n_geometry_total += n_geometry[ind];
    for (int ind = 0; ind <= 3; ++ind) // number of dof on ind dimensional geometry
        n_dof_geometry[ind] = calc_binomial_coefficient(M-1, ind);
    for (int ind = 0; ind <= 3; ++ind)
        std::cerr << "n_dof_geometry[" << ind << "] = " << n_dof_geometry[ind] << ",\t";
    std::cerr << '\n';
    n_dof_total = 0; // number of total degree of freedom
    for (int ind = 0; ind <= 3; ++ind)
        n_dof_total += n_geometry[ind] * n_dof_geometry[ind]; // sum up all dof on each dimensional geometry
    std::cerr << "n_dof_total = " << n_dof_total << '\n';
    // assign point_ref
    for (int ind = 0; ind <= 3; ++ind){
        point_ref_mesh[ind].resize(mesh.n_geometry(ind));
        for (int i = 0; i < mesh.n_geometry(ind); ++i){
            AFEPack::Point<3> point_tmp = mesh.point(mesh.geometry(ind, i).vertex(0));
            for (int indt = 1; indt <= ind; ++indt)
                point_tmp += mesh.point(mesh.geometry(ind, i).vertex(indt));
            point_tmp /= (ind + 1);
            point_ref_mesh[ind][i] = point_tmp;
        }
    }


    // record order of nodes on each dimensional geometry
    // number_node[ind+1=1:2][order+1=1:ind+1][]: the node number of ind+1-dimensional geometry
    number_node.resize(2);
    number_edge.resize(n_geometry[2]);
    for (int ind = 1; ind <= 2; ++ind){
        number_node[ind-1].resize(n_geometry[ind]);
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo)
            number_node[ind-1][ind_geo].resize(ind+1, -1);
    }
    for (int ind_geo = 0; ind_geo < n_geometry[2]; ++ind_geo)
        number_edge[ind_geo].resize(3, -1);
    // flag that whether order of start/end point on global face is the same as that on global edge
    flag_sameorder_edgeonface.resize(n_geometry[2]);
    for (int ind_face = 0; ind_face < n_geometry[2]; ++ind_face)
	flag_sameorder_edgeonface[ind_face].resize(3, true);
    // type of projection for edge and face on each element
    type_projection.resize(n_element);
    type_projection_inv.resize(n_element); // inverse projection type, for face on each element
    // traverse finite element space, record number_node
    // the order of node number determine others, i.e., correspond to nodes on reference as P0, P1, P2, P3
    std::vector<bool> flag_assign_edge(n_geometry[1], false);
    std::vector<bool> flag_assign_face(n_geometry[2], false);
    
    index_geometry_onelement.resize(n_element);
    int index_start_point[] = {0, 0, 0, 1, 1, 2};
    int index_end_point[]   = {1, 2, 3, 2, 3, 3};
    int index_face_point[] = {1, 2, 3,
                              0, 2, 3,
                              0, 1, 3,
                              0, 1, 2};
    int index_face_edge[] = {5, 4, 3,
                             5, 2, 1,
                             4, 2, 0,
                             3, 1, 0};
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // initialize
        index_geometry_onelement[ind_ele].resize(3);
        index_geometry_onelement[ind_ele][0].resize(4); // index for nodes
        index_geometry_onelement[ind_ele][1].resize(6); // index for edges
        index_geometry_onelement[ind_ele][2].resize(4); // index for faces
        type_projection[ind_ele].resize(2);
        type_projection[ind_ele][0].resize(6, -1); // 6 edge, 0-the same order; 1-inverse order
        type_projection[ind_ele][1].resize(4, -1); // 4 face, 0-the same order; 1-5-corresponding variation
        type_projection_inv[ind_ele].resize(4, -1); // 4 face, 0-the same order; 1-5-corresponding variation
        // read node number of this element, correspond to P0, P1, P2, P3
        // the node number determine the order of edge and face on this element
        for (int ind_node = 0; ind_node < 4; ++ind_node)
            index_geometry_onelement[ind_ele][0][ind_node] = mesh.geometry(3, ind_ele).vertex(ind_node);

        // read edge number of this element, correspond to E01, E02, E03, E12, E13, E23
        int ind_face_this;
        // consider face in front of P3, which is P0P1P2
        ind_face_this = mesh.geometry(3, ind_ele).boundary(3);
        for (int ind_e = 0; ind_e < 3; ++ind_e){
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][0]) // in front of P0
                index_geometry_onelement[ind_ele][1][3] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E12
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][1]) // in front of P1
                index_geometry_onelement[ind_ele][1][1] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E02
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][2]) // in front of P2
                index_geometry_onelement[ind_ele][1][0] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E01
        }
        // consider face in front of P2, which is P0P1P3
        ind_face_this = mesh.geometry(3, ind_ele).boundary(2);
        for (int ind_e = 0; ind_e < 3; ++ind_e){
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][0]) // in front of P0
                index_geometry_onelement[ind_ele][1][4] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E13
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][1]) // in front of P1
                index_geometry_onelement[ind_ele][1][2] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E03
        }
        // consider face in front of P1, which is P0P2P3
        ind_face_this = mesh.geometry(3, ind_ele).boundary(1);
        for (int ind_e = 0; ind_e < 3; ++ind_e)
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][0]) // in front of P0
                index_geometry_onelement[ind_ele][1][5] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E23

        // read face number of this element, correspond to F0=P1P2P3, F1=P0P2P3, F2=P0P1P3, F3=P0P1P2
        for (int ind_face = 0; ind_face < 4; ++ind_face)
            index_geometry_onelement[ind_ele][2][ind_face] = mesh.geometry(3, ind_ele).boundary(ind_face);
       
        // traverse 6 edges of this element
        for (int ind_edge = 0; ind_edge < 6; ++ind_edge){
            int &index_e = index_geometry_onelement[ind_ele][1][ind_edge];
            if (!flag_assign_edge[index_e]){ // if this edge hasn't been traversed
                // record start point and end point
                number_node[0][index_e][0] = index_geometry_onelement[ind_ele][0][index_start_point[ind_edge]];
                number_node[0][index_e][1] = index_geometry_onelement[ind_ele][0][index_end_point[ind_edge]];
                // record projection type
                type_projection[ind_ele][0][ind_edge] = 0;
                // update flag
                flag_assign_edge[index_e] = true;
            }
            else{
                // record projection type
                if (number_node[0][index_e][0] == index_geometry_onelement[ind_ele][0][index_start_point[ind_edge]]) // the same start point
                    type_projection[ind_ele][0][ind_edge] = 0;
                else
                    type_projection[ind_ele][0][ind_edge] = 1;
            }
        }
        
        // traveres 4 faces of this element
        for (int ind_face = 0; ind_face < 4; ++ind_face){
            int &index_f = index_geometry_onelement[ind_ele][2][ind_face];
            if (!flag_assign_face[index_f]){ // if this face hasn't been traversed
                // record order of nodes
                number_node[1][index_f][0] = index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 0]];
                number_node[1][index_f][1] = index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 1]];
                number_node[1][index_f][2] = index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 2]];
                number_edge[index_f][0] = index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 0]];
                number_edge[index_f][1] = index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 1]];
                number_edge[index_f][2] = index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 2]];
                // record whether order of endpoints of edge the same as those on global edge, denote face as 012
                if (number_node[0][index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 0]]][0]
                    != index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 1]]) // startpoint of 12 != 1
                    flag_sameorder_edgeonface[index_f][0] = false;
                if (number_node[0][index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 1]]][0]
                    != index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 0]]) // startpoint of 02 != 0
                    flag_sameorder_edgeonface[index_f][1] = false;
                if (number_node[0][index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 2]]][0]
                    != index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 0]]) // startpoint of 01 != 0
                    flag_sameorder_edgeonface[index_f][2] = false;
                // record projection type
                type_projection[ind_ele][1][ind_face] = 0;
                type_projection_inv[ind_ele][ind_face] = 0;
                // update flag
                flag_assign_face[index_f] = true;
            }
            else{
                // record projection type: 012, 021, 102, 120, 201, 210
                std::vector<int> order_node_on_ref(3, -1); // order of nodes on reference element
                for (int ind_n = 0; ind_n < 3; ++ind_n)
                    for (int ind_n_r = 0; ind_n_r < 3; ++ind_n_r)
                        if (index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + ind_n]] == number_node[1][index_f][ind_n_r])
                            order_node_on_ref[ind_n] = ind_n_r;
                int flag_number = order_node_on_ref[0] * 100 + order_node_on_ref[1] * 10 + order_node_on_ref[2];
                switch (flag_number){
                case 12: // 0 1 2, 0 1 2
                    type_projection[ind_ele][1][ind_face] = 0;
                    type_projection_inv[ind_ele][ind_face] = 0;
                    break;
                case 21: // 0 2 1, 0 2 1
                    type_projection[ind_ele][1][ind_face] = 1;
                    type_projection_inv[ind_ele][ind_face] = 1;
                    break;
                case 102: // 1 0 2, 1 0 2
                    type_projection[ind_ele][1][ind_face] = 2;
                    type_projection_inv[ind_ele][ind_face] = 2;
                    break;
                case 120: // 1 2 0, 2 0 1
                    type_projection[ind_ele][1][ind_face] = 3;
                    type_projection_inv[ind_ele][ind_face] = 4;
                    break;
                case 201: // 2 0 1, 1 2 0
                    type_projection[ind_ele][1][ind_face] = 4;
                    type_projection_inv[ind_ele][ind_face] = 3;
                    break;
                case 210: // 2 1 0, 2 1 0
                    type_projection[ind_ele][1][ind_face] = 5;
                    type_projection_inv[ind_ele][ind_face] = 5;
                    break;
                default:
                    std::cout << "error: assign type of projection!\n";
                }
            }
        }
    }
    std::cerr << "assign type_projection\n";

    // calculate volume for each element
    val_volume.resize(fem_space.n_element());
    AFEPack::Point<3> point_tmp;
    point_tmp[0] = point_tmp[1] = point_tmp[2] = 1.0/3.0;
    for (unsigned int ind_ele = 0; ind_ele < fem_space.n_element(); ++ind_ele){
    	valuetype volume = fem_space.element(ind_ele).templateElement().volume();
    	valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed
    	val_volume[ind_ele] = fabs(volume * jacobian);
    }
    // val_volume.resize(mesh.n_geometry(3));
    // for (unsigned int ind_ele = 0; ind_ele < mesh.n_geometry(3); ++ind_ele){
    // 	AFEPack::Point<3> &p0 = mesh.point(mesh.geometry(3, ind_ele).vertex(0));
    // 	AFEPack::Point<3> &p1 = mesh.point(mesh.geometry(3, ind_ele).vertex(1));
    // 	AFEPack::Point<3> &p2 = mesh.point(mesh.geometry(3, ind_ele).vertex(2));
    // 	AFEPack::Point<3> &p3 = mesh.point(mesh.geometry(3, ind_ele).vertex(3));
    // 	val_volume[ind_ele] = calc_volume_tetrahedron(p0, p1, p2, p3);
    // }
}

TEMPLATE_TSEM
void THIS_TSEM::build_geometry_info(Mesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    // setup, basic variables
    point_ref_mesh.resize(4); // [i = 0:3], barycenter of i dimensional geometry
    for (int ind = 0; ind <= 3; ++ind)
        n_geometry[ind] = mesh.n_geometry(ind);
    n_element = n_geometry[3];
    std::cerr << "read number of geometries, n_geometry[0:3] = [" << n_geometry[0] << ", " << n_geometry[1] << ", " << n_geometry[2] << ", " << n_geometry[3] << "]\n";
    std::cerr << "n_element = " << n_element << '\n';
    n_geometry_total = 0;
    for (int ind = 0; ind <= 3; ++ind)
        n_geometry_total += n_geometry[ind];
    for (int ind = 0; ind <= 3; ++ind) // number of dof on ind dimensional geometry
        n_dof_geometry[ind] = calc_binomial_coefficient(M-1, ind);
    for (int ind = 0; ind <= 3; ++ind)
        std::cerr << "n_dof_geometry[" << ind << "] = " << n_dof_geometry[ind] << ",\t";
    std::cerr << '\n';
    n_dof_total = 0; // number of total degree of freedom
    for (int ind = 0; ind <= 3; ++ind)
        n_dof_total += n_geometry[ind] * n_dof_geometry[ind]; // sum up all dof on each dimensional geometry
    std::cerr << "n_dof_total = " << n_dof_total << '\n';
    // assign point_ref
    for (int ind = 0; ind <= 3; ++ind){
        point_ref_mesh[ind].resize(mesh.n_geometry(ind));
        for (int i = 0; i < mesh.n_geometry(ind); ++i){
            AFEPack::Point<3> point_tmp = mesh.point(mesh.geometry(ind, i).vertex(0));
            for (int indt = 1; indt <= ind; ++indt)
                point_tmp += mesh.point(mesh.geometry(ind, i).vertex(indt));
            point_tmp /= (ind + 1);
            point_ref_mesh[ind][i] = point_tmp;
        }
    }


    // record order of nodes on each dimensional geometry
    // number_node[ind+1=1:2][order+1=1:ind+1][]: the node number of ind+1-dimensional geometry
    number_node.resize(2);
    number_edge.resize(n_geometry[2]);
    for (int ind = 1; ind <= 2; ++ind){
        number_node[ind-1].resize(n_geometry[ind]);
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo)
            number_node[ind-1][ind_geo].resize(ind+1, -1);
    }
    for (int ind_geo = 0; ind_geo < n_geometry[2]; ++ind_geo)
        number_edge[ind_geo].resize(3, -1);
    // flag that whether order of start/end point on global face is the same as that on global edge
    flag_sameorder_edgeonface.resize(n_geometry[2]);
    for (int ind_face = 0; ind_face < n_geometry[2]; ++ind_face)
	flag_sameorder_edgeonface[ind_face].resize(3, true);
    // type of projection for edge and face on each element
    type_projection.resize(n_element);
    type_projection_inv.resize(n_element); // inverse projection type, for face on each element
    // traverse finite element space, record number_node
    // the order of node number determine others, i.e., correspond to nodes on reference as P0, P1, P2, P3
    std::vector<bool> flag_assign_edge(n_geometry[1], false);
    std::vector<bool> flag_assign_face(n_geometry[2], false);
    
    index_geometry_onelement.resize(n_element);
    int index_start_point[] = {0, 0, 0, 1, 1, 2};
    int index_end_point[]   = {1, 2, 3, 2, 3, 3};
    int index_face_point[] = {1, 2, 3,
                              0, 2, 3,
                              0, 1, 3,
                              0, 1, 2};
    int index_face_edge[] = {5, 4, 3,
                             5, 2, 1,
                             4, 2, 0,
                             3, 1, 0};
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // initialize
        index_geometry_onelement[ind_ele].resize(3);
        index_geometry_onelement[ind_ele][0].resize(4); // index for nodes
        index_geometry_onelement[ind_ele][1].resize(6); // index for edges
        index_geometry_onelement[ind_ele][2].resize(4); // index for faces
        type_projection[ind_ele].resize(2);
        type_projection[ind_ele][0].resize(6, -1); // 6 edge, 0-the same order; 1-inverse order
        type_projection[ind_ele][1].resize(4, -1); // 4 face, 0-the same order; 1-5-corresponding variation
        type_projection_inv[ind_ele].resize(4, -1); // 4 face, 0-the same order; 1-5-corresponding variation
        // read node number of this element, correspond to P0, P1, P2, P3
        // the node number determine the order of edge and face on this element
        for (int ind_node = 0; ind_node < 4; ++ind_node)
            index_geometry_onelement[ind_ele][0][ind_node] = mesh.geometry(3, ind_ele).vertex(ind_node);

        // read edge number of this element, correspond to E01, E02, E03, E12, E13, E23
        int ind_face_this;
        // consider face in front of P3, which is P0P1P2
        ind_face_this = mesh.geometry(3, ind_ele).boundary(3);
        for (int ind_e = 0; ind_e < 3; ++ind_e){
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][0]) // in front of P0
                index_geometry_onelement[ind_ele][1][3] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E12
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][1]) // in front of P1
                index_geometry_onelement[ind_ele][1][1] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E02
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][2]) // in front of P2
                index_geometry_onelement[ind_ele][1][0] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E01
        }
        // consider face in front of P2, which is P0P1P3
        ind_face_this = mesh.geometry(3, ind_ele).boundary(2);
        for (int ind_e = 0; ind_e < 3; ++ind_e){
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][0]) // in front of P0
                index_geometry_onelement[ind_ele][1][4] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E13
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][1]) // in front of P1
                index_geometry_onelement[ind_ele][1][2] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E03
        }
        // consider face in front of P1, which is P0P2P3
        ind_face_this = mesh.geometry(3, ind_ele).boundary(1);
        for (int ind_e = 0; ind_e < 3; ++ind_e)
            if (mesh.geometry(2, ind_face_this).vertex(ind_e) == index_geometry_onelement[ind_ele][0][0]) // in front of P0
                index_geometry_onelement[ind_ele][1][5] = mesh.geometry(2, ind_face_this).boundary(ind_e); // correspond to E23

        // read face number of this element, correspond to F0=P1P2P3, F1=P0P2P3, F2=P0P1P3, F3=P0P1P2
        for (int ind_face = 0; ind_face < 4; ++ind_face)
            index_geometry_onelement[ind_ele][2][ind_face] = mesh.geometry(3, ind_ele).boundary(ind_face);
       
        // traverse 6 edges of this element
        for (int ind_edge = 0; ind_edge < 6; ++ind_edge){
            int &index_e = index_geometry_onelement[ind_ele][1][ind_edge];
            if (!flag_assign_edge[index_e]){ // if this edge hasn't been traversed
                // record start point and end point
                number_node[0][index_e][0] = index_geometry_onelement[ind_ele][0][index_start_point[ind_edge]];
                number_node[0][index_e][1] = index_geometry_onelement[ind_ele][0][index_end_point[ind_edge]];
                // record projection type
                type_projection[ind_ele][0][ind_edge] = 0;
                // update flag
                flag_assign_edge[index_e] = true;
            }
            else{
                // record projection type
                if (number_node[0][index_e][0] == index_geometry_onelement[ind_ele][0][index_start_point[ind_edge]]) // the same start point
                    type_projection[ind_ele][0][ind_edge] = 0;
                else
                    type_projection[ind_ele][0][ind_edge] = 1;
            }
        }
        
        // traveres 4 faces of this element
        for (int ind_face = 0; ind_face < 4; ++ind_face){
            int &index_f = index_geometry_onelement[ind_ele][2][ind_face];
            if (!flag_assign_face[index_f]){ // if this face hasn't been traversed
                // record order of nodes
                number_node[1][index_f][0] = index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 0]];
                number_node[1][index_f][1] = index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 1]];
                number_node[1][index_f][2] = index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 2]];
                number_edge[index_f][0] = index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 0]];
                number_edge[index_f][1] = index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 1]];
                number_edge[index_f][2] = index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 2]];
                // record whether order of endpoints of edge the same as those on global edge, denote face as 012
                if (number_node[0][index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 0]]][0]
                    != index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 1]]) // startpoint of 12 != 1
                    flag_sameorder_edgeonface[index_f][0] = false;
                if (number_node[0][index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 1]]][0]
                    != index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 0]]) // startpoint of 02 != 0
                    flag_sameorder_edgeonface[index_f][1] = false;
                if (number_node[0][index_geometry_onelement[ind_ele][1][index_face_edge[3*ind_face + 2]]][0]
                    != index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + 0]]) // startpoint of 01 != 0
                    flag_sameorder_edgeonface[index_f][2] = false;
                // record projection type
                type_projection[ind_ele][1][ind_face] = 0;
                type_projection_inv[ind_ele][ind_face] = 0;
                // update flag
                flag_assign_face[index_f] = true;
            }
            else{
                // record projection type: 012, 021, 102, 120, 201, 210
                std::vector<int> order_node_on_ref(3, -1); // order of nodes on reference element
                for (int ind_n = 0; ind_n < 3; ++ind_n)
                    for (int ind_n_r = 0; ind_n_r < 3; ++ind_n_r)
                        if (index_geometry_onelement[ind_ele][0][index_face_point[3*ind_face + ind_n]] == number_node[1][index_f][ind_n_r])
                            order_node_on_ref[ind_n] = ind_n_r;
                int flag_number = order_node_on_ref[0] * 100 + order_node_on_ref[1] * 10 + order_node_on_ref[2];
                switch (flag_number){
                case 12: // 0 1 2, 0 1 2
                    type_projection[ind_ele][1][ind_face] = 0;
                    type_projection_inv[ind_ele][ind_face] = 0;
                    break;
                case 21: // 0 2 1, 0 2 1
                    type_projection[ind_ele][1][ind_face] = 1;
                    type_projection_inv[ind_ele][ind_face] = 1;
                    break;
                case 102: // 1 0 2, 1 0 2
                    type_projection[ind_ele][1][ind_face] = 2;
                    type_projection_inv[ind_ele][ind_face] = 2;
                    break;
                case 120: // 1 2 0, 2 0 1
                    type_projection[ind_ele][1][ind_face] = 3;
                    type_projection_inv[ind_ele][ind_face] = 4;
                    break;
                case 201: // 2 0 1, 1 2 0
                    type_projection[ind_ele][1][ind_face] = 4;
                    type_projection_inv[ind_ele][ind_face] = 3;
                    break;
                case 210: // 2 1 0, 2 1 0
                    type_projection[ind_ele][1][ind_face] = 5;
                    type_projection_inv[ind_ele][ind_face] = 5;
                    break;
                default:
                    std::cout << "error: assign type of projection!\n";
                }
            }
        }
    }
    std::cerr << "assign type_projection\n";

    // calculate volume for each element
    val_volume.resize(fem_space.n_element());
    AFEPack::Point<3> point_tmp;
    point_tmp[0] = point_tmp[1] = point_tmp[2] = 1.0/3.0;
    for (unsigned int ind_ele = 0; ind_ele < fem_space.n_element(); ++ind_ele){
    	valuetype volume = fem_space.element(ind_ele).templateElement().volume();
    	valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed
    	val_volume[ind_ele] = fabs(volume * jacobian);
    }
}

TEMPLATE_TSEM
void THIS_TSEM::build_global_transform(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    // assign the weight_location, geometry_dimension and geometry_order of dof, which determine the order of these geometry in discretized matrix
    /* {weight_local, geometry_dimension, geometry_order}:
     *    weight_local:       a number for sorting, read from fem_space, determine the order of all dimensional geomtry
     *    geometry_dimension: the dimension of geometry corresponds to this set
     *    geometry_order:     the order of geometry in the same dimensional ones corresponds to this set
     */
    weight_location.resize(n_geometry_total); // the weight_location of 0: 3 dimensional geometry in turns
    geometry_dimension.resize(n_geometry_total);
    geometry_order.resize(n_geometry_total); // [i = 0:n_geometry_total-1] the order of i-th entry in weight_location, whose dimension is geometry_dimension[i]
    // correspondence between fem dof/element and geometry
    // transform_femdof2geometry: index from fem dof to 0 & 1 dimensional geometry, for 1 dimensional geometry, its index plus mesh.n_geometry(0)
    int n_dof = fem_space.n_dof();
    transform_femdof2geometry.resize(n_dof, -1);
    // location_geometry: location of all geometry (0:3 dimensional) according to increasing order of weight_location
    location_geometry.resize(n_geometry_total);
    // location_actualdof: start index of geometry in actual discretized matrix
    location_actualdof.resize(n_geometry_total);
    /* geometry in finite element space                       -> total index of all geometries     -> position in actual discretized system
     * transform_femdof2geometry (0 & 1 dimensional geometry) -> total index = 0: n_geometry_total -> location_actualdof[location_geometry[total index]]
     * transform_femele2geometry (2 & 3 dimensional geometry)
     */
    // assign geometry_dimension and geometry_order
    for (int ind = 0; ind <= 3; ++ind){
        int index_start = 0;
        for (int indt = 0; indt < ind; ++indt)
            index_start += n_geometry[indt];
        for (int i = 0; i < n_geometry[ind]; ++i){
            geometry_dimension[index_start + i] = ind;
            geometry_order[index_start + i] = i;
        }
    }
    // construct correspondence between fem dof and geometry (I): find weight_location for 0 & 1 dimensional geometry according to the order in fem_space
    for (int i = 0; i < n_dof; ++i)
        for (int ind = 0; ind <= 1 && transform_femdof2geometry[i] < 0; ++ind) // as all order in fem_space is natural number, so >= 0
            for (int j = 0; j < n_geometry[ind]; ++j)
                if (distance(point_ref_mesh[ind][j], fem_space.dofInfo(i).interp_point) < tol_zero){
                    transform_femdof2geometry[i] = ind * n_geometry[0] + j; // if ind == 1, add n_geometry[0], total index of point or edge
                    weight_location[transform_femdof2geometry[i]] = i;
                    break;
                }
    // calculate weight_location for 2 & 3 dimensional geometry
    for (int ind = 2; ind <= 3; ++ind){
        int index_start = n_geometry[0] + n_geometry[1] + (ind-2) * n_geometry[2]; // if ind == 3, add n_geometry[2]
        for (int i = 0; i < n_geometry[ind]; ++i){
            valuetype weight_tmp = 0.0;
            for (int indt = 0; indt <= ind; ++indt)
                weight_tmp += weight_location[mesh.geometry(ind, i).vertex(indt)];
            weight_tmp /= ind + 1;
            weight_location[index_start + i] = weight_tmp;
        }
    }
    // sort according to weight_location
    for (int i = 1; i < n_geometry_total; ++i)
        for (int j = i; j > 0; --j)
            if (weight_location[j] < weight_location[j-1]){
                valuetype tmp = weight_location[j];
                weight_location[j] = weight_location[j-1];
                weight_location[j-1] = tmp;
                int tmpi = geometry_dimension[j];
                geometry_dimension[j] = geometry_dimension[j-1];
                geometry_dimension[j-1] = tmpi;
                tmpi = geometry_order[j];
                geometry_order[j] = geometry_order[j-1];
                geometry_order[j-1] = tmpi;
            }
    for (int i = 0; i < n_geometry_total; ++i){ // assign location for geometry_order[i]-th geometry_dimension[i] dimensional geometry
        int index = geometry_order[i]; // recover the index of geometry in whole order
        for (int ind = 0; ind < geometry_dimension[i]; ++ind)
            index += n_geometry[ind];
        location_geometry[index] = i;
    }
    // insert dof into each location_geometry
    location_actualdof[0] = 0;
    for (int i = 1; i < n_geometry_total; ++i)
        location_actualdof[i] = location_actualdof[i-1] + n_dof_geometry[geometry_dimension[i-1]];
    std::cerr << "assign location_actualdof, location_geometry\n";
    
    
    // traverse element, construct transformation between local and global
    // initialize
    // assign conversion between lexicography order and the order in jia2022
    conversion.resize(2);
    conversion[0].resize(n_dof_face);
    conversion[1].resize(n_dof_geometry[3]);
    for (unsigned int l1 = 2; l1 <= M; ++l1)
	for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
	    unsigned int ind_lex = (l1+l2+1-3) * (l1+l2-3) / 2 + l2-1; // location of (l1-2, l2-1) in lexicography order
	    // unsigned int ind_jia = (M+1)*(l1-2) - (l1-1-2)*(l1-2)/2 + l2-1;
	    // conversion[0][ind_lex] = ind_jia;
	    conversion[0][ind_lex] = ind_lex;
	}
    for (unsigned int ind_index = 0; ind_index < n_dof_geometry[3]; ++ind_index){
	Multiindex<3> index_now = correspondence.number2index(ind_index);
	unsigned int l1 = index_now.index[0], l2 = index_now.index[1], l3 = index_now.index[2];
	// unsigned int ind_jia = ((l1-1)*l1*(2*l1-1) + l1*(M+1)*(M+2)*6 - (2*M+3)*(l1-1)*l1*3) / 12
	//     + (M-l1+1)*l2 - (l2-1)*l2/2
	//     + l3;
	// conversion[1][ind_index] = ind_jia;
	conversion[1][ind_index] = ind_index;
    }
    std::cerr << "assign conversion\n";
    // expression of local coefficient by linear summation of global ones
    transform_n_global2local.resize(n_element);
    transform_ind_global2local.resize(n_element);
    transform_val_global2local.resize(n_element);
    // expression of global coefficient on each-dimensional geometry by local ones, equivalent to column compress of transformation matrix
    transform_n_local2global.resize(n_element);
    transform_ind_local2global.resize(n_element);
    transform_val_local2global.resize(n_element);
    // assign transform_local2global and transform_global2local
    std::vector<int> index_edgedof_local(n_dof_edge * 6); // index of multiindex corresponds to edge dof on local element
    for (int l = 2; l <= M; ++l){
        index_edgedof_local[n_dof_edge*0 + l-2] = correspondence.index2number(Unitary_Multiindex[0]*l); // (l, 0, 0)
        index_edgedof_local[n_dof_edge*1 + l-2] = correspondence.index2number(Unitary_Multiindex[1]*l); // (0, l, 0)
        index_edgedof_local[n_dof_edge*2 + l-2] = correspondence.index2number(Unitary_Multiindex[2]*l); // (0, 0, l)
        index_edgedof_local[n_dof_edge*4 + l-2] = correspondence.index2number(Unitary_Multiindex[0] + Unitary_Multiindex[2]*(l-1)); // (1,   0, l-1)
        index_edgedof_local[n_dof_edge*3 + l-2] = correspondence.index2number(Unitary_Multiindex[0] + Unitary_Multiindex[1]*(l-1)); // (1, l-1,   0)
        index_edgedof_local[n_dof_edge*5 + l-2] = correspondence.index2number(Unitary_Multiindex[1] + Unitary_Multiindex[2]*(l-1)); // (0,   1, l-1)
    }
    std::vector<int> index_facedof_local(n_dof_face * 4); // index of multiindex corresponds to face dof on local element
    for (int l1 = 2; l1 <= M; ++l1)
        for (int l2 = 1; l2 <= M-l1; ++l2){
            int ind_index = (l1+l2-3)*(l1+l2-2)/2 + l2-1;
            index_facedof_local[n_dof_face*0 + ind_index] =
                correspondence.index2number(Unitary_Multiindex[0] + Unitary_Multiindex[1]*(l1-1) + Unitary_Multiindex[2]*l2); // (1, l1-1, l2)
            index_facedof_local[n_dof_face*1 + ind_index] =
                correspondence.index2number(Unitary_Multiindex[1]*l1 + Unitary_Multiindex[2]*l2); // ( 0, l1, l2)
            index_facedof_local[n_dof_face*2 + ind_index] =
                correspondence.index2number(Unitary_Multiindex[0]*l1 + Unitary_Multiindex[2]*l2); // (l1,  0, l2)
            index_facedof_local[n_dof_face*3 + ind_index] =
                correspondence.index2number(Unitary_Multiindex[0]*l1 + Unitary_Multiindex[1]*l2); // (l1, l2,  0)
        }
    int index_edge_local[] = {5, 4, 3, // edge number for each face on local element
			      5, 2, 1,
			      4, 2, 0,
			      3, 1, 0};
    int order_edge_global[] = {0, 1, 2, // order of edge on global face in each variation
                               0, 2, 1,
                               1, 0, 2,
                               1, 2, 0,
                               2, 0, 1,
                               2, 1, 0};
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // initialize
        transform_n_global2local[ind_ele].resize(n_index);
        transform_ind_global2local[ind_ele].resize(n_index);
        transform_val_global2local[ind_ele].resize(n_index);
        transform_n_local2global[ind_ele].resize(n_index, 0);
        transform_ind_local2global[ind_ele].resize(n_index);
        transform_val_local2global[ind_ele].resize(n_index);
        
        // construct transformation node dof
        // P0, P1, P2, P3 excatly correspond first 4 multiindex (0,0,0), (1,0,0), (0,1,0), (0,0,1)
        for (int ind_node = 0; ind_node < 4; ++ind_node){ // traverse node: P0, P1, P2, P3
            int ind_index = ind_node;
            int ind_geo = index_geometry_onelement[ind_ele][0][ind_node];
            int ind_dof_global = location_actualdof[location_geometry[ind_geo]];
            int n_dep = 1; // number of dependency
            transform_n_global2local[ind_ele][ind_index] = n_dep;
            transform_ind_global2local[ind_ele][ind_index].resize(n_dep);
            transform_val_global2local[ind_ele][ind_index].resize(n_dep);
            transform_ind_global2local[ind_ele][ind_index][0] = ind_dof_global;
            transform_val_global2local[ind_ele][ind_index][0] = 1;
        }
        
        // construct transformation for edge dof
        for (int ind_edge = 0; ind_edge < 6; ++ind_edge){ // traverse edge: E01, E02, E03, E12, E13, E23
            int ind_geo = index_geometry_onelement[ind_ele][1][ind_edge];
            int ind_dof_start = location_actualdof[location_geometry[ind_geo + n_geometry[0]]];
            for (int l = 2; l <= M; ++l){
                int ind_index = index_edgedof_local[n_dof_edge*ind_edge + l-2];
                int ind_dof_global = ind_dof_start + l-2;
                int n_dep = 1;
                transform_n_global2local[ind_ele][ind_index] = n_dep;
                transform_ind_global2local[ind_ele][ind_index].resize(n_dep);
                transform_val_global2local[ind_ele][ind_index].resize(n_dep);
                transform_ind_global2local[ind_ele][ind_index][0] = ind_dof_global;
                transform_val_global2local[ind_ele][ind_index][0] = 1;
                if (type_projection[ind_ele][0][ind_edge] == 1 && l%2 == 1)
                    transform_val_global2local[ind_ele][ind_index][0] = -1;
            }
        }
        
        // construct transformation for face dof
        for (int ind_face = 0; ind_face < 4; ++ind_face){
            int ind_geo = index_geometry_onelement[ind_ele][2][ind_face];
            int ind_facedof_start = location_actualdof[location_geometry[ind_geo + n_geometry[0] + n_geometry[1]]];
            for (int ind_dof = 0; ind_dof < n_dof_geometry[2]; ++ind_dof){
                int ind_index = index_facedof_local[n_dof_face*ind_face + ind_dof];
                int ind_var = type_projection_inv[ind_ele][ind_face];
                int n_dep = 0;
                for (int ind_e = 0; ind_e < 3; ++ind_e)
                    n_dep += n_expr_edge[ind_var][ind_dof][ind_e];
                n_dep += n_expr_face[ind_var][ind_dof];
                transform_n_global2local[ind_ele][ind_index] = n_dep;
                transform_ind_global2local[ind_ele][ind_index].resize(n_dep);
                transform_val_global2local[ind_ele][ind_index].resize(n_dep);
                // assign edge contribution
                int pos_start = 0;
                for (int ind_e = 0; ind_e < 3; ++ind_e){
                    int ind_edge = number_edge[index_geometry_onelement[ind_ele][2][ind_face]][ind_e];
                    int ind_edgedof_start = location_actualdof[location_geometry[ind_edge + n_geometry[0]]];
                    for (int ind_nnz = 0; ind_nnz < n_expr_edge[ind_var][ind_dof][ind_e]; ++ind_nnz){
                        transform_ind_global2local[ind_ele][ind_index][pos_start + ind_nnz]
                            = ind_edgedof_start + ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        transform_val_global2local[ind_ele][ind_index][pos_start + ind_nnz]
                            = val_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        if (!flag_sameorder_edgeonface[index_geometry_onelement[ind_ele][2][ind_face]][ind_e]
                            && ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz] % 2 == 1)
                            transform_val_global2local[ind_ele][ind_index][pos_start + ind_nnz] *= -1;
                    }
                    pos_start += n_expr_edge[ind_var][ind_dof][ind_e];
                }
                // assign face contribution
                for (int ind_nnz = 0; ind_nnz < n_expr_face[ind_var][ind_dof]; ++ind_nnz){
                    transform_ind_global2local[ind_ele][ind_index][pos_start + ind_nnz]
                        = ind_facedof_start + conversion[0][ind_expr_face[ind_var][ind_dof][ind_nnz]];
                    transform_val_global2local[ind_ele][ind_index][pos_start + ind_nnz]
                        = val_expr_face[ind_var][ind_dof][ind_nnz];
                }
                if (pos_start + n_expr_face[ind_var][ind_dof] != n_dep)
                    std::cerr << "error: assign transform_global2local for face dof\n";
            }
        }
        
        // construct transformation for interior dof
        int ind_interiordof_start = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]];
        for (int ind_dof = 0; ind_dof < n_dof_geometry[3]; ++ind_dof){
            Multiindex<3> index_now = correspondence.number2index(ind_dof);
            index_now = index_now + Unitary_Multiindex[0]*2 + Unitary_Multiindex[1] + Unitary_Multiindex[2];
            int ind_index = correspondence.index2number(index_now);
            int n_dep = 1;
            transform_n_global2local[ind_ele][ind_index] = n_dep;
            transform_ind_global2local[ind_ele][ind_index].resize(n_dep);
            transform_val_global2local[ind_ele][ind_index].resize(n_dep);
            transform_ind_global2local[ind_ele][ind_index][0] = ind_interiordof_start + conversion[1][ind_dof];
            transform_val_global2local[ind_ele][ind_index][0] = 1;
        }

        // traverse all dimensional geometries, count contribution for transform_local2global
        for (int ind_node = 0; ind_node < 4; ++ind_node) // ind_node -> ind_index -> 1-1 global vertex coefficient
            transform_n_local2global[ind_ele][ind_node]++;
        for (int ind_edge = 0; ind_edge < 6; ++ind_edge) // ind_edge, ind_dof -> ind_index -> 1-1 global edge coefficient
            for (int ind_dof = 0; ind_dof < n_dof_edge; ++ind_dof)
                transform_n_local2global[ind_ele][index_edgedof_local[n_dof_edge*ind_edge + ind_dof]]++;
        for (int ind_face = 0; ind_face < 4; ++ind_face){
            int& ind_var = type_projection[ind_ele][1][ind_face];
            for (int ind_dof = 0; ind_dof < n_dof_face; ++ind_dof){
                // count edge contribution
                for (int ind_e = 0; ind_e < 3; ++ind_e){ // ind_dof, ind_e -> edge index -> corresponding multiindex
                    int& edge_index_local = index_edge_local[3*ind_face + ind_e]; // local index of this edge, between 0 and 5
                    for (int ind_nnz = 0; ind_nnz < n_expr_edge[ind_var][ind_dof][ind_e]; ++ind_nnz){
                        int& edgedof_index = ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        int& multiindex_index = index_edgedof_local[n_dof_edge*edge_index_local + edgedof_index];
                        transform_n_local2global[ind_ele][multiindex_index]++;
                    }
                }
                // count face contribution
                for (int ind_nnz = 0; ind_nnz < n_expr_face[ind_var][ind_dof]; ++ind_nnz){
                    int& facedof_index = ind_expr_face[ind_var][ind_dof][ind_nnz];
                    int& multiindex_index = index_facedof_local[n_dof_face*ind_face + facedof_index];
                    transform_n_local2global[ind_ele][multiindex_index]++;
                }
            }
        }
        for (int ind_interiordof = 0; ind_interiordof < n_dof_geometry[3]; ++ind_interiordof){ // index of interior dof -> corresponding multiindex 
            Multiindex<3> index_tmp = correspondence.number2index(ind_interiordof); //                                     -> 1-1 global interior coefficient
            index_tmp = index_tmp + Unitary_Multiindex[0]*2 + Unitary_Multiindex[1] + Unitary_Multiindex[2];
            int multiindex_index = correspondence.index2number(index_tmp);
            transform_n_local2global[ind_ele][multiindex_index]++;
        }

        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            int& n_dep = transform_n_local2global[ind_ele][ind_index];
            transform_ind_local2global[ind_ele][ind_index].resize(n_dep);
            transform_val_local2global[ind_ele][ind_index].resize(n_dep);
        }
        
        // traverse all dimensional geometries again, assign transform_local2global
        std::vector<int> pointer(n_index, 0);
        for (int ind_node = 0; ind_node < 4; ++ind_node){
            int& multiindex_index = ind_node;
            int& ind_geo = index_geometry_onelement[ind_ele][0][ind_node];
            int& ind_dof_global = location_actualdof[location_geometry[ind_geo]];
            transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_dof_global;
            transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = 1;
            pointer[multiindex_index]++;
        }
        for (int ind_edge = 0; ind_edge < 6; ++ind_edge){
            int& ind_geo = index_geometry_onelement[ind_ele][1][ind_edge];
            int& ind_startdof_global = location_actualdof[location_geometry[ind_geo + n_geometry[0]]];
            for (int ind_dof = 0; ind_dof < n_dof_edge; ++ind_dof){
                int& multiindex_index = index_edgedof_local[n_dof_edge*ind_edge + ind_dof];
                transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_startdof_global + ind_dof;
                transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = 1;
                if (type_projection[ind_ele][0][ind_edge] == 1 && ind_dof % 2 == 1)
                    transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = -1;
                pointer[multiindex_index]++;
            }
        }
        for (int ind_face = 0; ind_face < 4; ++ind_face){
            int& ind_var = type_projection[ind_ele][1][ind_face];
            int& ind_geo = index_geometry_onelement[ind_ele][2][ind_face];
            int& ind_startdof_global = location_actualdof[location_geometry[ind_geo + n_geometry[0] + n_geometry[1]]];
            for (int ind_dof = 0; ind_dof < n_dof_face; ++ind_dof){
                // int& multiindex_index = index_facedof_local[n_dof_face*ind_face + ind_dof];
                // add edge contribution
                for (int ind_e = 0; ind_e < 3; ++ind_e){
                    int& edge_index_local = index_edge_local[3*ind_face + ind_e]; // local index of this edge, between 0 and 5
                    for (int ind_nnz = 0; ind_nnz < n_expr_edge[ind_var][ind_dof][ind_e]; ++ind_nnz){
                        int& edgedof_index = ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        int& multiindex_index = index_edgedof_local[n_dof_edge*edge_index_local + edgedof_index];
                        transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_startdof_global + conversion[0][ind_dof];
                        transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = val_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        pointer[multiindex_index]++;
                    }
                }
                // add face contribution
                for (int ind_nnz = 0; ind_nnz < n_expr_face[ind_var][ind_dof]; ++ind_nnz){
                    int& facedof_index = ind_expr_face[ind_var][ind_dof][ind_nnz];
                    int& multiindex_index = index_facedof_local[n_dof_face*ind_face + facedof_index];
                    transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_startdof_global + conversion[0][ind_dof];
                    transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = val_expr_face[ind_var][ind_dof][ind_nnz];
                    pointer[multiindex_index]++;
                }
            }
        }
        int& ind_startdof_interior_global = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]];
        for (int ind_interiordof = 0; ind_interiordof < n_dof_geometry[3]; ++ind_interiordof){
            Multiindex<3> index_tmp = correspondence.number2index(ind_interiordof);
            index_tmp = index_tmp + Unitary_Multiindex[0]*2 + Unitary_Multiindex[1] + Unitary_Multiindex[2];
            int multiindex_index = correspondence.index2number(index_tmp);
            transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_startdof_interior_global + conversion[1][ind_interiordof];
            transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = 1;
            pointer[multiindex_index]++;
        }

        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            if (pointer[ind_index] != transform_n_local2global[ind_ele][ind_index])
                std::cerr << "error: assign transform_n_local2global\n";
    }
    std::cerr << "assign transform_global2local and transform_local2global\n";
    // for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
    // 	std::cout << "ind_ele = " << ind_ele << '\n';
    // 	for (int ind_index = 0; ind_index < n_index; ++ind_index){
    // 	    std::cout << "\tind_index = " << ind_index << ", transform_n_global2local = " << transform_n_global2local[ind_ele][ind_index]
    // 		      << ", (index, weight) =";
    // 	    for (int i = 0; i < transform_n_global2local[ind_ele][ind_index]; ++i)
    // 		std::cout << " (" << transform_ind_global2local[ind_ele][ind_index][i] << ',' << transform_val_global2local[ind_ele][ind_index][i] << ')';
    // 	    std::cout << '\n';
    // 	}
    // }

    // given correct transform_global2local, test transform_local2global
    Vector<valuetype> u_g(n_dof_total); // g for global
    for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
        u_g(ind_dof) = rand() * 1.0 / RAND_MAX;
    Vector<valuetype> u_l(n_index); // l for local
    for (int ind_index = 0; ind_index < n_index; ++ind_index)
        u_l(ind_index) = rand() * 1.0 / RAND_MAX;
    valuetype dif_max_global_local = -1;
    valuetype dif_max_local_global = -1;
    int n_geometry_local[4] = {4, 6, 4, 1};
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // initialize global u
        Vector<valuetype> u_g_t(n_dof_total);
        for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
            u_g_t(ind_dof) = 0;
        // assign global u by transform_local2global
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_local2global[ind_ele][ind_index]; ++ind_nnz)
                u_g_t(transform_ind_local2global[ind_ele][ind_index][ind_nnz]) += u_l(ind_index) * transform_val_local2global[ind_ele][ind_index][ind_nnz];
        // initialize local u
        Vector<valuetype> u_l_t(n_index);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            u_l_t(ind_index) = 0;
        // assign tmp local u by transform_global2local
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                u_l_t(ind_index) += u_g_t(transform_ind_global2local[ind_ele][ind_index][ind_nnz]) * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        u_l_t -= u_l;
        if (dif_max_global_local < u_l_t.linfty_norm()) dif_max_global_local = u_l_t.linfty_norm();

        // test global -> local -> global
        std::vector<bool> flag_involved(n_dof_total, false);
        for (int ind = 0; ind < 3; ++ind)
            for (int ind_geo_local = 0; ind_geo_local < n_geometry_local[ind]; ++ind_geo_local){
                int ind_dof_start = index_geometry_onelement[ind_ele][ind][ind_geo_local];
                for (int indt = 0; indt < ind; ++indt)
                    ind_dof_start += n_geometry[indt];
                int& ind_dof_start_global = location_actualdof[location_geometry[ind_dof_start]];
                for (int ind_dof = 0; ind_dof < n_dof_geometry[ind]; ++ind_dof)
                    flag_involved[ind_dof_start_global + ind_dof] = true;
        }
        for (int ind_dof = 0; ind_dof < n_dof_geometry[3]; ++ind_dof)
            flag_involved[location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]] + ind_dof] = true;
        Vector<valuetype> v_g(n_dof_total);
        for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
            if (flag_involved[ind_dof]) v_g(ind_dof) = u_g(ind_dof);
            else v_g(ind_dof) = 0;
        Vector<valuetype> v_l(n_index);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            v_l(ind_index) = 0;
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                v_l(ind_index) += v_g(transform_ind_global2local[ind_ele][ind_index][ind_nnz]) * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        Vector<valuetype> v_g_t(n_dof_total);
        for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
            v_g_t(ind_dof) = 0;
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_local2global[ind_ele][ind_index]; ++ind_nnz)
                v_g_t(transform_ind_local2global[ind_ele][ind_index][ind_nnz]) += v_l(ind_index) * transform_val_local2global[ind_ele][ind_index][ind_nnz];
        v_g_t -= v_g;
        if (dif_max_local_global < v_g_t.linfty_norm()) dif_max_local_global = v_g_t.linfty_norm();
    }
    std::cerr << "test transform_global2local and transform_local2global with rand vector,\n\tdif_max_global_local = " << dif_max_global_local << ", dif_max_local_global = " << dif_max_local_global << '\n';
}

TEMPLATE_TSEM
void THIS_TSEM::build_global_transform(Mesh<3> &mesh, FEMSpace<double, 3> &fem_space)
{
    // assign the weight_location, geometry_dimension and geometry_order of dof, which determine the order of these geometry in discretized matrix
    /* {weight_local, geometry_dimension, geometry_order}:
     *    weight_local:       a number for sorting, read from fem_space, determine the order of all dimensional geomtry
     *    geometry_dimension: the dimension of geometry corresponds to this set
     *    geometry_order:     the order of geometry in the same dimensional ones corresponds to this set
     */
    weight_location.resize(n_geometry_total); // the weight_location of 0: 3 dimensional geometry in turns
    geometry_dimension.resize(n_geometry_total);
    geometry_order.resize(n_geometry_total); // [i = 0:n_geometry_total-1] the order of i-th entry in weight_location, whose dimension is geometry_dimension[i]
    // correspondence between fem dof/element and geometry
    // transform_femdof2geometry: index from fem dof to 0 & 1 dimensional geometry, for 1 dimensional geometry, its index plus mesh.n_geometry(0)
    int n_dof = fem_space.n_dof();
    transform_femdof2geometry.resize(n_dof, -1);
    // location_geometry: location of all geometry (0:3 dimensional) according to increasing order of weight_location
    location_geometry.resize(n_geometry_total);
    // location_actualdof: start index of geometry in actual discretized matrix
    location_actualdof.resize(n_geometry_total);
    /* geometry in finite element space                       -> total index of all geometries     -> position in actual discretized system
     * transform_femdof2geometry (0 & 1 dimensional geometry) -> total index = 0: n_geometry_total -> location_actualdof[location_geometry[total index]]
     * transform_femele2geometry (2 & 3 dimensional geometry)
     */
    // assign geometry_dimension and geometry_order
    for (int ind = 0; ind <= 3; ++ind){
        int index_start = 0;
        for (int indt = 0; indt < ind; ++indt)
            index_start += n_geometry[indt];
        for (int i = 0; i < n_geometry[ind]; ++i){
            geometry_dimension[index_start + i] = ind;
            geometry_order[index_start + i] = i;
        }
    }
    // construct correspondence between fem dof and geometry (I): find weight_location for 0 & 1 dimensional geometry according to the order in fem_space
    for (int i = 0; i < n_dof; ++i)
        for (int ind = 0; ind <= 1 && transform_femdof2geometry[i] < 0; ++ind) // as all order in fem_space is natural number, so >= 0
            for (int j = 0; j < n_geometry[ind]; ++j)
                if (distance(point_ref_mesh[ind][j], fem_space.dofInfo(i).interp_point) < tol_zero){
                    transform_femdof2geometry[i] = ind * n_geometry[0] + j; // if ind == 1, add n_geometry[0], total index of point or edge
                    weight_location[transform_femdof2geometry[i]] = i;
                    break;
                }
    // calculate weight_location for 2 & 3 dimensional geometry
    for (int ind = 2; ind <= 3; ++ind){
        int index_start = n_geometry[0] + n_geometry[1] + (ind-2) * n_geometry[2]; // if ind == 3, add n_geometry[2]
        for (int i = 0; i < n_geometry[ind]; ++i){
            valuetype weight_tmp = 0.0;
            for (int indt = 0; indt <= ind; ++indt)
                weight_tmp += weight_location[mesh.geometry(ind, i).vertex(indt)];
            weight_tmp /= ind + 1;
            weight_location[index_start + i] = weight_tmp;
        }
    }
    // sort according to weight_location
    for (int i = 1; i < n_geometry_total; ++i)
        for (int j = i; j > 0; --j)
            if (weight_location[j] < weight_location[j-1]){
                valuetype tmp = weight_location[j];
                weight_location[j] = weight_location[j-1];
                weight_location[j-1] = tmp;
                int tmpi = geometry_dimension[j];
                geometry_dimension[j] = geometry_dimension[j-1];
                geometry_dimension[j-1] = tmpi;
                tmpi = geometry_order[j];
                geometry_order[j] = geometry_order[j-1];
                geometry_order[j-1] = tmpi;
            }
    for (int i = 0; i < n_geometry_total; ++i){ // assign location for geometry_order[i]-th geometry_dimension[i] dimensional geometry
        int index = geometry_order[i]; // recover the index of geometry in whole order
        for (int ind = 0; ind < geometry_dimension[i]; ++ind)
            index += n_geometry[ind];
        location_geometry[index] = i;
    }
    // insert dof into each location_geometry
    location_actualdof[0] = 0;
    for (int i = 1; i < n_geometry_total; ++i)
        location_actualdof[i] = location_actualdof[i-1] + n_dof_geometry[geometry_dimension[i-1]];
    std::cerr << "assign location_actualdof, location_geometry\n";
    
    
    // traverse element, construct transformation between local and global
    // initialize
    // assign conversion between lexicography order and the order in jia2022
    conversion.resize(2);
    conversion[0].resize(n_dof_face);
    conversion[1].resize(n_dof_geometry[3]);
    for (unsigned int l1 = 2; l1 <= M; ++l1)
	for (unsigned int l2 = 1; l2 <= M-l1; ++l2){
	    unsigned int ind_lex = (l1+l2+1-3) * (l1+l2-3) / 2 + l2-1; // location of (l1-2, l2-1) in lexicography order
	    // unsigned int ind_jia = (M+1)*(l1-2) - (l1-1-2)*(l1-2)/2 + l2-1;
	    // conversion[0][ind_lex] = ind_jia;
	    conversion[0][ind_lex] = ind_lex;
	}
    for (unsigned int ind_index = 0; ind_index < n_dof_geometry[3]; ++ind_index){
	Multiindex<3> index_now = correspondence.number2index(ind_index);
	unsigned int l1 = index_now.index[0], l2 = index_now.index[1], l3 = index_now.index[2];
	// unsigned int ind_jia = ((l1-1)*l1*(2*l1-1) + l1*(M+1)*(M+2)*6 - (2*M+3)*(l1-1)*l1*3) / 12
	//     + (M-l1+1)*l2 - (l2-1)*l2/2
	//     + l3;
	// conversion[1][ind_index] = ind_jia;
	conversion[1][ind_index] = ind_index;
    }
    std::cerr << "assign conversion\n";
    // expression of local coefficient by linear summation of global ones
    transform_n_global2local.resize(n_element);
    transform_ind_global2local.resize(n_element);
    transform_val_global2local.resize(n_element);
    // expression of global coefficient on each-dimensional geometry by local ones, equivalent to column compress of transformation matrix
    transform_n_local2global.resize(n_element);
    transform_ind_local2global.resize(n_element);
    transform_val_local2global.resize(n_element);
    // assign transform_local2global and transform_global2local
    std::vector<int> index_edgedof_local(n_dof_edge * 6); // index of multiindex corresponds to edge dof on local element
    for (int l = 2; l <= M; ++l){
        index_edgedof_local[n_dof_edge*0 + l-2] = correspondence.index2number(Unitary_Multiindex[0]*l); // (l, 0, 0)
        index_edgedof_local[n_dof_edge*1 + l-2] = correspondence.index2number(Unitary_Multiindex[1]*l); // (0, l, 0)
        index_edgedof_local[n_dof_edge*2 + l-2] = correspondence.index2number(Unitary_Multiindex[2]*l); // (0, 0, l)
        index_edgedof_local[n_dof_edge*4 + l-2] = correspondence.index2number(Unitary_Multiindex[0] + Unitary_Multiindex[2]*(l-1)); // (1,   0, l-1)
        index_edgedof_local[n_dof_edge*3 + l-2] = correspondence.index2number(Unitary_Multiindex[0] + Unitary_Multiindex[1]*(l-1)); // (1, l-1,   0)
        index_edgedof_local[n_dof_edge*5 + l-2] = correspondence.index2number(Unitary_Multiindex[1] + Unitary_Multiindex[2]*(l-1)); // (0,   1, l-1)
    }
    std::vector<int> index_facedof_local(n_dof_face * 4); // index of multiindex corresponds to face dof on local element
    for (int l1 = 2; l1 <= M; ++l1)
        for (int l2 = 1; l2 <= M-l1; ++l2){
            int ind_index = (l1+l2-3)*(l1+l2-2)/2 + l2-1;
            index_facedof_local[n_dof_face*0 + ind_index] =
                correspondence.index2number(Unitary_Multiindex[0] + Unitary_Multiindex[1]*(l1-1) + Unitary_Multiindex[2]*l2); // (1, l1-1, l2)
            index_facedof_local[n_dof_face*1 + ind_index] =
                correspondence.index2number(Unitary_Multiindex[1]*l1 + Unitary_Multiindex[2]*l2); // ( 0, l1, l2)
            index_facedof_local[n_dof_face*2 + ind_index] =
                correspondence.index2number(Unitary_Multiindex[0]*l1 + Unitary_Multiindex[2]*l2); // (l1,  0, l2)
            index_facedof_local[n_dof_face*3 + ind_index] =
                correspondence.index2number(Unitary_Multiindex[0]*l1 + Unitary_Multiindex[1]*l2); // (l1, l2,  0)
        }
    int index_edge_local[] = {5, 4, 3, // edge number for each face on local element
			      5, 2, 1,
			      4, 2, 0,
			      3, 1, 0};
    int order_edge_global[] = {0, 1, 2, // order of edge on global face in each variation
                               0, 2, 1,
                               1, 0, 2,
                               1, 2, 0,
                               2, 0, 1,
                               2, 1, 0};
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // initialize
        transform_n_global2local[ind_ele].resize(n_index);
        transform_ind_global2local[ind_ele].resize(n_index);
        transform_val_global2local[ind_ele].resize(n_index);
        transform_n_local2global[ind_ele].resize(n_index, 0);
        transform_ind_local2global[ind_ele].resize(n_index);
        transform_val_local2global[ind_ele].resize(n_index);
        
        // construct transformation node dof
        // P0, P1, P2, P3 excatly correspond first 4 multiindex (0,0,0), (1,0,0), (0,1,0), (0,0,1)
        for (int ind_node = 0; ind_node < 4; ++ind_node){ // traverse node: P0, P1, P2, P3
            int ind_index = ind_node;
            int ind_geo = index_geometry_onelement[ind_ele][0][ind_node];
            int ind_dof_global = location_actualdof[location_geometry[ind_geo]];
            int n_dep = 1; // number of dependency
            transform_n_global2local[ind_ele][ind_index] = n_dep;
            transform_ind_global2local[ind_ele][ind_index].resize(n_dep);
            transform_val_global2local[ind_ele][ind_index].resize(n_dep);
            transform_ind_global2local[ind_ele][ind_index][0] = ind_dof_global;
            transform_val_global2local[ind_ele][ind_index][0] = 1;
        }
        
        // construct transformation for edge dof
        for (int ind_edge = 0; ind_edge < 6; ++ind_edge){ // traverse edge: E01, E02, E03, E12, E13, E23
            int ind_geo = index_geometry_onelement[ind_ele][1][ind_edge];
            int ind_dof_start = location_actualdof[location_geometry[ind_geo + n_geometry[0]]];
            for (int l = 2; l <= M; ++l){
                int ind_index = index_edgedof_local[n_dof_edge*ind_edge + l-2];
                int ind_dof_global = ind_dof_start + l-2;
                int n_dep = 1;
                transform_n_global2local[ind_ele][ind_index] = n_dep;
                transform_ind_global2local[ind_ele][ind_index].resize(n_dep);
                transform_val_global2local[ind_ele][ind_index].resize(n_dep);
                transform_ind_global2local[ind_ele][ind_index][0] = ind_dof_global;
                transform_val_global2local[ind_ele][ind_index][0] = 1;
                if (type_projection[ind_ele][0][ind_edge] == 1 && l%2 == 1)
                    transform_val_global2local[ind_ele][ind_index][0] = -1;
            }
        }
        
        // construct transformation for face dof
        for (int ind_face = 0; ind_face < 4; ++ind_face){
            int ind_geo = index_geometry_onelement[ind_ele][2][ind_face];
            int ind_facedof_start = location_actualdof[location_geometry[ind_geo + n_geometry[0] + n_geometry[1]]];
            for (int ind_dof = 0; ind_dof < n_dof_geometry[2]; ++ind_dof){
                int ind_index = index_facedof_local[n_dof_face*ind_face + ind_dof];
                int ind_var = type_projection_inv[ind_ele][ind_face];
                int n_dep = 0;
                for (int ind_e = 0; ind_e < 3; ++ind_e)
                    n_dep += n_expr_edge[ind_var][ind_dof][ind_e];
                n_dep += n_expr_face[ind_var][ind_dof];
                transform_n_global2local[ind_ele][ind_index] = n_dep;
                transform_ind_global2local[ind_ele][ind_index].resize(n_dep);
                transform_val_global2local[ind_ele][ind_index].resize(n_dep);
                // assign edge contribution
                int pos_start = 0;
                for (int ind_e = 0; ind_e < 3; ++ind_e){
                    int ind_edge = number_edge[index_geometry_onelement[ind_ele][2][ind_face]][ind_e];
                    int ind_edgedof_start = location_actualdof[location_geometry[ind_edge + n_geometry[0]]];
                    for (int ind_nnz = 0; ind_nnz < n_expr_edge[ind_var][ind_dof][ind_e]; ++ind_nnz){
                        transform_ind_global2local[ind_ele][ind_index][pos_start + ind_nnz]
                            = ind_edgedof_start + ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        transform_val_global2local[ind_ele][ind_index][pos_start + ind_nnz]
                            = val_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        if (!flag_sameorder_edgeonface[index_geometry_onelement[ind_ele][2][ind_face]][ind_e]
                            && ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz] % 2 == 1)
                            transform_val_global2local[ind_ele][ind_index][pos_start + ind_nnz] *= -1;
                    }
                    pos_start += n_expr_edge[ind_var][ind_dof][ind_e];
                }
                // assign face contribution
                for (int ind_nnz = 0; ind_nnz < n_expr_face[ind_var][ind_dof]; ++ind_nnz){
                    transform_ind_global2local[ind_ele][ind_index][pos_start + ind_nnz]
                        = ind_facedof_start + conversion[0][ind_expr_face[ind_var][ind_dof][ind_nnz]];
                    transform_val_global2local[ind_ele][ind_index][pos_start + ind_nnz]
                        = val_expr_face[ind_var][ind_dof][ind_nnz];
                }
                if (pos_start + n_expr_face[ind_var][ind_dof] != n_dep)
                    std::cerr << "error: assign transform_global2local for face dof\n";
            }
        }
        
        // construct transformation for interior dof
        int ind_interiordof_start = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]];
        for (int ind_dof = 0; ind_dof < n_dof_geometry[3]; ++ind_dof){
            Multiindex<3> index_now = correspondence.number2index(ind_dof);
            index_now = index_now + Unitary_Multiindex[0]*2 + Unitary_Multiindex[1] + Unitary_Multiindex[2];
            int ind_index = correspondence.index2number(index_now);
            int n_dep = 1;
            transform_n_global2local[ind_ele][ind_index] = n_dep;
            transform_ind_global2local[ind_ele][ind_index].resize(n_dep);
            transform_val_global2local[ind_ele][ind_index].resize(n_dep);
            transform_ind_global2local[ind_ele][ind_index][0] = ind_interiordof_start + conversion[1][ind_dof];
            transform_val_global2local[ind_ele][ind_index][0] = 1;
        }

        // traverse all dimensional geometries, count contribution for transform_local2global
        for (int ind_node = 0; ind_node < 4; ++ind_node) // ind_node -> ind_index -> 1-1 global vertex coefficient
            transform_n_local2global[ind_ele][ind_node]++;
        for (int ind_edge = 0; ind_edge < 6; ++ind_edge) // ind_edge, ind_dof -> ind_index -> 1-1 global edge coefficient
            for (int ind_dof = 0; ind_dof < n_dof_edge; ++ind_dof)
                transform_n_local2global[ind_ele][index_edgedof_local[n_dof_edge*ind_edge + ind_dof]]++;
        for (int ind_face = 0; ind_face < 4; ++ind_face){
            int& ind_var = type_projection[ind_ele][1][ind_face];
            for (int ind_dof = 0; ind_dof < n_dof_face; ++ind_dof){
                // count edge contribution
                for (int ind_e = 0; ind_e < 3; ++ind_e){ // ind_dof, ind_e -> edge index -> corresponding multiindex
                    int& edge_index_local = index_edge_local[3*ind_face + ind_e]; // local index of this edge, between 0 and 5
                    for (int ind_nnz = 0; ind_nnz < n_expr_edge[ind_var][ind_dof][ind_e]; ++ind_nnz){
                        int& edgedof_index = ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        int& multiindex_index = index_edgedof_local[n_dof_edge*edge_index_local + edgedof_index];
                        transform_n_local2global[ind_ele][multiindex_index]++;
                    }
                }
                // count face contribution
                for (int ind_nnz = 0; ind_nnz < n_expr_face[ind_var][ind_dof]; ++ind_nnz){
                    int& facedof_index = ind_expr_face[ind_var][ind_dof][ind_nnz];
                    int& multiindex_index = index_facedof_local[n_dof_face*ind_face + facedof_index];
                    transform_n_local2global[ind_ele][multiindex_index]++;
                }
            }
        }
        for (int ind_interiordof = 0; ind_interiordof < n_dof_geometry[3]; ++ind_interiordof){ // index of interior dof -> corresponding multiindex 
            Multiindex<3> index_tmp = correspondence.number2index(ind_interiordof); //                                     -> 1-1 global interior coefficient
            index_tmp = index_tmp + Unitary_Multiindex[0]*2 + Unitary_Multiindex[1] + Unitary_Multiindex[2];
            int multiindex_index = correspondence.index2number(index_tmp);
            transform_n_local2global[ind_ele][multiindex_index]++;
        }

        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            int& n_dep = transform_n_local2global[ind_ele][ind_index];
            transform_ind_local2global[ind_ele][ind_index].resize(n_dep);
            transform_val_local2global[ind_ele][ind_index].resize(n_dep);
        }
        
        // traverse all dimensional geometries again, assign transform_local2global
        std::vector<int> pointer(n_index, 0);
        for (int ind_node = 0; ind_node < 4; ++ind_node){
            int& multiindex_index = ind_node;
            int& ind_geo = index_geometry_onelement[ind_ele][0][ind_node];
            int& ind_dof_global = location_actualdof[location_geometry[ind_geo]];
            transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_dof_global;
            transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = 1;
            pointer[multiindex_index]++;
        }
        for (int ind_edge = 0; ind_edge < 6; ++ind_edge){
            int& ind_geo = index_geometry_onelement[ind_ele][1][ind_edge];
            int& ind_startdof_global = location_actualdof[location_geometry[ind_geo + n_geometry[0]]];
            for (int ind_dof = 0; ind_dof < n_dof_edge; ++ind_dof){
                int& multiindex_index = index_edgedof_local[n_dof_edge*ind_edge + ind_dof];
                transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_startdof_global + ind_dof;
                transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = 1;
                if (type_projection[ind_ele][0][ind_edge] == 1 && ind_dof % 2 == 1)
                    transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = -1;
                pointer[multiindex_index]++;
            }
        }
        for (int ind_face = 0; ind_face < 4; ++ind_face){
            int& ind_var = type_projection[ind_ele][1][ind_face];
            int& ind_geo = index_geometry_onelement[ind_ele][2][ind_face];
            int& ind_startdof_global = location_actualdof[location_geometry[ind_geo + n_geometry[0] + n_geometry[1]]];
            for (int ind_dof = 0; ind_dof < n_dof_face; ++ind_dof){
                // int& multiindex_index = index_facedof_local[n_dof_face*ind_face + ind_dof];
                // add edge contribution
                for (int ind_e = 0; ind_e < 3; ++ind_e){
                    int& edge_index_local = index_edge_local[3*ind_face + ind_e]; // local index of this edge, between 0 and 5
                    for (int ind_nnz = 0; ind_nnz < n_expr_edge[ind_var][ind_dof][ind_e]; ++ind_nnz){
                        int& edgedof_index = ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        int& multiindex_index = index_edgedof_local[n_dof_edge*edge_index_local + edgedof_index];
                        transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_startdof_global + conversion[0][ind_dof];
                        transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = val_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                        pointer[multiindex_index]++;
                    }
                }
                // add face contribution
                for (int ind_nnz = 0; ind_nnz < n_expr_face[ind_var][ind_dof]; ++ind_nnz){
                    int& facedof_index = ind_expr_face[ind_var][ind_dof][ind_nnz];
                    int& multiindex_index = index_facedof_local[n_dof_face*ind_face + facedof_index];
                    transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_startdof_global + conversion[0][ind_dof];
                    transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = val_expr_face[ind_var][ind_dof][ind_nnz];
                    pointer[multiindex_index]++;
                }
            }
        }
        int& ind_startdof_interior_global = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]];
        for (int ind_interiordof = 0; ind_interiordof < n_dof_geometry[3]; ++ind_interiordof){
            Multiindex<3> index_tmp = correspondence.number2index(ind_interiordof);
            index_tmp = index_tmp + Unitary_Multiindex[0]*2 + Unitary_Multiindex[1] + Unitary_Multiindex[2];
            int multiindex_index = correspondence.index2number(index_tmp);
            transform_ind_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = ind_startdof_interior_global + conversion[1][ind_interiordof];
            transform_val_local2global[ind_ele][multiindex_index][pointer[multiindex_index]] = 1;
            pointer[multiindex_index]++;
        }

        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            if (pointer[ind_index] != transform_n_local2global[ind_ele][ind_index])
                std::cerr << "error: assign transform_n_local2global\n";
    }
    std::cerr << "assign transform_global2local and transform_local2global\n";

    // given correct transform_global2local, test transform_local2global
    Vector<valuetype> u_g(n_dof_total); // g for global
    for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
        u_g(ind_dof) = rand() * 1.0 / RAND_MAX;
    Vector<valuetype> u_l(n_index); // l for local
    for (int ind_index = 0; ind_index < n_index; ++ind_index)
        u_l(ind_index) = rand() * 1.0 / RAND_MAX;
    valuetype dif_max_global_local = -1;
    valuetype dif_max_local_global = -1;
    int n_geometry_local[4] = {4, 6, 4, 1};
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // initialize global u
        Vector<valuetype> u_g_t(n_dof_total);
        for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
            u_g_t(ind_dof) = 0;
        // assign global u by transform_local2global
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_local2global[ind_ele][ind_index]; ++ind_nnz)
                u_g_t(transform_ind_local2global[ind_ele][ind_index][ind_nnz]) += u_l(ind_index) * transform_val_local2global[ind_ele][ind_index][ind_nnz];
        // initialize local u
        Vector<valuetype> u_l_t(n_index);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            u_l_t(ind_index) = 0;
        // assign tmp local u by transform_global2local
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                u_l_t(ind_index) += u_g_t(transform_ind_global2local[ind_ele][ind_index][ind_nnz]) * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        u_l_t -= u_l;
        if (dif_max_global_local < u_l_t.linfty_norm()) dif_max_global_local = u_l_t.linfty_norm();

        // test global -> local -> global
        std::vector<bool> flag_involved(n_dof_total, false);
        for (int ind = 0; ind < 3; ++ind)
            for (int ind_geo_local = 0; ind_geo_local < n_geometry_local[ind]; ++ind_geo_local){
                int ind_dof_start = index_geometry_onelement[ind_ele][ind][ind_geo_local];
                for (int indt = 0; indt < ind; ++indt)
                    ind_dof_start += n_geometry[indt];
                int& ind_dof_start_global = location_actualdof[location_geometry[ind_dof_start]];
                for (int ind_dof = 0; ind_dof < n_dof_geometry[ind]; ++ind_dof)
                    flag_involved[ind_dof_start_global + ind_dof] = true;
        }
        for (int ind_dof = 0; ind_dof < n_dof_geometry[3]; ++ind_dof)
            flag_involved[location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]] + ind_dof] = true;
        Vector<valuetype> v_g(n_dof_total);
        for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
            if (flag_involved[ind_dof]) v_g(ind_dof) = u_g(ind_dof);
            else v_g(ind_dof) = 0;
        Vector<valuetype> v_l(n_index);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            v_l(ind_index) = 0;
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                v_l(ind_index) += v_g(transform_ind_global2local[ind_ele][ind_index][ind_nnz]) * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        Vector<valuetype> v_g_t(n_dof_total);
        for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
            v_g_t(ind_dof) = 0;
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_local2global[ind_ele][ind_index]; ++ind_nnz)
                v_g_t(transform_ind_local2global[ind_ele][ind_index][ind_nnz]) += v_l(ind_index) * transform_val_local2global[ind_ele][ind_index][ind_nnz];
        v_g_t -= v_g;
        if (dif_max_local_global < v_g_t.linfty_norm()) dif_max_local_global = v_g_t.linfty_norm();
    }
    std::cerr << "test transform_global2local and transform_local2global with rand vector,\n\tdif_max_global_local = " << dif_max_global_local << ", dif_max_local_global = " << dif_max_local_global << '\n';
}

TEMPLATE_TSEM
void THIS_TSEM::read_interp_info(const std::string &interp_filename)
{
    // number of nonzero contribution on edges
    n_expr_edge.resize(6);;
    // index of nonzero contribution on edges
    ind_expr_edge.resize(6);
    // value of nonzero contribution on edges
    val_expr_edge.resize(6);
    // number of nonzero contribution on face
    n_expr_face.resize(6);
    // index of nonzero contribution on face
    ind_expr_face.resize(6);
    // value of nonzero contribution on face
    val_expr_face.resize(6);
    
    // read projection info from file Info_Interpolate
    std::ifstream input_interp(interp_filename);
    int M_tmp;
    input_interp >> M_tmp;
    n_dof_edge = M - 1;
    n_dof_face = (M-2) * (M-3) / 2 + M-2;
    int n_dof_edge_read = M_tmp - 1;
    int n_dof_face_read = (M_tmp-2) * (M_tmp-3) / 2 + M_tmp-2;
    
    // read info
    for (int ind_var = 0; ind_var < 6; ++ind_var){
        // initialize
        n_expr_edge[ind_var].resize(n_dof_face_read);
        ind_expr_edge[ind_var].resize(n_dof_face_read);
        val_expr_edge[ind_var].resize(n_dof_face_read);
        n_expr_face[ind_var].resize(n_dof_face_read);
        ind_expr_face[ind_var].resize(n_dof_face_read);
        val_expr_face[ind_var].resize(n_dof_face_read);
        // traverse dof, read edge interpolation info
        for (int ind_dof = 0; ind_dof < n_dof_face_read; ++ind_dof){
            n_expr_edge[ind_var][ind_dof].resize(3);
            ind_expr_edge[ind_var][ind_dof].resize(3);
            val_expr_edge[ind_var][ind_dof].resize(3);
            for (int ind_e = 0; ind_e < 3; ++ind_e){
                int cnt;
                input_interp >> cnt;
                n_expr_edge[ind_var][ind_dof][ind_e] = cnt;
                ind_expr_edge[ind_var][ind_dof][ind_e].resize(cnt);
                val_expr_edge[ind_var][ind_dof][ind_e].resize(cnt);
                for (int ind_nnz = 0; ind_nnz < cnt; ++ind_nnz){
                    input_interp >> ind_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                    input_interp >> val_expr_edge[ind_var][ind_dof][ind_e][ind_nnz];
                }
            }
        }
        // traverse dof, read face interpolation info
        for (int ind_dof = 0; ind_dof < n_dof_face_read; ++ind_dof){
            int count;
            input_interp >> count;
            n_expr_face[ind_var][ind_dof] = count;
            ind_expr_face[ind_var][ind_dof].resize(count);
            val_expr_face[ind_var][ind_dof].resize(count);
            for (int ind_nnz = 0; ind_nnz < count; ++ind_nnz){
                input_interp >> ind_expr_face[ind_var][ind_dof][ind_nnz];
                input_interp >> val_expr_face[ind_var][ind_dof][ind_nnz];
            }
        }
    }
    input_interp.close();
    std::cerr << "read projection info, M_tmp = " << M_tmp
              << ", n_dof_edge_read = " << n_dof_edge_read << ", n_dof_face_read = " << n_dof_face_read
              << '\n';
}

#undef TEMPLATE_TSEM
#undef THIS_TSEM

