#include "../include/TetrahedralSEM.h"
// This file contains functions calculating error or difference in kinds of NORM.
// Functions:
// calc_l2_error();
// calc_l2_error_gradient();
// calc_h1_error();
// calc_l2_error_gradient_difference();
// calc_l2_density_difference();
// calc_l2_difference();
// calc_h1_difference();
// calc_l2_error_component();
// calc_density_l2_difference();


#define TEMPLATE_TSEM template<typename valuetype>
#define THIS_TSEM TSEM<valuetype>


TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_error(Vector<valuetype> &sol, std::vector<std::vector<valuetype> > &val_exc)
{
    // for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
    // 	std::cout << sol(ind_dof) << ' ';
    // std::cout << '\n';
    valuetype err_l2 = 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	// recover local coefficient on this element
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += sol(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
	// for (int ind_index = 0; ind_index < n_index; ++ind_index)
	//     std::cout << val_coef[ind_index] << ' ';
	// std::cout << '\n';

	valuetype count = 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
            valuetype val_sol = 0;
            for (int ind_index = 0; ind_index < n_index; ++ind_index)
                val_sol += val_coef[ind_index] * basis_value_actual[p][ind_index];
            // count the contribution
            count += Weight[2][p] * pow(val_sol - val_exc[ind_ele][p], 2);
	    // valuetype dif = fabs(val_sol - val_exc[ind_ele][p]);
	    // if (dif > 1.0e-8)
	    // 	std::cout << "val_sol = " << val_sol << ", val_exc = " << val_exc[ind_ele][p] << ", dif = " << dif << ", count = " << count << ", Weight[2][p] = " << Weight[2][p] << '\n';
        }
        err_l2 += count * val_volume[ind_ele];
	// std::cout << "ind_ele = " << ind_ele << ", count = " << count << '\n';
    }
    err_l2 = sqrt(err_l2);
    // std::cout << "err_l2 = " << err_l2 << '\n';
    
    return err_l2;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_error(Vector<valuetype> &sol, FEMSpace<double, 3> &fem_space, valuetype (*func)(double*))
{
    valuetype err_l2 = 0;
    std::vector<valuetype> value_exact(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(ind_ele).local_to_global(QPoint);
	for (int p = 0; p < n_q_point[2]; ++p)
	    value_exact[p] = func(q_point[p]);
	// recover local coefficient on this element
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += sol(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
	valuetype count = 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
            valuetype val_sol = 0;
            for (int ind_index = 0; ind_index < n_index; ++ind_index)
                val_sol += val_coef[ind_index] * basis_value_actual[p][ind_index];
            // count the contribution
            count += Weight[2][p] * pow(val_sol - value_exact[p], 2);
        }
        err_l2 += count * val_volume[ind_ele];
    }
    err_l2 = sqrt(err_l2);
    
    return err_l2;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_error(Vector<valuetype> &sol, FEMSpace<double, 3> &fem_space, valuetype (*func)(double*, double), valuetype t)
{
    valuetype err_l2 = 0;
    std::vector<valuetype> value_exact(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(ind_ele).local_to_global(QPoint);
	for (int p = 0; p < n_q_point[2]; ++p)
	    value_exact[p] = func(q_point[p], t);
	// recover local coefficient on this element
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += sol(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
	valuetype count = 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
            valuetype val_sol = 0;
            for (int ind_index = 0; ind_index < n_index; ++ind_index)
                val_sol += val_coef[ind_index] * basis_value_actual[p][ind_index];
            // count the contribution
            count += Weight[2][p] * pow(val_sol - value_exact[p], 2);
        }
        err_l2 += count * val_volume[ind_ele];
    }
    err_l2 = sqrt(err_l2);
    
    return err_l2;
}

TEMPLATE_TSEM
valuetype THIS_TSEM:: calc_l2_error_gradient(RegularMesh<3> &mesh,
					     Vector<valuetype> &sol, std::vector<std::vector<std::vector<valuetype> > > &val_g_exc)
{
    valuetype err = (valuetype) 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	// valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	// recover local coefficient on this element
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += sol(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	// prepare coordinates of element vertices
	valuetype xj_k[3][4];
	for (int k = 0; k < 4; ++k){
	    int ind_node = index_geometry_onelement[ind_ele][0][k];
	    for (int j = 0; j < 3; ++j){
		xj_k[j][k] = mesh.point(ind_node)[j];
		// std::cout << "xj_k[" << j << "][" << k << "] = " << xj_k[j][k] << '\n';
	    }
	}
	valuetype d_jk[3][3];
	for (int j = 0; j < 3; ++j)
	    for (int k = 0; k < 3; ++k){
		d_jk[j][k] = xj_k[j][k+1] - xj_k[j][0];
		// std::cout << "d_jk[" << j << "][" << k << "] = xj_k[" << j << "][" << k+1 << "] - xj_k[" << j << "][0] = " << xj_k[j][k+1] << "\t- " << xj_k[j][0] << " = " << d_jk[j][k] << '\n';
	    }
	valuetype val_V = d_jk[0][0]*d_jk[1][1]*d_jk[2][2] + d_jk[0][1]*d_jk[1][2]*d_jk[2][0] + d_jk[0][2]*d_jk[1][0]*d_jk[2][1]
	                - d_jk[0][0]*d_jk[1][2]*d_jk[2][1] - d_jk[0][1]*d_jk[1][0]*d_jk[2][2] - d_jk[0][2]*d_jk[1][1]*d_jk[2][0];
	// std::cout << "ind_ele = " << ind_ele << ", val_V = " << val_V << '\n';
	// break;
	if (fabs(fabs(val_V) - 6*val_volume[ind_ele])/fabs(val_V) > tol_zero) std::cout << "find different volume, ind_ele = " << ind_ele << ", volume = " << val_volume[ind_ele] << ", val_V = " << fabs(val_V/6) << ", dif = " << fabs(fabs(val_volume[ind_ele] - fabs(val_V/6))) << '\n';
	
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
	    
	valuetype count = (valuetype) 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
	    valuetype val_g_sol[3] = {0, 0, 0};
	    for (int ind = 0; ind < 3; ++ind)
		for (int ind_index = 0; ind_index < n_index; ++ind_index)
		    for (int indt = 0; indt < 3; ++indt)
			val_g_sol[ind] += parxj_parxk[indt][ind] * val_coef[ind_index] * basis_gradient_actual[p][ind_index][indt];
            // count the contribution
	    valuetype cnt_tmp = (valuetype) 0;
	    for (int ind = 0; ind < 3; ++ind)
		cnt_tmp += pow(val_g_sol[ind] - val_g_exc[ind_ele][p][ind], 2);
            count += Weight[2][p] * cnt_tmp;
	}
	err += count * val_volume[ind_ele];
    }
    
    err = sqrt(err);
    return err;
}


TEMPLATE_TSEM
valuetype THIS_TSEM:: calc_l2_error_gradient(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space,
					     Vector<valuetype> &sol, void (*func)(double*, valuetype*))
{
    valuetype err = (valuetype) 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	// recover local coefficient on this element
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += sol(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	// prepare coordinates of element vertices
	valuetype xj_k[3][4];
	for (int k = 0; k < 4; ++k){
	    int ind_node = index_geometry_onelement[ind_ele][0][k];
	    for (int j = 0; j < 3; ++j){
		xj_k[j][k] = mesh.point(ind_node)[j];
		// std::cout << "xj_k[" << j << "][" << k << "] = " << xj_k[j][k] << '\n';
	    }
	}
	valuetype d_jk[3][3];
	for (int j = 0; j < 3; ++j)
	    for (int k = 0; k < 3; ++k){
		d_jk[j][k] = xj_k[j][k+1] - xj_k[j][0];
		// std::cout << "d_jk[" << j << "][" << k << "] = xj_k[" << j << "][" << k+1 << "] - xj_k[" << j << "][0] = " << xj_k[j][k+1] << "\t- " << xj_k[j][0] << " = " << d_jk[j][k] << '\n';
	    }
	valuetype val_V = d_jk[0][0]*d_jk[1][1]*d_jk[2][2] + d_jk[0][1]*d_jk[1][2]*d_jk[2][0] + d_jk[0][2]*d_jk[1][0]*d_jk[2][1]
	                - d_jk[0][0]*d_jk[1][2]*d_jk[2][1] - d_jk[0][1]*d_jk[1][0]*d_jk[2][2] - d_jk[0][2]*d_jk[1][1]*d_jk[2][0];
	// std::cout << "ind_ele = " << ind_ele << ", val_V = " << val_V << '\n';
	// break;
	if (fabs(fabs(val_V) - 6*val_volume[ind_ele])/fabs(val_V) > tol_zero) std::cout << "find different volume, ind_ele = " << ind_ele << ", volume = " << val_volume[ind_ele] << ", val_V = " << fabs(val_V/6) << ", dif = " << fabs(fabs(val_volume[ind_ele] - fabs(val_V/6))) << '\n';
	
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
	    
	valuetype count = (valuetype) 0;
        for (int p = 0; p < n_q_point[2]; ++p){
	    // calculate gradient value
	    AFEPack::Point<3> qp = fem_space.element(ind_ele).local_to_global(QPoint[p]);
	    valuetype g_ref[3]; func(qp, g_ref);
            // recover the value of solution
	    valuetype val_g_sol[3] = {0, 0, 0};
	    for (int ind = 0; ind < 3; ++ind)
		for (int ind_index = 0; ind_index < n_index; ++ind_index)
		    for (int indt = 0; indt < 3; ++indt)
			val_g_sol[ind] += parxj_parxk[indt][ind] * val_coef[ind_index] * basis_gradient_actual[p][ind_index][indt];
            // count the contribution
	    valuetype cnt_tmp = (valuetype) 0;
	    for (int ind = 0; ind < 3; ++ind)
		cnt_tmp += pow(val_g_sol[ind] - g_ref[ind], 2);
            count += Weight[2][p] * cnt_tmp;
	}
	err += count * val_volume[ind_ele];
    }
    
    err = sqrt(err);
    return err;
}

TEMPLATE_TSEM
valuetype THIS_TSEM:: calc_l2_error_gradient(Mesh<3> &mesh,
					     Vector<valuetype> &sol, std::vector<std::vector<std::vector<valuetype> > > &val_g_exc)
{
    valuetype err = (valuetype) 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	// valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	// recover local coefficient on this element
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += sol(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	// prepare coordinates of element vertices
	valuetype xj_k[3][4];
	for (int k = 0; k < 4; ++k){
	    int ind_node = index_geometry_onelement[ind_ele][0][k];
	    for (int j = 0; j < 3; ++j){
		xj_k[j][k] = mesh.point(ind_node)[j];
		// std::cout << "xj_k[" << j << "][" << k << "] = " << xj_k[j][k] << '\n';
	    }
	}
	valuetype d_jk[3][3];
	for (int j = 0; j < 3; ++j)
	    for (int k = 0; k < 3; ++k){
		d_jk[j][k] = xj_k[j][k+1] - xj_k[j][0];
		// std::cout << "d_jk[" << j << "][" << k << "] = xj_k[" << j << "][" << k+1 << "] - xj_k[" << j << "][0] = " << xj_k[j][k+1] << "\t- " << xj_k[j][0] << " = " << d_jk[j][k] << '\n';
	    }
	valuetype val_V = d_jk[0][0]*d_jk[1][1]*d_jk[2][2] + d_jk[0][1]*d_jk[1][2]*d_jk[2][0] + d_jk[0][2]*d_jk[1][0]*d_jk[2][1]
	                - d_jk[0][0]*d_jk[1][2]*d_jk[2][1] - d_jk[0][1]*d_jk[1][0]*d_jk[2][2] - d_jk[0][2]*d_jk[1][1]*d_jk[2][0];
	// std::cout << "ind_ele = " << ind_ele << ", val_V = " << val_V << '\n';
	// break;
	if (fabs(fabs(val_V) - 6*val_volume[ind_ele])/fabs(val_V) > tol_zero) std::cout << "find different volume, ind_ele = " << ind_ele << ", volume = " << val_volume[ind_ele] << ", val_V = " << fabs(val_V/6) << ", dif = " << fabs(fabs(val_volume[ind_ele] - fabs(val_V/6))) << '\n';
	
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
	    
	valuetype count = (valuetype) 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
	    valuetype val_g_sol[3] = {0, 0, 0};
	    for (int ind = 0; ind < 3; ++ind)
		for (int ind_index = 0; ind_index < n_index; ++ind_index)
		    for (int indt = 0; indt < 3; ++indt)
			val_g_sol[ind] += parxj_parxk[indt][ind] * val_coef[ind_index] * basis_gradient_actual[p][ind_index][indt];
            // count the contribution
	    valuetype cnt_tmp = (valuetype) 0;
	    for (int ind = 0; ind < 3; ++ind)
		cnt_tmp += pow(val_g_sol[ind] - val_g_exc[ind_ele][p][ind], 2);
            count += Weight[2][p] * cnt_tmp;
	}
	err += count * val_volume[ind_ele];
    }
    
    err = sqrt(err);
    return err;
}

TEMPLATE_TSEM
valuetype THIS_TSEM:: calc_l2_error_gradient(Mesh<3> &mesh,
					     Vector<valuetype> &sol, std::vector<std::vector<valuetype> > &val_g_exc, unsigned int idx_ele)
{
    valuetype err = (valuetype) 0;

    // recover local coefficient on this element
    std::vector<valuetype> val_coef(n_index);
    calc_coef_onElement(sol, val_coef, idx_ele);
    
    // prepare coordinates of element vertices
    valuetype xj_k[3][4];
    for (int k = 0; k < 4; ++k){
	int ind_node = index_geometry_onelement[idx_ele][0][k];
	for (int j = 0; j < 3; ++j){
	    xj_k[j][k] = mesh.point(ind_node)[j];
	    // std::cout << "xj_k[" << j << "][" << k << "] = " << xj_k[j][k] << '\n';
	}
    }
    valuetype d_jk[3][3];
    for (int j = 0; j < 3; ++j)
	for (int k = 0; k < 3; ++k){
	    d_jk[j][k] = xj_k[j][k+1] - xj_k[j][0];
	    // std::cout << "d_jk[" << j << "][" << k << "] = xj_k[" << j << "][" << k+1 << "] - xj_k[" << j << "][0] = " << xj_k[j][k+1] << "\t- " << xj_k[j][0] << " = " << d_jk[j][k] << '\n';
	}
    valuetype val_V = d_jk[0][0]*d_jk[1][1]*d_jk[2][2] + d_jk[0][1]*d_jk[1][2]*d_jk[2][0] + d_jk[0][2]*d_jk[1][0]*d_jk[2][1]
	- d_jk[0][0]*d_jk[1][2]*d_jk[2][1] - d_jk[0][1]*d_jk[1][0]*d_jk[2][2] - d_jk[0][2]*d_jk[1][1]*d_jk[2][0];
    // std::cout << "ind_ele = " << ind_ele << ", val_V = " << val_V << '\n';
    // break;
    if (fabs(fabs(val_V) - 6*val_volume[idx_ele])/fabs(val_V) > tol_zero) std::cout << "find different volume, idx_ele = " << idx_ele << ", volume = " << val_volume[idx_ele] << ", val_V = " << fabs(val_V/6) << ", dif = " << fabs(fabs(val_volume[idx_ele] - fabs(val_V/6))) << '\n';
	
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
	    
    valuetype count = (valuetype) 0;
    for (int p = 0; p < n_q_point[2]; ++p){
	// recover the value of solution
	valuetype val_g_sol[3] = {0, 0, 0};
	for (int ind = 0; ind < 3; ++ind)
	    for (int ind_index = 0; ind_index < n_index; ++ind_index)
		for (int indt = 0; indt < 3; ++indt)
		    val_g_sol[ind] += parxj_parxk[indt][ind] * val_coef[ind_index] * basis_gradient_actual[p][ind_index][indt];
	// count the contribution
	valuetype cnt_tmp = (valuetype) 0;
	for (int ind = 0; ind < 3; ++ind)
	    cnt_tmp += pow(val_g_sol[ind] - val_g_exc[p][ind], 2);
	count += Weight[2][p] * cnt_tmp;
    }
    err += count * val_volume[idx_ele];
    
    err = sqrt(err);
    return err;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_h1_error(RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space,
				   Vector<valuetype> &sol, valuetype (*func)(double*), void (*func_g)(double*, valuetype*))
{
    double err = sqrt(pow(calc_l2_error(sol, fem_space, func), 2) +
		      pow(calc_l2_error_gradient(mesh, fem_space, sol, func_g), 2));
    return err;
}

TEMPLATE_TSEM
valuetype THIS_TSEM:: calc_l2_error_gradient_difference(RegularMesh<3> &mesh,
							Vector<valuetype> &u, Vector<valuetype> &v)
{
    valuetype err = (valuetype) 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	// recover local coefficient on this element
        std::vector<valuetype> val_coef_dif(n_index, 0), val_coef_v(n_index, 0);
	calc_coef_onElement(u, val_coef_dif, ind_ele);
	calc_coef_onElement(v, val_coef_v, ind_ele);
	for (unsigned int idx_index = 0; idx_index < n_index; ++idx_index) val_coef_dif[idx_index] -= val_coef_v[idx_index];

	// prepare coordinates of element vertices
	valuetype xj_k[3][4];
	for (int k = 0; k < 4; ++k){
	    int &ind_node = index_geometry_onelement[ind_ele][0][k];
	    for (int j = 0; j < 3; ++j) xj_k[j][k] = mesh.point(ind_node)[j];
	}
	valuetype d_jk[3][3];
	for (int j = 0; j < 3; ++j)
	    for (int k = 0; k < 3; ++k) d_jk[j][k] = xj_k[j][k+1] - xj_k[j][0];
	valuetype val_V = d_jk[0][0]*d_jk[1][1]*d_jk[2][2] + d_jk[0][1]*d_jk[1][2]*d_jk[2][0] + d_jk[0][2]*d_jk[1][0]*d_jk[2][1]
	    - d_jk[0][0]*d_jk[1][2]*d_jk[2][1] - d_jk[0][1]*d_jk[1][0]*d_jk[2][2] - d_jk[0][2]*d_jk[1][1]*d_jk[2][0];
	    
	if (fabs(fabs(val_V) - 6*val_volume[ind_ele])/fabs(val_V) > tol_zero) std::cout << "find different volume, ind_ele = " << ind_ele << ", volume = " << val_volume[ind_ele] << ", val_V = " << fabs(val_V/6) << ", dif = " << fabs(fabs(val_volume[ind_ele] - fabs(val_V/6))) << '\n';
	
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
		
	valuetype count = (valuetype) 0;
        for (unsigned int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
	    valuetype val_g_sol[3] = {0, 0, 0};
	    for (unsigned int ind = 0; ind < 3; ++ind)
		for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index)
		    for (unsigned int indt = 0; indt < 3; ++indt)
			val_g_sol[ind] += parxj_parxk[indt][ind] * val_coef_dif[ind_index] * basis_gradient_actual[p][ind_index][indt];
            // count the contribution
	    valuetype cnt_tmp = (valuetype) 0;
	    for (unsigned int ind = 0; ind < 3; ++ind)
		cnt_tmp += pow(val_g_sol[ind], 2);
            count += Weight[2][p] * cnt_tmp;
	}
	err += count * val_volume[ind_ele];
    }
    
    err = sqrt(err);
    return err;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_density_difference(Vector<valuetype> &u, std::vector<std::vector<valuetype> > &v)
{
    valuetype dif_l2 = 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	// recover local coefficient on this element
        std::vector<valuetype> val_coef_u(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef_u[ind_index] += u(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
       
	valuetype count = 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
            valuetype val_u = (valuetype) 0;
            for (int ind_index = 0; ind_index < n_index; ++ind_index)
                val_u += val_coef_u[ind_index] * basis_value_actual[p][ind_index];
            // count the contribution
            count += Weight[2][p] * pow(pow(val_u,2) - v[ind_ele][p], 2);
        }
        dif_l2 += count * val_volume[ind_ele];
    }
    dif_l2 = sqrt(dif_l2);
    
    return dif_l2;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_density_difference(Vector<valuetype> &u, Vector<valuetype> &v)
{
    valuetype dif_l2 = 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	// recover local coefficient on this element
        std::vector<valuetype> val_coef_u(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef_u[ind_index] += u(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        std::vector<valuetype> val_coef_v(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef_v[ind_index] += v(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	valuetype count = 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
            valuetype val_u = (valuetype) 0, val_v = (valuetype) 0;
            for (int ind_index = 0; ind_index < n_index; ++ind_index){
                val_u += val_coef_u[ind_index] * basis_value_actual[p][ind_index];
                val_v += val_coef_v[ind_index] * basis_value_actual[p][ind_index];
	    }
            // count the contribution
            count += Weight[2][p] * pow(pow(val_u,2) - pow(val_v,2), 2);
        }
        dif_l2 += count * val_volume[ind_ele];
    }
    dif_l2 = sqrt(dif_l2);
    
    return dif_l2;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_density_difference(std::vector<Vector<valuetype> > &psi, FEMSpace<double, 3> &fem_space, double(*psi_re)(double *, double), double(*psi_im)(double *, double), double t)
{
    valuetype dif_l2 = 0, count;
    std::vector<valuetype> val_psi_re(n_q_point[2]), val_psi_im(n_q_point[2]);
    for (unsigned int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(ind_ele).local_to_global(QPoint);
	calc_val_qp_onElement(psi[0], val_psi_re, ind_ele);
	calc_val_qp_onElement(psi[1], val_psi_im, ind_ele);
	count = 0;
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    count += Weight[2][p] * pow(pow(val_psi_re[p],2)        + pow(val_psi_im[p],2) -
					pow(psi_re(q_point[p],t),2) - pow(psi_im(q_point[p],t),2), 2);
        dif_l2 += count * val_volume[ind_ele];
    }
    dif_l2 = sqrt(dif_l2);
    
    return dif_l2;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_difference(Vector<valuetype> &u, Vector<valuetype> &v)
{
    valuetype dif_l2 = 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	// recover local coefficient on this element
        std::vector<valuetype> val_coef_u(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef_u[ind_index] += u(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        std::vector<valuetype> val_coef_v(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef_v[ind_index] += v(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	valuetype count = 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
            valuetype val_u = (valuetype) 0, val_v = (valuetype) 0;
            for (int ind_index = 0; ind_index < n_index; ++ind_index){
                val_u += val_coef_u[ind_index] * basis_value_actual[p][ind_index];
                val_v += val_coef_v[ind_index] * basis_value_actual[p][ind_index];
	    }
            // count the contribution
            count += Weight[2][p] * pow(val_u - val_v, 2);
        }
        dif_l2 += count * val_volume[ind_ele];
    }
    dif_l2 = sqrt(dif_l2);
    
    return dif_l2;
}


TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_difference(Vector<valuetype> &u, Vector<valuetype> &v, int polynomial_order, int flag)
{ // flag: 0 - each, 1 - <= polynomial order
    valuetype dif_l2 = 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	// recover local coefficient on this element
        std::vector<valuetype> val_coef_u(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef_u[ind_index] += u(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        std::vector<valuetype> val_coef_v(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef_v[ind_index] += v(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
	for (int ind_index = 0; ind_index < n_index; ++ind_index){
	    Multiindex<3> index_now = correspondence.number2index(ind_index);
	    if (flag == 0 && index_now.sum() == polynomial_order) continue;
	    if (flag == 1 && index_now.sum() <=  polynomial_order) continue;
	    val_coef_u[ind_index] = 0;
	    val_coef_v[ind_index] = 0;
	}

	valuetype count = 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
            valuetype val_u = (valuetype) 0, val_v = (valuetype) 0;
            for (int ind_index = 0; ind_index < n_index; ++ind_index){
                val_u += val_coef_u[ind_index] * basis_value_actual[p][ind_index];
                val_v += val_coef_v[ind_index] * basis_value_actual[p][ind_index];
	    }
            // count the contribution
            count += Weight[2][p] * pow(val_u - val_v, 2);
        }
        dif_l2 += count * val_volume[ind_ele];
    }
    dif_l2 = sqrt(dif_l2);
    
    return dif_l2;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_h1_difference(RegularMesh<3> &mesh, Vector<valuetype> &u, Vector<valuetype> &v)
{
    valuetype err_l2 = calc_l2_difference(u, v);
    valuetype err_l2_gradient = calc_l2_error_gradient_difference(mesh, u, v);
    return err_l2 + err_l2_gradient;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_l2_error_component(Vector<valuetype> &u, int polynomial_order_less, int polynomial_order_more)
{ // calcultate L2 error between polynomial_order_less and polynomial_order_large,
  // suppose polynomial_order_less < polynomial_order_more <= M
    valuetype err_l2 = 0;
    std::vector<bool> flag_less(n_index, true);
    std::vector<bool> flag_more(n_index, true);
    for (int ind_index = 0; ind_index < n_index; ++ind_index){
	Multiindex<3> index_now = correspondence.number2index(ind_index);
	if (index_now.sum() > polynomial_order_less) flag_less[ind_index] = false;
	if (index_now.sum() > polynomial_order_more) flag_more[ind_index] = false;
    }
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	// recover local coefficient on this element
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index)
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += u(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];

	valuetype count = 0;
        for (int p = 0; p < n_q_point[2]; ++p){
            // recover the value of solution
            valuetype val_less = (valuetype) 0, val_more = (valuetype) 0;
            for (int ind_index = 0; ind_index < n_index; ++ind_index){
		if (flag_less[ind_index]) val_less += val_coef[ind_index] * basis_value_actual[p][ind_index];
		if (flag_more[ind_index]) val_more += val_coef[ind_index] * basis_value_actual[p][ind_index];
	    }
            // count the contribution
            count += Weight[2][p] * pow(val_less - val_more, 2);
        }
        err_l2 += count * val_volume[ind_ele];
    }
    err_l2 = sqrt(err_l2);
    
    return err_l2;
}

TEMPLATE_TSEM
valuetype THIS_TSEM::calc_density_l2_difference(std::vector<Vector<valuetype> > &u, std::vector<Vector<valuetype> > &v, std::vector<valuetype> &n_occupation)
{
    unsigned int n_orbital = u.size();
    valuetype dif_l2 = 0, count;
    std::vector<valuetype> val_qp_u(n_q_point[2]);
    std::vector<valuetype> val_qp_v(n_q_point[2]);
    std::vector<valuetype> val_dif(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed

	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    val_dif[p] = 0;

	for (unsigned int ind_orbital = 0; ind_orbital < n_orbital; ++ind_orbital){
	    calc_val_qp_onElement(u[ind_orbital], val_qp_u, ind_ele);
	    calc_val_qp_onElement(v[ind_orbital], val_qp_v, ind_ele);
	    for (unsigned int p = 0; p < n_q_point[2]; ++p)
		val_dif[p] += n_occupation[ind_orbital] * (pow(val_qp_u[p],2) - pow(val_qp_v[p],2));
	}

	count = 0;
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    count += Weight[2][p] * pow(val_dif[p], 2);

	dif_l2 += count * val_volume[ind_ele];
    }
    dif_l2 = sqrt(dif_l2);
    
    return dif_l2;
}

#undef TEMPLATE_TSEM
#undef THIS_TSEM
