#include "../include/TetrahedralSEM.h"
// This file contains functions for TSEM discretization, including:
//     1. building mass matrix, mass_V and stiff matrix corresponding to TSEM discretization of (u_h,v_h), (u_h,v_h)_V, (\nabla u_h,\nabla v_h)
//     2. calculating right-hand-side
//     3. building boundary mark and imposing boundary condition
// Functions:
// build_flag_bm()
// build_mass_matrix_init()
// build_mass_matrix()
// build_mass_V_matrix()
// build_stiff_matrix_init()
// build_stiff_matrix()
// calc_rhs()
// impose_zero_boundary_condition()
// impose_zero_boundary_condition_complex()
// impose_boundary_condition_rowOnly()
// impose_boundary_condition()
// is_on_boundary()


#define TEMPLATE_TSEM template<typename valuetype>
#define THIS_TSEM TSEM<valuetype>


TEMPLATE_TSEM
void THIS_TSEM::build_flag_bm(valuetype bnd_left, valuetype bnd_right)
{ // build boundary mark for cubic region [Bnd_Left, Bnd_Right]^3, inner-0, outer-1

    flag_bm.resize(4); // boundary flag for all dimensional geometry
    flag_bm_dof.resize(n_dof_total, 0);
    for (int ind = 0; ind <= 3; ++ind){
	flag_bm[ind].resize(n_geometry[ind], 0);
        if (ind == 3) break;
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
            flag_bm[ind][ind_geo] = is_on_boundary(point_ref_mesh[ind][ind_geo], bnd_left, bnd_right) ? 1 : 0;
	    if (flag_bm[ind][ind_geo] == 0) continue;
	    int location = ind_geo + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
	    int row_start = location_actualdof[location_geometry[location]];
	    for (unsigned int ind_dof = 0; ind_dof < n_dof_geometry[ind]; ++ind_dof)
		flag_bm_dof[row_start + ind_dof] = 1;
	}
    }
    std::cerr << "build boundary flag for geometrys with boundary info: Bnd_Left = " << bnd_left << ", Bnd_Right = " << bnd_right << '\n';
}

TEMPLATE_TSEM
void THIS_TSEM::build_flag_bm(valuetype radius)
{ // build boundary mark for ball region centering at (0, 0, 0) with radius, inner-0, outer-1

    flag_bm.resize(4); // boundary flag for all dimensional geometry
    flag_bm_dof.resize(n_dof_total, 0);
    for (int ind = 0; ind <= 3; ++ind){
        flag_bm[ind].resize(n_geometry[ind], 0);
        if (ind == 3) break;
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
	    switch (ind){
	    case 0:
		flag_bm[0][ind_geo] = is_on_boundary(point_ref_mesh[0][ind_geo], radius) ? 1 : 0;
		break;
	    case 1:
		flag_bm[1][ind_geo] = flag_bm[0][number_node[0][ind_geo][0]]!=0 && flag_bm[0][number_node[0][ind_geo][1]]!=0 ? 1 : 0;
		break;
	    case 2:
		flag_bm[2][ind_geo] = flag_bm[0][number_node[1][ind_geo][0]]!=0 && flag_bm[0][number_node[1][ind_geo][1]]!=0 && flag_bm[0][number_node[1][ind_geo][2]]!=0 ? 1 : 0;
		break;
	    default:
		break;
	    }
	    if (flag_bm[ind][ind_geo] == 0) continue;
	    int location = ind_geo + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
	    int row_start = location_actualdof[location_geometry[location]];
	    for (unsigned int ind_dof = 0; ind_dof < n_dof_geometry[ind]; ++ind_dof)
		flag_bm_dof[row_start + ind_dof] = 1;
	}
    }
    std::cerr << "build boundary flag for geometrys with boundary info: center = (0, 0, 0), radius = " << radius << '\n';
}

TEMPLATE_TSEM
void THIS_TSEM::build_flag_bm(valuetype bnd_left, valuetype bnd_right, valuetype axis[])
{ // build boundary mark for cubic region [Bnd_Left, Bnd_Right]^3, inner-0, outer-1

    flag_bm.resize(4); // boundary flag for all dimensional geometry
    flag_bm_dof.resize(n_dof_total, 0);
    std::vector<AFEPack::Point<3> > axis_coord(3);
    for (unsigned int ind_axis = 0; ind_axis < 3; ++ind_axis)
	for (unsigned int ind = 0; ind < 3; ++ind)
	    axis_coord[ind_axis][ind] = axis[ind_axis*3 + ind];
    for (int ind = 0; ind <= 3; ++ind){
	flag_bm[ind].resize(n_geometry[ind], 0);
        if (ind == 3) break;
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
	    AFEPack::Point<3> point_now;
	    for (unsigned int ind_c = 0; ind_c < 3; ++ind_c)
		point_now[ind_c] = calc_inner_product(axis_coord[ind_c], point_ref_mesh[ind][ind_geo]);
            flag_bm[ind][ind_geo] = is_on_boundary(point_now, bnd_left, bnd_right) ? 1 : 0;
	    if (flag_bm[ind][ind_geo] == 0) continue;
	    int location = ind_geo + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
	    int row_start = location_actualdof[location_geometry[location]];
	    for (unsigned int ind_dof = 0; ind_dof < n_dof_geometry[ind]; ++ind_dof)
		flag_bm_dof[row_start + ind_dof] = 1;
	}
    }
    std::cerr << "build boundary flag for geometrys with boundary info: Bnd_Left = " << bnd_left << ", Bnd_Right = " << bnd_right << '\n';
}

TEMPLATE_TSEM
void THIS_TSEM::build_flag_bm(RegularMesh<3> &mesh)
{
    flag_bm.resize(4); // boundary flag for all dimensional geometry
    flag_bm_dof.resize(n_dof_total);
    for (unsigned int ind = 0; ind <= 3; ++ind){
	flag_bm[ind].resize(n_geometry[ind]);
	if (ind == 3) break;
	for (unsigned int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo)
	    flag_bm[ind][ind_geo] = mesh.geometry(ind, ind_geo).boundaryMark();
    }
}

TEMPLATE_TSEM
void THIS_TSEM::build_mass_matrix_init()
{
    // construct sparsity pattern for mass matrix element
    std::vector<unsigned int> n_nnz_per_row_mass(n_index, 5 * 7 * 9); // the largest one
					   SparsityPattern sp_mass_matrix_element(n_index, n_index, n_nnz_per_row_mass);
    for (int row = 0; row < n_index; ++row){
        Multiindex<3> index_row = correspondence.number2index(row);
        for (int d1 = -2; d1 <= 2; ++d1)
            for (int d2 = -3; d2 <= 3; ++d2)
                for (int d3 = -4; d3 <= 4; ++d3){
                    Multiindex<3> index_col = index_row + Unitary_Multiindex[0] * d1 + Unitary_Multiindex[1] * (d2-d1) + Unitary_Multiindex[2] * (d3-d2);
                    if (!(Zero_Multiindex <= index_col) || index_col.sum() > M) continue;
		    // std::cout << "index_col = [" << index_col.index[0] << ", " << index_col.index[1] << ", " << index_col.index[2] <<  "]\n";
                    int col = correspondence.index2number(index_col);
                    sp_mass_matrix_element.add(row, col);
                }
    }
    sp_mass_matrix_element.compress();

    // construct sparse matrix of stiff matrix element
    SparseMatrix<valuetype> mass_matrix_element(sp_mass_matrix_element);
    valuetype C[3]; // coefficient [C1, C2, C3]
    for (int row = 0; row < n_index; ++row){
        Multiindex<3> index_row = correspondence.number2index(row);
        int l1 = index_row.index[0], l2 = index_row.index[1], l3 = index_row.index[2];
        for (int d1 = -2; d1 <= 2; ++d1){
            int k1 = l1 + d1;
            if (k1 < 0) continue; // continue if multiindex has negative component
            C[0] = (valuetype) 0;
            if (l1 == 0)
                switch (k1){
                case 0:  C[0] = (valuetype) 2;       break;
                case 2:  C[0] = ((valuetype) -1)/3;  break;
                default: C[0] = (valuetype) 0;
                }
            else if (l1 == 1)
                switch (k1){
                case 1:  C[0] = ((valuetype) 2)/3;   break;
                case 3:  C[0] = ((valuetype) -2)/15; break;
                default: C[0] = (valuetype) 0;
                }
            else{
                valuetype sum = (valuetype) 0;
                for (int i1 = ((d1-1 > -1) ? d1-1 : -1); i1 <= ((d1+1 < 1) ? d1+1 : 1); ++i1)
                    sum += (((valuetype) calc_delta(i1, 0)) + calc_coefficient_a(-1, -1, i1+2, l1+i1)) * (((valuetype) calc_delta(d1, i1)) - calc_coefficient_a(-1, -1, d1-i1+2, l1+d1));
                C[0] = ((valuetype) (l1-1)) / ((2*l1-1)*2*l1) * sum;
            }
            for (int d2 = -3; d2 <= 3; ++d2){
                int d2a = d2 - d1; // actual d2
                int k2 = l2 + d2a;
                if (k2 < 0) continue; // continue if multiindex index has negative component
                int aph2 = 2 * l1 - 1;
                C[1] = (valuetype) 0;
                if (l1 == 0 && l2 == 0)
                    switch (d1){
                    case 0: switch (k2){
                        case 0: C[1] = (valuetype) 2; break;
                        case 1: C[1] = ((valuetype) -2)/3; break;
                        case 2: C[1] = ((valuetype) -1)/3; break;
                        case 3: C[1] = ((valuetype) 2)/15; break;
                        }
                        break;
                    case 1: switch (k2){
                        case 0: C[1] = ((valuetype) 8)/3; break;
                        case 1: C[1] = ((valuetype) 4)/3; break;
                        case 2: C[1] = ((valuetype) -2)/5; break;
                        }
                        break;
                    case 2: switch (k2){
                        case 0: C[1] = (valuetype) 4; break;
                        case 1: C[1] = ((valuetype) 16)/5; break;
                        }
                    }
                else if (l1 > 0 && l2 == 0){
                    valuetype numerator = (valuetype) pow(2, 2*l1+d1+2);
                    if (d2a != 0){
                        numerator *= (valuetype) (2*l1 + 2*d1 - 1 + d2a);
                        for (int i = 0; i < d2a-1; ++i)
                            numerator *= (valuetype) d1-2 + i;
                    }
                    valuetype denominator = (valuetype) 1;
                    for (int i = 0; i < d2a+1; ++i)
                        denominator *= (valuetype) 2*l1+d1+2 + i;
                    C[1] = numerator / denominator;
                }
                else if (l1 == 0 && k1 == 0 && l2 == 1)
                    switch (k2){
                    case 0: C[1] = ((valuetype) -2)/3; break;
                    case 1: C[1] = ((valuetype) 2)/3; break;
                    case 2: C[1] = ((valuetype) 1)/15; break;
                    case 3: C[1] = ((valuetype) -2)/15; break;
                    case 4: C[1] = ((valuetype) 2)/35;
                    }
                else if (l1 == 0 && k1 > 0 && l2 == 1)
                    for (int j = 0; j <= 1; ++j){
                        valuetype numerator = (valuetype) pow(2, d1+2+j);
                        if (j == 1) numerator *= (valuetype) -1;
                        if (k2 != 0){
                            numerator *= (valuetype) 2*d1 - 1 + k2;
                            for (int i = 0; i < k2-1; ++i)
                                numerator *= (valuetype) d1-2-j + i;
                        }
                        valuetype denominator = (valuetype) 1;
                        for (int i = 0; i < k2+1; ++i)
                            denominator *= (valuetype) d1+2+j + i;
                        C[1] += numerator / denominator;
                    }
                else{
                    int ds = d2a + d1;
                    for (int i21 = ((ds-2 > -1) ? ds-2 : -1); i21 <= ((ds+2 < 1) ? ds+2 : 1); ++i21){
                        valuetype tmp = ((valuetype) calc_delta(i21, 0)) + calc_coefficient_a(aph2, -1, i21+2, l2+i21);
                        for (int i22 = ((ds-i21-1 > -1) ? ds-i21-1 : -1); i22 <= ((ds-i21+1 < 1) ? ds-i21+1 : 1); ++i22){
                            int i2 = i21 + i22;
                            valuetype tmp2;
                            switch (d1){
                            case -2: tmp2 =                       calc_coefficient_c(aph2-2, -1, i22+2, l2+i2+1); break;
                            case 2:  tmp2 =                     4*calc_coefficient_g(aph2,   -1, i22+2, l2+i2-1); break;
                            default: tmp2 = (((valuetype) calc_delta(i22, 0)) - calc_coefficient_a(aph2,   -1, i22+2, l2+i2));
                            }
                            switch (d1){
                            case -2: tmp2 *=                       calc_coefficient_c(aph2-4, -1, d2a-i2,   k2); break;
                            case -1: tmp2 *=                       calc_coefficient_c(aph2-2, -1, d2a-i2+1, k2); break;
                            case 0:  tmp2 *= ((valuetype) calc_delta(d2a, i2)) - calc_coefficient_a(aph2,   -1, d2a-i2+2, k2); break;
                            case 1:  tmp2 *=                     4*calc_coefficient_g(aph2,   -1, d2a-i2+3, k2); break;
                            case 2:  tmp2 *=                     4*calc_coefficient_g(aph2+2, -1, d2a-i2+4, k2);
                            }
                            C[1] += tmp * tmp2;
                        }
                    }
                    C[1] *= ((valuetype) pow(2, 2*l1-1) * (2*l1 + l2 - 1)) / (l2 * (2*l1 + 2*l2 - 1)); // the gamma
                }
                if (l1 >= 1 && k2 == 0){
                    valuetype numerator = (valuetype) pow(2, l1+k1+2);
                    if (l2 != 0){
                        numerator *= (valuetype) 2*l1 + l2 - 1;
                        for (int i = 0; i < l2-1; ++i)
                            numerator *= (valuetype) l1-k1-2 + i;
                    }
                    valuetype denominator = (valuetype) 1;
                    for (int i = 0; i < l2+1; ++i)
                        denominator *= (valuetype) l1+k1+2 + i;
                    C[1] = numerator / denominator;
                }
                if (l1 >= 1 && k1 == 0 && k2 == 1){
                    C[1] = (valuetype) 0;
                    for (int j = 0; j <= 1; ++j){
                        valuetype numerator = (valuetype) pow(2, l1+k1+j+2);
                        if (j == 1) numerator *= (valuetype) -1;
                        if (l2 != 0){
                            numerator *= (valuetype) 2*l1 + l2 - 1;
                            for (int i = 0; i < l2-1; ++i)
                                numerator *= (valuetype) l1-k1-j-2 + i;
                        }
                        valuetype denominator = (valuetype) 1;
                        for (int i = 0; i < l2+1; ++i)
                            denominator *= (valuetype) l1+k1+j+2 + i;
                        C[1] += numerator / denominator;
                    }
                }
                for (int d3 = -4; d3 <= 4; ++d3){
                    int d = d1 + d2a, l = l1 + l2, k = k1 + k2; // d1 is exactly actual d1
                    int d3a = d3 - d; // actual d3
                    int k3 = l3 + d3a;
                    if (k3 < 0 || k1+k2+k3 > M) continue; // continue if exist negative component or the summation exceed M
                    int aph3 = 2 * l1 + 2 * l2 - 1;
                    C[2] = (valuetype) 0;
                    if (l1 == 0 && l2 == 0 && l3 == 0)
                        switch (d){
                        case 0: switch (k3){
                            case 0: C[2] = ((valuetype) 8)/3; break;
                            case 1: C[2] = ((valuetype) -4)/3; break;
                            case 2: C[2] = ((valuetype) -2)/5; break;
                            case 3: C[2] = ((valuetype) 4)/15; break;
                            case 4: C[2] = ((valuetype) -2)/35; break;
                            }
                            break;
                        case 1: switch (k3){
                            case 0: C[2] = (valuetype) 4; break;
                            case 1: C[2] = ((valuetype) 8)/5; break;
                            case 2: C[2] = ((valuetype) -4)/5; break;
                            case 3: C[2] = ((valuetype) 16)/105;
                            }
                            break;
                        case 2: switch (k3){
                            case 0: C[2] = ((valuetype) 32)/5; break;
                            case 1: C[2] = ((valuetype) 64)/15; break;
                            case 2: C[2] = ((valuetype) -16)/21;
                            }
                            break;
                        case 3: switch (k3){
                            case 0: C[2] = ((valuetype) 32)/3; break;
                            case 1: C[2] = ((valuetype) 64)/7;
                            }
                        }
                    else if (l > 0 && l3 == 0){
                        valuetype numerator = (valuetype) pow(2, 2*l+d+3);
                        if (d3a != 0){
                            numerator *= (valuetype) (2*l + 2*d - 1 + d3a);
                            for (int i = 0; i < d3a-1; ++i)
                                numerator *= (valuetype) d-3 + i;
                        }
                        valuetype denominator = (valuetype) 1;
                        for (int i = 0; i < d3a+1; ++i)
                            denominator *= (valuetype) 2*l+d+3 + i;
                        C[2] = numerator / denominator;
                    }
                    else if (l1 == 0 && l2 == 0 && k1 == 0 && k2 == 0 && l3 == 1)
                        switch (k3){
                        case 0: C[2] = ((valuetype) -4)/3; break;
                        case 1: C[2] = ((valuetype) 16)/15; break;
                        case 2: C[2] = ((valuetype) 2)/15; break;
                        case 3: C[2] = ((valuetype) -4)/21; break;
                        case 4: C[2] = ((valuetype) 4)/35; break;
                        case 5: C[2] = ((valuetype) -8)/315;
                        }
                    else if (l1 == 0 && l2 == 0 && d > 0 && l3 == 1)
                        for (int j = 0; j <= 1; ++j){
                            valuetype numerator = (valuetype) pow(2, d+3+j);
                            if (j == 1) numerator *= (valuetype) -1;
                            if (k3 != 0){
                                numerator *= (valuetype) 2*d - 1 + k3;
                                for (int i = 0; i < k3-1; ++i)
                                    numerator *= (valuetype) d-3-j + i;
                            }
                            valuetype denominator = (valuetype) 1;
                            for (int i = 0; i < k3+1; ++i)
                                denominator *= (valuetype) d+3+j + i;
                            C[2] += numerator / denominator;
                        }
                    else{
                        int ds = d3a + d;
                        for (int i31 = ((ds-3 > -1) ? ds-3 : -1); i31 <= ((ds+3 < 1) ? ds+3 : 1); ++i31){
                            valuetype tmp = ((valuetype) calc_delta(i31, 0)) + calc_coefficient_a(aph3, -1, i31+2, l3+i31);
                            for (int i32 = ((ds-i31-2 > -1) ? ds-i31-2 : -1); i32 <= ((ds-i31+2 < 1) ? ds-i31+2 : 1); ++i32){
                                valuetype tmp2;
                                switch (d){
                                case -3: tmp2 =                      calc_coefficient_c(aph3-2, -1, i32+2, l3+i31+i32+1); break;
                                case 3:  tmp2 =                    4*calc_coefficient_g(aph3,   -1, i32+2, l3+i31+i32-1); break;
                                default: tmp2 = ((valuetype) calc_delta(i32, 0)) - calc_coefficient_a(aph3,   -1, i32+2, l3+i31+i32);
                                }
                                for (int i33 = ((ds-i31-i32-1 > -1) ? ds-i31-i32-1 : -1); i33 <= ((ds-i31-i32+1 < 1) ? ds-i31-i32+1 : 1); ++i33){
                                    int i3 = i31 + i32 + i33;
                                    valuetype tmp3;
                                    switch (d){
                                    case -3: tmp3 =                      calc_coefficient_c(aph3-4, -1, i33+2, l3+i3+2); break;
                                    case -2: tmp3 =                      calc_coefficient_c(aph3-2, -1, i33+2, l3+i3+1); break;
                                    case 2:  tmp3 =                    4*calc_coefficient_g(aph3,   -1, i33+2, l3+i3-1); break;
                                    case 3:  tmp3 =                    4*calc_coefficient_g(aph3+2, -1, i33+2, l3+i3-2); break;
                                    default: tmp3 = ((valuetype) calc_delta(i33, 0)) - calc_coefficient_a(aph3,   -1, i33+2, l3+i3);
                                    }
                                    switch (d){
                                    case -3: tmp3 *=                       calc_coefficient_c(aph3-6, -1, d3a-i3-1, k3); break;
                                    case -2: tmp3 *=                       calc_coefficient_c(aph3-4, -1, d3a-i3,   k3); break;
                                    case -1: tmp3 *=                       calc_coefficient_c(aph3-2, -1, d3a-i3+1, k3); break;
                                    case 0:  tmp3 *= ((valuetype) calc_delta(d3a, i3)) - calc_coefficient_a(aph3,   -1, d3a-i3+2, k3); break;
                                    case 1:  tmp3 *=                     4*calc_coefficient_g(aph3,   -1, d3a-i3+3, k3); break;
                                    case 2:  tmp3 *=                     4*calc_coefficient_g(aph3+2, -1, d3a-i3+4, k3); break;
                                    case 3:  tmp3 *=                     4*calc_coefficient_g(aph3+4, -1, d3a-i3+5, k3);
                                    }
                                    C[2] += tmp * tmp2 * tmp3;
                                }
                            }
                        }
                        C[2] *= ((valuetype) pow(2, 2*l-1) * (2*l + l3 - 1)) / (l3 * (2*l + 2*l3 - 1)); // the gamma
                    }
                    if (l >= 1 && k3 == 0){
                        valuetype numerator = (valuetype) pow(2, l+k+3);
                        if (l3 != 0){
                            numerator *= (valuetype) 2*l + l3 - 1;
                            for (int i = 0; i < l3-1; ++i)
                                numerator *= (valuetype) l-k-3 + i;
                        }
                        valuetype denominator = (valuetype) 1;
                        for (int i = 0; i < l3+1; ++i)
                            denominator *= (valuetype) l+k+3 + i;
                        C[2] = numerator / denominator;
                    }
                    if (l >= 1 && k == 0 && k3 == 1){
                        C[2] = (valuetype) 0;
                        for (int j = 0; j <= 1; ++j){
                            valuetype numerator = (valuetype) pow(2, l+k+j+3);
                            if (j == 1) numerator *= (valuetype) -1;
                            if (l3 != 0){
                                numerator *= (valuetype) 2*l + l3 - 1;
                                for (int i = 0; i < l3-1; ++i)
                                    numerator *= (valuetype) l-k-j-3 + i;
                            }
                            valuetype denominator = (valuetype) 1;
                            for (int i = 0; i < l3+1; ++i)
                                denominator *= (valuetype) l+k+j+3 + i;
                            C[2] += numerator / denominator;
                        }
                    }
                    valuetype val = C[0] * C[1] * C[2];
                    if (fabs(val) < tol_zero) continue; // continue if one of coefficient is 0
                    val *= pow((valuetype) 0.5, 2*(l1+k1)+l2+k2+6);
                    Multiindex<3> index_col = Unitary_Multiindex[0] * k1 + Unitary_Multiindex[1] * k2 + Unitary_Multiindex[2] * k3;
                    int col = correspondence.index2number(index_col);
                    mass_matrix_element.add(row, col, val);
                }
            }
        }
    }
    std::cerr << "generate mass matrix element\n";

    // std::ofstream sparsity("./info_mass_element_sparsity.m");
    // std::vector<unsigned int> ind_actual(n_index, 0);
    // unsigned int n_inner_index = 1;
    // Multiindex<3> index_cmp;
    // index_cmp.index[0] = 2; index_cmp.index[1] = 1; index_cmp.index[2] = 1;
    // for (unsigned int ind_index = 0; ind_index < n_index; ++ind_index){
    // 	Multiindex<3> index_now = correspondence.number2index(ind_index);
    // 	if (index_cmp <= index_now){
    // 	    unsigned int l1 = index_now.index[0]-2, l2 = index_now.index[1]-1, l3 = index_now.index[2]-1;
    // 	    unsigned int tmp;
    // 	    tmp = ((l1-1)*l1*(2*l1-1) + l1*(M+1)*(M+2)*6 - (2*M+3)*(l1-1)*l1*3) / 12
    // 		+ (M-l1+1)*l2 - (l2-1)*l2/2
    // 		+ l3;
    // 	    ind_actual[ind_index] = tmp + 1;
    // 	}
    // }
    // SparseMatrix<double>::iterator spm_ite = mass_matrix_element.begin(0);
    // SparseMatrix<double>::iterator spm_end = mass_matrix_element.end(mass_matrix_element.m()-1);
    // sparsity << "row_massEle = [";
    // for (; spm_ite != spm_end; ++spm_ite)
    // 	if (ind_actual[spm_ite->row()] > 0 && ind_actual[spm_ite->column()] > 0)
    // 	    sparsity << ind_actual[spm_ite->row()] << ' ';
    // sparsity << "];\n";
    // spm_ite = mass_matrix_element.begin(0);
    // sparsity << "col_massEle = [";
    // for (; spm_ite != spm_end; ++spm_ite)
    // 	if (ind_actual[spm_ite->row()] > 0 && ind_actual[spm_ite->column()] > 0)
    // 	    sparsity << ind_actual[spm_ite->column()] << ' ';
    // sparsity << "];\n";
    // spm_ite = mass_matrix_element.begin(0);
    // sparsity << "val_massEle = [";
    // for (; spm_ite != spm_end; ++spm_ite)
    // 	if (ind_actual[spm_ite->row()] > 0 && ind_actual[spm_ite->column()] > 0)
    // 	    sparsity << spm_ite->value() << ' ';
    // sparsity << "];\n";
    // sparsity << "spm_massEle = sparse(row_massEle, col_massEle, val_massEle);";
    // sparsity.close();

    // assign stiff matrix actual with n_transform_local, transform_local, weight_transform_local
    std::vector<unsigned int> n_nnz_per_row_actual_mass(n_index, 5 * 7 * 9 * 16 * 16);
    //     the largest one, h^{(3,2)} * max(n_transform_local)^2 * max(n_index_variation)^2
    sp_mass_matrix_element_actual.reinit(n_index, n_index, n_nnz_per_row_actual_mass);
    // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = mass_matrix_element.begin(0);
    // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = mass_matrix_element.end(mass_matrix_element.m()-1);
    // for (; spm_ite != spm_end; ++spm_ite){
    // 	int row = spm_ite->row();
    // 	int col = spm_ite->column();
    const SparsityPattern &sp_pattern = mass_matrix_element.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    for (unsigned int row = 0; row < mass_matrix_element.m(); ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_row_tl = 0; ind_row_tl < n_transform_local[row]; ++ind_row_tl) // tr for transform local
		for (int ind_col_tl = 0; ind_col_tl < n_transform_local[col]; ++ind_col_tl)
		    sp_mass_matrix_element_actual.add(transform_local[row][ind_row_tl],
						      transform_local[col][ind_col_tl]);
	}
    sp_mass_matrix_element_actual.compress();

    mass_matrix_element_actual.reinit(sp_mass_matrix_element_actual);
    // spm_ite = mass_matrix_element.begin(0);
    // for (; spm_ite != spm_end; ++spm_ite){
    // 	int row = spm_ite->row();
    // 	int col = spm_ite->column();
    // 	valuetype val = spm_ite->value();
    for (unsigned int row = 0; row < mass_matrix_element.m(); ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    valuetype val = mass_matrix_element.global_entry(pos_whole);
	    for (int ind_row_tl = 0; ind_row_tl < n_transform_local[row]; ++ind_row_tl) // tr for transform local
		for (int ind_col_tl = 0; ind_col_tl < n_transform_local[col]; ++ind_col_tl)
		    mass_matrix_element_actual.add(transform_local[row][ind_row_tl],
						   transform_local[col][ind_col_tl],
						   val
						   * weight_transform_local[row][ind_row_tl] * weight_transform_local[col][ind_col_tl]);
	}
    std::cerr << "given transform_local, assign actual mass matrix element\n";

    
    // const SparsityPattern &sp_pattern_t = mass_matrix_element_actual.get_sparsity_pattern();
    // const std::size_t *row_start_t = sp_pattern_t.get_rowstart_indices();
    // const unsigned int *col_nums_t = sp_pattern_t.get_column_numbers();
    // std::ofstream output("./mass_element.m");
    // unsigned int index = 1;
    // for (unsigned int row = 0; row < mass_matrix_element_actual.m(); ++row)
    // 	for (unsigned int pos_whole = row_start_t[row]; pos_whole < row_start_t[row+1]; ++pos_whole){
    // 	    const unsigned int &col = col_nums_t[pos_whole];
    // 	    output << "row(" << index << ") = " << row << "; "
    // 		   << "col(" << index << ") = " << col << "; "
    // 		   << "val(" << index << ") = " << mass_matrix_element_actual.global_entry(pos_whole) << ";\n";
    // 	    index++;
    // 	}
}

TEMPLATE_TSEM
void THIS_TSEM::build_mass_matrix(FEMSpace<double, 3> &fem_space)
{
    // assemble mass matrix
    n_nonzero_per_row_mass.resize(n_dof_total, 0);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	const SparsityPattern &sp_pattern = mass_matrix_element_actual.get_sparsity_pattern();
	const std::size_t *row_start = sp_pattern.get_rowstart_indices();
	const unsigned int *col_nums = sp_pattern.get_column_numbers();
	for (unsigned int row = 0; row < mass_matrix_element_actual.m(); ++row){
	    int cnt = 0;
	    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		const unsigned int &col = col_nums[pos_whole];
        // int m_matrix = mass_matrix_element_actual.m();
        // for (int row = 0; row < m_matrix; ++row){
        //     SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = mass_matrix_element_actual.begin(row);
        //     SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = mass_matrix_element_actual.end(row);
		// for (; spm_ite != spm_end; ++spm_ite)
		    // cnt += transform_n_global2local[ind_ele][spm_ite->column()];
		cnt += transform_n_global2local[ind_ele][col];
	    }
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][row]; ++ind_nnz)
                n_nonzero_per_row_mass[transform_ind_global2local[ind_ele][row][ind_nnz]] += cnt;
        }
    }

    sp_mass_matrix.reinit(n_dof_total, n_dof_total, n_nonzero_per_row_mass);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // int m_matrix = mass_matrix_element_actual.m();
        // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = mass_matrix_element_actual.begin(0);
        // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = mass_matrix_element_actual.end(m_matrix-1);
        // for (; spm_ite != spm_end; ++spm_ite){
        //     int row = spm_ite->row();
        //     int col = spm_ite->column();
	const SparsityPattern &sp_pattern = mass_matrix_element_actual.get_sparsity_pattern();
	const std::size_t *row_start = sp_pattern.get_rowstart_indices();
	const unsigned int *col_nums = sp_pattern.get_column_numbers();
	for (unsigned int row = 0; row < mass_matrix_element_actual.m(); ++row)
	    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		const unsigned int &col = col_nums[pos_whole];
		for (int ind_nnz_row = 0; ind_nnz_row < transform_n_global2local[ind_ele][row]; ++ind_nnz_row)
		    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][col]; ++ind_nnz_col)
			sp_mass_matrix.add(transform_ind_global2local[ind_ele][row][ind_nnz_row],
					   transform_ind_global2local[ind_ele][col][ind_nnz_col]);
	    }
    }
    sp_mass_matrix.compress();
    
    mass_matrix.reinit(sp_mass_matrix);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        const std::vector<int>& element_dof = fem_space.element(ind_ele).dof();
        valuetype coef = 6 * val_volume[ind_ele]; // the actual volume is fabs(volume * jacobian)
        // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = mass_matrix_element_actual.begin(0);
        // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = mass_matrix_element_actual.end(mass_matrix_element_actual.m()-1);
        // for (; spm_ite != spm_end; ++spm_ite){
        //     int row = spm_ite->row();
        //     int col = spm_ite->column();
        //     valuetype val = spm_ite->value();
	const SparsityPattern &sp_pattern = mass_matrix_element_actual.get_sparsity_pattern();
	const std::size_t *row_start = sp_pattern.get_rowstart_indices();
	const unsigned int *col_nums = sp_pattern.get_column_numbers();
	for (unsigned int row = 0; row < mass_matrix_element_actual.m(); ++row)
	    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		const unsigned int &col = col_nums[pos_whole];
		valuetype val = mass_matrix_element_actual.global_entry(pos_whole);
		for (int ind_nnz_row = 0; ind_nnz_row < transform_n_global2local[ind_ele][row]; ++ind_nnz_row)
		    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][col]; ++ind_nnz_col)
			mass_matrix.add(transform_ind_global2local[ind_ele][row][ind_nnz_row],
					transform_ind_global2local[ind_ele][col][ind_nnz_col],
					val * coef
					* transform_val_global2local[ind_ele][row][ind_nnz_row] * transform_val_global2local[ind_ele][col][ind_nnz_col]);
	    }
    }
    std::cerr << "construct mass matrix\n";

    
    // // assemble essential mass matrix
    // std::vector<unsigned int> n_nonzero_per_row_mass_essential;
    // n_nonzero_per_row_mass_essential.resize(n_dof_total, 0);
    // for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
    // 	const SparsityPattern &sp_pattern = mass_matrix_element_actual.get_sparsity_pattern();
    // 	const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    // 	const unsigned int *col_nums = sp_pattern.get_column_numbers();
    // 	for (unsigned int row = 0; row < mass_matrix_element_actual.m(); ++row){
    // 	    int cnt = 0;
    // 	    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
    // 		const unsigned int &col = col_nums[pos_whole];
    // 		// cnt += transform_n_global2local[ind_ele][col];
    // 		cnt += transform_n_local2global[ind_ele][col];
    // 	    }
    //         for (int ind_nnz = 0; ind_nnz < transform_n_local2global[ind_ele][row]; ++ind_nnz)
    //             n_nonzero_per_row_mass_essential[transform_ind_local2global[ind_ele][row][ind_nnz]] += cnt;
    //     }
    // }

    // sp_mass_matrix_essential.reinit(n_dof_total, n_dof_total, n_nonzero_per_row_mass_essential);
    // for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
    // 	const SparsityPattern &sp_pattern = mass_matrix_element_actual.get_sparsity_pattern();
    // 	const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    // 	const unsigned int *col_nums = sp_pattern.get_column_numbers();
    // 	for (unsigned int row = 0; row < mass_matrix_element_actual.m(); ++row)
    // 	    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
    // 		const unsigned int &col = col_nums[pos_whole];
    // 		// for (int ind_nnz_row = 0; ind_nnz_row < transform_n_local2global[ind_ele][row]; ++ind_nnz_row)
    // 		//     for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][col]; ++ind_nnz_col)
    // 		// 	sp_mass_matrix_essential.add(transform_ind_local2global[ind_ele][row][ind_nnz_row],
    // 		// 				     transform_ind_global2local[ind_ele][col][ind_nnz_col]);
    // 		for (int ind_nnz_row = 0; ind_nnz_row < transform_n_local2global[ind_ele][row]; ++ind_nnz_row)
    // 		    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_local2global[ind_ele][col]; ++ind_nnz_col)
    // 			sp_mass_matrix_essential.add(transform_ind_local2global[ind_ele][row][ind_nnz_row],
    // 						     transform_ind_local2global[ind_ele][col][ind_nnz_col]);
    // 	    }
    // }
    // sp_mass_matrix_essential.compress();
    
    // mass_matrix_essential.reinit(sp_mass_matrix_essential);
    // for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
    //     valuetype volume = fem_space.element(ind_ele).templateElement().volume();
    //     AFEPack::Point<3> point_tmp;
    //     for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
    //     valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed
    //     const std::vector<int>& element_dof = fem_space.element(ind_ele).dof();
    //     valuetype coef = 6 * fabs(volume * jacobian); // the actual volume is fabs(volume * jacobian)
    // 	const SparsityPattern &sp_pattern = mass_matrix_element_actual.get_sparsity_pattern();
    // 	const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    // 	const unsigned int *col_nums = sp_pattern.get_column_numbers();
    // 	for (unsigned int row = 0; row < mass_matrix_element_actual.m(); ++row)
    // 	    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
    // 		const unsigned int &col = col_nums[pos_whole];
    // 		valuetype val = mass_matrix_element_actual.global_entry(pos_whole);
    // 		// for (int ind_nnz_row = 0; ind_nnz_row < transform_n_local2global[ind_ele][row]; ++ind_nnz_row)
    // 		//     for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][col]; ++ind_nnz_col)
    // 		// 	mass_matrix_essential.add(transform_ind_local2global[ind_ele][row][ind_nnz_row],
    // 		// 				  transform_ind_global2local[ind_ele][col][ind_nnz_col],
    // 		// 				  val * coef
    // 		// 				  * transform_val_local2global[ind_ele][row][ind_nnz_row] * transform_val_global2local[ind_ele][col][ind_nnz_col]);
    // 		for (int ind_nnz_row = 0; ind_nnz_row < transform_n_local2global[ind_ele][row]; ++ind_nnz_row)
    // 		    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_local2global[ind_ele][col]; ++ind_nnz_col)
    // 			mass_matrix_essential.add(transform_ind_local2global[ind_ele][row][ind_nnz_row],
    // 						  transform_ind_local2global[ind_ele][col][ind_nnz_col],
    // 						  val * coef
    // 						  * transform_val_local2global[ind_ele][row][ind_nnz_row] * transform_val_local2global[ind_ele][col][ind_nnz_col]);
    // 	    }
    // }
}

TEMPLATE_TSEM
void THIS_TSEM::build_mass_V_matrix(std::vector<unsigned int> &n_nonzero_per_row_mass_V, SparsityPattern &sp_pattern, SparseMatrix<valuetype> &sp_matrix,
				    std::vector<std::vector<valuetype> > &value_V)
{ // function V is define in main cpp file
    
    std::vector<std::vector<std::vector<int> > > col_local_contribution;
    std::vector<std::vector<std::vector<valuetype> > > val_local_contribution;
    std::vector<std::vector<int> > n_local_contribution;
    col_local_contribution.resize(n_element);
    val_local_contribution.resize(n_element);
    n_local_contribution.resize(n_element);
    for (int i = 0; i < n_element; ++i){
        col_local_contribution[i].resize(n_index);
        val_local_contribution[i].resize(n_index);
        n_local_contribution[i].resize(n_index);
    }
    // std::cout << "setup for calculation of mass_V\n";

    std::vector<int>       col_nonzero_entry(n_index);
    std::vector<valuetype> val_nonzero_entry(n_index);
    int n_nonzero_entry;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	// valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
	// valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp); // the determinant of jacobian is fixed
        for (int ind_index_row = 0; ind_index_row < n_index; ++ind_index_row){
            n_nonzero_entry = 0;
            for (int ind_index_col = 0; ind_index_col < n_index; ++ind_index_col){
                valuetype count = 0; // the integral between row index and col index
                for (int p = 0; p < n_q_point[2]; ++p)
                    count += Weight[2][p] * basis_value_actual[p][ind_index_row] * basis_value_actual[p][ind_index_col] * value_V[ind_ele][p];
                if (fabs(count) > tol_zero){
                    col_nonzero_entry[n_nonzero_entry] = ind_index_col;
                    val_nonzero_entry[n_nonzero_entry] = count * val_volume[ind_ele];
                    n_nonzero_entry++;
                }
            }
            col_local_contribution[ind_ele][ind_index_row].resize(n_nonzero_entry);
            val_local_contribution[ind_ele][ind_index_row].resize(n_nonzero_entry);
            for (int j = 0; j < n_nonzero_entry; ++j){
                col_local_contribution[ind_ele][ind_index_row][j] = col_nonzero_entry[j];
                val_local_contribution[ind_ele][ind_index_row][j] = val_nonzero_entry[j];
            }
            n_local_contribution[ind_ele][ind_index_row] = n_nonzero_entry;
        }
    }
    // std::cout << "calculate nonzero contribution on each elements for mass_V\n";
    
    n_nonzero_per_row_mass_V.resize(n_dof_total, 0);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele)
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            // count nonzero contribution for each col
            int cnt = 0;
            for (int ind_nnz = 0; ind_nnz < n_local_contribution[ind_ele][ind_index]; ++ind_nnz)
                cnt += transform_n_global2local[ind_ele][col_local_contribution[ind_ele][ind_index][ind_nnz]];
            // add number to nonzero array
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                n_nonzero_per_row_mass_V[transform_ind_global2local[ind_ele][ind_index][ind_nnz]] += cnt;
        }
    
    sp_pattern.reinit(n_dof_total, n_dof_total, n_nonzero_per_row_mass_V);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele)
        for (int ind_index_row = 0; ind_index_row < n_index; ++ind_index_row)
            for (int ind_nnz_contribution = 0; ind_nnz_contribution < n_local_contribution[ind_ele][ind_index_row]; ++ind_nnz_contribution){
                int& ind_index_col = col_local_contribution[ind_ele][ind_index_row][ind_nnz_contribution];
                for (int ind_nnz_row = 0; ind_nnz_row < transform_n_global2local[ind_ele][ind_index_row]; ++ind_nnz_row)
                    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][ind_index_col]; ++ind_nnz_col)
                        sp_pattern.add(transform_ind_global2local[ind_ele][ind_index_row][ind_nnz_row],
				       transform_ind_global2local[ind_ele][ind_index_col][ind_nnz_col]);
            }
    sp_pattern.compress();
    
    sp_matrix.reinit(sp_pattern);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele)
        for (int ind_index_row = 0; ind_index_row < n_index; ++ind_index_row)
            for (int ind_nnz_contribution = 0; ind_nnz_contribution < n_local_contribution[ind_ele][ind_index_row]; ++ind_nnz_contribution){
                int& ind_index_col = col_local_contribution[ind_ele][ind_index_row][ind_nnz_contribution];
                for (int ind_nnz_row = 0; ind_nnz_row < transform_n_global2local[ind_ele][ind_index_row]; ++ind_nnz_row)
                    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][ind_index_col]; ++ind_nnz_col)
                        sp_matrix.add(transform_ind_global2local[ind_ele][ind_index_row][ind_nnz_row],
				      transform_ind_global2local[ind_ele][ind_index_col][ind_nnz_col],
				      val_local_contribution[ind_ele][ind_index_row][ind_nnz_contribution]
				      * transform_val_global2local[ind_ele][ind_index_row][ind_nnz_row]
				      * transform_val_global2local[ind_ele][ind_index_col][ind_nnz_col]);
            }    
    std::cerr << "construct mass type matrix for V\n";
}

TEMPLATE_TSEM
void THIS_TSEM::build_mass_V_matrix(std::vector<unsigned int> &n_nonzero_per_row_mass_V, SparsityPattern &sp_pattern, SparseMatrix<valuetype> &sp_matrix,
				    FEMSpace<double, 3> &fem_space, valuetype(*func)(double*))
{ // function V is define in main cpp file
    
    std::vector<std::vector<std::vector<int> > > col_local_contribution;
    std::vector<std::vector<std::vector<valuetype> > > val_local_contribution;
    std::vector<std::vector<int> > n_local_contribution;
    col_local_contribution.resize(n_element);
    val_local_contribution.resize(n_element);
    n_local_contribution.resize(n_element);
    for (int i = 0; i < n_element; ++i){
        col_local_contribution[i].resize(n_index);
        val_local_contribution[i].resize(n_index);
        n_local_contribution[i].resize(n_index);
    }
    // std::cout << "setup for calculation of mass_V\n";

    std::vector<int>       col_nonzero_entry(n_index);
    std::vector<valuetype> val_nonzero_entry(n_index);
    int n_nonzero_entry;
    std::vector<double> value_V(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(ind_ele).local_to_global(QPoint);
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    value_V[p] = func(q_point[p]);
        for (int ind_index_row = 0; ind_index_row < n_index; ++ind_index_row){
            n_nonzero_entry = 0;
            for (int ind_index_col = 0; ind_index_col < n_index; ++ind_index_col){
                valuetype count = 0; // the integral between row index and col index
                for (int p = 0; p < n_q_point[2]; ++p)
                    count += Weight[2][p] * basis_value_actual[p][ind_index_row] * basis_value_actual[p][ind_index_col] * value_V[p];
                if (fabs(count) > tol_zero){
                    col_nonzero_entry[n_nonzero_entry] = ind_index_col;
                    val_nonzero_entry[n_nonzero_entry] = count * val_volume[ind_ele];
                    n_nonzero_entry++;
                }
            }
            col_local_contribution[ind_ele][ind_index_row].resize(n_nonzero_entry);
            val_local_contribution[ind_ele][ind_index_row].resize(n_nonzero_entry);
            for (int j = 0; j < n_nonzero_entry; ++j){
                col_local_contribution[ind_ele][ind_index_row][j] = col_nonzero_entry[j];
                val_local_contribution[ind_ele][ind_index_row][j] = val_nonzero_entry[j];
            }
            n_local_contribution[ind_ele][ind_index_row] = n_nonzero_entry;
        }
    }
    // std::cout << "calculate nonzero contribution on each elements for mass_V\n";
    
    n_nonzero_per_row_mass_V.resize(n_dof_total, 0);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele)
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            // count nonzero contribution for each col
            int cnt = 0;
            for (int ind_nnz = 0; ind_nnz < n_local_contribution[ind_ele][ind_index]; ++ind_nnz)
                cnt += transform_n_global2local[ind_ele][col_local_contribution[ind_ele][ind_index][ind_nnz]];
            // add number to nonzero array
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                n_nonzero_per_row_mass_V[transform_ind_global2local[ind_ele][ind_index][ind_nnz]] += cnt;
        }
    
    sp_pattern.reinit(n_dof_total, n_dof_total, n_nonzero_per_row_mass_V);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele)
        for (int ind_index_row = 0; ind_index_row < n_index; ++ind_index_row)
            for (int ind_nnz_contribution = 0; ind_nnz_contribution < n_local_contribution[ind_ele][ind_index_row]; ++ind_nnz_contribution){
                int& ind_index_col = col_local_contribution[ind_ele][ind_index_row][ind_nnz_contribution];
                for (int ind_nnz_row = 0; ind_nnz_row < transform_n_global2local[ind_ele][ind_index_row]; ++ind_nnz_row)
                    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][ind_index_col]; ++ind_nnz_col)
                        sp_pattern.add(transform_ind_global2local[ind_ele][ind_index_row][ind_nnz_row],
				       transform_ind_global2local[ind_ele][ind_index_col][ind_nnz_col]);
            }
    sp_pattern.compress();
    
    sp_matrix.reinit(sp_pattern);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele)
        for (int ind_index_row = 0; ind_index_row < n_index; ++ind_index_row)
            for (int ind_nnz_contribution = 0; ind_nnz_contribution < n_local_contribution[ind_ele][ind_index_row]; ++ind_nnz_contribution){
                int& ind_index_col = col_local_contribution[ind_ele][ind_index_row][ind_nnz_contribution];
                for (int ind_nnz_row = 0; ind_nnz_row < transform_n_global2local[ind_ele][ind_index_row]; ++ind_nnz_row)
                    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][ind_index_col]; ++ind_nnz_col)
                        sp_matrix.add(transform_ind_global2local[ind_ele][ind_index_row][ind_nnz_row],
				      transform_ind_global2local[ind_ele][ind_index_col][ind_nnz_col],
				      val_local_contribution[ind_ele][ind_index_row][ind_nnz_contribution]
				      * transform_val_global2local[ind_ele][ind_index_row][ind_nnz_row]
				      * transform_val_global2local[ind_ele][ind_index_col][ind_nnz_col]);
            }    
    std::cerr << "construct mass type matrix for V\n";
}
TEMPLATE_TSEM
void THIS_TSEM::build_stiff_matrix_init()
{ // prepare coefficients and sparse matrix for assignment of stiff matrix

    // construct sparsity pattern for stiff matrix element
    int n_index_stiff = correspondence.n_index_end(M-1); // as the order of polynomial for stiff matrix is M - 1
    // std::cout << "stiff matrix element has size " << n_index_stiff << '\n';
    std::vector<unsigned int> n_nnz_per_row(n_index_stiff, 5 * 5 * 5); // the largest one, h^{(3,2)}
    std::vector<SparsityPattern> sp_stiff_matrix_element(6);
    for (int i = 0; i < 6; ++i)
        sp_stiff_matrix_element[i].reinit(n_index_stiff, n_index_stiff, n_nnz_per_row);
    for (int row = 0; row < n_index_stiff; ++row){
        Multiindex<3> index_now = correspondence.number2index(row);
        for (int d1 = -2; d1 <= 2; ++d1)
            for (int d2 = -2; d2 <= 2; ++d2)
                for (int d3 = -2; d3 <= 2; ++d3){
                    Multiindex<3> index_tmp = index_now + Unitary_Multiindex[0] * d1 + Unitary_Multiindex[1] * (d2-d1) + Unitary_Multiindex[2] * (d3-d2);
		    // std::cout << "index_tmp = [" << index_tmp.index[0] << ", " << index_tmp.index[1] << ", " << index_tmp.index[2] << "]\n";
                    if (!(Zero_Multiindex <= index_tmp) || index_tmp.sum() >= M) continue;
                    int col = correspondence.index2number(index_tmp);
		    // std::cout << "index_tmp = [" << index_tmp.index[0] << ", " << index_tmp.index[1] << ", " << index_tmp.index[2] <<  "]\n";
		    // if (col >= 20) std::cout << "find col = " << col << " exceeds size, index_tmp = [" << index_tmp.index[0] << ", " << index_tmp.index[1] << ", " << index_tmp.index[2] << "], sum = " << index_tmp.sum() << "\n";
                    // the sparsity pattern of stiff matrix element (1, 1), which corresponds to multiindex (2, 0, 0)
                    if (d1 == 0 && abs(d2) <= 1)
                        sp_stiff_matrix_element[0].add(row, col);
                    // the sparsity pattern of stiff matrix element (2, 1), which corresponds to multiindex (1, 1, 0)
                    if (abs(d1) <= 1 && abs(d2) <= 1)
                        sp_stiff_matrix_element[1].add(row, col);
                    // the sparsity pattern of stiff matrix element (3, 1), which corresponds to multiindex (1, 0, 1)
                    if (abs(d1) <= 1)
                        sp_stiff_matrix_element[2].add(row, col);
                    // the sparsity pattern of stiff matrix element (2, 2), which corresponds to multiindex (0, 2, 0)
                    if (abs(d1) <= 1 && abs(d2) <= 1)
                        sp_stiff_matrix_element[3].add(row, col);
                    // the sparsity pattern of stiff matrix element (3, 2), which corresponds to multiindex (0, 1, 1)
                    sp_stiff_matrix_element[4].add(row, col);
                    // the sparsity pattern of stiff matrix element (3, 3), which corresponds to multiindex (0, 0, 2)
                    if (abs(d1) <= 1)
                        sp_stiff_matrix_element[5].add(row, col);
                }
    }
    for (int i = 0; i < 6; ++i)
        sp_stiff_matrix_element[i].compress();
    
    // construct sparse matrix of stiff matrix element
    std::vector<SparseMatrix<valuetype> > stiff_matrix_element(6);
    for (int i = 0; i < 6; ++i)
        stiff_matrix_element[i].reinit(sp_stiff_matrix_element[i]);
    valuetype C1[6]; // coefficient C1
    valuetype C2[6]; // coefficient C2 for each choice of d1
    valuetype C3[6]; // coefficient C3 for each choice of d1+d2a = d2
    for (int row = 0; row < n_index_stiff; ++row){
        Multiindex<3> index_row = correspondence.number2index(row);
        int l1 = index_row.index[0], l2 = index_row.index[1], l3 = index_row.index[2];
        for (int d1 = -2; d1 <= 2; ++d1){
            // initialize
            for (int i = 0; i < 6; ++i)
                C1[i] = (valuetype) 0;
            // continue if multiindex has negative component
            int k1 = l1 + d1;
            if (k1 < 0) continue;
            // corresponds to (1, 1), or multiindex (2, 0, 0)
            if (d1 == 0) C1[0] = ((valuetype) 2) / (2*l1 + 1); // d1 = 0
            if (l1 == 0){
                if (0 <= d1 && d1 <= 1){
                    // corresponds to (2, 1), or multiindex (1, 1, 0)
                    C1[1] = 2 - 3 * d1; // d1 = 0: 1
                    // corresponds to (3, 1), or multiindex (1, 0, 1)
                    C1[2] = 2 - 3 * d1; // d1 = 0: 1
                    // corresponds to (2, 2), or multiindex (0, 2, 0)
                    C1[3] = 2 - d1; // d1 = 0: 1
                    // corresponds to (3, 3), or multiindex (0, 0, 2)
                    C1[5] = 2 - d1; // d1 = 0: 1
                }
                // corresponds to (3, 2), or mulltindex (0, 1, 1)
                switch (d1){ // d1 = 0: 2
                case 0: C1[4] = (valuetype) 2; break;
                case 2: C1[4] = ((valuetype) -1)/3; break;
                default: C1[4] = (valuetype) 0;
                }
            }
            else{
                if (abs(d1) <= 1){
                    // corresponds to (2, 1), or multiindex (1, 1, 0)
                    C1[1] = ((valuetype) 1) / (2*l1) * (calc_delta(d1, 0) - calc_coefficient_a(-1, 0, d1+2, l1+d1)); // d1 = -1: 1
                    // corresponds to (3, 1), or multiindex (1, 0, 1)
                    C1[2] = ((valuetype) 1) / (2*l1) * (calc_delta(d1, 0) - calc_coefficient_a(-1, 0, d1+2, l1+d1)); // d1 = -1: 1
                    // corresponds to (2, 2), or multiindex (0, 2, 0)
                    C1[3] = ((valuetype) 1) / (2*l1) * (calc_delta(d1, 0) + calc_coefficient_a(0, -1, d1+2, l1+d1)); // d1 = -1: 1
                    // correspondes to (3, 3), or multiindex (0, 0, 2)
                    C1[5] = ((valuetype) 1) / (2*l1) * (calc_delta(d1, 0) + calc_coefficient_a(0, -1, d1+2, l1+d1)); // d1 = -1: 1
                }
                // corresponds to (3, 2), or multiindex (0, 1, 1)
                if (l1 == 1)
                    switch (k1){
                    case 1: C1[4] = ((valuetype) 2)/3; break;
                    case 3: C1[4] = ((valuetype) -2)/15; break;
                    default: C1[4] = (valuetype) 0;
                    }
                else{
                    valuetype sum = 0;
                    for (int i1 = ((d1-1 > -1) ? d1-1 : -1);
                         i1 <= ((d1+1 < 1) ? d1+1 : 1); ++i1)
                        sum += (calc_delta(i1, 0) + calc_coefficient_a(-1, -1, i1+2, l1+i1)) * (calc_delta(d1, i1) - calc_coefficient_a(-1, -1, d1-i1+2, l1+d1));
                    C1[4] = ((valuetype) (l1-1)) / ((2*l1-1)*2*l1) * sum; // d1 = -2: 2
                }
            }
            for (int d2 = -2; d2 <= 2; ++d2){
                int d2a = d2 - d1; // actual d2
                // continue if multiindex index has negative component
                int k2 = l2 + d2a;
                if (k2 < 0) continue;
                int aph2 = 2 * l1;
                for (int i = 0; i < 6; ++i) // initialize
                    C2[i] = (valuetype) 0;
                if (abs(d1) <= 1 && abs(d2) <= 1){
                    // corresponds to (2, 1), or multiindex (1, 1, 0)
                    C2[1] = ((valuetype) pow(2, aph2+1)) / (aph2+2*l2+1);
                    switch (d1){ // d2 = -1: 1
                    case -1: C2[1] *= calc_coefficient_c(aph2-2, 0, d2a+1, l2+d2a);
                        break;
                    case  0: C2[1] *= ((valuetype) calc_delta(d2a, 0)) - calc_coefficient_a(aph2, 0, d2a+2, l2+d2a);
                        break;
                    case  1: C2[1] *= 4 * calc_coefficient_g(aph2, 0, d2a+3, l2+d2a);
		    }
                    // corresponds to (2, 2), or multiindex (0, 2, 0)
                    C2[3] = C2[1];
		}
                if (l2 == 0){
                    if (d1 == 0 && 0 <= d2a && d2a <= 1)
                        // corresponds to (1, 1), or multiindex (2, 0, 0)
                        C2[0] = ((valuetype) pow(2, aph2+2)) / (aph2+2+d2a); // d1 = 0, d2a = 0: 1
                    if (abs(d1) <= 1 && 0 <= d2a && d2a <= 3){
                        // corresponds to (3, 1), or multiindex (1, 0, 1)
                        C2[2] = (valuetype) pow(2, aph2+d1+2);
                        valuetype numerator = 1;
                        if (d2a != 0){
                            for (int i = 0; i < d2a-1; ++i)
                                numerator *= (valuetype) d1-1 + i;
                            numerator *= (valuetype) aph2 + 2*d1 + d2a;
                        }
                        valuetype denominator = 1;
                        for (int i = 0; i < d2a+1; ++i)
                            denominator *= (valuetype) aph2+d1+2 + i;
                        C2[2] *= numerator / denominator; // d2a = 0: 3
                        // corresponds to (3, 3), or multiindex (0, 0, 2)
                        C2[5] = C2[2];// d2a = 0: 3
                    }
                }
                else{
                    if (d1 == 0 && abs(d2) <= 1)
                        // corresponds to (1, 1), or multiindex (2, 0, 0)
                        C2[0] = ((valuetype) pow(2, aph2+1) * (aph2+l2+1)) / (l2 * (aph2+2*l2+1))
                            * (((valuetype) calc_delta(d2a, 0)) + calc_coefficient_a(aph2+1, -1, d2a+2, l2+d2a)); // d1 = 0, d2 = -1: 1
                    if (abs(d1) <= 1){
                        // corresponds to (3, 1), or multiindex (1, 0, 1)
                        C2[2] = ((valuetype) pow(2, aph2) * (aph2+l2)) / (l2 * (aph2+2*l2));
                        valuetype sum = 0;
                        for (int i2 = ((d2a-1 + d1 > -1) ? d2a-1 + d1 : -1);
                             i2 <= ((d2a+1 + d1 < 1) ? d2a+1 + d1 : 1); ++i2){
                            // if (l2 + i2 < 0) continue;
                            valuetype tmp = ((valuetype) calc_delta(i2, 0)) + calc_coefficient_a(aph2, -1, i2+2, l2+i2);
                            switch (d1){
                            case -1: sum += tmp * calc_coefficient_c(aph2-2, -1, d2a-i2+1, l2+d2a);
                                break;
                            case  0: sum += tmp * (((valuetype) calc_delta(d2a, i2)) - calc_coefficient_a(aph2, -1, d2a-i2+2, l2+d2a));
                                break;
                            case  1: sum += tmp * 4 * calc_coefficient_g(aph2, -1, d2a-i2+3, l2+d2a);
                            }
                        }
                        C2[2] *= sum; // d1 = -1: 1, d2 = -2: 2
                        // corresponds to (3, 3), or multiindex (0, 0, 2)
                        C2[5] = C2[2]; // d1 = -1: 1, d2 = -2: 2
                    }
                }
                // corresponds to (3, 2), or multiindex (0, 1, 1)
                if (l1 == 0 && l2 == 0)
                    switch (d1){
                    case 0:
                        switch (d2a){
                        case 0: C2[4] = (valuetype) 2; break;
                        case 1: C2[4] = ((valuetype) -4)/3; break;
                        case 2: C2[4] = ((valuetype) 1)/3; break;
                        default: C2[4] = (valuetype) 0;
                        }
                        break;
                    case 1:
                        switch (d2a){
                        case 0: C2[4] = ((valuetype) 8)/3; break;
                        case 1: C2[4] = ((valuetype) -2)/3; break;
                        default: C2[4] = (valuetype) 0;
                        }
                        break;
                    case 2:
                        if (d2a == 0) C2[4] = (valuetype) 4;
                        break;
                    default: C2[4] = (valuetype) 0;
                    }
                else{
                    C2[4] = ((valuetype) pow(2, aph2)) / (aph2+2*l2);
                    valuetype sum = (valuetype) 0;
                    if (d1 == -2)
                        for (int i2 = ((d2a-2 > 0) ? d2a-2 : 0);
                             i2 <= ((d2a < 2) ? d2a : 2); ++i2){
                            sum += calc_coefficient_c(aph2-3, 0, i2+1, l2+i2) * calc_coefficient_c(aph2-5, 0, d2a-i2+1, l2+d2a);
                        }
                    else if (d1 == 2)
                        for (int i2 = ((d2a > -2) ? d2a : -2);
                             i2 <= ((d2a+2 < 0) ? d2a+2 : 0); ++i2){
                            sum += 16 * calc_coefficient_g(aph2-1, 0, i2+3, l2+i2) * calc_coefficient_g(aph2+1, 0, d2a-i2+3, l2+d2a);
                        }
                    else
                        for (int i2 = ((d2a-1 + d1 > -1) ? d2a-1 + d1 : -1);
                             i2 <= ((d2a+1 + d1 < 1) ? d2a+1 + d1 : 1); ++i2){
                            valuetype tmp = 1-abs(i2) - calc_coefficient_a(aph2-1, 0, i2+2, l2+i2);
                            switch (d1){
                            case -1: sum += tmp * calc_coefficient_c(aph2-3, 0, d2a-i2+1, l2+d2a);
                                break;
                            case  0: sum += tmp * (((valuetype) calc_delta(d2a, i2)) - calc_coefficient_a(aph2-1, 0, d2a-i2+2, l2+d2a));
                                break;
                            case  1: sum += tmp * 4 * calc_coefficient_g(aph2-1, 0, d2a-i2+3, l2+d2a);
                            }
                        }
                    C2[4] *= sum; // d1 = -2: 2, d2 = -2: 2
                }
                // special case for (3, 2), or multiindex (0, 1, 1)
                if (l1 == 0 && k1 == 2 && l2 >= 1 && k2 == 0)
                    switch (l2){
                    case 1: C2[4] = ((valuetype) -16)/5;  break;
                    case 2: C2[4] = ((valuetype) 8)/5;    break;
                    case 3: C2[4] = ((valuetype) -16)/35; break;
                    case 4: C2[4] = ((valuetype) 2)/35;   break;
                    }
                for (int d3 = -2; d3 <= 2; ++d3){
                    int d = d1 + d2a; // d1 is exactly actual d1
                    int d3a = d3 - d; // actual d3
                    int k3 = l3 + d3a;
                    // continue if multiindex has negative component or has sum larger than M
                    Multiindex<3> index_col = Unitary_Multiindex[0] * k1 + Unitary_Multiindex[1] * k2 + Unitary_Multiindex[2] * k3;
                    if (!(Zero_Multiindex <= index_col) || index_col.sum() >= M) continue;
                    int aph3 = 2 * l1 + 2 * l2;
                    // corresponds to (3, 2), or multiindex (0, 1, 1)
                    for (int i = 0; i < 6; ++i) // initialize
                        C3[i] = (valuetype) 0;
                    C3[4] = pow(2, aph3+1) / (aph3+2*l3+1.0);
                    valuetype sum = 0;
                    if (d == -2)
                        for (int i3 = ((d3a-2 > 0) ? d3a-2 : 0);
                             i3 <= ((d3a < 2) ? d3a : 2); ++i3){
                                sum += calc_coefficient_c(aph3-2, 0, i3+1, l3+i3) * calc_coefficient_c(aph3-4, 0, d3a-i3+1, l3+d3a);
                        }
                    if (d == 2)
                        for (int i3 = ((d3a > -2) ? d3a : -2);
                             i3 <= ((d3a+2 < 0) ? d3a+2 : 0); ++i3){
                                sum += 16 * calc_coefficient_g(aph3, 0, i3+3, l3+i3) * calc_coefficient_g(aph3+2, 0, d3a-i3+3, l3+d3a);
                        }
                    if (-1 <= d && d <= 1)
                        for (int i3 = ((d3a-1 + d > -1) ? d3a-1 + d : -1);
                             i3 <= ((d3a+1 + d < 1) ? d3a+1 + d : 1); ++i3){
                            valuetype tmp = calc_delta(i3, 0) - calc_coefficient_a(aph3, 0, i3+2, l3+i3);
                            switch (d){
                            case -1: sum += tmp * calc_coefficient_c(aph3-2, 0, d3a-i3+1, l3+d3a);
                                break;
                            case  0: sum += tmp * (((valuetype) calc_delta(d3a, i3)) - calc_coefficient_a(aph3, 0, d3a-i3+2, l3+d3a));
                                break;
                            case  1: sum += tmp * 4 * calc_coefficient_g(aph3, 0, d3a-i3+3, l3+d3a);
                            }
                        }
                    C3[4] *= sum; // d1 = -2: 2, d2 = -2: 2, d3 = -2: 2
                    if (abs(d1) <= 1){
                        // corresponds to (3, 1), or multiindex (1, 0, 1)
                        C3[2] = C3[4]; // d1 = -1: 1, d2 = -2: 2, d3 = -2: 2
                        // corresponds to (3, 3), or multiindex (0, 0, 2)
                        C3[5] = C3[4]; // d1 = -1: 1, d2 = -2: 2, d3 = -2: 2
                    }
                    if (l3 == 0){
                        if (abs(d1) <= 1 && abs(d2) <= 1 && 0 <= d3a && d3a <= 3){
                            // corresponds to (2, 1), or multiindex (1, 1, 0)
                            C3[1] = (valuetype) pow(2, aph3+d+3);
                            valuetype numerator = 1;
                            if (d3a != 0){
                                for (int i = 0; i < d3a-1; ++i)
                                    numerator *= (valuetype) d-1 + i;
                                numerator *= (valuetype) aph3 + 2*d + 1 + d3a;
                            }
                            valuetype denominator = (valuetype) 1;
                            for (int i = 0; i < d3a+1; ++i)
                                denominator *= (valuetype) aph3+d+3 + i;
                            C3[1] *= numerator / denominator; // d1 = -1: 1, d2 = -1: 1, d3a = 0: 3
                            if (d1 == 0)
                                // corresponds to (1, 1), or multiindex (2, 0, 0)
                                C3[0] = C3[1]; // d1 = 0, d2 = -1: 1, d3a = 0: 3
                            // corresponds to (2, 2), or multiindex (0, 2, 0)
                            C3[3] = C3[1]; // d1 = -1: 1, d2 = -1: 1, d3a = 0: 3
                        }
                    }
                    else{
                        if (abs(d1) <= 1 && abs(d2) <= 1){
                            // corresponds to (2, 1), or multiindex (1, 1, 0)
                            C3[1] = ((valuetype) pow(2, aph3+1) * (aph3+l3+1)) / (l3 * (aph3+2*l3+1));
                            sum = (valuetype) 0;
                            for (int i3 = ((d3a-1 + d > -1) ? d3a-1 + d : -1);
                                 i3 <= ((d3a+1 + d < 1) ? d3a+1 + d : 1); ++i3){
                                valuetype tmp = ((valuetype) calc_delta(i3, 0)) + calc_coefficient_a(aph3+1, -1, i3+2, l3+i3);
                                switch (d){
                                case -1: sum += tmp * calc_coefficient_c(aph3-1, -1, d3a-i3+1, l3+d3a);
                                    break;
                                case  0: sum += tmp * (((valuetype) calc_delta(d3a, i3)) - calc_coefficient_a(aph3+1, -1, d3a-i3+2, l3+d3a));
                                    break;
                                case  1: sum += tmp * 4 * calc_coefficient_g(aph3+1, -1, d3a-i3+3, l3+d3a);
                                }
                            }
                            C3[1] *= sum; // d1 = -1: 1, d2 = -1: 1, d3 = -2: 2
                            if (d1 == 0)
                                // corresponds to (1, 1), or multiindex (2, 0, 0)
                                C3[0] = C3[1]; // d1 = 0, d2 = -1: 1, d3 = -2: 2
                            // corresponds to (2, 2), or multiindex (0, 2, 0)
                            C3[3] = C3[1]; // d1 = -1: 1, d2 = -1: 1, d3 = -2: 2
                        }
                    }
                    int col = correspondence.index2number(index_col);
                    for (int i = 0; i < 6; ++i){
                        valuetype ttmp = C1[i] * C2[i] * C3[i];
                        if (fabs(ttmp) < tol_zero) continue;
                        ttmp *= pow((valuetype) 0.5, 2*(l1+k1)+l2+k2+6);
                        stiff_matrix_element[i].add(row, col, ttmp);
                    }
                }
	    }
	}
    }
    std::cerr << "generate 6 stiff matrix element\n";
    
    // calculate coefficient and multiindex for the assignment of discretized matrix
    n_index_variation.resize(6); // number of the coefficient for the derivative of generalized jacobi polynomial
    index_variation.resize(6); // variation of multiindex corresponding to the coefficients
    // assign n_index_variation, intialize index_variation
    n_index_variation[0] = 1; // $\partial x1$,               or multiindex (2, 0, 0)
    n_index_variation[1] = 2; // $\partial x2 - \partial x1$, or multiindex (1, 1, 0)
    n_index_variation[2] = 4; // $\partial x1 - \partial x3$, or multiindex (1, 0, 1)
    n_index_variation[3] = 2; // $\partial x2$,               or multiindex (0, 2, 0)
    n_index_variation[4] = 2; // $\partial x3 - \partial x2$, or multiindex (0, 1, 1)
    n_index_variation[5] = 4; // $\partial x3$,               or multiindex (0, 0, 2)
    for (int ind = 0; ind < 6; ++ind)
        index_variation[ind].resize(n_index_variation[ind]);
    // assign index_variation respectively
    // $\partial x1$,               or multiindex (2, 0, 0)
    index_variation[0][0] = Unitary_Multiindex[0];
    // $\partial x2 - \partial x1$, or multiindex (1, 1, 0), p = 0, 1 in turns
    index_variation[1][0] = Unitary_Multiindex[1];
    index_variation[1][1] = Unitary_Multiindex[0];
    // $\partial x1 - \partial x3$, or multiindex (1, 0, 1), (p, q) = (0, 0), (1, 0), (0, 1), (1, 1) in turns
    index_variation[2][0] = Unitary_Multiindex[2];
    index_variation[2][1] = Unitary_Multiindex[0] - Unitary_Multiindex[1] + Unitary_Multiindex[2];
    index_variation[2][2] = Unitary_Multiindex[1];
    index_variation[2][3] = Unitary_Multiindex[0];
    // $\partial x2$,               or multiindex (0, 2, 0), p = 0, 1 in turns
    for (int i = 0; i < n_index_variation[3]; ++i)
        index_variation[3][i] = index_variation[1][i];
    // $\partial x3 - \partial x2$, or multiindex (0, 1, 1), q = 0, 1 in turns, correspond to (0, q) for (p, q)
    index_variation[4][0] = Unitary_Multiindex[2];
    index_variation[4][1] = Unitary_Multiindex[1];
    // $\partial x3$,               or multiindex (0, 0, 2), (p, q) = (0, 0), (1, 0), (0, 1), (1, 1) in turns
    for (int i = 0; i < n_index_variation[5]; ++i)
        index_variation[5][i] = index_variation[2][i];
    // calculate the corresponding coefficient
    coefficient_derivative.resize(6);
    for (int ind = 0; ind < 6; ++ind){
        coefficient_derivative[ind].resize(n_index_variation[ind]);
        for (int indt = 0; indt < n_index_variation[ind]; ++indt){
            coefficient_derivative[ind][indt].resize(n_index);
            for (int i = 0; i < n_index; ++i)
                if (index_variation[ind][indt] <= correspondence.number2index(i))
                    coefficient_derivative[ind][indt][i] = calc_coefficient_D(ind, indt, correspondence.number2index(i));
                else coefficient_derivative[ind][indt][i] = 0;
        }
    }
    std::cerr << "calculate coefficients for stiff matrix\n";
    
    // assign stiff matrix actual with n_transform_local, transform_local, weight_transform_local
    std::vector<unsigned int> n_nnz_per_row_actual(n_index, 5 * 5 * 5 * 16 * 16);
    //     the largest one, h^{(3,2)} * max(n_transform_local)^2 * max(n_index_variation)^2
    sp_stiff_matrix_element_actual.resize(6);
    stiff_matrix_element_actual.resize(6);
    for (int ind_var = 0; ind_var < 6; ++ind_var){
	// std::cout << "ind_var = " << ind_var << '\n';
        sp_stiff_matrix_element_actual[ind_var].reinit(n_index, n_index, n_nnz_per_row_actual);
        // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = stiff_matrix_element[ind_var].begin(0);
        // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = stiff_matrix_element[ind_var].end(stiff_matrix_element[ind_var].m()-1);
	const SparsityPattern &sp_pattern = stiff_matrix_element[ind_var].get_sparsity_pattern();
	const std::size_t *row_start = sp_pattern.get_rowstart_indices();
	const unsigned int *col_nums = sp_pattern.get_column_numbers();
	for (unsigned int row = 0; row < stiff_matrix_element[ind_var].m(); ++row)
        // for (; spm_ite != spm_end; ++spm_ite){
        //     int row = spm_ite->row();
	    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
            // int col = spm_ite->column();
		const unsigned int &col = col_nums[pos_whole];
	    // std::cout << "row = " << row << ", col = " << col << '\n';
		for (int ind_row_v = 0; ind_row_v < n_index_variation[ind_var]; ++ind_row_v) // v for variation
		    for (int ind_col_v = 0; ind_col_v < n_index_variation[ind_var]; ++ind_col_v){
			Multiindex<3> index_row = correspondence.number2index(row) + index_variation[ind_var][ind_row_v];
			Multiindex<3> index_col = correspondence.number2index(col) + index_variation[ind_var][ind_col_v];
			if (!(Zero_Multiindex <= index_row) || !(Zero_Multiindex <= index_col)) continue;
			int row_actual = correspondence.index2number(index_row);
			int col_actual = correspondence.index2number(index_col);
			for (int ind_row_tl = 0; ind_row_tl < n_transform_local[row_actual]; ++ind_row_tl) // tr for transform local
			    for (int ind_col_tl = 0; ind_col_tl < n_transform_local[col_actual]; ++ind_col_tl){
				sp_stiff_matrix_element_actual[ind_var].add(transform_local[row_actual][ind_row_tl],
									    transform_local[col_actual][ind_col_tl]);
			    }
		    }
	    }
        sp_stiff_matrix_element_actual[ind_var].compress();
	// std::cout << "build sparsity pattern for " << ind_var << "-th actual stiff matrix element\n";

	stiff_matrix_element_actual[ind_var].reinit(sp_stiff_matrix_element_actual[ind_var]);
        // spm_ite = stiff_matrix_element[ind_var].begin(0);
        // for (; spm_ite != spm_end; ++spm_ite){
        //     int row = spm_ite->row();
        //     int col = spm_ite->column();
        //     valuetype val = spm_ite->value();
	for (unsigned int row = 0; row < stiff_matrix_element[ind_var].m(); ++row)
        // for (; spm_ite != spm_end; ++spm_ite){
        //     int row = spm_ite->row();
	    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
            // int col = spm_ite->column();
		const unsigned int &col = col_nums[pos_whole];
		valuetype val = stiff_matrix_element[ind_var].global_entry(pos_whole);
		for (int ind_row_v = 0; ind_row_v < n_index_variation[ind_var]; ++ind_row_v) // v for variation
		    for (int ind_col_v = 0; ind_col_v < n_index_variation[ind_var]; ++ind_col_v){
			Multiindex<3> index_row = correspondence.number2index(row) + index_variation[ind_var][ind_row_v];
			Multiindex<3> index_col = correspondence.number2index(col) + index_variation[ind_var][ind_col_v];
			if (!(Zero_Multiindex <= index_row) || !(Zero_Multiindex <= index_col)) continue;
			int row_actual = correspondence.index2number(index_row);
			int col_actual = correspondence.index2number(index_col);
			for (int ind_row_tl = 0; ind_row_tl < n_transform_local[row_actual]; ++ind_row_tl) // tr for transform local
			    for (int ind_col_tl = 0; ind_col_tl < n_transform_local[col_actual]; ++ind_col_tl)
				stiff_matrix_element_actual[ind_var].add(transform_local[row_actual][ind_row_tl],
									 transform_local[col_actual][ind_col_tl],
									 val
									 * coefficient_derivative[ind_var][ind_row_v][row_actual]
									 * coefficient_derivative[ind_var][ind_col_v][col_actual]
									 * weight_transform_local[row_actual][ind_row_tl]
									 * weight_transform_local[col_actual][ind_col_tl]);
		    
		    }
	    }
    }
    std::cerr << "given transform_local, generate 6 actual stiff matrix element\n";
}

TEMPLATE_TSEM
void THIS_TSEM::build_stiff_matrix(FEMSpace<double, 3> &fem_space)
{ // given fem_space, construct stiff matrix

    // construct discretized matrix and right-hand-side
    n_nonzero_per_row_stiff.resize(n_dof_total, 0);
    // std::cout << n_nonzero_per_row_stiff.size() << '\n';
    // std::cout << "assign n_nonzero_per_row_stiff\n";
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele)
        for (int ind_var = 0; ind_var < 6; ++ind_var){
	    // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = stiff_matrix_element_actual[ind_var].begin(0);
	    // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = stiff_matrix_element_actual[ind_var].end(stiff_matrix_element_actual[ind_var].m()-1);
	    const SparsityPattern &sp_pattern = stiff_matrix_element_actual[ind_var].get_sparsity_pattern();
	    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
	    const unsigned int *col_nums = sp_pattern.get_column_numbers();
	    for (unsigned int row = 0; row < stiff_matrix_element_actual[ind_var].m(); ++row)
		for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		    const unsigned int &col = col_nums[pos_whole];
	    // for (; spm_ite != spm_end; ++spm_ite){
	    // 	int row = spm_ite->row();
	    // 	int col = spm_ite->column();
		    for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][row]; ++ind_nnz)
			n_nonzero_per_row_stiff[transform_ind_global2local[ind_ele][row][ind_nnz]] += transform_n_global2local[ind_ele][col];
		}   
        }
    // std::cout << "count nonzero elements for stiff matrix\n";
    // std::cout << n_nonzero_per_row_stiff.size() << '\n';

    sp_stiff_matrix.reinit(n_dof_total, n_dof_total, n_nonzero_per_row_stiff);
    // std::cout << "reinit sparsity pattern for stiff matrix\n";
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele)
        for (int ind_var = 0; ind_var < 6; ++ind_var){
            // int m_matrix = stiff_matrix_element_actual[ind_var].m();
            // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = stiff_matrix_element_actual[ind_var].begin(0);
            // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = stiff_matrix_element_actual[ind_var].end(stiff_matrix_element_actual[ind_var].m()-1);
	    const SparsityPattern &sp_pattern = stiff_matrix_element_actual[ind_var].get_sparsity_pattern();
	    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
	    const unsigned int *col_nums = sp_pattern.get_column_numbers();
	    for (unsigned int row = 0; row < stiff_matrix_element_actual[ind_var].m(); ++row)
		for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		    const unsigned int &col = col_nums[pos_whole];
            // for (; spm_ite != spm_end; ++spm_ite){
            //     int row = spm_ite->row();
            //     int col = spm_ite->column();
		    if (col >= n_index) std::cout << "find illegal col = " << col << '\n';
		    // std::cout << "row = " << row << ", col = " << col << '\n';
		    for (int ind_nnz_row = 0; ind_nnz_row < transform_n_global2local[ind_ele][row]; ++ind_nnz_row)
			for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][col]; ++ind_nnz_col)
			    sp_stiff_matrix.add(transform_ind_global2local[ind_ele][row][ind_nnz_row],
						transform_ind_global2local[ind_ele][col][ind_nnz_col]);
		}
        }
    sp_stiff_matrix.compress();
    // std::cout << "build sparsity pattern for stiff matrix\n";
    
    stiff_matrix.reinit(sp_stiff_matrix);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	// std::cout << "ind_ele = " << ind_ele << '\n';
        const std::vector<int>& element_dof = fem_space.element(ind_ele).dof();
        valuetype coef = ((valuetype) 2) / (3 * val_volume[ind_ele]);
        // traverse each stiff_matrix_element, assign their contribution to stiff matrix
        for (int ind1 = 0; ind1 <= 2; ++ind1) // traverse index of face
            for (int ind2 = ind1+1; ind2 <= 3; ++ind2){
                int ind_p1, ind_p2; // endpoints of the common line between face $F_{ind1}$ and $F_{ind2}$
                if (ind1 == 0)
                    ind_p1 = ind2 % 2 + 1; // ind2 = 1: 2, 3; 2: 1, 3; 3: 1, 2
                else
                    ind_p1 = 0;
                ind_p2 = 6 - ind1 - ind2 - ind_p1;
                AFEPack::Point<3> p1 = fem_space.dofInfo(element_dof[ind1]).interp_point,   p2 = fem_space.dofInfo(element_dof[ind2]).interp_point;
                AFEPack::Point<3> ps = fem_space.dofInfo(element_dof[ind_p1]).interp_point, pe = fem_space.dofInfo(element_dof[ind_p2]).interp_point;
                valuetype coef_dihedral_angle = 0.25 * calc_inner_product(calc_cross_product(pe - ps, p1 - ps), calc_cross_product(pe - ps, p2 - ps));
                if (fabs(coef_dihedral_angle) < tol_zero) continue;
                int ind = ((ind1 == 0) ? correspondence.index2number(Unitary_Multiindex[ind2-1] * 2) - 4 // corresponds to $x_{ind2}^2$, 4=correspondence.n_begin(2)
                           : correspondence.index2number(Unitary_Multiindex[ind1-1] + Unitary_Multiindex[ind2-1]) - 4); // corresponds to $x_{ind2} - x_{ind1}$
                // add contribution of stiff_matrix_element
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = stiff_matrix_element_actual[ind].begin(0);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = stiff_matrix_element_actual[ind].end(stiff_matrix_element_actual[ind].m()-1);
                // for (; spm_ite != spm_end; ++spm_ite){
                //     int row = spm_ite->row();
                //     int col = spm_ite->column();
                //     valuetype val = spm_ite->value();
		const SparsityPattern &sp_pattern = stiff_matrix_element_actual[ind].get_sparsity_pattern();
		const std::size_t *row_start = sp_pattern.get_rowstart_indices();
		const unsigned int *col_nums = sp_pattern.get_column_numbers();
		for (unsigned int row = 0; row < stiff_matrix_element_actual[ind].m(); ++row)
		    for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
			const unsigned int &col = col_nums[pos_whole];
			valuetype val = stiff_matrix_element_actual[ind].global_entry(pos_whole);
			for (int ind_nnz_row = 0; ind_nnz_row < transform_n_global2local[ind_ele][row]; ++ind_nnz_row)
			    for (int ind_nnz_col = 0; ind_nnz_col < transform_n_global2local[ind_ele][col]; ++ind_nnz_col)
				stiff_matrix.add(transform_ind_global2local[ind_ele][row][ind_nnz_row],
						 transform_ind_global2local[ind_ele][col][ind_nnz_col],
						 val * coef * coef_dihedral_angle
						 * transform_val_global2local[ind_ele][row][ind_nnz_row] * transform_val_global2local[ind_ele][col][ind_nnz_col]);
		    }
            }
    }
    std::cerr << "construct stiff matrix\n";
}
TEMPLATE_TSEM
void THIS_TSEM::calc_rhs(Vector<valuetype> &rhs, std::vector<std::vector<valuetype> > &value_f)
{
    rhs.reinit(n_dof_total);
    for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof) rhs(ind_dof) = 0;
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // valuetype volume = fem_space.element(ind_ele).templateElement().volume();
        // AFEPack::Point<3> point_tmp;
        // for (int ind = 0; ind < 3; ++ind) point_tmp[ind] = 1.0/3;
        // valuetype jacobian = fem_space.element(ind_ele).local_to_global_jacobian(point_tmp);
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            valuetype cnt = 0;
            for (int p = 0; p < n_q_point[2]; ++p)
                cnt += Weight[2][p] * value_f[ind_ele][p] * basis_value_actual[p][ind_index];
            cnt *= val_volume[ind_ele];
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                rhs(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    += cnt * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        }
    }
    std::cerr << "calculate rhs\n";
}

TEMPLATE_TSEM
void THIS_TSEM::calc_rhs(Vector<valuetype> &rhs, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*))
{
    rhs.reinit(n_dof_total);
    for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof) rhs(ind_dof) = 0;
    std::vector<double> value_f(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(ind_ele).local_to_global(QPoint);
	for (int p = 0; p < n_q_point[2]; ++p)
	    value_f[p] = func(q_point[p]);
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            valuetype cnt = 0;
            for (int p = 0; p < n_q_point[2]; ++p)
                cnt += Weight[2][p] * value_f[p] * basis_value_actual[p][ind_index];
            cnt *= val_volume[ind_ele];
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                rhs(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    += cnt * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        }
    }
    std::cerr << "calculate rhs\n";
}

TEMPLATE_TSEM
void THIS_TSEM::calc_rhs(Vector<valuetype> &rhs, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*, double), double t)
{
    rhs.reinit(n_dof_total);
    for (int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof) rhs(ind_dof) = 0;
    std::vector<double> value_f(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	std::vector<AFEPack::Point<3> > q_point = fem_space.element(ind_ele).local_to_global(QPoint);
	for (int p = 0; p < n_q_point[2]; ++p)
	    value_f[p] = func(q_point[p], t);
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            valuetype cnt = 0;
            for (int p = 0; p < n_q_point[2]; ++p)
                cnt += Weight[2][p] * value_f[p] * basis_value_actual[p][ind_index];
            cnt *= val_volume[ind_ele];
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                rhs(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    += cnt * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        }
    }
    std::cerr << "calculate rhs\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_zero_boundary_condition(SparseMatrix<valuetype> &sp_matrix, unsigned int flag_boundary)
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    
    // impose boundary condition
    for (int ind = 0; ind <= 2; ++ind)
        for (int ind_geometry = 0; ind_geometry < n_geometry[ind]; ++ind_geometry){
            // if (mesh.boundaryMark(ind, ind_geometry) != 1) continue;
            if (flag_bm[ind][ind_geometry] != flag_boundary) continue;
            int location = ind_geometry + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
            unsigned int row_begin = location_actualdof[location_geometry[location]];
            for (int i = 0; i < n_dof_geometry[ind]; ++i){
                unsigned int row = row_begin + i;
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = sp_matrix.begin(row);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = sp_matrix.end(row);
                // for (; spm_ite != spm_end; ++spm_ite){
                //     int col = spm_ite->column();
		for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		    const unsigned int &col = col_nums[pos_whole];
                    if (col != row){
                        sp_matrix.set(row, col, 0);
                        sp_matrix.set(col, row, 0);
                    }
                }
            }
        }
    std::cerr << "impose boundary condition to matrix\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_zero_boundary_condition_complex(SparseMatrix<valuetype> &sp_matrix, unsigned int flag_boundary)
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    
    // impose boundary condition
    for (int ind = 0; ind <= 2; ++ind)
        for (int ind_geometry = 0; ind_geometry < n_geometry[ind]; ++ind_geometry){
            // if (mesh.boundaryMark(ind, ind_geometry) != 1) continue;
            if (flag_bm[ind][ind_geometry] != flag_boundary) continue;
            int location = ind_geometry + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
            unsigned int row_begin = location_actualdof[location_geometry[location]];
            for (int i = 0; i < n_dof_geometry[ind]; ++i){
                unsigned int row = row_begin + i;
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = sp_matrix.begin(row);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = sp_matrix.end(row);
                // for (; spm_ite != spm_end; ++spm_ite){
                //     int col = spm_ite->column();
		for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		    const unsigned int &col = col_nums[pos_whole];
                    if (col != row){
                        sp_matrix.set(row, col, 0);
                        sp_matrix.set(col, row, 0);
                    }
                }

		row += n_dof_total;
		for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		    const unsigned int &col = col_nums[pos_whole];
                    if (col != row){
                        sp_matrix.set(row, col, 0);
                        sp_matrix.set(col, row, 0);
                    }
                }
            }
        }
    std::cerr << "impose boundary condition to matrix\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_zero_boundary_condition(Vector<valuetype> &rhs, unsigned int flag_boundary)
{
    // impose zero boundary condition
    for (int ind = 0; ind <= 2; ++ind)
	for (int ind_geometry = 0; ind_geometry < n_geometry[ind]; ++ind_geometry){
	    if (flag_bm[ind][ind_geometry] != flag_boundary) continue;
	    int location = ind_geometry + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
	    int row_start = location_actualdof[location_geometry[location]];
	    for (int i = 0; i < n_dof_geometry[ind]; ++i)
		rhs(row_start + i) = 0;
	}
    std::cerr << "impose boundary condition to rhs\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_zero_boundary_condition_complex(Vector<valuetype> &rhs, unsigned int flag_boundary)
{
    // impose zero boundary condition
    for (int ind = 0; ind <= 2; ++ind)
	for (int ind_geometry = 0; ind_geometry < n_geometry[ind]; ++ind_geometry){
	    if (flag_bm[ind][ind_geometry] != flag_boundary) continue;
	    int location = ind_geometry + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
	    int row_start = location_actualdof[location_geometry[location]];
	    for (int i = 0; i < n_dof_geometry[ind]; ++i)
		rhs(row_start + i) = 0;
	    row_start += n_dof_total;
	    for (int i = 0; i < n_dof_geometry[ind]; ++i)
		rhs(row_start + i) = 0;
	}
    std::cerr << "impose boundary condition to rhs\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_boundary_condition_rowOnly(SparseMatrix<valuetype> &sp_matrix, unsigned int flag_boundary)
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    
    // impose boundary condition
    for (int ind = 0; ind <= 2; ++ind)
        for (int ind_geometry = 0; ind_geometry < n_geometry[ind]; ++ind_geometry){
            if (flag_bm[ind][ind_geometry] != flag_boundary) continue;
            int location = ind_geometry + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
            unsigned int row_begin = location_actualdof[location_geometry[location]];
            for (int i = 0; i < n_dof_geometry[ind]; ++i){
                unsigned int row = row_begin + i;
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = sp_matrix.begin(row);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = sp_matrix.end(row);
                // for (; spm_ite != spm_end; ++spm_ite){
                //     int col = spm_ite->column();
		for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		    const unsigned int &col = col_nums[pos_whole];
                    if (col != row)
                        sp_matrix.set(row, col, 0);
                }
            }
        }
    std::cerr << "impose boundary condition to matrix\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_boundary_condition_rowOnly(SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &rhs, std::vector<std::vector<std::vector<valuetype> > > &value_bnd, unsigned int flag_boundary)
{
    
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    
    // impose boundary condition
    for (int ind = 0; ind <= 2; ++ind)
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
            // if (mesh.boundaryMark(ind, ind_geo) != 1) continue;
            if (flag_bm[ind][ind_geo] != flag_boundary) continue;
            int location = ind_geo + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
            unsigned int row_begin = location_actualdof[location_geometry[location]];
            std::vector<double> val(n_dof_geometry[ind], 0);
            // evaluate val in the same way to interpolation
            if (ind == 0) // 0 dimensional geometry, only 1 point
                val[0] = value_bnd[ind][ind_geo][0];
            if (ind == 1){ // 1 dimensional geometry
                int ind_point_s = number_node[0][ind_geo][0], ind_point_e = number_node[0][ind_geo][1];
                int location_s = location_actualdof[location_geometry[ind_point_s]], location_e = location_actualdof[location_geometry[ind_point_e]];
                double c_s = rhs(location_s) / sp_matrix.diag_element(location_s), c_e = rhs(location_e) / sp_matrix.diag_element(location_e);
                for (int l = 2; l <= M; ++l){ // dof locate on this edge
                    double count = 0;
                    for (int p = 0; p < n_q_point[0]; ++p)
                        count += Weight[0][p]
			    * (value_bnd[ind][ind_geo][p] - c_s*QPoint_Barycentric[0][p][1] - c_e*QPoint_Barycentric[0][p][0]) * basis_value_interp[0][p][l-2];
                    val[l-2] = count * -l * (2*l - 1) / (2 * (l-1));
                }
            }
            if (ind == 2){ // 2 dimensional geometry
                // traverse 2 dimensional geometry
                for (int l1 = 2; l1 <= M; ++l1)
                    for (int l2 = 1; l2 <= M-l1; ++l2){
                        int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                        double count = 0;
                        for (int p = 0; p < n_q_point[1]; ++p){
                            double xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1], rp = QPoint_Barycentric[1][p][2];
                            double count_p = value_bnd[ind][ind_geo][p]; // the contribution from function u
                            for (int ind = 0; ind < 3; ++ind){ // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
                                int& pos_vertex_global = location_actualdof[location_geometry[number_node[1][ind_geo][ind]]];
                                // double val_interp_vertex = rhs(pos_vertex_global) / stiff_matrix.diag_element(pos_vertex_global);
                                double val_interp_vertex = rhs(pos_vertex_global) / sp_matrix.diag_element(pos_vertex_global);
                                count_p -= val_interp_vertex * QPoint_Barycentric[1][p][(ind+2)%3];
                            }
                            int location_dof_edge[3]; // start location of dof on each edge
                            for (int ind_edge = 0; ind_edge < 3; ++ind_edge)
                                location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_geo][ind_edge] + n_geometry[0]]];
                            std::vector<std::vector<double> > val_edgedof(3);
                            for (int ind_e = 0; ind_e < 3; ++ind_e){
                                val_edgedof[ind_e].resize(n_dof_edge);
                                for (int l = 2; l <= M; ++l){
                                    val_edgedof[ind_e][l-2] = rhs(location_dof_edge[ind_e] + l-2) / sp_matrix.diag_element(location_dof_edge[ind_e] + l-2);
                                    if (!flag_sameorder_edgeonface[ind_geo][ind_e] && l % 2 == 1)
                                        val_edgedof[ind_e][l-2] *= -1;
                                }
                            }
                            for (int l = 2; l <= M; ++l){ // substract the contribution from edge
                                count_p += 2 * (xp*yp * val_edgedof[0][l-2] * basis_value_addition[p][l-2][1] +
                                                yp*rp * val_edgedof[1][l-2] * basis_value_addition[p][l-2][1] +
                                                xp*rp * val_edgedof[2][l-2] * basis_value_addition[p][l-2][0]);
                            }
                            count += Weight[1][p] * count_p * basis_value_interp[1][p][ind_index];
                        }
                        val[ind_index] = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
                    }
            }
	    if (ind == 2)
		 for (int i = 0; i < n_dof_geometry[ind]; ++i){
		     unsigned int row = row_begin + conversion[0][i];
		     rhs(row) = sp_matrix.diag_element(row) * val[i];
		 }
	    else
		for (int i = 0; i < n_dof_geometry[ind]; ++i){
		    unsigned int row = row_begin + i;
		    rhs(row) = sp_matrix.diag_element(row) * val[i];
		}
        }
    std::cerr << "impose boundary condition\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_boundary_condition(SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &rhs, std::vector<std::vector<std::vector<valuetype> > > &value_bnd,
					  bool flag_modify_matrix, unsigned int flag_boundary)
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    
    // impose boundary condition
    for (int ind = 0; ind <= 2; ++ind)
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
            // if (mesh.boundaryMark(ind, ind_geo) != 1) continue;
            // if (flag_bm[ind][ind_geo] != flag_boundary) continue;
	    if (flag_boundary == -1 && flag_bm[ind][ind_geo] == 0) continue;
	    if (flag_boundary != -1 && flag_bm[ind][ind_geo] != flag_boundary) continue;
            int location = ind_geo + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
            unsigned int row_begin = location_actualdof[location_geometry[location]];
            std::vector<double> val(n_dof_geometry[ind], 0);
            // evaluate val in the same way to interpolation
            if (ind == 0) // 0 dimensional geometry, only 1 point
                val[0] = value_bnd[ind][ind_geo][0];
            if (ind == 1){ // 1 dimensional geometry
                int ind_point_s = number_node[0][ind_geo][0], ind_point_e = number_node[0][ind_geo][1];
                int location_s = location_actualdof[location_geometry[ind_point_s]], location_e = location_actualdof[location_geometry[ind_point_e]];
                double c_s = rhs(location_s) / sp_matrix.diag_element(location_s), c_e = rhs(location_e) / sp_matrix.diag_element(location_e);
                for (int l = 2; l <= M; ++l){ // dof locate on this edge
                    double count = 0;
                    for (int p = 0; p < n_q_point[0]; ++p)
                        count += Weight[0][p]
			    * (value_bnd[ind][ind_geo][p] - c_s*QPoint_Barycentric[0][p][1] - c_e*QPoint_Barycentric[0][p][0]) * basis_value_interp[0][p][l-2];
                    val[l-2] = count * -l * (2*l - 1) / (2 * (l-1));
                }
            }
            if (ind == 2){ // 2 dimensional geometry
                // traverse 2 dimensional geometry
                for (int l1 = 2; l1 <= M; ++l1)
                    for (int l2 = 1; l2 <= M-l1; ++l2){
                        int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                        double count = 0;
                        for (int p = 0; p < n_q_point[1]; ++p){
                            double xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1], rp = QPoint_Barycentric[1][p][2];
                            double count_p = value_bnd[ind][ind_geo][p]; // the contribution from function u
                            for (int ind = 0; ind < 3; ++ind){ // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
                                int& pos_vertex_global = location_actualdof[location_geometry[number_node[1][ind_geo][ind]]];
                                // double val_interp_vertex = rhs(pos_vertex_global) / stiff_matrix.diag_element(pos_vertex_global);
                                double val_interp_vertex = rhs(pos_vertex_global) / sp_matrix.diag_element(pos_vertex_global);
                                count_p -= val_interp_vertex * QPoint_Barycentric[1][p][(ind+2)%3];
                            }
                            int location_dof_edge[3]; // start location of dof on each edge
                            for (int ind_edge = 0; ind_edge < 3; ++ind_edge)
                                location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_geo][ind_edge] + n_geometry[0]]];
                            std::vector<std::vector<double> > val_edgedof(3);
                            for (int ind_e = 0; ind_e < 3; ++ind_e){
                                val_edgedof[ind_e].resize(n_dof_edge);
                                for (int l = 2; l <= M; ++l){
                                    val_edgedof[ind_e][l-2] = rhs(location_dof_edge[ind_e] + l-2) / sp_matrix.diag_element(location_dof_edge[ind_e] + l-2);
                                    if (!flag_sameorder_edgeonface[ind_geo][ind_e] && l % 2 == 1)
                                        val_edgedof[ind_e][l-2] *= -1;
                                }
                            }
                            for (int l = 2; l <= M; ++l){ // substract the contribution from edge
                                count_p += 2 * (xp*yp * val_edgedof[0][l-2] * basis_value_addition[p][l-2][1] +
                                                yp*rp * val_edgedof[1][l-2] * basis_value_addition[p][l-2][1] +
                                                xp*rp * val_edgedof[2][l-2] * basis_value_addition[p][l-2][0]);
                            }
                            count += Weight[1][p] * count_p * basis_value_interp[1][p][ind_index];
                        }
                        val[ind_index] = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
                    }
            }
            for (int i = 0; i < n_dof_geometry[ind]; ++i){
                unsigned int row = row_begin + (ind == 2 ? conversion[0][i] : i);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = sp_matrix.begin(row);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = sp_matrix.end(row);
                // for (; spm_ite != spm_end; ++spm_ite){
                //     int col = spm_ite->column();
		for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
		    const unsigned int &col = col_nums[pos_whole];
                    if (col == row)
                        // rhs(row) = spm_ite->value() * val[i];
			// rhs(row) = sp_matrix.global_entry(pos_whole) * val[i];
			rhs(row) = sp_matrix.diag_element(row) * val[i];
                    else{
			// if (!flag_bm_dof[col])
			    rhs(col) -= sp_matrix.el(col, row) * val[i];
			if (flag_modify_matrix){
			    sp_matrix.set(row, col, 0);
			    sp_matrix.set(col, row, 0);
			}
                    }
                }
            }
        }
    std::cerr << "impose boundary condition\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_boundary_condition(SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &rhs,
					  RegularMesh<3> &mesh, valuetype(*func)(double*),
					  bool flag_modify_matrix, unsigned int flag_boundary)
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();

    // impose boundary condition
    for (int ind = 0; ind <= 2; ++ind)
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
            // if (mesh.boundaryMark(ind, ind_geo) != 1) continue;
            // if (flag_bm[ind][ind_geo] != flag_boundary) continue;
	    if (flag_boundary == -1 && flag_bm[ind][ind_geo] == 0) continue;
	    if (flag_boundary != -1 && flag_bm[ind][ind_geo] != flag_boundary) continue;
            int location = ind_geo + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
            unsigned int row_begin = location_actualdof[location_geometry[location]];
            std::vector<valuetype> val(n_dof_geometry[ind], 0);
            // evaluate val in the same way to interpolation
            if (ind == 0) // 0 dimensional geometry, only 1 point
                val[0] = func(mesh.point(ind_geo));
            if (ind == 1){ // 1 dimensional geometry
                int ind_point_s = number_node[0][ind_geo][0], ind_point_e = number_node[0][ind_geo][1];
    	    	std::vector<valuetype> value_bnd(n_q_point[0]);
    	    	for (int p = 0; p < n_q_point[0]; ++p){
    	    	    AFEPack::Point<3> p_tmp = mesh.point(ind_point_s), p_ttmp = mesh.point(ind_point_e);
    	    	    p_tmp  *= QPoint_Barycentric[0][p][1];
    	    	    p_ttmp *= QPoint_Barycentric[0][p][0];
    	    	    p_tmp  += p_ttmp;
    	    	    value_bnd[p] = func(p_tmp);
    	    	}
                int location_s = location_actualdof[location_geometry[ind_point_s]], location_e = location_actualdof[location_geometry[ind_point_e]];
                valuetype c_s = rhs(location_s) / sp_matrix.diag_element(location_s), c_e = rhs(location_e) / sp_matrix.diag_element(location_e);
                for (int l = 2; l <= M; ++l){ // dof locate on this edge
                    valuetype count = 0;
                    for (int p = 0; p < n_q_point[0]; ++p)
                        count += Weight[0][p]
    	    		    * (value_bnd[p] - c_s*QPoint_Barycentric[0][p][1] - c_e*QPoint_Barycentric[0][p][0]) * basis_value_interp[0][p][l-2];
                    val[l-2] = count * -l * (2*l - 1) / (2 * (l-1));
                }
            }
            if (ind == 2){ // 2 dimensional geometry
    	    	std::vector<valuetype> value_bnd(n_q_point[1]);
    	    	for (int p = 0; p < n_q_point[1]; ++p){
    	    	    AFEPack::Point<3> p_tmp;
    	    	    for (int indt = 0; indt < 3; ++indt) p_tmp[indt] = 0;
    	    	    for (int ind_p = 0; ind_p <= 2; ++ind_p){ // vertex[0]=(0,0)<->[2]=1-x-y; vertex[1]=(1,0)<->[0]=x; vertex[2]=(0,1)<->[1]=y
    	    		int ind_vertex = number_node[1][ind_geo][ind_p];
    	    		AFEPack::Point<3> p_ttmp = mesh.point(ind_vertex);
    	    		p_ttmp *= QPoint_Barycentric[1][p][(ind_p+2)%3];
    	    		p_tmp  += p_ttmp;
    	    	    }
    	    	    value_bnd[p] = func(p_tmp);
    	    	}
                // traverse 2 dimensional geometry
                for (int l1 = 2; l1 <= M; ++l1)
                    for (int l2 = 1; l2 <= M-l1; ++l2){
                        int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                        valuetype count = 0;
                        for (int p = 0; p < n_q_point[1]; ++p){
                            valuetype xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1], rp = QPoint_Barycentric[1][p][2];
                            valuetype count_p = value_bnd[p]; // the contribution from function u
                            for (int ind = 0; ind < 3; ++ind){ // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
                                int& pos_vertex_global = location_actualdof[location_geometry[number_node[1][ind_geo][ind]]];
                                // valuetype val_interp_vertex = rhs(pos_vertex_global) / stiff_matrix.diag_element(pos_vertex_global);
                                valuetype val_interp_vertex = rhs(pos_vertex_global) / sp_matrix.diag_element(pos_vertex_global);
                                count_p -= val_interp_vertex * QPoint_Barycentric[1][p][(ind+2)%3];
                            }
                            int location_dof_edge[3]; // start location of dof on each edge
                            for (int ind_edge = 0; ind_edge < 3; ++ind_edge)
                                location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_geo][ind_edge] + n_geometry[0]]];
                            std::vector<std::vector<valuetype> > val_edgedof(3);
                            for (int ind_e = 0; ind_e < 3; ++ind_e){
                                val_edgedof[ind_e].resize(n_dof_edge);
                                for (int l = 2; l <= M; ++l){
                                    val_edgedof[ind_e][l-2] = rhs(location_dof_edge[ind_e] + l-2) / sp_matrix.diag_element(location_dof_edge[ind_e] + l-2);
                                    if (!flag_sameorder_edgeonface[ind_geo][ind_e] && l % 2 == 1)
                                        val_edgedof[ind_e][l-2] *= -1;
                                }
                            }
                            for (int l = 2; l <= M; ++l){ // substract the contribution from edge
                                count_p += 2 * (xp*yp * val_edgedof[0][l-2] * basis_value_addition[p][l-2][1] +
                                                yp*rp * val_edgedof[1][l-2] * basis_value_addition[p][l-2][1] +
                                                xp*rp * val_edgedof[2][l-2] * basis_value_addition[p][l-2][0]);
                            }
                            count += Weight[1][p] * count_p * basis_value_interp[1][p][ind_index];
                        }
                        val[ind_index] = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
                    }
	    }
            for (int i = 0; i < n_dof_geometry[ind]; ++i){
                unsigned int row = row_begin + (ind == 2 ? conversion[0][i] : i);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = sp_matrix.begin(row);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = sp_matrix.end(row);
                // for (; spm_ite != spm_end; ++spm_ite){
                //     int col = spm_ite->column();
    	    	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
    	    	    const unsigned int &col = col_nums[pos_whole];
                    if (col == row)
                        // rhs(row) = spm_ite->value() * val[i];
    	    		// rhs(row) = sp_matrix.global_entry(pos_whole) * val[i];
    	    		rhs(row) = sp_matrix.diag_element(row) * val[i];
                    else{
    	    		// if (!flag_bm_dof[col])
    	    		    rhs(col) -= sp_matrix.el(col, row) * val[i];
    	    		if (flag_modify_matrix){
    	    		    sp_matrix.set(row, col, 0);
    	    		    sp_matrix.set(col, row, 0);
    	    		}
                    }
                }
            }
        }
    std::cerr << "impose boundary condition\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_boundary_condition(SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &rhs,
					  Mesh<3>& mesh, valuetype(*func)(double*),
					  bool flag_modify_matrix, unsigned int flag_boundary)
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();

    // impose boundary condition
    for (int ind = 0; ind <= 2; ++ind)
        for (int ind_geo = 0; ind_geo < n_geometry[ind]; ++ind_geo){
	    if (flag_boundary == -1 && flag_bm[ind][ind_geo] == 0) continue;
	    if (flag_boundary != -1 && flag_bm[ind][ind_geo] != flag_boundary) continue;
            int location = ind_geo + ((ind > 0) ? n_geometry[0] : 0) + ((ind > 1) ? n_geometry[1] : 0);
            unsigned int row_begin = location_actualdof[location_geometry[location]];
            std::vector<valuetype> val(n_dof_geometry[ind], 0);
            // evaluate val in the same way to interpolation
            if (ind == 0) // 0 dimensional geometry, only 1 point
                val[0] = func(mesh.point(ind_geo));
            if (ind == 1){ // 1 dimensional geometry
                int ind_point_s = number_node[0][ind_geo][0], ind_point_e = number_node[0][ind_geo][1];
    	    	std::vector<valuetype> value_bnd(n_q_point[0]);
    	    	for (int p = 0; p < n_q_point[0]; ++p){
    	    	    AFEPack::Point<3> p_tmp = mesh.point(ind_point_s), p_ttmp = mesh.point(ind_point_e);
    	    	    p_tmp  *= QPoint_Barycentric[0][p][1];
    	    	    p_ttmp *= QPoint_Barycentric[0][p][0];
    	    	    p_tmp  += p_ttmp;
    	    	    value_bnd[p] = func(p_tmp);
    	    	}
                int location_s = location_actualdof[location_geometry[ind_point_s]], location_e = location_actualdof[location_geometry[ind_point_e]];
                valuetype c_s = rhs(location_s) / sp_matrix.diag_element(location_s), c_e = rhs(location_e) / sp_matrix.diag_element(location_e);
                for (int l = 2; l <= M; ++l){ // dof locate on this edge
                    valuetype count = 0;
                    for (int p = 0; p < n_q_point[0]; ++p)
                        count += Weight[0][p]
    	    		    * (value_bnd[p] - c_s*QPoint_Barycentric[0][p][1] - c_e*QPoint_Barycentric[0][p][0]) * basis_value_interp[0][p][l-2];
                    val[l-2] = count * -l * (2*l - 1) / (2 * (l-1));
                }
            }
            if (ind == 2){ // 2 dimensional geometry
    	    	std::vector<valuetype> value_bnd(n_q_point[1]);
    	    	for (int p = 0; p < n_q_point[1]; ++p){
    	    	    AFEPack::Point<3> p_tmp;
    	    	    for (int indt = 0; indt < 3; ++indt) p_tmp[indt] = 0;
    	    	    for (int ind_p = 0; ind_p <= 2; ++ind_p){ // vertex[0]=(0,0)<->[2]=1-x-y; vertex[1]=(1,0)<->[0]=x; vertex[2]=(0,1)<->[1]=y
    	    		int ind_vertex = number_node[1][ind_geo][ind_p];
    	    		AFEPack::Point<3> p_ttmp = mesh.point(ind_vertex);
    	    		p_ttmp *= QPoint_Barycentric[1][p][(ind_p+2)%3];
    	    		p_tmp  += p_ttmp;
    	    	    }
    	    	    value_bnd[p] = func(p_tmp);
    	    	}
                // traverse 2 dimensional geometry
                for (int l1 = 2; l1 <= M; ++l1)
                    for (int l2 = 1; l2 <= M-l1; ++l2){
                        int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                        valuetype count = 0;
                        for (int p = 0; p < n_q_point[1]; ++p){
                            valuetype xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1], rp = QPoint_Barycentric[1][p][2];
                            valuetype count_p = value_bnd[p]; // the contribution from function u
                            for (int ind = 0; ind < 3; ++ind){ // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
                                int& pos_vertex_global = location_actualdof[location_geometry[number_node[1][ind_geo][ind]]];
                                // valuetype val_interp_vertex = rhs(pos_vertex_global) / stiff_matrix.diag_element(pos_vertex_global);
                                valuetype val_interp_vertex = rhs(pos_vertex_global) / sp_matrix.diag_element(pos_vertex_global);
                                count_p -= val_interp_vertex * QPoint_Barycentric[1][p][(ind+2)%3];
                            }
                            int location_dof_edge[3]; // start location of dof on each edge
                            for (int ind_edge = 0; ind_edge < 3; ++ind_edge)
                                location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_geo][ind_edge] + n_geometry[0]]];
                            std::vector<std::vector<valuetype> > val_edgedof(3);
                            for (int ind_e = 0; ind_e < 3; ++ind_e){
                                val_edgedof[ind_e].resize(n_dof_edge);
                                for (int l = 2; l <= M; ++l){
                                    val_edgedof[ind_e][l-2] = rhs(location_dof_edge[ind_e] + l-2) / sp_matrix.diag_element(location_dof_edge[ind_e] + l-2);
                                    if (!flag_sameorder_edgeonface[ind_geo][ind_e] && l % 2 == 1)
                                        val_edgedof[ind_e][l-2] *= -1;
                                }
                            }
                            for (int l = 2; l <= M; ++l){ // substract the contribution from edge
                                count_p += 2 * (xp*yp * val_edgedof[0][l-2] * basis_value_addition[p][l-2][1] +
                                                yp*rp * val_edgedof[1][l-2] * basis_value_addition[p][l-2][1] +
                                                xp*rp * val_edgedof[2][l-2] * basis_value_addition[p][l-2][0]);
                            }
                            count += Weight[1][p] * count_p * basis_value_interp[1][p][ind_index];
                        }
                        val[ind_index] = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
                    }
	    }
            for (int i = 0; i < n_dof_geometry[ind]; ++i){
                unsigned int row = row_begin + (ind == 2 ? conversion[0][i] : i);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_ite = sp_matrix.begin(row);
                // SparseMatrix<VALUETYPE_ITERATOR_TSEM>::iterator spm_end = sp_matrix.end(row);
                // for (; spm_ite != spm_end; ++spm_ite){
                //     int col = spm_ite->column();
    	    	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
    	    	    const unsigned int &col = col_nums[pos_whole];
                    if (col == row)
                        // rhs(row) = spm_ite->value() * val[i];
    	    		// rhs(row) = sp_matrix.global_entry(pos_whole) * val[i];
    	    		rhs(row) = sp_matrix.diag_element(row) * val[i];
                    else{
    	    		// if (!flag_bm_dof[col])
    	    		    rhs(col) -= sp_matrix.el(col, row) * val[i];
    	    		if (flag_modify_matrix){
    	    		    sp_matrix.set(row, col, 0);
    	    		    sp_matrix.set(col, row, 0);
    	    		}
                    }
                }
            }
        }
    std::cerr << "impose boundary condition\n";
}

TEMPLATE_TSEM
void THIS_TSEM::impose_boundary_condition(Vector<double>& rhs, RegularMesh<3>& mesh, valuetype(*func)(double*), unsigned int flag_boundary)
{// impose nuemann boundary condition
    for (unsigned int idx_face = 0; idx_face < n_geometry[2]; ++idx_face){
	if (flag_bm[2][idx_face] != flag_boundary) continue;
	// prepare contribution vector
	std::vector<valuetype> cnt_v(3, 0);
	std::vector<std::vector<valuetype> > cnt_e(3, std::vector<valuetype> (n_dof_geometry[1], 0));
	std::vector<valuetype> cnt_f(n_dof_geometry[2], 0);
	std::vector<int>& node = number_node[1][idx_face], edge = number_edge[idx_face];
	// attain vertex coordinates
	AFEPack::Point<3>& p0 = mesh.point(node[0]), p1 = mesh.point(node[1]), p2 = mesh.point(node[2]);
	// calculate volume, weight * volume
	valuetype volume = calc_area_triangle(p0, p1, p2);
	// traverse quadrature points, count contributions
	for (unsigned int idx_qp = 0; idx_qp < n_q_point[1]; ++idx_qp){
	    // attain physical coordinates of quadrature point
	    std::vector<valuetype>& qb = QPoint_Barycentric[1][idx_qp];
	    AFEPack::Point<3> p_qp(p0[0]*qb[2] + p1[0]*qb[0] + p2[0]*qb[1],
				   p0[1]*qb[2] + p1[1]*qb[0] + p2[1]*qb[1],
				   p0[2]*qb[2] + p1[2]*qb[0] + p2[2]*qb[1]);
	    // calculate value of neumann boundary
	    valuetype val_g = func(p_qp);
	    valuetype wgv = Weight[1][idx_qp] * val_g * volume;
	    // count contribution to each dof
	    for (unsigned int idx_v = 0; idx_v < 3; ++idx_v) cnt_v[idx_v] += wgv * qb[(idx_v+2)%3];
	    valuetype& xp = qb[0], yp = qb[1], rp = qb[2];
	    for (unsigned int idx_dof = 0; idx_dof < n_dof_geometry[1]; ++idx_dof){
		cnt_e[0][idx_dof] += wgv * -2*xp*yp*basis_value_addition[idx_qp][idx_dof][1];
		cnt_e[1][idx_dof] += wgv * -2*yp*rp*basis_value_addition[idx_qp][idx_dof][1];
		cnt_e[2][idx_dof] += wgv * -2*xp*rp*basis_value_addition[idx_qp][idx_dof][0];
	    }
	    for (unsigned int idx_dof = 0; idx_dof < n_dof_geometry[2]; ++idx_dof)
		cnt_f[idx_dof] += wgv * basis_value[1][idx_qp][idx_dof];
	}
	// add contribution to each dof
	std::vector<int> idx_dof_v(3), idx_dof_e(3);
	for (unsigned int idx = 0; idx < 3; ++idx){
	    idx_dof_v[idx] = location_actualdof[location_geometry[node[idx]]];
	    idx_dof_e[idx] = location_actualdof[location_geometry[edge[idx]+n_geometry[0]]];
	}
	int idx_dof_f = location_actualdof[location_geometry[n_geometry[0]+n_geometry[1]+idx_face]];
	for (unsigned int idx_v = 0; idx_v < 3; ++idx_v) rhs(idx_dof_v[idx_v]) += cnt_v[idx_v];
	for (unsigned int idx_e = 0; idx_e < 3; ++idx_e){
	    if (flag_sameorder_edgeonface[idx_face][idx_e]) continue;
	    for (unsigned int idx_dof = 0; idx_dof < n_dof_geometry[1]; ++idx_dof) if (idx_dof % 2 == 1) cnt_e[idx_e][idx_dof] *= -1;
	}
	for (unsigned int idx_e = 0; idx_e < 3; ++idx_e)
	    for (unsigned int idx_dof = 0; idx_dof < n_dof_geometry[1]; ++idx_dof)
		rhs(idx_dof_e[idx_e] + idx_dof) += cnt_e[idx_e][idx_dof];
	for (unsigned int idx_dof = 0; idx_dof < n_dof_geometry[2]; ++idx_dof)
	    rhs(idx_dof_f + idx_dof) += cnt_f[idx_dof];
    }
}

TEMPLATE_TSEM
bool THIS_TSEM::is_on_boundary(AFEPack::Point<3> &p, valuetype bnd_left, valuetype bnd_right)
{// judge whether point p is on boundary of [Bnd_Left, Bnd_Right]^3
    bool flag = false;
    for (int ind = 0; ind < 3; ++ind){
        bool flag_onleftbnd  = fabs(p[ind] - bnd_left)  < tol_zero;
        bool flag_onrightbnd = fabs(p[ind] - bnd_right) < tol_zero;
        if (flag_onleftbnd || flag_onrightbnd) flag = true;
    }
    // if (fabs(p[0] - 1.0/3) < 1.0e-8 && fabs(p[1] - 1.0/3) < 1.0e-8 && fabs(p[2] - 1.0/3) < 1.0e-8)
    //     flag = true; // for tetrahedron
    return flag;
}

TEMPLATE_TSEM
bool THIS_TSEM::is_on_boundary(AFEPack::Point<3> &p, valuetype radius)
{// judge whether point p is on boundary of ball centering at (0,0,0) with radius
    // std::cout << "p = (" << p[0] << ", " << p[1] << ", " << p[2] << "), dis = " << fabs(radius-sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])) << '\n';
    // return fabs(sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]) - radius) < radius * 0.1;
    return fabs(sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]) - radius) < tol_zero;
}


#undef TEMPLATE_TSEM
#undef THIS_TSEM
