#include "../include/TetrahedralSEM.h"
// This file contains functions for interpolation
// Functions:
// calc_interpolation()


#define TEMPLATE_TSEM template<typename valuetype>
#define THIS_TSEM TSEM<valuetype>


TEMPLATE_TSEM
void THIS_TSEM::calc_interpolation(Vector<valuetype> &interp, std::vector<std::vector<std::vector<valuetype> > > &val_interp)
{
    interp.reinit(n_dof_total);
    
    // traverse 0 dimensional geometry, correspond to vertex in fem dof
    for (int ind_point = 0; ind_point < n_geometry[0]; ++ind_point)
        interp(location_actualdof[location_geometry[ind_point]]) = val_interp[0][ind_point][0];

    // traverse 1 dimensional geometry
    for (int ind_edge = 0; ind_edge < n_geometry[1]; ++ind_edge){
	int ind_point_s = number_node[0][ind_edge][0], ind_point_e = number_node[0][ind_edge][1];
        double c_s = interp(location_actualdof[location_geometry[ind_point_s]]), c_e = interp(location_actualdof[location_geometry[ind_point_e]]);
        for (int l = 2; l <= M; ++l){ // dof locate on this edge
            int location_dof = location_actualdof[location_geometry[ind_edge+n_geometry[0]]] + l-2; // position of this dof
            valuetype count = 0;
            for (int p = 0; p < n_q_point[0]; ++p)
                count += Weight[0][p] * (val_interp[1][ind_edge][p] - c_s*QPoint_Barycentric[0][p][1] - c_e*QPoint_Barycentric[0][p][0]) * basis_value_interp[0][p][l-2];
            interp(location_dof) = count * -l * (2*l - 1) / (2 * (l-1));
        }
    }
    
    // traverse 2 dimensional geometry
    for (int ind_face = 0; ind_face < n_geometry[2]; ++ind_face){
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2){
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                int location_dof = location_actualdof[location_geometry[ind_face+n_geometry[0]+n_geometry[1]]] + ind_index;
                valuetype count = 0;
                for (int p = 0; p < n_q_point[1]; ++p){
                    valuetype xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1], rp = QPoint_Barycentric[1][p][2];
                    valuetype count_p = val_interp[2][ind_face][p]; // the contribution from function u
                    for (int ind = 0; ind < 3; ++ind) // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
                        count_p -= interp(location_actualdof[location_geometry[number_node[1][ind_face][ind]]]) * QPoint_Barycentric[1][p][(ind+2)%3];
                    int location_dof_edge[3]; // start location of dof on each edge
                    for (int ind_edge = 0; ind_edge < 3; ++ind_edge)
                        location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_face][ind_edge] + n_geometry[0]]];
                    std::vector<std::vector<valuetype> > val_edgedof(3);
                    for (int ind_e = 0; ind_e < 3; ++ind_e){
                        val_edgedof[ind_e].resize(n_dof_edge);
                        for (int l = 2; l <= M; ++l){
                            val_edgedof[ind_e][l-2] = interp(location_dof_edge[ind_e] + l-2);
                            if (!flag_sameorder_edgeonface[ind_face][ind_e] && l % 2 == 1)
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
                interp(location_dof) = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
            }
    }
    
    // traverse fem element, assign value for interior dof
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
        // get local coefficients for vertex, edge and face
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            Multiindex<3> index_now = correspondence.number2index(ind_index);
            if (index_now.index[0] >= 2 && index_now.index[1] >= 1 && index_now.index[2] >= 1)
                continue;
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += interp(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        }
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2)
                for (int l3 = 1; l3 <= M-l1-l2; ++l3){
                    int ind_index = correspondence.index2number(Unitary_Multiindex[0] * (l1-2) + Unitary_Multiindex[1] * (l2-1) + Unitary_Multiindex[2] * (l3-1));
                    int location_dof = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]] + ind_index;
                    valuetype count = 0;
                    for (int p = 0; p < n_q_point[2]; ++p){
                        valuetype count_p = val_interp[3][ind_ele][p];
                        for (int j = 0; j < n_index; ++j){ // traverse multiindex for vertex, edge and face
                            Multiindex<3> index_tmp = correspondence.number2index(j);
                            if (index_tmp.index[0] >= 2 && index_tmp.index[1] >= 1 && index_tmp.index[2] >= 1) // the interior model function are orthogonal
                                continue;
                            count_p -= val_coef[j] * basis_value_actual[p][j];
                        }
                        count += Weight[2][p] * count_p *  basis_value_interp[2][p][ind_index];
                    }
                    // for interior dof, global coefficient is exactly local one
                    interp(location_dof) = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) * (2*l1+2*l2+2*l3-1) / (pow(2, 2*l1+l2-5) * (l1-1) * 6);
                }
    }
    std::cerr << "interpolation exact solution\n";
}

TEMPLATE_TSEM
void THIS_TSEM::calc_interpolation(Vector<valuetype> &interp, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*))
{
    interp.reinit(n_dof_total);

    // traverse 0 dimensional geometry, correspond to vertex in fem dof
    for (int ind_point = 0; ind_point < n_geometry[0]; ++ind_point)
        interp(location_actualdof[location_geometry[ind_point]]) = func(mesh.point(ind_point));

    // traverse 1 dimensional geometry
    std::vector<double> val_interp1(n_q_point[0]);
    for (int ind_edge = 0; ind_edge < n_geometry[1]; ++ind_edge){
	int ind_point_s = number_node[0][ind_edge][0], ind_point_e = number_node[0][ind_edge][1];
        double c_s = interp(location_actualdof[location_geometry[ind_point_s]]), c_e = interp(location_actualdof[location_geometry[ind_point_e]]);
	for (int p = 0; p < n_q_point[0]; ++p){
	    AFEPack::Point<3> p_tmp = mesh.point(ind_point_s), p_ttmp = mesh.point(ind_point_e);
	    p_tmp  *= QPoint_Barycentric[0][p][1];
	    p_ttmp *= QPoint_Barycentric[0][p][0];
	    p_tmp  += p_ttmp;
	    val_interp1[p] = func(p_tmp);
	}
        for (int l = 2; l <= M; ++l){ // dof locate on this edge
            int location_dof = location_actualdof[location_geometry[ind_edge+n_geometry[0]]] + l-2; // position of this dof
            valuetype count = 0;
            for (int p = 0; p < n_q_point[0]; ++p)
                count += Weight[0][p] * (val_interp1[p] - c_s*QPoint_Barycentric[0][p][1] - c_e*QPoint_Barycentric[0][p][0]) * basis_value_interp[0][p][l-2];
            interp(location_dof) = count * -l * (2*l - 1) / (2 * (l-1));
        }
    }
    
    // traverse 2 dimensional geometry
    std::vector<double> val_interp2(n_q_point[1]);
    for (int ind_face = 0; ind_face < n_geometry[2]; ++ind_face){
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2){
		for (int p = 0; p < n_q_point[1]; ++p){
    		    AFEPack::Point<3> p_tmp;
    		    for (int ind = 0; ind < 3; ++ind) p_tmp[ind] = 0;
    		    for (int ind_p = 0; ind_p <= 2; ++ind_p){ // vertex[0]=(0,0)<->[2]=1-x-y; vertex[1]=(1,0)<->[0]=x; vertex[2]=(0,1)<->[1]=y
    			int ind_vertex = number_node[1][ind_face][ind_p];
    			AFEPack::Point<3> p_ttmp = mesh.point(ind_vertex);
    			p_ttmp *= QPoint_Barycentric[1][p][(ind_p+2)%3];
    			p_tmp  += p_ttmp;
    		    }
    		    val_interp2[p] = func(p_tmp);
    		}
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                int location_dof = location_actualdof[location_geometry[ind_face+n_geometry[0]+n_geometry[1]]] + ind_index;
                valuetype count = 0;
                for (int p = 0; p < n_q_point[1]; ++p){
                    valuetype xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1], rp = QPoint_Barycentric[1][p][2];
                    valuetype count_p = val_interp2[p]; // the contribution from function u
                    for (int ind = 0; ind < 3; ++ind) // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
                        count_p -= interp(location_actualdof[location_geometry[number_node[1][ind_face][ind]]]) * QPoint_Barycentric[1][p][(ind+2)%3];
                    int location_dof_edge[3]; // start location of dof on each edge
                    for (int ind_edge = 0; ind_edge < 3; ++ind_edge)
                        location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_face][ind_edge] + n_geometry[0]]];
                    std::vector<std::vector<valuetype> > val_edgedof(3);
                    for (int ind_e = 0; ind_e < 3; ++ind_e){
                        val_edgedof[ind_e].resize(n_dof_edge);
                        for (int l = 2; l <= M; ++l){
                            val_edgedof[ind_e][l-2] = interp(location_dof_edge[ind_e] + l-2);
                            if (!flag_sameorder_edgeonface[ind_face][ind_e] && l % 2 == 1)
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
                interp(location_dof) = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
            }
    }
    
    // traverse fem element, assign value for interior dof
    std::vector<double> val_interp3(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    val_interp3[p] = func(fem_space.element(ind_ele).local_to_global(QPoint[p]));
        // get local coefficients for vertex, edge and face
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            Multiindex<3> index_now = correspondence.number2index(ind_index);
            if (index_now.index[0] >= 2 && index_now.index[1] >= 1 && index_now.index[2] >= 1)
                continue;
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += interp(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        }
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2)
                for (int l3 = 1; l3 <= M-l1-l2; ++l3){
                    int ind_index = correspondence.index2number(Unitary_Multiindex[0] * (l1-2) + Unitary_Multiindex[1] * (l2-1) + Unitary_Multiindex[2] * (l3-1));
                    int location_dof = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]] + ind_index;
                    valuetype count = 0;
                    for (int p = 0; p < n_q_point[2]; ++p){
                        valuetype count_p = val_interp3[p];
                        for (int j = 0; j < n_index; ++j){ // traverse multiindex for vertex, edge and face
                            Multiindex<3> index_tmp = correspondence.number2index(j);
                            if (index_tmp.index[0] >= 2 && index_tmp.index[1] >= 1 && index_tmp.index[2] >= 1) // the interior model function are orthogonal
                                continue;
                            count_p -= val_coef[j] * basis_value_actual[p][j];
                        }
                        count += Weight[2][p] * count_p *  basis_value_interp[2][p][ind_index];
                    }
                    // for interior dof, global coefficient is exactly local one
                    interp(location_dof) = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) * (2*l1+2*l2+2*l3-1) / (pow(2, 2*l1+l2-5) * (l1-1) * 6);
                }
    }
    std::cerr << "interpolation exact solution\n";
}

TEMPLATE_TSEM
void THIS_TSEM::calc_interpolation(Vector<valuetype> &interp, Mesh<3> &mesh, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*))
{
    interp.reinit(n_dof_total);

    // traverse 0 dimensional geometry, correspond to vertex in fem dof
    for (int ind_point = 0; ind_point < n_geometry[0]; ++ind_point)
        interp(location_actualdof[location_geometry[ind_point]]) = func(mesh.point(ind_point));

    // traverse 1 dimensional geometry
    std::vector<double> val_interp1(n_q_point[0]);
    for (int ind_edge = 0; ind_edge < n_geometry[1]; ++ind_edge){
	int ind_point_s = number_node[0][ind_edge][0], ind_point_e = number_node[0][ind_edge][1];
        double c_s = interp(location_actualdof[location_geometry[ind_point_s]]), c_e = interp(location_actualdof[location_geometry[ind_point_e]]);
	for (int p = 0; p < n_q_point[0]; ++p){
	    AFEPack::Point<3> p_tmp = mesh.point(ind_point_s), p_ttmp = mesh.point(ind_point_e);
	    p_tmp  *= QPoint_Barycentric[0][p][1];
	    p_ttmp *= QPoint_Barycentric[0][p][0];
	    p_tmp  += p_ttmp;
	    val_interp1[p] = func(p_tmp);
	}
        for (int l = 2; l <= M; ++l){ // dof locate on this edge
            int location_dof = location_actualdof[location_geometry[ind_edge+n_geometry[0]]] + l-2; // position of this dof
            valuetype count = 0;
            for (int p = 0; p < n_q_point[0]; ++p)
                count += Weight[0][p] * (val_interp1[p] - c_s*QPoint_Barycentric[0][p][1] - c_e*QPoint_Barycentric[0][p][0]) * basis_value_interp[0][p][l-2];
            interp(location_dof) = count * -l * (2*l - 1) / (2 * (l-1));
        }
    }
    
    // traverse 2 dimensional geometry
    std::vector<double> val_interp2(n_q_point[1]);
    for (int ind_face = 0; ind_face < n_geometry[2]; ++ind_face){
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2){
		for (int p = 0; p < n_q_point[1]; ++p){
    		    AFEPack::Point<3> p_tmp;
    		    for (int ind = 0; ind < 3; ++ind) p_tmp[ind] = 0;
    		    for (int ind_p = 0; ind_p <= 2; ++ind_p){ // vertex[0]=(0,0)<->[2]=1-x-y; vertex[1]=(1,0)<->[0]=x; vertex[2]=(0,1)<->[1]=y
    			int ind_vertex = number_node[1][ind_face][ind_p];
    			AFEPack::Point<3> p_ttmp = mesh.point(ind_vertex);
    			p_ttmp *= QPoint_Barycentric[1][p][(ind_p+2)%3];
    			p_tmp  += p_ttmp;
    		    }
    		    val_interp2[p] = func(p_tmp);
    		}
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                int location_dof = location_actualdof[location_geometry[ind_face+n_geometry[0]+n_geometry[1]]] + ind_index;
                valuetype count = 0;
                for (int p = 0; p < n_q_point[1]; ++p){
                    valuetype xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1], rp = QPoint_Barycentric[1][p][2];
                    valuetype count_p = val_interp2[p]; // the contribution from function u
                    for (int ind = 0; ind < 3; ++ind) // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
                        count_p -= interp(location_actualdof[location_geometry[number_node[1][ind_face][ind]]]) * QPoint_Barycentric[1][p][(ind+2)%3];
                    int location_dof_edge[3]; // start location of dof on each edge
                    for (int ind_edge = 0; ind_edge < 3; ++ind_edge)
                        location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_face][ind_edge] + n_geometry[0]]];
                    std::vector<std::vector<valuetype> > val_edgedof(3);
                    for (int ind_e = 0; ind_e < 3; ++ind_e){
                        val_edgedof[ind_e].resize(n_dof_edge);
                        for (int l = 2; l <= M; ++l){
                            val_edgedof[ind_e][l-2] = interp(location_dof_edge[ind_e] + l-2);
                            if (!flag_sameorder_edgeonface[ind_face][ind_e] && l % 2 == 1)
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
                interp(location_dof) = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
            }
    }
    
    // traverse fem element, assign value for interior dof
    std::vector<double> val_interp3(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    val_interp3[p] = func(fem_space.element(ind_ele).local_to_global(QPoint[p]));
        // get local coefficients for vertex, edge and face
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            Multiindex<3> index_now = correspondence.number2index(ind_index);
            if (index_now.index[0] >= 2 && index_now.index[1] >= 1 && index_now.index[2] >= 1)
                continue;
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += interp(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        }
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2)
                for (int l3 = 1; l3 <= M-l1-l2; ++l3){
                    int ind_index = correspondence.index2number(Unitary_Multiindex[0] * (l1-2) + Unitary_Multiindex[1] * (l2-1) + Unitary_Multiindex[2] * (l3-1));
                    int location_dof = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]] + ind_index;
                    valuetype count = 0;
                    for (int p = 0; p < n_q_point[2]; ++p){
                        valuetype count_p = val_interp3[p];
                        for (int j = 0; j < n_index; ++j){ // traverse multiindex for vertex, edge and face
                            Multiindex<3> index_tmp = correspondence.number2index(j);
                            if (index_tmp.index[0] >= 2 && index_tmp.index[1] >= 1 && index_tmp.index[2] >= 1) // the interior model function are orthogonal
                                continue;
                            count_p -= val_coef[j] * basis_value_actual[p][j];
                        }
                        count += Weight[2][p] * count_p *  basis_value_interp[2][p][ind_index];
                    }
                    // for interior dof, global coefficient is exactly local one
                    interp(location_dof) = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) * (2*l1+2*l2+2*l3-1) / (pow(2, 2*l1+l2-5) * (l1-1) * 6);
                }
    }
    std::cerr << "interpolation exact solution\n";
}

TEMPLATE_TSEM
void THIS_TSEM::calc_interpolation(Vector<valuetype> &interp, RegularMesh<3> &mesh, FEMSpace<double, 3> &fem_space, valuetype(*func)(double*, double), double t)
{
    interp.reinit(n_dof_total);

    // traverse 0 dimensional geometry, correspond to vertex in fem dof
    for (int ind_point = 0; ind_point < n_geometry[0]; ++ind_point)
        interp(location_actualdof[location_geometry[ind_point]]) = func(mesh.point(ind_point), t);

    // traverse 1 dimensional geometry
    std::vector<double> val_interp1(n_q_point[0]);
    for (int ind_edge = 0; ind_edge < n_geometry[1]; ++ind_edge){
	int ind_point_s = number_node[0][ind_edge][0], ind_point_e = number_node[0][ind_edge][1];
        double c_s = interp(location_actualdof[location_geometry[ind_point_s]]), c_e = interp(location_actualdof[location_geometry[ind_point_e]]);
	for (int p = 0; p < n_q_point[0]; ++p){
	    AFEPack::Point<3> p_tmp = mesh.point(ind_point_s), p_ttmp = mesh.point(ind_point_e);
	    p_tmp  *= QPoint_Barycentric[0][p][1];
	    p_ttmp *= QPoint_Barycentric[0][p][0];
	    p_tmp  += p_ttmp;
	    val_interp1[p] = func(p_tmp, t);
	}
        for (int l = 2; l <= M; ++l){ // dof locate on this edge
            int location_dof = location_actualdof[location_geometry[ind_edge+n_geometry[0]]] + l-2; // position of this dof
            valuetype count = 0;
            for (int p = 0; p < n_q_point[0]; ++p)
                count += Weight[0][p] * (val_interp1[p] - c_s*QPoint_Barycentric[0][p][1] - c_e*QPoint_Barycentric[0][p][0]) * basis_value_interp[0][p][l-2];
            interp(location_dof) = count * -l * (2*l - 1) / (2 * (l-1));
        }
    }
    
    // traverse 2 dimensional geometry
    std::vector<double> val_interp2(n_q_point[1]);
    for (int ind_face = 0; ind_face < n_geometry[2]; ++ind_face){
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2){
		for (int p = 0; p < n_q_point[1]; ++p){
    		    AFEPack::Point<3> p_tmp;
    		    for (int ind = 0; ind < 3; ++ind) p_tmp[ind] = 0;
    		    for (int ind_p = 0; ind_p <= 2; ++ind_p){ // vertex[0]=(0,0)<->[2]=1-x-y; vertex[1]=(1,0)<->[0]=x; vertex[2]=(0,1)<->[1]=y
    			int ind_vertex = number_node[1][ind_face][ind_p];
    			AFEPack::Point<3> p_ttmp = mesh.point(ind_vertex);
    			p_ttmp *= QPoint_Barycentric[1][p][(ind_p+2)%3];
    			p_tmp  += p_ttmp;
    		    }
    		    val_interp2[p] = func(p_tmp, t);
    		}
                int ind_index = (1+l1+l2-3) * (l1+l2-3) / 2 + l2-1; // order of multiindex (l1-2, l2-1)
                int location_dof = location_actualdof[location_geometry[ind_face+n_geometry[0]+n_geometry[1]]] + ind_index;
                valuetype count = 0;
                for (int p = 0; p < n_q_point[1]; ++p){
                    valuetype xp = QPoint_Barycentric[1][p][0], yp = QPoint_Barycentric[1][p][1], rp = QPoint_Barycentric[1][p][2];
                    valuetype count_p = val_interp2[p]; // the contribution from function u
                    for (int ind = 0; ind < 3; ++ind) // substract the contribution from vertex: 0 -> 2, 1 -> 0, 2 -> 1 ((x + 2) % 3)
                        count_p -= interp(location_actualdof[location_geometry[number_node[1][ind_face][ind]]]) * QPoint_Barycentric[1][p][(ind+2)%3];
                    int location_dof_edge[3]; // start location of dof on each edge
                    for (int ind_edge = 0; ind_edge < 3; ++ind_edge)
                        location_dof_edge[ind_edge] = location_actualdof[location_geometry[number_edge[ind_face][ind_edge] + n_geometry[0]]];
                    std::vector<std::vector<valuetype> > val_edgedof(3);
                    for (int ind_e = 0; ind_e < 3; ++ind_e){
                        val_edgedof[ind_e].resize(n_dof_edge);
                        for (int l = 2; l <= M; ++l){
                            val_edgedof[ind_e][l-2] = interp(location_dof_edge[ind_e] + l-2);
                            if (!flag_sameorder_edgeonface[ind_face][ind_e] && l % 2 == 1)
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
                interp(location_dof) = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) / (pow(2, l1) * (l1-1));
            }
    }
    
    // traverse fem element, assign value for interior dof
    std::vector<double> val_interp3(n_q_point[2]);
    for (int ind_ele = 0; ind_ele < n_element; ++ind_ele){
	for (unsigned int p = 0; p < n_q_point[2]; ++p)
	    val_interp3[p] = func(fem_space.element(ind_ele).local_to_global(QPoint[p]), t);
        // get local coefficients for vertex, edge and face
        std::vector<valuetype> val_coef(n_index, 0);
        for (int ind_index = 0; ind_index < n_index; ++ind_index){
            Multiindex<3> index_now = correspondence.number2index(ind_index);
            if (index_now.index[0] >= 2 && index_now.index[1] >= 1 && index_now.index[2] >= 1)
                continue;
            for (int ind_nnz = 0; ind_nnz < transform_n_global2local[ind_ele][ind_index]; ++ind_nnz)
                val_coef[ind_index] += interp(transform_ind_global2local[ind_ele][ind_index][ind_nnz])
                    * transform_val_global2local[ind_ele][ind_index][ind_nnz];
        }
        for (int l1 = 2; l1 <= M; ++l1)
            for (int l2 = 1; l2 <= M-l1; ++l2)
                for (int l3 = 1; l3 <= M-l1-l2; ++l3){
                    int ind_index = correspondence.index2number(Unitary_Multiindex[0] * (l1-2) + Unitary_Multiindex[1] * (l2-1) + Unitary_Multiindex[2] * (l3-1));
                    int location_dof = location_actualdof[location_geometry[ind_ele + n_geometry[0] + n_geometry[1] + n_geometry[2]]] + ind_index;
                    valuetype count = 0;
                    for (int p = 0; p < n_q_point[2]; ++p){
                        valuetype count_p = val_interp3[p];
                        for (int j = 0; j < n_index; ++j){ // traverse multiindex for vertex, edge and face
                            Multiindex<3> index_tmp = correspondence.number2index(j);
                            if (index_tmp.index[0] >= 2 && index_tmp.index[1] >= 1 && index_tmp.index[2] >= 1) // the interior model function are orthogonal
                                continue;
                            count_p -= val_coef[j] * basis_value_actual[p][j];
                        }
                        count += Weight[2][p] * count_p *  basis_value_interp[2][p][ind_index];
                    }
                    // for interior dof, global coefficient is exactly local one
                    interp(location_dof) = count * -l1 * (2*l1-1) * (2*l1+2*l2-1) * (2*l1+2*l2+2*l3-1) / (pow(2, 2*l1+l2-5) * (l1-1) * 6);
                }
    }
    std::cerr << "interpolation exact solution\n";
}


TEMPLATE_TSEM
void THIS_TSEM::calc_interpolation(std::vector<Vector<valuetype> > &psi, std::vector<valuetype> &n_occ,  RegularMesh<3> &mesh, const std::string &filename_interpMesh, std::string &filename_out)
{
    Mesh<3> mesh_out; mesh_out.readData(filename_interpMesh);

    TemplateGeometry<3> template_geometry;
    CoordTransform<3, 3> coord_transform;
    TemplateDOF<3> template_dof;
    BasisFunctionAdmin<double, 3, 3> basis_function;
    std::vector<TemplateElement<double, 3, 3> > template_element;
    template_geometry.readData("tetrahedron.tmp_geo");
    coord_transform.readData("tetrahedron.crd_trs");
    template_dof.reinit(template_geometry);     template_dof.readData("tetrahedron.1.tmp_dof");
    basis_function.reinit(template_dof);        basis_function.readData("tetrahedron.1.bas_fun");
    template_element.resize(1);
    template_element[0].reinit(template_geometry, template_dof, coord_transform, basis_function);

    FEMSpace<double, 3> fem_space; fem_space.reinit(mesh_out, template_element);
    int n_ele = mesh_out.n_geometry(3);
    fem_space.element().resize(n_ele);
    for (unsigned int idx_ele = 0; idx_ele < n_ele; ++idx_ele) fem_space.element(idx_ele).reinit(fem_space, idx_ele, 0);
    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();
    int n_dof = fem_space.n_dof();

    FEMFunction<double, 3> output(fem_space);
    for (unsigned int idx_dof = 0; idx_dof < n_dof; ++idx_dof){
	AFEPack::Point<3>& p = fem_space.dofInfo(idx_dof).interp_point;
	output(idx_dof) = calc_val_point(mesh, p, psi, n_occ);
    }
    output.writeOpenDXData(filename_out);
}


#undef TEMPLATE_TSEM
#undef THIS_TSEM
