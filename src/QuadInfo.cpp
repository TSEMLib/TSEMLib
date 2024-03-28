#include "../include/TetrahedralSEM.h"
// This file contains functions for the configuration of the quadrature from the input and building the quadrature information for calculations
// Functions:
// read_quad_info()
// build_quad_info()


#define TEMPLATE_TSEM template<typename valuetype>
#define THIS_TSEM TSEM<valuetype>


TEMPLATE_TSEM
void THIS_TSEM::read_quad_info(std::vector<int> &n_qp, std::vector<std::string> &q_fn)
{ // read quadrature information, initialize corresponding variables

    n_q_point.resize(3);
    // std::cout << "n_q_point.size() = " << n_q_point.size() << '\n';
    // record number of quadrature points
    for (int ind = 0; ind < 3; ++ind)
	n_q_point[ind] = n_qp[ind];
    // initialize QPoint, QPoint_Barycentric
    QPoint.resize(n_q_point[2]);
    QPoint_Barycentric.resize(2);
    for (int ind = 0; ind < 2; ++ind){
        QPoint_Barycentric[ind].resize(n_q_point[ind]);
        for (int i = 0; i < n_q_point[ind]; ++i)
            QPoint_Barycentric[ind][i].resize(ind+2);
    }
    Weight.resize(3);
    for (int ind = 0; ind < 3; ++ind)
	Weight[ind].resize(n_q_point[ind]);
    // read quadrature info from q_fn
    for (int ind = 0; ind < 3; ++ind){
	std::ifstream infile;
	infile.open(q_fn[ind]);
	valuetype sum_weight = 0;
        for (int p = 0; p < n_q_point[ind]; ++p){
            for (int indt = 0; indt <= ind; ++indt)
                if (ind == 2)
                    infile >> QPoint[p][indt];
                else
                    infile >> QPoint_Barycentric[ind][p][indt];
            if (ind == 0) // modify coordinate to be barycenter one
                QPoint_Barycentric[ind][p][0] = (QPoint_Barycentric[ind][p][0] + 1) / 2;
            if (ind < 2){ // calculate the barycenter coordinate of the last point
                QPoint_Barycentric[ind][p][ind+1] = 1;
                for (int indt = 0; indt <= ind; ++indt)
                    QPoint_Barycentric[ind][p][ind+1] -= QPoint_Barycentric[ind][p][indt];
            }
            infile >> Weight[ind][p];
	    sum_weight += Weight[ind][p];
        }
        std::cerr << "Read " << ind+1 << "d quadrature info, find " << n_q_point[ind] << " pairs of points and weights, with weight summation " << sum_weight << '\n';
        infile.close();
    }
    // std::cout << "QPoint.size() = " << QPoint.size() << '\n';
    // for (int p = 0; p < QPoint.size(); ++p)
    // 	std::cout << "([" << QPoint[p][0] << ',' << QPoint[p][1] << ',' << QPoint[p][2] << "], " << Weight[2][p] << ") ";
    // std::cout << '\n';
}

TEMPLATE_TSEM
void THIS_TSEM::build_quad_info(int qinfo_order)
{
    n_q_point.resize(3);
    for (unsigned int idx = 0; idx < 3; ++idx) n_q_point[idx] = pow(qinfo_order, idx+1);
    QPoint.resize(n_q_point[2]);
    QPoint_Barycentric.resize(2);
    for (unsigned int idx = 0; idx < 2; ++idx) QPoint_Barycentric[idx].resize(n_q_point[idx], std::vector<valuetype> (idx+2));
    Weight.resize(3);
    for (unsigned int idx = 0; idx < 3; ++idx) Weight[idx].resize(n_q_point[idx]);
    
    std::ostringstream oss;
    oss << path_1DQuadInfo << "qinfo_1d_order" << qinfo_order;
    std::ifstream in(oss.str());
    double sum = 0;
    for (unsigned int idx_qp = 0; idx_qp < n_q_point[0]; ++idx_qp){
	in >> QPoint_Barycentric[0][idx_qp][0] >> Weight[0][idx_qp];
	QPoint_Barycentric[0][idx_qp][0] = (1+QPoint_Barycentric[0][idx_qp][0])*0.5; // [-1, 1] -> [0, 1]
	QPoint_Barycentric[0][idx_qp][1] = 1 - QPoint_Barycentric[0][idx_qp][0];
	sum += Weight[0][idx_qp];
	// std::cout << QPoint_Barycentric[0][idx_qp][0] << '\t' << QPoint_Barycentric[0][idx_qp][1] << '\t' << Weight[0][idx_qp] << '\n';
    }
    in.close();
    std::cerr << "Read 1d qinfo, order = " << qinfo_order << ", find " << n_q_point[0] << " points, with sum of weight " << sum << '\n';

    unsigned int n_qp_1d = n_q_point[0];
    sum = 0;
    for (unsigned int idx1_qp = 0; idx1_qp < n_qp_1d; ++idx1_qp)
	for (unsigned int idx2_qp = 0; idx2_qp < n_qp_1d; ++idx2_qp){
	    unsigned int idx_qp_2d = idx1_qp*n_qp_1d + idx2_qp;
	    std::vector<valuetype> &p1_1d = QPoint_Barycentric[0][idx1_qp];
	    std::vector<valuetype> &p2_1d = QPoint_Barycentric[0][idx2_qp];
	    std::vector<valuetype> &p_2d  = QPoint_Barycentric[1][idx_qp_2d];
	    p_2d[0] = p1_1d[0] * p2_1d[1];
	    p_2d[1] = p2_1d[0];
	    p_2d[2] = 1 - p_2d[0] - p_2d[1];
	    Weight[1][idx_qp_2d] = 2*p2_1d[1] * Weight[0][idx1_qp]*Weight[0][idx2_qp];
	    sum += Weight[1][idx_qp_2d];
	}
    std::cerr << "Build 2d quadrature info from 1d qinfo, generate " << n_q_point[1] << " points, with sum of weight " << sum << '\n';

    sum = 0;
    for (unsigned int idx1_qp = 0; idx1_qp < n_qp_1d; ++idx1_qp)
	for (unsigned int idx2_qp = 0; idx2_qp < n_qp_1d; ++idx2_qp)
	    for (unsigned int idx3_qp = 0; idx3_qp < n_qp_1d; ++idx3_qp){
		unsigned int idx_qp_3d = (idx1_qp*n_qp_1d + idx2_qp)*n_qp_1d + idx3_qp;
		std::vector<valuetype> &p1_1d = QPoint_Barycentric[0][idx1_qp];
		std::vector<valuetype> &p2_1d = QPoint_Barycentric[0][idx2_qp];
		std::vector<valuetype> &p3_1d = QPoint_Barycentric[0][idx3_qp];
		AFEPack::Point<3> &p_3d  = QPoint[idx_qp_3d];
		p_3d[0] = p1_1d[0] * p2_1d[1] * p3_1d[1];
		p_3d[1] = p2_1d[0] * p3_1d[1];
		p_3d[2] = p3_1d[0];
		Weight[2][idx_qp_3d] = 6 * p2_1d[1] * pow(p3_1d[1],2) * Weight[0][idx1_qp]*Weight[0][idx2_qp]*Weight[0][idx3_qp];
		sum += Weight[2][idx_qp_3d];
	    }
    std::cerr << "Build 3d quadrature info from 1d qinfo, generate " << n_q_point[2] << " points, with sum of weight " << sum << '\n';
}

#undef TEMPLATE_TSEM
#undef THIS_TSEM
