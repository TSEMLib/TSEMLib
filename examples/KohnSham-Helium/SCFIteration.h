#ifndef __SCFITERATION_H__
#define __SCFITERATION_H__

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>

#include "AFEPack/Geometry.h"
#include "AFEPack/FEMSpace.h"

#include "lac/sparsity_pattern.h"
#include "lac/sparse_matrix.h"
#include "lac/vector.h"

#include "TetrahedralSEM.h"
#include "TSEMSolver.h"
#include "TSEMVisualize.h"

#define DIM 3
#define PI (4*atan(1.0))

double wave_init(double*);

class SCFIteration{

private:
    double Tol_Zero = 1.0e-12;
    
    // filename of info
    std::string info_evolve;
    // mesh
    HGeometryTree<DIM> h_tree;
    IrregularMesh<DIM> *irregular_mesh;
    // fem space
    TemplateGeometry<DIM> template_geometry;
    CoordTransform<DIM, DIM> coord_transform;
    TemplateDOF<DIM> template_dof;
    BasisFunctionAdmin<double, DIM, DIM> basis_function;
    std::vector<TemplateElement<double, DIM, DIM> > template_element;
    int n_element;
    FEMSpace<double, DIM> fem_space;
    // tsem
    int order_tsem;
    TSEM<double> tsem;
    int n_dof_total, n_qp;
    // evolution info
    double coef_mix;
    double tolerance_scf, tolerance_lobpcg, tolerance_solver;
    int max_step_scf, max_step_lobpcg, max_step_solver;
    // nucleus info
    unsigned int n_nu = 1, charge_tot = 2, n_orbital = 1;
    std::vector<int> charge_nu = {2};
    std::vector<std::vector<double> > pos_nu = {{0,0,0}};
    std::vector<double> n_occupation = {2};
    // io info
    bool flag_read, flag_write, flag_write_geometry, flag_output;
    std::string filename_read, filename_write, filename_output_energy, filename_output_orthonormal;

    // psi and density, energy
    std::vector<Vector<double> > psi;
    std::vector<std::vector<double> > density;

    // value of external potential on element
    std::vector<std::vector<double> > val_V_Ext;
    // Vector for Hartree potential
    Vector<double> v_Har;
    // coefficient for the calculation of multipole expansion
    double coef_me[10]; // me for multipole expansion, :={<1>, <x1>, <x2>, <x3>, <x1^2>, <x1x2>, <x1x3>, <x2^2>, <x2x3>, <x3^2>}
    std::vector<std::vector<std::vector<double> > > coef_vH_me;

    double ene_kin = 0;

private:
    // functions for calculating Hartree potential
    void init_coef_vH_me();
    void calc_coef_me();
    double calc_multipole_expansion(AFEPack::Point<DIM>&);
    void calc_V_Har_bnd(std::vector<std::vector<std::vector<double> > >&);
    void calc_V_Har();

    // functions for calculating potential and energy
    double calc_V_Ext(AFEPack::Point<DIM>&, std::vector<int>&, std::vector<std::vector<double> >&);
    double calc_V_X(double);
    double calc_E_X(double);
    double calc_V_C(double);
    double calc_E_C(double);
    double calc_V_XC(double);
    double calc_E_XC(double);
    double calc_energy(std::vector<Vector<double> >&, bool);


    // function for SCF iteration
    void normalize();
    void assemble_matrix(std::vector<Vector<double> >&, std::vector<unsigned int>&, SparsityPattern&, SparseMatrix<double>&);
    void update_density();

public:
    SCFIteration(const std::string& filename_evolve);
    ~SCFIteration();

    void init();
    void run();
    void solve();
    void outputSolution();
};

class Visualize : public TSEMVisualize<double>
{
public:
    Visualize();
    ~Visualize();

    virtual void get_indicator();
};

#endif
