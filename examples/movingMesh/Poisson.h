#ifndef __POISSON_H_
#define __POISSON_H_

#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <algorithm>

#include "AFEPack/Geometry.h"
#include "AFEPack/FEMSpace.h"
#include "AFEPack/MovingMesh3D.h"

#include "lac/sparsity_pattern.h"
#include "lac/sparse_matrix.h"
#include "lac/vector.h"

#include "TetrahedralSEM.h"
#include "TSEMSolver.h"

#define DIM 3
#define PI (4*atan(1.0))

class Poisson : public MovingMesh3D
{
    
private:
    double Tol_Zero = 1.0e-12;
    double Coef_Scale = 1;
    
    // mesh
    std::string path_mesh;
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
    int n_dof_total;

    // r-adaptive info
    int n_move, n_ite;
    SparseMatrix<double> sp_matrix;
    Vector<double> solution, rhs;
    
    

public:
    Poisson(const std::string&, int, int);
    virtual ~Poisson();

    void init();
    void run();
    void solve();
    void output_info(int);
    void output_err(int, std::string);
    virtual void getMonitor();
    virtual void updateSolution();
    virtual void outputSolution();
};



#endif
