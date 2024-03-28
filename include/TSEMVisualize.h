/*
 * visualization of tetrahedral spectral element method
 */
#ifndef __TSEMVISUALIZE_
#define __TSEMVISUALIZE_

#include <string>
#include <vector>
#include <algorithm>

#include "AFEPack/FEMSpace.h"
#include "AFEPack/Geometry.h"
#include "AFEPack/HGeometry.h"

#include "lac/vector.h"

#include "TetrahedralSEM.h"


template<typename valuetype> class TSEMVisualize
{
public:
    TSEM<valuetype>* tsem;
protected:
    // coefficient for AFEPack local refinement
    valuetype convergenceCoefficient = 8, refine_threshold = 1.33333;
    int n_global_refine, n_local_refine;
    valuetype rate_refine;
    valuetype tolerance;
    int n_smooth_indicator;
    // for building interpolation
    HGeometryTree<3>* h_tree;
    Vector<valuetype>* sol;
    std::vector<Vector<valuetype> >* psi;
    std::vector<valuetype>* n_occ;
    unsigned int flag_interp = 0; // 1: single TSEM solution; 2: TSEM wave functions with corresponding occupation numbers
public:
    // for mesh adaption
    IrregularMesh<3>* irregular_mesh_original;
    IrregularMesh<3>* irregular_mesh;
    std::vector<valuetype> val_indicator, val_interp;

public:
    TSEMVisualize();
    ~TSEMVisualize();

public:
    int  n_globalRefine() const { return n_global_refine; };
    int& n_globalRefine(){ return n_global_refine; };
    int  n_localRefine() const { return n_local_refine; };
    int& n_localRefine(){ return n_local_refine; };
    valuetype  rateRefine() const { return rate_refine; };
    valuetype& rateRefine(){return rate_refine; };
    valuetype  Tolerance() const { return tolerance; };
    valuetype& Tolerance(){ return tolerance; };
    int  n_smoothIndicator() const { return n_smooth_indicator; };
    int& n_smoothIndicator(){ return n_smooth_indicator; };

    void read(TSEM<valuetype>* p_tsem, HGeometryTree<3>* hg_tree, Vector<valuetype>* solution);
    void read(TSEM<valuetype>* p_tsem, HGeometryTree<3>* hg_tree, std::vector<Vector<valuetype> >* wave_function, std::vector<valuetype>* n_occupation);
    void build();
    virtual void interp();
    virtual void get_indicator();
    virtual void smooth_indicator();
    void write(std::string filename_output);
    void write(std::string filename_mesh, std::string filename_output);
};



#undef __TSEMVISUALIZE_
#endif
