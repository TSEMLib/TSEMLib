#include "../include/TSEMVisualize.h"
// This contains functions for visualization of TSEM solution
// Functions:
// read()
// build()
// interp()
// get_indicator()
// smooth_indicator()
// write()


template class TSEMVisualize<double>;
template class TSEMVisualize<float>;

#define TEMPLATE_TSEMVISUALIZE template<typename valuetype>
#define THIS_TSEMVISUALIZE TSEMVisualize<valuetype>


TEMPLATE_TSEMVISUALIZE
THIS_TSEMVISUALIZE::TSEMVisualize() :
    n_global_refine(2), n_local_refine(3), rate_refine(0.5),
    n_smooth_indicator(2), tolerance(1.0e-8)
{}

TEMPLATE_TSEMVISUALIZE
THIS_TSEMVISUALIZE::~TSEMVisualize()
{
    delete irregular_mesh_original;
    delete irregular_mesh;
}

TEMPLATE_TSEMVISUALIZE
void THIS_TSEMVISUALIZE::read(TSEM<valuetype>* p_tsem, HGeometryTree<3>* hg_tree, Vector<valuetype>* solution)
{
    tsem = p_tsem; h_tree = hg_tree; sol = solution;
    flag_interp = 1;
}
TEMPLATE_TSEMVISUALIZE
void THIS_TSEMVISUALIZE::read(TSEM<valuetype>* p_tsem, HGeometryTree<3>* hg_tree, std::vector<Vector<valuetype> >* wave_function, std::vector<valuetype>* n_occupation)
{
    tsem = p_tsem; h_tree = hg_tree; psi = wave_function; n_occ = n_occupation;
    flag_interp = 2;
}

TEMPLATE_TSEMVISUALIZE
void THIS_TSEMVISUALIZE::build()
{
    irregular_mesh_original = new IrregularMesh<3>;
    irregular_mesh_original->reinit(*h_tree);
    irregular_mesh_original->semiregularize();
    irregular_mesh_original->regularize(false);
    irregular_mesh = new IrregularMesh<3>;
    irregular_mesh->reinit(*h_tree);
    if (n_global_refine > 0) irregular_mesh->globalRefine(n_global_refine);

    // obtain interpolation to linear fem space by local refinement
    for (int n_refine = 0; n_refine <= n_local_refine; ++n_refine){
	std::cerr << "n_refine = " << n_refine << '\n';
	// interpolate tsem solution to fem space on current mesh (represented by value on nodes)
	interp();
	std::cout << "interpolation for " << irregular_mesh->regularMesh().n_geometry(0) << " points done\n";
	if (n_refine == n_local_refine) break;
	// calculate error indicator
	get_indicator();
	std::cerr << "get indicator\n";
	// smooth indicator
	if (n_smooth_indicator > 0){
	    smooth_indicator();	std::cerr << "smooth indicator\n"; }
	valuetype max_indicator = *std::max_element(val_indicator.begin(), val_indicator.end());
	if (max_indicator < tolerance) break;
	// adapt mesh for better interpolation
	RegularMesh<3>& mesh = irregular_mesh->regularMesh();
	Indicator<3> indicator(mesh);
	for (unsigned int idx_ele = 0; idx_ele < mesh.n_geometry(3); ++idx_ele) indicator[idx_ele] = val_indicator[idx_ele];
	MeshAdaptor<3> mesh_adaptor(*irregular_mesh);
	mesh_adaptor.setIndicator(indicator);
	mesh_adaptor.tolerence() = max_indicator * rate_refine / (convergenceCoefficient * refine_threshold);
	mesh_adaptor.is_refine_only() = true;
	mesh_adaptor.adapt();
    }
    // std::cout << "get mesh\n";
}

TEMPLATE_TSEMVISUALIZE
void THIS_TSEMVISUALIZE::interp()
{
    irregular_mesh->semiregularize(); irregular_mesh->regularize(false);
    RegularMesh<3>& mesh = irregular_mesh->regularMesh();
    val_interp.resize(mesh.n_geometry(0));
    std::vector<bool> flag_traversed(mesh.n_geometry(0), false);
    
    // IrregularMeshPair<3> mesh_pair(irregular_mesh, irregular_mesh_original);
    // ActiveElementPairIterator<3> the_pair = mesh_pair.beginActiveElementPair();
    // ActiveElementPairIterator<3> end_pair = mesh_pair.endActiveElementPair();
    // for (; the_pair != end_pair; ++the_pair){
    // 	// std::cout << "index = " << the_pair(0).index << ' ' << the_pair(1).index
    // 	// 	  << ", state: " << the_pair.state() << ", less_than = " << ActiveElementPairIterator<3>::LESS_THAN
    // 	// 	  << '\n';
    // 	if (the_pair.state() != ActiveElementPairIterator<3>::LESS_THAN) continue;
    // 	for (unsigned int idx_p = 0; idx_p < mesh.geometry(3, the_pair(0).index).n_vertex(); ++idx_p){
    // 	    int idx_point = mesh.geometry(3, the_pair(0).index).vertex(idx_p);
    // 	    if (flag_traversed[idx_point]) continue;
    // 	    switch (flag_interp){
    // 	    case 1:
    // 		val_interp[idx_point] = tsem->calc_val_inElement(mesh, the_pair(1).index, mesh.point(idx_point), *sol);
    // 		break;
    // 	    case 2:
    // 		val_interp[idx_point] = tsem->calc_val_inElement(mesh, the_pair(1).index, mesh.point(idx_point), *psi, *n_occ);
    // 		break;
    // 	    }
    // 	    flag_traversed[idx_point] = true;
    // 	}	
    // }

    // for (int i = 0; i < mesh.n_geometry(0); ++i) std::cout << i << ' ' << mesh.point(i) << ' ';
    RegularMesh<3>& mesh_original = irregular_mesh_original->regularMesh();
    for (unsigned int idx_p = 0; idx_p < mesh.n_geometry(0); ++idx_p){
    	switch (flag_interp){
    	case 1:
    	    val_interp[idx_p] = tsem->calc_val_point(mesh_original, mesh.point(idx_p), *sol);
    	    break;
    	case 2:
    	    val_interp[idx_p] = tsem->calc_val_point(mesh_original, mesh.point(idx_p), *psi, *n_occ);
    	    break;
    	}
    }
}

TEMPLATE_TSEMVISUALIZE
void THIS_TSEMVISUALIZE::get_indicator()
{
    val_indicator.resize(irregular_mesh->regularMesh().n_geometry(3), 1);
}

TEMPLATE_TSEMVISUALIZE
void THIS_TSEMVISUALIZE::smooth_indicator()
{
    RegularMesh<3>& mesh = irregular_mesh->regularMesh();
    std::vector<std::vector<int> > idx_vtx = {{0,1,2,3}, {0,1,3,4}};
    std::vector<valuetype> indicator_node(mesh.n_geometry(0));
    std::vector<valuetype> volume(mesh.n_geometry(3));
    std::vector<valuetype> mass_lumping(mesh.n_geometry(0), 0);
    for (unsigned int idx_ele = 0; idx_ele < mesh.n_geometry(3); ++idx_ele){
	int idx = mesh.geometry(3,idx_ele).n_vertex() == 5 ? 1: 0;
	AFEPack::Point<3>& x0 = mesh.point(mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][0]));
	AFEPack::Point<3>& x1 = mesh.point(mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][1]));
	AFEPack::Point<3>& x2 = mesh.point(mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][2]));
	AFEPack::Point<3>& x3 = mesh.point(mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][3]));
	volume[idx_ele] = fabs(((x1[0] - x0[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
				(x1[1] - x0[1])*(x2[2] - x0[2])*(x3[0] - x0[0]) +
				(x1[2] - x0[2])*(x2[0] - x0[0])*(x3[1] - x0[1]) -
				(x1[0] - x0[0])*(x2[2] - x0[2])*(x3[1] - x0[1]) -
				(x1[1] - x0[1])*(x2[0] - x0[0])*(x3[2] - x0[2]) -
				(x1[2] - x0[2])*(x2[1] - x0[1])*(x3[0] - x0[0])));
	for (unsigned int idx_p = 0; idx_p < 4; ++idx_p) mass_lumping[mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][idx_p])] += volume[idx_ele];
    }
    for (unsigned int idx_smooth = 0; idx_smooth < n_smooth_indicator; ++idx_smooth){
	std::fill(indicator_node.begin(), indicator_node.end(), 0);
	for (unsigned int idx_ele = 0; idx_ele < mesh.n_geometry(3); ++idx_ele){
	    int idx = mesh.geometry(3,idx_ele).n_vertex() == 5 ? 1: 0;
	    for (unsigned int idx_p = 0; idx_p < 4; ++idx_p) indicator_node[mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][idx_p])] += val_indicator[idx_ele] * volume[idx_ele];
	}
	for (unsigned int idx_p = 0; idx_p < mesh.n_geometry(0); ++idx_p) indicator_node[idx_p] *= 0.25;
	std::fill(val_indicator.begin(), val_indicator.end(), 0);
	for (unsigned int idx_ele = 0; idx_ele < mesh.n_geometry(3); ++idx_ele){
	    int idx = mesh.geometry(3,idx_ele).n_vertex() == 5 ? 1: 0;
	    for (unsigned int idx_p = 0; idx_p < 4; ++idx_p) val_indicator[idx_ele] += indicator_node[mesh.geometry(3,idx_ele).vertex(idx_vtx[idx][idx_p])];
	}
    }
}

TEMPLATE_TSEMVISUALIZE
void THIS_TSEMVISUALIZE::write(std::string filename_output)
{
    irregular_mesh->semiregularize(); irregular_mesh->regularize(false);
    RegularMesh<3>& mesh = irregular_mesh->regularMesh();
    
    // generate fem_space for output
    std::vector<TemplateElement<double, 3, 3> > template_element(3);
    TemplateGeometry<3> template_geometry_tetrahedron;
    CoordTransform<3, 3> coord_transform_tetrahedron;
    TemplateDOF<3> template_dof_tetrahedron;
    BasisFunctionAdmin<double, 3, 3> basis_function_tetrahedron;
    template_geometry_tetrahedron.readData("tetrahedron.tmp_geo");
    coord_transform_tetrahedron.readData("tetrahedron.crd_trs");
    template_dof_tetrahedron.reinit(template_geometry_tetrahedron);
    template_dof_tetrahedron.readData("tetrahedron.1.tmp_dof");
    basis_function_tetrahedron.reinit(template_dof_tetrahedron);
    basis_function_tetrahedron.readData("tetrahedron.1.bas_fun");
    template_element[0].reinit(template_geometry_tetrahedron,
			       template_dof_tetrahedron,
			       coord_transform_tetrahedron,
			       basis_function_tetrahedron);
    TemplateGeometry<3> template_geometry_twin_tetrahedron;
    CoordTransform<3, 3> coord_transform_twin_tetrahedron;
    TemplateDOF<3> template_dof_twin_tetrahedron;
    BasisFunctionAdmin<double, 3, 3> basis_function_twin_tetrahedron;
    template_geometry_twin_tetrahedron.readData("twin_tetrahedron.tmp_geo");
    coord_transform_twin_tetrahedron.readData("twin_tetrahedron.crd_trs");
    template_dof_twin_tetrahedron.reinit(template_geometry_twin_tetrahedron);
    template_dof_twin_tetrahedron.readData("twin_tetrahedron.1.tmp_dof");
    basis_function_twin_tetrahedron.reinit(template_dof_twin_tetrahedron);
    basis_function_twin_tetrahedron.readData("twin_tetrahedron.1.bas_fun");
    template_element[1].reinit(template_geometry_twin_tetrahedron,
			       template_dof_twin_tetrahedron,
			       coord_transform_twin_tetrahedron,
			       basis_function_twin_tetrahedron);
    TemplateGeometry<3> template_geometry_four_tetrahedron;
    CoordTransform<3, 3> coord_transform_four_tetrahedron;
    TemplateDOF<3> template_dof_four_tetrahedron;
    BasisFunctionAdmin<double, 3, 3> basis_function_four_tetrahedron;
    template_geometry_four_tetrahedron.readData("four_tetrahedron.tmp_geo");
    coord_transform_four_tetrahedron.readData("four_tetrahedron.crd_trs");
    template_dof_four_tetrahedron.reinit(template_geometry_four_tetrahedron);
    template_dof_four_tetrahedron.readData("four_tetrahedron.1.tmp_dof");
    basis_function_four_tetrahedron.reinit(template_dof_four_tetrahedron);
    basis_function_four_tetrahedron.readData("four_tetrahedron.1.bas_fun");
    template_element[2].reinit(template_geometry_four_tetrahedron,
			       template_dof_four_tetrahedron,
			       coord_transform_four_tetrahedron,
			       basis_function_four_tetrahedron);

    // construct finite element space
    FEMSpace<double, 3> fem_space(mesh, template_element);
    int n_ele = mesh.n_geometry(3);
    fem_space.element().resize(n_ele);
    for (unsigned int idx_ele = 0; idx_ele < n_ele; ++idx_ele){
	int idx = mesh.geometry(3,idx_ele).n_vertex()-4;
	if (idx > 1) idx = 2;
    	fem_space.element(idx_ele).reinit(fem_space, idx_ele, idx);
    }
    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();
    std::cerr << "form fem_space for output, with n_dof " << fem_space.n_dof() << '\n';

    FEMFunction<double, 3> output(fem_space);
    double tz = 1.0e-8; // tolerance_zero
    for (unsigned int idx_p = 0; idx_p < fem_space.n_dof(); ++idx_p){
        AFEPack::Point<3>& p = fem_space.dofInfo(idx_p).interp_point;
        if (distance(p, mesh.point(idx_p)) < tz){
            output(idx_p) = val_interp[idx_p]; continue;
        }
        for (unsigned int idx_pm = 0; idx_pm < mesh.n_geometry(0); ++idx_pm)
            if (distance(p, mesh.point(idx_pm)) < tz){
                output(idx_p) = val_interp[idx_pm]; break;
            }
    }
    output.writeOpenDXData(filename_output);
    std::cerr << "write dx data to " << filename_output << '\n';
}

TEMPLATE_TSEMVISUALIZE
void THIS_TSEMVISUALIZE::write(std::string filename_mesh, std::string filename_output)
{
    Mesh<3> mesh;
    mesh.readData(filename_mesh);
    
    TemplateGeometry<3> template_geometry;
    CoordTransform<3, 3> coord_transform;
    TemplateDOF<3> template_dof;
    BasisFunctionAdmin<double, 3, 3> basis_function;
    std::vector<TemplateElement<double, 3, 3> > template_element;
    template_geometry.readData("tetrahedron.tmp_geo");
    coord_transform.readData("tetrahedron.crd_trs");
    template_dof.reinit(template_geometry);     template_dof.readData("tetrahedron.2.tmp_dof");
    basis_function.reinit(template_dof);        basis_function.readData("tetrahedron.2.bas_fun");
    template_element.resize(1);
    template_element[0].reinit(template_geometry, template_dof, coord_transform, basis_function);

    // fem space
    FEMSpace<double, 3> fem_space;
    fem_space.reinit(mesh, template_element);
    int n_ele = mesh.n_geometry(3);
    fem_space.element().resize(n_ele);
    for (unsigned int idx_ele = 0; idx_ele < n_ele; ++idx_ele) fem_space.element(idx_ele).reinit(fem_space, idx_ele, 0);
    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();
    int n_dof = fem_space.n_dof();

    // interpolation
    FEMFunction<double, 3> output(fem_space);
    for (unsigned int idx_dof = 0; idx_dof < n_dof; ++idx_dof){
	AFEPack::Point<3>& p = fem_space.dofInfo(idx_dof).interp_point;
	switch (flag_interp){
	case 1:
	    output(idx_dof) = tsem->calc_val_point(mesh, p, *sol);
	    break;
	case 2:
	    output(idx_dof) = tsem->calc_val_point(mesh, p, *psi, *n_occ);
	    break;
	}
    }
    output.writeOpenDXData(filename_output);
}


#undef TEMPLATE_TSEMVISUALIZE
#undef THIS_TSEMVISUALIZE
