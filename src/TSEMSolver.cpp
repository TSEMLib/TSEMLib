#include "../include/TSEMSolver.h"
// This file contains functions for solving the linear system from TSEM discretizations
// Functions:
// init_PMG()
// build_restrict_identity()
// build_restrict_approximateL2Interpolation()
// build_interpolate_identity()
// build_interpolate_approximateL2Interpolation()
// build_interpolate_approximateL2InterpolationTransfer()
// build_M_exactPartial()
// build_M_PMPT()
// get_PMPT()
// ite_GaussSeidel()
// ite_BackwardGaussSeidel()
// ite_SymmetricGaussSeidel()
// ite_PMG()
// ite_Jacobi()
// ite_SOR()
// ite_SSOR()
// solve_PMG()
// solve_GaussSeidel()
// solve_CG()
// solve_PCG()
// solve_InversePower()
// solve_LOBPCG()
// calc_inv()


template class TSEMSolver<double>;

#define TEMPLATE_TSEMSOLVER template<typename valuetype>
#define THIS_TSEMSOLVER TSEMSolver<valuetype>


TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::init_PMG(const SparseMatrix<valuetype> &sp_matrix,
				     const SparseMatrix<valuetype> &mass, int polynomial_order, std::vector<unsigned int> & sum_dof)
{
  std::cerr << "PMGSolver initializing ..." << '\n';
    
  // setup
  int n_dof_total = sp_matrix.m();
  int n_dof_this_layer = sp_matrix.m();
  int n_order_each_layer = 1;

  // n_project = polynomial_order - 1;
  n_project = (polynomial_order+n_order_each_layer-1)/n_order_each_layer - 1
    + ((polynomial_order-1)%n_order_each_layer == 0 ? 0 : 1);
  std::vector<int> polynomial_order_each_layer(n_project+1);
  polynomial_order_each_layer[0] = polynomial_order;
  polynomial_order_each_layer[n_project] = 1;
  // for (unsigned int i = n_project-1; i >= 1; --i)
  // 	polynomial_order_each_layer[i] = polynomial_order_each_layer[i+1] + n_order_each_layer;
  for (unsigned int i = 1; i < n_project; ++i)
    polynomial_order_each_layer[i] = polynomial_order_each_layer[i-1] - n_order_each_layer;

  // std::vector<unsigned int> sum_dof_tmp(sum_dof.size());
  // for (unsigned int i = 0; i < sum_dof.size(); ++i)
  // 	sum_dof_tmp[i] = sum_dof[i];
  // unsigned int tmp_po = polynomial_order - 1;
  // n_project = 0;
  // while (tmp_po > 0){
  // 	tmp_po /= 2;
  // 	n_project++;
  // }
    
  location_dof.resize(n_project);
  location_dof_inv.resize(n_project);
  flag_fine_only.resize(n_project);
  n_dof_layer.resize(polynomial_order);
  n_dof_layer[0] = n_dof_total;
    
  // calculate basic info
  core_at.resize(n_project);
  for (unsigned int ind_project = 0; ind_project < n_project; ++ind_project){ // except the last layer
    // int polynomial_order_fine   = polynomial_order - ind_project; // select out the dof with multiindex whose sum <= polynomial_order_fine/coarse
    // int polynomial_order_coarse = polynomial_order - ind_project - 1; // summation of multiindex on coarse layer
    // int polynomial_order_fine = ind_project == 0 ? 0 : ind_project + 1; 
    // int polynomial_order_coarse   = ind_project + 2;
    // int polynomial_order_fine   = polynomial_order - ind_project*n_order_each_layer; // select out the dof with multiindex whose sum <= polynomial_order_fine/coarse
    // int polynomial_order_coarse = polynomial_order - ind_project*n_order_each_layer - n_order_each_layer; // summation of multiindex on coarse layer
    int polynomial_order_fine   = polynomial_order_each_layer[ind_project];
    int polynomial_order_coarse = polynomial_order_each_layer[ind_project+1];
	
    int n_dof_next_layer = 0;
    int p_core_at = -1; // pointer for core_at, as we assign its value according to total dof
    core_at[ind_project].resize(n_dof_this_layer, -1);

    location_dof[ind_project].resize(n_dof_this_layer);
    location_dof_inv[ind_project].resize(n_dof_total, -1);
    flag_fine_only[ind_project].resize(n_dof_this_layer, false);
        
    for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof){
      bool flag_fine   = sum_dof[ind_dof] <= polynomial_order_fine;
      bool flag_coarse = sum_dof[ind_dof] <= polynomial_order_coarse;
      // bool flag_fine   = sum_dof[ind_dof] >= polynomial_order_fine;
      // bool flag_coarse = sum_dof[ind_dof] >= polynomial_order_coarse;
      // bool flag_fine   = sum_dof_tmp[ind_dof] > 0;
      // bool flag_coarse = sum_dof_tmp[ind_dof] > 0 && sum_dof_tmp[ind_dof]%2 == 1;
	    
      if (flag_fine)   p_core_at++;
      if (flag_coarse) core_at[ind_project][p_core_at] = n_dof_next_layer++;
            
      if (flag_fine){
	location_dof[ind_project][p_core_at] = ind_dof;
	location_dof_inv[ind_project][ind_dof] = p_core_at;
      }
      if (flag_fine && !flag_coarse)
	flag_fine_only[ind_project][p_core_at] = true;
    }
    n_dof_this_layer = n_dof_next_layer;
    n_dof_layer[ind_project+1] = n_dof_this_layer;
    std::cerr << "\tfind " << n_dof_next_layer << " points in coarse grid\n";

    // for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
    //     if (sum_dof_tmp[ind_dof] % 2 == 0)
    // 	sum_dof_tmp[ind_dof] = 0;
    //     else
    // 	sum_dof_tmp[ind_dof] = (sum_dof_tmp[ind_dof]+1) / 2;
  }

    
  // build restrict operator P, interpolate operator PT and matrix M
  build_restrict_identity();
  // build_restrict_approximateL2Interpolation(mass);
  build_interpolate_identity();
  // build_interpolate_approximateL2Interpolation(mass);
  // build_interpolate_approximateL2InterpolationTransfer(mass);
  build_M_exactPartial(sp_matrix);
  // build_M_PMPT(sp_matrix);
    
  // solver.reinit(*projected_matrix[n_project]);
    
  std::cerr << "\tOK! grid levels: " << n_project << '\n';
}


TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::init_PMG(const SparseMatrix<valuetype> &D, const SparseMatrix<valuetype> &L,
			       const SparseMatrix<valuetype> &mass, int polynomial_order, std::vector<unsigned int> &sum_dof)
{
    std::cerr << "PMGSolver initializing ..." << '\n';
    
    // setup
    int n_dof_total = D.m();
    int n_dof_this_layer = D.m();
    int n_order_each_layer = 1;

    // n_project = polynomial_order - 1;
    n_project = (polynomial_order+n_order_each_layer-1)/n_order_each_layer - 1
    	+ ((polynomial_order-1)%n_order_each_layer == 0 ? 0 : 1);
    std::vector<int> polynomial_order_each_layer(n_project+1);
    polynomial_order_each_layer[0] = polynomial_order;
    polynomial_order_each_layer[n_project] = 1;
    for (unsigned int i = 1; i < n_project; ++i)
    	polynomial_order_each_layer[i] = polynomial_order_each_layer[i-1] - n_order_each_layer;
    
    location_dof.resize(n_project);
    location_dof_inv.resize(n_project);
    flag_fine_only.resize(n_project);
    n_dof_layer.resize(polynomial_order);
    n_dof_layer[0] = n_dof_total;
    
    // calculate basic info
    core_at.resize(n_project);
    for (unsigned int ind_project = 0; ind_project < n_project; ++ind_project){ // except the last layer
	int polynomial_order_fine   = polynomial_order_each_layer[ind_project];
	int polynomial_order_coarse = polynomial_order_each_layer[ind_project+1];
	
        int n_dof_next_layer = 0;
        int p_core_at = -1; // pointer for core_at, as we assign its value according to total dof
        core_at[ind_project].resize(n_dof_this_layer, -1);

        location_dof[ind_project].resize(n_dof_this_layer);
        location_dof_inv[ind_project].resize(n_dof_total, -1);
        flag_fine_only[ind_project].resize(n_dof_this_layer, false);
        
        for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof){
            bool flag_fine   = sum_dof[ind_dof] <= polynomial_order_fine;
            bool flag_coarse = sum_dof[ind_dof] <= polynomial_order_coarse;
	    
            if (flag_fine)   p_core_at++;
            if (flag_coarse) core_at[ind_project][p_core_at] = n_dof_next_layer++;
            
            if (flag_fine){
                location_dof[ind_project][p_core_at] = ind_dof;
                location_dof_inv[ind_project][ind_dof] = p_core_at;
            }
            if (flag_fine && !flag_coarse)
                flag_fine_only[ind_project][p_core_at] = true;
        }
        n_dof_this_layer = n_dof_next_layer;
	n_dof_layer[ind_project+1] = n_dof_this_layer;
        std::cerr << "\tfind " << n_dof_next_layer << " points in coarse grid\n";
    }

    
    // build restrict operator P, interpolate operator PT and matrix M
    build_restrict_identity();
    build_interpolate_identity();
    build_M_exactPartial(D, L);
    
    std::cerr << "\tOK! grid levels: " << n_project << '\n';
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::build_restrict_identity()
{
    
    // count nonzero entries
    std::vector<std::vector<unsigned int> > n_nonzero(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	n_nonzero[ind_project].resize(n_dof_layer[ind_project+1], 2);


    // build sparsity pattern
    // std::vector<SparsityPattern&> sp_pattern(n_project);
    sp_pattern_restrict.resize(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	// sp_pattern[ind_project] = *(new SparsityPattern(n_dof_layer[ind_project], n_dof_layer[ind_project+1], n_nonzero[ind_project]));
	sp_pattern_restrict[ind_project].reinit(n_dof_layer[ind_project+1], n_dof_layer[ind_project], n_nonzero[ind_project]);
    // add coarse dof on coarse -> fine
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_pattern_restrict[ind_project].add(core_at[ind_project][ind_dof_fine], ind_dof_fine);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
        sp_pattern_restrict[ind_project].compress();


    // assign sparse matrix
    std::vector<SparseMatrix<valuetype>*> sp_matrix(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_matrix[ind_project] = new SparseMatrix<valuetype>(sp_pattern_restrict[ind_project]);
    // add coarse dof on fine -> coarse
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_matrix[ind_project]->add(core_at[ind_project][ind_dof_fine], ind_dof_fine, (valuetype) 1);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	restrict_matrix.push_back(sp_matrix[ind_project]);


    std::cerr << "build restrict matrix\n";
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::build_restrict_approximateL2Interpolation(const SparseMatrix<valuetype> &mass)
{
    const SparsityPattern &spM = mass.get_sparsity_pattern();
    const std::size_t *row_start = spM.get_rowstart_indices();
    const unsigned int *col_nums = spM.get_column_numbers();
    int n_dof_total = mass.m();
    // std::cout << "mass.m() = " << n_dof_total << '\n';

    
    // count nonzero entries
    std::vector<std::vector<unsigned int> > n_nonzero(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	n_nonzero[ind_project].resize(n_dof_layer[ind_project+1], 2);
    for (unsigned int row = 0; row < mass.m(); ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    // if row = find only, col = coarse, n_nnz++
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] != -1) continue; // fine grid only
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution
		n_nonzero[ind_project][core_at[ind_project][location_dof_inv[ind_project][col]]]++;
	    }
	}
    // std::cerr << "count nonzero entries\n";


    // build sparsity pattern
    // std::vector<SparsityPattern&> sp_pattern(n_project);
    sp_pattern_restrict.resize(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	// sp_pattern[ind_project] = *(new SparsityPattern(n_dof_layer[ind_project+1], n_dof_layer[ind_project], n_nonzero[ind_project]));
	sp_pattern_restrict[ind_project].reinit(n_dof_layer[ind_project+1], n_dof_layer[ind_project], n_nonzero[ind_project]);
    // add coarse dof on fine -> coarse
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_pattern_restrict[ind_project].add(core_at[ind_project][ind_dof_fine], ind_dof_fine);
    // add approximate L2 interpolation of DoFs with higher polynomial order
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] != -1) continue; // fine grid only
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_pattern_restrict[ind_project].add(core_at[ind_project][location_dof_inv[ind_project][col]],
						     location_dof_inv[ind_project][row]);
	    }
	}
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_pattern_restrict[ind_project].compress();
    // std::cerr << "build sparse pattern\n";


    // assign sparse matrix
    std::vector<SparseMatrix<valuetype>*> sp_matrix(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_matrix[ind_project] = new SparseMatrix<valuetype>(sp_pattern_restrict[ind_project]);
    // add coarse dof on fine -> coarse
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_matrix[ind_project]->add(core_at[ind_project][ind_dof_fine], ind_dof_fine, (valuetype) 1);
    // add approximate L2 interpolation of DoFs with higher polynomial order
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] != -1) continue; // fine grid only
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_matrix[ind_project]->add(core_at[ind_project][location_dof_inv[ind_project][col]],
					    location_dof_inv[ind_project][row],
					    mass.global_entry(pos_whole) / mass.diag_element(row));
	    }
	}
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	restrict_matrix.push_back(sp_matrix[ind_project]);


    std::cerr << "build restrict matrix\n";
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::build_interpolate_identity()
{
    // count nonzero entries
    std::vector<std::vector<unsigned int> > n_nonzero(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	n_nonzero[ind_project].resize(n_dof_layer[ind_project], 2);


    // build sparsity pattern
    // std::vector<SparsityPattern&> sp_pattern(n_project);
    sp_pattern_interpolate.resize(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	// sp_pattern[ind_project] = *(new SparsityPattern(n_dof_layer[ind_project], n_dof_layer[ind_project+1], n_nonzero[ind_project]));
	sp_pattern_interpolate[ind_project].reinit(n_dof_layer[ind_project], n_dof_layer[ind_project+1], n_nonzero[ind_project]);
    // add coarse dof on fine -> coarse
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_pattern_interpolate[ind_project].add(ind_dof_fine, core_at[ind_project][ind_dof_fine]);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_pattern_interpolate[ind_project].compress();


    // assign sparse matrix
    std::vector<SparseMatrix<valuetype>*> sp_matrix(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_matrix[ind_project] = new SparseMatrix<valuetype>(sp_pattern_interpolate[ind_project]);
    // add coarse dof on fine -> coarse
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_matrix[ind_project]->add(ind_dof_fine, core_at[ind_project][ind_dof_fine], (valuetype) 1);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	interpolate_matrix.push_back(sp_matrix[ind_project]);


    std::cerr << "build interpolate matrix\n";
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::build_interpolate_approximateL2Interpolation(const SparseMatrix<valuetype> &mass)
{
    const SparsityPattern &spM = mass.get_sparsity_pattern();
    const std::size_t *row_start = spM.get_rowstart_indices();
    const unsigned int *col_nums = spM.get_column_numbers();
    int n_dof_total = mass.m();

    
    // count nonzero entries
    std::vector<std::vector<unsigned int> > n_nonzero(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	n_nonzero[ind_project].resize(n_dof_layer[ind_project], 2);
    for (unsigned int row = 0; row < mass.m(); ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    // if row = coarse, col = find only, n_nnz++
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] != -1) continue; // fine grid only
		// add contribution
		n_nonzero[ind_project][location_dof_inv[ind_project][col]]++;
	    }
	}


    // build sparsity pattern
    // std::vector<SparsityPattern&> sp_pattern(n_project);
    sp_pattern_interpolate.resize(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	// sp_pattern[ind_project] = *(new SparsityPattern(n_dof_layer[ind_project+1], n_dof_layer[ind_project], n_nonzero[ind_project]));
	sp_pattern_interpolate[ind_project].reinit(n_dof_layer[ind_project], n_dof_layer[ind_project+1], n_nonzero[ind_project]);
    // add coarse dof on coarse -> fine
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_pattern_interpolate[ind_project].add(ind_dof_fine, core_at[ind_project][ind_dof_fine]);
    // add approximate L2 interpolation of DoFs with higher polynomial order
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] != -1) continue; // fine grid only
		// add contribution of n_nonzero
		sp_pattern_interpolate[ind_project].add(location_dof_inv[ind_project][col],
							core_at[ind_project][location_dof_inv[ind_project][row]]);
	    }
	}
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_pattern_interpolate[ind_project].compress();


    // assign sparse matrix
    std::vector<SparseMatrix<valuetype>*> sp_matrix(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_matrix[ind_project] = new SparseMatrix<valuetype>(sp_pattern_interpolate[ind_project]);
    // add coarse dof on coarse -> fine
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_matrix[ind_project]->add(ind_dof_fine, core_at[ind_project][ind_dof_fine], (valuetype) 1);
    // add approximate L2 interpolation of DoFs with higher polynomial order
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] != -1) continue; // fine grid only
		// add contribution of n_nonzero
		sp_matrix[ind_project]->add(location_dof_inv[ind_project][col],
					    core_at[ind_project][location_dof_inv[ind_project][row]],
					    mass.global_entry(pos_whole) / mass.diag_element(row));
	    }
	}
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	interpolate_matrix.push_back(sp_matrix[ind_project]);


    std::cerr << "build interpolate matrix\n";
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::build_interpolate_approximateL2InterpolationTransfer(const SparseMatrix<valuetype> &mass)
{
    const SparsityPattern &spM = mass.get_sparsity_pattern();
    const std::size_t *row_start = spM.get_rowstart_indices();
    const unsigned int *col_nums = spM.get_column_numbers();
    int n_dof_total = mass.m();

    
    // count nonzero entries
    std::vector<std::vector<unsigned int> > n_nonzero(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	n_nonzero[ind_project].resize(n_dof_layer[ind_project], 2);
    for (unsigned int row = 0; row < mass.m(); ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    // if row = find only, col = coarse, n_nnz++
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] != -1) continue; // fine grid only
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution
		n_nonzero[ind_project][location_dof_inv[ind_project][row]]++;
	    }
	}


    // build sparsity pattern
    // std::vector<SparsityPattern&> sp_pattern(n_project);
    sp_pattern_interpolate.resize(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	// sp_pattern[ind_project] = *(new SparsityPattern(n_dof_layer[ind_project+1], n_dof_layer[ind_project], n_nonzero[ind_project]));
	sp_pattern_interpolate[ind_project].reinit(n_dof_layer[ind_project], n_dof_layer[ind_project+1], n_nonzero[ind_project]);
    // add coarse dof on fine -> coarse
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_pattern_interpolate[ind_project].add(ind_dof_fine, core_at[ind_project][ind_dof_fine]);
    // add approximate L2 interpolation of DoFs with higher polynomial order
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] != -1) continue; // fine grid only
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_pattern_interpolate[ind_project].add(location_dof_inv[ind_project][row],
							core_at[ind_project][location_dof_inv[ind_project][col]]);
	    }
	}
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_pattern_interpolate[ind_project].compress();


    // assign sparse matrix
    std::vector<SparseMatrix<valuetype>*> sp_matrix(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_matrix[ind_project] = new SparseMatrix<valuetype>(sp_pattern_interpolate[ind_project]);
    // add coarse dof on fine -> coarse
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	for (int ind_dof_fine = 0; ind_dof_fine < n_dof_layer[ind_project]; ++ind_dof_fine)
	    if (core_at[ind_project][ind_dof_fine] != -1)
		sp_matrix[ind_project]->add(ind_dof_fine, core_at[ind_project][ind_dof_fine], (valuetype) 1);
    // add approximate L2 interpolation of DoFs with higher polynomial order
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] != -1) continue; // fine grid only
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_matrix[ind_project]->add(location_dof_inv[ind_project][row],
					    core_at[ind_project][location_dof_inv[ind_project][col]],
					    mass.global_entry(pos_whole) / mass.diag_element(row));
	    }
	}
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	interpolate_matrix.push_back(sp_matrix[ind_project]);


    std::cerr << "build interpolate matrix\n";
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::build_M_exactPartial(const SparseMatrix<valuetype> &sp_matrix_finest)
{// assign all matrix simutaneously
 // index [ind = 0: n_project-1] correspond to matrix for DoFs with polynomial order <= polynomial_order - (ind+1)
 // e.g., sp_matrix[1] is exactly the matrix for DoFs with polynomial order <= polynomial_order - 2
    
    int n_dof_total = sp_matrix_finest.m();
    projected_matrix.push_back(&sp_matrix_finest);
    const SparsityPattern &sp_pattern_finest = sp_matrix_finest.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern_finest.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern_finest.get_column_numbers();
    
    
    // count nonzero entries
    std::vector<std::vector<unsigned int> > n_nonzero(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	n_nonzero[ind_project].resize(n_dof_layer[ind_project+1]);
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // at least on fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue;
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // at least on fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue;
		// add contribution of n_nonzero
		n_nonzero[ind_project][core_at[ind_project][location_dof_inv[ind_project][row]]]++;
	    }
	}


    // build sparsity pattern
    // std::vector<SparsityPattern&> sp_pattern(n_project);
    sp_pattern_projected_matrix.resize(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	// sp_pattern[ind_project] = *(new SparsityPattern(n_dof_layer[ind_project+1], n_dof_layer[ind_project+1], n_nonzero[ind_project]));
	sp_pattern_projected_matrix[ind_project].reinit(n_dof_layer[ind_project+1], n_dof_layer[ind_project+1], n_nonzero[ind_project]);
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_pattern_projected_matrix[ind_project].add(core_at[ind_project][location_dof_inv[ind_project][row]],
							     core_at[ind_project][location_dof_inv[ind_project][col]]);
	    }
	}
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_pattern_projected_matrix[ind_project].compress();
    

    // assign sparse matrix
    std::vector<SparseMatrix<valuetype>*> sp_matrix(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	sp_matrix[ind_project] = new SparseMatrix<valuetype>(sp_pattern_projected_matrix[ind_project]);
    for (unsigned int row = 0; row < n_dof_total; ++row)
	for (unsigned int pos_whole = row_start[row]; pos_whole < row_start[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_matrix[ind_project]->add(core_at[ind_project][location_dof_inv[ind_project][row]],
					    core_at[ind_project][location_dof_inv[ind_project][col]],
					    sp_matrix_finest.global_entry(pos_whole));
	    }
	}
    for (int ind_project = 0; ind_project < n_project; ++ind_project)
	projected_matrix.push_back(sp_matrix[ind_project]);


    std::cerr << "build matrices on each grid\n";
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::build_M_exactPartial(const SparseMatrix<valuetype> &sp_matrix_D_finest, const SparseMatrix<valuetype> &sp_matrix_L_finest)
{// assign all matrix simutaneously
 // index [ind = 0: n_project-1] correspond to matrix for DoFs with polynomial order <= polynomial_order - (ind+1)
 // e.g., sp_matrix[1] is exactly the matrix for DoFs with polynomial order <= polynomial_order - 2
    
    int n_dof_total = sp_matrix_D_finest.m();
    projected_matrix_D.push_back(&sp_matrix_D_finest);
    projected_matrix_L.push_back(&sp_matrix_L_finest);
    const SparsityPattern &sp_pattern_D_finest = sp_matrix_D_finest.get_sparsity_pattern();
    const std::size_t *row_start_D = sp_pattern_D_finest.get_rowstart_indices();
    const unsigned int *col_nums_D = sp_pattern_D_finest.get_column_numbers();
    const SparsityPattern &sp_pattern_L_finest = sp_matrix_L_finest.get_sparsity_pattern();
    const std::size_t *row_start_L = sp_pattern_L_finest.get_rowstart_indices();
    const unsigned int *col_nums_L = sp_pattern_L_finest.get_column_numbers();
    
    
    // count nonzero entries
    std::vector<std::vector<unsigned int> > n_nonzero_D(n_project), n_nonzero_L(n_project);
    for (unsigned int ind_project = 0; ind_project < n_project; ++ind_project){
	n_nonzero_D[ind_project].resize(n_dof_layer[ind_project+1]);
	n_nonzero_L[ind_project].resize(n_dof_layer[ind_project+1]);
    }
    for (unsigned int row = 0; row < n_dof_total; ++row){
	for (unsigned int pos_whole = row_start_D[row]; pos_whole < row_start_D[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums_D[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // at least on fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue;
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // at least on fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue;
		// add contribution of n_nonzero
		n_nonzero_D[ind_project][core_at[ind_project][location_dof_inv[ind_project][row]]]++;
	    }
	}
	for (unsigned int pos_whole = row_start_L[row]; pos_whole < row_start_L[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums_L[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // at least on fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue;
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // at least on fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue;
		// add contribution of n_nonzero
		n_nonzero_L[ind_project][core_at[ind_project][location_dof_inv[ind_project][row]]]++;
	    }
	}
    }

    // build sparsity pattern
    sp_pattern_projected_matrix_D.resize(n_project);
    sp_pattern_projected_matrix_L.resize(n_project);
    for (unsigned int ind_project = 0; ind_project < n_project; ++ind_project){
	sp_pattern_projected_matrix_D[ind_project].reinit(n_dof_layer[ind_project+1], n_dof_layer[ind_project+1], n_nonzero_D[ind_project]);
	sp_pattern_projected_matrix_L[ind_project].reinit(n_dof_layer[ind_project+1], n_dof_layer[ind_project+1], n_nonzero_L[ind_project]);
    }
    for (unsigned int row = 0; row < n_dof_total; ++row){
	for (unsigned int pos_whole = row_start_D[row]; pos_whole < row_start_D[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums_D[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_pattern_projected_matrix_D[ind_project].add(core_at[ind_project][location_dof_inv[ind_project][row]],
							       core_at[ind_project][location_dof_inv[ind_project][col]]);
	    }
	}
	for (unsigned int pos_whole = row_start_L[row]; pos_whole < row_start_L[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums_L[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_pattern_projected_matrix_L[ind_project].add(core_at[ind_project][location_dof_inv[ind_project][row]],
							       core_at[ind_project][location_dof_inv[ind_project][col]]);
	    }
	}
    }
    for (int ind_project = 0; ind_project < n_project; ++ind_project){
	sp_pattern_projected_matrix_D[ind_project].compress();
	sp_pattern_projected_matrix_L[ind_project].compress();
    }
    

    // assign sparse matrix
    std::vector<SparseMatrix<valuetype>*> sp_matrix_D(n_project);
    std::vector<SparseMatrix<valuetype>*> sp_matrix_L(n_project);
    for (int ind_project = 0; ind_project < n_project; ++ind_project){
	sp_matrix_D[ind_project] = new SparseMatrix<valuetype>(sp_pattern_projected_matrix_D[ind_project]);
	sp_matrix_L[ind_project] = new SparseMatrix<valuetype>(sp_pattern_projected_matrix_L[ind_project]);
    }
    for (unsigned int row = 0; row < n_dof_total; ++row){
	for (unsigned int pos_whole = row_start_D[row]; pos_whole < row_start_D[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums_D[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_matrix_D[ind_project]->add(core_at[ind_project][location_dof_inv[ind_project][row]],
					      core_at[ind_project][location_dof_inv[ind_project][col]],
					      sp_matrix_D_finest.global_entry(pos_whole));
	    }
	}
	for (unsigned int pos_whole = row_start_L[row]; pos_whole < row_start_L[row+1]; ++pos_whole){
	    const unsigned int &col = col_nums_L[pos_whole];
	    for (int ind_project = 0; ind_project < n_project; ++ind_project){
		// judge whether row is on this layer
		if (location_dof_inv[ind_project][row] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][row]] == -1) continue; // fine grid -> coarse grid
		// judge whether col is on this layer
		if (location_dof_inv[ind_project][col] == -1) continue; // whole -> fine grid
		if (core_at[ind_project][location_dof_inv[ind_project][col]] == -1) continue; // fine grid -> coarse grid
		// add contribution of n_nonzero
		sp_matrix_L[ind_project]->add(core_at[ind_project][location_dof_inv[ind_project][row]],
					      core_at[ind_project][location_dof_inv[ind_project][col]],
					      sp_matrix_L_finest.global_entry(pos_whole));
	    }
	}
    }
    for (int ind_project = 0; ind_project < n_project; ++ind_project){
	projected_matrix_D.push_back(sp_matrix_D[ind_project]);
	projected_matrix_L.push_back(sp_matrix_L[ind_project]);
    }


    std::cerr << "build matrices on each grid\n";
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::build_M_PMPT(const SparseMatrix<valuetype> &sp_matrix_finest)
{// construct projected_matrix according to restrict_matrix and interpolate_matrix
    projected_matrix.push_back(&sp_matrix_finest);
    std::cerr << "initial matrix has size " << sp_matrix_finest.m() << " and sparsity " << sp_matrix_finest.n_nonzero_elements() / pow(sp_matrix_finest.m(), 2) << '\n';

    for (unsigned int ind_project = 0; ind_project < n_project; ++ind_project){
	SparseMatrix<valuetype> *PMPT;
	PMPT = get_PMPT(*restrict_matrix[ind_project], *projected_matrix[ind_project], *interpolate_matrix[ind_project]);
	projected_matrix.push_back(PMPT);
	std::cerr << "build matrix on " << ind_project << "-th grid, with size " << PMPT->m() << " and sparsity " << PMPT->n_nonzero_elements() / pow(PMPT->m(), 2) << "\n";

	// valuetype max_rate = (valuetype) -1;
	// const std::size_t *row_start = PMPT->get_sparsity_pattern().get_rowstart_indices();
	// for (unsigned int row = 0; row < PMPT->m(); ++row){
	//     valuetype sum_row = (valuetype) 0;
	//     for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
	// 	sum_row += fabs(PMPT->global_entry(pos_whole));
	//     valuetype rat_tmp = sum_row / fabs(PMPT->diag_element(row));
	//     if (max_rate < rat_tmp) max_rate = rat_tmp;
	// }
	// valuetype min_rate = max_rate;
	// for (unsigned int row = 0; row < PMPT->m(); ++row){
	//     valuetype sum_row = (valuetype) 0;
	//     for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
	// 	sum_row += fabs(PMPT->global_entry(pos_whole));
	//     valuetype rat_tmp = sum_row / fabs(PMPT->diag_element(row));
	//     if (fabs(rat_tmp) < tol_zero) continue;
	//     if (min_rate > rat_tmp) min_rate = rat_tmp;
	// }
	// valuetype total_rat = (valuetype) 0;
	// unsigned int count = 0;
	// for (unsigned int row = 0; row < PMPT->m(); ++row){
	//     valuetype sum_row = (valuetype) 0;
	//     for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
	// 	sum_row += fabs(PMPT->global_entry(pos_whole));
	//     valuetype rat_tmp = sum_row / fabs(PMPT->diag_element(row));
	//     if (fabs(rat_tmp) < tol_zero) continue;
	//     total_rat += rat_tmp;
	//     count++;
	// }
	// valuetype avg_rat = total_rat / count;
	// unsigned int count_less = 0;
	// for (unsigned int row = 0; row < PMPT->m(); ++row){
	//     valuetype sum_row = (valuetype) 0;
	//     for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
	// 	sum_row += fabs(PMPT->global_entry(pos_whole));
	//     valuetype rat_tmp = sum_row / fabs(PMPT->diag_element(row));
	//     if (fabs(rat_tmp) < tol_zero) continue;
	//     if (rat_tmp < avg_rat) count_less++;
	// }
	// std::cerr << "\tradio between sum of |un-diangonal entries| and diagonal one has minimum " << min_rate << " and maximum " << max_rate << ", average : " << avg_rat << " with proportion " << count_less * 1.0 / count  << '\n';
    }
    
    std::cerr << "build matrices on each grid\n";
}

TEMPLATE_TSEMSOLVER
SparseMatrix<valuetype> *THIS_TSEMSOLVER::get_PMPT(const SparseMatrix<valuetype> &P, const SparseMatrix<valuetype> &M, const SparseMatrix<valuetype> &PT) const
{    
    const SparsityPattern& spP = P.get_sparsity_pattern();
    const SparsityPattern& spM = M.get_sparsity_pattern();
    const SparsityPattern& spPT = PT.get_sparsity_pattern();
    const std::size_t * P_rowstart = spP.get_rowstart_indices();
    const unsigned int * P_colnums = spP.get_column_numbers();
    const std::size_t * M_rowstart = spM.get_rowstart_indices();
    const unsigned int * M_colnums = spM.get_column_numbers();
    const std::size_t * PT_rowstart = spPT.get_rowstart_indices();
    const unsigned int * PT_colnums = spPT.get_column_numbers();
    std::vector<unsigned int> row_length(P.m(), 0);
    std::vector<bool> flag(P.m(), true);
    std::vector<unsigned int> index(P.m());
    std::vector<std::vector<unsigned int> > col_index(P.m());
    for (unsigned int i = 0; i < P.m(); ++i){
        row_length[i] = 1; /**< add the diagonal entry at first */
        flag[i] = false;
        index[0] = i;
        for (unsigned int j = P_rowstart[i]; j < P_rowstart[i+1]; ++j) {
            const unsigned int& a = P_colnums[j];
            for (unsigned int k = M_rowstart[a]; k < M_rowstart[a+1]; ++k) {
                const unsigned int& b = M_colnums[k];
                for (unsigned int l = PT_rowstart[b]; l < PT_rowstart[b+1]; ++l) {
                    const unsigned int&  c = PT_colnums[l];
                    if (flag[c]) {
                        index[row_length[i]++] = c;
                        flag[c] = false;
                    }
                }
            }
        }
        col_index[i].resize(row_length[i]);
        for (unsigned int j = 0; j < row_length[i]; ++j){
            col_index[i][j] = index[j];
            flag[index[j]] = true;
        }
    }
    SparsityPattern& spA = *(new SparsityPattern(P.m(), row_length));
    for (unsigned int i = 0; i < P.m(); ++i)
        for (unsigned int j = 0; j < row_length[i]; ++j)
            spA.add(i, col_index[i][j]);
    spA.compress();
    SparseMatrix<valuetype> *A = new SparseMatrix<valuetype>(spA);
    std::vector<valuetype> row_entry(P.m(), 0.0);
    for (unsigned int i = 0; i < P.m(); ++i){
        for (unsigned int j = P_rowstart[i]; j < P_rowstart[i+1]; ++j){
            const unsigned int& a = P_colnums[j];
            for (unsigned int k = M_rowstart[a]; k < M_rowstart[a+1]; ++k){
                const unsigned int& b = M_colnums[k];
                for (unsigned int l = PT_rowstart[b]; l < PT_rowstart[b+1]; ++l){
                    const unsigned int& c = PT_colnums[l];
                    row_entry[c] += P.global_entry(j) * M.global_entry(k) * PT.global_entry(l);
                }
            }
        }
        for (unsigned int j = 0; j < row_length[i]; ++j){
            A->add(i, col_index[i][j], row_entry[col_index[i][j]]);
            row_entry[col_index[i][j]] = 0.0;
        }
    }
    
    return A;
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_GaussSeidel(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, const int &step) const
{
    const SparsityPattern& sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    for (unsigned int n_ite = 0; n_ite < step; ++n_ite){
        for (unsigned int row = 0; row < sp_matrix.m(); ++row){
            valuetype r0 = r(row);
            for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
                r0 -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
            x(row) = r0 / sp_matrix.global_entry(row_start[row]);
        }
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_GaussSeidel(const SparseMatrix<valuetype> &D, const SparseMatrix<valuetype> &L, std::vector<Vector<valuetype> > &x, const std::vector<Vector<valuetype> > &r, const int &step) const
{
    int n_dof = D.m();
    const SparsityPattern& sp_pattern_D = D.get_sparsity_pattern();
    const std::size_t *row_start_D = sp_pattern_D.get_rowstart_indices();
    const unsigned int *col_nums_D = sp_pattern_D.get_column_numbers();
    const SparsityPattern& sp_pattern_L = L.get_sparsity_pattern();
    const std::size_t *row_start_L = sp_pattern_L.get_rowstart_indices();
    const unsigned int *col_nums_L = sp_pattern_L.get_column_numbers();
    for (unsigned int n_ite = 0; n_ite < step; ++n_ite){
        for (unsigned int row = 0; row < n_dof; ++row){
	    std::vector<valuetype> r0{r[0](row), r[1](row)};
            for (unsigned int pos_whole = row_start_D[row]+1; pos_whole < row_start_D[row+1]; ++pos_whole){
                r0[0] -= x[0](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
                r0[1] -= x[1](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
	    }
            for (unsigned int pos_whole = row_start_L[row]+1; pos_whole < row_start_L[row+1]; ++pos_whole){
                r0[0] += x[1](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
                r0[1] -= x[0](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
	    }
	    calc_inv(D.diag_element(row), L.diag_element(row), r0);
	    x[0](row) = r0[0]; x[1](row) = r0[1];
        }
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_GaussSeidel(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, const int &step, unsigned int ind_project) const
{ // smooth DoFs on coarse layer first, then the ones only on fine layer
    const SparsityPattern& sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    for (unsigned int n_ite = 0; n_ite < step; ++n_ite){
        for (unsigned int row = 0; row < sp_matrix.m(); ++row){
	    if (core_at[ind_project][row] == -1) continue;
            valuetype r0 = r(row);
            for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
                r0 -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
            x(row) = r0 / sp_matrix.global_entry(row_start[row]);
        }
        for (unsigned int row = 0; row < sp_matrix.m(); ++row){
	    if (core_at[ind_project][row] != -1) continue;
            valuetype r0 = r(row);
            for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
                r0 -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
            x(row) = r0 / sp_matrix.global_entry(row_start[row]);
        }
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_BackwardGaussSeidel(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, const int &step) const
{
    const SparsityPattern& sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    for (unsigned int n_ite = 0; n_ite < step; ++n_ite){
        for (int row = sp_matrix.m()-1; row >= 0; --row){
            valuetype r0 = r(row);
            for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
                r0 -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
            x(row) = r0 / sp_matrix.global_entry(row_start[row]);
        }
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_BackwardGaussSeidel(const SparseMatrix<valuetype> &D, const SparseMatrix<valuetype> &L, std::vector<Vector<valuetype> > &x, const std::vector<Vector<valuetype> > &r, const int &step) const
{
    int n_dof = D.m();
    const SparsityPattern& sp_pattern_D = D.get_sparsity_pattern();
    const std::size_t *row_start_D = sp_pattern_D.get_rowstart_indices();
    const unsigned int *col_nums_D = sp_pattern_D.get_column_numbers();
    const SparsityPattern& sp_pattern_L = L.get_sparsity_pattern();
    const std::size_t *row_start_L = sp_pattern_L.get_rowstart_indices();
    const unsigned int *col_nums_L = sp_pattern_L.get_column_numbers();
    for (unsigned int n_ite = 0; n_ite < step; ++n_ite){
        for (int row = n_dof-1; row >= 0; --row){
	    std::vector<valuetype> r0{r[0](row), r[1](row)};
            for (unsigned int pos_whole = row_start_D[row]+1; pos_whole < row_start_D[row+1]; ++pos_whole){
                r0[0] -= x[0](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
                r0[1] -= x[1](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
	    }
            for (unsigned int pos_whole = row_start_L[row]+1; pos_whole < row_start_L[row+1]; ++pos_whole){
                r0[0] += x[1](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
                r0[1] -= x[0](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
	    }
	    calc_inv(D.diag_element(row), L.diag_element(row), r0);
	    x[0](row) = r0[0]; x[1](row) = r0[1];
        }
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_BackwardGaussSeidel(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, const int &step, unsigned int ind_project) const
{ // smooth DoFs only on fine layer first, then the ones on coarse layer
    const SparsityPattern& sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    for (unsigned int n_ite = 0; n_ite < step; ++n_ite){
        for (unsigned int row = 0; row < sp_matrix.m(); ++row){
	    if (core_at[ind_project][row] != -1) continue;
            valuetype r0 = r(row);
            for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
                r0 -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
            x(row) = r0 / sp_matrix.global_entry(row_start[row]);
        }
        for (unsigned int row = 0; row < sp_matrix.m(); ++row){
	    if (core_at[ind_project][row] == -1) continue;
            valuetype r0 = r(row);
            for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
                r0 -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
            x(row) = r0 / sp_matrix.global_entry(row_start[row]);
        }
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_SymmetricGaussSeidel(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, const int &step) const
{
    const SparsityPattern& sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    for (unsigned int n_ite = 0; n_ite < step; ++n_ite){
        for (unsigned int row = 0; row < sp_matrix.m(); ++row){
            valuetype r0 = r(row);
            for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
                r0 -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
            x(row) = r0 / sp_matrix.global_entry(row_start[row]);
        }
	for (int row = sp_matrix.m()-1; row >= 0; --row){
	    valuetype r0 = r(row);
            for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
                r0 -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
            x(row) = r0 / sp_matrix.global_entry(row_start[row]);
	}
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_SymmetricGaussSeidel(const SparseMatrix<valuetype> &D, const SparseMatrix<valuetype> &L, std::vector<Vector<valuetype> > &x, const std::vector<Vector<valuetype> > &r, const int &step) const
{
    int n_dof = D.m();
    const SparsityPattern& sp_pattern_D = D.get_sparsity_pattern();
    const std::size_t *row_start_D = sp_pattern_D.get_rowstart_indices();
    const unsigned int *col_nums_D = sp_pattern_D.get_column_numbers();
    const SparsityPattern& sp_pattern_L = L.get_sparsity_pattern();
    const std::size_t *row_start_L = sp_pattern_L.get_rowstart_indices();
    const unsigned int *col_nums_L = sp_pattern_L.get_column_numbers();
    for (unsigned int n_ite = 0; n_ite < step; ++n_ite){
        for (unsigned int row = 0; row < n_dof; ++row){
	    std::vector<valuetype> r0{r[0](row), r[1](row)};
            for (unsigned int pos_whole = row_start_D[row]+1; pos_whole < row_start_D[row+1]; ++pos_whole){
                r0[0] -= x[0](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
                r0[1] -= x[1](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
	    }
            for (unsigned int pos_whole = row_start_L[row]+1; pos_whole < row_start_L[row+1]; ++pos_whole){
                r0[0] += x[1](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
                r0[1] -= x[0](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
	    }
	    calc_inv(D.diag_element(row), L.diag_element(row), r0);
	    x[0](row) = r0[0]; x[1](row) = r0[1];
        }
        for (int row = n_dof-1; row >= 0; --row){
	    std::vector<valuetype> r0{r[0](row), r[1](row)};
            for (unsigned int pos_whole = row_start_D[row]+1; pos_whole < row_start_D[row+1]; ++pos_whole){
                r0[0] -= x[0](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
                r0[1] -= x[1](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
	    }
            for (unsigned int pos_whole = row_start_L[row]+1; pos_whole < row_start_L[row+1]; ++pos_whole){
                r0[0] += x[1](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
                r0[1] -= x[0](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
	    }
	    calc_inv(D.diag_element(row), L.diag_element(row), r0);
	    x[0](row) = r0[0]; x[1](row) = r0[1];
	}
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_PMG(Vector<valuetype> &x, const Vector<valuetype> &r,
			      unsigned int ite_step, unsigned int smooth_step) const
{
    Vector<valuetype> r1(r), v0, v1, v2(r.size());
    std::vector<Vector<valuetype>*> projected_r(n_project + 1);
    std::vector<Vector<valuetype>*> projected_x(n_project + 1);
    projected_x[0] = &x;
    projected_r[0] = &r1;
    for (unsigned int i = 1; i <= n_project; ++i) {
        projected_r[i] = new Vector<valuetype>(projected_matrix[i]->m());
        projected_x[i] = new Vector<valuetype>(projected_matrix[i]->m());
    }

    
    for (unsigned int n_ite = 0; n_ite < ite_step; ++n_ite){
        for (unsigned int i = 0; i < n_project; ++i){
	    // smooth
	    // ite_Jacobi(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], smooth_step);
            ite_GaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
            // ite_SymmetricGaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
	    // ite_SSOR(*projected_matrix[i], *projected_x[i], *projected_r[i], 1.0, smooth_step);
            // ite_GaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step, i);
	    // calculate residual
            v0.reinit(projected_x[i]->size());
            projected_matrix[i]->vmult(v0, *projected_x[i]);
            v1 = *projected_r[i];
            v1 -= v0;
            // restrict 
            restrict_matrix[i]->vmult(*projected_r[i+1], v1);
            (*projected_x[i+1]) = 0;
        }

	// solve on the most coarse grid, by Gauss Seidel iterator
	// ite_Jacobi(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], smooth_step);
        // ite_GaussSeidel(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], smooth_step);
	ite_SymmetricGaussSeidel(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], (smooth_step+1)/2);
	// solver.solve(*projected_x[n_project], *projected_r[n_project]);
        // ite_GaussSeidel(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], smooth_step/2);
        // ite_BackwardGaussSeidel(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], smooth_step/2);
	// ite_SSOR(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], 1.0, smooth_step);
	
        for (int i = n_project-1; i >= 0; --i){
	    // interpolate
            v0.reinit(projected_x[i]->size());
            interpolate_matrix[i]->vmult(v0, *projected_x[i+1]);
	    // add contribution from coarse grid
            (*projected_x[i]) += v0;
	    // smooth
	    // ite_Jacobi(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], smooth_step);
            // ite_GaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
            ite_BackwardGaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
            // ite_SymmetricGaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
	    // ite_SSOR(*projected_matrix[i], *projected_x[i], *projected_r[i], 1.0, smooth_step);
            // ite_BackwardGaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step, i);
        }
    }

    
    // delete intermedia variables
    for (unsigned int i = 1; i <= n_project; ++i){
        delete projected_r[i];
        delete projected_x[i];
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_PMG(std::vector<Vector<valuetype> > &x, const std::vector<Vector<valuetype> > &r,
			      unsigned int ite_step, unsigned int smooth_step) const
{
    int n_dof = x[0].size();
    Vector<valuetype> v0, v1;
    std::vector<std::vector<Vector<valuetype> > > projected_r(n_project+1, std::vector<Vector<valuetype> > (2));
    std::vector<std::vector<Vector<valuetype> > > projected_x(n_project+1, std::vector<Vector<valuetype> > (2));
    projected_x[0][0] = x[0]; projected_x[0][1] = x[1];
    projected_r[0][0] = r[0]; projected_r[0][1] = r[1];
    for (unsigned int i = 1; i <= n_project; ++i)
	for (unsigned int idx_ri = 0; idx_ri <= 1; ++idx_ri){
	    projected_r[i][idx_ri].reinit(projected_matrix_D[i]->m());
	    projected_x[i][idx_ri].reinit(projected_matrix_D[i]->m());
	}

    for (unsigned int n_ite = 0; n_ite < ite_step; ++n_ite){
        for (unsigned int i = 0; i < n_project; ++i){
	    // smooth
	    ite_GaussSeidel(*projected_matrix_D[i], *projected_matrix_L[i], projected_x[i], projected_r[i], smooth_step);
	    // calculate residual, restrict
	    // real part
            v0.reinit(projected_x[i][0].size());
            v1 = projected_r[i][0];
            projected_matrix_D[i]->vmult(v0, projected_x[i][0]); v1 -= v0;
	    projected_matrix_L[i]->vmult(v0, projected_x[i][1]); v1 += v0;
            restrict_matrix[i]->vmult(projected_r[i+1][0], v1);
            projected_x[i+1][0] = 0;
	    // imaginary part
            v1 = projected_r[i][1];
            projected_matrix_D[i]->vmult(v0, projected_x[i][1]); v1 -= v0;
	    projected_matrix_L[i]->vmult(v0, projected_x[i][0]); v1 -= v0;
            restrict_matrix[i]->vmult(projected_r[i+1][1], v1);
            projected_x[i+1][1] = 0;
        }
	// solve on the most coarse grid, by Gauss Seidel iterator
	ite_SymmetricGaussSeidel(*projected_matrix_D[n_project], *projected_matrix_L[n_project], projected_x[n_project], projected_r[n_project], (smooth_step+1)/2);
        for (int i = n_project-1; i >= 0; --i){
	    // interpolate, add contribution to fine grid
            v0.reinit(projected_x[i][0].size());
            interpolate_matrix[i]->vmult(v0, projected_x[i+1][0]); projected_x[i][0] += v0;
            interpolate_matrix[i]->vmult(v0, projected_x[i+1][1]); projected_x[i][1] += v0;
	    // smooth
            ite_BackwardGaussSeidel(*projected_matrix_D[i], *projected_matrix_L[i], projected_x[i], projected_r[i], smooth_step);
        }
    }
    x[0] = projected_x[0][0]; x[1] = projected_x[0][1];
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_Jacobi(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, unsigned int ite_step) const
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    Vector<valuetype> y(x.size());
    
    for (int n_ite = 0; n_ite < ite_step; ++n_ite){
    	y = x;
    	for (unsigned int row = 0; row < r.size(); ++row){
    	    valuetype count = (valuetype) r(row);
    	    for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
    		count -= y(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    	    x(row) = count / sp_matrix.global_entry(row_start[row]);
    	}
    }

    // for (unsigned int row = 0; row < r.size(); ++row)
    // 	x(row) = r(row) / sp_matrix.global_entry(row_start[row]);
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_Jacobi(const SparseMatrix<valuetype> &D, const SparseMatrix<valuetype> &L,
					 std::vector<Vector<valuetype> > &x, const std::vector<Vector<valuetype> > &r, unsigned int ite_step) const
{
    const SparsityPattern &sp_pattern_D = D.get_sparsity_pattern();
    const std::size_t *row_start_D = sp_pattern_D.get_rowstart_indices();
    const unsigned int *col_nums_D = sp_pattern_D.get_column_numbers();
    const SparsityPattern &sp_pattern_L = L.get_sparsity_pattern();
    const std::size_t *row_start_L = sp_pattern_L.get_rowstart_indices();
    const unsigned int *col_nums_L = sp_pattern_L.get_column_numbers();
    int n_dof = x[0].size();
    std::vector<Vector<valuetype> > y(2, Vector<valuetype> (n_dof));
    
    for (unsigned int n_ite = 0; n_ite < ite_step; ++n_ite){
    	y[0] = x[0]; y[1] = x[1];
    	for (unsigned int row = 0; row < n_dof; ++row){
    	    std::vector<valuetype> count{r[0](row), r[1](row)};
    	    for (unsigned int pos_whole = row_start_D[row]+1; pos_whole < row_start_D[row+1]; ++pos_whole){
		count[0] -= y[0](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
		count[1] -= y[1](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
	    }
    	    for (unsigned int pos_whole = row_start_L[row]+1; pos_whole < row_start_L[row+1]; ++pos_whole){
		count[0] += y[1](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
		count[1] -= y[0](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
	    }
	    calc_inv(D.diag_element(row), L.diag_element(row), count);
    	    x[0](row) = count[0]; x[1](row) = count[1];
    	}
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_SOR(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, double omega, unsigned int ite_step) const
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    double omega_ = 1 - omega;


    for (int n_ite = 0; n_ite < ite_step; ++n_ite){
	for (unsigned int row = 0; row < r.size(); ++row){
	    valuetype count = (valuetype) r(row);
	    for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
	    x(row) = omega_ * x(row) + omega * count / sp_matrix.global_entry(row_start[row]);
	}
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_SSOR(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, double omega, unsigned int ite_step) const
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    double omega_ = 1 - omega, omega__ = 2 - omega;


    for (int n_ite = 0; n_ite < ite_step; ++n_ite){
    	for (unsigned int row = 0; row < r.size(); ++row){
    	    valuetype count = (valuetype) r(row);
    	    for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
    		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    	    x(row) = omega_ * x(row) + omega * count / sp_matrix.global_entry(row_start[row]);
    	}
    	for (int row = r.size()-1; row >= 0; --row){
    	    valuetype count = (valuetype) r(row);
    	    for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
    		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    	    x(row) = omega_ * x(row) + omega * count / sp_matrix.global_entry(row_start[row]);
    	}
    }


    // for (unsigned int row = 0; row < r.size(); ++row){
    // 	valuetype count = (valuetype) r(row) * omega__;
    // 	for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
    // 	    if (col_nums[pos_whole] < row)
    // 		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    // 	x(row) = omega * count / sp_matrix.global_entry(row_start[row]);
    // }
    // for (int row = r.size()-1; row >= 0; --row){
    // 	valuetype count = x(row) * sp_matrix.global_entry(row_start[row]) / omega;
    // 	for (unsigned int pos_whole = row_start[row+1]; pos_whole < row_start[row+1]; ++pos_whole)
    // 	    if (col_nums[pos_whole] > row)
    // 		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    // 	x(row) = omega * count / sp_matrix.global_entry(row_start[row]);
    // }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_SSOR(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r, double omega,
			       unsigned int max_step, double tol) const
{
    const SparsityPattern &sp_pattern = sp_matrix.get_sparsity_pattern();
    const std::size_t *row_start = sp_pattern.get_rowstart_indices();
    const unsigned int *col_nums = sp_pattern.get_column_numbers();
    double omega_ = 1 - omega;
    Vector<valuetype> res(r), tmp(r.size());
    sp_matrix.vmult(tmp, x);
    res.add(-1, tmp);
    double residual_init = res.l1_norm();
    double residual = residual_init;


    for (int n_ite = 0; n_ite < max_step; ++n_ite){
    	for (unsigned int row = 0; row < r.size(); ++row){
    	    valuetype count = (valuetype) r(row);
    	    for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
    		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    	    x(row) = omega_ * x(row) + omega * count / sp_matrix.global_entry(row_start[row]);
    	}
    	for (int row = r.size()-1; row >= 0; --row){
    	    valuetype count = (valuetype) r(row);
    	    for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
    		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    	    x(row) = omega_ * x(row) + omega * count / sp_matrix.global_entry(row_start[row]);
    	}
	
	res = r;
	sp_matrix.vmult(tmp, x);
	res.add(-1, tmp);
	residual = res.l1_norm();
	if (residual < tol * residual_init) break;
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::ite_SSOR(const SparseMatrix<valuetype> &D, const SparseMatrix<valuetype> &L,
				       std::vector<Vector<valuetype> > &x, const std::vector<Vector<valuetype> > &r, double omega, unsigned int ite_step) const
{
    const SparsityPattern &sp_pattern_D = D.get_sparsity_pattern();
    const std::size_t *row_start_D = sp_pattern_D.get_rowstart_indices();
    const unsigned int *col_nums_D = sp_pattern_D.get_column_numbers();
    const SparsityPattern &sp_pattern_L = L.get_sparsity_pattern();
    const std::size_t *row_start_L = sp_pattern_L.get_rowstart_indices();
    const unsigned int *col_nums_L = sp_pattern_L.get_column_numbers();
    double omega_ = 1 - omega, omega__ = 2 - omega;


    for (unsigned int n_ite = 0; n_ite < ite_step; ++n_ite){
    	for (unsigned int row = 0; row < r.size(); ++row){
	    std::vector<valuetype> count{r[0](row), r[1](row)};
    	    for (unsigned int pos_whole = row_start_D[row]+1; pos_whole < row_start_D[row+1]; ++pos_whole){
    		count[0] -= x[0](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
    		count[1] -= x[1](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
	    }
    	    for (unsigned int pos_whole = row_start_L[row]+1; pos_whole < row_start_L[row+1]; ++pos_whole){
    		count[0] += x[1](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
    		count[1] -= x[0](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
	    }
	    calc_inv(D.diag_element(row), L.diag_element(row), count);
	    x[0](row) = omega_*x[0](row) + omega*count[0];
	    x[1](row) = omega_*x[1](row) + omega*count[1];
    	}
    	for (int row = r.size()-1; row >= 0; --row){
	    std::vector<valuetype> count{r[0](row), r[1](row)};
    	    for (unsigned int pos_whole = row_start_D[row]+1; pos_whole < row_start_D[row+1]; ++pos_whole){
    		count[0] -= x[0](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
    		count[1] -= x[1](col_nums_D[pos_whole]) * D.global_entry(pos_whole);
	    }
    	    for (unsigned int pos_whole = row_start_L[row]+1; pos_whole < row_start_L[row+1]; ++pos_whole){
    		count[0] += x[1](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
    		count[1] -= x[0](col_nums_L[pos_whole]) * L.global_entry(pos_whole);
	    }
	    calc_inv(D.diag_element(row), L.diag_element(row), count);
	    x[0](row) = omega_*x[0](row) + omega*count[0];
	    x[1](row) = omega_*x[1](row) + omega*count[1];
    	}
    }


    // for (unsigned int row = 0; row < r.size(); ++row){
    // 	valuetype count = (valuetype) r(row) * omega__;
    // 	for (unsigned int pos_whole = row_start[row]+1; pos_whole < row_start[row+1]; ++pos_whole)
    // 	    if (col_nums[pos_whole] < row)
    // 		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    // 	x(row) = omega * count / sp_matrix.global_entry(row_start[row]);
    // }
    // for (int row = r.size()-1; row >= 0; --row){
    // 	valuetype count = x(row) * sp_matrix.global_entry(row_start[row]) / omega;
    // 	for (unsigned int pos_whole = row_start[row+1]; pos_whole < row_start[row+1]; ++pos_whole)
    // 	    if (col_nums[pos_whole] > row)
    // 		count -= x(col_nums[pos_whole]) * sp_matrix.global_entry(pos_whole);
    // 	x(row) = omega * count / sp_matrix.global_entry(row_start[row]);
    // }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_PMG(Vector<valuetype> &x, const Vector<valuetype> &r,
				double tol, unsigned int max_step, unsigned int smooth_step) const
{
    Vector<valuetype> r1(r), v0, v1, v2(r.size());
    projected_matrix[0]->vmult(v2, x);
    v2.add(-1.0, r);
    valuetype residual = v2.l1_norm();
    valuetype init_residual = residual;
    std::cerr << "\tinitial residual " << init_residual << '\n';
    unsigned int iteration_step = 0;
    std::vector<Vector<valuetype>*> projected_r(n_project + 1);
    std::vector<Vector<valuetype>*> projected_x(n_project + 1);
    projected_x[0] = &x;
    projected_r[0] = &r1;
    for (unsigned int i = 1; i <= n_project; ++i) {
        projected_r[i] = new Vector<valuetype>(projected_matrix[i]->m());
        projected_x[i] = new Vector<valuetype>(projected_matrix[i]->m());
    }


    std::vector<valuetype> err(max_step+2);
    err[0] = residual;
    
    // while(residual >= tol*init_residual){
    while (residual >= tol){
        for (int i = 0; i < n_project; ++i){
	    // smooth
            ite_GaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
	    // calculate residual
            v0.reinit(projected_x[i]->size());
            projected_matrix[i]->vmult(v0, *projected_x[i]);
            v1 = *projected_r[i];
            v1 -= v0;
            // restrict 
            restrict_matrix[i]->vmult(*projected_r[i+1], v1);
            (*projected_x[i+1]) = 0;
        }

	// solve on the most coarse grid, by Gauss Seidel iterator
        // ite_GaussSeidel(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], smooth_step);
	ite_SymmetricGaussSeidel(*projected_matrix[n_project], *projected_x[n_project], *projected_r[n_project], smooth_step/2);

	// solver.solve(*projected_x[n_project], *projected_r[n_project]);
	
        for (int i = n_project-1; i >= 0; --i){
	    // interpolate
            v0.reinit(projected_x[i]->size());
            interpolate_matrix[i]->vmult(v0, *projected_x[i+1]);
	    // add contribution from coarse grid
            (*projected_x[i]) += v0;
	    // smooth
            // ite_GaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
            ite_BackwardGaussSeidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
        }

	// calculate residual and corresponding l1 norm
        projected_matrix[0]->vmult(v2, x);
        v2.add(-1.0, r);
        residual = v2.l1_norm();

	// update iteration step
        iteration_step++;
	err[iteration_step] = residual;
        if (iteration_step > max_step) break;
    }

    if (flag_output_err){
	std::ofstream output("./PMG_Data.m");
	for (int i = 0; i <= iteration_step; ++i)
	    output << "err(" << i+1 << ") = " << err[i] << ";\n";
	output.close();
    }

    // delete intermedia variables
    for (unsigned int i = 1; i <= n_project; ++i){
        delete projected_r[i];
        delete projected_x[i];
    }


    // output residual and iteration info
    // if (residual < tol*init_residual){
    if (residual < tol){
        std::cerr << "\r\tconverge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
    else {
        std::cerr << "\r\tfailed to converge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_PMG(std::vector<Vector<valuetype> > &x, const std::vector<Vector<valuetype> > &r,
				double tol, unsigned int max_step, unsigned int smooth_step) const
{
    int n_dof = x[0].size();
    Vector<valuetype> v0, v1;
    std::vector<std::vector<Vector<valuetype> > > projected_r(n_project+1, std::vector<Vector<valuetype> > (2));
    std::vector<std::vector<Vector<valuetype> > > projected_x(n_project+1, std::vector<Vector<valuetype> > (2));
    projected_x[0][0] = x[0]; projected_x[0][1] = x[1];
    projected_r[0][0] = r[0]; projected_r[0][1] = r[1];
    for (unsigned int i = 1; i <= n_project; ++i)
    	for (unsigned int idx_ri = 0; idx_ri <= 1; ++idx_ri){
    	    projected_r[i][idx_ri].reinit(projected_matrix_D[i]->m());
    	    projected_x[i][idx_ri].reinit(projected_matrix_D[i]->m());
    	}
    v0.reinit(n_dof);
    v1 = r[0];
    projected_matrix_D[0]->vmult(v0, x[0]); v1 -= v0;
    projected_matrix_L[0]->vmult(v0, x[1]); v1 += v0;
    valuetype residual = v1.l1_norm();
    v1 = r[1];
    projected_matrix_D[0]->vmult(v0, x[1]); v1 -= v0;
    projected_matrix_L[0]->vmult(v0, x[0]); v1 -= v0;
    residual += v1.l1_norm();
    std::cerr << "\tinitial residual " << residual << '\n';
    unsigned int iteration_step = 0;

    std::vector<valuetype> err(max_step+2);
    err[0] = residual;
    
    // while(residual >= tol*init_residual){
    while (residual >= tol){
        for (unsigned int i = 0; i < n_project; ++i){
    	    // smooth
    	    ite_GaussSeidel(*projected_matrix_D[i], *projected_matrix_L[i], projected_x[i], projected_r[i], smooth_step);
    	    // calculate residual, restrict
    	    // real part
            v0.reinit(projected_x[i][0].size());
            v1 = projected_r[i][0];
            projected_matrix_D[i]->vmult(v0, projected_x[i][0]); v1 -= v0;
    	    projected_matrix_L[i]->vmult(v0, projected_x[i][1]); v1 += v0;
            restrict_matrix[i]->vmult(projected_r[i+1][0], v1);
            projected_x[i+1][0] = 0;
    	    // imaginary part
            v1 = projected_r[i][1];
            projected_matrix_D[i]->vmult(v0, projected_x[i][1]); v1 -= v0;
    	    projected_matrix_L[i]->vmult(v0, projected_x[i][0]); v1 -= v0;
            restrict_matrix[i]->vmult(projected_r[i+1][1], v1);
            projected_x[i+1][1] = 0;
        }
    	// solve on the most coarse grid, by Gauss Seidel iterator
    	ite_SymmetricGaussSeidel(*projected_matrix_D[n_project], *projected_matrix_L[n_project], projected_x[n_project], projected_r[n_project], (smooth_step+1)/2);
        for (int i = n_project-1; i >= 0; --i){
    	    // interpolate, add contribution to fine grid
            v0.reinit(projected_x[i][0].size());
            interpolate_matrix[i]->vmult(v0, projected_x[i+1][0]); projected_x[i][0] += v0;
            interpolate_matrix[i]->vmult(v0, projected_x[i+1][1]); projected_x[i][1] += v0;
    	    // smooth
            ite_BackwardGaussSeidel(*projected_matrix_D[i], *projected_matrix_L[i], projected_x[i], projected_r[i], smooth_step);
        }

    	// calculate residual and corresponding l1 norm
    	v1 = projected_r[0][0];
    	projected_matrix_D[0]->vmult(v0, projected_x[0][0]); v1 -= v0;
    	projected_matrix_L[0]->vmult(v0, projected_x[0][1]); v1 += v0;
    	residual = v1.l1_norm();
    	v1 = projected_r[0][1];
    	projected_matrix_D[0]->vmult(v0, projected_x[0][1]); v1 -= v0;
    	projected_matrix_L[0]->vmult(v0, projected_x[0][0]); v1 -= v0;
    	residual += v1.l1_norm();
    	// if (iteration_step % 200 == 0)
    	//     std::cout << "iteration_step = " << iteration_step << ", residual = " << residual << '\n';
	
    	// update iteration step
        iteration_step++;
    	err[iteration_step] = residual;
        if (iteration_step > max_step) break;
    }
    x[0] = projected_x[0][0]; x[1] = projected_x[0][1];

    if (flag_output_err){
    	std::ofstream output("./PMG_Data.m");
    	for (int i = 0; i <= iteration_step; ++i)
    	    output << "err(" << i+1 << ") = " << err[i] << ";\n";
    	output.close();
    }


    // output residual and iteration info
    // if (residual < tol*init_residual){
    if (residual < tol){
        std::cerr << "\r\tconverge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
    else {
        std::cerr << "\r\tfailed to converge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_GaussSeidel(Vector<valuetype> &x, const Vector<valuetype> &r,
					double tol, unsigned int max_step, unsigned int smooth_step) const
{
    Vector<valuetype> r1(r), v0, v1, v2(r.size());
    // calculate initial residual
    projected_matrix[0]->vmult(v2, x);
    v2.add(-1.0, r);
    valuetype residual = v2.l1_norm();
    valuetype init_residual = residual;
    std::cerr << "\tinitial residual " << init_residual << '\n';
    unsigned int iteration_step = 0;
    std::vector<Vector<valuetype>*> projected_r(n_project + 1);
    std::vector<Vector<valuetype>*> projected_x(n_project + 1);
    projected_x[0] = &x;
    projected_r[0] = &r1;

    
    std::vector<valuetype> err(max_step+2);
    err[0] = residual;
    
    // while (residual >= tol*init_residual){
    while (residual >= tol){
	// smooth
	ite_SymmetricGaussSeidel(*projected_matrix[0], *projected_x[0], *projected_r[0], smooth_step/2);

	// calculate l1 norm of residual
	projected_matrix[0]->vmult(v2, x);
	v2.add(-1.0, r);
	residual = v2.l1_norm();
        
	// update iteration step
	iteration_step++;
	err[iteration_step] = residual;
	if (iteration_step > max_step) break;
    }


    if (flag_output_err){
	std::ofstream output("./GS_Data.m");
	for (int i = 0; i <= iteration_step; ++i)
	    output << "err(" << i+1 << ") = " << err[i] << ";\n";
	output.close();
    }
    

    // output residual and iteration info
    // if (residual < tol*init_residual){
    if (residual < tol){
        std::cerr << "\r\tconverge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
    else {
        std::cerr << "\r\tfailed to converge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
}


TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_CG(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r,
			       double tol, unsigned int max_step) const
{
    Vector<valuetype> tmp(r.size()), rk(r.size()), pk(r.size());
    // calculate initial residual
    sp_matrix.vmult(rk, x);
    rk.add(-1.0, r); // r0 = sp_matrix * x0 - rhs
    pk = rk; pk *= -1;
    valuetype residual = rk.l1_norm();
    valuetype init_residual = residual;
    std::cerr << "\tinitial residual " << init_residual << '\n';
    unsigned int iteration_step = 0;


    std::vector<valuetype> err(max_step+2);
    err[0] = residual;


    valuetype rkTrk = (valuetype) 0, pkTApk, alphak, rkTrk_next, betak;
    for (unsigned int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	rkTrk += rk(ind_dof) * rk(ind_dof);
    // while (residual >= tol*init_residual){
    while (residual >= tol){
	// calculate step size
        pkTApk = (valuetype) 0;
	sp_matrix.vmult(tmp, pk);
	for (unsigned int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	    pkTApk += pk(ind_dof) * tmp(ind_dof);
	alphak = rkTrk / pkTApk;

	// update xk+1
	x.add(alphak, pk);

	// update rk+1
	rk.add(alphak, tmp);

	// calcualte step size for updating direction pk
	rkTrk_next = (valuetype) 0;
	for (unsigned int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	    rkTrk_next += rk(ind_dof) *  rk(ind_dof);
	betak = rkTrk_next / rkTrk;

	// update direction pk+1
	pk *= betak;
	pk.add(-1, rk);

	// record residual on this step
	residual = rk.l1_norm();

	// update iteration step and rkTrk
	iteration_step++;
	err[iteration_step] = residual;
	if (iteration_step > max_step) break;
	rkTrk = rkTrk_next;
    }


    if (flag_output_err){
	std::ofstream output("./CG_Data.m");
	for (int i = 0; i <= iteration_step; ++i)
	    output << "err(" << i+1 << ") = " << err[i] << ";\n";
	output.close();
    }
    

    // output residual and iteration info
    // if (residual < tol*init_residual){
    if (residual < tol){
        std::cerr << "\r\tconverge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
    else {
        std::cerr << "\r\tfailed to converge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_PCG(const SparseMatrix<valuetype> &sp_matrix, Vector<valuetype> &x, const Vector<valuetype> &r,
				double tol, unsigned int max_step, unsigned int flag_preconditioner, unsigned int ite_step, double var) const
{/* preconditioned conjugate gradient method
  * flag_preconditioner: 1 - diagonal matrix (Jacobi preconditioner); 2 - sor; 3 - ssor; 4 - PMG
  * var1, var2: extra variable for PCG, ssor: var1 - omega; PMG: var1 - PMG iteration step, var2 - PMG smooth step
  */
    // setup
    double omega = var;
    unsigned int smooth_step;
    if (flag_preconditioner == 4)
	smooth_step = (unsigned int) var;

    
    Vector<valuetype> tmp(r.size()), rk(r.size()), yk(x), pk(r);
    // calculate initial residual
    sp_matrix.vmult(rk, x);
    rk.add(-1.0, r); // r0 = sp_matrix * x0 - rhs
    switch (flag_preconditioner){
    case 0:
	yk = rk;
	break;
    case 1:
	yk = rk;
	for (unsigned int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	    yk(ind_dof) /= sp_matrix.diag_element(ind_dof);
	break;
    case 2:
	if (omega <= 0 || omega >= 2){
	    std::cerr << "invalid sor parameter!\n";
	    return;
	}
	ite_SOR(sp_matrix, yk, rk, omega, ite_step);
	break;
    case 3:
	if (omega <= 0 || omega >= 2){
	    std::cerr << "invalid ssor parameter!\n";
	    return;
	}
	ite_SSOR(sp_matrix, yk, rk, omega, ite_step);
	// ite_SSOR(sp_matrix, yk, rk, omega, 1000, 1.0e-1);
	break;
    case 4:
	ite_PMG(yk, rk, ite_step, smooth_step);
	break;
    default:
	std::cerr << "undefined preconditioner type!\n";
	return;
    }
    pk = yk; pk *= -1;
    valuetype residual = rk.l1_norm();
    valuetype init_residual = residual;
    std::cerr << "\tinitial residual " << init_residual << '\n';
    unsigned int iteration_step = 0;


    std::vector<valuetype> err(max_step+2);
    err[0] = residual;


    valuetype rkTyk = (valuetype) 0, pkTApk, alphak, rkTyk_next, betak;
    for (int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	rkTyk += rk(ind_dof) * yk(ind_dof);
    // while (residual >= tol*init_residual){
    while (residual >= tol){
	// calculate step size
        pkTApk = (valuetype) 0;
	sp_matrix.vmult(tmp, pk);
	for (unsigned int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	    pkTApk += pk(ind_dof) * tmp(ind_dof);
	alphak = rkTyk / pkTApk;

	// update xk
	x.add(alphak, pk);

	// update rk
	rk.add(alphak, tmp);

	// calculate yk+1
	// yk = rk;
	yk *= 0;
	switch (flag_preconditioner){
	case 0:
	    break;
	case 1:
	    yk = rk;
	    for (unsigned int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
		yk(ind_dof) /= sp_matrix.diag_element(ind_dof);
	    break;
	case 2:
	    ite_SOR(sp_matrix, yk, rk, omega, ite_step);
	    break;
	case 3:
	    // yk = rk;
	    // for (int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	    // 	yk(ind_dof) /= sp_matrix.diag_element(ind_dof);
	    ite_SSOR(sp_matrix, yk, rk, omega, ite_step);
	    // ite_SSOR(sp_matrix, yk, rk, omega, 1000, 1.0e-1);
	    break;
	case 4:
	    // yk = rk;
	    // for (int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	    // 	yk(ind_dof) /= sp_matrix.diag_element(ind_dof);
	    ite_PMG(yk, rk, ite_step, smooth_step);
	    // yk = x;
	    // ite_PMG(yk, r, ite_step, smooth_step);
	    // yk.add(-1, x);
	    break;
	}


	// calcualte step size for updating direction pk
	rkTyk_next = (valuetype) 0;
	for (unsigned int ind_dof = 0; ind_dof < r.size(); ++ind_dof)
	    rkTyk_next += rk(ind_dof) *  yk(ind_dof);
	betak = rkTyk_next / rkTyk;

	// update direction pk+1
	pk *= betak;
	pk.add(-1, yk);

	// record residual on this step
	residual = rk.l1_norm();

	// update iteration step and rkTrk
	iteration_step++;
	err[iteration_step] = residual;
	// if (iteration_step % 200 == 0) std::cout << "iteration_step = " << iteration_step << ", err = " << residual << '\n';
	if (iteration_step > max_step) break;
	rkTyk = rkTyk_next;

	if (flag_output_intermedia)
	    std::cout << "iteration_step = " << iteration_step << ", residual = " << residual << '\n';
    }


    if (flag_output_err){
	std::ofstream output("./PCG_Data.m");
	for (int i = 0; i <= iteration_step; ++i)
	    output << "err(" << i+1 << ") = " << err[i] << ";\n";
	output.close();
    }
    

    // output residual and iteration info
    // if (residual < tol*init_residual){
    if (residual < tol){
        std::cerr << "\r\tconverge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
    else {
        std::cerr << "\r\tfailed to converge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_PCG(const SparseMatrix<valuetype> &D, const SparseMatrix<valuetype> &L, 
				std::vector<Vector<valuetype> > &x, const std::vector<Vector<valuetype> > &r,
				double tol, unsigned int max_step) const
{/* preconditioned conjugate gradient method for complex system [D U; -U D], use Jacobi preconditioner
  */
    int n_dof = r[0].size();
    Vector<valuetype> tmp(n_dof), tmp1(n_dof);
    std::vector<Vector<valuetype> > rk(2, Vector<valuetype> (n_dof));
    std::vector<Vector<valuetype> > yk(2, Vector<valuetype> (n_dof));
    std::vector<Vector<valuetype> > pk(2, Vector<valuetype> (n_dof));
    // calculate initial residual, r0 = sp_matrix * x0 - rhs
    D.vmult(rk[0], x[0]); L.vmult(tmp, x[1]); rk[0].add(-1.0, tmp); rk[0].add(-1.0, r[0]);
    D.vmult(rk[1], x[1]); L.vmult(tmp, x[0]); rk[1].add( 1.0, tmp); rk[1].add(-1.0, r[1]);
    std::vector<valuetype> cnt(2);
    for (unsigned int idx_dof = 0; idx_dof < n_dof; ++idx_dof){
    	cnt[0] = rk[0](idx_dof); cnt[1] = rk[1](idx_dof);
    	calc_inv(D.diag_element(idx_dof), L.diag_element(idx_dof), cnt);
    	yk[0](idx_dof) = cnt[0]; yk[1](idx_dof) = cnt[1];
    }
    // yk[0] *= 0; yk[1] *= 0;
    // ite_Jacobi(D, L, yk, rk, 1);
    pk[0] = yk[0]; pk[1] = yk[1]; pk[0] *= -1; pk[1] *= -1;
    valuetype residual = rk[0].l1_norm() + rk[1].l1_norm();
    valuetype init_residual = residual;
    std::cerr << "\tinitial residual " << init_residual << '\n';
    unsigned int iteration_step = 0;


    std::vector<valuetype> err(max_step+2);
    err[0] = residual;


    valuetype rkTyk = (valuetype) 0, pkTApk, alphak, rkTyk_next, betak;
    for (unsigned int idx_dof = 0; idx_dof < n_dof; ++idx_dof)
	rkTyk += rk[0](idx_dof)*yk[0](idx_dof) + rk[1](idx_dof)*yk[1](idx_dof);
    // while (residual >= tol*init_residual){
    while (residual >= tol){
	// calculate step size
        pkTApk = (valuetype) 0;
	D.vmult(tmp, pk[0]); D.vmult(tmp1, pk[1]);
	for (unsigned int idx_dof = 0; idx_dof < n_dof; ++idx_dof)
	    pkTApk += pk[0](idx_dof)*tmp(idx_dof) + pk[1](idx_dof)*tmp1(idx_dof);
	alphak = rkTyk / pkTApk;

	// update xk
	x[0].add(alphak, pk[0]); x[1].add(alphak, pk[1]);

	// update rk
	rk[0].add(alphak, tmp);  L.vmult(tmp,  pk[1]); rk[0].add(-alphak, tmp);
	rk[1].add(alphak, tmp1); L.vmult(tmp1, pk[0]); rk[1].add( alphak, tmp1);

	// calculate yk+1
	// for (unsigned int idx_dof = 0; idx_dof < n_dof; ++idx_dof){
	//     cnt[0] = rk[0](idx_dof); cnt[1] = rk[1](idx_dof);
	//     calc_inv(D.diag_element(idx_dof), L.diag_element(idx_dof), cnt);
	//     yk[0](idx_dof) = cnt[0]; yk[1](idx_dof) = cnt[1];
	// }
	yk[0] = 0; yk[1] = 0;
	// yk[0] = rk[0]; yk[1] = rk[1];
	// ite_Jacobi(D, L, yk, rk, 1);
	// ite_SSOR(D, L, yk, rk, 1.2, 1);
	ite_PMG(yk, rk, 1, 1);


	// calcualte step size for updating direction pk
	rkTyk_next = (valuetype) 0;
	for (unsigned int idx_dof = 0; idx_dof < n_dof; ++idx_dof)
	    rkTyk_next += rk[0](idx_dof)*yk[0](idx_dof) + rk[1](idx_dof)*yk[1](idx_dof);
	betak = rkTyk_next / rkTyk;

	// update direction pk+1
	pk[0] *= betak;	pk[0].add(-1.0, yk[0]);
	pk[1] *= betak;	pk[1].add(-1.0, yk[1]);

	// record residual on this step
	residual = rk[0].l1_norm() + rk[1].l1_norm();

	// update iteration step and rkTrk
	iteration_step++;
	err[iteration_step] = residual;
	if (iteration_step > max_step) break;
	rkTyk = rkTyk_next;

	// if (flag_output_intermedia)
	if (iteration_step % 100 == 0)
	    std::cout << "iteration_step = " << iteration_step << ", residual = " << residual << '\n';
    }


    if (flag_output_err){
	std::ofstream output("./PCG_Data.m");
	for (int i = 0; i <= iteration_step; ++i)
	    output << "err(" << i+1 << ") = " << err[i] << ";\n";
	output.close();
    }
    

    // output residual and iteration info
    // if (residual < tol*init_residual){
    if (residual < tol){
        std::cerr << "\r\tconverge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
    else {
        std::cerr << "\r\tfailed to converge with residual " << residual
                  << " at step " << iteration_step << "." << std::endl;
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_InversePower(Vector<valuetype> &x, const SparseMatrix<valuetype> &A, const SparseMatrix<valuetype> &B,
					 bool Flag_Output_Intermediate, double tol, unsigned int max_step) const
{// Ax = \lambda_min Bx, x^TBx = I, break criterion: e = x^TAx, |e_last-e_now| < tol
    Vector<valuetype> rhs(x.size()), tmp(x.size());
    valuetype e_last, e_now, cnt;
    B.vmult(rhs, x); // preparation
    for (unsigned int n_ite = 0; n_ite < max_step; ++n_ite){
	// iterate
	solve_PCG(A, x, rhs, tol*1.0e-2, max_step*100);
	// normalize
	cnt = 0;
	B.vmult(rhs, x);
	for (unsigned int ind_dof = 0; ind_dof < x.size(); ++ind_dof)
	    cnt += rhs(ind_dof) * x(ind_dof);
	x /= sqrt(cnt); rhs /= sqrt(cnt);
	// calculate e_now
	e_now = 0;
	A.vmult(tmp, x);
	for (unsigned int ind_dof = 0; ind_dof < x.size(); ++ind_dof)
	    e_now += tmp(ind_dof) * x(ind_dof);
	// determine whether break, update e_last
	if (n_ite > 0 && fabs(e_now-e_last) < tol) break;
	if (Flag_Output_Intermediate && n_ite > 0)
	    std::cout << "\tn_ite = " << n_ite << ", e_now = " << e_now << ", dif = " << fabs(e_now-e_last) << '\n';
	e_last = e_now;
    }
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_LOBPCG(TSEM<valuetype> &tsem, Vector<valuetype> &x, const SparseMatrix<valuetype> &A, const SparseMatrix<valuetype> &B,
				   bool Flag_Output_Intermediate, valuetype tol, unsigned int max_step) const
{// Ax = \lambda_min Bx, x^TBx = I, break criterion: e = x^TAx, |e_last-e_now| < tol
     // lapack parameters
    valuetype a[9], b[9];
    valuetype vl, vu;
    int m;
    valuetype w[3];
    int ldz = 3;
    valuetype z[3];
    int ifail[3], info;

    
    // initialize
    // vector consist of [x, p, w]
    unsigned int n_dof_total = x.size();
    std::vector<Vector<valuetype> >  xpw(3, Vector<valuetype> (n_dof_total));
    std::vector<Vector<valuetype> > Axpw(3, Vector<valuetype> (n_dof_total));
    std::vector<Vector<valuetype> > Bxpw(3, Vector<valuetype> (n_dof_total));
    valuetype eig, eig_last;
    valuetype tmp_n = 0, tmp_d = 0; // n for numerator, d for denominator
    xpw[0] = x;
    A.vmult(Axpw[0], xpw[0]);
    B.vmult(Bxpw[0], xpw[0]);
    for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof){
	tmp_n += xpw[0](ind_dof) * Axpw[0](ind_dof);
	tmp_d += xpw[0](ind_dof) * Bxpw[0](ind_dof);
    }
    eig_last = eig = tmp_n / tmp_d;
    // the first search direction is set to be residual
    xpw[1] = Bxpw[0];    xpw[1] *= eig;    xpw[1].add(-1.0, Axpw[0]);
    A.vmult(Axpw[1], xpw[1]);
    B.vmult(Bxpw[1], xpw[1]);
    // vec = [xpw, Axpw, Bxpw] on last step
    std::vector<Vector<valuetype> > vec_last(9, Vector<valuetype> (n_dof_total));
    // x, Ax, Bx
    vec_last[0] = xpw[0];    vec_last[3] = Axpw[0];    vec_last[6] = Bxpw[0];
    // p, Ap, Bp
    vec_last[1] = xpw[1];    vec_last[4] = Axpw[1];    vec_last[7] = Bxpw[1];
    // TSEMSolver<valuetype> solver;
    Vector<valuetype> resi(n_dof_total);
    for (int n_ite = 0; n_ite < max_step; ++n_ite){
	// get residual
	if (n_ite == 0)
	    resi = xpw[1];
	else{
	    resi = Bxpw[0];
	    resi *= eig;
	    resi.add(-1.0, Axpw[0]);
	}

	// assign sparse matrix for preconditioner
	std::vector<unsigned int> nnz(n_dof_total, 0);
	SparsityPattern spp;
	SparseMatrix<valuetype> spm;
	const SparsityPattern &spp_M = tsem.mass_matrix.get_sparsity_pattern();
	const std::size_t *row_start_M = spp_M.get_rowstart_indices();
	const unsigned int *col_nums_M = spp_M.get_column_numbers();
	const SparsityPattern &spp_S = tsem.stiff_matrix.get_sparsity_pattern();
	const std::size_t *row_start_S = spp_S.get_rowstart_indices();
	const unsigned int *col_nums_S = spp_S.get_column_numbers();
	for (unsigned int row = 0; row < n_dof_total; ++row)
	    nnz[row] += row_start_M[row+1]-row_start_M[row] + row_start_S[row+1]-row_start_S[row];
	spp.reinit(n_dof_total, n_dof_total, nnz);
	for (unsigned int row = 0; row < n_dof_total; ++row){
	    for (unsigned int pos_whole = row_start_M[row]; pos_whole < row_start_M[row+1]; ++pos_whole)
		spp.add(row, col_nums_M[pos_whole]);
	    for (unsigned int pos_whole = row_start_S[row]; pos_whole < row_start_S[row+1]; ++pos_whole)
		spp.add(row, col_nums_S[pos_whole]);
	}
	spp.compress();
	spm.reinit(spp);
	for (unsigned int row = 0; row < n_dof_total; ++row){
	    for (unsigned int pos_whole = row_start_M[row]; pos_whole < row_start_M[row+1]; ++pos_whole)
		spm.add(row,  col_nums_M[pos_whole],  tsem.mass_matrix.global_entry(pos_whole)*-eig);
	    for (unsigned int pos_whole = row_start_S[row]; pos_whole < row_start_S[row+1]; ++pos_whole)
		spm.add(row, col_nums_S[pos_whole], tsem.stiff_matrix.global_entry(pos_whole)*0.5);
	}

	// solve
	tsem.impose_zero_boundary_condition(spm);
	solve_PCG(spm, xpw[2], resi, tol*1.0e-2, max_step*100);


	// assign Axpw, Bxpw for residual
	A.vmult(Axpw[2], xpw[2]);
	B.vmult(Bxpw[2], xpw[2]);

	// copy info about w to vec_last
	vec_last[2] =  xpw[2];
	vec_last[5] = Axpw[2];
	vec_last[8] = Bxpw[2];
    
	
	// assign local matrix, calculate its n_pair eigenvectors for corresponding n_pair smallest eigenvalue
	for (unsigned int row = 0; row < 3; ++row)
	    for (unsigned int col = 0; col < 3; ++col){
		tmp_n = tmp_d = 0;
		for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof){
		    tmp_n += xpw[row](ind_dof) * Axpw[col](ind_dof);
		    tmp_d += xpw[row](ind_dof) * Bxpw[col](ind_dof);
		}
		a[row*3 + col] = tmp_n;
		b[row*3 + col] = tmp_d;
	    }


	// solve local eigenvalue problem
	info = LAPACKE_dsygvx(102, 1, 'V', 'I', 'U', 3, a, 3, b, 3, vl, vu, 1, 1, -1, &m, w, z, ldz, ifail);
	if (info > 0){
	    std::cout << "n_ite_lobpcg = " << n_ite << ", the algorithm failed to compute eigenvalues."
		      << ", info = " << info
		      << "\n";
	    break;
	}


	// update xpw, Axpw and Bxpw
	// x, Ax, Bx
	for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
	    xpw[0](ind_dof) = Axpw[0](ind_dof) = Bxpw[0](ind_dof) = 0;
	for (unsigned int ind_p = 0; ind_p < 3; ++ind_p){ // row index of z
	    xpw[ 0].add(z[ind_p], vec_last[ind_p]); // x, from xpw
	    Axpw[0].add(z[ind_p], vec_last[ind_p + 3]); // Ax, from Axpw
	    Bxpw[0].add(z[ind_p], vec_last[ind_p + 6]); // Bx, from Axpw
	}
	// p, Ap, Bp
	for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
	    xpw[1](ind_dof) = Axpw[1](ind_dof) = Bxpw[1](ind_dof) = 0;
	for (unsigned int ind_p = 1; ind_p < 3; ++ind_p){ // the update of p only relates to p and w
	    xpw[ 1].add(z[ind_p], vec_last[ind_p]); // p, from xpw
	    Axpw[1].add(z[ind_p], vec_last[ind_p + 3]); // Ap, from Axpw
	    Bxpw[1].add(z[ind_p], vec_last[ind_p + 6]); // Bp, from Bxpw
	}

	
	// calculate err, break if dif met criterion
	valuetype err = fabs(w[0] - eig_last);
	if (Flag_Output_Intermediate)
	    std::cout << "\tn_ite_lobpcg = " << n_ite << ", err_ene = " << err << '\n';
	if (err < tol) break;
	
	// update vec_last, eig
	// vec_last
	vec_last[0] =  xpw[0]; vec_last[1] =  xpw[1]; vec_last[3] = Axpw[0];
	vec_last[4] = Axpw[1]; vec_last[6] = Bxpw[0]; vec_last[7] = Bxpw[1];
	// eig
	eig_last = eig = w[0];
    }

    x = xpw[0];
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::solve_LOBPCG(TSEM<valuetype> &tsem, std::vector<Vector<valuetype> > &x, const SparseMatrix<valuetype> &A, const SparseMatrix<valuetype> &B,
				   bool Flag_Output_Intermediate, valuetype tol, unsigned int max_step) const
{// Ax = \lambda_min Bx, x^TBx = I, break criterion: e = x^TAx, |e_last-e_now| < tol
     // lapack parameters
    int n_pair = x.size();
    int n = n_pair * 3;
    valuetype a[n*n], b[n*n];
    valuetype vl, vu;
    int m;
    valuetype w[n];
    int ldz = n;
    valuetype z[ldz*n_pair];
    int ifail[n], info;

    
    // initialize
    // vector consist of [x, p, w]
    unsigned int n_dof_total = x[0].size();
    std::vector<Vector<valuetype> >  xpw(n_pair*3, Vector<valuetype> (n_dof_total));
    std::vector<Vector<valuetype> > Axpw(n_pair*3, Vector<valuetype> (n_dof_total));
    std::vector<Vector<valuetype> > Bxpw(n_pair*3, Vector<valuetype> (n_dof_total));
    valuetype eig[n_pair], eig_last[n_pair];
    valuetype tmp_n = 0, tmp_d = 0; // n for numerator, d for denominator
    for (unsigned int ind_pair = 0; ind_pair < n_pair; ++ind_pair){
	xpw[ind_pair] = x[ind_pair];
	A.vmult(Axpw[ind_pair], xpw[ind_pair]);
	B.vmult(Bxpw[ind_pair], xpw[ind_pair]);
	tmp_n = tmp_d = 0;
	for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof){
	    tmp_n += xpw[ind_pair](ind_dof) * Axpw[ind_pair](ind_dof);
	    tmp_d += xpw[ind_pair](ind_dof) * Bxpw[ind_pair](ind_dof);
	}
	eig_last[ind_pair] = eig[ind_pair] = tmp_n / tmp_d;
	// the first search direction is set to be residual
	xpw[ind_pair + n_pair] = Bxpw[ind_pair];
	xpw[ind_pair + n_pair] *= eig[ind_pair];
	xpw[ind_pair + n_pair].add(-1.0, Axpw[ind_pair]);
	A.vmult(Axpw[ind_pair + n_pair], xpw[ind_pair + n_pair]);
	B.vmult(Bxpw[ind_pair + n_pair], xpw[ind_pair + n_pair]);
    }
    // vec = [xpw, Axpw, Bxpw] on last step
    std::vector<Vector<valuetype> > vec_last(n_pair*9, Vector<valuetype> (n_dof_total));
    for (unsigned int ind_pair = 0; ind_pair < n_pair; ++ind_pair){
	// x, Ax, Bx
	vec_last[ind_pair] = xpw[ind_pair];
	vec_last[ind_pair + n_pair*3] = Axpw[ind_pair];
	vec_last[ind_pair + n_pair*6] = Bxpw[ind_pair];
	// p, Ap, Bp
	vec_last[ind_pair + n_pair] = xpw[ind_pair + n_pair];
	vec_last[ind_pair + n_pair*4] = Axpw[ind_pair + n_pair];
	vec_last[ind_pair + n_pair*7] = Bxpw[ind_pair + n_pair];
    }
    for (int n_ite = 0; n_ite < max_step; ++n_ite){
	// attain w by preconditioner, free particle: (S/2-eig*M) * w = r
#pragma omp parallel for
	for (unsigned int ind_pair = 0; ind_pair < n_pair; ++ind_pair){
	    // get residual
	    Vector<valuetype> resi(n_dof_total);
	    if (n_ite == 0)
		resi = xpw[ind_pair + n_pair];
	    else{
		resi = Bxpw[ind_pair];
		resi *= eig[ind_pair];
		resi.add(-1.0, Axpw[ind_pair]);
	    }

	    // assign sparse matrix for preconditioner
	    std::vector<unsigned int> nnz(n_dof_total, 0);
	    SparsityPattern spp;
	    SparseMatrix<valuetype> spm;
	    const SparsityPattern &spp_M = tsem.mass_matrix.get_sparsity_pattern();
	    const std::size_t *row_start_M = spp_M.get_rowstart_indices();
	    const unsigned int *col_nums_M = spp_M.get_column_numbers();
	    const SparsityPattern &spp_S = tsem.stiff_matrix.get_sparsity_pattern();
	    const std::size_t *row_start_S = spp_S.get_rowstart_indices();
	    const unsigned int *col_nums_S = spp_S.get_column_numbers();
	    for (unsigned int row = 0; row < n_dof_total; ++row)
		nnz[row] += (row_start_M[row+1]-row_start_M[row]) + (row_start_S[row+1]-row_start_S[row]);
	    spp.reinit(n_dof_total, n_dof_total, nnz);
	    for (unsigned int row = 0; row < n_dof_total; ++row){
		for (unsigned int pos_whole = row_start_M[row]; pos_whole < row_start_M[row+1]; ++pos_whole)
		    spp.add(row, col_nums_M[pos_whole]);
		for (unsigned int pos_whole = row_start_S[row]; pos_whole < row_start_S[row+1]; ++pos_whole)
		    spp.add(row, col_nums_S[pos_whole]);
	    }
	    spp.compress();
	    spm.reinit(spp);
	    for (unsigned int row = 0; row < n_dof_total; ++row){
		for (unsigned int pos_whole = row_start_M[row]; pos_whole < row_start_M[row+1]; ++pos_whole)
		    spm.add(row, col_nums_M[pos_whole],  tsem.mass_matrix.global_entry(pos_whole)*-eig[ind_pair]);
		for (unsigned int pos_whole = row_start_S[row]; pos_whole < row_start_S[row+1]; ++pos_whole)
		    spm.add(row, col_nums_S[pos_whole], tsem.stiff_matrix.global_entry(pos_whole)*0.5);
	    }

	    // solve
	    tsem.impose_zero_boundary_condition(spm);
	    tsem.impose_zero_boundary_condition(resi);
	    solve_PCG(spm, xpw[ind_pair + n_pair*2], resi, tol*1.0e-2, max_step*100);

	    // assign Axpw, Bxpw for residual
	    A.vmult(Axpw[ind_pair + n_pair*2], xpw[ind_pair + n_pair*2]);
	    B.vmult(Bxpw[ind_pair + n_pair*2], xpw[ind_pair + n_pair*2]);

	    // copy info about w to vec_last
	    vec_last[ind_pair + n_pair*2] =  xpw[ind_pair + n_pair*2];
	    vec_last[ind_pair + n_pair*5] = Axpw[ind_pair + n_pair*2];
	    vec_last[ind_pair + n_pair*8] = Bxpw[ind_pair + n_pair*2];
	}
    
	
	// assign local matrix, calculate its n_pair eigenvectors for corresponding n_pair smallest eigenvalue
	for (unsigned int row = 0; row < n_pair*3; ++row)
	    for (unsigned int col = 0; col < n_pair*3; ++col){
		tmp_n = tmp_d = 0;
		for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof){
		    tmp_n += xpw[row](ind_dof) * Axpw[col](ind_dof);
		    tmp_d += xpw[row](ind_dof) * Bxpw[col](ind_dof);
		}
		a[row*n_pair*3 + col] = tmp_n;
		b[row*n_pair*3 + col] = tmp_d;
	    }


	// solve local eigenvalue problem
	info = LAPACKE_dsygvx(102, 1, 'V', 'I', 'U', 3, a, 3, b, 3, vl, vu, 1, 1, -1, &m, w, z, ldz, ifail);
	if (info > 0){
	    std::cout << "n_ite_lobpcg = " << n_ite << ", the algorithm failed to compute eigenvalues."
		      << ", info = " << info
		      << "\n";
	    break;
	}


	// update xpw, Axpw and Bxpw
	for (unsigned int ind_pair = 0; ind_pair < n_pair; ++ind_pair){ // column index of z
	    // x, Ax, Bx
	    for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
		xpw[ind_pair](ind_dof) = Axpw[ind_pair](ind_dof) = Bxpw[ind_pair](ind_dof) = 0;
	    for (unsigned int ind_p = 0; ind_p < n_pair*3; ++ind_p){ // row index of z
		xpw[ ind_pair].add(z[ind_p + ind_pair*ldz], vec_last[ind_p]); // x, from xpw
		Axpw[ind_pair].add(z[ind_p + ind_pair*ldz], vec_last[ind_p + n_pair*3]); // Ax, from Axpw
		Bxpw[ind_pair].add(z[ind_p + ind_pair*ldz], vec_last[ind_p + n_pair*6]); // Bx, from Axpw
	    }
	    // p, Ap, Bp
	    for (unsigned int ind_dof = 0; ind_dof < n_dof_total; ++ind_dof)
		xpw[ind_pair + n_pair](ind_dof) = Axpw[ind_pair + n_pair](ind_dof) = Bxpw[ind_pair + n_pair](ind_dof) = 0;
	    for (unsigned int ind_p = n_pair; ind_p < n_pair*3; ++ind_p){ // the update of p only relates to p and w
		xpw[ ind_pair + n_pair].add(z[ind_p + ind_pair*ldz], vec_last[ind_p]); // p, from xpw[n_pair*(1:2)]
		Axpw[ind_pair + n_pair].add(z[ind_p + ind_pair*ldz], vec_last[ind_p + n_pair*3]); // Ap, from Axpw[n_pair*(1:2)]
		Bxpw[ind_pair + n_pair].add(z[ind_p + ind_pair*ldz], vec_last[ind_p + n_pair*6]); // Bp, from Bxpw[n_pair*(1:2)]
	    }
	}

	
	// calculate err, break if dif met criterion
	valuetype err = 0, errt;
	for (unsigned int ind_pair = 0; ind_pair < n_pair; ++ind_pair){
	    errt = fabs(w[ind_pair] - eig_last[ind_pair]);
	    if (errt > err) err = errt;
	}
	if (Flag_Output_Intermediate)
	    std::cout << "\tn_ite_lobpcg = " << n_ite << ", err_ene = " << err << '\n';
	if (err < tol) break;
	
	// update vec_last, eig
	for (unsigned int ind_pair = 0; ind_pair < n_pair; ++ind_pair){
	    // vec_last
	    vec_last[ind_pair]            =  xpw[ind_pair];
	    vec_last[ind_pair + n_pair]   =  xpw[ind_pair + n_pair];
	    vec_last[ind_pair + n_pair*3] = Axpw[ind_pair];
	    vec_last[ind_pair + n_pair*4] = Axpw[ind_pair + n_pair];
	    vec_last[ind_pair + n_pair*6] = Bxpw[ind_pair];
	    vec_last[ind_pair + n_pair*7] = Bxpw[ind_pair + n_pair];
	    // eig
	    eig_last[ind_pair] = eig[ind_pair] = w[ind_pair];
	}
    }

    // copy x from xpw
    for (unsigned int ind_pair = 0; ind_pair < n_pair; ++ind_pair)
	x[ind_pair] = xpw[ind_pair];
}

TEMPLATE_TSEMSOLVER
void THIS_TSEMSOLVER::calc_inv(const valuetype d, const valuetype l, std::vector<valuetype> &r) const
{ // calculate r = [d -l; l d]^{-1}r
	if (r.size() != 2) std::cerr << "wrong input for calc_inv!\n";
	valuetype r1 = r[0], r2 = r[1];
	valuetype dn = d*d + l*l;
	// if (dn < 1.0e-8) std::cerr << "near singular local matrix! d = " << d << ", l = " << l << "\n";
	// r[0] = ( d*r1 + l*r2) / dn;
	// r[1] = (-l*r1 + d*r2) / dn;
	valuetype t = l/d;
	r[0] = r1 / (d + l*t) + r2 / (d/t + l);
	r[1] =-r1 / (d/t + l) + r2 / (d + l*t);
}

#undef TEMPLATE_TSEMSOLVER
#undef THIS_TSEMSOLVER
