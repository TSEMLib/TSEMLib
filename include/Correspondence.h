/*
 * realization of correspondence
 * declare Correspondence<DIM> correspondence; at the end of this file
 */
#ifndef __Correspondence_
#define __Correspondence_

#include "Multiindex.h"


template <int dim> class Correspondence
{/* correpondence between
  *   multiindex and number of multiindex in lexicographical order
  */

private:
    int M; // truncation order
    int number_index; // number of multiindex $|\alpha|\leq M$
    std::vector<Multiindex<dim> > number_to_index; // number_to_index[i] = i-th multiindex in lexicographical order
    std::vector<std::vector<int> > number_multiindex; /* number_multiindex[i=0:dim-1]: number of multiindices with $(i+1)$ components and summation less than some number
						       *     number_multiindex[i=0:dim-1][j=0:M]: number of multiindices $\beta=(\beta_k)_{k=0}^i\leq j$
						       */
    int calc_binomial_coefficient(int n, int m);
public:
    // Multiindex unitary_multiindex[dim];
    Correspondence();
    Correspondence(int order_truncate);
    void init(int order_truncate); // initialize with truncation order
    int index2number(Multiindex<dim> multiindex); // return the number of given multiindex in lexicographical order
    Multiindex<dim> number2index(int number_index){ // return the (number_index)-th multiindex
        return number_to_index[number_index];}
    int n_index(){ // return number of multiindex $|\alpha|\leq M$, where M is truncation order
        return number_multiindex[dim-1][M];}
    int n_index_begin(int order){ // return begin number of multiindex for $|\alpha|=order$ in lexicographical order
        if (order == 0) return 0;
        else return number_multiindex[dim-1][order-1];}
    int n_index_end(int order){ // return begin number of multiindex for $|\alpha|=order+1$ in lexicographical order, which is also the number of multiindex $|\alpha|\leq order$
        return number_multiindex[dim-1][order];}
};


#define TEMPLATE_CRR template<int dim>
#define THIS_CRR Correspondence<dim>

TEMPLATE_CRR
THIS_CRR::Correspondence(){};

TEMPLATE_CRR
THIS_CRR::Correspondence(int order_truncate){
    init(order_truncate);
}

TEMPLATE_CRR
int THIS_CRR::calc_binomial_coefficient(int n, int m)
{// calculate binomial coefficient $\binom{n}{m}$
    int val = 1;
    if (n > 0 && n >= m){
        if (n-m < m) m = n-m;
        for (int i = 0; i < m; ++i){val *= (n-i); val /= (i+1);}}
    return val;
}

TEMPLATE_CRR
void THIS_CRR::init(int order_truncate)
{// initialize correspondence with given truncation order
    
    M = order_truncate;
    // calculate number_multiindex
    number_multiindex.resize(dim); // number_multiindex[i=0:dim-1]: number of multiindices with $(i+1)$ components and summation less than some number
    for (int i = 0; i < dim; ++i)
        number_multiindex[i].resize(M+1);
    for (int i = 0; i < dim; ++i) // number_multiindex[i=0:dim-1][j=0:M]: number of multiindices $\beta=(\beta_k)_{k=0}^i\leq j$
        for (int j = 0; j <= M; ++j){
            number_multiindex[i][j] = calc_binomial_coefficient(i+j, i);
            if (j > 0) number_multiindex[i][j] += number_multiindex[i][j-1];
        }
    number_index = number_multiindex[dim-1][M];
    // calculate i-th multiindex
    number_to_index.resize(number_index);
    for (int i = 0; i < number_index; ++i){
        int tmp_i = i+1;
        int sum = M;
        for (int j = 0; j < dim; ++j){ // traverse summation of multiindex $|(\alpha_k)_{k=j}^{dim-1}|$
            int k = sum - 1;
            for (; k >= 0; k -= 1)
                if (number_multiindex[dim-1-j][k] < tmp_i)
                    break;
            if (j > 0) number_to_index[i].index[j-1] = sum - (k+1); // now step = 1
            sum = k + 1;
            if (k >= 0) tmp_i -= number_multiindex[dim-1-j][k];
        }
        number_to_index[i].index[dim-1] = sum;
    }
}

TEMPLATE_CRR
int THIS_CRR::index2number(Multiindex<dim> multiindex)
{// calculate the number of given multiindex
    int num = 0; // default for zero multiindex
    for (int i = 0; i < dim; ++i){
        int sum = 0; // summation of last several components of multiindex $|(\alpha_k)_{k=i}^{dim-1}|$
        for (int j = i; j < dim; ++j)
            sum += multiindex.index[j];
        if (sum == 0) break;
        else num += number_multiindex[dim-1-i][sum-1];
    }
    
    return num;
}



#endif
