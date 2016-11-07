#include <RcppArmadillo.h>

//#include <cmath.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

// [[Rcpp::export]]
//Function to calculate the p3 part of the FRK variance.
//Called within the FRK function.
vec p3Cpp(vec& index, 
			vec& p3, 
			double sigxi,
			mat& SigInv1,
			mat& SigInv2,
			mat& DInv,
			mat& Sp,
			mat& KSSigInv){
				
	//Initialize variables.			
	int n = index.n_rows;
	double ESigInvE;
	int ind = 0;
	
	//Loop through indices.
	for(int i=0; i<n; i++){
		ind = index(i);
        ESigInvE = (DInv(ind,ind) - sum(SigInv1.row(ind) * SigInv2.col(ind)));
		p3(i) +=  2 * sigxi * sum(Sp.row(i) * KSSigInv.col(i)) + pow(sigxi,2) * ESigInvE;
	}
	
	return p3;
				
} //end p3Cpp function.



