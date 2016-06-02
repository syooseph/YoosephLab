#include "gem.h"

GEM::GEM()
{
}

GEM::GEM( MatrixType *m,
            size_t nr,
            size_t nc,
            size_t cpu,
            size_t num,
            double tol,
            double psu,
            bool eq,
            bool v )
{
    matrix    = m;
    nrow      = nr;
    ncol      = nc;
    ncpu      = cpu;
    max_iter  = num;
    tolerance = tol;
    pseudo    = psu;
    even      = eq;    
    verbose   = v;

	t_exp = t_max = t_con = 0.0;
}

GEM::~GEM()
{
    delete[] mixing;
    delete[] mixold;
    
    //if ( ncpu > 1 ) delete[] LogLiks;
}

void GEM::run()
{
    allocate();

    performEM();
}

void GEM::allocate()
{
	double t0 = mytime();
	
    //-------------------
    // Mixing coefficient
    //-------------------    
    mixing = new double[ncol];    
    mixold = new double[ncol];
    
    //-----------------------------------------------
    // Matrices for joint and posterior probabilities
    //-----------------------------------------------
	Q = MatrixType(*matrix);
	R = MatrixType(*matrix);
    // Q = MatrixType(nrow, ncol);
    // R = MatrixType(nrow, ncol);
    // for ( size_t i = 0; i < nrow; i++ ) {
    //     for ( size_t j = 0; j < ncol; j++ ) {
	// 		double v = (*matrix)(i,j);
    //         if ( v == 0 ) continue;
	// 		Q(i,j) = 0.0;
	// 		R(i,j) = 0.0;
    //         // Q[i].insert( std::pair<uint32_t, double>( it->first, 0.0) );
    //         // R[i].insert( std::pair<uint32_t, double>( it->first, 0.0) );
    //     }
    // }

    // //------------------------------------------------------------
    // // Auxilary array for multi-threaded loglikelihood calculation
    // //------------------------------------------------------------
    // if ( ncpu > 1 )
    //     LogLiks = new double[nrow];

	//empties = std::vector<std::vector<bool>>( ncpu, std::vector<bool>(ncol, true));

	if ( verbose ) std::cerr << "Q & R matrix allocated:" << mytime()-t0 << "\n";
}

void GEM::performEM()
{
    //-------------------------------
    // Initialize mixing coefficients
    //-------------------------------
    initMixingCoefficient();
    //addPseudoCount();
    
    converge = false;
    likelihood = likelidiff = exp(10); //Initial likeliehood
    // error = exp(10);  // sum of errors of mixing coefficients
	merr = exp(10);  // sum of errors of mixing coefficients
    rmse = exp(10);

    iter = 0;
    while ( iter < max_iter && !converge ) {
        double t0 = mytime();
        iter++;

        doExpectation();
        doMaximization();
        checkConvergency();
		
		//---------------------------
        // Update mixing coefficients
        //---------------------------
        std::copy(mixing, mixing+ncol, mixold);        

        if ( verbose )
			showIterationSummary(t0);
    }
}

void GEM::initMixingCoefficient()
{
    //---------------------------------
    // Simple equal weight coefficients
    //---------------------------------
    if ( even ) {
        for ( size_t c = 0; c < ncol; c++ )
            mixing[c] = 1.0/ncol;
        return;
    }

    //--------------------------------------------------
    // Intitalize coefficient based on number column sum
    //--------------------------------------------------
    double total = 0;
    size_t c,r;
    std::vector<double> colsum(ncol,0.0);            
#pragma omp parallel for schedule(dynamic, 1) if (ncpu) private(c,r) num_threads(ncpu)                    
    
    for ( c = 0; c < ncol; c++ ) {
        size_t sum = 0;
        for ( r = 0; r < nrow; r++ ) {
			double v = (*matrix)(r,c);
            if ( v == 0 ) continue;
            sum += v;
        }
        colsum[c] = sum;
#pragma omp critical
        total += sum;
    }

    for ( size_t c = 0; c < ncol; c++ ) 
        mixing[c] = colsum[c]/total;
}

// void GEM::addPseudoCount()
// {
//     if ( pseudo == 0 ) return;

//     size_t i,j;
// #pragma omp parallel for if (ncpu) private(i,j) num_threads(ncpu)
//     for ( i = 0; i < nrow; i++ ) 
//         for ( j = 0; j < ncol; j++ )
//             matrix[i][j] += pseudo;
// }


// void GEM::getEmptyColumns( size_t w, size_t r )
// {
// 	for ( size_t j = 0; j < ncol; j++ ) {
// 		double v = (*matrix)(r,j);
// 		empties[w][j] = (v==0);
// 	}
// }

void GEM::doExpectation()
{
	double t0 = mytime();
	
    size_t i,j;

#pragma omp parallel for schedule(static, 1) if (ncpu) private(i,j) num_threads(ncpu)        
    for ( i = 0; i < nrow; i++ ) {

// 		size_t w =  omp_get_thread_num();

// // 		if ( verbose )
// // #pragma omp critical
// // 			std::cerr << "row:" << i << "\twhich:" << w << "\n";
// 		getEmptyColumns(w, i);
		
		double rsum = 0;
        for ( j = 0; j < ncol; j++ ) {
			//if ( empties[w][j] ) continue;
			double v = (*matrix)(i,j);
            if ( v == 0 ) continue;
			double n = mixing[j]*v;
            //Q(i,j) = mixing[j]*v;
			Q(i,j) = n;
			//Q(i,j) = mixing[j]*(*matrix)(i,j);
			rsum += n;//Q(i,j);
        }

        // double rsum = 0;
        // for ( j = 0; j < ncol; j++ ) {
		// 	double v = (*matrix)(i,j);
		// 	if ( v == 0 ) continue;
        //     rsum += Q(i,j);
        // }

        for ( size_t j = 0; j < ncol; j++ ) {
			//if ( empties[w][j] ) continue;
			double v = (*matrix)(i,j);
			if ( v == 0 ) continue;
            R(i,j) = Q(i,j)/rsum;
        }
    }

    // if ( verbose ) {
    //     std::cerr << "Q matrix\n";
    //     for ( i = 0; i < nrow; i++ ) {
    //         std::cerr << i << "\t";
    //         for ( j = 0; j < ncol; j++ ) {
    //             auto it = Q[i].find(j);
    //             if ( it == Q[i].end() )
    //                 std::cerr << 0 << " ";
    //             else std::cerr << it->second << " ";
    //         }
    //         std::cerr << "\n";
    //     }
    //     std::cerr << "R matrix\n";
    //     for ( i = 0; i < nrow; i++ ) {
    //         std::cerr << i << "\t";
    //         for ( j = 0; j < ncol; j++ ) {
    //             auto it = R[i].find(j);
    //             if ( it == R[i].end() )
    //                 std::cerr << 0 << " ";
    //             else std::cerr << it->second << " ";
    //         }
    //         std::cerr << "\n";
    //     }
    // }

	t_exp += (mytime()-t0);
}

void GEM::doMaximization()
{
	double t0 = mytime();
	// Need to test whether multi-threading is faster here
	// column based operation
	size_t i,j;
#pragma omp parallel for if (ncpu) private(i,j) num_threads(ncpu)                
	for ( j = 0; j < ncol; j++ ) {
		double csum = 0.0;
		for ( i = 0; i < nrow; i++ ) {
			double v = (*matrix)(i,j);
			if ( v == 0 ) continue;
			csum += R(i,j);
		}
		mixing[j] = csum/nrow;
	}
	// size_t i,j;
    // std::vector<double> colsum(ncol,0.0);
    // for ( i = 0; i < nrow; i++ ) 
    //     for ( j = 0; j < ncol; j++ )
    //         colsum[j] += R(i,j);
	
    // for ( j = 0; j < ncol; j++ )
    //     mixing[j] = colsum[j]/nrow;

	t_max += (mytime()-t0);
}

void GEM::checkConvergency()
{
	double t0 = mytime();
	
    double l = computeLikelihood();
    likelidiff = std::abs(likelihood-l);
    //if ( iter >1 && likelidiff < tolerance ) converge = true;
    likelihood = l;
	
	merr = computeMeanAbsError();
    rmse = computeRMSE();
	if ( iter >1 && merr < tolerance ) converge = true;

	t_con += (mytime()-t0);
}

double GEM::computeLikelihood()
{
    double loglik = 0;
    size_t i,j;

    double* LogLiks = new double[nrow];    
    std::fill(LogLiks, LogLiks+nrow, 0.0);
#pragma omp parallel for private(i,j) num_threads(ncpu)
    for ( i = 0; i < nrow; i++ ) {
        for ( j = 0; j < ncol; j++ ) {
			double q = Q(i,j);
            if ( q == 0 ) continue;
            //if ( isinf( log(q) ) ) continue;
			if ( ! std::isfinite( log(q) ) ) continue;
            LogLiks[i] += ( R(i,j) * log(q) );
        }
    }
    
    for ( i = 0; i < nrow; i++ )
        loglik += LogLiks[i];

    delete[] LogLiks;
    return loglik;
}

double GEM::computeMeanAbsError()
{
    if ( iter == 1 ) return exp(10);
    double e = 0;
    for ( size_t c = 0; c < ncol; c++ ) 
        e += std::abs(mixold[c]-mixing[c]);
    e /= ncol;
    return e;
}

double GEM::computeRMSE()
{
    if ( iter == 1 ) return exp(10);
    double e = 0;
    for ( size_t c = 0; c < ncol; c++ ) 
        e += std::abs(mixold[c]-mixing[c]);
    return e;
}
void GEM::showIterationSummary( double t0 )
{
	std::cerr << "Iteration:" << iter << "\telapsed:" << mytime()-t0 << "\n";
	std::cerr << "Likelihood change:" << std::setprecision(10) << likelidiff << "\n";
	std::cerr << "Absolute error:" << merr << "\n";
	std::cerr << "mixing coefficient:" << "\n";
	for ( size_t j = 0; j < ncol; j++ )
		std::cerr << mixing[j] << " ";
	std::cerr << "\n";
	printf("Timing (Expecation:%.2f, Maximization:%.2f, Convergence:%.2f)\n", t_exp, t_max, t_con);
}

// void GEM::computeAbundance( double *abundance, uint32_t *lengths )
// {
//     size_t j;
//     double total_ratio = 0;
// #pragma omp parallel for private(j) num_threads(ncpu)    
//     for ( j = 0; j < ncol; j++ ) {
//         double ratio = mixing[j]/lengths[j];
// #pragma omp critical
//         total_ratio += ratio;
//     }

// #pragma omp parallel for private(j) num_threads(ncpu)        
//     for ( j = 0; j < ncol; j++ )
//         abundance[j] = mixing[j]/(lengths[j]*total_ratio);
// }

// double GEM::at(MatrixType &m, size_t i, size_t j)
// {
//     auto it = m[i].find(j);
//     if ( it == m[i].end() )
//         return 0;
//     return it->second;
// }
