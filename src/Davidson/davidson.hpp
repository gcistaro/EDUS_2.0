#include <complex>
#include <Geometry/Vector.hpp>

class Davidson
{
    private:
        std::complex<double>* A_ = nullptr;
        int dim_ = 0;
        int num_solutions_ = 0;
        int max_iter_ = 0;

    public:
        Davidson(std::complex<double>* A__, int dim__, int num_solutions__, int max_iter__) 
         : A_(A__), dim_(dim__), num_solutions_(num_solutions__), max_iter_(max_iter__){};

        void solve()
        {
            /* Create matrix A from the pointera */
            Matrix<std::complex<double>> A(A_, {dim_, dim_});
            /* create the subspace vectors */
            std::vector<Vector<std::complex<double>>> v(max_iter_);
            std::vector<Vector<std::complex<double>>> Av(max_iter_);
            std::vector<Vector<std::complex<double>>> r(max_iter_);
            std::vector<Vector<std::complex<double>>> t_perp(dim_);
            std::vector<Vector<std::complex<double>>> t_parall(dim_);
            std::vector<Vector<std::complex<double>>> t(dim_);
            Matrix<std::complex<double>> M(500, 500);

            int dim_subspace     = 0;   
            int dim_old_subspace = 0;
            int new_vectors      = 2.*num_solutions_; 

            srand((unsigned)time(0)); 
        
            /* initialize t vectors */
            for( int it = 0; it < new_vectors; ++it ) {
                t[it].initialize(dim_);
                for( int ix = 0; ix < dim_; ++ix ) {
                    t[it](ix) = rand();
                }
                t[it] = t[it]/t[it].norm();
            }

            output::print("Starting iterations.");
            /* start davidson iterations */
            for( int iter = 0; iter < max_iter_; ++iter ) {
                output::print("iter = ", iter);

                dim_old_subspace = dim_subspace; 
                dim_subspace     = dim_subspace + new_vectors;
                output::print("dim_old_subspace:", dim_old_subspace);
                output::print("dim_subspace:", dim_subspace);
                output::print("new_vectors:", new_vectors);

                Matrix<std::complex<double>> M_iter(dim_subspace, dim_subspace);

                /* project all new vectors t onto orthogonal complement of V_{k-1} */
                for( int it = 0; it < new_vectors; ++it ) {
                    t_parall[it].initialize(dim_);
                    t_perp[it].initialize(dim_);
                    for( int i = 0; i < dim_old_subspace + it; ++i ) {
                        t_parall[it] = t_parall[it] + ( v[i].dot(t[it]) )*v[i];
                    }
                    t_perp[it] = t[it] - t_parall[it]; 
                    /* Calculate new subspace vectors */

                    v[dim_old_subspace + it] = t_perp[it] / t_perp[it].norm();

                    /* Calculate the matrix-vector product */
                    Av[dim_old_subspace + it] = A*v[dim_old_subspace + it];
                }

                /* Calculate subspace matrix --Warning this can be improved 
                   Here we define the new blocks of the matrix M, here depicted:
                        ( x  x  x  o  o)
                        ( x  x  x  o  o)
                   M =  ( x  x  x  o  o) 
                        ( o  o  o  o  o) 
                        ( o  o  o  o  o) 

                    x-> old step 
                    o-> new step
                */

                for( int irow = dim_old_subspace; irow < dim_subspace; ++irow ) {
                    for( int icol = 0; icol < dim_subspace; ++icol ) {
                        M( irow, icol ) =  v[irow].dot(Av[icol]);
                        M( icol, irow ) = std::conj(M( iter, icol ));
                    }
                }

                /* copy the subspace matrix */
                for( int irow = 0; irow < dim_subspace; ++irow ){
                    for( int icol = 0; icol < dim_subspace ; ++icol ) {
                        M_iter(irow, icol) = M(irow, icol);
                    }
                }
                
                /* Solve subspace eigenvalue problem */
                Matrix<std::complex<double>> s(dim_subspace, dim_subspace);
                mdarray<double,1> theta({dim_subspace});
                M_iter.diagonalize( s, theta );

                /* Calculate Ritz vectors */
                Matrix<std::complex<double>> u(dim_, dim_subspace);
                Matrix<std::complex<double>> Au(dim_, dim_subspace);
                std::vector<Vector<std::complex<double>>> r(dim_subspace);

                /* Compute approximate eigenvectors u_i = V*s_i */
                for( int j = 0; j < dim_; ++j ) {
                    for( int i = 0; i < dim_subspace; ++i ) {
                        r[i].initialize(dim_);
                        for( int k = 0; k < dim_subspace; ++k ) {
                            u(j,i) += v[i](k)*s(k,i);
                            Au(j,i) += Av[i](k)*s(k,i);
                        }
                        /* Compute residual */
                        r[i](j) = Au(j,i) - u(j,i)*theta(i);
                    }
                }

                /* check if eigenvectors are converged */
                new_vectors = 0;
                for( int i = 0; i < dim_subspace; ++i ) {
                    auto norm = r[i].norm();
                    if( norm < 1.e-08 ) {
                        output::print("eigenvector", i, "converged:", norm );
                    }
                    else {
                        output::print("eigenvector", i, "not converged:", norm );
                        t[new_vectors].initialize(dim_);
                        for( int ix = 0; ix < dim_; ++ix ) {
                            t[new_vectors](ix) = r[i](ix)/(theta[i] - A(i,i)); 
                        }
                        new_vectors++;
                    }
                }

                output::print("# of unconverged eigenvectors: ", new_vectors);
            }
            
        }
};