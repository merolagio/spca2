#ifndef SUPPORT_SHARED_H
#define SUPPORT_SHARED_H

#include <RcppEigen.h>

// Spectra settings shared by all truncated eigendecompositions.
// Matrices smaller than SPECTRA_MIN_DIM are solved with dense Eigen.
static const int SPECTRA_MIN_DIM = 25;
static const double SPECTRA_TOL = 1e-10;
static const int SPECTRA_MAXIT = 1000;

int spectra_ncv(int nev, int n);

Eigen::MatrixXd makeCorComp_int(const Eigen::Ref<const Eigen::MatrixXd>& A,
                                const Eigen::Ref<const Eigen::MatrixXd>& S,
                                int ncomps);

Eigen::MatrixXd makeSubS(const Eigen::Ref<const Eigen::MatrixXd>& S,
                         Eigen::VectorXi& e);

void maked_loopColF(const Eigen::Ref<const Eigen::MatrixXd>& S,
                    Eigen::MatrixXd& M,
                    Eigen::VectorXi& e);

void makeXdC(const Eigen::Ref<const Eigen::MatrixXd>& X,
             const Eigen::Ref<const Eigen::VectorXi>& ind,
             Eigen::MatrixXd& Xd);

Eigen::MatrixXd CorCompC(const Eigen::Ref<const Eigen::MatrixXd>& A,
                         const Eigen::Ref<const Eigen::MatrixXd>& S,
                         int d = 0);
Eigen::MatrixXd cor_int(const Eigen::Ref<const Eigen::MatrixXd>& X,
                        bool center = true,
                        bool scale = true);

Eigen::VectorXd eigvecPMC(const Eigen::MatrixXd& X,
                          double& val,
                          double eps = 1e-5,
                          int maxiter = 100);

Eigen::VectorXd GeigvecPMC(const Eigen::MatrixXd& A,
                           const Eigen::MatrixXd& B,
                           double& val,
                           double eps = 1e-5,
                           int maxiter = 100);

#endif
