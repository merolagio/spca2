#include <RcppEigen.h>
#include <cmath>
#include <stdexcept>
#include <string>

#include "support_tall.h"
#include "support_fat.h"
// Spectra headers are provided by RSpectra (LinkingTo: RSpectra)
#include <Spectra/SymEigsSolver.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;
using namespace std;

// Basic validation shared by the full-eigen and power-method PCA paths.
// More detailed data handling is done by the R wrapper.
static void validate_pcaC_inputs(const Eigen::Ref<const Eigen::MatrixXd>& M,
                                 int ncomps,
                                 bool data_matrix,
                                 bool fat_matrix)
{
  const int n = M.rows();
  const int p = M.cols();

  if (n < 1 || p < 1)
    Rcpp::stop("M must have at least one row and one column");

  if (!M.allFinite())
    Rcpp::stop("M contains non-finite values");

  if (!data_matrix) {
    if (n != p)
      Rcpp::stop("covariance/correlation input must be square");

    const double nrm = M.norm();
    if (nrm > 0.0 && (M - M.transpose()).norm() > 1e-6 * nrm)
      Rcpp::stop("covariance/correlation input must be symmetric");
  }

  if (fat_matrix) {
    if (!data_matrix)
      Rcpp::stop("fat_matrix = TRUE requires a data matrix");

    if (n >= p)
      Rcpp::stop("fat_matrix = TRUE requires n < p");
  }

  const int max_comps = fat_matrix ? n : p;
  if (ncomps < 1 || ncomps > max_comps)
    Rcpp::stop("ncomps is outside the allowed range for this backend");
}

// Exported low-level helper for computing leading eigenpairs by power method.
// For fat_matrix = TRUE, M is a data matrix; otherwise M is a covariance matrix.
// [[Rcpp::export]]
Rcpp::List PMnEigenpairs(const Eigen::Map<Eigen::MatrixXd>& M,
                      int ncomps,
                      bool fat_matrix = false,
                      double epsPM = 1e-5,
                      int maxiterPM = 100)
{
  try {
    if (ncomps < 1)
      Rcpp::stop("ncomps must be positive");
    if (epsPM <= 0.0 || !std::isfinite(epsPM))
      Rcpp::stop("epsPM must be positive and finite");
    if (maxiterPM < 1)
      Rcpp::stop("maxiterPM must be positive");
    if (!M.allFinite())
      Rcpp::stop("M contains non-finite values");

    if (fat_matrix) {
      if (M.rows() >= M.cols())
        Rcpp::stop("fat_matrix = TRUE requires n < p");
      if (ncomps > M.rows())
        Rcpp::stop("ncomps cannot exceed n in the fat backend");
      return PMnEigenpairs_fat(M, ncomps, epsPM, maxiterPM);
    }

    if (M.rows() != M.cols())
      Rcpp::stop("tall PMnEigenpairs expects a square covariance/correlation matrix");
    if (ncomps > M.cols())
      Rcpp::stop("ncomps cannot exceed matrix dimension");
    return PMnEigenpairs_tall(M, ncomps, epsPM, maxiterPM);
  } catch (std::exception& ex) {
    Rcpp::stop(std::string("PMnEigenpairs: ") + ex.what());
  } catch (...) {
    Rcpp::stop("PMnEigenpairs: unknown C++ exception");
  }
}

// Cross-product operator for Spectra, without forming the cross product.
// Tall (fat = false): y = X'(X x) / denom, dimension p.
// Fat  (fat = true):  y = X(X' x) / denom, dimension n.
class CrossProdOp
{
private:
  Eigen::Ref<const Eigen::MatrixXd> X;
  const bool fat;
  const double denom;
  Eigen::VectorXd work;

public:
  CrossProdOp(const Eigen::Ref<const Eigen::MatrixXd>& X_, bool fat_, double denom_) :
    X(X_), fat(fat_), denom(denom_), work(fat_ ? X_.cols() : X_.rows())
  {}

  int rows() const { return fat ? X.rows() : X.cols(); }
  int cols() const { return rows(); }

  void perform_op(const double* x_in, double* y_out)
  {
    const int d = rows();
    Eigen::Map<const Eigen::VectorXd> x(x_in, d);
    Eigen::Map<Eigen::VectorXd> y(y_out, d);
    if (fat) {
      work.noalias() = X.transpose() * x;
      y.noalias() = X * work;
    } else {
      work.noalias() = X * x;
      y.noalias() = X.transpose() * work;
    }
    y /= denom;
  }
};

// Leading nev eigenpairs by Spectra. Returns false if Spectra did not converge.
template <typename OpType>
static bool spectra_leading(OpType& op, int nev, int dim,
                            Eigen::VectorXd& values, Eigen::MatrixXd& vectors)
{
  Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, OpType>
    eigs(&op, nev, spectra_ncv(nev, dim));
  eigs.init();
  eigs.compute(SPECTRA_MAXIT, SPECTRA_TOL, Spectra::LARGEST_ALGE);
  if (eigs.info() != Spectra::SUCCESSFUL)
    return false;
  values = eigs.eigenvalues();
  vectors = eigs.eigenvectors(nev);
  return true;
}

// pcaC
//
// Main C++ PCA backend used by the R wrapper.
// Computes ordinary PCA for tall/covariance inputs and row-space PCA for fat
// data matrices. Data matrices are assumed to be centred (and scaled if
// required) by the R wrapper; the cross products are divided by n - 1, so
// the eigenvalues are those of the covariance (or correlation) matrix.
//
// If PM = true (default in R), only the requested leading eigenpairs are
// computed with Spectra (Lanczos). For data matrices the cross product is
// never formed: Spectra works with the operator X'X or XX'. Dense Eigen is
// used instead when PM = false, when ncomps is at least half the dimension,
// when the dimension is smaller than SPECTRA_MIN_DIM, or (with a warning)
// when Spectra does not converge.
//
// @param M Numeric matrix. If data_matrix = true, M is a data matrix. If
//   data_matrix = false, M is a covariance or correlation matrix.
// @param ncomps Number of principal components/eigenpairs to return.
// @param data_matrix Logical. TRUE when M is a data matrix; FALSE when M is a
//   covariance or correlation matrix.
// @param fat_matrix Logical. TRUE uses the fat matrix backend based on M M';
//   FALSE uses the tall/covariance backend.
// @param PM Logical. TRUE uses Spectra truncated eigendecomposition; FALSE
//   uses full dense eigendecomposition.
// @param epsPM,maxiterPM Ignored, kept for backward compatibility.
//
// @return List with:
//   weights: p x ncomps matrix of unit L2 weight vectors.
//   scores: n x ncomps matrix of PC scores when M is a data matrix, otherwise
//     an empty matrix.
//   vexpPC: variance explained by the returned PCs, scaled by total variance.
//   eigenvalues: eigenvalues returned by the selected backend. With Spectra,
//     only ncomps eigenvalues are returned.
//   totvar: total variance, equal to the trace of the covariance matrix.
//   ncomps: number of returned components.
//   fat_matrix: backend flag used by the computation.
//   PM: TRUE if Spectra was used.
// [[Rcpp::export]]
Rcpp::List pcaC(const Eigen::Map<Eigen::MatrixXd>& M,
                int ncomps,
                bool data_matrix = true,
                bool fat_matrix = false,
                bool PM = true,
                double epsPM = 1e-5,
                int maxiterPM = 100)
{
  (void) epsPM;
  (void) maxiterPM;
  try {
    validate_pcaC_inputs(M, ncomps, data_matrix, fat_matrix);

    const int n = M.rows();
    if (data_matrix && n < 2)
      Rcpp::stop("a data matrix needs at least two rows");
    const double denom = data_matrix ? double(n - 1) : 1.0;

    const int dim = fat_matrix ? n : M.cols();
    bool use_spectra = PM && (dim >= SPECTRA_MIN_DIM) && (2 * ncomps < dim);

    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd vec;   // eigenvectors of the matrix actually decomposed
    const double totvar = data_matrix ? M.squaredNorm() / denom : M.trace();

    if (use_spectra) {
      bool ok = false;
      if (data_matrix) {
        CrossProdOp op(M, fat_matrix, denom);
        ok = spectra_leading(op, ncomps, dim, eigenvalues, vec);
      } else {
        Spectra::DenseSymMatProd<double> op(M);
        ok = spectra_leading(op, ncomps, dim, eigenvalues, vec);
      }
      if (!ok) {
        Rcpp::warning("Spectra did not converge; using dense eigendecomposition");
        use_spectra = false;
      }
    }

    if (!use_spectra) {
      Eigen::MatrixXd K;
      if (!data_matrix)
        K = M;
      else if (fat_matrix)
        K = (M * M.transpose()) / denom;
      else
        K = (M.transpose() * M) / denom;

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
      if (es.info() != Eigen::Success)
        Rcpp::stop("pcaC: dense eigendecomposition failed");
      eigenvalues = Eigen::VectorXd(es.eigenvalues()).reverse();
      vec = Eigen::MatrixXd(es.eigenvectors()).rowwise().reverse().leftCols(ncomps);
    }

    Eigen::MatrixXd weights;
    Eigen::MatrixXd scores(0, 0);

    if (fat_matrix) {
      // Row-space eigenvectors u: weights a = X'u / ||X'u||,
      // scores X a = sqrt((n - 1) * lambda) u.
      weights = M.transpose() * vec;
      scores = vec;
      for (int j = 0; j < ncomps; ++j) {
        const double nrm = weights.col(j).norm();
        if (!std::isfinite(nrm) || nrm <= 0.0)
          Rcpp::stop("pcaC: zero or non-finite weight norm in fat backend");
        weights.col(j) /= nrm;
        if (weights(0, j) < 0.0) {
          weights.col(j) = -weights.col(j);
          scores.col(j) = -scores.col(j);
        }
        scores.col(j) *= std::sqrt(std::max(0.0, eigenvalues(j)) * denom);
      }
    } else {
      weights = vec;
      for (int j = 0; j < ncomps; ++j) {
        if (weights(0, j) < 0.0)
          weights.col(j) = -weights.col(j);
      }
      if (data_matrix)
        scores = M * weights;
    }

    return Rcpp::List::create(
      Rcpp::Named("weights") = weights,
      Rcpp::Named("scores") = scores,
      Rcpp::Named("vexpPC") = eigenvalues.head(ncomps) / totvar,
      Rcpp::Named("eigenvalues") = eigenvalues,
      Rcpp::Named("totvar") = totvar,
      Rcpp::Named("ncomps") = ncomps,
      Rcpp::Named("fat_matrix") = fat_matrix,
      Rcpp::Named("PM") = use_spectra
    );
  } catch (std::exception& ex) {
    Rcpp::stop(std::string("pcaC: ") + ex.what());
  } catch (...) {
    Rcpp::stop("pcaC: unknown C++ exception");
  }
}

// eigvalsC
//
// All eigenvalues of the covariance (or correlation) matrix, without
// eigenvectors. Used by eigvals().
// For a data matrix (assumed centred, and scaled if required, by the R
// wrapper) the eigenvalues are those of X'X / (n - 1), computed from the
// smaller of X'X and XX' (the nonzero eigenvalues are the same).
// For a covariance or correlation matrix they are those of M.
//
// @param M Numeric matrix: data matrix or covariance/correlation matrix.
// @param data_matrix Logical. TRUE when M is a data matrix.
//
// @return List with eigenvalues (decreasing) and totvar (the trace).
// [[Rcpp::export]]
Rcpp::List eigvalsC(const Eigen::Map<Eigen::MatrixXd>& M,
                    bool data_matrix = true)
{
  try {
    const int n = M.rows();
    const int p = M.cols();
    if (n < 1 || p < 1)
      Rcpp::stop("M must have at least one row and one column");
    if (!M.allFinite())
      Rcpp::stop("M contains non-finite values");

    Eigen::MatrixXd K;
    if (data_matrix) {
      if (n < 2)
        Rcpp::stop("a data matrix needs at least two rows");
      const double denom = double(n - 1);
      if (n < p)
        K = (M * M.transpose()) / denom;
      else
        K = (M.transpose() * M) / denom;
    } else {
      if (n != p)
        Rcpp::stop("covariance/correlation input must be square");
      K = M;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K, Eigen::EigenvaluesOnly);
    if (es.info() != Eigen::Success)
      Rcpp::stop("eigvalsC: eigendecomposition failed");

    Eigen::VectorXd ev = Eigen::VectorXd(es.eigenvalues()).reverse();
    return Rcpp::List::create(
      Rcpp::Named("eigenvalues") = ev,
      Rcpp::Named("totvar") = K.trace()
    );
  } catch (std::exception& ex) {
    Rcpp::stop(std::string("eigvalsC: ") + ex.what());
  } catch (...) {
    Rcpp::stop("eigvalsC: unknown C++ exception");
  }
}
