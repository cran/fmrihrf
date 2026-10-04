#define ARMA_DONT_USE_WRAPPER
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// Single-threaded by design; we compile as C++17 (set in src/Makevars).

using arma::uword;
using namespace Rcpp;


// [[Rcpp::export]]
List neural_input_rcpp(List x, double from, double to, double resolution) {
  int n = (to - from) / resolution;
  NumericVector time(n);
  NumericVector out(n);
  NumericVector ons = x["onsets"];
  NumericVector dur = x["duration"];
  NumericVector amp = x["amplitude"];
  
  for (int i = 0; i < ons.length(); i++) {
    double on = ons[i];
    double d = dur[i];
    int startbin = (int) ((on - from) / resolution) + 1;
    if (d > 0) {
      int endbin = (int) ((on - from) / resolution + d / resolution) + 1;
      for (int j = startbin; j <= endbin; j++) {
        out[j-1] += amp[i];
      }
    } else {
      out[startbin-1] += amp[i];
    }
  }
  
  for (int i = 0; i < n; i++) {
    time[i] = from + (i + 0.5) * resolution;
  }
  
  List result;
  result["time"] = time;
  result["neural_input"] = out;
  return result;
}


/*  Accumulate coef * A_j(x) into a difference array, where A_j is the running
 *  integral of the linear hat (tent) basis function centred on bin j.
 *
 *  A block over [s, e] projected onto the hat basis has bin weights
 *  A_j(e) - A_j(s), which is trapezoid quadrature with the two endpoints
 *  placed at their exact sub-bin positions. Both properties matter: bin-aligned
 *  trapezoid weights quantised the block edges, and exact edges with plain
 *  rectangle weights dropped the quadrature to first order.
 */
static inline void addHatStep(arma::vec& diff, double pos, double coef,
                              double dt, uword nBins)
{
    if (pos <= 0.0) return;                        // nothing accumulated yet
    const double maxPos = (double) (nBins - 1);
    if (pos > maxPos) pos = maxPos;

    const uword a = (uword) std::floor(pos);
    const double f = pos - (double) a;
    const double v = coef * dt;

    // Full weight on every bin strictly left of a.
    if (a > 0) { diff[0] += v; diff[a] -= v; }

    // Partial weight on the two bins bracketing pos.
    const double wa  = v * (0.5 + f - 0.5 * f * f);
    const double wa1 = v * (0.5 * f * f);
    diff[a]     += wa;   diff[a + 1] -= wa;
    diff[a + 1] += wa1;  diff[a + 2] -= wa1;
}

/*──────────────────────────────────────────────────────────────────────
  Fine-grid neural input, built in O(E + N) with a difference array.

  Point events (duration <= 0) contribute a unit-mass impulse, so that
  convolution with the sampled HRF reproduces amp * h(t - onset)
  independently of `dt`.

  Block events (duration > 0) contribute the hat-basis projection of the
  boxcar, so that convolution approximates
  amp * \int_0^duration h(t - onset - u) du and converges as dt -> 0.
  Previously the boxcar had unit height and no dt factor at all, so a block
  response was a bare sample count that grew as ~1/dt: the amplitude of
  every epoch regressor was a function of the `precision` argument.

  Neither onsets nor block edges are snapped to the grid; both are placed at
  their exact sub-bin positions.

  With `summate = false` each block is divided by its total mass
  (`duration`), giving the duration-averaged response rather than the
  accumulated one.
──────────────────────────────────────────────────────────────────────*/
static arma::vec buildImpulseTrain(const arma::vec& on,
                                   const arma::vec& dur,
                                   const arma::vec& amp,
                                   double t0, double t1,
                                   double dt,
                                   bool summate)
{
    if (dt <= 0.0) Rcpp::stop("dt must be > 0 in buildImpulseTrain");

    const uword nBins = (uword) std::floor((t1 - t0) / dt) + 1;
    arma::vec diff(nBins + 2, arma::fill::zeros);      // 2 guard slots
    const double maxPos = (double) (nBins - 1);

    for (uword i = 0; i < on.n_elem; ++i) {
        const double d = dur[i];
        double pos_s = (on[i] - t0) / dt;
        if (pos_s > maxPos) continue;                  // starts past the window

        if (!(d > 0.0)) {
            // Impulse of unit mass, split linearly across the two bracketing
            // bins. Snapping it to the nearest bin instead cost up to dt/2 of
            // onset error, which at the default precision of 0.33 s was ~12%
            // of peak for an onset that fell between bins.
            if (pos_s < 0.0) pos_s = 0.0;
            const uword a = (uword) std::floor(pos_s);
            const double f = pos_s - (double) a;
            const double w0 = amp[i] * (1.0 - f);
            const double w1 = amp[i] * f;
            diff[a]     += w0;  diff[a + 1] -= w0;
            diff[a + 1] += w1;  diff[a + 2] -= w1;
            continue;
        }

        const double pos_e = (on[i] + d - t0) / dt;
        if (pos_e <= 0.0) continue;                    // ends before the window

        // Block weights are the difference of the two running hat integrals,
        // so the block carries mass amp * duration for any sub-bin alignment.
        const double c = amp[i] * (summate ? 1.0 : 1.0 / d);
        addHatStep(diff, pos_e,  c, dt, nBins);
        addHatStep(diff, pos_s, -c, dt, nBins);
    }
    // head(nBins) drops the guard slots.
    return arma::cumsum(diff.head(nBins));
}

// [[Rcpp::export]]
NumericMatrix evaluate_regressor_convolution(NumericVector grid,
                                             NumericVector onsets,
                                             NumericVector durations,
                                             NumericVector amplitudes,
                                             NumericMatrix hrf_values,
                                             double hrf_span,
                                             double start,
                                             double end,
                                             double precision,
                                             bool summate = true) {
  int ngrid = grid.size();
  int nfine = int((end - start) / precision) + 1;
  NumericVector finegrid(nfine);
  for (int i = 0; i < nfine; i++) {
    finegrid[i] = start + i * precision;
  }
  arma::vec neural_input = buildImpulseTrain(as<arma::vec>(onsets),
                                             as<arma::vec>(durations),
                                             as<arma::vec>(amplitudes),
                                             start, end, precision, summate);
  // Convolve neural input with HRF for each basis
  int nbasis = hrf_values.ncol();
  arma::mat conv_result(nfine, nbasis);
  for (int b = 0; b < nbasis; b++) {
    arma::vec hrf_b = hrf_values(_, b);
    arma::vec conv_b = arma::conv(neural_input, hrf_b);
    // Trim the convolution result to match the fine grid size
    conv_result.col(b) = conv_b.subvec(0, nfine - 1);
  }
  // Interpolate conv_result to grid
  NumericMatrix outmat(ngrid, nbasis);
  for (int b = 0; b < nbasis; b++) {
    arma::vec conv_b = conv_result.col(b);
    for (int i = 0; i < ngrid; i++) {
      double t = grid[i];
      if (t <= finegrid[0]) {
        outmat(i, b) = conv_b[0];
      } else if (t >= finegrid[nfine - 1]) {
        outmat(i, b) = conv_b[nfine - 1];
      } else {
        // Linear interpolation
        int idx = int((t - start) / precision);
        if (idx >= nfine - 1) idx = nfine - 2;
        double t1 = finegrid[idx];
        double t2 = finegrid[idx + 1];
        double y1 = conv_b[idx];
        double y2 = conv_b[idx + 1];
        double alpha = (t - t1) / (t2 - t1);
        outmat(i, b) = y1 + alpha * (y2 - y1);
      }
    }
  }
  return outmat;
}


/*──────────────────────────────────────────────────────────────────────
  Unified Rcpp Wrapper for Regressor Evaluation
──────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
SEXP evaluate_regressor_cpp(const arma::vec& grid,
                              const arma::vec& onsets,
                              const arma::vec& durations,
                              const arma::vec& amplitudes,
                              const arma::mat& hrf_matrix, // Use a common name for HRF input
                              double hrf_span,
                              double precision,
                              std::string method = "conv",
                              bool summate = true) {
    // The generated Rcpp export wrapper translates exceptions after C++
    // stack unwinding; do not longjmp through live C++ objects here.
    if (method != "conv") {
        Rcpp::stop("Invalid method for evaluate_regressor_cpp; only 'conv' is supported.");
    }
    // Fine convolution window: back off by the HRF span so events preceding
    // the grid still contribute, and run past the last block so its tail fits.
    double start = grid.min() - hrf_span;
    double max_dur = (durations.n_elem > 0) ? durations.max() : 0.0;
    double last_onset = (onsets.n_elem > 0) ? onsets.max() : grid.max();
    double end = std::max(grid.max(), last_onset + max_dur) + hrf_span;

    NumericMatrix result_mat = evaluate_regressor_convolution(
                                    Rcpp::wrap(grid),
                                    Rcpp::wrap(onsets),
                                    Rcpp::wrap(durations),
                                    Rcpp::wrap(amplitudes),
                                    Rcpp::wrap(hrf_matrix),
                                    hrf_span,
                                    start,
                                    end,
                                    precision,
                                    summate);
    return Rcpp::wrap(result_mat);
}
