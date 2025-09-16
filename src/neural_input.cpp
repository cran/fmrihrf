#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#define ARMA_DONT_USE_WRAPPER
#include <RcppArmadillo.h>
// #include <RcppParallel.h> // REMOVED
#include <unordered_map>
#include <mutex> // For HRFworker error_mutex, can be removed if HRFworker is fully removed
#include <sstream> // For HRFworker error messages, can be removed if HRFworker is fully removed

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]] // REMOVED
// We compile as C++17 (set in src/Makevars)

using arma::uword;
// using namespace RcppParallel; // REMOVED
using namespace Rcpp;


/*──── helper: next power of two ────────────────────────────────────*/
/* Armadillo ≥ 13.0 supposedly exposes arma::next_power_of_two(), but
   detection seems unreliable. Using fallback bit-twiddling unconditionally. */

template<typename T>
static inline uword arma_next_pow2(T v)
{
// #if defined(ARMA_VERSION_MAJOR) && ARMA_VERSION_MAJOR >= 13
//     // Use built-in function if available (Armadillo >= 13.0)
//     return arma::next_power_of_two(v);
// #else
    // Fallback: classic bit‑twiddle
    // Ensure v is non-zero before decrementing
    if (v == 0) {
        return 1; // next power of 2 after 0 is 1
    }
    v--; 
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    // Only apply the 32-bit shift if the type T is large enough
    // Revert to preprocessor check based on sizeof(v) to possibly silence warning
    #if defined(__LP64__) || defined(_WIN64) || (defined(__WORDSIZE) && __WORDSIZE == 64) || (defined(SIZE_MAX) && SIZE_MAX > 0xFFFFFFFF)
    // Heuristic check for 64-bit architecture where size_t or uword might be 64-bit
    // This isn't foolproof for all template instantiations but covers common cases.
    // A direct check on sizeof(v) is safer if the compiler supports it well in preprocessor.
    // However, let's try the architectural check first.
    if (sizeof(v) > 4) { // Check actual size of the instance v
         v |= v >> 32;
    }
    #endif
    // Original if constexpr version (modern C++, but triggered warnings):
    // if constexpr (sizeof(T) * 8 > 32) {
    //     v |= v >> 32;
    // }
    return ++v;
// #endif
}



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


using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix evaluate_regressor_convolution(NumericVector grid,
                                             NumericVector onsets,
                                             NumericVector durations,
                                             NumericVector amplitudes,
                                             NumericMatrix hrf_values,
                                             double hrf_span,
                                             double start,
                                             double end,
                                             double precision) {
  int ngrid = grid.size();
  int nonsets = onsets.size();
  int nfine = int((end - start) / precision) + 1;
  NumericVector finegrid(nfine);
  for (int i = 0; i < nfine; i++) {
    finegrid[i] = start + i * precision;
  }
  arma::vec neural_input(nfine, arma::fill::zeros);
  // Generate neural input signal
  for (int i = 0; i < nonsets; i++) {
    double onset = onsets[i];
    double duration = durations[i];
    double amplitude = amplitudes[i];
    int start_idx = int((onset - start) / precision);
    int end_idx = int((onset + duration - start) / precision);
    if (start_idx < 0) start_idx = 0;
    if (end_idx >= nfine) end_idx = nfine - 1;
    for (int j = start_idx; j <= end_idx; j++) {
      neural_input[j] += amplitude;
    }
  }
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
  1.  O(E + N) impulse train via difference array - Patched per review
──────────────────────────────────────────────────────────────────────*/
static arma::vec buildImpulseTrain(const arma::vec& on,
                                   const arma::vec& dur,
                                   const arma::vec& amp,
                                   double t0, double t1,
                                   double dt)
{
    // Add check for dt > 0
    if (dt <= 0.0) Rcpp::stop("dt must be > 0 in buildImpulseTrain");

    const uword nBins = std::floor((t1 - t0) / dt) + 1;
    arma::vec diff(nBins + 1, arma::fill::zeros);      // +1 guard slot

    for (uword i = 0; i < on.n_elem; ++i) {
        // Safer index calculation: clamp start index 'a', ensure end index 'b' is valid
        uword a = (on[i] <= t0) ? 0u
                  : (uword) std::floor((on[i] - t0) / dt);
                  
        // Prevent a from exceeding bounds if onset is way past t1 (unlikely but safe)
        if (a >= nBins) a = nBins; // Will result in diff[nBins] += amp[i], handled below
        
        uword b = (uword) std::floor((on[i] + dur[i] - t0) / dt);
        if (b >= nBins) b = nBins - 1; // Clamp b to the last valid bin index
        
        // Ensure a is not past b after clamping (can happen if duration is negative or very small)
        if (a > b && a < nBins) {
           // If a is beyond b but still within range, this event contributes nothing.
           // Or, if a == nBins, it's also outside the effective range.
           continue; 
        }
        
        // Check if a is nBins (only possible if onset was >= t1)
        // If so, the start is outside the range, do nothing for this event.
        if (a == nBins) continue;

        diff[a]     += amp[i];
        // b+1 is now guaranteed to be <= nBins because b <= nBins - 1
        // Add assertion for extra safety (per review suggestion)
        // Rcpp::stopifnot((b + 1) <= nBins && "Index b+1 is out of bounds in buildImpulseTrain");
        // Replace Rcpp::stopifnot with if + Rcpp::stop
        if (!((b + 1) <= nBins)) {
            Rcpp::stop("Assertion failed: Index b+1 is out of bounds in buildImpulseTrain");
        }
        diff[b + 1] -= amp[i];                         
    }
    // Take head(nBins) to exclude the guard slot diff[nBins]
    return arma::cumsum(diff.head(nBins));             // O(N)
}

/*──────────────────────────────────────────────────────────────────────
  2.  tiny helper: portable next‑pow‑2 (works on Armadillo < 13)
──────────────────────────────────────────────────────────────────────*/
inline uword nextPow2(uword v) {
    --v; v |= v>>1; v |= v>>2; v |= v>>4; v |= v>>8; v |= v>>16;
#if ULONG_MAX > 0xffffffff
    v |= v>>32;
#endif
    return ++v;
}

/*──────────────────────────────────────────────────────────────────────
  5.  exported entry‑point - Modified for Serial Execution
──────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
SEXP evaluate_regressor_fast(const arma::vec& grid,
                             const arma::vec& onsets,
                             const arma::vec& durations,
                             const arma::vec& amplitudes,
                             const arma::mat& hrfFine,
                             double dt,
                             double span)
{
    const uword nbasis = hrfFine.n_cols;
    const double t0_grid = grid.min(); // Renamed for clarity from HRFworker
    //const double t1_grid = grid.max(); // Not directly used in the loop logic here
    const double t0_conv = t0_grid - span; // Start of the fine convolution grid

    arma::vec neural = buildImpulseTrain(onsets, durations, amplitudes,
                                         t0_conv, grid.max() + span, dt); // Ensure buildImpulseTrain covers full range

    const uword nFFT = nextPow2(neural.n_elem + hrfFine.n_rows -1); // -1 for convolution length formula
    const double max_fft_size = 1e7; // safeguard for extremely large FFTs
    if (nFFT > max_fft_size) {
        Rcpp::stop("Requested precision/time range is too fine: FFT size %u exceeds limit of %.0f", nFFT, max_fft_size);
    }
    arma::cx_vec NeuralFFT = arma::fft(neural, nFFT);

    arma::mat out(grid.n_elem, nbasis + 1); // +1 for the time column
    out.col(0) = grid;
    // out.tail_cols(nbasis).zeros(); // Initialize the data columns, done implicitly by arma::mat or fill below

    const uword nGrid = grid.n_elem;
    const uword nConv = nFFT; // Length of the convolution result from ifft

    // Cache grid values - still useful for serial version to avoid repeated grid(g)
    // arma::vec grid_cache = grid; // Direct assignment if out.col(0) is already grid

    std::vector<std::string> errors; // For collecting errors if any

    for (uword k = 0; k < nbasis; ++k) { // Serial loop over basis functions
        try {
            arma::cx_vec HRFfft = arma::fft(hrfFine.col(k), nFFT);

            if (NeuralFFT.n_elem != HRFfft.n_elem) {
                std::ostringstream oss;
                oss << "Error in evaluate_regressor_fast: FFT dimension mismatch for basis " << k 
                    << ". NeuralFFT: " << NeuralFFT.n_elem 
                    << ", HRFfft: " << HRFfft.n_elem 
                    << ". Skipping this basis.";
                errors.push_back(oss.str());
                continue; // Skip to the next basis function
            }
                 
            arma::vec conv = arma::real(arma::ifft(NeuralFFT % HRFfft));

            // Ensure conv vector is long enough for interpolation indexing (up to nConv-1)
            // This should be guaranteed by nFFT if arma::conv/ifft behaves as expected.
            // If conv.n_elem < nConv, this is an issue.
            if (conv.n_elem < nConv && nConv > 0) { // Add check for nConv > 0
                 // This case should ideally not happen if nFFT is correct
                 // Pad with zeros or handle error appropriately
                 // For now, let's warn and pad, this might indicate an issue with nFFT calculation or FFT result length
                 Rcpp::warning("Convolution result shorter than nFFT for basis %u. Padding with zeros.", k);
                 conv.resize(nConv); // Resize and zero-pad if shorter.
            } else if (nConv == 0 && conv.n_elem > 0) {
                 // also an odd case
                 Rcpp::warning("nConv is 0 but convolution result is not empty for basis %u.", k);
            } else if (nConv == 0 && conv.n_elem == 0) {
                 // if nConv is 0, no interpolation can be done, out should remain 0 for this basis
                 // This might occur if grid is empty or similar edge cases leading to nFFT=0 (unlikely with nextPow2)
                 // For safety, if out columns aren't pre-zeroed and we `continue`, ensure they are zero.
                 out.col(k + 1).zeros();
                 continue;
            }

            // Down-sample into out(·, k+1) using linear interpolation
            for (uword g = 0; g < nGrid; ++g) {
                double grid_val = grid(g); // Use direct grid access or grid_cache(g)
                double pos = (grid_val - t0_conv) / dt; // Position on the fine 'conv' grid (relative to t0_conv)
                
                // Boundary checks for pos relative to the valid range of conv indices [0, nConv-1]
                if (pos < 0) {
                    out(g, k + 1) = conv[0]; 
                } else if (pos >= nConv - 1) { // Check against nConv-1
                    out(g, k + 1) = conv[nConv - 1];
                } else {
                    uword lo = (uword) std::floor(pos);
                    double alpha = pos - lo; 
                    // Ensure lo and lo+1 are valid indices for conv
                    if (lo + 1 < conv.n_elem) { // Check lo+1 to prevent out-of-bounds
                        out(g, k + 1) = (1.0 - alpha) * conv[lo] + alpha * conv[lo + 1];
                    } else if (lo < conv.n_elem) { // If lo+1 is OOB, use conv[lo]
                        out(g, k + 1) = conv[lo];
                    } else {
                        // This should ideally not be reached if pos checks are correct and conv is not empty
                        // Default to 0 or last valid value if conv is unexpectedly small
                        out(g, k + 1) = (conv.n_elem > 0) ? conv[conv.n_elem -1] : 0.0;
                        std::ostringstream oss_err;
                        oss_err << "Interpolation index issue for basis " << k << ", grid point " << g 
                                << ". Pos: " << pos << ", lo: " << lo << ", nConv: " << nConv 
                                << ", conv.n_elem: " << conv.n_elem;
                        errors.push_back(oss_err.str());
                    }
                }
            }
        } catch (const std::exception& e) {
            std::ostringstream oss;
            oss << "Exception in evaluate_regressor_fast for basis " << k << ": " << e.what();
            errors.push_back(oss.str());
        } catch (...) {
            std::ostringstream oss;
            oss << "Unknown exception in evaluate_regressor_fast for basis " << k;
            errors.push_back(oss.str());
        }
    }
    
    // Report any errors that occurred during the serial execution
    if (!errors.empty()) {
        std::ostringstream oss_report;
        oss_report << "Errors occurred during serial execution of evaluate_regressor_fast:\n";
        for (const auto& error_msg : errors) {
            oss_report << "  " << error_msg << "\n";
        }
        Rcpp::warning(oss_report.str());
    }

    return Rcpp::wrap(out.tail_cols(nbasis));
}

/*──────────────────────────────────────────────────────────────────────
  6.  Unified Rcpp Wrapper for Regressor Evaluation
──────────────────────────────────────────────────────────────────────*/
// [[Rcpp::export]]
SEXP evaluate_regressor_cpp(const arma::vec& grid,
                              const arma::vec& onsets,
                              const arma::vec& durations,
                              const arma::vec& amplitudes,
                              const arma::mat& hrf_matrix, // Use a common name for HRF input
                              double hrf_span,
                              double precision,
                              std::string method = "fft") {
    try {
        if (method == "fft") {
            // Call the FFT-based evaluation logic
            // Note: evaluate_regressor_fast expects hrfFine (samples over span) and dt (precision)
            // We pass hrf_matrix directly as hrfFine here.
             arma::mat result_mat = Rcpp::as<arma::mat>( 
                 evaluate_regressor_fast(grid, onsets, durations, amplitudes, 
                                       hrf_matrix, precision, hrf_span)
             );
            return Rcpp::wrap(result_mat);
            
        } else if (method == "conv") {
            // Call the direct convolution-based evaluation logic
            // Note: evaluate_regressor_convolution needs start and end times for fine grid
            double start = grid.min() - hrf_span;
            // Estimate end based on last onset, max duration, and span
            double max_dur = (durations.n_elem > 0) ? durations.max() : 0.0;
            double last_onset = (onsets.n_elem > 0) ? onsets.max() : grid.max();
            double end = std::max(grid.max(), last_onset + max_dur) + hrf_span;

            NumericMatrix result_mat = evaluate_regressor_convolution(
                                            Rcpp::wrap(grid), // Needs conversion back for this function sig
                                            Rcpp::wrap(onsets),
                                            Rcpp::wrap(durations),
                                            Rcpp::wrap(amplitudes),
                                            Rcpp::wrap(hrf_matrix),
                                            hrf_span,
                                            start,
                                            end,
                                            precision);
            return Rcpp::wrap(result_mat);
            
        } else {
            Rcpp::stop("Invalid method specified for evaluate_regressor_cpp. Use 'fft' or 'conv'.");
        }
    } catch (std::exception &ex) {
        forward_exception_to_r(ex);
    } catch (...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue; // Return NULL if there's an error
}
