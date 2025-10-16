// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace paratoric::statistics {

// next power of two
inline std::size_t next_pow2(std::size_t n) {
    if (n == 0) return 1;
    --n;
    n |= n >> 1;  n |= n >> 2;  n |= n >> 4;
    n |= n >> 8;  n |= n >> 16;
#if SIZE_MAX > 0xFFFFFFFFu
    n |= n >> 32;
#endif
    return n + 1;
}

// In-place radix-2 FFT. Could use faster external libraries like MKL but not worth the trouble for this
inline void fft_inplace(std::vector<std::complex<double>>& a, bool inverse) {
    const std::size_t n = a.size();
    if (n == 0) return;
    // n must be power of two
#ifndef NDEBUG
    if ((n & (n - 1)) != 0) throw std::runtime_error("fft_inplace: size not power of two");
#endif

    // bit-reversal permutation
    std::size_t j = 0;
    for (std::size_t i = 1; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // Cooley–Tukey stages
    const double pi = 3.141592653589793238462643383279502884;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = 2.0 * pi / static_cast<double>(len) * (inverse ? -1.0 : 1.0);
        const std::complex<double> wlen = { std::cos(ang), std::sin(ang) };

        for (std::size_t i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            const std::size_t half = len >> 1;
            for (std::size_t j = 0; j < half; ++j) {
                const std::complex<double> u = a[i + j];
                const std::complex<double> v = a[i + j + half] * w;
                a[i + j]         = u + v;
                a[i + j + half]  = u - v;
                w *= wlen;
            }
        }
    }

    if (inverse) {
        const double inv_n = 1.0 / static_cast<double>(n);
        for (auto& x : a) x *= inv_n; // scale only on IFFT
    }
}

// unbiased FFT-based version O(N log N))
inline std::vector<double> get_autocorrelation_function_fft(std::span<const double> data) {
    const std::size_t N = data.size();
    if (N == 0) return {0.0};

    // mean-center
    const double mean = std::accumulate(data.begin(), data.end(), 0.0) / static_cast<double>(N);

    std::vector<double> y(N);
    for (std::size_t i = 0; i < N; ++i) y[i] = data[i] - mean;

    // variance check (constant sequence)
    double c0_sum = 0.0;
    for (double v : y) c0_sum += v * v;
    if (c0_sum == 0.0) {
        // Define rho(0)=1 and no higher lags (undefined). Keep API compact.
        return std::vector<double>{1.0};
    }

    // zero-padding
    const std::size_t M = next_pow2(N * 2);
    std::vector<std::complex<double>> a(M, std::complex<double>(0.0, 0.0));
    for (std::size_t i = 0; i < N; ++i) a[i] = std::complex<double>(y[i], 0.0);

    // FFT --> power spectrum --> IFFT to get linear autocovariance sums
    fft_inplace(a, /*inverse=*/false);
    for (std::size_t k = 0; k < M; ++k) a[k] *= std::conj(a[k]); // |Y|^2
    fft_inplace(a, /*inverse=*/true);

    // a[k].real() now holds sum_{i=0}^{N-1-k} y[i]*y[i+k] for k < N (thanks to zero padding)

    const double C0 = c0_sum / static_cast<double>(N); // unbiased denominator for normalization
    std::vector<double> acf;
    acf.reserve(N);

    // lag 0
    acf.push_back(1.0);

    // lags 1..N-1: unbiased covariance Ck = sum / (N - k), then normalize by C0
    for (std::size_t k = 1; k < N; ++k) {
        const double sum = std::real(a[k]);                  // raw sum at lag k
        const double Ck  = sum / static_cast<double>(N - k); // unbiased
        acf.push_back(Ck / C0);
    }

    return acf;
}

// Naive version O(n^2) 
inline std::vector<double> get_autocorrelation_function_naive(std::span<const double> data) {
    const size_t N = data.size();
    if (N == 0) return {0.0};

    // Mean
    const double mean = std::accumulate(data.begin(), data.end(), 0.0) / static_cast<double>(N);

    // Centered variance sum (C(0) * N)
    double c0_sum = 0.0;
    for (double x : data) {
        const double y = x - mean;
        c0_sum += y * y;
    }

    if (c0_sum == 0.0) {
        std::vector<double> acf(1, 1.0);
        return acf;
    }

    std::vector<double> acf(N);
    acf[0] = 1.0;

    for (size_t k = 1; k < N; ++k) {
        double num = 0.0;
        for (size_t i = 0; i + k < N; ++i) {
            num += (data[i] - mean) * (data[i + k] - mean);
        }
        // Unbiased covariance at lag k
        const double Ck = num / static_cast<double>(N - k);
        const double C0 = c0_sum / static_cast<double>(N);
        acf[k] = Ck / C0;
    }
    return acf;
}

/**
 * @brief This function returns the autocorrelation function of the input data.
 * 
 * @param data the data of which the autocorrelation function will be calculated
 * 
 * @return the autocorrelation function
 */
inline std::vector<double> get_autocorrelation_function(std::span<const double> data, std::string mode = "fft") { 
    if (mode == "fft") {
        return get_autocorrelation_function_fft(data);
    } else if (mode == "naive") {
        return get_autocorrelation_function_naive(data);
    } else {
        throw std::invalid_argument(std::format("Mode {} not found.", mode));
        return std::vector<double>{};
    }
    

}

/**
 * @brief This function returns the integrated autocorrelation time of the input autocorrelation function acf.
 * 
 * @param acf the autocorrelation function which will be integrated
 * 
 * @return the integrated autocorrelation time
 */
inline double get_autocorrelation_time(const std::vector<double>& acf) {
    const size_t N = acf.size();
    if (N == 0 || std::abs(acf[0]) == 0.0) return 0.5;

    // Find window W via Geyer's IPS: stop at first m with rho(2m)+rho(2m+1) <= 0
    size_t W = 0;
    for (size_t m = 1; 2*m + 1 < N; ++m) {
        const double pair_sum = acf[2*m] + acf[2*m + 1];
        if (pair_sum <= 0.0) {
            W = 2*m;               // include up to 2m-1
            break;
        }
    }
    if (W == 0) {
        // Fallback: first non-positive rho(k), else all
        size_t cutoff = N;
        for (size_t k = 1; k < N; ++k) {
            if (acf[k] <= 0.0) { cutoff = k; break; }
        }
        W = (cutoff == N ? N - 1 : cutoff - 1);
    }

    double tau = 0.5;
    for (size_t k = 1; k <= W; ++k) {
        tau += acf[k];
    }
    if (!std::isfinite(tau)) return 0.5;
    return std::max(0.5, tau);
}


} // namespace paratoric::statistics



