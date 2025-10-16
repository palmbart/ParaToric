// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#pragma once

#include "rng/rng.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

namespace paratoric::statistics {

using RNG = paratoric::rng::RNG;

/**
 * Based on the following publications:
 *
 * Politis, D. N., & White, H. (2004). 
 * Automatic Block-Length Selection for the Dependent Bootstrap. 
 * Econometric Reviews, 23(1), 53-70. https://doi.org/10.1081/ETC-120028836
 *
 * Patton, A., Politis, D. N., & White, H. (2009). 
 * Correction to “Automatic Block-Length Selection for the Dependent Bootstrap” by D. Politis and H. White. 
 * Econometric Reviews, 28(4), 372-375. https://doi.org/10.1080/07474930802459016
 * 
 * @param data the data which will be bootstrapped using the stationary bootstrap
 * 
 * @return the optimal block length for the stationary bootstrap
 */
inline double opt_block_length(const std::vector<double>& data) {
    // We use the convention 1/N for the autocorrelation here, not 1/(N-k)
    const size_t N = data.size();
    if (N < 2) return 1.0;

    double b_max = std::ceil(std::min(3 * std::sqrt((double)N), N / 3.0));
    int K_n = std::max(5, static_cast<int>(std::log(N)));
    int m_max = static_cast<int>(std::ceil(std::sqrt((double)N))) + K_n;
    const double c = 2.0;

    std::vector<double> R_k(m_max + 1);
    std::vector<double> ac_abs(m_max + 1);

    double mean = std::accumulate(data.begin(), data.end(), 0.0) / N;
    std::vector<double> delta(N);
    for (size_t i = 0; i < N; ++i) {
        delta[i] = data[i] - mean;
    }

    int opt_m = -1;

    for (int k = 0; k <= m_max; ++k) {
        double d1_sq = 0.0, d2_sq = 0.0, summand = 0.0;
        for (size_t i = 0; i + k + 1 < N; ++i) {
            d1_sq += delta[i + k + 1] * delta[i + k + 1];
            d2_sq += delta[i] * delta[i];
        }
        for (size_t i = 0; i + k < N; ++i) {
            summand += delta[i + k] * delta[i];
        }

        R_k[k] = summand / N;
        if (d1_sq * d2_sq == 0.0) {
            return 1.0;
        }

        ac_abs[k] = std::abs(summand) / std::sqrt(d1_sq * d2_sq);

        if (k >= K_n && opt_m < 0) {
            bool small = true;
            double threshold = c * std::sqrt(std::log(N) / N);
            for (int j = k - K_n; j < k; ++j) {
                if (ac_abs[j] >= threshold) {
                    small = false;
                    break;
                }
            }
            if (small) opt_m = k - K_n;
        }
    }

    int m = (opt_m >= 0 ? 2 * std::max(opt_m, 1) : m_max);
    m = std::min(m, m_max);

    double G = 0.0;
    double g0_hat = R_k[0];
    for (int k = 1; k <= m; ++k) {
        double h = std::min(1.0, 2.0 * (1.0 - double(k) / m));
        G += 2.0 * h * k * R_k[k];
        g0_hat += 2.0 * h * R_k[k];
    }

    const double det_c = 2.0; // stationary bootstrap correction
    double D_sb = det_c * g0_hat * g0_hat;
    if (D_sb == 0.0) return 1.0;

    double b_opt = std::cbrt((2.0 * G * G) / D_sb) * std::cbrt((double)N);
    return std::min(b_opt, b_max);
}

/**
 * Based on the following publication:
 *
 * Politis, D. N., & Romano, J. P. (1994). 
 * The Stationary Bootstrap. Journal of the American Statistical Association, 
 *  89(428), 1303-1313. https://doi.org/10.1080/01621459.1994.10476870
 * 
 * @param data the data which will be bootstrapped using the stationary bootstrap
 * @param average_block_length the block length (parameter of the stationary bootstrap)
 * @param sample_size the sample size of the bootstrap
 * @param rng the Mersenne Twister random number generator
 * 
 * @return a vector with the the bootstrapped indices of the original data
 */
inline std::vector<size_t> stationary_bootstrap(
    const std::vector<double>& data, double average_block_length, 
    size_t sample_size, std::shared_ptr<RNG> rng
) {
    size_t N = data.size();
    std::uniform_real_distribution<double> udist(0.0, 1.0);
    std::uniform_int_distribution<size_t> idist(0, N - 1);
    std::vector<size_t> result(sample_size);
    double p = 1.0 / average_block_length;
    size_t start = idist(*rng), offset = 0;
    for (size_t i = 0; i < sample_size; ++i) {
        if (i > 0 && udist(*rng) > p) {
            ++offset;
        } else {
            start = idist(*rng);
            offset = 0;
        }
        result[i] = (start + offset) % N;
    }
    return result;
}

/**
 * @brief This function will perform n_iter stationary bootstraps of input data and return the mean, the Binder ratio and the standard error, respectively.
 * 
 * Autocorrelation effects are included in the stationary bootstrap.
 * 
 * @param data the data which will be bootstrapped using the stationary bootstrap
 * @param n_iter (optional) - the number of bootstraps, defaults to 1000
 * 
 * @return a tuple of the mean, the standard error of the mean, the Binder ratio and the standard error of the Binder ratio
 */
inline std::tuple<double, double, double, double> get_bootstrap_statistics(
    const std::vector<double>& data, std::shared_ptr<RNG> rng, 
    size_t n_iter = 1000
) {
    const size_t N = data.size();
    double blk_len = opt_block_length(data);
    // Naive point estimate (original sample mean)
    const double m_hat = std::accumulate(data.begin(), data.end(), 0.0) / N;
    double sum2_raw = 0.0, sum4_raw = 0.0; 
    for (double v : data) { sum2_raw += v*v; sum4_raw += v*v*v*v; }
    double e2_raw = sum2_raw / N, e4_raw = sum4_raw / N;
    if (e2_raw == 0.0) { e2_raw = 1e-5; e4_raw = 0.0; }
    const double binder_hat = e4_raw / (e2_raw * e2_raw);

    std::vector<double> means, binders;
    means.reserve(n_iter); binders.reserve(n_iter);
    for (size_t it = 0; it < n_iter; ++it) {
        auto idx = stationary_bootstrap(data, blk_len, N, rng);
        double sum = 0, sum2 = 0, sum4 = 0;
        for (auto i: idx) {
            double v = data[i]; sum += v;
            sum2 += v*v; sum4 += v*v*v*v;
        }
        double m = sum / N;
        means.push_back(m);
        double e2 = sum2 / N, e4 = sum4 / N;
        if (e2 == 0.0) { e2 = 1e-5; e4 = 0; }
        binders.push_back(e4 / (e2*e2));
    }
    auto stats = [&](const std::vector<double>& v) {
        double mu = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double var = 0;
        for (auto x: v) var += (x-mu)*(x-mu);
        var /= (v.size() - 1);
        return std::make_pair(mu, std::sqrt(var));
    };
    auto [m_mean, m_std] = stats(means);
    auto [b_mean, b_std] = stats(binders);
    // Bias-corrected mean: 2*m_hat - mean(bootstrap means)
    const double m_mean_bias_corrected = 2.0 * m_hat - m_mean;
    const double b_mean_bias_corrected = 2.0 * binder_hat  - b_mean;
    return {m_mean_bias_corrected, m_std, b_mean_bias_corrected, b_std};
}

/**
 * @brief This function will perform n_iter stationary bootstraps of input half and full representing the half/full Wilson/'t Hooft loops of the Fredenhagen-Marcu. It returns the mean, the Binder ratio and the standard error, respectively.
 * 
 * Autocorrelation effects are included in the stationary bootstrap. Due to the non-trivial structure of the Fredenhagen-Marcu fraction, a custom function is required.
 * 
 * @param half the half-loop (numerator) of the Fredenhagen-Marcu which will be bootstrapped using the stationary bootstrap
 * @param full the full-loop (denominator) of the Fredenhagen-Marcu which will be bootstrapped using the stationary bootstrap
 * @param n_iter (optional) - the number of bootstraps, defaults to 1000
 * 
 * @return a tuple of the mean, the standard error of the mean, the Binder ratio and the standard error of the Binder ratio
 */
inline std::tuple<double,double,double,double> get_bootstrap_statistics_fm(
    const std::vector<double>& half, const std::vector<double>& full, 
    std::shared_ptr<RNG> rng, size_t n_iter = 1000
) {
    const size_t N = half.size();
    if (full.size() != N) throw std::invalid_argument("half/full length mismatch");

    double blk_len = std::max(opt_block_length(half), opt_block_length(full));

    double sum_numerator_hat = 0;
    double sum_denominator_hat = 0;
    for (size_t i = 1; i < half.size(); ++i) {
        sum_numerator_hat += half[i];
        sum_denominator_hat += full[i];
    }

    double m_numerator_hat = sum_numerator_hat / N;
    double m_denominator_hat = std::sqrt(std::abs(sum_denominator_hat / N));
    double m_hat = m_denominator_hat > 0 ? m_numerator_hat/m_denominator_hat : 0.;

    std::vector<double> means, binders;
    means.reserve(n_iter); binders.reserve(n_iter);

    for (size_t it = 0; it < n_iter; ++it) {
        auto idx = stationary_bootstrap(half, blk_len, N, rng);
        double sum_numerator = 0;
        double sum_denominator = 0;
        for (auto i: idx) {
            sum_numerator += half[i];
            sum_denominator += full[i];
        }

        double m_numerator = sum_numerator / N;
        double m_denominator = std::sqrt(std::abs(sum_denominator / N));
        double m = m_denominator > 0 ? m_numerator/m_denominator : 0.;
        means.push_back(m);
        binders.push_back(0.);
    }
    auto stats = [&](const std::vector<double>& v) {
        double mu = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double var = 0;
        for (auto x: v) var += (x-mu)*(x-mu);
        var /= (v.size() - 1);
        return std::make_pair(mu, std::sqrt(var));
    };
    auto [m_mean, m_std] = stats(means);
    auto [b_mean, b_std] = stats(binders);
    // Bias-corrected mean: 2*m_hat - mean(bootstrap means)
    const double m_mean_bias_corrected = 2.0 * m_hat - m_mean;
    return {m_mean_bias_corrected, m_std, b_mean, b_std};
}

/**
 * @brief This function will perform n_iter stationary bootstraps of input half and full representing the connected and disconnected part of the susceptibility. It returns the mean, the Binder ratio and the standard error, respectively.
 * 
 * Autocorrelation effects are included in the stationary bootstrap. Due to the non-trivial structure of the susceptibility fraction, a custom function is required.
 * 
 * @param half the half-loop (numerator) of the Fredenhagen-Marcu which will be bootstrapped using the stationary bootstrap
 * @param full the full-loop (denominator) of the Fredenhagen-Marcu which will be bootstrapped using the stationary bootstrap
 * @param n_iter (optional) - the number of bootstraps, defaults to 1000
 * 
 * @return a tuple of the mean, the standard error of the mean, the Binder ratio and the standard error of the Binder ratio
 */
inline std::tuple<double,double,double,double> get_bootstrap_statistics_susceptibility(
    const std::vector<double>& kL, const std::vector<double>& kR, double beta, 
    double h, double Nsites, std::shared_ptr<RNG> rng, size_t n_iter = 1000
) {
    const size_t N = kL.size();
    if (kR.size() != N)
        throw std::invalid_argument("kL/kR length mismatch");
    if (N < 2)
        throw std::invalid_argument("need >=2 samples for bootstrap");

    double blk_len = std::max(opt_block_length(kL), opt_block_length(kR));

    // --- original-sample targets for bias correction ---
    double sumL_raw = 0.0, sumR_raw = 0.0, sumLR_raw = 0.0;
    double sumLR2_raw = 0.0, sumLR4_raw = 0.0;
    for (size_t i = 0; i < N; ++i) {
        const double L = kL[i];
        const double R = kR[i];
        const double LR = L * R;
        sumL_raw  += L;
        sumR_raw  += R;
        sumLR_raw += LR;
        sumLR2_raw += LR * LR;
        sumLR4_raw += LR * LR * LR * LR;
    }
    const double meanL_raw  = sumL_raw  / N;
    const double meanR_raw  = sumR_raw  / N;
    const double meanLR_raw = sumLR_raw / N;
    const double cov_hat = (meanLR_raw - meanL_raw * meanR_raw);
    double e2_raw = sumLR2_raw / N, e4_raw = sumLR4_raw / N;
    if (e2_raw == 0.0) { e2_raw = 1e-12; e4_raw = 0.0; }
    const double binder_hat = e4_raw / (e2_raw * e2_raw);

    std::vector<double> cov_samples;
    std::vector<double> binder_samples;
    cov_samples.reserve(n_iter);
    binder_samples.reserve(n_iter);

    for (size_t it = 0; it < n_iter; ++it) {
        auto idx = stationary_bootstrap(kL, blk_len, N, rng);

        double sumL = 0.0, sumR = 0.0, sumLR = 0.0;
        double sumLR2 = 0.0, sumLR4 = 0.0;

        for (auto i : idx) {
            const double L = kL[i];
            const double R = kR[i];
            const double LR = L * R;

            sumL   += L;
            sumR   += R;
            sumLR  += LR;
            sumLR2 += LR * LR;
            sumLR4 += LR * LR * LR * LR;
        }

        const double meanL  = sumL  / N;
        const double meanR  = sumR  / N;
        const double meanLR = sumLR / N;

        const double covLR = (meanLR - meanL * meanR);
        cov_samples.push_back(covLR);

        double e2 = sumLR2 / N;
        double e4 = sumLR4 / N;
        if (e2 == 0.0) { e2 = 1e-12; e4 = 0.0; }
        binder_samples.push_back(e4 / (e2 * e2));
    }

    auto summarize = [](const std::vector<double>& v) {
        const double mu = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double var = 0.0;
        for (double x : v) {
            const double d = x - mu;
            var += d * d;
        }
        var /= (v.size() > 1 ? (v.size() - 1) : 1);
        return std::make_pair(mu, std::sqrt(var));
    };

    auto [cov_mean,    cov_std]    = summarize(cov_samples);
    auto [binder_mean, binder_std] = summarize(binder_samples);
    // Bias-corrected mean: 2*m_hat - mean(bootstrap means)
    const double cov_mean_bias_corrected    = 2.0 * cov_hat    - cov_mean;
    const double binder_mean_bias_corrected = 2.0 * binder_hat - binder_mean;

    return {cov_mean_bias_corrected, cov_std, binder_mean_bias_corrected, binder_std};
}


// TODO
inline std::tuple<double,double,double,double> bootstrap_offdiag_susceptibility(
    const std::vector<double>& kvec, double beta, double h, 
    double Nsites, std::shared_ptr<RNG> rng, size_t n_iter = 1000
) {
    if (kvec.size() < 2) throw std::invalid_argument("need >=2 samples");

    double blk_len = opt_block_length(kvec);

    std::vector<double> chi_samples;
    std::vector<double> binder_samples;
    chi_samples.reserve(n_iter);
    binder_samples.reserve(n_iter);

    const size_t N = kvec.size();

    // --- original-sample chi for bias correction ---
    double sumk_raw = 0.0, sumk2_raw = 0.0;
    for (double k : kvec) { sumk_raw += k; sumk2_raw += k*k; }
    const double meank_raw  = sumk_raw  / N;
    const double meank2_raw = sumk2_raw / N;
    const double varK_raw   = meank2_raw - meank_raw * meank_raw;
    const double chi_hat    = (varK_raw - meank_raw) / (beta * Nsites * h * h);

    for (size_t it = 0; it < n_iter; ++it) {
        auto idx = stationary_bootstrap(kvec, blk_len, N, rng);

        double sumk = 0.0;
        double sumk2 = 0.0;
        //double sumk4 = 0.0;

        for (auto i : idx) {
            const double k = kvec[i];
            sumk  += k;
            sumk2 += k * k;
            //const double k2 = k * k;
            //sumk4 += k2 * k2;
        }

        const double meank  = sumk  / N;
        const double meank2 = sumk2 / N;
        //const double meank4 = sumk4 / N;

        // Moments  
        const double varK = (meank2 - meank * meank);
        // Convert to intensive susceptibility per site:
        const double chi = (varK - meank) / (beta * Nsites * h * h);

        chi_samples.push_back(chi);
        // TODO
        binder_samples.push_back( 0. );
    }

    auto summarize = [](const std::vector<double>& v) {
        const double mu = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double var = 0.0;
        for (double x : v) {
            const double d = x - mu;
            var += d * d;
        }
        var /= (v.size() > 1 ? (v.size() - 1) : 1);
        return std::make_pair(mu, std::sqrt(var));
    };

    auto [chi_mean, chi_std]       = summarize(chi_samples);
    auto [binder_mean, binder_std] = summarize(binder_samples);
    const double chi_mean_bias_corrected = 2.0 * chi_hat - chi_mean;
    return {chi_mean_bias_corrected, chi_std, binder_mean, binder_std};
}

} // namespace paratoric::statistics



