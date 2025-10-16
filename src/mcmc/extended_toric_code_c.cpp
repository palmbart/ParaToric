// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#include "paratoric/mcmc/extended_toric_code_c.h"

#include "paratoric/mcmc/extended_toric_code.hpp"
#include "paratoric/types/types.hpp"

#include <new>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <mutex>

using paratoric::ExtendedToricCode;
using paratoric::Config;
using paratoric::LatSpec;
using paratoric::ParamSpec;
using paratoric::SimSpec;
using paratoric::OutSpec;
using paratoric::Result;

/* ---------- thread-local error string ---------- */

static thread_local std::string g_last_error;

extern "C" const char* ptc_last_error(void) {
    return g_last_error.c_str();
}

static void set_error(const char* what) {
    g_last_error = (what ? what : "unknown error");
}

/* ---------- Opaque handle ---------- */

struct ptc_handle_t {
    std::unique_ptr<ExtendedToricCode> impl;
};

/* ---------- Small helpers: C<->C++ conversions ---------- */

static LatSpec to_cpp(const ptc_lat_spec_t& c) {
    LatSpec s{};
    s.basis        = c.basis;
    s.lattice_type = c.lattice_type ? std::string(c.lattice_type) : std::string{};
    s.system_size  = c.system_size;
    s.beta         = c.beta;
    s.boundaries   = c.boundaries ? std::string(c.boundaries) : std::string{};
    s.default_spin = c.default_spin;
    return s;
}

static ParamSpec to_cpp(const ptc_param_spec_t& c) {
    ParamSpec p{};
    p.mu    = c.mu;
    p.h     = c.h;
    p.J     = c.J;
    p.lmbda = c.lmbda;
    p.h_therm     = c.h_therm;
    p.lmbda_therm = c.lmbda_therm;

    if (c.h_hys && c.h_hys_len) p.h_hys.assign(c.h_hys, c.h_hys + c.h_hys_len);
    if (c.lmbda_hys && c.lmbda_hys_len) p.lmbda_hys.assign(c.lmbda_hys, c.lmbda_hys + c.lmbda_hys_len);
    return p;
}

static SimSpec to_cpp(const ptc_sim_spec_t& c) {
    SimSpec s{};
    s.N_samples         = c.N_samples;
    s.N_thermalization  = c.N_thermalization;
    s.N_between_samples = c.N_between_samples;
    s.N_resamples       = c.N_resamples;
    s.custom_therm      = c.custom_therm;
    s.seed              = c.seed;

    if (c.observables && c.N_observables) {
        s.observables.reserve(c.N_observables);
        for (size_t i = 0; i < c.N_observables; ++i) {
            s.observables.emplace_back(c.observables[i] ? c.observables[i] : "");
        }
    }
    return s;
}

static OutSpec to_cpp(const ptc_out_spec_t& c) {
    OutSpec o{};
    if (c.path_out) o.path_out = c.path_out;

    if (c.paths_out && c.N_paths_out) {
        o.paths_out.reserve(c.N_paths_out);
        for (size_t i = 0; i < c.N_paths_out; ++i) {
            if (c.paths_out[i]) o.paths_out.emplace_back(c.paths_out[i]);
        }
    }

    o.save_snapshots = c.save_snapshots;
    return o;
}

static Config to_cpp(const ptc_config_t& c) {
    Config cfg{};
    cfg.lat_spec       = to_cpp(c.lat);
    cfg.param_spec     = to_cpp(c.params);
    cfg.sim_spec       = to_cpp(c.sim);
    cfg.out_spec       = to_cpp(c.out);
    return cfg;
}

/* ----- Result (C++) -> (C) ----- */

static void free_series(ptc_series_t* s) {
    if (!s) return;
    for (size_t i = 0; i < s->nrows; ++i) {
        free(s->rows[i].data);
        s->rows[i].data = nullptr;
        s->rows[i].len = 0;
    }
    free(s->rows);
    s->rows = nullptr;
    s->nrows = 0;
}

static ptc_status_t fill_from_cpp(const Result& r, ptc_result_t* out) {
    if (!out) return PTC_STATUS_INVALID_ARGUMENT;

    /* Zero the output first. */
    *out = ptc_result_t{};

    auto copy_dvec = [](const std::vector<double>& v)->ptc_dvec_t {
        ptc_dvec_t dv{nullptr, 0};
        if (!v.empty()) {
            dv.data = (double*)malloc(sizeof(double) * v.size());
            if (!dv.data) return dv;
            for (size_t i = 0; i < v.size(); ++i) dv.data[i] = v[i];
            dv.len = v.size();
        }
        return dv;
    };

    /* series: vector<vector<variant<complex<double>, double>>> */
    {
        const auto& S = r.series;
        out->series.nrows = S.size();
        if (out->series.nrows) {
            out->series.rows = (ptc_scalar_vec_t*)calloc(out->series.nrows, sizeof(ptc_scalar_vec_t));
            if (!out->series.rows) return PTC_STATUS_NO_MEMORY;
            for (size_t i = 0; i < S.size(); ++i) {
                const auto& row = S[i];
                ptc_scalar_vec_t sv{nullptr, 0};
                if (!row.empty()) {
                    sv.data = (ptc_scalar_t*)malloc(sizeof(ptc_scalar_t) * row.size());
                    if (!sv.data) return PTC_STATUS_NO_MEMORY;
                    sv.len = row.size();
                    for (size_t j = 0; j < row.size(); ++j) {
                        const auto& val = row[j];
                        ptc_scalar_t s{};
                        if (std::holds_alternative<std::complex<double>>(val)) {
                            auto z = std::get<std::complex<double>>(val);
                            s.kind = PTC_SCALAR_COMPLEX;
                            s.re = z.real(); s.im = z.imag();
                        } else {
                            s.kind = PTC_SCALAR_REAL;
                            s.re = std::get<double>(val);
                            s.im = 0.0;
                        }
                        sv.data[j] = s;
                    }
                }
                out->series.rows[i] = sv;
            }
        }
    }

    out->acc_ratio    = copy_dvec(r.acc_ratio);
    out->mean      = copy_dvec(r.mean);
    out->mean_std  = copy_dvec(r.mean_std);
    out->binder    = copy_dvec(r.binder);
    out->binder_std= copy_dvec(r.binder_std);
    out->tau_int   = copy_dvec(r.tau_int);

    /* series_hys: vector<vector<vector<variant<...>>>> -> blocks of series */
    {
        const auto& H = r.series_hys;
        out->series_hys.nblocks = H.size();
        if (out->series_hys.nblocks) {
            out->series_hys.blocks = (ptc_series_t*)calloc(out->series_hys.nblocks, sizeof(ptc_series_t));
            if (!out->series_hys.blocks) return PTC_STATUS_NO_MEMORY;
            for (size_t b = 0; b < H.size(); ++b) {
                const auto& block = H[b]; /* vector<vector<variant<...>>> */
                ptc_series_t sblk{};
                sblk.nrows = block.size();
                if (sblk.nrows) {
                    sblk.rows = (ptc_scalar_vec_t*)calloc(sblk.nrows, sizeof(ptc_scalar_vec_t));
                    if (!sblk.rows) return PTC_STATUS_NO_MEMORY;
                    for (size_t i = 0; i < block.size(); ++i) {
                        const auto& row = block[i];
                        ptc_scalar_vec_t sv{nullptr, 0};
                        if (!row.empty()) {
                            sv.data = (ptc_scalar_t*)malloc(sizeof(ptc_scalar_t) * row.size());
                            if (!sv.data) return PTC_STATUS_NO_MEMORY;
                            sv.len = row.size();
                            for (size_t j = 0; j < row.size(); ++j) {
                                const auto& val = row[j];
                                ptc_scalar_t s{};
                                if (std::holds_alternative<std::complex<double>>(val)) {
                                    auto z = std::get<std::complex<double>>(val);
                                    s.kind = PTC_SCALAR_COMPLEX;
                                    s.re = z.real(); s.im = z.imag();
                                } else {
                                    s.kind = PTC_SCALAR_REAL;
                                    s.re = std::get<double>(val);
                                    s.im = 0.0;
                                }
                                sv.data[j] = s;
                            }
                        }
                        sblk.rows[i] = sv;
                    }
                }
                out->series_hys.blocks[b] = sblk;
            }
        }
    }

    /* double matrices for hysteresis aggregates */
    auto copy_mat = [&](const std::vector<std::vector<double>>& M)->ptc_dmat_t {
        ptc_dmat_t dm{nullptr, 0};
        if (!M.empty()) {
            dm.nrows = M.size();
            dm.rows = (ptc_dvec_t*)calloc(dm.nrows, sizeof(ptc_dvec_t));
            if (!dm.rows) { dm.nrows = 0; return dm; }
            for (size_t i = 0; i < M.size(); ++i) {
                dm.rows[i] = copy_dvec(M[i]);
            }
        }
        return dm;
    };

    out->mean_hys      = copy_mat(r.mean_hys);
    out->mean_std_hys  = copy_mat(r.mean_std_hys);
    out->binder_hys    = copy_mat(r.binder_hys);
    out->binder_std_hys= copy_mat(r.binder_std_hys);
    out->tau_int_hys   = copy_mat(r.tau_int_hys);

    return PTC_STATUS_OK;
}

/* ---------- public API ---------- */

extern "C" ptc_status_t ptc_create(ptc_handle_t** out_handle) {
    if (!out_handle) return PTC_STATUS_INVALID_ARGUMENT;
    try {
        auto h = new ptc_handle_t;
        h->impl = std::make_unique<ExtendedToricCode>();
        *out_handle = h;
        return PTC_STATUS_OK;
    } catch (const std::bad_alloc&) {
        set_error("allocation failed");
        return PTC_STATUS_NO_MEMORY;
    } catch (const std::exception& e) {
        set_error(e.what());
        return PTC_STATUS_INTERNAL_ERROR;
    }
}

extern "C" void ptc_destroy(ptc_handle_t* handle) {
    delete handle;
}

extern "C" void ptc_result_destroy(ptc_result_t* r) {
    if (!r) return;

    free_series(&r->series);

    free(r->acc_ratio.data);
    free(r->mean.data);
    free(r->mean_std.data);
    free(r->binder.data);
    free(r->binder_std.data);
    free(r->tau_int.data);

    for (size_t b = 0; b < r->series_hys.nblocks; ++b) {
        free_series(&r->series_hys.blocks[b]);
    }
    free(r->series_hys.blocks);

    if (r->mean_hys.rows)      { for (size_t i=0;i<r->mean_hys.nrows;i++)      free(r->mean_hys.rows[i].data);      free(r->mean_hys.rows); }
    if (r->mean_std_hys.rows)  { for (size_t i=0;i<r->mean_std_hys.nrows;i++)  free(r->mean_std_hys.rows[i].data);  free(r->mean_std_hys.rows); }
    if (r->binder_hys.rows)    { for (size_t i=0;i<r->binder_hys.nrows;i++)    free(r->binder_hys.rows[i].data);    free(r->binder_hys.rows); }
    if (r->binder_std_hys.rows){ for (size_t i=0;i<r->binder_std_hys.nrows;i++)free(r->binder_std_hys.rows[i].data);free(r->binder_std_hys.rows); }
    if (r->tau_int_hys.rows)   { for (size_t i=0;i<r->tau_int_hys.nrows;i++)   free(r->tau_int_hys.rows[i].data);   free(r->tau_int_hys.rows); }

    *r = ptc_result_t{};
}

static ptc_status_t call_cpp(const ptc_config_t* c,
                             ptc_result_t* out,
                             Result (*fn)(const Config&))
{
    if (!c || !out) return PTC_STATUS_INVALID_ARGUMENT;

    *out = ptc_result_t{};

    auto fail = [&](ptc_status_t st) {
        ptc_result_destroy(out);
        return st;
    };

    try {
        Config cfg = to_cpp(*c);
        Result r = fn(cfg);  // call the static method

        const ptc_status_t st = fill_from_cpp(r, out);
        if (st != PTC_STATUS_OK) return fail(st);
        return PTC_STATUS_OK;
    } catch (const std::invalid_argument& e) {
        set_error(e.what());
        return fail(PTC_STATUS_INVALID_ARGUMENT);
    } catch (const std::bad_alloc&) {
        set_error("allocation failed");
        return fail(PTC_STATUS_NO_MEMORY);
    } catch (const std::exception& e) {
        set_error(e.what());
        return fail(PTC_STATUS_RUNTIME_ERROR);
    } catch (...) {
        set_error("unknown exception");
        return fail(PTC_STATUS_INTERNAL_ERROR);
    }
}

extern "C" ptc_status_t ptc_get_thermalization(ptc_handle_t* handle,
                                               const ptc_config_t* config,
                                               ptc_result_t* out)
{
    if (!handle || !handle->impl) return PTC_STATUS_INVALID_ARGUMENT;
    return call_cpp(config, out, &ExtendedToricCode::get_thermalization);
}

extern "C" ptc_status_t ptc_get_sample(ptc_handle_t* handle,
                                       const ptc_config_t* config,
                                       ptc_result_t* out)
{
    if (!handle || !handle->impl) return PTC_STATUS_INVALID_ARGUMENT;
    return call_cpp(config, out, &ExtendedToricCode::get_sample);
}

extern "C" ptc_status_t ptc_get_hysteresis(ptc_handle_t* handle,
                                           const ptc_config_t* config,
                                           ptc_result_t* out)
{
    if (!handle || !handle->impl) return PTC_STATUS_INVALID_ARGUMENT;
    return call_cpp(config, out, &ExtendedToricCode::get_hysteresis);
}
