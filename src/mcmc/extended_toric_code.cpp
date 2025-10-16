// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#include "mcmc/extended_toric_code_qmc.hpp"
#include "paratoric/mcmc/extended_toric_code.hpp"
#include "paratoric/types/types.hpp"

#include <complex>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace paratoric {

struct ExtendedToricCode::Impl {
    template<char B, typename Fn>
    static auto with_backend(Fn&& fn) {
        ExtendedToricCodeQMC<B> q{};
        return fn(q);
    }

    template<typename Fn>
    static auto dispatch(char basis, Fn&& fn) {
        switch (basis) {
            case 'x': return with_backend<'x'>(std::forward<Fn>(fn));
            case 'z': return with_backend<'z'>(std::forward<Fn>(fn));
            default:  throw std::invalid_argument("basis must be 'x' or 'z'");
        }
    }
};

ExtendedToricCode::ExtendedToricCode()  : impl_(std::make_unique<Impl>()) {}
ExtendedToricCode::~ExtendedToricCode() = default;
ExtendedToricCode::ExtendedToricCode(ExtendedToricCode&&) noexcept = default;
ExtendedToricCode& ExtendedToricCode::operator=(ExtendedToricCode&&) noexcept = default;

auto ExtendedToricCode::get_thermalization(const Config& config)
-> Result
{

    return Impl::dispatch(config.lat_spec.basis, [&](auto& q) {
        return q.get_thermalization(config);
    });
}

auto ExtendedToricCode::get_sample(const Config& config)
-> Result
{

    return Impl::dispatch(config.lat_spec.basis, [&](auto& q) {
        return q.get_sample(config);
    });
}

auto ExtendedToricCode::get_hysteresis(const Config& config)
-> Result
{

    return Impl::dispatch(config.lat_spec.basis, [&](auto& q) {
        return q.get_hysteresis(config);
    });
}

} // namespace paratoric
