// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#define BOOST_TEST_MODULE TestTimeSearch

#include "lattice/time_search.hpp"

#include <boost/test/unit_test.hpp>

#include <array>
#include <random>

namespace paratoric::detail {
namespace {

// Independent reference: partition the interval at every event and evaluate
// each constant segment at its midpoint, counting flips with a linear scan.
double reference_integral(std::span<const WorldlineView> lines, double left, double right) {
    if (left >= right) return 0.0;
    std::vector<double> cuts{left, right};
    for (const auto& line : lines) {
        for (double tau : line.flips) {
            if (left < tau && tau < right) cuts.push_back(tau);
        }
    }
    std::sort(cuts.begin(), cuts.end());
    double result = 0.0;
    for (std::size_t i = 1; i < cuts.size(); ++i) {
        if (cuts[i - 1] == cuts[i]) continue;
        const double midpoint = 0.5 * (cuts[i - 1] + cuts[i]);
        int product = 1;
        for (const auto& line : lines) {
            int spin = line.spin;
            for (double tau : line.flips) {
                if (tau <= midpoint) spin = -spin;
            }
            product *= spin;
        }
        result += product * (cuts[i] - cuts[i - 1]);
    }
    return result;
}

void check_times(const std::vector<double>& actual, std::initializer_list<double> expected) {
    BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(), expected.end());
}

} // namespace

BOOST_AUTO_TEST_CASE(bounds_match_standard_searches) {
    // Exercise tiny, scan, and binary-search paths, including duplicate runs.
    for (int size : {0, 1, 2, 7, 8, 31, 64, 65, 129}) {
        std::vector<double> times;
        for (int i = 0; i < size; ++i) times.push_back(i / 3);
        for (double tau = -0.5; tau <= size + 0.5; tau += 0.5) {
            BOOST_CHECK(time_lower_bound(times.begin(), times.end(), tau)
                        == std::lower_bound(times.begin(), times.end(), tau));
            BOOST_CHECK(time_upper_bound(times.begin(), times.end(), tau)
                        == std::upper_bound(times.begin(), times.end(), tau));
        }
    }
}

BOOST_AUTO_TEST_CASE(exact_and_periodic_neighbor_searches) {
    std::vector<double> times{1., 2., 2., 4.};
    BOOST_CHECK_EQUAL(time_index(times, 2.), 1);
    BOOST_CHECK_THROW(time_index(times, 3.), std::runtime_error);
    BOOST_CHECK_THROW(time_index({}, 1.), std::runtime_error);
    BOOST_CHECK_EQUAL(next_event_time(times, 2.), 4.);
    BOOST_CHECK_EQUAL(previous_event_time(times, 2.), 1.);
    BOOST_CHECK_EQUAL(next_event_time(times, 4.), 1.);
    BOOST_CHECK_EQUAL(previous_event_time(times, 1.), 4.);
    BOOST_CHECK_EQUAL(next_event_time({}, 3.), 3.);
    BOOST_CHECK_EQUAL(previous_event_time({}, 3.), 3.);
    times = {2., 2.};
    BOOST_CHECK_EQUAL(next_event_time(times, 2.), 2.);
    BOOST_CHECK_EQUAL(previous_event_time(times, 2.), 2.);
}

BOOST_AUTO_TEST_CASE(history_insert_erase_and_move) {
    std::vector<double> times;
    insert_time_pair(times, 1., 3.);
    insert_time(times, 2.);
    insert_time_pair(times, 1., 3.);
    check_times(times, {1., 1., 2., 3., 3.});
    erase_time_pair(times, 1., 3.);
    check_times(times, {1., 2., 3.});
    erase_time(times, 2.);
    check_times(times, {1., 3.});
    BOOST_CHECK_THROW(erase_time(times, 2.), std::runtime_error);
    BOOST_CHECK_THROW(erase_time_pair(times, 1., 4.), std::runtime_error);
    BOOST_CHECK_THROW(insert_time_pair(times, 3., 1.), std::invalid_argument);
    BOOST_CHECK_THROW(erase_time_pair(times, 3., 1.), std::invalid_argument);
    check_times(times, {1., 3.});
    move_time(times, 0, 2., true);
    check_times(times, {2., 3.});
    move_time(times, 0, 4., false);
    check_times(times, {3., 4.});
    move_time(times, 1, 1., false);
    check_times(times, {1., 3.});
    times = {2.};
    move_time(times, 0, 4., false);
    check_times(times, {4.});
    times.clear();
    insert_time_pair(times, 2., 2.);
    check_times(times, {2., 2.});
    erase_time(times, 2.);
    check_times(times, {2.});
}

BOOST_AUTO_TEST_CASE(spin_and_integral_endpoint_conventions) {
    const std::vector<double> times{0., 1., 1., 3., 4.};
    const WorldlineView line{1, times};
    BOOST_CHECK_EQUAL(spin_at_time(line, 0.), -1);
    BOOST_CHECK_EQUAL(spin_at_time(line, 1.), -1);
    BOOST_CHECK_EQUAL(spin_at_time(line, 3.), 1);
    BOOST_CHECK_EQUAL(spin_at_time(line, 4.), -1);
    BOOST_CHECK_EQUAL(spin_integral(line, 0., 4.), -2.);
    BOOST_CHECK_EQUAL(spin_integral(line, 1., 3.), -2.);
    BOOST_CHECK_EQUAL(spin_integral(line, 3., 4.), 1.);
    BOOST_CHECK_EQUAL(spin_integral(line, 1., 1.), 0.);
    BOOST_CHECK_EQUAL(spin_integral(line, 3., 1.), 0.);
    BOOST_CHECK_EQUAL(spin_integral({-1, {}}, 0., 4.), -4.);
    BOOST_CHECK_EQUAL(spin_integral_no_inner_flips(line, 1., 3.), -2.);
    BOOST_CHECK_EQUAL(spin_integral_no_inner_flips(line, 1., 3., 3, 3.), -2.);
    BOOST_CHECK_EQUAL(spin_integral_no_inner_flips(line, 1., 3., 1, 1.), -2.);
    BOOST_CHECK_EQUAL(spin_integral_no_inner_flips(line, 3., 4., 3, 3.), 1.);
}

BOOST_AUTO_TEST_CASE(product_integrals_match_independent_segment_reference) {
    std::mt19937 rng(1927);
    for (int count : {0, 1, 2, 4, 9, 17}) {
        for (int trial = 0; trial < 30; ++trial) {
            std::vector<std::vector<double>> histories(count);
            std::vector<WorldlineView> lines;
            for (auto& history : histories) {
                const int n = rng() % 90;
                for (int i = 0; i < n; ++i) history.push_back((rng() % 65) * 0.25);
                std::sort(history.begin(), history.end());
                lines.push_back({rng() % 2 ? 1 : -1, history});
            }
            for (auto [left, right] : {std::pair{0., 16.}, {1., 7.}, {3.25, 3.5}, {4., 4.}, {7., 1.}}) {
                BOOST_CHECK_EQUAL(spin_product_integral(lines, left, right), reference_integral(lines, left, right));
                for (const auto& line : lines) {
                    BOOST_CHECK_EQUAL(spin_integral(line, left, right),
                                      reference_integral(std::span(&line, 1), left, right));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(site_storage_adapts_without_lattice_descriptors) {
    struct Site { int spin; std::vector<double> flips; };
    const std::array<Site, 2> sites{{{1, {1., 3.}}, {-1, {2., 3.}}}};
    const auto view = [](const Site& site) -> WorldlineView { return {site.spin, site.flips}; };
    // Product is -1 on [0,1], +1 on [1,2], -1 on [2,4]. The simultaneous
    // endpoint flips at 3 cancel in the bond product.
    BOOST_CHECK_EQUAL(spin_product_integral(sites, 0., 4., view), -2.);
    BOOST_CHECK_EQUAL(spin_integral(view(sites[0]), 0., 4.), 0.);
    const std::array<WorldlineView, 2> repeated{{view(sites[0]), view(sites[0])}};
    BOOST_CHECK_EQUAL(spin_product_integral(repeated, 0., 4.), 4.);
}

BOOST_AUTO_TEST_CASE(weighted_integrals) {
    const auto weight = [](double a, double b) {
        return triangular_weight_primitive(b, 4.) - triangular_weight_primitive(a, 4.);
    };
    BOOST_CHECK_EQUAL(weighted_spin_integral({1, {}}, 0., 4., weight), 4.);
    const std::vector<double> times{1., 3.};
    BOOST_CHECK_EQUAL(weighted_spin_integral({1, times}, 0., 4., weight), -2.);
    BOOST_CHECK_EQUAL(weighted_spin_integral({1, times}, 1., 3., weight), -3.);
    BOOST_CHECK_EQUAL(weighted_spin_integral({-1, times}, 1., 3., weight), 3.);
    BOOST_CHECK_EQUAL(weighted_spin_integral({1, times}, 1., 1., weight), 0.);
    // A different weight confirms the integration helper is not tied to an estimator.
    BOOST_CHECK_SMALL(weighted_spin_integral({1, times}, 0., 4.,
        [](double a, double b) { return (b*b*b - a*a*a) / 3.; }) - 4., 1.e-13);
}

BOOST_AUTO_TEST_CASE(origin_rotation_preserves_periodic_worldlines) {
    const std::vector<double> original{0.5, 1.5, 2.5, 3.5};
    const WorldlineView line{-1, original};
    std::vector<double> scratch;
    for (double cut : {0., 0.5, 1., 2., 3.75}) {
        auto times = original;
        const auto crossed = rotate_times(times, cut, 4., scratch);
        const int spin = crossed & 1 ? -line.spin : line.spin;
        BOOST_CHECK(std::is_sorted(times.begin(), times.end()));
        BOOST_CHECK_EQUAL(times.size(), original.size());
        BOOST_CHECK_SMALL(spin_integral({spin, times}, 0., 4.) - spin_integral(line, 0., 4.), 1.e-14);
        for (double tau : {0.125, 0.875, 1.125, 2.875, 3.875}) {
            double old_tau = tau + cut;
            if (old_tau >= 4.) old_tau -= 4.;
            BOOST_CHECK_EQUAL(spin_at_time({spin, times}, tau), spin_at_time(line, old_tau));
        }
    }
    std::vector<double> times{1., 1., 2., 3.};
    BOOST_CHECK_EQUAL(rotate_times(times, 1., 4., scratch), 0);
    check_times(times, {std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::epsilon(), 1., 2.});
    const auto buffer_size = scratch.size();
    times.clear();
    BOOST_CHECK_EQUAL(rotate_times(times, 1., 4., scratch), 0);
    BOOST_CHECK_EQUAL(scratch.size(), buffer_size);
}

} // namespace paratoric::detail
