// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#define BOOST_TEST_MODULE TestTupleParity

#include "lattice/lattice.hpp"
#include "lattice/lattice_geometry.hpp"
#include "lattice/tuple_parity.hpp"
#include "mcmc/extended_toric_code_qmc.hpp"
#include "paratoric/mcmc/extended_toric_code.hpp"
#include "rng/rng.hpp"
#include "tuple_parity_reference.hpp"

#include <boost/log/core.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace {

struct QuietLogging {
    bool enabled = boost::log::core::get()->get_logging_enabled();
    QuietLogging() { boost::log::core::get()->set_logging_enabled(false); }
    ~QuietLogging() { boost::log::core::get()->set_logging_enabled(enabled); }
};
BOOST_GLOBAL_FIXTURE(QuietLogging);

auto incidence(const paratoric::LatticeGeometry& geometry, char basis) {
    std::map<std::pair<int, int>, int> indices;
    for (std::size_t i = 0; i < geometry.edges.size(); ++i) {
        const auto& edge = geometry.edges[i];
        indices.emplace(std::minmax(edge.source, edge.target), static_cast<int>(i));
    }
    const auto count = basis == 'x' ? geometry.plaquettes.size() : geometry.vertices.size();
    std::vector<std::vector<int>> columns(count);
    if (basis == 'x') {
        for (std::size_t i = 0; i < count; ++i) {
            for (const auto& [u, v] : geometry.plaquettes[i]) {
                columns[i].push_back(indices.at(std::minmax(u, v)));
            }
        }
    } else {
        for (std::size_t i = 0; i < geometry.edges.size(); ++i) {
            const auto& edge = geometry.edges[i];
            columns[edge.source].push_back(static_cast<int>(i));
            columns[edge.target].push_back(static_cast<int>(i));
        }
    }
    return columns;
}

auto check_kernel(const paratoric::LatticeGeometry& geometry, char basis, std::size_t dimension) {
    const auto columns = incidence(geometry, basis);
    const auto kernel = paratoric::detail::tuple_parity_kernel(columns, geometry.edges.size());
    BOOST_REQUIRE_EQUAL(kernel.size(), dimension);
    std::set<int> last_indices;
    for (const auto& relation : kernel) {
        BOOST_REQUIRE(!relation.empty());
        BOOST_CHECK(std::is_sorted(relation.begin(), relation.end()));
        // Distinct last nonzero coordinates prove linear independence.
        BOOST_CHECK(last_indices.insert(relation.back()).second);
        std::vector<bool> edge_parity(geometry.edges.size(), false);
        for (const int tuple : relation) {
            for (const int edge : columns[tuple]) edge_parity[edge] = !edge_parity[edge];
        }
        BOOST_CHECK(std::none_of(edge_parity.begin(), edge_parity.end(), [](bool odd) { return odd; }));
    }
    return kernel;
}

double series_mean(const paratoric::Result& result, std::size_t observable = 0) {
    double sum = 0.;
    for (const auto& value : result.series[observable]) sum += std::get<double>(value);
    return sum / result.series[observable].size();
}

paratoric::Config commuting_config(char basis, const std::string& geometry, int size,
                                   const std::string& boundaries, double beta) {
    paratoric::Config config;
    config.lat_spec = {basis, geometry, size, beta, boundaries, 1};
    config.param_spec = {.mu = basis == 'z' ? 1. : 0., .h = 0.,
                         .J = basis == 'x' ? 1. : 0., .lmbda = 0.};
    config.sim_spec.N_thermalization = 10000;
    config.sim_spec.N_samples = 30000;
    config.sim_spec.N_between_samples = 30;
    config.sim_spec.N_resamples = 4;
    config.sim_spec.seed = 123;
    config.sim_spec.observables = {basis == 'x' ? "plaquette_z" : "star_x"};
    return config;
}

// For P commuting involutions with their only relation product_i B_i = I,
// Z is proportional to cosh(beta)^P + sinh(beta)^P. This includes both
// tuple parities; an even-only chain instead gives tanh(beta).
double commuting_mean(double beta, int count) {
    const double t = std::tanh(beta);
    return (t + std::pow(t, count - 1)) / (1. + std::pow(t, count));
}

} // namespace

BOOST_AUTO_TEST_CASE(kernel_covers_supported_geometries_and_boundaries) {
    for (const auto* type : {"square", "triangular", "honeycomb", "kagome", "cubic"}) {
        for (const auto* boundary : {"open", "periodic"}) {
            BOOST_TEST_CONTEXT(type << " / " << boundary) {
                const auto geometry = paratoric::make_lattice_geometry(type, 6, boundary);
                const bool periodic = std::string(boundary) == "periodic";
                const std::size_t x_dimension = std::string(type) == "cubic"
                    ? (periodic ? 218 : 125) : (periodic ? 1 : 0);
                check_kernel(geometry, 'x', x_dimension);
                check_kernel(geometry, 'z', 1);
            }
        }
    }
    paratoric::Lattice lat({'x', "square", 6, 1., "periodic", 1});
    const auto& relations = lat.get_tuple_parity_relations();
    BOOST_REQUIRE_EQUAL(relations.size(), 1);
    BOOST_CHECK_EQUAL(relations.front().size(), 36);
    BOOST_CHECK(&relations == &lat.get_tuple_parity_relations());
}

BOOST_AUTO_TEST_CASE(kernel_spans_local_cube_surfaces_and_global_relations) {
    const auto cube = paratoric::make_lattice_geometry("cubic", 2, "open");
    const auto single_relation = check_kernel(cube, 'x', 1);
    BOOST_REQUIRE_EQUAL(single_relation.front().size(), 6);

    const auto torus = paratoric::make_lattice_geometry("cubic", 3, "periodic");
    const auto kernel = check_kernel(torus, 'x', 29); // 26 independent cubes + 3 wrapping surfaces
    for (auto surface : torus.cube_plaquettes) {
        BOOST_REQUIRE_EQUAL(surface.size(), 6);
        std::sort(surface.begin(), surface.end());
        for (auto relation = kernel.rbegin(); relation != kernel.rend(); ++relation) {
            if (std::binary_search(surface.begin(), surface.end(), relation->back())) {
                std::vector<int> reduced;
                std::set_symmetric_difference(surface.begin(), surface.end(), relation->begin(), relation->end(),
                                             std::back_inserter(reduced));
                surface = std::move(reduced);
            }
        }
        BOOST_CHECK(surface.empty());
    }
}

BOOST_AUTO_TEST_CASE(conditional_poisson_counts_have_the_correct_parity_and_mean) {
    paratoric::rng::RNG rng(123);
    constexpr int draws = 20000;
    for (const double mean : {1e-12, 0.01, 0.7, 1.5, 40.}) {
        for (const bool odd : {false, true}) {
            double sum = 0.;
            for (int i = 0; i < draws; ++i) {
                const auto count = paratoric::rng::poisson_with_parity(rng, mean, odd);
                BOOST_REQUIRE_EQUAL(count % 2, odd ? 1 : 0);
                sum += count;
            }
            const double expected = odd ? mean / std::tanh(mean) : mean * std::tanh(mean);
            const double variance = std::max(0., mean * mean + expected - expected * expected);
            BOOST_CHECK_SMALL(sum / draws - expected, 6. * std::sqrt(variance / draws) + 1e-4);
        }
    }
    BOOST_CHECK_EQUAL(paratoric::rng::poisson_with_parity(rng, 0., false), 0);
    BOOST_CHECK_EQUAL(paratoric::rng::poisson_with_parity(rng, 0., true), 1);
}

BOOST_AUTO_TEST_CASE(single_event_parity_proposals_obey_detailed_balance) {
    for (const double mean : {1e-12, 0.1, 1.5, 10., 10000.}) {
        const paratoric::detail::PoissonParityProposal proposal(1., mean);
        BOOST_CHECK_EQUAL(proposal.birth_probability(0), 1.);
        for (std::size_t count = 0; count < 64; ++count) {
            const double log_forward = proposal.log_ratio(count, true);
            const double log_reverse = proposal.log_ratio(count + 1, false);
            BOOST_CHECK_SMALL(log_forward + log_reverse, 1e-14);
            // Adjacent Poisson count weights have ratio x/(k+1). Include both
            // proposal directions and Metropolis acceptance in the flux check.
            const double forward = proposal.birth_probability(count) * std::min(1., std::exp(log_forward));
            const double reverse = mean / (count + 1)
                * ((count + 1) / (mean + count + 1)) * std::min(1., std::exp(log_reverse));
            BOOST_CHECK_CLOSE_FRACTION(forward, reverse, 1e-12);
        }
    }
    const paratoric::detail::PoissonParityProposal tiny(1e-200, 1e-200);
    BOOST_CHECK(std::isfinite(tiny.log_ratio(0, true)));
    BOOST_CHECK_SMALL(tiny.log_ratio(0, true) + tiny.log_ratio(1, false), 1e-14);
}

namespace {

template<char Basis>
void compare_parity_transitions() {
    // Reuse the backends across geometries and coupling changes to exercise
    // scratch resizing and coupling changes as in hysteresis.
    auto actual_rng = std::make_shared<paratoric::rng::RNG>(321);
    auto reference_rng = std::make_shared<paratoric::rng::RNG>(321);
    paratoric::ExtendedToricCodeQMC<Basis> actual(actual_rng), reference(reference_rng);
    for (const auto* geometry : {"square", "triangular", "honeycomb", "kagome", "cubic"}) {
#ifndef NDEBUG
        // Debug construction requests unimplemented Kagome FM loop paths.
        // Release and sanitizer builds cover its parity transitions below.
        if (std::string(geometry) == "kagome") continue;
#endif
        for (const auto* boundary : {"open", "periodic"}) {
            const paratoric::LatSpec spec{Basis, geometry, 6, 3., boundary, -1};
            paratoric::Lattice lat(spec, actual_rng), reference_lat(spec, reference_rng);
            const int count = Basis == 'x' ? lat.get_plaquette_count() : lat.get_vertex_count();
            for (int tuple = 0; tuple < count; ++tuple) {
                const auto edges = Basis == 'x' ? lat.get_plaquette_edges(tuple) : lat.get_star_edges(tuple);
                const auto reference_edges = Basis == 'x' ? reference_lat.get_plaquette_edges(tuple)
                                                          : reference_lat.get_star_edges(tuple);
                // Include empty histories, coincident events, and events at zero.
                // Even tuple counts give a closed initial worldline on every edge.
                for (int event = 0; event < 2 * (tuple % 3); ++event) {
                    const double time = event / 2;
                    lat.insert_tuple_flip(tuple, edges, time);
                    reference_lat.insert_tuple_flip(tuple, reference_edges, time);
                }
            }
            lat.init_potential_energy();
            reference_lat.init_potential_energy();
            for (const double field : {0., 0.2, -0.2, 0.}) {
                // Zero-field moves intentionally leave inactive edge caches
                // stale. Rebuild them when changing the diagonal coupling.
                lat.init_potential_energy();
                reference_lat.init_potential_energy();
                for (const double coupling : {0., 0.7, 1.4}) {
                    const double h = Basis == 'x' ? field : 0.;
                    const double lmbda = Basis == 'z' ? field : 0.;
                    const double mu = Basis == 'z' ? coupling : 1.2;
                    const double J = Basis == 'x' ? coupling : 1.2;
                    double energy = paratoric::TupleParityTestAccess::energy(actual, lat, h, mu, J, lmbda);
                    double reference_energy = paratoric::TupleParityTestAccess::energy(reference, reference_lat, h, mu, J, lmbda);
                    for (int step = 0; step < 40; ++step) {
                        BOOST_TEST_CONTEXT(Basis << " / " << geometry << " / " << boundary
                                           << " / field=" << field << " / coupling=" << coupling
                                           << " / step=" << step) {
                            double ratio = -1., reference_ratio = -1.;
                            paratoric::TupleParityTestAccess::step(actual, lat, energy, ratio, spec.beta, h, mu, J, lmbda);
                            paratoric::TupleParityTestAccess::reference(reference, reference_lat, reference_energy,
                                                                      reference_ratio, spec.beta, h, mu, J, lmbda);
                            BOOST_REQUIRE_EQUAL(ratio, reference_ratio);
                            BOOST_REQUIRE_EQUAL(energy, reference_energy);
                            BOOST_REQUIRE_EQUAL(std::bit_cast<std::uint64_t>(ratio), std::bit_cast<std::uint64_t>(reference_ratio));
                            BOOST_REQUIRE_EQUAL(std::bit_cast<std::uint64_t>(energy), std::bit_cast<std::uint64_t>(reference_energy));
                            BOOST_REQUIRE(actual_rng->rng == reference_rng->rng);
                            const double integrated = paratoric::TupleParityTestAccess::energy(actual, lat, h, mu, J, lmbda);
                            BOOST_CHECK_SMALL(energy - integrated, 1e-9 * (1. + std::abs(integrated)));
                            bool tuple_histories_match = true;
                            for (int tuple = 0; tuple < count; ++tuple) {
                                const auto times = lat.get_tuple_spin_flips(tuple);
                                const auto reference_times = reference_lat.get_tuple_spin_flips(tuple);
                                tuple_histories_match &= std::equal(times.begin(), times.end(), reference_times.begin(), reference_times.end());
                            }
                            BOOST_REQUIRE(tuple_histories_match);
                            bool edge_histories_match = true, sorted = true, closed = true, spins_match = true, caches_match = true;
                            double max_cache_error = 0.;
                            // Each edge appears in a star; compare by its endpoint
                            // indices across independently constructed Boost graphs.
                            for (int vertex = 0; vertex < lat.get_vertex_count(); ++vertex) {
                                const auto edges = lat.get_star_edges(vertex);
                                const auto reference_edges = reference_lat.get_star_edges(vertex);
                                for (std::size_t i = 0; i < edges.size(); ++i) {
                                    const int events = lat.get_spin_flip_count(edges[i]);
                                    const bool same_count = events == reference_lat.get_spin_flip_count(reference_edges[i]);
                                    edge_histories_match &= same_count;
                                    closed &= events % 2 == 0;
                                    if (!same_count) continue;
                                    double previous = 0.;
                                    for (int event = 0; event < events; ++event) {
                                        const double time = lat.get_spin_flip_imag_time(edges[i], event);
                                        edge_histories_match &= time == reference_lat.get_spin_flip_imag_time(reference_edges[i], event);
                                        sorted &= time >= previous;
                                        previous = time;
                                    }
                                    spins_match &= lat.get_spin(edges[i]) == reference_lat.get_spin(reference_edges[i]);
                                    caches_match &= std::bit_cast<std::uint64_t>(lat.get_potential_edge_energy(edges[i]))
                                        == std::bit_cast<std::uint64_t>(reference_lat.get_potential_edge_energy(reference_edges[i]));
                                    if (field != 0.) {
                                        const double direct = lat.integrated_edge_energy(edges[i], 0., spec.beta);
                                        max_cache_error = std::max(max_cache_error, std::abs(lat.get_potential_edge_energy(edges[i]) - direct));
                                    }
                                }
                            }
                            BOOST_REQUIRE(edge_histories_match);
                            BOOST_CHECK(sorted);
                            BOOST_CHECK(closed);
                            BOOST_CHECK(spins_match);
                            BOOST_CHECK(caches_match);
                            BOOST_CHECK_SMALL(max_cache_error, 1e-10);
                        }
                    }
                }
            }
        }
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(optimized_transitions_match_the_original_in_both_bases) {
    compare_parity_transitions<'x'>();
    compare_parity_transitions<'z'>();
}

BOOST_AUTO_TEST_CASE(four_toggle_edges_match_the_original_transition) {
    const paratoric::LatSpec spec{'x', "cubic", 6, 3., "periodic", 1};
    auto actual_rng = std::make_shared<paratoric::rng::RNG>(941);
    auto reference_rng = std::make_shared<paratoric::rng::RNG>(941);
    paratoric::Lattice actual(spec, actual_rng), reference(spec, reference_rng);
    const auto& relations = actual.get_tuple_parity_relations();
    std::vector<int> relation;
    // Two cube surfaces meeting at an edge give four toggles there. This
    // identity relation exercises spill storage and the generic sorting path.
    for (std::size_t i = 1; i < relations.size(); ++i) {
        auto first = relations.front(), second = relations[i];
        std::sort(first.begin(), first.end());
        std::sort(second.begin(), second.end());
        relation.clear();
        std::set_symmetric_difference(first.begin(), first.end(), second.begin(), second.end(),
                                     std::back_inserter(relation));
        std::map<paratoric::Lattice::Edge, int> occurrences;
        for (const int tuple : relation)
            for (const auto& edge : actual.get_plaquette_edges(tuple)) ++occurrences[edge];
        if (std::any_of(occurrences.begin(), occurrences.end(), [](const auto& entry) { return entry.second == 4; })) break;
        relation.clear();
    }
    BOOST_REQUIRE(!relation.empty());
    for (auto* lattice : {&actual, &reference}) {
        // Replace the test lattice's relation before its edge cache is built.
        // The underlying lattice is mutable; select only this identity product
        // so every proposal exercises it instead of another kernel direction.
        auto& selected = const_cast<std::vector<std::vector<int>>&>(lattice->get_tuple_parity_relations());
        selected = {relation};
        for (const int tuple : relation) {
            for (const double time : {0., 0., 1.5, 1.5})
                lattice->insert_tuple_flip(tuple, lattice->get_plaquette_edges(tuple), time);
        }
        lattice->init_potential_energy();
    }
    const auto& edges = actual.get_tuple_parity_edges(0);
    BOOST_REQUIRE(std::any_of(edges.begin(), edges.end(), [](const auto& edge) { return edge.tuple_positions.size() == 4; }));
    paratoric::ExtendedToricCodeQMC<'x'> mc(actual_rng), reference_mc(reference_rng);
    double energy = paratoric::TupleParityTestAccess::energy(mc, actual, 0.2, 1.2, 1., 0.);
    double reference_energy = energy;
    for (int step = 0; step < 500; ++step) {
        double ratio = 0., reference_ratio = 0.;
        paratoric::TupleParityTestAccess::step(mc, actual, energy, ratio, spec.beta, 0.2, 1.2, 1., 0.);
        paratoric::TupleParityTestAccess::reference(reference_mc, reference, reference_energy,
                                                  reference_ratio, spec.beta, 0.2, 1.2, 1., 0.);
        BOOST_REQUIRE_EQUAL(ratio, reference_ratio);
        BOOST_REQUIRE_EQUAL(energy, reference_energy);
        BOOST_REQUIRE_EQUAL(std::bit_cast<std::uint64_t>(ratio), std::bit_cast<std::uint64_t>(reference_ratio));
        BOOST_REQUIRE_EQUAL(std::bit_cast<std::uint64_t>(energy), std::bit_cast<std::uint64_t>(reference_energy));
        BOOST_REQUIRE(actual_rng->rng == reference_rng->rng);
        for (const int tuple : relation) {
            const auto times = actual.get_tuple_spin_flips(tuple), reference_times = reference.get_tuple_spin_flips(tuple);
            BOOST_REQUIRE_EQUAL_COLLECTIONS(times.begin(), times.end(), reference_times.begin(), reference_times.end());
        }
        for (const auto& affected : edges) {
            const auto [u, v] = actual.vertices_of_edge(affected.edge);
            const auto reference_edge = reference.edge_in_between(u, v);
            const int count = actual.get_spin_flip_count(affected.edge);
            BOOST_REQUIRE_EQUAL(count, reference.get_spin_flip_count(reference_edge));
            for (int event = 0; event < count; ++event)
                BOOST_REQUIRE_EQUAL(actual.get_spin_flip_imag_time(affected.edge, event), reference.get_spin_flip_imag_time(reference_edge, event));
            const double direct = actual.integrated_edge_energy(affected.edge, 0., spec.beta);
            BOOST_CHECK_SMALL(actual.get_potential_edge_energy(affected.edge) - direct, 1e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(cubic_moves_use_local_surfaces_and_cached_edge_incidence) {
    for (const auto* boundary : {"periodic", "open"}) {
        paratoric::Lattice lat({'x', "cubic", 6, 1., boundary, 1});
        const auto& relations = lat.get_tuple_parity_relations();
        const bool periodic = std::string(boundary) == "periodic";
        BOOST_REQUIRE_EQUAL(relations.size(), periodic ? 218 : 125);
        std::size_t cubes = 0;
        for (std::size_t i = 0; i < relations.size(); ++i) {
            const auto& relation = relations[i];
            cubes += relation.size() == 6;
            BOOST_CHECK_LE(relation.size(), periodic ? 36 : 6);
            const auto& edges = lat.get_tuple_parity_edges(i);
            BOOST_CHECK(&edges == &lat.get_tuple_parity_edges(i));
            std::vector<int> occurrences(relation.size(), 0);
            for (const auto& affected : edges) {
                BOOST_CHECK_EQUAL(affected.tuple_positions.size() % 2, 0);
                for (const int position : affected.tuple_positions) {
                    BOOST_REQUIRE_LT(position, relation.size());
                    const auto tuple_edges = lat.get_plaquette_edges(relation[position]);
                    BOOST_CHECK(std::find(tuple_edges.begin(), tuple_edges.end(), affected.edge) != tuple_edges.end());
                    ++occurrences[position];
                }
            }
            BOOST_CHECK(std::all_of(occurrences.begin(), occurrences.end(), [](int count) { return count == 4; }));
        }
        BOOST_CHECK_EQUAL(cubes, periodic ? 215 : 125);
    }
}

BOOST_AUTO_TEST_CASE(zero_field_workflows_keep_active_caches_consistent) {
    for (const char basis : {'x', 'z'}) {
        for (const auto* geometry : {"square", "cubic"}) {
            auto config = commuting_config(basis, geometry, 6, "open", 2.);
            config.param_spec = {.mu = 1.2, .h = basis == 'x' ? 0.2 : 0.,
                                 .J = 1.2, .lmbda = basis == 'z' ? 0.2 : 0.};
            config.sim_spec.N_samples = 20;
            config.sim_spec.N_thermalization = 500;
            config.sim_spec.N_between_samples = 100;
            // Workflow end checks compare the running action with direct history
            // integration. Open cubic plaquettes include local parity relations.
            BOOST_CHECK_NO_THROW(paratoric::ExtendedToricCode::get_sample(config));
            BOOST_CHECK_NO_THROW(paratoric::ExtendedToricCode::get_thermalization(config));
            config.param_spec.h_hys = basis == 'x' ? std::vector{0.1, 0.2} : std::vector{0., 0.};
            config.param_spec.lmbda_hys = basis == 'z' ? std::vector{0.1, 0.2} : std::vector{0., 0.};
            BOOST_CHECK_NO_THROW(paratoric::ExtendedToricCode::get_hysteresis(config));
        }
    }
}

// Long statistical comparisons run in optimized builds. Debug checks above
// exercise topology, conditional counts and all workflows with bounded counts.
// Existing Debug construction also rejects x-basis extents below six while
// preparing unused FM paths (audit finding 8).
#ifdef NDEBUG
BOOST_AUTO_TEST_CASE(commuting_square_torus_matches_full_trace_in_both_bases) {
    for (const char basis : {'x', 'z'}) {
        const auto config = commuting_config(basis, "square", 6, "periodic", 2.);
        const auto result = paratoric::ExtendedToricCode::get_sample(config);
        const double measured = series_mean(result);
        BOOST_TEST_MESSAGE("basis=" << basis << " torus mean=" << measured);
        BOOST_CHECK_SMALL(measured - commuting_mean(2., 36), 0.006);
    }
}

BOOST_AUTO_TEST_CASE(audit_square_torus_and_single_cube_match_exact_traces) {
    for (const auto* geometry : {"square", "cubic"}) {
        const bool cube = std::string(geometry) == "cubic";
        const auto config = commuting_config('x', geometry, cube ? 2 : 3,
                                             cube ? "open" : "periodic", 1.5);
        const auto result = paratoric::ExtendedToricCode::get_sample(config);
        const double measured = series_mean(result);
        const double expected = commuting_mean(1.5, cube ? 6 : 9);
        BOOST_TEST_MESSAGE(geometry << " mean=" << measured << " exact=" << expected);
        BOOST_CHECK_SMALL(measured - expected, 0.008);
    }
}
BOOST_AUTO_TEST_CASE(diagonal_action_matches_the_other_basis) {
    // The reference basis has a positive single-spin field and uses the original
    // combination update. Compare both dual zero-field cases, exercising cache
    // updates on accepted relations and history restoration on rejections.
    for (const char parity_basis : {'x', 'z'}) {
        auto config = commuting_config(parity_basis, "square", 6, "periodic", 2.);
        config.param_spec = {.mu = 1.2, .h = parity_basis == 'x' ? 0.2 : 0.,
                             .J = 1.2, .lmbda = parity_basis == 'z' ? 0.2 : 0.};
        config.sim_spec.N_between_samples = 100;
        config.sim_spec.observables = {parity_basis == 'x' ? "plaquette_z" : "star_x",
                                      parity_basis == 'x' ? "sigma_x" : "sigma_z"};
        const auto zero_field = paratoric::ExtendedToricCode::get_sample(config);
        config.lat_spec.basis = parity_basis == 'x' ? 'z' : 'x';
        const auto reference = paratoric::ExtendedToricCode::get_sample(config);
        for (std::size_t i = 0; i < 2; ++i) {
            BOOST_CHECK_SMALL(series_mean(zero_field, i) - series_mean(reference, i), 0.025);
        }
    }
}
#endif
