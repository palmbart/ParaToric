// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#include "lattice/lattice.hpp"
#include "lattice/lattice_geometry.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>

#include <algorithm> 
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <span>
#include <sstream>
#include <stack>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

#define UNUSED(expr) do { (void)(expr); } while (0)

namespace paratoric {

void Lattice::SnapshotSpoolState::reset() {
    if (stream.is_open()) {
        stream.close();
    }
    if (!path.empty()) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        path.clear();
    }
    edge_count = 0;
    sample_count = 0;
}

Lattice::SnapshotSpoolState::~SnapshotSpoolState() {
    reset();
}

void Lattice::check_input_validity() const {
    if (BASIS != 'x' && BASIS != 'z') {
        throw std::invalid_argument("Basis has to be 'x' or 'z'.");
    }

    if (LATTICE_TYPE != "square"
        && LATTICE_TYPE != "cubic"
        && LATTICE_TYPE != "honeycomb"
        && LATTICE_TYPE != "triangular"
        && LATTICE_TYPE != "kagome") {
        throw std::invalid_argument(std::format(
            "Lattice type \"{}\" is not supported.", LATTICE_TYPE
        ));
    }

    if (SYSTEM_SIZE <= 0) {
        throw std::invalid_argument("System size has to be strictly positive.");
    }

    if (DEFAULT_SPIN != -1 && DEFAULT_SPIN != 1) {
        throw std::invalid_argument("Default spin has to be 1 or -1.");
    }

    if (!(BETA > 0)) {
        throw std::invalid_argument("Beta has to be strictly positive.");
    }

    if (BOUNDARIES != "periodic" && BOUNDARIES != "open") {
        throw std::invalid_argument("Boundaries can be \"periodic\" or \"open\".");
    }
}

void Lattice::set_rng(std::shared_ptr<RNG>& rng_inp) {
    rng = rng_inp;
}

std::shared_ptr<paratoric::rng::RNG> Lattice::get_rng() {
    return rng;
}

void Lattice::init_lattice_graph() {
    auto geometry = make_lattice_geometry(LATTICE_TYPE, SYSTEM_SIZE, BOUNDARIES);
    std::tie(half_path_vector, full_path_vector) = make_fredenhagen_marcu_paths(
        geometry, LATTICE_TYPE, SYSTEM_SIZE, BOUNDARIES, BASIS
    );
    LATTICE_DIMENSIONALITY = geometry.dimensionality;
    g = LatticeGraph(geometry.vertices.size());
    for (std::size_t v = 0; v < geometry.vertices.size(); ++v) {
        const auto& position = geometry.vertices[v];
        g[v].x = position[0];
        g[v].y = position[1];
        g[v].z = position[2];
    }
    edge_vector.reserve(geometry.edges.size());
    for (const auto& edge : geometry.edges) {
        EdgeData data(edge.orientation);
        data.spin = DEFAULT_SPIN;
        data.source_vertex = edge.source;
        data.target_vertex = edge.target;
        boost::add_edge(edge.source, edge.target, data, g);
        edge_vector.emplace_back(edge.source, edge.target);
    }

    plaquette_vector = std::move(geometry.plaquettes);
    cube_vector = std::move(geometry.cubes);
    plaquette_part_of_cube_lookup = std::move(geometry.plaquette_cubes);
    cube_has_plaquettes_lookup = std::move(geometry.cube_plaquettes);
    MAX_COORDINATES = std::move(geometry.max_coordinates);
    MAX_PLAQUETTE_COORDINATES = std::move(geometry.max_plaquette_coordinates);
    for (const auto& position : geometry.plaquette_coordinates) {
        plaquette_x_vector.push_back(position[0]);
        plaquette_y_vector.push_back(position[1]);
        if (LATTICE_DIMENSIONALITY == 3) plaquette_z_vector.push_back(position[2]);
    }
    for (const auto& position : geometry.cube_coordinates) {
        cube_x_vector.push_back(position[0]);
        cube_y_vector.push_back(position[1]);
        cube_z_vector.push_back(position[2]);
    }
    plaquette_flip_vector.resize(plaquette_vector.size());
    integrated_plaquette_energy_vector.resize(plaquette_vector.size(), 0.);
    edge_dist = std::uniform_int_distribution<int>(0, get_edge_count() - 1);
    vertex_dist = std::uniform_int_distribution<int>(0, get_vertex_count() - 1);
    plaquette_dist = std::uniform_int_distribution<int>(0, static_cast<int>(plaquette_vector.size()) - 1);
}

// Resolve geometry once so proposals can use cached adjacency and descriptors.
void Lattice::build_caches_() {
    egde_cache_.clear();
    for (auto e : boost::make_iterator_range(boost::edges(g))) {
        g[e].part_of_plaquette_lookup.clear();
        egde_cache_.emplace_back(e);
    }

    // Plaquette -> edges in construction order and sorted unique vertices.
    plaquette_edges_cache_.clear();
    plaquette_edges_cache_.resize(plaquette_vector.size());
    plaquette_vertices_cache_.clear();
    plaquette_vertices_cache_.resize(plaquette_vector.size());

    for (size_t p = 0; p < plaquette_vector.size(); ++p) {
        const auto vpairs = get_plaquette_vertex_pairs(static_cast<int>(p));
        auto& pedges = plaquette_edges_cache_[p];
        auto& pverts = plaquette_vertices_cache_[p];
        pedges.clear();
        pedges.reserve(vpairs.size());
        pverts.clear();
        pverts.reserve(vpairs.size() * 2);

        for (const auto& pr : vpairs) {
            auto uv = boost::edge(pr.first, pr.second, g);
            if (!uv.second) {
                // Retry the vertex pair in reverse order.
                uv = boost::edge(pr.second, pr.first, g);
            }
            if (!uv.second) {
                // Reject missing geometry before caching an invalid descriptor.
                throw std::runtime_error("build_caches_: missing edge between plaquette vertices");
            }
            pedges.emplace_back(uv.first);
            pverts.emplace_back(pr.first);
            pverts.emplace_back(pr.second);
            g[uv.first].part_of_plaquette_lookup.emplace_back(static_cast<int>(p));
        }

        std::sort(pverts.begin(), pverts.end());
        pverts.erase(std::unique(pverts.begin(), pverts.end()), pverts.end());
    }

    // Star center -> incident edges in graph iteration order.
    auto vindex = boost::get(boost::vertex_index, g);
    const auto V = static_cast<size_t>(boost::num_vertices(g));

    star_edges_cache_.clear();
    star_edges_cache_.resize(V);

    for (auto v : boost::make_iterator_range(boost::vertices(g))) {
        const size_t idx = static_cast<size_t>(get(vindex, v));
        auto& lst = star_edges_cache_[idx];
        lst.clear();
        lst.reserve(static_cast<size_t>(boost::out_degree(v, g)));

        for (auto e : boost::make_iterator_range(boost::out_edges(v, g))) {
            lst.emplace_back(e);
        }
    }

    // Star center -> sorted unique adjacent plaquette indices.
    star_plaquettes_cache_.clear();
    star_plaquettes_cache_.resize(V);
    for (size_t v = 0; v < V; ++v) {
        auto& plist = star_plaquettes_cache_[v];
        plist.clear();
        for (const auto& e : star_edges_cache_[v]) {
            const auto& p_lookup = g[e].part_of_plaquette_lookup;
            plist.insert(plist.end(), p_lookup.begin(), p_lookup.end());
        }
        std::sort(plist.begin(), plist.end());
        plist.erase(std::unique(plist.begin(), plist.end()), plist.end());
    }
}

int Lattice::get_anyon_count() {
    int anyon_count = 0;
    if (BASIS == 'x') {
        for (const auto& v : boost::make_iterator_range(boost::vertices(g))) {
            if (get_vertex_nn_spins_prod(v) == -1) {
                anyon_count += 1;
            } 
        } 
    } else {
        int sum = 0;
        const int P = get_plaquette_count();
        for (int p = 0; p < P; ++p) {
            const auto& pedges = get_plaquette_edges(p);
            const int prod = get_tuple_prod(pedges);  
            sum += (prod == -1);                   
        }
        anyon_count = sum;
    }
    return anyon_count;
}

int Lattice::get_spin_flip_index(const Edge& edg, double tau) {
    return static_cast<int>(detail::time_index(g[edg].spin_flips, tau));
}

int Lattice::get_single_spin_flip_index(const Edge& edg, double tau) {
    return static_cast<int>(detail::time_index(g[edg].single_spin_flips, tau));
}

std::span<const double> Lattice::get_tuple_spin_flips(int t_index) {
    if (BASIS == 'x') {
        const auto& v = plaquette_flip_vector[t_index];
        return {v.data(), v.size()};
    } else {
        const auto& v = g[t_index].star_flips;
        return {v.data(), v.size()};
    }
}

std::span<const double> Lattice::get_single_spin_flips(const Edge& edg) {
    const auto& v = g[edg].single_spin_flips;
    return {v.data(), v.size()};
}

void Lattice::delete_double_single_spin_flip(const Edge& edg, 
    double imag_time_single_spin_flip, 
    double imag_time_next_single_spin_flip
) {
    detail::erase_time_pair(g[edg].spin_flips, imag_time_single_spin_flip, imag_time_next_single_spin_flip);
    detail::erase_time_pair(g[edg].single_spin_flips, imag_time_single_spin_flip, imag_time_next_single_spin_flip);
}

void Lattice::delete_single_spin_flip(const Edge& edg, int spin_flip_index) {
    auto& spin_flips = g[edg].spin_flips;
    const auto imag_time = spin_flips[spin_flip_index];
    spin_flips.erase(spin_flips.begin() + spin_flip_index);
    detail::erase_time(g[edg].single_spin_flips, imag_time);
}

void Lattice::delete_double_tuple_flip(
    int tuple_index, 
    std::span<const Edge> tuple_edges, 
    double imag_time_tuple_flip, 
    double imag_time_next_tuple_flip
) {
    if (imag_time_next_tuple_flip < imag_time_tuple_flip) [[unlikely]] {
        throw std::invalid_argument("delete_double_tuple_flip: imag_time_next_tuple_flip has to be larger than imag_time_tuple_flip.");
    }
    for (const Edge& edg : tuple_edges) {
        detail::erase_time_pair(g[edg].spin_flips, imag_time_tuple_flip, imag_time_next_tuple_flip);
    }
    auto& tuple_flips = BASIS == 'x' ? plaquette_flip_vector[tuple_index] : g[tuple_index].star_flips;
    detail::erase_time_pair(tuple_flips, imag_time_tuple_flip, imag_time_next_tuple_flip);
}

void Lattice::delete_tuple_flip(int tuple_index, std::span<const Edge> tuple_edges, double imag_time_tuple_flip) {
    for (const Edge& edg : tuple_edges) {
        detail::erase_time(g[edg].spin_flips, imag_time_tuple_flip);
    }
    auto& tuple_flips = BASIS == 'x' ? plaquette_flip_vector[tuple_index] : g[tuple_index].star_flips;
    detail::erase_time(tuple_flips, imag_time_tuple_flip);
}

bool Lattice::check_tuple_flip_present_tuple(std::span<const Edge> tuple_edges, double tau) {
    for (const Edge& edg : tuple_edges) {
        const auto& spin_flips = g[edg].spin_flips;
        if (!std::binary_search(spin_flips.begin(), spin_flips.end(), tau))
            return false;
    }
    return true;
}

bool Lattice::check_plaquette_flip_at_edge(const Edge& edg, double tau) {
    const auto& plist = g[edg].part_of_plaquette_lookup; // usually size ≤ 2
    for (int p_index : plist) {
        const auto& pedges = get_plaquette_edges(p_index);
        if (check_tuple_flip_present_tuple(pedges, tau)) {
            return true; // early exit
        }
    }
    return false;
}

bool Lattice::check_star_flip_at_edge(const Edge& edg, double tau) {
    auto [source_v, target_v] = vertices_of_edge(edg);
    std::array<int, 2> centers = {source_v, target_v};
    return std::any_of(centers.begin(), centers.end(), [this, tau](int center_index) {
        const auto& sedges = get_star_edges(center_index);
        return check_tuple_flip_present_tuple(sedges, tau);
    });
}

bool Lattice::check_spin_flips_present_tuple(std::span<const Edge> tuple_edges, double tau) {
    return std::any_of(tuple_edges.begin(), tuple_edges.end(), [this, tau](const Edge& edg) {
        return std::binary_search(g[edg].spin_flips.begin(), g[edg].spin_flips.end(), tau);
    });
}

double Lattice::flip_next_imag_time(const Edge& edg, double tau) {
    return detail::next_event_time(g[edg].spin_flips, tau);
}

std::vector<double> Lattice::flip_next_imag_times_tuple(std::span<const Edge> tuple_edges, double tau) {
    auto small = flip_next_imag_times_tuple_small(tuple_edges, tau);
    return {small.begin(), small.end()};
}

Lattice::SmallEnergyVector Lattice::flip_next_imag_times_tuple_small(std::span<const Edge> tuple_edges, double tau) {
    SmallEnergyVector imag_times;
    imag_times.reserve(tuple_edges.size());
    for (const Edge& edg : tuple_edges) {
        imag_times.emplace_back(flip_next_imag_time(edg, tau));
    }
    return imag_times;
}

double Lattice::flip_prev_imag_time(const Edge& edg, double tau) {
    return detail::previous_event_time(g[edg].spin_flips, tau);
}

std::vector<double> Lattice::flip_prev_imag_times_tuple(std::span<const Edge> tuple_edges, double tau) {
    auto small = flip_prev_imag_times_tuple_small(tuple_edges, tau);
    return {small.begin(), small.end()};
}

Lattice::SmallEnergyVector Lattice::flip_prev_imag_times_tuple_small(std::span<const Edge> tuple_edges, double tau) {
    SmallEnergyVector imag_times;
    imag_times.reserve(tuple_edges.size());
    for (const Edge& edg : tuple_edges) {
        imag_times.emplace_back(flip_prev_imag_time(edg, tau));
    }
    return imag_times;
}

std::pair<double, double> Lattice::tuple_flip_window(
    std::span<const Edge> tuple_edges,
    double tau,
    SmallIndexVector* flip_indices
) {
    bool initialized = false;
    double tau_left = tau;
    double tau_right = tau;
    if (flip_indices != nullptr) {
        flip_indices->clear();
        flip_indices->reserve(tuple_edges.size());
    }

    for (const Edge& edg : tuple_edges) {
        const auto& spin_flips = g[edg].spin_flips;
        double next = tau;
        double prev = tau;

        if (!spin_flips.empty()) [[likely]] {
            const auto lower = detail::time_lower_bound(spin_flips.begin(), spin_flips.end(), tau);
            // The lower bound is already at the event (normally unique).
            // Skip equal-time events without another search of the entire tail.
            auto upper = lower;
            while (upper != spin_flips.end() && *upper == tau) ++upper;

            if (flip_indices != nullptr) {
                flip_indices->emplace_back(static_cast<int>(lower - spin_flips.begin()));
            }

            next = (upper != spin_flips.end()) ? *upper : spin_flips.front();
            if (next < tau) {
                next += BETA;
            }

            prev = (lower != spin_flips.begin()) ? *(lower - 1) : spin_flips.back();
            if (prev > tau) {
                prev -= BETA;
            }
        }

        if (!initialized) [[unlikely]] {
            tau_left = prev;
            tau_right = next;
            initialized = true;
        } else {
            if (next < tau_right) {
                tau_right = next;
            }
            if (prev > tau_left) {
                tau_left = prev;
            }
        }
    }

    return {modulo(tau_left, BETA), modulo(tau_right, BETA)};
}

void Lattice::insert_double_single_spin_flip(const Edge& edg, double tau_left, double tau_right) {
    detail::insert_time_pair(g[edg].spin_flips, tau_left, tau_right);
    detail::insert_time_pair(g[edg].single_spin_flips, tau_left, tau_right);
}

void Lattice::insert_single_spin_flip(const Edge& edg, double tau) {
    detail::insert_time(g[edg].spin_flips, tau);
    detail::insert_time(g[edg].single_spin_flips, tau);
}

void Lattice::insert_double_tuple_flip(
    int tuple_index, 
    std::span<const Edge> tuple_edges, 
    double tau_left, 
    double tau_right) {
    if (tau_right < tau_left) [[unlikely]] {
        throw std::invalid_argument("insert_double_tuple_flip: tau_right has to be larger than tau_left.");
    }
    for (const Edge& edg : tuple_edges) {
        detail::insert_time_pair(g[edg].spin_flips, tau_left, tau_right);
    }
    auto& tuple_flips = BASIS == 'x' ? plaquette_flip_vector[tuple_index] : g[tuple_index].star_flips;
    detail::insert_time_pair(tuple_flips, tau_left, tau_right);
}

void Lattice::insert_tuple_flip(int tuple_index, std::span<const Edge> tuple_edges, double tau) {
    for (const Edge& edg : tuple_edges) {
        detail::insert_time(g[edg].spin_flips, tau);
    }
    auto& tuple_flips = BASIS == 'x' ? plaquette_flip_vector[tuple_index] : g[tuple_index].star_flips;
    detail::insert_time(tuple_flips, tau);
}

void Lattice::move_spin_flip(
    const Edge& edg, 
    int spin_flip_index, 
    double tau_new, 
    bool no_move_over_beta
) {
    auto& spin_flips = g[edg].spin_flips;
    auto& single_spin_flips = g[edg].single_spin_flips;
    const auto single_index = no_move_over_beta
        ? detail::time_index(single_spin_flips, spin_flips[spin_flip_index])
        : (spin_flip_index == 0 ? 0 : single_spin_flips.size() - 1);
    detail::move_time(single_spin_flips, single_index, tau_new, no_move_over_beta);
    detail::move_time(spin_flips, spin_flip_index, tau_new, no_move_over_beta);
}

void Lattice::move_tuple_flip(
    int tuple_index, 
    std::span<const Edge> tuple_edges, 
    double tau_old, 
    double tau_new, 
    bool no_move_over_beta,
    std::span<const int> edge_flip_indices,
    int tuple_flip_index
) {
    if (!edge_flip_indices.empty() && edge_flip_indices.size() != tuple_edges.size()) [[unlikely]] {
        throw std::invalid_argument(
            "move_tuple_flip: edge_flip_indices and tuple_edges must have equal sizes."
        );
    }
    for (size_t i = 0; i < tuple_edges.size(); ++i) {
        auto& flips = g[tuple_edges[i]].spin_flips;
        const auto index = edge_flip_indices.empty()
            ? detail::time_index(flips, tau_old)
            : static_cast<std::size_t>(edge_flip_indices[i]);
        detail::move_time(flips, index, tau_new, no_move_over_beta);
    }
    auto& tuple_flips = BASIS == 'x' ? plaquette_flip_vector[tuple_index] : g[tuple_index].star_flips;
    const auto index = tuple_flip_index >= 0
        ? static_cast<std::size_t>(tuple_flip_index)
        : detail::time_index(tuple_flips, tau_old);
    detail::move_time(tuple_flips, index, tau_new, no_move_over_beta);
}

void Lattice::print_spins() {
    for (const auto& vertex : boost::make_iterator_range(boost::vertices(g))) {
        for (const auto& neighbor : boost::make_iterator_range(boost::adjacent_vertices(vertex, g))) { 
            if (get_spin(edge_in_between(vertex, neighbor)) == -1) {
                std::println("Edge between {0} and {1} has spin {2}", vertex, neighbor, get_spin(edge_in_between(vertex, neighbor)));
            }
        }
    }
}

void Lattice::print_spin_flip_imag_times(const Edge& edg) {
    std::println("Total number of spin flips: {}", get_spin_flip_count(edg));
    std::println("Spin at imaginary time 0/beta: {}", get_spin(edg));
    std::println("Spin flip imaginary times: {}", g[edg].spin_flips);
}

void Lattice::print_tuple_flip_imag_times(std::span<const Edge> tuple_edges) {
    std::println("Printing spin flips of tuple\n");
    for (const Edge& edg : tuple_edges) {
        const auto [source_v, target_v] = vertices_of_edge(edg);
        std::println("Edge between {0} and {1}:", source_v, target_v);
        print_spin_flip_imag_times(edg);
    }
}

void Lattice::flip_spin(const Edge& edg) {
    g[edg].spin *= -1;
}

void Lattice::flip_star(int v) {
    for (const auto& edg : boost::make_iterator_range(boost::out_edges(v, g))) { 
        flip_spin(edg);
    }
}

std::tuple<Lattice::Edge, int, int> Lattice::get_random_edge() {
    const size_t edg_index = static_cast<size_t>(
        paratoric::rng::uniform_index(*rng, static_cast<std::uint64_t>(egde_cache_.size()))
    );
    auto [u, v] = edge_vector[edg_index];
    auto edg = egde_cache_[edg_index];
    return {edg, u, v};
}

Lattice::Edge Lattice::get_random_edge_descriptor() {
    const size_t edg_index = static_cast<size_t>(
        paratoric::rng::uniform_index(*rng, static_cast<std::uint64_t>(egde_cache_.size()))
    );
    return egde_cache_[edg_index];
}

int Lattice::get_random_vertex() {
    return static_cast<int>(
        paratoric::rng::uniform_index(*rng, static_cast<std::uint64_t>(get_vertex_count()))
    );
}

int Lattice::get_random_plaquette_index() {
    return static_cast<int>(
        paratoric::rng::uniform_index(*rng, static_cast<std::uint64_t>(plaquette_vector.size()))
    );
}

std::vector<std::pair<int,int>> Lattice::get_plaquette_vertex_pairs(int p_index) {
    return plaquette_vector[p_index];
}

std::vector<int> Lattice::get_cube_vertices(int c_index) {
    return cube_vector[c_index];
}

int Lattice::get_tuple_sum(std::span<const Edge> tuple_edges) {
    const int result = std::accumulate(
        tuple_edges.begin(), 
        tuple_edges.end(), 
        0, 
        [&](int lhs, const Edge& rhs) {return lhs + get_spin(rhs);}
    );
    return result;
}

int Lattice::get_tuple_prod(std::span<const Edge> tuple_edges) {
    const int result = std::accumulate(
        tuple_edges.begin(), 
        tuple_edges.end(), 
        1, 
        [&](int lhs, const Edge& rhs) {return lhs * get_spin(rhs);}
    );
    return result;
}

void Lattice::flip_tuple(std::span<const Edge> tuple_edges) {
    for (const Edge& edg : tuple_edges) {
        flip_spin(edg);
    }
}

int Lattice::get_vertex_nn_spins_prod(int v) {
    int prod = 1;
    for (auto e : boost::make_iterator_range(boost::out_edges(v, g))) {
        prod *= get_spin(e);
    }
    return prod;
}

double Lattice::integrated_tuple_energy_diff(
    std::span<const Edge> tuple_edges, 
    double imag_time_1, 
    double imag_time_2
) {
    return -2. * integrated_tuple_energy(tuple_edges, imag_time_1, imag_time_2);
}

std::tuple<double, Lattice::SmallIndexVector, Lattice::SmallEnergyVector> 
Lattice::integrated_star_energy_diff(
    const Edge& edg, 
    double imag_time_1, 
    double imag_time_2, 
    bool total_cache
) {
    if (imag_time_1 == imag_time_2) [[unlikely]] {
        throw std::invalid_argument("integrated_star_energy_diff: Time interval has to be non-zero.");
    } 

    SmallIndexVector star_centers;
    SmallEnergyVector star_potential_energy_diffs;

    auto [source_v, target_v] = vertices_of_edge(edg);
    double energy = 0.;
    double source_energy = 0.;
    double target_energy = 0.;
    if (total_cache) {
        source_energy = -2 * get_potential_star_energy(source_v);
        target_energy = -2 * get_potential_star_energy(target_v);
    } else {
        const auto& sedges_source = get_star_edges(source_v);
        source_energy = (BASIS == 'x')
            ? integrated_tuple_energy_diff_single_flips(sedges_source, imag_time_1, imag_time_2)
            : integrated_tuple_energy_diff(sedges_source, imag_time_1, imag_time_2);
        const auto& sedges_target = get_star_edges(target_v);
        target_energy = (BASIS == 'x')
            ? integrated_tuple_energy_diff_single_flips(sedges_target, imag_time_1, imag_time_2)
            : integrated_tuple_energy_diff(sedges_target, imag_time_1, imag_time_2);
    }
    star_centers.emplace_back(source_v);
    star_potential_energy_diffs.emplace_back(source_energy);
    energy += source_energy;
    star_centers.emplace_back(target_v);
    star_potential_energy_diffs.emplace_back(target_energy);
    energy += target_energy;

    return {energy, std::move(star_centers), std::move(star_potential_energy_diffs)};
}

std::tuple<double, Lattice::SmallIndexVector, Lattice::SmallEnergyVector> 
Lattice::integrated_plaquette_energy_diff(
    const Edge& edg, 
    double imag_time_1, 
    double imag_time_2, 
    bool total_cache
) {
    if (imag_time_1 == imag_time_2) [[unlikely]] {
        throw std::invalid_argument("integrated_plaquette_energy_diff: Time interval has to be non-zero.");
    } 
    SmallIndexVector plaquette_indices;
    SmallEnergyVector plaquette_potential_energy_diffs;
    plaquette_indices.reserve(g[edg].part_of_plaquette_lookup.size());
    plaquette_potential_energy_diffs.reserve(g[edg].part_of_plaquette_lookup.size());

    double energy = 0.;
    for (int p_index : g[edg].part_of_plaquette_lookup ) {
        double energy_p = 0.;
        if (total_cache) {
            energy_p = -2*get_potential_plaquette_energy(p_index);
        } else {
            const auto& pedges = get_plaquette_edges(p_index);
            energy_p = (BASIS == 'z')
                ? integrated_tuple_energy_diff_single_flips(pedges, imag_time_1, imag_time_2)
                : integrated_tuple_energy_diff(pedges, imag_time_1, imag_time_2);
        }
        plaquette_indices.emplace_back(p_index);
        plaquette_potential_energy_diffs.emplace_back(energy_p);
        energy += energy_p;
    }
    return {energy, std::move(plaquette_indices), std::move(plaquette_potential_energy_diffs)};
}

double Lattice::total_integrated_star_energy() {
    double total_integrated_star_energy = 0.;
    for (size_t star_center = 0; star_center < (size_t)get_vertex_count(); ++star_center) {
        const auto& sedges = get_star_edges(star_center);
        total_integrated_star_energy += (BASIS == 'x')
            ? integrated_tuple_energy_single_flips(sedges, 0, BETA)
            : integrated_tuple_energy(sedges, 0, BETA);
    }
    return total_integrated_star_energy;
}

double Lattice::total_integrated_plaquette_energy() {
    double total_integrated_plaquette_energy = 0.;
    for (size_t plaquette_index = 0; plaquette_index < plaquette_vector.size(); ++plaquette_index) {
        const auto& pedges = get_plaquette_edges(plaquette_index);
        total_integrated_plaquette_energy += (BASIS == 'z')
            ? integrated_tuple_energy_single_flips(pedges, 0, BETA)
            : integrated_tuple_energy(pedges, 0, BETA);
    }
    return total_integrated_plaquette_energy;
}

namespace {
    // Insert (time, tag) into a short sorted local schedule.
    inline void insert_sorted_by_time(std::vector<std::pair<double,int>>& v,
                                      double t, int tag) {
        auto it = v.begin();
        // Preserve insertion order for equal times.
        for (; it != v.end() && it->first <= t; ++it) {}
        v.insert(it, {t, tag});
    }
} // namespace

[[gnu::hot]]
std::tuple<double, Lattice::SmallIndexVector, Lattice::SmallEnergyVector> 
Lattice::integrated_star_energy_diff_combination(
    int plaquette_index,
    double imag_time_1,
    double imag_time_2,
    std::span<const double> spin_flip_lookup,   // aligned with plaquette_edges
    double imag_time_tuple_flip)
{
    if (imag_time_1 == imag_time_2) [[unlikely]] {
        throw std::invalid_argument("integrated_star_energy_diff_combination: Time interval has to be non-zero.");
    } 

    const auto plaquette_edges = get_plaquette_edges(plaquette_index);
    const auto& cached_vertices = plaquette_vertices_cache_[static_cast<size_t>(plaquette_index)];

    SmallIndexVector star_centers;
    star_centers.assign(cached_vertices.begin(), cached_vertices.end());
    SmallEnergyVector star_potential_energy_diffs;
    star_potential_energy_diffs.reserve(star_centers.size());

    thread_local std::vector<std::pair<int,int>> edge_vertices;
    edge_vertices.clear();
    edge_vertices.reserve(plaquette_edges.size());
    for (const auto& e : plaquette_edges) {
        edge_vertices.emplace_back(vertices_of_edge(e));
    }

    double energy = 0.;
    thread_local std::vector<std::pair<double,int>> local_spin_flips;
    local_spin_flips.clear();
    local_spin_flips.reserve(1 + plaquette_edges.size());

    for (const auto& v : star_centers) {
        const auto star_edges = get_star_edges(v);
        local_spin_flips.clear();
        insert_sorted_by_time(local_spin_flips, imag_time_tuple_flip, 1);

        // Singles are only on plaquette edges incident to this star center.
        for (size_t idx = 0; idx < edge_vertices.size(); ++idx) {
            const auto& [source_v, target_v] = edge_vertices[idx];
            if (source_v == v || target_v == v) {
                insert_sorted_by_time(local_spin_flips, spin_flip_lookup[idx], 0);
            }
        }

        double energy_v = integrated_tuple_energy_diff_combination_from_flips(
            star_edges, imag_time_1, imag_time_2, local_spin_flips, BASIS == 'x'
        );
        star_potential_energy_diffs.emplace_back(energy_v);
        energy += energy_v;
    }
    return {energy, std::move(star_centers), std::move(star_potential_energy_diffs)};
}

[[gnu::hot]]
std::tuple<double, Lattice::SmallIndexVector, Lattice::SmallEnergyVector> 
Lattice::integrated_plaquette_energy_diff_combination(
    int star_index,
    double imag_time_1,
    double imag_time_2,
    std::span<const double> spin_flip_lookup,   // aligned with star_edges
    double imag_time_tuple_flip)
{
    if (imag_time_1 == imag_time_2) [[unlikely]] {
        throw std::invalid_argument(
            "integrated_plaquette_energy_diff_combination: Time interval has to be non-zero.");
    }

    const auto star_edges = get_star_edges(star_index);
    const auto& cached_plaquettes = star_plaquettes_cache_[static_cast<size_t>(star_index)];
    SmallIndexVector unique_plaquettes;
    unique_plaquettes.assign(cached_plaquettes.begin(), cached_plaquettes.end());

    thread_local std::vector<std::pair<double,int>> flips;
    flips.reserve(8);

    SmallEnergyVector plaquette_potential_energy_diffs;
    plaquette_potential_energy_diffs.reserve(unique_plaquettes.size());

    double energy_sum = 0.0;

    for (int p_index : unique_plaquettes) {
        const auto plaq = get_plaquette_edges(p_index); // span<const Edge>
        flips.clear();
        insert_sorted_by_time(flips, imag_time_tuple_flip, 1);

        // Add singles on star-edges that belong to this plaquette.
        for (size_t i = 0; i < star_edges.size(); ++i) {
            const auto& p_lookup = g[star_edges[i]].part_of_plaquette_lookup;
            if (std::find(p_lookup.begin(), p_lookup.end(), p_index) != p_lookup.end()) {
                insert_sorted_by_time(flips, spin_flip_lookup[i], 0);
            }
        }

        const double ep = integrated_tuple_energy_diff_combination_from_flips(
            plaq, imag_time_1, imag_time_2, flips, BASIS == 'z');
        plaquette_potential_energy_diffs.emplace_back(ep);
        energy_sum += ep;
    }

    return { energy_sum, std::move(unique_plaquettes), std::move(plaquette_potential_energy_diffs) };
}

double Lattice::integrated_edge_energy_diff(const Edge& edg, double imag_time_1, double imag_time_2) {
    return - 2. * integrated_edge_energy(edg, imag_time_1, imag_time_2);
}

double Lattice::total_integrated_edge_energy() {
    double total_integrated_edge_energy = 0.;
    for (const auto& edg : egde_cache_) {
        total_integrated_edge_energy += integrated_edge_energy(edg, 0, BETA);
    }
    return total_integrated_edge_energy;
}

double Lattice::total_integrated_edge_energy_weighted() {
    double total_integrated_edge_energy = 0.;
    for (const auto& edg : egde_cache_) {
        total_integrated_edge_energy += integrated_edge_energy_weighted(edg, 0, BETA);
    }
    return total_integrated_edge_energy;
}

void Lattice::init_potential_energy() {
    for (const auto& edg : egde_cache_) {
        set_potential_edge_energy(edg, integrated_edge_energy(edg, 0, BETA));
    }
    if (BASIS == 'x') {
        for (size_t star_center = 0; star_center < (size_t)get_vertex_count(); ++star_center) {
            const auto& sedges = get_star_edges(star_center);
            set_potential_star_energy(star_center, integrated_tuple_energy_single_flips(sedges, 0, BETA));
        }
    } else {
        for (size_t plaquette_index = 0; plaquette_index < plaquette_vector.size(); ++plaquette_index) {
            const auto& pedges = get_plaquette_edges(plaquette_index);
            set_potential_plaquette_energy(plaquette_index, integrated_tuple_energy_single_flips(pedges, 0, BETA));
        }
    }
}

double Lattice::get_diag_single_energy() {
    const int energy = std::accumulate(
        egde_cache_.begin(), 
        egde_cache_.end(), 
        0, 
        [&](int lhs, const Edge& rhs) {return lhs + get_spin(rhs);}
    );
    return static_cast<double>(energy); 
}

std::complex<double> Lattice::get_diag_M_M() {
    const int magnetization = std::accumulate(
        egde_cache_.begin(), 
        egde_cache_.end(), 
        0, 
        [&](int lhs, const Edge& rhs) {return lhs + get_spin(rhs);}
    );

    double integrated_magnetization = total_integrated_edge_energy();
    return {static_cast<double>(integrated_magnetization / (double)(get_edge_count()) ), static_cast<double>(magnetization / (double)(get_edge_count()))}; 
}

std::complex<double> Lattice::get_diag_dynamical_M_M() {
    const int magnetization = std::accumulate(
        egde_cache_.begin(), 
        egde_cache_.end(), 
        0, 
        [&](int lhs, const Edge& rhs) {return lhs + get_spin(rhs);}
    );

    double integrated_magnetization = total_integrated_edge_energy_weighted();
    // Eq. (9) style estimators correspond to 1/2 * \int_0^\beta min(tau, beta-tau) C(tau) dtau.
    // integrated_edge_energy_weighted() returns the full triangular-kernel integral, so we apply 1/2 here.
    integrated_magnetization *= 0.5;
    return {static_cast<double>(integrated_magnetization / (double)(get_edge_count()) ), static_cast<double>(magnetization) / (double)(get_edge_count())}; 
}

std::complex<double> Lattice::get_non_diag_M_M() {
    double k_total = 0.0;
    for (const auto& edg : egde_cache_) {
        k_total += g[edg].single_spin_flips.size();
    }
    // Both components carry the raw count; the off-diagonal reducer uses real().
    return {k_total, k_total};
}

// Integrate sigma_e(tau) * w(tau) over [imag_time_1, imag_time_2]
// using w(tau)=min(tau, beta-tau) on [0, beta].
double Lattice::integrated_edge_energy_weighted(
    const Edge& edg,
    double imag_time_1,
    double imag_time_2)
{
    // Handle wrap-around by splitting into two monotone segments
    if (imag_time_1 == imag_time_2) {
        throw std::invalid_argument("integrated_edge_energy_weighted: time interval must be non-zero");
    }
    if (imag_time_1 > imag_time_2) {
        return integrated_edge_energy_weighted(edg, imag_time_1, BETA)
             + integrated_edge_energy_weighted(edg, 0.0,        imag_time_2);
    }

    return detail::weighted_spin_integral(
        {get_spin(edg), g[edg].spin_flips}, imag_time_1, imag_time_2,
        [this](double a, double b) {
            return detail::triangular_weight_primitive(b, BETA)
                 - detail::triangular_weight_primitive(a, BETA);
        }
    );
}

std::complex<double> Lattice::get_kL_kR_single() {
    // Counts refer to the current cut at zero and beta/2. The caller controls
    // when the time origin is randomized with rotate_imag_time().
    double kL = 0.0, kR = 0.0;
    const double half = 0.5 * BETA;

    for (const auto& edg : egde_cache_) {
        const auto& flips = g[edg].single_spin_flips; // sorted, in [0, beta)
        // Split the single-spin event count at the midpoint of the period.
        for (double t : flips) {
            if (t < half) ++kL; else ++kR;
        }
    }
    // Pack (k_L, k_R)
    return {kL, kR};
}

double Lattice::get_non_diag_single_energy_x() {
    double energy_beta_lmbda = 0.;
    for (const auto& edg : egde_cache_) {
        energy_beta_lmbda += g[edg].single_spin_flips.size(); 
    }
    return energy_beta_lmbda / BETA;
}

double Lattice::get_non_diag_single_energy_z() {
    double energy_beta_h = 0.;
    for (const auto& edg : egde_cache_) {
        energy_beta_h += g[edg].single_spin_flips.size();
    }
    return energy_beta_h / BETA;
}

double Lattice::get_non_diag_tuple_energy_x() {
    double energy_beta_J = 0.;
    for (size_t plaquette_index = 0; plaquette_index < plaquette_vector.size(); ++plaquette_index) {
        energy_beta_J += plaquette_flip_vector[plaquette_index].size();
    }
    return energy_beta_J / BETA;
}

double Lattice::get_non_diag_tuple_energy_z() {
    double energy_beta_mu = 0.;
    for (size_t star_center = 0; star_center < (size_t)get_vertex_count(); ++star_center) {
        energy_beta_mu += g[star_center].star_flips.size();
    }
    return energy_beta_mu / BETA;
}

double Lattice::get_diag_tuple_energy_x() {
    const auto vertices_it = boost::make_iterator_range(boost::vertices(g));
    const int energy = std::accumulate(
        vertices_it.begin(), 
        vertices_it.end(), 
        0, 
        [&](int lhs, int rhs) {return lhs + get_vertex_nn_spins_prod(rhs);}
    );
    return static_cast<double>(energy);
}

double Lattice::get_diag_tuple_energy_z() {
    int sum = 0;
    const int P = get_plaquette_count();
    for (int p = 0; p < P; ++p) {
        const auto& pedges = get_plaquette_edges(p);       
        sum += get_tuple_prod(pedges);                  
    }
    return static_cast<double>(sum);
}

std::complex<double> Lattice::fredenhagen_marcu() {
    const int half_prod = std::accumulate(
        half_path_vector.begin(), 
        half_path_vector.end(), 
        1, 
        [&](int lhs, const std::pair<int, int>& rhs) {return lhs * get_spin(edge_in_between(rhs.first, rhs.second));}
    );
    const int full_prod = std::accumulate(
        full_path_vector.begin(), 
        full_path_vector.end(), 
        1, 
        [&](int lhs, const std::pair<int, int>& rhs) {return lhs * get_spin(edge_in_between(rhs.first, rhs.second));}
    );
    return std::complex<double> {static_cast<double>(half_prod), static_cast<double>(full_prod)};
}

double Lattice::get_staggered_imaginary_times_plaquette() {
    double energy = 0.;
    const int plaquette_index = get_random_plaquette_index();
    
    std::vector<double>& plaquette_flips = plaquette_flip_vector[plaquette_index];
    int plaquette_flip_count = plaquette_flips.size();

    if (plaquette_flip_count > 0) [[likely]] {
        for (int i = 0; i < plaquette_flip_count; ++i) {
            if (i == 0) [[unlikely]] {
                energy += (i%2?-1.:1.) * (plaquette_flips[i]-0);
            } else [[likely]] {
                energy += (i%2?-1.:1.) * (plaquette_flips[i]-plaquette_flips[i-1]);
            }
        }
        energy += ((plaquette_flip_count)%2?-1.:1.) * (BETA - plaquette_flips[plaquette_flip_count-1]);
    } else [[unlikely]] {
        energy += BETA;
    }
    
    return energy / BETA;
}

double Lattice::get_staggered_imaginary_times_star() {
    double energy = 0.;
    const int center_index = get_random_vertex();
    
    std::vector<double>& star_flips = g[center_index].star_flips;
    int star_flip_count = star_flips.size();

    if (star_flip_count > 0) [[likely]] {
        for (int i = 0; i < star_flip_count; ++i) {
            if (i == 0) [[unlikely]] {
                energy += (i%2?-1.:1.) * (star_flips[i]-0);
            } else [[likely]] {
                energy += (i%2?-1.:1.) * (star_flips[i]-star_flips[i-1]);
            }
        }
        energy += ((star_flip_count)%2?-1.:1.) * (BETA - star_flips[star_flip_count-1]);
    } else {
        energy += BETA;
    }

    return energy / BETA;
}

bool Lattice::is_winding_percolating() { 
    // Set boundary values based on BOUNDARIES and dimensionality.
    const int x_left = (BOUNDARIES == "periodic") ? int(0.5 * MAX_COORDINATES[0]) : 0;
    const int y_top  = (BOUNDARIES == "periodic") ? int(0.5 * MAX_COORDINATES[1]) : 0;
    const int z_shallow = (LATTICE_DIMENSIONALITY == 3 && BOUNDARIES == "periodic")
                              ? int(0.5 * MAX_COORDINATES[2])
                              : 0;

    // Get all vertices from the graph.
    const auto vs = boost::make_iterator_range(boost::vertices(g));

    // A lambda that returns the coordinate of a vertex.
    auto get_coord = [this](int vertex, char coord) -> int {
        switch (coord) {
            case 'x': return int(g[vertex].x);
            case 'y': return int(g[vertex].y);
            case 'z': return int(g[vertex].z);
            default:  throw std::invalid_argument("Invalid coordinate");
        }
    };

    // Lambda to perform DFS along a specified direction.
    auto check_direction = [this, &vs, get_coord](char coord, int boundary_value) -> bool {
        std::vector<int> filtered;
        std::copy_if(vs.begin(), vs.end(), std::back_inserter(filtered),
                     [=](int v) { return get_coord(v, coord) == boundary_value; });
        
        // Create a discovered vector for DFS.
        std::vector<bool> discovered(get_vertex_count(), false);

        // For each starting vertex, perform DFS.
        for (int start : filtered) {
            try {
                // Storage for winding numbers.
                std::vector<int> winding(get_vertex_count(), INT_MAX);
                // Stack holds pairs: (vertex, current winding number).
                std::stack<std::pair<int, int>> dfs_stack;
                dfs_stack.push({start, 0});
                
                while (!dfs_stack.empty()) {
                    auto [v, wn] = dfs_stack.top();
                    dfs_stack.pop();

                    if (discovered[v])
                        continue;
                    discovered[v] = true;
                    winding[v] = wn;

                    // Traverse neighbors.
                    for (int neighbor : boost::make_iterator_range(boost::adjacent_vertices(v, g))) {
                        // Check if edge between v and neighbor is active.
                        if (get_spin(edge_in_between(v, neighbor)) == -1) {
                            int new_wn = wn;
                            int coord_v = get_coord(v, coord);
                            int coord_n = get_coord(neighbor, coord);
                            if (coord_n > boundary_value && coord_v == boundary_value)
                                new_wn = wn + 1;
                            else if (coord_v > boundary_value && coord_n == boundary_value)
                                new_wn = wn - 1;
                            
                            if (!discovered[neighbor])
                                dfs_stack.push({neighbor, new_wn});
                            else if (winding[neighbor] != new_wn)
                                throw FoundPercolation();
                        }
                    }
                }
            } catch (const FoundPercolation&) {
                return true;
            }
        }
        return false;
    };

    // Check percolation along x, y and, if in 3D, z directions.
    if (check_direction('x', x_left)) return true;
    if (check_direction('y', y_top))  return true;
    if (LATTICE_DIMENSIONALITY == 3 && check_direction('z', z_shallow))
        return true;
    
    return false;
}

bool Lattice::is_percolating() {
    int x_left = static_cast<int>(0 * MAX_COORDINATES[0]);
    int x_right = static_cast<int>(1 * MAX_COORDINATES[0]);

    int y_top = static_cast<int>(0 * MAX_COORDINATES[1]);
    int y_bottom = static_cast<int>(1 * MAX_COORDINATES[1]);

    int z_shallow{}, z_deep{};

    if (LATTICE_DIMENSIONALITY == 3) {
        z_shallow = static_cast<int>(0 * MAX_COORDINATES[2]);
        z_deep = static_cast<int>(1 * MAX_COORDINATES[2]);
    }

    // Get all vertices using Boost.
    const auto vs = boost::make_iterator_range(boost::vertices(g));

    // Helper lambda to fetch a vertex's coordinate.
    auto get_coord = [this](int v, char axis) -> int {
        switch (axis) {
            case 'x': return int(g[v].x);
            case 'y': return int(g[v].y);
            case 'z': return int(g[v].z);
            default:  throw std::invalid_argument("Invalid axis");
        }
    };

    // Lambda to check percolation along a given axis.
    // It performs a DFS from vertices on the 'start_bound' side and checks for any connection
    // to vertices on the 'target_bound' side through active edges (spin == -1).
    auto check_direction = [this, &vs, &get_coord](char axis, int start_bound, int target_bound) -> bool {
        std::vector<int> start_vertices;
        std::copy_if(vs.begin(), vs.end(), std::back_inserter(start_vertices),
                     [=](int v) { return get_coord(v, axis) == start_bound; });
        
        // Build a set of target vertices for fast lookup.
        std::unordered_set<int> target_set;
        for (int v : vs) {
            if (get_coord(v, axis) == target_bound)
                target_set.insert(v);
        }
        
        std::vector<bool> discovered(get_vertex_count(), false);
        std::stack<int> dfs;
        for (int start : start_vertices) {
            if (discovered[start])
                continue;
            dfs.push(start);
            while (!dfs.empty()) {
                int v = dfs.top();
                dfs.pop();
                if (discovered[v])
                    continue;
                discovered[v] = true;
                if (target_set.find(v) != target_set.end())
                    return true;
                for (int neighbor : boost::make_iterator_range(boost::adjacent_vertices(v, g))) {
                    if (get_spin(edge_in_between(v, neighbor)) != -1 || discovered[neighbor])
                        continue;

                    // Open-style percolation must not traverse periodic seam links.
                    // On periodic lattices those are exactly start-boundary <-> target-boundary hops.
                    const int coord_v = get_coord(v, axis);
                    const int coord_n = get_coord(neighbor, axis);
                    const bool crosses_seam =
                        (coord_v == start_bound && coord_n == target_bound) ||
                        (coord_v == target_bound && coord_n == start_bound);
                    if (crosses_seam)
                        continue;

                    dfs.push(neighbor);
                }
            }
        }
        return false;
    };

    if (check_direction('x', x_left, x_right))
        return true;
    if (check_direction('y', y_top, y_bottom))
        return true;
    if (LATTICE_DIMENSIONALITY == 3 && check_direction('z', z_shallow, z_deep))
        return true;
    
    return false;
}

bool Lattice::is_winding_cube_percolating() {
    if (LATTICE_DIMENSIONALITY != 3) {
        throw std::invalid_argument("Lattice dimensionality has to be 3.");
    }

    int x_left = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[0]);

    int y_top = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[1]);

    int z_shallow = 0;

    if (BOUNDARIES == "periodic") {
        x_left = static_cast<int>(0.5 * MAX_PLAQUETTE_COORDINATES[0]);

        y_top = static_cast<int>(0.5 * MAX_PLAQUETTE_COORDINATES[1]);
    }

    if (LATTICE_DIMENSIONALITY == 3) {
        z_shallow = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[2]);
        if (BOUNDARIES == "periodic") {
            z_shallow = static_cast<int>(0.5 * MAX_PLAQUETTE_COORDINATES[2]);
        }
    }

    // Create a list of all cube indices.
    std::vector<int> cubes(get_cube_count());
    std::iota(cubes.begin(), cubes.end(), 0);

    // Lambda to perform DFS for cube percolation along a given axis.
    // cube_coord: the coordinate vector (cube_x_vector, cube_y_vector, or cube_z_vector).
    // boundary: the starting boundary value for that axis.
    auto check_cube_direction = [this, &cubes](const std::vector<double>& cube_coord, int boundary) -> bool {
        // Filter cubes on the start boundary.
        std::vector<int> start;
        std::copy_if(cubes.begin(), cubes.end(), std::back_inserter(start),
                     [&](int i) { return static_cast<int>(cube_coord[i]) == boundary; });
        
        std::vector<bool> discovered(get_cube_count(), false);
        std::vector<int> winding(get_cube_count(), INT_MAX);
        std::stack<std::pair<int, int>> dfs;  // Pair: {cube index, winding number}

        // For each starting cube, do a DFS.
        for (int c : start) {
            if (!discovered[c])
                dfs.push({c, 0});
            while (!dfs.empty()) {
                auto [cur, wn] = dfs.top();
                dfs.pop();
                if (discovered[cur])
                    continue;
                discovered[cur] = true;
                winding[cur] = wn;
                
                // Iterate over all plaquettes adjacent to cube 'cur'.
                for (int p_index : cube_has_plaquettes_lookup[cur]) {
                    // A shared face connects cubes when its spin product is +1.
                    const auto& pedges = get_plaquette_edges(p_index);
                    if (get_tuple_prod(pedges) == 1) {
                        // For each neighboring cube via this plaquette.
                        for (int neighbor : plaquette_part_of_cube_lookup[p_index]) {
                            if (neighbor == cur)
                                continue;
                            int new_wn = wn;
                            // Update winding based on crossing the boundary.
                            if (static_cast<int>(cube_coord[neighbor]) > boundary 
                            && static_cast<int>(cube_coord[cur]) == boundary)
                                new_wn = wn + 1;
                            else if (static_cast<int>(cube_coord[cur]) > boundary 
                            && static_cast<int>(cube_coord[neighbor]) == boundary)
                                new_wn = wn - 1;
                            
                            if (!discovered[neighbor])
                                dfs.push({neighbor, new_wn});
                            else if (winding[neighbor] != new_wn)
                                return true; // Conflict detected: percolation!
                        }
                    }
                }
            }
        }
        return false;
    };

    // Check percolation in the x, y, and z directions.
    if (check_cube_direction(cube_x_vector, x_left))
        return true;
    if (check_cube_direction(cube_y_vector, y_top))
        return true;
    if (check_cube_direction(cube_z_vector, z_shallow))
        return true;
    
    return false;
} 

bool Lattice::is_winding_plaquette_percolating() {
    int x_left = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[0]);

    int y_top = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[1]);

    int z_shallow = 0;

    if (BOUNDARIES == "periodic") {
        x_left = static_cast<int>(0.5 * MAX_PLAQUETTE_COORDINATES[0]);

        y_top = static_cast<int>(0.5 * MAX_PLAQUETTE_COORDINATES[1]);
    }

    if (LATTICE_DIMENSIONALITY == 3) {
        z_shallow = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[2]);
        if (BOUNDARIES == "periodic") {
            z_shallow = static_cast<int>(0.5 * MAX_PLAQUETTE_COORDINATES[2]);
        }
    }

    // Build a list of all plaquette indices.
    std::vector<int> ps(get_plaquette_count());
    std::iota(ps.begin(), ps.end(), 0);

    // Lambda for DFS on plaquettes along a given axis.
    // coord_vector is the coordinate vector for the axis (e.g. plaquette_x_vector),
    // boundary is the "start" boundary value.
    auto check_direction = [this, &ps](const std::vector<double>& coord_vector, int boundary) -> bool {
        // Filter plaquettes on the start boundary.
        std::vector<int> start;
        std::copy_if(ps.begin(), ps.end(), std::back_inserter(start),
                     [&](int i) { return static_cast<int>(coord_vector[i]) == boundary; });
        
        std::vector<bool> discovered(get_plaquette_count(), false);
        std::vector<int> winding(get_plaquette_count(), INT_MAX);
        std::stack<std::pair<int, int>> stack;  // {plaquette index, winding number}
        
        for (int p : start) {
            if (!discovered[p])
                stack.push({p, 0});
            while (!stack.empty()) {
                auto [cur, wn] = stack.top();
                stack.pop();
                if (discovered[cur])
                    continue;
                discovered[cur] = true;
                winding[cur] = wn;
                // Traverse neighboring plaquettes via each edge of current plaquette.
                const auto& pedges = get_plaquette_edges(cur);
                for (auto edg : pedges) {
                    if (get_spin(edg) == -1) {
                        for (int p_neighbor : g[edg].part_of_plaquette_lookup) {
                            if (p_neighbor == cur)
                                continue;
                            int new_wn = wn;
                            // Update winding number based on crossing the boundary.
                            if (static_cast<int>(coord_vector[p_neighbor]) > boundary 
                            && static_cast<int>(coord_vector[cur]) == boundary)
                                new_wn = wn + 1;
                            else if (static_cast<int>(coord_vector[cur]) > boundary 
                            && static_cast<int>(coord_vector[p_neighbor]) == boundary)
                                new_wn = wn - 1;
                            
                            if (!discovered[p_neighbor])
                                stack.push({p_neighbor, new_wn});
                            else if (winding[p_neighbor] != new_wn)
                                return true;  // Found percolation.
                        }
                    }
                }
            }
        }
        return false;
    };

    if (check_direction(plaquette_x_vector, x_left))
        return true;
    if (check_direction(plaquette_y_vector, y_top))
        return true;
    if (LATTICE_DIMENSIONALITY == 3 && check_direction(plaquette_z_vector, z_shallow))
        return true;

    return false;
}

bool Lattice::is_plaquette_percolating() {
    int x_left = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[0]);
    int x_right = static_cast<int>(1 * MAX_PLAQUETTE_COORDINATES[0]);

    int y_top = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[1]);
    int y_bottom = static_cast<int>(1 * MAX_PLAQUETTE_COORDINATES[1]);

    int z_shallow{}, z_deep{};
    if (LATTICE_DIMENSIONALITY == 3) {
        z_shallow = static_cast<int>(0 * MAX_PLAQUETTE_COORDINATES[2]);
        z_deep = static_cast<int>(1 * MAX_PLAQUETTE_COORDINATES[2]);
    }

    // Build a list of all plaquette indices.
    std::vector<int> ps(get_plaquette_count());
    std::iota(ps.begin(), ps.end(), 0);

    // Lambda to check percolation along a given axis.
    // It performs a DFS from plaquettes on the 'start_bound' side and checks
    // for any connection to plaquettes on the 'target_bound' side through active
    // links (edges with spin == -1).
    auto check_direction = [this, &ps](const std::vector<double>& coord_vector, int start_bound, int target_bound) -> bool {
        std::vector<int> start_plaquettes;
        std::copy_if(ps.begin(), ps.end(), std::back_inserter(start_plaquettes),
                     [&](int p) { return static_cast<int>(coord_vector[p]) == start_bound; });

        // Build a set of target plaquettes for fast lookup.
        std::unordered_set<int> target_set;
        for (int p : ps) {
            if (static_cast<int>(coord_vector[p]) == target_bound) {
                target_set.insert(p);
            }
        }

        std::vector<bool> discovered(get_plaquette_count(), false);
        std::stack<int> dfs;
        for (int start : start_plaquettes) {
            if (discovered[start]) {
                continue;
            }
            dfs.push(start);
            while (!dfs.empty()) {
                const int cur = dfs.top();
                dfs.pop();

                if (discovered[cur]) {
                    continue;
                }
                discovered[cur] = true;

                if (target_set.find(cur) != target_set.end()) {
                    return true;
                }

                const auto& pedges = get_plaquette_edges(cur);
                for (auto edg : pedges) {
                    if (get_spin(edg) == -1) {
                        for (int p_neighbor : g[edg].part_of_plaquette_lookup) {
                            if (p_neighbor == cur || discovered[p_neighbor]) {
                                continue;
                            }

                            // Open-style percolation must not traverse periodic seam links.
                            // On periodic lattices those are exactly start-boundary <-> target-boundary hops.
                            const int cur_coord = static_cast<int>(coord_vector[cur]);
                            const int neighbor_coord = static_cast<int>(coord_vector[p_neighbor]);
                            const bool crosses_seam =
                                (cur_coord == start_bound && neighbor_coord == target_bound) ||
                                (cur_coord == target_bound && neighbor_coord == start_bound);
                            if (crosses_seam) {
                                continue;
                            }

                            if (!discovered[p_neighbor]) {
                                dfs.push(p_neighbor);
                            }
                        }
                    }
                }
            }
        }
        return false;
    };

    if (check_direction(plaquette_x_vector, x_left, x_right)) {
        return true;
    }
    if (check_direction(plaquette_y_vector, y_top, y_bottom)) {
        return true;
    }
    if (LATTICE_DIMENSIONALITY == 3 && check_direction(plaquette_z_vector, z_shallow, z_deep)) {
        return true;
    }

    return false;
}

int Lattice::largest_plaquette_cluster() {
    const size_t P = static_cast<size_t>(get_plaquette_count());
    if (P == 0) return 0;

    thread_local std::vector<int> parent, rankv, comp_size;
    if (parent.size() < P) {
        parent.resize(P);
        rankv.resize(P);
        comp_size.resize(P);
    }

    std::iota(parent.begin(), parent.begin() + P, 0);
    std::fill(rankv.begin(), rankv.begin() + P, 0);
    std::fill(comp_size.begin(), comp_size.begin() + P, 1);

    // Path-halving 
    auto find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    auto unite = [&](int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) return;
        if (rankv[ra] < rankv[rb]) std::swap(ra, rb);
        parent[rb] = ra;
        if (rankv[ra] == rankv[rb]) ++rankv[ra];
        comp_size[ra] += comp_size[rb];
    };

    for (auto e : egde_cache_) {
        if (g[e].spin != -1) continue;
        const auto& pls = g[e].part_of_plaquette_lookup; 
        if (pls.size() >= 2) {
            const int p0 = pls[0];
            for (size_t i = 1; i < pls.size(); ++i) {
                unite(p0, pls[i]);
            }
        }
    }

    int best = 0;
    for (int p = 0; p < static_cast<int>(P); ++p) {
        if (parent[p] == p) best = std::max(best, comp_size[p]);
    }
    return best;
}


int Lattice::largest_cluster() {
    const size_t N = get_vertex_count();

    thread_local std::vector<int> parent, rankv, edge_count;
    if (parent.size() < N) {
        parent.resize(N);
        rankv.resize(N);
        edge_count.resize(N);
    }

    std::iota(parent.begin(), parent.begin() + N, 0);
    std::fill(rankv.begin(), rankv.begin() + N, 0);
    std::fill(edge_count.begin(), edge_count.begin() + N, 0);

    // Path-halving
    auto find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    int best = 0;

    for (auto e : egde_cache_) {
        if (g[e].spin != -1) continue;

        int u = boost::source(e, g);
        int v = boost::target(e, g);

        int ru = find(u);
        int rv = find(v);

        if (ru == rv) {
            int c = ++edge_count[ru];
            if (c > best) best = c;
        }
        else {
            if (rankv[ru] < rankv[rv]) std::swap(ru, rv);
            parent[rv] = ru;
            if (rankv[ru] == rankv[rv]) ++rankv[ru];

            edge_count[ru] += edge_count[rv] + 1;
            if (edge_count[ru] > best) best = edge_count[ru];
        }
    }
    return best;
}

double Lattice::percolation_strength() {
    bool percolating = 0;
    if (BOUNDARIES == "periodic") {
        percolating = is_winding_percolating();
    } else {
        percolating = is_percolating();
    }

    if (percolating) {
        return largest_cluster() / static_cast<double>(get_edge_count());
    } else {
        return 0.;
    }
}

double Lattice::percolation_probability() {
    if (BOUNDARIES == "periodic") {
        return is_winding_percolating();
    } else {
        return is_percolating();
    }
}

double Lattice::plaquette_percolation_strength() {
    bool percolating = 0;
    if (BOUNDARIES == "periodic") {
        percolating = is_winding_plaquette_percolating();
    } else {
        percolating = is_plaquette_percolating();
    }

    if (percolating) {
        return largest_plaquette_cluster() / static_cast<double>(get_plaquette_count());
    } else {
        return 0.;
    }
}

double Lattice::plaquette_percolation_probability() {
    if (BOUNDARIES == "periodic") {
        return is_winding_plaquette_percolating();
    } else {
        return is_plaquette_percolating();
    }
}

double Lattice::cube_percolation_strength() {
    return -1.;
}

double Lattice::cube_percolation_probability() {
    if (BOUNDARIES == "periodic") {
        return is_winding_cube_percolating();
    } else {
        // Cube percolation with open boundaries is not implemented.
        return -1.;
    }
}

void Lattice::rotate_imag_time() {
    std::uniform_real_distribution<double> new_times_dist(0, BETA);
    const double tau_0 = new_times_dist(*rng);
    // Reuse one buffer across all histories. Copy the two sorted segments in
    // their new order while shifting them, avoiding an in-place rotation plus
    // a separate modulo pass over every event.
    std::vector<double> shifted_times;
    const auto rotate_times = [&](std::vector<double>& times) {
        return detail::rotate_times(times, tau_0, BETA, shifted_times);
    };

    for (const auto& edg : egde_cache_) {
        rotate_times(g[edg].single_spin_flips);
        const size_t pivot_index = rotate_times(g[edg].spin_flips);
        // flips crossing the cut = pivot_index
        if (pivot_index % 2 == 1) {
            g[edg].spin *= -1;
        }
    }

    if (BASIS == 'x') {
        for (int p_index = 0; p_index < get_plaquette_count(); ++p_index) {
            rotate_times(plaquette_flip_vector[p_index]);
        }
    } else {
        for (int s_index = 0; s_index < get_vertex_count(); ++s_index) {
            rotate_times(g[s_index].star_flips);
        } 
    }

}

void Lattice::ensure_snapshot_spool_() {
    if (snapshot_spool_.active()) return;

    const auto dir = std::filesystem::temp_directory_path();
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto self = reinterpret_cast<std::uintptr_t>(this);

    for (int attempt = 0; attempt < 100; ++attempt) {
        auto path = dir / (
            "paratoric_snapshots_" + std::to_string(stamp) + "_" +
            std::to_string(self) + "_" + std::to_string(attempt) + ".bin"
        );
        if (std::filesystem::exists(path)) continue;

        snapshot_spool_.stream.open(path, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!snapshot_spool_.stream) {
            snapshot_spool_.stream.close();
            continue;
        }

        snapshot_spool_.path = path;
        snapshot_spool_.edge_count = egde_cache_.size();
        snapshot_spool_.sample_count = 0;

        for (const auto& edg : egde_cache_) {
            g[edg].spin_string.clear();
            g[edg].spin_string.shrink_to_fit();
        }
        return;
    }

    throw std::runtime_error("Could not create temporary snapshot spool file.");
}

void Lattice::update_spin_string() {
    ensure_snapshot_spool_();

    for (const auto& edg : egde_cache_) {
        const char spin = g[edg].spin > 0 ? '1' : '0';
        snapshot_spool_.stream.write(&spin, 1);
    }
    if (!snapshot_spool_.stream) {
        throw std::runtime_error("Could not write snapshot to temporary spool file.");
    }
    ++snapshot_spool_.sample_count;
}

void Lattice::write_snapshot_graphml_from_spool_(
    const std::string& file_name,
    const std::filesystem::path& output_directory
) {
    if (snapshot_spool_.stream.is_open()) {
        snapshot_spool_.stream.flush();
        snapshot_spool_.stream.close();
    }

    if (snapshot_spool_.edge_count != egde_cache_.size()) {
        throw std::runtime_error("Snapshot spool edge count does not match current lattice.");
    }

    struct ScopedTempFile {
        std::filesystem::path path{};
        ~ScopedTempFile() {
            if (!path.empty()) {
                std::error_code ec;
                std::filesystem::remove(path, ec);
            }
        }
    };

    const std::size_t edge_count = snapshot_spool_.edge_count;
    const std::size_t sample_count = snapshot_spool_.sample_count;
    const auto max_stream_offset = static_cast<std::uintmax_t>(std::numeric_limits<std::streamoff>::max());
    if (edge_count != 0 && sample_count > max_stream_offset / edge_count) {
        throw std::runtime_error("Snapshot spool is too large for stream offsets.");
    }

    ScopedTempFile column_spool;
    {
        const auto dir = std::filesystem::temp_directory_path();
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        const auto self = reinterpret_cast<std::uintptr_t>(this);

        for (int attempt = 0; attempt < 100; ++attempt) {
            auto path = dir / (
                "paratoric_snapshots_columns_" + std::to_string(stamp) + "_" +
                std::to_string(self) + "_" + std::to_string(attempt) + ".bin"
            );
            if (std::filesystem::exists(path)) continue;

            std::fstream column_stream(path, std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
            if (!column_stream) {
                column_stream.close();
                continue;
            }

            column_spool.path = path;
            const std::size_t target_tile_bytes = 8ULL * 1024ULL * 1024ULL;
            const std::size_t samples_per_tile = edge_count == 0
                ? 1
                : std::max<std::size_t>(1, target_tile_bytes / edge_count);
            std::vector<char> row_tile(samples_per_tile * edge_count);
            std::vector<char> edge_tile(samples_per_tile);

            std::ifstream row_stream(snapshot_spool_.path, std::ios::binary | std::ios::in);
            if (!row_stream) {
                throw std::runtime_error("Could not reopen temporary snapshot spool file.");
            }

            for (std::size_t sample_begin = 0; sample_begin < sample_count; sample_begin += samples_per_tile) {
                const std::size_t tile_samples = std::min(samples_per_tile, sample_count - sample_begin);
                const std::size_t tile_bytes = tile_samples * edge_count;
                row_stream.read(row_tile.data(), static_cast<std::streamsize>(tile_bytes));
                if (row_stream.gcount() != static_cast<std::streamsize>(tile_bytes)) {
                    throw std::runtime_error("Snapshot spool file ended unexpectedly.");
                }

                for (std::size_t edge_index = 0; edge_index < edge_count; ++edge_index) {
                    for (std::size_t sample = 0; sample < tile_samples; ++sample) {
                        edge_tile[sample] = row_tile[sample * edge_count + edge_index];
                    }

                    const auto offset = static_cast<std::streamoff>(
                        edge_index * sample_count + sample_begin
                    );
                    column_stream.seekp(offset);
                    column_stream.write(edge_tile.data(), static_cast<std::streamsize>(tile_samples));
                    if (!column_stream) {
                        throw std::runtime_error("Could not write transposed snapshot spool file.");
                    }
                }
            }
            column_stream.close();
            break;
        }

        if (column_spool.path.empty()) {
            throw std::runtime_error("Could not create temporary transposed snapshot spool file.");
        }
    }

    std::filesystem::path path_file(file_name + ".xml");
    std::filesystem::path path_out = output_directory / path_file;
    std::ofstream graphml_file(path_out);
    if (!graphml_file) {
        throw std::runtime_error(std::format("Could not open GraphML output file \"{}\".", path_out.string()));
    }

    graphml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    graphml_file << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n";
    graphml_file << "  <key id=\"key0\" for=\"edge\" attr.name=\"spin\" attr.type=\"string\" />\n";
    graphml_file << "  <key id=\"key1\" for=\"node\" attr.name=\"x\" attr.type=\"double\" />\n";
    graphml_file << "  <key id=\"key2\" for=\"node\" attr.name=\"y\" attr.type=\"double\" />\n";
    graphml_file << "  <key id=\"key3\" for=\"node\" attr.name=\"z\" attr.type=\"double\" />\n";
    graphml_file << "  <graph id=\"G\" edgedefault=\"undirected\" parse.nodeids=\"canonical\" parse.edgeids=\"canonical\" parse.order=\"nodesfirst\">\n";

    for (const auto& v : boost::make_iterator_range(boost::vertices(g))) {
        graphml_file << "    <node id=\"n" << v << "\">\n";
        graphml_file << "      <data key=\"key1\">" << g[v].x << "</data>\n";
        graphml_file << "      <data key=\"key2\">" << g[v].y << "</data>\n";
        graphml_file << "      <data key=\"key3\">" << g[v].z << "</data>\n";
        graphml_file << "    </node>\n";
    }

    std::ifstream column_stream(column_spool.path, std::ios::binary | std::ios::in);
    if (!column_stream) {
        throw std::runtime_error("Could not reopen transposed snapshot spool file.");
    }
    const std::size_t xml_chunk_samples = 64ULL * 1024ULL;
    std::vector<char> edge_chunk(std::min(xml_chunk_samples, std::max<std::size_t>(sample_count, 1)));

    for (std::size_t edge_index = 0; edge_index < edge_count; ++edge_index) {
        const auto& edg = egde_cache_[edge_index];
        graphml_file << "    <edge id=\"e" << edge_index << "\" source=\"n"
                     << boost::source(edg, g) << "\" target=\"n" << boost::target(edg, g) << "\">\n";
        graphml_file << "      <data key=\"key0\">";

        const auto offset = static_cast<std::streamoff>(edge_index * sample_count);
        column_stream.seekg(offset);
        bool first_sample = true;
        for (std::size_t sample_begin = 0; sample_begin < sample_count; sample_begin += xml_chunk_samples) {
            const std::size_t chunk_samples = std::min(xml_chunk_samples, sample_count - sample_begin);
            column_stream.read(edge_chunk.data(), static_cast<std::streamsize>(chunk_samples));
            if (column_stream.gcount() != static_cast<std::streamsize>(chunk_samples)) {
                throw std::runtime_error("Transposed snapshot spool file ended unexpectedly.");
            }

            for (std::size_t sample = 0; sample < chunk_samples; ++sample) {
                if (!first_sample) graphml_file << ' ';
                first_sample = false;
                if (edge_chunk[sample] == '1') graphml_file << '1';
                else graphml_file << "-1";
            }
        }
        graphml_file << "</data>\n";
        graphml_file << "    </edge>\n";
    }

    graphml_file << "  </graph>\n";
    graphml_file << "</graphml>\n";
    if (!graphml_file) {
        throw std::runtime_error(std::format("Could not finish GraphML output file \"{}\".", path_out.string()));
    }

    snapshot_spool_.reset();
}

void Lattice::write_graph(const std::string& file_name, const std::filesystem::path& output_directory) {
    if (snapshot_spool_.active()) {
        write_snapshot_graphml_from_spool_(file_name, output_directory);
        return;
    }

    // Assemble ofstream
    std::filesystem::path path_file(file_name + ".xml");
    std::filesystem::path path_out = output_directory / path_file;
    std::ofstream graphml_file(path_out);

    // Assemble dynamic properties
    boost::dynamic_properties dp;
    dp.property("spin", boost::get(&EdgeData::spin_string, g));
    dp.property("x", boost::get(&VertexData::x, g));
    dp.property("y", boost::get(&VertexData::y, g));
    dp.property("z", boost::get(&VertexData::z, g));

    // Write the output file
    boost::write_graphml(graphml_file, g, dp, true);
}

} // namespace paratoric
