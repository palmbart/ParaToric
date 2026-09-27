// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#pragma once

#include "paratoric/types/types.hpp"
#include "rng/rng.hpp"
#include "lattice/time_search.hpp"

#include <boost/container/small_vector.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <concepts>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <string>
#include <tuple>
#include <vector>

namespace paratoric {

/**
 * @brief Graph geometry, spin histories, and bare integrals for continuous-time QMC.
 *
 * Spins live on edges. Each edge stores its spin at the time origin and a
 * sorted history of all flips, with single-spin events also stored separately.
 * Tuple events are plaquette flips in the x-basis and star flips in the z-basis.
 * A completed QMC update must preserve periodic spin histories and keep the
 * edge histories consistent with the separate single and tuple histories.
 *
 * Event times lie within one period beta. Neighbor searches wrap around the
 * period; ordinary energy integrals take ordered, non-wrapping intervals.
 * Bare integrals contain spin or spin-product factors only, without Hamiltonian
 * minus signs or couplings. QMC applies those factors and maintains the active
 * energy caches after accepting an update.
 *
 * Geometry and adjacency caches are built once at construction. Edge views
 * borrow that storage; event views must be reacquired after changing the
 * corresponding history. Indices and edge descriptors must belong to this
 * lattice. Low-level setters do not validate these invariants.
 */
class Lattice {    
public:
    /** @brief Control-flow exception used to stop a percolation search early. */
    class FoundPercolation : public std::exception {
        public:
            /** @brief Return the diagnostic message for a completed percolation search. */
            const char* what() const noexcept override {
                return "Found percolation";
            }
    };

    /** @brief Star events and coordinates at one graph vertex. */
    struct VertexData {
        // Sorted star-event times (off-diagonal in the z-basis).
        std::vector<double> star_flips;
        // Cached integral of the star spin product over [0, beta], used in x.
        double integrated_star_energy;
        // Real-space coordinates used by geometry-dependent observables.
        double x = 0.;
        double y = 0.;
        double z = 0.;
    };

    /** @brief Spin history, bare edge integral, and geometry of one graph edge. */
    struct EdgeData {
        // Spin (+1 or -1) at the time origin; beta is identified with zero.
        int spin = 1;
        // Legacy GraphML property. Production snapshot histories are spooled to disk.
        std::string spin_string = "";
        // Sorted times of all spin flips, including events contributed by tuples.
        std::vector<double> spin_flips; 
        // Sorted single-spin events, also present in spin_flips.
        std::vector<double> single_spin_flips; 
        // Cached integral of the spin over [0, beta], without sign or coupling.
        double integrated_edge_energy;
        // Indices of adjacent plaquettes, populated by build_caches_().
        std::vector<int> part_of_plaquette_lookup;
        int source_vertex = -1;
        int target_vertex = -1;
        // Cubic edge direction (x, y, or z), used by cube percolation.
        std::string orientation;

        /** @brief Set the edge's direction label, defaulting to x. */
        EdgeData(std::string o = "x") : orientation(o){ }
    };

    using RNG = paratoric::rng::RNG;
    using LatticeGraph = boost::adjacency_list<
        boost::vecS, boost::vecS,
        boost::undirectedS,
        VertexData,
        EdgeData
    >;
    using VertexWhitelist = std::set<int>; // TODO use vertex descriptor instead of int
    using VertexPair = std::pair<int, int>; // TODO use vertex descriptor instead of int
    using Vertex = boost::graph_traits<LatticeGraph>::vertex_descriptor;
    using Edge = boost::graph_traits<LatticeGraph>::edge_descriptor;
    using SmallIndexVector = boost::container::small_vector<int, 8>;
    using SmallEnergyVector = boost::container::small_vector<double, 8>;
    
    /**
     * @brief Create data structure for storing and modifying a lattice for continuous QMC of any graph geometry.
     * 
     * @param lat_spec the lattice specification object
     * @param lat_spec.basis the spin basis (x or z)
     * @param lat_spec.lattice_type the lattice type 
     * @param lat_spec.system_size the system size of the lattice (length in 1D, width/height in 2D)
     * @param lat_spec.beta the inverse temperature beta (needed for the imaginary time dimension)
     * @param lat_spec.boundaries the boundary condition of the lattice (periodic or open)
     * @param lat_spec.default_spin the default spin on edges (-1 or 1)
     * @param rng (optional) shared pointer to RNG (MT19937). If null, a default RNG is created.
     * 
     */
    Lattice(const LatSpec& lat_spec, std::shared_ptr<RNG> rng = nullptr) : rng(rng ? std::move(rng) : std::make_shared<RNG>()) {
        BASIS = lat_spec.basis;
        // Geometry construction determines LATTICE_DIMENSIONALITY.
        LATTICE_TYPE = lat_spec.lattice_type;
        SYSTEM_SIZE = lat_spec.system_size;
        BETA = lat_spec.beta;
        BOUNDARIES = lat_spec.boundaries;
        DEFAULT_SPIN = lat_spec.default_spin;

        check_input_validity();

        init_lattice_graph();

        build_caches_();
        init_potential_energy();

    }  

    /** @brief Create an uninitialized placeholder; assign a constructed lattice before use. */
    Lattice() = default; 
    /** @brief Release lattice storage and clean up any temporary snapshot spool. */
    ~Lattice() = default;

    /** @brief Copy lattice state and share its RNG; start with an empty snapshot spool. */
    Lattice(Lattice const&) = default;
    /** @brief Copy lattice state and share its RNG; discard the destination's snapshot spool. */
    Lattice& operator=(Lattice const&) = default;

    /** @brief Count edges with spin +1 at the time origin. */
    inline int get_non_string_count();

    /** @brief Count edges with spin -1 at the time origin. */
    inline int get_string_count();

    /** @brief Return the number of graph vertices. */
    inline int get_vertex_count();

    /** @brief Return the number of graph edges. */
    inline int get_edge_count();

    /** @brief Return the number of elementary plaquettes. */
    inline int get_plaquette_count();

    /** @brief Return the number of elementary cubes. */
    inline int get_cube_count();

    /** @brief Count star defects in x or plaquette defects in z at the time origin. */
    int get_anyon_count();

    /**
     * @brief Replace the shared RNG used for lattice proposals and time rotations.
     * @param rng_inp Non-null shared generator; ownership is shared, not cloned.
     */
    void set_rng(std::shared_ptr<RNG>& rng_inp);

    /** @brief Return shared ownership of the lattice's RNG. */
    std::shared_ptr<RNG> get_rng();

    /** @brief Return the edge spin at the time origin. */
    inline int get_spin(const Edge& edg);
    /** @brief Read the cached bare edge integral over [0, beta]. */
    inline double get_potential_edge_energy(const Edge& edg);
    /** @brief Replace the cached bare edge integral. */
    inline void set_potential_edge_energy(const Edge& edg, double potential_energy);
    /** @brief Add a bare integral change to the edge cache. */
    inline void add_potential_edge_energy(const Edge& edg, double diff);
    /** @brief Read the cached bare star integral (maintained in the x-basis). */
    inline double get_potential_star_energy(int star_index);
    /** @brief Replace the cached bare star integral. */
    inline void set_potential_star_energy(int star_index, double potential_energy);
    /** @brief Add a bare integral change to the star cache. */
    inline void add_potential_star_energy(int star_index, double diff);
    /** @brief Read the cached bare plaquette integral (maintained in the z-basis). */
    inline double get_potential_plaquette_energy(int plaquette_index);
    /** @brief Replace the cached bare plaquette integral. */
    inline void set_potential_plaquette_energy(int plaquette_index, double potential_energy);
    /** @brief Add a bare integral change to the plaquette cache. */
    inline void add_potential_plaquette_energy(int plaquette_index, double diff);

    /** @brief Return the edge's geometry-specific direction label. */
    inline std::string get_orientation(const Edge& edg);

    /** @brief Return a time by index in the edge's full flip history. */
    inline double get_spin_flip_imag_time(const Edge& edg, int spin_flip_index);

    /**
     * @brief Find the first exact occurrence of tau in the edge's full history.
     * @return Zero-based index in that history.
     * @throws std::runtime_error If no event has exactly the requested time.
     */
    int get_spin_flip_index(const Edge& edg, double tau);

    /**
     * @brief Find the first exact occurrence of tau in the edge's single-spin history.
     * @return Zero-based index in that history.
     * @throws std::runtime_error If no event has exactly the requested time.
     */
    int get_single_spin_flip_index(const Edge& edg, double tau);

    /**
     * @brief Borrow the sorted single-spin event times on edg.
     * @return Read-only view; reacquire it after mutating this edge's event history.
     */
    std::span<const double> get_single_spin_flips(const Edge& edg);

    /**
     * @brief Borrow the sorted event times of an update tuple.
     * @param t_index Plaquette index in the x-basis, or star center in the z-basis.
     * @return Read-only view; reacquire it after mutating this tuple's history.
     */
    std::span<const double> get_tuple_spin_flips(int t_index);

    /**
     * @brief Overwrite one time in the edge's full history.
     * @pre The index is valid and the replacement preserves sorted order.
     * @note Updates only this history; the caller must synchronize related histories.
     */
    inline void set_spin_flip_imag_time(
        const Edge& edg, int spin_flip_index, double imag_time
    );

    /**
     * @brief Overwrite one time in the edge's single-spin history.
     * @pre The index is valid and the replacement preserves sorted order.
     * @note Updates only this history; the caller must synchronize related histories.
     */
    inline void set_single_spin_flip_imag_time(
        const Edge& edg, int spin_flip_index, double imag_time
    );

    /**
     * @brief Remove two ordered single-spin times from both edge histories.
     * @pre Both events exist and imag_time_single_spin_flip < imag_time_next_single_spin_flip.
     * @throws std::runtime_error If a requested event is absent.
     */
    void delete_double_single_spin_flip(
        const Edge& edg, double imag_time_single_spin_flip, double imag_time_next_single_spin_flip
    );

    /**
     * @brief Remove a single-spin event from both of an edge's histories.
     * @param spin_flip_index Index in the full history, not the single-spin history.
     * @pre The indexed event is a single-spin event.
     * @throws std::runtime_error If the corresponding single-spin event is absent.
     */
    void delete_single_spin_flip(const Edge& edg, int spin_flip_index);

    /**
     * @brief Remove two tuple events from the tuple and its edge histories.
     * @param tuple_index Plaquette index in x, or star center in z.
     * @param tuple_edges Edges belonging to tuple_index.
     * @pre The requested times exist in increasing order in every affected history.
     * @throws std::runtime_error If a requested event is absent.
     */
    void delete_double_tuple_flip(
        int tuple_index, 
        std::span<const Edge> tuple_edges, 
        double imag_time_tuple_flip, 
        double imag_time_next_tuple_flip
    );

    /**
     * @brief Remove one tuple event from the tuple and its edge histories.
     * @param tuple_index Plaquette index in x, or star center in z.
     * @param tuple_edges Edges belonging to tuple_index.
     * @pre The requested time exists in every affected history.
     * @throws std::runtime_error If a requested event is absent.
     */
    void delete_tuple_flip(
        int tuple_index, std::span<const Edge> tuple_edges, double imag_time_tuple_flip
    );

    /**
     * @brief Test whether every tuple edge has an event at exactly tau.
     * @note Checks full edge histories, not the tuple's separate event list.
     *       Coincident single-spin events also satisfy this test.
     */
    bool check_tuple_flip_present_tuple(std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Test adjacent plaquettes for coincident events on all their edges.
     * @see check_tuple_flip_present_tuple() for the event-presence criterion.
     */
    bool check_plaquette_flip_at_edge(const Edge& edg, double tau);

    /**
     * @brief Test both endpoint stars for coincident events on all their edges.
     * @see check_tuple_flip_present_tuple() for the event-presence criterion.
     */
    bool check_star_flip_at_edge(const Edge& edg, double tau);

    /** @brief Test whether any tuple edge has an event at exactly tau. */
    bool check_spin_flips_present_tuple(std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Find the neighboring event strictly after tau, wrapping at beta.
     * @return A time in the stored period; tau if the edge has no events.
     *         A history containing only tau also returns tau after wrapping.
     */
    double flip_next_imag_time(const Edge& edg, double tau);

    /** @brief Apply flip_next_imag_time() to each edge, preserving tuple_edges order. */
    std::vector<double> flip_next_imag_times_tuple(std::span<const Edge> tuple_edges, double tau);
    /** @brief Small-vector version of flip_next_imag_times_tuple(). */
    SmallEnergyVector flip_next_imag_times_tuple_small(std::span<const Edge> tuple_edges, double tau);
    /**
     * @brief Intersect neighboring-event windows around a tuple event at tau.
     * @param tuple_edges Edges of the tuple.
     * @param tau Event time present on every edge for a tuple-move proposal.
     * @param flip_indices Optional output ranks in full edge histories, in input order.
     * @return (left, right) modulo beta; a wrapped interval can have left > right.
     * @note Edges without events constrain the window to tau and add no rank.
     *       The optional ranks are complete only when every edge has a history.
     */
    std::pair<double, double> tuple_flip_window(
        std::span<const Edge> tuple_edges,
        double tau,
        SmallIndexVector* flip_indices = nullptr
    );

    /**
     * @brief Find the neighboring event strictly before tau, wrapping at beta.
     * @return A time in the stored period; tau if the edge has no events.
     *         A history containing only tau also returns tau after wrapping.
     */
    double flip_prev_imag_time(const Edge& edg, double tau);

    /** @brief Apply flip_prev_imag_time() to each edge, preserving tuple_edges order. */
    std::vector<double> flip_prev_imag_times_tuple(std::span<const Edge> tuple_edges, double tau);
    /** @brief Small-vector version of flip_prev_imag_times_tuple(). */
    SmallEnergyVector flip_prev_imag_times_tuple_small(std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Insert two times into both an edge's full and single-spin histories.
     * @pre 0 <= tau_left < tau_right <= beta; the caller prevents unwanted coincidences.
     */
    void insert_double_single_spin_flip(const Edge& edg, double tau_left, double tau_right);

    /**
     * @brief Insert tau into both an edge's full and single-spin histories.
     * @pre tau is in the stored time period; the caller prevents unwanted coincidences.
     */
    void insert_single_spin_flip(const Edge& edg, double tau);

    /**
     * @brief Insert two events into the tuple history and every tuple edge's full history.
     * @param tuple_index Plaquette index in x, or star center in z.
     * @param tuple_edges Edges belonging to tuple_index.
     * @pre 0 <= tau_left < tau_right <= beta; the caller prevents unwanted coincidences.
     * @throws std::invalid_argument If tau_right < tau_left.
     */
    void insert_double_tuple_flip(
        int tuple_index, std::span<const Edge> tuple_edges, double tau_left, double tau_right
    );

    /**
     * @brief Insert an event into the tuple history and every tuple edge's full history.
     * @param tuple_index Plaquette index in x, or star center in z.
     * @param tuple_edges Edges belonging to tuple_index.
     * @pre tau lies in the stored time period; the caller prevents unwanted coincidences.
     */
    void insert_tuple_flip(int tuple_index, std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Move a single-spin event in both edge histories.
     * @param spin_flip_index Index of the single event in the full edge history.
     * @param tau_new New time within the neighboring-event window.
     * @param no_move_over_beta False when the event crosses the time origin.
     * @pre No other event is crossed, except through the periodic time boundary.
     * @note The caller reverses the time-origin spin for a boundary crossing and
     *       applies the corresponding energy-cache changes.
     */
    void move_spin_flip(const Edge& edg, int spin_flip_index, double tau_new, bool no_move_over_beta);

    /**
     * @brief Move a tuple event in its history and all affected edge histories.
     * @param tuple_index Plaquette index in x, or star center in z.
     * @param tuple_edges Edges belonging to tuple_index.
     * @param tau_old Existing tuple-event time.
     * @param tau_new New time within the common neighboring-event window.
     * @param no_move_over_beta False when the event crosses the time origin.
     * @param edge_flip_indices Optional full-history ranks, aligned with tuple_edges.
     * @param tuple_flip_index Optional tuple-history rank; -1 requests a search.
     * @pre The event exists and crosses no other event except through the time boundary.
     * @throws std::invalid_argument If a nonempty edge_flip_indices has the wrong size.
     * @note The caller reverses time-origin spins for boundary crossings and updates caches.
     */
    void move_tuple_flip(
        int tuple_index, 
        std::span<const Edge> tuple_edges, 
        double tau_old, 
        double tau_new, 
        bool no_move_over_beta,
        std::span<const int> edge_flip_indices = {},
        int tuple_flip_index = -1
    );

    /** @brief Print every edge's spin at the time origin to standard output. */
    void print_spins();

    /** @brief Print the full event history of one edge to standard output. */
    void print_spin_flip_imag_times(const Edge& edg);

    /** @brief Print the full event history of each tuple edge to standard output. */
    void print_tuple_flip_imag_times(std::span<const Edge> tuple_edges);

    /** @brief Count all single and tuple events in an edge's full history. */
    inline int get_spin_flip_count(const Edge& edg);

    /** @brief Reverse one edge's spin at the time origin; histories and caches are unchanged. */
    void flip_spin(const Edge& edg);

    /** @brief Reverse the time-origin spins incident to vertex v; caches are unchanged. */
    void flip_star(int v);
    /**
     * @brief Return the edge connecting two vertices.
     * @pre The edge exists; only debug builds throw std::runtime_error if it is absent.
     */
    inline Edge edge_in_between(int v_1, int v_2);
    /** @brief Test whether two vertices share an edge. */
    inline bool exists_edge(int v_1, int v_2);
    /** @brief Return the source and target vertex indices of an edge. */
    inline std::pair<int, int> vertices_of_edge(const Edge& edg);

    /** @brief Draw a uniform edge and return (descriptor, source vertex, target vertex). */
    std::tuple<Edge, int, int> get_random_edge();
    /** @brief Draw a uniform edge descriptor using the shared RNG. */
    Edge get_random_edge_descriptor();

    /** @brief Draw a vertex index uniformly using the shared RNG. */
    int get_random_vertex();

    /** @brief Draw a plaquette index uniformly using the shared RNG. */
    int get_random_plaquette_index();
    /** @brief Copy the vertex pairs defining a plaquette's edges, in construction order. */
    std::vector<std::pair<int,int>> get_plaquette_vertex_pairs(int p_index);
    /** @brief Copy the vertex indices defining an elementary cube. */
    std::vector<int> get_cube_vertices(int c_index);
    /** @brief Borrow cached plaquette edges in construction order; valid while the geometry lives. */
    inline std::span<const Edge> get_plaquette_edges(int p_index);
    /** @brief Borrow cached incident edges at a star center; valid while the geometry lives. */
    inline std::span<const Edge> get_star_edges(int center_index);
    /** @brief Sum the tuple's spins at the time origin. */
    int get_tuple_sum(std::span<const Edge> tuple_edges);
    /** @brief Multiply the tuple's spins at the time origin. */
    int get_tuple_prod(std::span<const Edge> tuple_edges);
    /** @brief Reverse the tuple's time-origin spins; histories and caches are unchanged. */
    void flip_tuple(std::span<const Edge> tuple_edges);
    /** @brief Multiply all time-origin spins incident to vertex v. */
    int get_vertex_nn_spins_prod(int v);

    /**
     * @brief Return -2 times the bare tuple integral for a spin-product reversal.
     * @pre imag_time_1 < imag_time_2; split wrapped intervals at beta.
     */
    double integrated_tuple_energy_diff(
        std::span<const Edge> tuple_edges, 
        double imag_time_1, 
        double imag_time_2
    );

    /**
     * @brief Integrate a tuple-product change from a local combination schedule.
     * @param spin_flip_lookup Sorted (time, type) events: exactly one tuple event
     *        tagged 1, plus one single event tagged 0 for each overlapping edge.
     * @note Equal-time events combine by parity. The tuple event toggles the local
     *       product only for odd overlap. imag_time_2 is the history-integration
     *       cutoff; later schedule intervals use the product at that cutoff.
     */
    inline double integrated_tuple_energy_diff_combination(
        std::span<const Edge> tuple_edges, 
        double imag_time_1, 
        double imag_time_2, 
        const std::vector<std::pair<double, int>>& spin_flip_lookup
    );

    /**
     * @brief Integrate the product of tuple spins without a minus sign or coupling.
     * @pre Use ordered, non-wrapping bounds within [0, beta].
     * @return Bare integral; zero if imag_time_1 >= imag_time_2.
     */
    inline double integrated_tuple_energy(
        std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2
    );

    /**
     * @brief Bare star-integral changes when one edge is reversed on a time interval.
     * @param total_cache Use -2 * cached integrals only for a full-period flip.
     * @return (sum, affected indices, aligned per-star changes), ordered by
     *         source and target vertices. No coupling or Hamiltonian minus sign is applied.
     * @pre imag_time_1 < imag_time_2; split wrapped intervals at beta.
     * @throws std::invalid_argument If the bounds are equal.
     */
    std::tuple<double, SmallIndexVector, SmallEnergyVector> 
    integrated_star_energy_diff(
        const Edge& edg, double imag_time_1, double imag_time_2, bool total_cache
    );

    /**
     * @brief Bare plaquette-integral changes when one edge is reversed on a time interval.
     * @param total_cache Use -2 * cached integrals only for a full-period flip.
     * @return (sum, affected indices, aligned per-plaquette changes), ordered by
     *         the edge's plaquette adjacency list. No coupling or Hamiltonian minus sign is applied.
     * @pre imag_time_1 < imag_time_2; split wrapped intervals at beta.
     * @throws std::invalid_argument If the bounds are equal.
     */
    std::tuple<double, SmallIndexVector, SmallEnergyVector> 
    integrated_plaquette_energy_diff(
        const Edge& edg, double imag_time_1, double imag_time_2, bool total_cache
    );

    /** @brief Recompute the sum of bare star integrals over [0, beta]. */
    double total_integrated_star_energy();

    /** @brief Recompute the sum of bare plaquette integrals over [0, beta]. */
    double total_integrated_plaquette_energy();

    /**
     * @brief Bare star-integral changes from one tuple event and per-edge single events.
     * @param plaquette_index Index of the update tuple.
     * @param spin_flip_lookup One single-event time per edge, aligned with get_plaquette_edges(plaquette_index).
     * @param imag_time_tuple_flip Time of the update tuple's event.
     * @return (sum, sorted unique star indices, aligned bare changes).
     * @pre imag_time_1 < imag_time_2 and the per-edge array matches the tuple size.
     * @throws std::invalid_argument If the bounds are equal.
     * @see integrated_tuple_energy_diff_combination() for schedule parity and cutoff rules.
     */
    std::tuple<double, SmallIndexVector, SmallEnergyVector> 
    integrated_star_energy_diff_combination(
        int plaquette_index, 
        double imag_time_1, 
        double imag_time_2, 
        std::span<const double> spin_flip_lookup, 
        double imag_time_tuple_flip
    );

    /**
     * @brief Bare plaquette-integral changes from one tuple event and per-edge single events.
     * @param star_index Index of the update tuple.
     * @param spin_flip_lookup One single-event time per edge, aligned with get_star_edges(star_index).
     * @param imag_time_tuple_flip Time of the update tuple's event.
     * @return (sum, sorted unique plaquette indices, aligned bare changes).
     * @pre imag_time_1 < imag_time_2 and the per-edge array matches the tuple size.
     * @throws std::invalid_argument If the bounds are equal.
     * @see integrated_tuple_energy_diff_combination() for schedule parity and cutoff rules.
     */
    std::tuple<double, SmallIndexVector, SmallEnergyVector> 
    integrated_plaquette_energy_diff_combination(
        int star_index, 
        double imag_time_1, 
        double imag_time_2, 
        std::span<const double> spin_flip_lookup, 
        double imag_time_tuple_flip
    );

    /**
     * @brief Return -2 times the bare edge integral for a spin reversal.
     * @pre imag_time_1 < imag_time_2; split wrapped intervals at beta.
     * @throws std::invalid_argument If the bounds are equal.
     */
    double integrated_edge_energy_diff(const Edge& edg, double imag_time_1, double imag_time_2);
    /**
     * @brief Compute a bare edge-integral change on an interval with no inner events.
     * @param known_flip_index Optional lower-bound rank of known_flip_time; -1 searches.
     * @param known_flip_time Event bordering the interval whose rank is already known.
     * @pre imag_time_1 < imag_time_2, with no event strictly inside the interval.
     * @throws std::invalid_argument If the bounds are equal.
     */
    inline double integrated_edge_energy_diff_no_inner_flips(
        const Edge& edg, double imag_time_1, double imag_time_2,
        int known_flip_index = -1, double known_flip_time = 0.
    );

    /**
     * @brief Integrate the bare edge change from a sorted proposed flip schedule.
     * @param spin_flip_lookup Sorted (time, type) pairs; only times are used here.
     * @pre imag_time_1 < imag_time_2; the change has even parity before the first
     *      schedule event in the interval. Events outside [imag_time_1, imag_time_2)
     *      are ignored. Equal-time events cancel by parity.
     * @throws std::invalid_argument If the bounds are equal.
     */
    inline double integrated_edge_energy_diff_combination(
        const Edge& edg, 
        double imag_time_1, 
        double imag_time_2, 
        std::vector<std::pair<double,int>>& spin_flip_lookup
    );
    /**
     * @brief Two-event overload; event times may be supplied in either order.
     * @see integrated_edge_energy_diff_combination() for interval and parity rules.
     */
    inline double integrated_edge_energy_diff_combination(
        const Edge& edg,
        double imag_time_1,
        double imag_time_2,
        double imag_time_flip_1,
        double imag_time_flip_2
    );

    /**
     * @brief Integrate one spin without a Hamiltonian minus sign or coupling.
     * @pre imag_time_1 < imag_time_2 within [0, beta]; split wrapped intervals.
     * @throws std::invalid_argument If the bounds are equal.
     */
    inline double integrated_edge_energy(const Edge& edg, double imag_time_1, double imag_time_2);

    /**
     * @brief Integrate an edge spin with weight min(tau, beta - tau).
     * @pre Bounds lie in [0, beta]; reversed bounds request a wrap through beta.
     * @throws std::invalid_argument If the bounds are equal.
     */
    double integrated_edge_energy_weighted(const Edge& edg, double imag_time_1, double imag_time_2);

    /** @brief Recompute the sum of bare edge integrals over [0, beta]. */
    double total_integrated_edge_energy();

    /** @brief Sum all edge integrals weighted by min(tau, beta - tau) over [0, beta]. */
    double total_integrated_edge_energy_weighted();

    /**
     * @brief Rebuild bare edge caches and the basis's diagonal tuple caches.
     * Stars are initialized in x and plaquettes in z. Call after direct history
     * edits or before reusing caches whose coupling was previously zero.
     */
    void init_potential_energy();

    /** @brief Sum spins at the time origin, without a minus sign or field coupling. */
    double get_diag_single_energy();

    /**
     * @brief Pack (integrated magnetization per edge, time-origin magnetization per edge).
     * @note The integral is not divided by beta. Both components are real estimators.
     */
    std::complex<double> get_diag_M_M();

    /**
     * @brief Pack (half the weighted magnetization integral per edge, magnetization per edge).
     * Uses the weight min(tau, beta - tau) and the time-origin magnetization.
     */
    std::complex<double> get_diag_dynamical_M_M();

    /** @brief Pack the total single-spin event count into both real and imaginary components. */
    std::complex<double> get_non_diag_M_M();

    /**
     * @brief Pack single-spin event counts in [0, beta/2) and [beta/2, beta).
     * @return Real part kL and imaginary part kR for the off-diagonal dynamical estimator.
     * @note The counts use the current time origin; this method does not rotate it.
     */
    std::complex<double> get_kL_kR_single();

    /** @brief Return the single-spin event count / beta; its negative estimates the lmbda energy term. */
    double get_non_diag_single_energy_x();

    /** @brief Return the single-spin event count / beta; its negative estimates the h energy term. */
    double get_non_diag_single_energy_z();

    /** @brief Return the plaquette event count / beta; its negative estimates the J energy term. */
    double get_non_diag_tuple_energy_x();

    /** @brief Return the star event count / beta; its negative estimates the mu energy term. */
    double get_non_diag_tuple_energy_z();

    /** @brief Sum star spin products at the time origin, without a minus sign or mu. */
    double get_diag_tuple_energy_x();

    /** @brief Sum plaquette spin products at the time origin, without a minus sign or J. */
    double get_diag_tuple_energy_z();

    /**
     * @brief Pack equal-time half and full Wilson/'t Hooft loop products.
     * @return Real part: half-loop product; imaginary part: full-loop product.
     *         The statistics layer forms the Fredenhagen-Marcu ratio from their means.
     * @see https://doi.org/10.1103/PhysRevLett.56.223.
     */
    std::complex<double> fredenhagen_marcu();

    /**
     * @brief Return the alternating gap sum / beta for a uniformly chosen plaquette.
     * @see https://doi.org/10.1103/PhysRevB.85.195104 for the related order parameter.
     */
    double get_staggered_imaginary_times_plaquette();

    /**
     * @brief Return the alternating gap sum / beta for a uniformly chosen star.
     * @see https://doi.org/10.1103/PhysRevB.85.195104 for the related order parameter.
     */
    double get_staggered_imaginary_times_star();

    /**
     * @brief Detect winding of a connected spin -1 edge cluster using lattice coordinates.
     * @note Depends on geometry-specific coordinates and boundary conventions.
     */
    bool is_winding_percolating();

    /**
     * @brief Detect a spin -1 edge cluster spanning opposite coordinate boundaries.
     * @note Depends on geometry-specific coordinates and boundary conventions.
     */
    bool is_percolating();

    /**
     * @brief Detect winding of plaquettes connected through spin -1 edges.
     * @note Depends on geometry-specific coordinates and boundary conventions.
     */
    bool is_winding_plaquette_percolating();

    /**
     * @brief Detect a plaquette cluster spanning opposite boundaries through spin -1 edges.
     * @note Depends on geometry-specific coordinates and boundary conventions.
     */
    bool is_plaquette_percolating();

    /**
     * @brief Detect winding of cubes connected through plaquettes with spin product +1.
     * @pre Use a cubic lattice with periodic boundaries.
     * @throws std::invalid_argument If lattice dimensionality is not three.
     */
    bool is_winding_cube_percolating();

    /** @brief Count spin -1 edges in the largest connected string cluster. */
    int largest_cluster();

    /** @brief Count plaquettes in the largest cluster connected through spin -1 edges. */
    int largest_plaquette_cluster();

    /**
     * @brief Return largest_cluster() / get_edge_count() if any string cluster percolates.
     * @return Fraction of all edges in the largest string cluster, or zero without percolation.
     * @note Periodic boundaries use winding; open boundaries use spanning.
     */
    double percolation_strength();

    /**
     * @brief Return the string percolation indicator (zero or one) for this configuration.
     * @note Uses winding for periodic boundaries and spanning for open boundaries.
     *       A probability is obtained by averaging this indicator over samples.
     */
    double percolation_probability();

    /**
     * @brief Return largest_plaquette_cluster() / get_plaquette_count() if percolating.
     * @return Fraction of all plaquettes in the largest cluster, or zero without percolation.
     * @note Periodic boundaries use winding; open boundaries use spanning.
     */
    double plaquette_percolation_strength();

    /**
     * @brief Return the plaquette percolation indicator (zero or one) for this configuration.
     * @note Uses winding for periodic boundaries and spanning for open boundaries.
     *       A probability is obtained by averaging this indicator over samples.
     */
    double plaquette_percolation_probability();

    /** @brief Return -1; cube percolation strength is not implemented. */
    double cube_percolation_strength();

    /**
     * @brief Return the cube winding indicator for periodic boundaries.
     * @return One or zero for periodic cubic lattices; -1 for unsupported open boundaries.
     * @throws std::invalid_argument For periodic lattices that are not three-dimensional.
     */
    double cube_percolation_probability();

    /**
     * @brief Shift all histories by one uniformly drawn time origin modulo beta.
     * Adjusts time-origin spins by the parity of events crossing the cut and keeps
     * histories sorted. Full-period bare integrals are invariant under this shift.
     */
    void rotate_imag_time();

    /**
     * @brief Append time-origin spins to a temporary binary snapshot spool.
     * @throws std::runtime_error If the temporary spool cannot be created or written.
     * @see write_graph() to export the accumulated snapshots.
     */
    void update_spin_string();

    /**
     * @brief Write geometry and serialized spin histories as GraphML in file_name.xml.
     * If snapshots have been spooled, export their per-edge spin strings and reset
     * the spool after a successful write. Otherwise export the graph's legacy
     * spin_string properties, which may be empty when no history has been recorded.
     * @param file_name Filename stem, without the .xml extension.
     * @param output_directory Existing writable destination directory.
     */
    void write_graph(const std::string& file_name, const std::filesystem::path& output_directory);

private:   
    char BASIS;
    int LATTICE_DIMENSIONALITY;
    std::string LATTICE_TYPE;   
    int SYSTEM_SIZE;  
    double BETA; 
    std::string BOUNDARIES;
    int DEFAULT_SPIN;

    // Graph vertices carry stars; graph edges carry spin histories.
    LatticeGraph g;
    // Plaquette index -> boundary vertex pairs, in construction order.
    std::vector<std::vector<std::pair<int,int>>> plaquette_vector;
    // Plaquette index -> sorted tuple events (off-diagonal in x).
    std::vector<std::vector<double>> plaquette_flip_vector;
    // Plaquette index -> cached bare spin-product integral (diagonal in z).
    std::vector<double> integrated_plaquette_energy_vector;
    // Plaquette coordinates used by percolation.
    std::vector<double> plaquette_x_vector;
    std::vector<double> plaquette_y_vector;
    std::vector<double> plaquette_z_vector;
    // Cube index -> vertices (cubic geometry only).
    std::vector<std::vector<int>> cube_vector;
    // Cube coordinates used by percolation.
    std::vector<double> cube_x_vector;
    std::vector<double> cube_y_vector;
    std::vector<double> cube_z_vector;
    // Plaquette index -> adjacent cubes.
    std::vector<std::vector<int>> plaquette_part_of_cube_lookup;
    // Cube index -> boundary plaquettes.
    std::vector<std::vector<int>> cube_has_plaquettes_lookup;

    // p -> edges of plaquette p (arbitrary length: 3/4/6/…)
    std::vector<std::vector<Edge>> plaquette_edges_cache_;
    // p -> unique star centers touched by plaquette p
    std::vector<std::vector<int>> plaquette_vertices_cache_;
    // v -> incident edges (star at vertex v)
    std::vector<std::vector<Edge>> star_edges_cache_;
    // v -> unique plaquettes touching star v
    std::vector<std::vector<int>> star_plaquettes_cache_;
    // Edge descriptors in graph iteration order.
    std::vector<Edge> egde_cache_;

    /**
     * @brief Own the temporary snapshot file; destruction removes it.
     * Copying creates an empty spool so lattice copies never share output streams.
     */
    struct SnapshotSpoolState {
        std::filesystem::path path{};
        std::ofstream stream{};
        std::size_t edge_count = 0;
        std::size_t sample_count = 0;

        /** @brief Create an empty spool without opening a temporary file. */
        SnapshotSpoolState() = default;
        /** @brief Create an empty spool when copying; recorded snapshots are not copied. */
        SnapshotSpoolState(const SnapshotSpoolState&)
            : path{}, stream{}, edge_count(0), sample_count(0) {}
        /** @brief Reset this spool without copying the source's file or snapshots. */
        SnapshotSpoolState& operator=(const SnapshotSpoolState&) {
            reset();
            return *this;
        }
        /** @brief Close the stream and attempt to remove its temporary file. */
        ~SnapshotSpoolState();

        /** @brief Test whether a temporary file path has been assigned to the spool. */
        bool active() const { return !path.empty(); }
        /**
         * @brief Close the stream, remove the temporary file, and clear the path and counts.
         * @note File-removal errors are ignored.
         */
        void reset();
    };
    SnapshotSpoolState snapshot_spool_;

    // These two vectors are used to calculate the Fredenhagen-Marcu order parameter
    std::vector<VertexPair> half_path_vector;
    std::vector<VertexPair> full_path_vector;

    std::vector<double> MAX_COORDINATES;
    std::vector<double> MAX_PLAQUETTE_COORDINATES;

    /** @brief Throw std::invalid_argument for invalid lattice specifications. */
    void check_input_validity() const;
    /** @brief Build geometry lookups and edge descriptors after graph construction. */
    void build_caches_();
    /** @brief Lazily create and open the temporary snapshot file. */
    void ensure_snapshot_spool_();
    /** @brief Transpose snapshot rows into GraphML edge histories and release the spool. */
    void write_snapshot_graphml_from_spool_(const std::string& file_name, const std::filesystem::path& output_directory);

    /** @brief Instantiate the graph and QMC storage from reusable geometry mappings. */
    void init_lattice_graph();
    /**
     * @brief Integrate a bare diagonal tuple product using only single-spin histories.
     * @pre Tuple events leave this product unchanged, as for stars in x or plaquettes in z.
     * @note Stars and plaquettes overlap on an even number of edges, so only single
     *       events need to be merged when evaluating a diagonal tuple integral.
     * @see integrated_tuple_energy() for interval and return-value conventions.
     */
    inline double integrated_tuple_energy_single_flips(
        std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2
    );
    /**
     * @brief Return -2 times the bare diagonal tuple integral for a product reversal.
     * @see integrated_tuple_energy_single_flips() for history and interval requirements.
     */
    inline double integrated_tuple_energy_diff_single_flips(
        std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2
    );
    /**
     * @brief Integrate a combination-update change using full or single-spin histories.
     * @param single_flips_only True selects single-spin histories; false selects all flips.
     * @pre When single_flips_only is true, tuple events leave the local product unchanged.
     * @return Bare tuple-integral change; zero for an empty proposal schedule.
     * @see integrated_tuple_energy_diff_combination() for schedule parity and cutoff rules.
     */
    inline double integrated_tuple_energy_diff_combination_from_flips(
        std::span<const Edge> tuple_edges,
        double imag_time_1,
        double imag_time_2,
        const std::vector<std::pair<double, int>>& spin_flip_lookup,
        bool single_flips_only
    );
    /**
     * @brief Integrate a bare tuple product by merging the selected edge histories.
     * @param single_flips_only True selects single-spin histories; false selects all flips.
     * @pre Use ordered, non-wrapping bounds within [0, beta]. When single_flips_only
     *      is true, tuple events must leave the product unchanged.
     * @return Bare integral; zero if imag_time_1 >= imag_time_2.
     * @note Equal-time flips combine by parity.
     */
    inline double integrated_tuple_energy_from_flips(
        std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2, bool single_flips_only
    );
    std::vector<std::pair<Vertex, Vertex>> edge_vector;
    mutable std::uniform_int_distribution<int> edge_dist;
    mutable std::uniform_int_distribution<int> vertex_dist;
    mutable std::uniform_int_distribution<int> plaquette_dist;

    // Lattice copies share this generator; copying the RNG object itself reseeds it.
    std::shared_ptr<RNG> rng;
    std::uniform_real_distribution<double> uniform_dist{0., 1.};

    /**
     * @brief Choose a uniform iterator from a nonempty range [start, end).
     * @param gen Generator used for the draw.
     * @return Iterator advanced to the selected element.
     */
    template<typename Iter, typename RandomGenerator>
    Iter random_element(Iter start, Iter end, RandomGenerator& gen) {
        std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
        std::advance(start, dis(gen));
        return start;
    }

    /**
     * @brief Choose an iterator using a lazily seeded generator local to this overload.
     * @pre [start, end) is nonempty.
     * @note This overload does not use the lattice's shared RNG.
     */
    template<typename Iter>
    Iter random_element(Iter start, Iter end) {
        static std::random_device rd;
        static RNG gen(rd());
        return random_element(start, end, gen);
    }

    /** @brief Wrap a into [0, b); b must be positive. */
    template<typename T>
    requires std::integral<T> || std::floating_point<T>
    constexpr T modulo(T a, T b) {
        if constexpr (std::integral<T>) {
            T result = a % b;
            return result >= 0 ? result : result + b;
        } else {
            T result = std::fmod(a, b);
            return result >= 0 ? result : result + b;
        }
    }
};

/** @brief Count edges with spin +1 at the time origin. */
inline int Lattice::get_non_string_count() {
    int result = 0;
    for (const auto& edg : egde_cache_) {
        if (g[edg].spin == 1) result += 1;
    }
    return result;
}

/** @brief Count edges with spin -1 at the time origin. */
inline int Lattice::get_string_count() {
    int result = 0;
    for (const auto& edg : egde_cache_) {
        if (g[edg].spin == -1) result += 1;
    }
    return result;
}

/** @brief Return the number of graph vertices. */
inline int Lattice::get_vertex_count() {
    return boost::num_vertices(g);
}

/** @brief Return the number of graph edges. */
inline int Lattice::get_edge_count() {
    return boost::num_edges(g);
}

/** @brief Return the number of elementary plaquettes. */
inline int Lattice::get_plaquette_count() {
    return plaquette_vector.size();
}

/** @brief Return the number of elementary cubes. */
inline int Lattice::get_cube_count() {
    return cube_vector.size();
}
/** @brief Return the edge spin at the time origin. */
inline int Lattice::get_spin(const Edge& edg) {
    return g[edg].spin;
}

/** @brief Read the cached bare edge integral over [0, beta]. */
inline double Lattice::get_potential_edge_energy(const Edge& edg) {
    return g[edg].integrated_edge_energy;
}

/** @brief Replace the cached bare edge integral. */
inline void Lattice::set_potential_edge_energy(const Edge& edg, double potential_energy) {
    g[edg].integrated_edge_energy = potential_energy;
}

/** @brief Add a bare integral change to the edge cache. */
inline void Lattice::add_potential_edge_energy(const Edge& edg, double diff) {
    g[edg].integrated_edge_energy += diff;
}

/** @brief Read the cached bare star integral (maintained in the x-basis). */
inline double Lattice::get_potential_star_energy(int star_index) {
    return g[star_index].integrated_star_energy;
}
    
/** @brief Replace the cached bare star integral. */
inline void Lattice::set_potential_star_energy(int star_index, double potential_energy) {
    g[star_index].integrated_star_energy = potential_energy;
}

/** @brief Add a bare integral change to the star cache. */
inline void Lattice::add_potential_star_energy(int star_index, double diff) {
    g[star_index].integrated_star_energy += diff;
}

/** @brief Read the cached bare plaquette integral (maintained in the z-basis). */
inline double Lattice::get_potential_plaquette_energy(int plaquette_index) {
    return integrated_plaquette_energy_vector[plaquette_index];
}
    
/** @brief Replace the cached bare plaquette integral. */
inline void Lattice::set_potential_plaquette_energy(int plaquette_index, double potential_energy) {
    integrated_plaquette_energy_vector[plaquette_index] = potential_energy;
}

/** @brief Add a bare integral change to the plaquette cache. */
inline void Lattice::add_potential_plaquette_energy(int plaquette_index, double diff) {
    integrated_plaquette_energy_vector[plaquette_index] += diff;
}

/** @brief Return the edge's geometry-specific direction label. */
inline std::string Lattice::get_orientation(const Edge& edg) {
    return g[edg].orientation;
}

/** @brief Return a time by index in the edge's full flip history. */
inline double Lattice::get_spin_flip_imag_time(const Edge& edg, int spin_flip_index) {
    return g[edg].spin_flips[spin_flip_index];
}

/** @brief Overwrite one time in the edge's full history. */
inline void Lattice::set_spin_flip_imag_time(const Edge& edg, int spin_flip_index, double imag_time) {
    g[edg].spin_flips[spin_flip_index] = imag_time;
}

/** @brief Overwrite one time in the edge's single-spin history. */
inline void Lattice::set_single_spin_flip_imag_time(const Edge& edg, int spin_flip_index, double imag_time) {
    g[edg].single_spin_flips[spin_flip_index] = imag_time;
}

/** @brief Count all single and tuple events in an edge's full history. */
inline int Lattice::get_spin_flip_count(const Edge& edg) {
    return g[edg].spin_flips.size();
}

/** @brief Return the edge connecting two vertices. */
inline Lattice::Edge Lattice::edge_in_between(int v_1, int v_2) {
    const auto edg_full = boost::edge(v_1, v_2, g);
#ifndef NDEBUG
    if (!edg_full.second) {
        throw std::runtime_error(std::format("There is no edge between vertex {} and {}.", v_1, v_2));
    }
#endif
    return edg_full.first;
}

/** @brief Borrow cached plaquette edges in construction order; valid while the geometry lives. */
inline std::span<const Lattice::Edge> Lattice::get_plaquette_edges(int p_index) {
    const auto& v = plaquette_edges_cache_[static_cast<size_t>(p_index)];
    return {v.data(), v.size()};
}

/** @brief Borrow cached incident edges at a star center; valid while the geometry lives. */
inline std::span<const Lattice::Edge> Lattice::get_star_edges(int center_index) {
    const auto& s = star_edges_cache_[static_cast<size_t>(center_index)];
    return {s.data(), s.size()};
}

/** @brief Test whether two vertices share an edge. */
inline bool Lattice::exists_edge(int v_1, int v_2) {
    return boost::edge(v_1, v_2, g).second;
}

/** @brief Return the source and target vertex indices of an edge. */
inline std::pair<int, int> Lattice::vertices_of_edge(const Edge& edg) {
    const Edge e = edg;
    const auto source_v = boost::source(e, g);
    const auto target_v = boost::target(e, g);
    return {source_v, target_v};
}

/** @brief Two-event overload; event times may be supplied in either order. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_edge_energy_diff_combination(
    const Edge& edg,
    double imag_time_1,
    double imag_time_2,
    double imag_time_flip_1,
    double imag_time_flip_2
) {
    if (imag_time_1 == imag_time_2) {
        throw std::invalid_argument(
            "integrated_edge_energy_diff_combination: time interval must be non-zero");
    }

    double t_first = imag_time_flip_1;
    double t_second = imag_time_flip_2;
    if (t_second < t_first) {
        std::swap(t_first, t_second);
    }

    const bool first_active = (t_first >= imag_time_1) && (t_first < imag_time_2);
    const bool second_active = (t_second >= imag_time_1) && (t_second < imag_time_2);
    if (!first_active && !second_active) {
        return 0.0;
    }
    if (first_active && second_active) {
        if (t_first == t_second) {
            return 0.0;
        }
        return -2.0 * integrated_edge_energy(edg, t_first, t_second);
    }

    const double t_start = first_active ? t_first : t_second;
    return -2.0 * integrated_edge_energy(edg, t_start, imag_time_2);
}

/** @brief Compute a bare edge-integral change on an interval with no inner events. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_edge_energy_diff_no_inner_flips(
    const Edge& edg, double imag_time_1, double imag_time_2,
    int known_flip_index, double known_flip_time
) {
    if (imag_time_1 == imag_time_2) {
        throw std::invalid_argument(
            "integrated_edge_energy_diff_no_inner_flips: time interval must be non-zero");
    }

    return -2.0 * detail::spin_integral_no_inner_flips(
        {get_spin(edg), g[edg].spin_flips}, imag_time_1, imag_time_2,
        known_flip_index, known_flip_time
    );
}

/** @brief Integrate the bare edge change from a sorted proposed flip schedule. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_edge_energy_diff_combination(
    const Edge& edg, 
    double imag_time_1, 
    double imag_time_2, 
    std::vector<std::pair<double,int>>& spin_flip_lookup
) {
    if (imag_time_1 == imag_time_2) {
        throw std::invalid_argument(
            "integrated_edge_energy_diff_combination: time interval must be non-zero");
    }

    if (spin_flip_lookup.empty()) {
        return 0.0;
    }
    if (spin_flip_lookup.size() == 2) {
        return integrated_edge_energy_diff_combination(
            edg,
            imag_time_1,
            imag_time_2,
            spin_flip_lookup[0].first,
            spin_flip_lookup[1].first
        );
    }

    const auto* delta_it = spin_flip_lookup.data();
    const auto* delta_end = delta_it + spin_flip_lookup.size();

    bool odd_parity = false;
    while (delta_it != delta_end && delta_it->first < imag_time_1) {
        ++delta_it;
    }

    double odd_integral = 0.0;
    double t_prev = imag_time_1;
    while (delta_it != delta_end && delta_it->first < imag_time_2) {
        const double t = delta_it->first;
        if (odd_parity && t_prev < t) {
            odd_integral += integrated_edge_energy(edg, t_prev, t);
        }

        bool toggles = false;
        do {
            toggles = !toggles;
            ++delta_it;
        } while (delta_it != delta_end && delta_it->first == t);
        if (toggles) {
            odd_parity = !odd_parity;
        }
        t_prev = t;
    }

    if (odd_parity && t_prev < imag_time_2) {
        odd_integral += integrated_edge_energy(edg, t_prev, imag_time_2);
    }

    return -2.0 * odd_integral;
}

/** @brief Integrate one spin without a Hamiltonian minus sign or coupling. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_edge_energy(
    const Edge& edg, double imag_time_1, double imag_time_2
) {
    if (imag_time_1 == imag_time_2) {
        throw std::invalid_argument(
            "integrated_edge_energy: time interval must be non-zero");
    }

    return detail::spin_integral({get_spin(edg), g[edg].spin_flips}, imag_time_1, imag_time_2);
}

/** @brief Integrate a tuple-product change from a local combination schedule. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_tuple_energy_diff_combination(
    std::span<const Edge> tuple_edges, 
    double imag_time_1, 
    double imag_time_2, 
    const std::vector<std::pair<double, int>>& spin_flip_lookup
) {
    return integrated_tuple_energy_diff_combination_from_flips(
        tuple_edges, imag_time_1, imag_time_2, spin_flip_lookup, false
    );
}

/** @brief Integrate a combination-update change using full or single-spin histories. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_tuple_energy_diff_combination_from_flips(
    std::span<const Edge> tuple_edges,
    double imag_time_1,
    double imag_time_2,
    const std::vector<std::pair<double, int>>& spin_flip_lookup,
    bool single_flips_only
) {
    if (spin_flip_lookup.empty()) {
        return 0.0;
    }

    // Tuple event toggles the local tuple-product iff the local overlap parity is odd.
    // In this code path, local overlap count == spin_flip_lookup.size() - 1.
    const bool tuple_flip_toggles = ((spin_flip_lookup.size() & 1) == 0);

    const auto* delta_it  = spin_flip_lookup.data();
    const auto* delta_end = delta_it + spin_flip_lookup.size();
    bool odd_parity = false;

    auto event_toggles = [&](const std::pair<double, int>& evt) {
        return (evt.second != 1) || tuple_flip_toggles;
    };

    while (delta_it != delta_end && delta_it->first < imag_time_1) {
        const double t = delta_it->first;
        bool toggles = false;
        do {
            if (event_toggles(*delta_it)) {
                toggles = !toggles;
            }
            ++delta_it;
        } while (delta_it != delta_end && delta_it->first == t);
        if (toggles) {
            odd_parity = !odd_parity;
        }
    }

    double odd_integral = 0.0;
    double t_prev = imag_time_1;

    bool spin_at_cutoff_ready = false;
    int spin_at_cutoff = 1;
    auto get_spin_at_cutoff = [&]() {
        if (!spin_at_cutoff_ready) {
            int spin_prod = 1;
            for (const Edge& e : tuple_edges) {
                const auto& flips = single_flips_only ? g[e].single_spin_flips : g[e].spin_flips;
                spin_prod *= detail::spin_at_time({get_spin(e), flips}, imag_time_2);
            }
            spin_at_cutoff = spin_prod;
            spin_at_cutoff_ready = true;
        }
        return spin_at_cutoff;
    };

    auto add_odd_interval = [&](double t_lo, double t_hi) {
        if (t_lo >= t_hi) {
            return;
        }
        if (t_lo < imag_time_2) {
            const double t_mid = std::min(t_hi, imag_time_2);
            if (t_lo < t_mid) {
                odd_integral += integrated_tuple_energy_from_flips(
                    tuple_edges, t_lo, t_mid, single_flips_only
                );
            }
            t_lo = t_mid;
        }
        if (t_lo < t_hi) {
            odd_integral += (t_hi - t_lo) * static_cast<double>(get_spin_at_cutoff());
        }
    };

    while (delta_it != delta_end) {
        const double t = delta_it->first;
        if (odd_parity && t_prev < t) {
            add_odd_interval(t_prev, t);
        }

        bool toggles = false;
        do {
            if (event_toggles(*delta_it)) {
                toggles = !toggles;
            }
            ++delta_it;
        } while (delta_it != delta_end && delta_it->first == t);
        if (toggles) {
            odd_parity = !odd_parity;
        }
        t_prev = t;
    }

    return -2.0 * odd_integral;
}

/** @brief Integrate the product of tuple spins without a minus sign or coupling. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_tuple_energy(
    std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2
) {
    return integrated_tuple_energy_from_flips(tuple_edges, imag_time_1, imag_time_2, false);
}

/** @brief Integrate a bare diagonal tuple product using only single-spin histories. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_tuple_energy_single_flips(
    std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2
) {
    return integrated_tuple_energy_from_flips(tuple_edges, imag_time_1, imag_time_2, true);
}

/** @brief Return -2 times the bare diagonal tuple integral for a product reversal. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_tuple_energy_diff_single_flips(
    std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2
) {
    return -2.0 * integrated_tuple_energy_single_flips(tuple_edges, imag_time_1, imag_time_2);
}

/** @brief Integrate a bare tuple product by merging the selected edge histories. */
[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_tuple_energy_from_flips(
    std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2, bool single_flips_only
) {
    return detail::spin_product_integral(
        tuple_edges, imag_time_1, imag_time_2,
        [&](const Edge& edg) -> detail::WorldlineView {
            return {get_spin(edg), single_flips_only ? g[edg].single_spin_flips : g[edg].spin_flips};
        }
    );
}

} // namespace paratoric
