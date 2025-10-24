// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#pragma once

#include "paratoric/types/types.hpp"
#include "rng/rng.hpp"

#include <boost/graph/adjacency_list.hpp>

#include <concepts>
#include <complex>
#include <filesystem>
#include <iostream>
#include <span>
#include <tuple>
#include <vector>

namespace paratoric {

/**
 * @class Lattice
 * @brief Data structure for storing and modifying a lattice for continuous QMC of any graph geometry.
 * 
 * TODO: This class is too large and too monolithic, split it up; externalize lattice construction
 */
class Lattice {    
public:
    // Class which is thrown in the search for percolation, see is_percolating()
    class FoundPercolation : public std::exception {
        public:
            const char* what() const noexcept override {
                return "Found percolation";
            }
    };

    /**
     * @brief Vertex properties for every vertex (i.e. lattice site) in the lattice.
     * 
     */
    struct VertexData {
        // In this vector the imaginary times of star flips are stored
        std::vector<double> star_flips;
        // In this vector we store the integrated potential energy of the star
        double integrated_star_energy;
        // These are the real space coordinates of the lattice site
        double x = 0.;
        double y = 0.;
        double z = 0.;
    };

    /**
     * @brief Edge properties for every edge (i.e. link between lattice sites) in the lattice.
     * 
     */
    struct EdgeData {
        // This int (+1 or -1) stores the local spin at imaginary time zero (beta)
        int spin = 1;
        // This string stores the spin for different snapshots
        std::string spin_string = "";
        // In this vector the imaginary times of any type of spin flip (single or tuple) are stored
        std::vector<double> spin_flips; 
        // In this vector the imaginary times of single spin flips are stored
        std::vector<double> single_spin_flips; 
        // In this vector we store the integrated potential energy of the edge
        double integrated_edge_energy;
        // In this vector we store all the plaquettes of which this edge is part of 
        std::vector<int> part_of_plaquette_lookup;
        // For the cubic lattice, this is x, y or z. Important for cube percolation
        std::string orientation;

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
        // Set LATTICE_DIMENSIONALITY later automatically
        LATTICE_TYPE = lat_spec.lattice_type;
        SYSTEM_SIZE = lat_spec.system_size;
        BETA = lat_spec.beta;
        BOUNDARIES = lat_spec.boundaries;
        DEFAULT_SPIN = lat_spec.default_spin;

        check_input_validity();

        g = init_lattice_graph(
            BASIS,
            LATTICE_TYPE, 
            SYSTEM_SIZE,
            BETA,
            BOUNDARIES,
            DEFAULT_SPIN
            );

        build_caches_();
        init_potential_energy();

        // Useful for debugging
        //write_graph("graph", "./");
    }  

    Lattice() = default; 
    ~Lattice() = default;

    Lattice(Lattice const&) = default;
    Lattice& operator=(Lattice const&) = default;

    /**
     * @brief Returns the total number of non-strings, i.e. links with spin +1.
     * 
     * @return total number of links with spin +1
     * 
     */
    inline int get_non_string_count();

    /**
     * @brief Returns the total number of strings, i.e. links with spin -1.
     * 
     * @return total number of links with spin -1
     * 
     */
    inline int get_string_count();

    /**
     * @brief Returns the total number of vertices in the graph.
     * 
     * @return number of vertices
     * 
     */
    inline int get_vertex_count();

    /**
     * @brief Returns the total number of edges in the graph.
     * 
     * @return number of edges
     * 
     */
    inline int get_edge_count();

    /**
     * @brief Returns the total number of plaquettes in the graph.
     * 
     * @return number of plaquettes
     * 
     */
    inline int get_plaquette_count();

    /**
     * @brief Returns the total number of cubes in the graph.
     * 
     * @return number of cubes
     * 
     */
    inline int get_cube_count();

    /**
     * @brief Returns the anyon count (e-anyons in the x-basis, m-anyons in the z-basis).
     * 
     * @return total anyon count
     * 
     */
    int get_anyon_count();

    /**
     * @brief Set the internal (pseudo-)random number generator.
     * 
     * @param rng_inp the random number generator of type RNG
     * 
     */
    void set_rng(std::shared_ptr<RNG>& rng_inp);

    /**
     * @brief Get the internal (pseudo-)random number generator.
     * 
     * @return random number generator stored in the Lattice instance
     * 
     */
    std::shared_ptr<RNG> get_rng();

    inline int get_spin(const Edge& edg);
    inline double get_potential_edge_energy(const Edge& edg);
    inline void set_potential_edge_energy(const Edge& edg, double potential_energy);
    inline void add_potential_edge_energy(const Edge& edg, double diff);
    inline double get_potential_star_energy(int star_index);
    inline void set_potential_star_energy(int star_index, double potential_energy);
    inline void add_potential_star_energy(int star_index, double diff);
    inline double get_potential_plaquette_energy(int plaquette_index);
    inline void set_potential_plaquette_energy(int plaquette_index, double potential_energy);
    inline void add_potential_plaquette_energy(int plaquette_index, double diff);

    /**
     * @brief Returns the orientation at the Edge edg.
     * 
     * @param edg the edge whose orientation is returned
     * 
     * @return orientation at edg
     * 
     */
    inline std::string get_orientation(const Edge& edg);

    /**
     * @brief Returns the imaginary time of the spin flip with index spin_flip_index at Edge edg.
     * 
     * @param edg the edge where the spin flip is located
     * @param spin_flip_index the index of the spin flip 
     * 
     * @return imaginary time of spin flip
     * 
     */
    inline double get_spin_flip_imag_time(const Edge& edg, int spin_flip_index);

    /**
     * @brief Returns the index of the spin flip with imaginary time tau at Edge edg.
     * 
     * @details Throw std::runtime_error when imaginary time is not found.
     * 
     * @param edg the edge where the spin flip is located
     * @param tau the imaginary time of the spin flip 
     * 
     * @return index of spin flip
     * 
     */
    int get_spin_flip_index(const Edge& edg, double tau);

    /**
     * @brief Returns the index of the single spin flip with imaginary time tau at Edge edg.
     * 
     * @details Throw std::runtime_error when imaginary time is not found.
     * 
     * @param edg the edge where the spin flip is located
     * @param tau the imaginary time of the spin flip 
     * 
     * @return index of spin flip
     * 
     */
    int get_single_spin_flip_index(const Edge& edg, double tau);

    /**
     * @brief Returns the single spin flip vector at Edge edg.
     * 
     * @param edg the edge on which the single spin flips should be returned
     * 
     * @return a vector containing the single spin flips at Edge edg
     * 
     */
    std::span<const double> get_single_spin_flips(const Edge& edg);

    /**
     * @brief Returns the tuple spin flip vector at tuple with tuple_index t_index.
     * 
     * @param t_index the index of the tuple
     * @return a vector containing the tuple spin flips at tuple t_index
     * 
     */
    std::span<const double> get_tuple_spin_flips(int t_index);

    /**
     * @brief Change the imaginary time of the spin flip with index spin_flip_index at Edge edg to imag_time.
     * 
     * @param edg the edge where the spin flip is located
     * @param spin_flip_index the index of the spin flip 
     * @param imag_time new imaginary time
     * 
     */
    inline void set_spin_flip_imag_time(
        const Edge& edg, int spin_flip_index, double imag_time
    );

    /**
     * @brief Change the imaginary time of the single spin flip with index spin_flip_index at Edge edg to imag_time.
     * 
     * @param edg the edge where the spin flip is located
     * @param spin_flip_index the index of the spin flip 
     * @param imag_time new imaginary time
     * 
     */
    inline void set_single_spin_flip_imag_time(
        const Edge& edg, int spin_flip_index, double imag_time
    );

    /**
     * @brief Deletes the imaginary time of the single spin flips with imaginary times imag_time_single_spin_flip and imag_time_next_single_spin_flip at Edge edg.
     * 
     * @param edg the edge where the single spin flips are deleted
     * @param imag_time_single_spin_flip the imaginary time of the single spin flip 
     * @param imag_time_next_single_spin_flip the imaginary time of the next single spin flip
     * 
     */
    void delete_double_single_spin_flip(
        const Edge& edg, double imag_time_single_spin_flip, double imag_time_next_single_spin_flip
    );

    /**
     * @brief Deletes the imaginary time of the single spin flip with index spin_flip_index_1 at Edge edg.
     * 
     * @param edg the edge where the single spin flip is deleted
     * @param spin_flip_index the index of the single spin flip 
     * 
     */
    void delete_single_spin_flip(const Edge& edg, int spin_flip_index);

    /**
     * @brief Deletes the imaginary times of the tuple flips with imaginary time imag_time_tuple_flip and imag_time_next_tuple_flip at the tuple with edges tuple_edges.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param imag_time_tuple_flip the imaginary time of the tuple flip 
     * @param imag_time_next_tuple_flip the imaginary time of the next tuple flip 
     * 
     */
    void delete_double_tuple_flip(
        int tuple_index, 
        std::span<const Edge> tuple_edges, 
        double imag_time_tuple_flip, 
        double imag_time_next_tuple_flip
    );

    /**
     * @brief Deletes the imaginary time of the tuple flip with imaginary time imag_time_tuple_flip at the tuple with edges tuple_edges.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param imag_time_tuple_flip the imaginary time of the tuple flip 
     * 
     */
    void delete_tuple_flip(
        int tuple_index, std::span<const Edge> tuple_edges, double imag_time_tuple_flip
    );

    /**
     * @brief Calculates whether there is a tuple flip at imaginary time tau on all edges of the tuple with edges tuple_edges.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param tau the imaginary time of the potential tuple flip 
     * 
     * @return true if there is a tuple flip at tau else false
     * 
     */
    bool check_tuple_flip_present_tuple(std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Calculates whether there is a plaquette flip at imaginary time tau at the plaquettes that contain the Edge edg.
     * 
     * @param edg the edge which is part of plaquettes
     * @param tau the imaginary time of the potential plaquette flip(s)
     * 
     * @return true if there is a plaquette flip at tau else false
     * 
     */
    bool check_plaquette_flip_at_edge(const Edge& edg, double tau);

    /**
     * @brief Calculates whether there is a star flip at imaginary time tau at the stars that contain the Edge edg.
     * 
     * @param edg the edge which is part of stars
     * @param tau the imaginary time of the potential star flip(s)
     * 
     * @return true if there is a star flip at tau else false
     * 
     */
    bool check_star_flip_at_edge(const Edge& edg, double tau);

    /**
     * @brief Calculates whether there is a spin flip at imaginary times tau_1 or tau_2 at the tuple with edges tuple_edges.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param tau the imaginary time of the potential spin flip
     * 
     * @return true if there is a tuple flip at tau_1 or tau_2 else false
     * 
     */
    bool check_spin_flips_present_tuple(std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Calculates the imaginary time of the next spin flip (of any type) after tau at the edge edg. 
     * If there is a spin_flip at tau, the next spin flip (possibly going over beta) is returned. 
     * If there is no spin flip, tau is returned. 
     * If there is at least one spin flip but not at tau, the next spin flip after tau (possibly going over beta) is returned.
     * 
     * @param edg the edge of interest
     * @param tau the imaginary time after which we look for a spin flip
     * 
     * @return imaginary time of next spin flip (of any type)
     * 
     */
    double flip_next_imag_time(const Edge& edg, double tau);

    /**
     * @brief Calculates the imaginary times of the next spin flips (of any type) after tau at the tuple with edges tuple_edges. 
     * If there is a spin_flip at tau, the next spin flip (possibly going over beta) is returned. 
     * If there is no spin flip, tau is returned. 
     * If there is at least one spin flip but not at tau, the next spin flip after tau (possibly going over beta) is returned.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param tau the imaginary time after which we look for spin flips
     * 
     * @return vector of imaginary times of next spin flips (of any type), see code for order of edges
     * 
     */
    std::vector<double> flip_next_imag_times_tuple(std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Calculates the imaginary time of the previous spin flip (of any type) before tau at the edge edg. 
     * If there is a spin_flip at tau, the previous spin flip (possibly going over beta) is returned. 
     * If there is no spin flip, tau is returned. 
     * If there is at least one spin flip but not at tau, the previous spin flip before tau (possibly going over beta) is returned.
     * 
     * @param edg the edge of interest
     * @param tau the imaginary time before which we look for a spin flip
     * 
     * @return imaginary time of previous spin flip (of any type)
     * 
     */
    double flip_prev_imag_time(const Edge& edg, double tau);

    /**
     * @brief Calculates the imaginary times of the previous spin flips (of any type) before tau at the tuple with edges tuple_edges. 
     * If there is a spin_flip at tau, the previous spin flip (possibly going over beta) is returned. 
     * If there is no spin flip, tau is returned. 
     * If there is at least one spin flip but not at tau, the previous spin flip before tau (possibly going over beta) is returned.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param tau the imaginary time before which we look for spin flips
     * 
     * @return vector of imaginary times of previous spin flips (of any type), see code for order of edges
     * 
     */
    std::vector<double> flip_prev_imag_times_tuple(std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Inserts the imaginary time of the spin flips with imaginary times tau_left and tau_right at Edge edg.
     * 
     * @param edg the edge where the spin flips are inserted
     * @param tau_left the imaginary time of the first single spin flip 
     * @param tau_right the imaginary time of the second single spin flip
     * 
     */
    void insert_double_single_spin_flip(const Edge& edg, double tau_left, double tau_right);

    /**
     * @brief Inserts the imaginary time of the spin flip with imaginary time tau at Edge edg.
     * 
     * @param edg the edge where the single spin flip is inserted
     * @param tau the imaginary time of the single spin flip 
     * 
     */
    void insert_single_spin_flip(const Edge& edg, double tau);

    /**
     * @brief Inserts the imaginary times of the tuple flips with imaginary time tau_left and tau_right at the tuple with edges tuple_edges.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param tau_left the imaginary time of the first tuple flip 
     * @param tau_right the imaginary time of the second tuple flip 
     * 
     */
    void insert_double_tuple_flip(
        int tuple_index, std::span<const Edge> tuple_edges, double tau_left, double tau_right
    );

    /**
     * @brief Inserts the imaginary time of the tuple flip with imaginary time tau at the tuple with edges tuple_edges.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param tau the imaginary time of the tuple flip 
     * 
     */
    void insert_tuple_flip(int tuple_index, std::span<const Edge> tuple_edges, double tau);

    /**
     * @brief Move the imaginary time of the spin flip with index spin_flip_index to imaginary time tau_new at Edge edg. 
     * no_move_over_beta specifies if the move is going over beta (false) or not (true).
     * 
     * @param edg the edge where the single spin flip is inserted
     * @param spin_flip_index the (old) spin flip index of the spin before the move
     * @param tau_new the new imaginary time
     * @param no_move_over_beta if the move is going over beta (false) or not (true) 
     * 
     */
    void move_spin_flip(const Edge& edg, int spin_flip_index, double tau_new, bool no_move_over_beta);

    /**
     * @brief Move the imaginary time of the tuple flip with old imaginary time tau_old to imaginary time tau_new at the tuple with edges tuple_edges.
     * no_move_over_beta specifies if the move is going over beta (false) or not (true).
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * @param tau_old the (old) imaginary time of the tuple flip before the move
     * @param tau_new the (new) imaginary time of the tuple flip after the move
     * @param no_move_over_beta if the move is going over beta (false) or not (true) 
     * 
     */
    void move_tuple_flip(
        int tuple_index, 
        std::span<const Edge> tuple_edges, 
        double tau_old, 
        double tau_new, 
        bool no_move_over_beta
    );

    /**
     * @brief Print all spins at imaginary time 0.
     * 
     */
    void print_spins();

    /**
     * @brief Print all imaginary times of spin flips at Edge edg.
     * 
     * @param edg the edge whose spin flips are printed 
     * 
     */
    void print_spin_flip_imag_times(const Edge& edg);

    /**
     * @brief Print all imaginary times of spin flips at the tuple with edges tuple_edges.
     * 
     * @param tuple_edges the vector which stores the edges of the tuple
     * 
     */
    void print_tuple_flip_imag_times(std::span<const Edge> tuple_edges);

    /**
     * @brief Return the number of spin flips at Edge edg.
     * 
     * @param edg the edge whose number of spin flips is returned 
     * 
     * @return number of spin flips on edg
     */
    inline int get_spin_flip_count(const Edge& edg);

    /**
     * @brief Flips the spin on edge edg.
     * 
     * @param edg the edge where the spin is flipped
     * 
     */
    void flip_spin(const Edge& edg);

    /**
     * @brief Flips the spins on the star with center v.
     * 
     * @param v the star center vertex 
     * 
     */
    void flip_star(int v);
    inline Edge edge_in_between(int v_1, int v_2);
    inline bool exists_edge(int v_1, int v_2);
    inline std::pair<int, int> vertices_of_edge(const Edge& edg);

    /**
     * @brief Return randomly selected edge from graph.
     * 
     * @return tuple of the edge descriptor, the source vertex and the target vertex of the selected edge
     * 
     */
    std::tuple<Edge, int, int> get_random_edge();

    /**
     * @brief Return randomly selected vertex from graph.
     * 
     * @return the index of the random vertex
     * 
     */
    int get_random_vertex();

    /**
     * @brief Returns the plaquette index of a random plaquette.
     * 
     * @return random plaquette index
     * 
     */
    int get_random_plaquette_index();
    std::vector<std::pair<int,int>> get_plaquette_vertex_pairs(int p_index);
    std::vector<int> get_cube_vertices(int c_index);
    inline std::span<const Edge> get_plaquette_edges(int p_index);
    inline std::span<const Edge> get_star_edges(int center_index);
    int get_tuple_sum(std::span<const Edge> tuple_edges);
    int get_tuple_prod(std::span<const Edge> tuple_edges);
    void flip_tuple(std::span<const Edge> tuple_edges);
    int get_vertex_nn_spins_prod(int v);

    /**
     * @brief Return the tuple energy DIFFERENCE (integrated over imaginary time) of the tuple with edges tuple_edges when flipping the spin 
     * on one edge in the tuple between the imaginary times imag_time_1 and imag_time_2.
     * 
     * @param tuple_edges the edges of the tuple
     * @param imag_time_1 the lower bound for imaginary time
     * @param imag_time_2 the upper bound for imaginary time
     * 
     * @return the tuple energy difference
     * 
     */
    double integrated_tuple_energy_diff(
        std::span<const Edge> tuple_edges, 
        double imag_time_1, 
        double imag_time_2
    );

    /**
     * @brief Return the tuple energy DIFFERENCE (integrated over imaginary time) of the tuple with edges tuple_edges when flipping the spins in spin_flip_lookup
     * in the relevant time interval from imag_time_1 to imag_time_2.
     * 
     * @param tuple_edges the edges of the tuple
     * @param imag_time_1 the lower bound for imaginary time
     * @param imag_time_2 the upper bound for imaginary time
     * @param spin_flip_lookup contains pairs with the imaginary time of the spin flip and the spin flip type (1: tuple, 0: single)
     * 
     * @return the tuple energy difference
     * 
     */
    inline double integrated_tuple_energy_diff_combination(
        std::span<const Edge> tuple_edges, 
        double imag_time_1, 
        double imag_time_2, 
        const std::vector<std::pair<double, int>>& spin_flip_lookup
    );

    /**
     * @brief Return the tuple energy (integrated over imaginary time) of the tuple with edges tuple_edges between the imaginary times imag_time_1 and imag_time_2.
     * 
     * @param tuple_edges the edges of the tuple
     * @param imag_time_1 the lower bound for imaginary time
     * @param imag_time_2 the upper bound for imaginary time
     * 
     * @return the tuple energy
     * 
     */
    inline double integrated_tuple_energy(
        std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2
    );

    /**
     * @brief Integrated STAR–energy difference for the two stars incident to an edge.
     *
     * Computes the (bare) change of each star’s energy, integrated over the
     * imaginary-time interval [@p imag_time_1, @p imag_time_2], when the spin on
     * @p edg is flipped across that interval. The two affected stars are the ones
     * centered at the edge’s source and target vertices. Returns their individual
     * contributions and their sum. Couplings are NOT applied here.
     *
     * @param edg            Edge whose spin is flipped; its endpoints define the two stars.
     * @param imag_time_1    Lower bound of imaginary time (inclusive).
     * @param imag_time_2    Upper bound of imaginary time (inclusive).
     * @param total_cache    If true, use the cached integrated star energies and
     *                       return @c -2 * cache for each of the two stars. This is
     *                       a constant-time fast path intended for full-period flips
     *                       (e.g., [0, β]). Do not set true for partial intervals.
     *
     * @return std::tuple<
     *           double,                 // sum of bare star-energy differences for the two stars
     *           std::vector<int>,       // star centers (vertex indices), order: [source_v, target_v]
     *           std::vector<double>     // per-star bare energy differences on [t1,t2],
     *                                   // aligned with the centers vector
     *         >
     *
     * @pre @p imag_time_2 > @p imag_time_1. The interval must be non-zero.
     * @throws std::invalid_argument if @p imag_time_1 == @p imag_time_2.
     *
     * @note “Bare” means no coupling strength factors are applied here.
     *       Apply couplings externally when forming acceptance ratios or totals.
     */
    std::tuple<double, std::vector<int>, std::vector<double>> 
    integrated_star_energy_diff(
        const Edge& edg, double imag_time_1, double imag_time_2, bool total_cache
    );

    /**
     * @brief Integrated plaquette-energy difference around an edge on a time interval.
     *
     * Computes, for every plaquette that touches @p edg, the bare energy change
     * when flipping the spin on @p edg over the imaginary-time interval 
     * [@p imag_time_1, @p imag_time_2]. Returns both the per-plaquette
     * contributions and their sum. Couplings strengths are NOT applied here.
     *
     * @param edg            Edge whose spin is flipped.
     * @param imag_time_1    Lower bound of imaginary time (inclusive).
     * @param imag_time_2    Upper bound of imaginary time (inclusive).
     * @param total_cache    If true, use the cached integrated plaquette energy and
     *                       return @c -2 * cache for each adjacent plaquette. This
     *                       is a constant-time fast path intended for full-period
     *                       flips (i.e. [0, β]). Do not set true for partial
     *                       intervals.
     *
     * @return std::tuple<
     *           double,                 // sum over adjacent plaquettes of their bare energy differences
     *           std::vector<int>,       // plaquette indices touching @p edg (same order as adjacency)
     *           std::vector<double>     // per-plaquette bare energy differences on [t1,t2]
     *         >
     *
     * @pre imag_time_2 > imag_time_1.  The interval must be non-zero.
     * @throws std::invalid_argument if @p imag_time_1 == @p imag_time_2.
     *
     * @note “Bare” means no coupling factors are applied. Apply coupling strength factors externally.
     * @note The two output vectors have the same length: the plaquette coordination of @p edg
     *       (typically 2 in 2D). Their order matches @c g[edg].part_of_plaquette_lookup.
     */
    std::tuple<double, std::vector<int>, std::vector<double>> 
    integrated_plaquette_energy_diff(
        const Edge& edg, double imag_time_1, double imag_time_2, bool total_cache
    );

    /**
     * @brief Return the total star energy (integrated over imaginary time) of all stars 
     * between the imaginary times 0 and beta.
     * 
     * @return the total integrated star energy
     * 
     */
    double total_integrated_star_energy();

    /**
     * @brief Return the total plaquette energy (integrated over imaginary time) of all plaquettes 
     * between the imaginary times 0 and beta.
     * 
     * @return the total integrated plaquette energy
     * 
     */
    double total_integrated_plaquette_energy();

    /**
     * @brief Integrated STAR–energy difference for a tuple–plus–single–flip combination update.
     *
     * For each star centered at a vertex touched by any edge in @p plaquette_edges,
     * this computes the (bare) change of the star energy integrated over the
     * imaginary-time interval [@p imag_time_1, @p imag_time_2], when a tuple flip
     * occurs at @p imag_time_tuple_flip and selected single-edge flips occur on the
     * edges of @p plaquette_edges at the times given in @p spin_flip_lookup.
     * The result returns the per-star contributions and their sum. Couplings are
     * NOT applied here.
     *
     * @param plaquette_index    The index of the plaquette under consideration.
     *                           The edges of this plaquette determine which stars (their incident
     *                           vertices) are affected.
     * @param imag_time_1        Lower bound of imaginary time (inclusive).
     * @param imag_time_2        Upper bound of imaginary time (inclusive).
     * @param spin_flip_lookup   Per-edge flip times aligned with @p plaquette_edges:
     *                           @c spin_flip_lookup[k] is the imaginary time of the
     *                           single spin flip on @c plaquette_edges[k]. Only those
     *                           edges that belong to a given star contribute to that
     *                           star’s local flip schedule.
     * @param imag_time_tuple_flip  Imaginary time of the tuple (plaquette) flip event.
     *
     * @return std::tuple<
     *           double,                 // sum of bare star-energy differences over all affected stars
     *           std::vector<int>,       // unique star centers (vertex indices), sorted ascending
     *           std::vector<double>     // per-star bare energy differences on [t1,t2],
     *                                   // in the same order as the vertex index vector
     *         >
     *
     * @pre @p imag_time_2 > @p imag_time_1. The interval must be non-zero.
     * @throws std::invalid_argument if @p imag_time_1 == @p imag_time_2.
     *
     * @note “Bare” means no coupling strength factors are applied.
     * @note For each star, the local flip schedule consists of the tuple flip
     *       at @p imag_time_tuple_flip and any single flips from @p plaquette_edges
     *       that are incident on the star. Equal-time flips are handled by parity,
     *       so their order does not affect the integral.
     * @note The vertex list is deduplicated and sorted; the per-star values are
     *       aligned with that list one-to-one.
     */
    std::tuple<double, std::vector<int>, std::vector<double>> 
    integrated_star_energy_diff_combination(
        int plaquette_index, 
        double imag_time_1, 
        double imag_time_2, 
        std::span<const double> spin_flip_lookup, 
        double imag_time_tuple_flip
    );

    /**
     * @brief Integrated PLAQUETTE–energy difference for a star–plus–single–flip combination update.
     *
     * For each plaquette touching any edge in @p star_edges, compute the (bare)
     * change of the plaquette energy integrated over
     * [@p imag_time_1, @p imag_time_2], when a tuple flip occurs at
     * @p imag_time_tuple_flip and selected single-edge flips occur on
     * @p star_edges at the times in @p spin_flip_lookup.
     * Returns the per-plaquette contributions and their sum.
     * Couplings are NOT applied here.
     *
     * @param star_index           Index of the star under consideration.
     *                             The edge of this star determine which plaquettes are affected.
     * @param imag_time_1          Lower bound of imaginary time (inclusive).
     * @param imag_time_2          Upper bound of imaginary time (inclusive).
     * @param spin_flip_lookup     Per-edge flip times aligned with @p star_edges:
     *                             @c spin_flip_lookup[k] is the imaginary time of the
     *                             single spin flip on @c star_edges[k]. Only those
     *                             edges that belong to a given plaquette contribute to
     *                             that plaquette’s local flip schedule.
     * @param imag_time_tuple_flip Imaginary time of the tuple (star) flip event.
     *
     * @return std::tuple<
     *           double,                 // sum of bare plaquette-energy differences over all affected plaquettes
     *           std::vector<int>,       // unique plaquette indices, sorted ascending
     *           std::vector<double>     // per-plaquette bare energy differences on [t1,t2],
     *                                   // in the same order as the index vector
     *         >
     *
     * @pre @p imag_time_2 > @p imag_time_1. The interval must be non-zero.
     * @throws std::invalid_argument if @p imag_time_1 == @p imag_time_2.
     *
     * @note “Bare” means no coupling factors (e.g., @c -J) are applied.
     * @note For each plaquette, the local flip schedule consists of the tuple flip
     *       at @p imag_time_tuple_flip and any single flips from @p star_edges
     *       that are incident on that plaquette. Equal-time flips are combined by
     *       parity; order does not affect the integral.
     * @note The plaquette list is deduplicated and sorted; per-plaquette values
     *       align one-to-one with that list.
     */
    std::tuple<double, std::vector<int>, std::vector<double>> 
    integrated_plaquette_energy_diff_combination(
        int star_index, 
        double imag_time_1, 
        double imag_time_2, 
        std::span<const double> spin_flip_lookup, 
        double imag_time_tuple_flip
    );

    /**
     * @brief Return the edge energy DIFFERENCE (integrated over imaginary time) at the edge edg when flipping the spin 
     * between the imaginary times imag_time_1 and imag_time_2.
     * 
     * @param edg the edge of interest
     * @param imag_time_1 the lower bound for imaginary time
     * @param imag_time_2 the upper bound for imaginary time
     * 
     * @return the edge energy difference
     * 
     */
    double integrated_edge_energy_diff(const Edge& edg, double imag_time_1, double imag_time_2);

    /**
     * @brief Return the edge energy DIFFERENCE (integrated over imaginary time) at the edge edg when flipping the spins in spin_flip_lookup in the relevant time interval
     * from imag_time_1 to imag_time_2.
     * 
     * @param edg the edge of interest
     * @param imag_time_1 the lower bound for imaginary time
     * @param imag_time_2 the upper bound for imaginary time
     * @param spin_flip_lookup contains pairs with the imaginary time of the spin flip and the spin flip type (1: tuple, 0: single)
     * 
     * @return the edge energy difference
     * 
     */
    inline double integrated_edge_energy_diff_combination(
        const Edge& edg, 
        double imag_time_1, 
        double imag_time_2, 
        std::vector<std::pair<double,int>>& spin_flip_lookup
    );

    /**
     * @brief Return the edge energy (integrated over imaginary time) at the edge edg 
     * between the imaginary times imag_time_1 and imag_time_2.
     * 
     * @param edg the edge of interest
     * @param imag_time_1 the lower bound for imaginary time
     * @param imag_time_2 the upper bound for imaginary time
     * 
     * @return the edge energy
     * 
     */
    inline double integrated_edge_energy(const Edge& edg, double imag_time_1, double imag_time_2);

    //TODO
    double integrated_edge_energy_weighted(const Edge& edg, double imag_time_1, double imag_time_2);

    /**
     * @brief Return the total edge energy (integrated over imaginary time) over all edges 
     * between the imaginary times 0 and beta.
     * 
     * @return the total integrated edge energy
     * 
     */
    double total_integrated_edge_energy();

    //TODO
    double total_integrated_edge_energy_weighted();

    //TODO
    void init_potential_energy();

    /**
     * @brief Returns the sum of the spin on all links at imaginary time beta (or equivalently zero) WITHOUT an overall negative sign WITHOUT multiplying it by h/lmbda.
     * 
     * @return spin energy 
     * 
     */
    double get_diag_single_energy();

    /**
     * @brief Returns the total magnetization in a format ready for the susceptibility calculation.
     * 
     * @return magnetization 
     * 
     */
    std::complex<double> get_diag_M_M();

    //TODO
    std::complex<double> get_diag_dynamical_M_M();

    /**
     * @brief Returns the total magnetization in a format ready for the susceptibility calculation.
     * 
     * @return magnetization 
     * 
     */
    std::complex<double> get_non_diag_M_M();

    /**
     * @brief Returns kL and kR in Eq. (9) of https://doi.org/10.1103/PhysRevX.5.031007.
     * 
     * @return kL and kR
     * 
     * @author Simon Mathias Linsel
     */
    std::complex<double> get_kL_kR_single();

    /**
     * @brief Returns the sum of all single spin flips in the imaginary time axis of the lattice divided by beta WITHOUT an overall negative sign. The result is the gauge field energy term multiplied by lmbda!
     * 
     * @return gauge field energy multiplied by lmbda 
     * 
     */
    double get_non_diag_single_energy_x();

    /**
     * @brief Returns the sum of all single spin flips in the imaginary time axis of the lattice divided by beta WITHOUT an overall negative sign. The result is the electric field energy term multiplied by h!
     * 
     * @return electric field energy multiplied by h 
     * 
     */
    double get_non_diag_single_energy_z();

    /**
     * @brief Returns the sum of all plaquette flips in the imaginary time axis of the lattice divided by beta WITHOUT an overall negative sign. The result is the plaquette energy term multiplied by J!
     * 
     * @return plaquette energy multiplied by J 
     * 
     */
    double get_non_diag_tuple_energy_x();

    /**
     * @brief Returns the sum of all star flips in the imaginary time axis of the lattice divided by beta WITHOUT an overall negative sign. The result is the star energy term multiplied by mu!
     * 
     * @return star energy multiplied by mu 
     * 
     */
    double get_non_diag_tuple_energy_z();

    /**
     * @brief Returns the sum of all star terms at imaginary time beta (or equivalently zero) WITHOUT an overall negative sign WITHOUT multiplying it by mu.
     * 
     * @return star energy with a minus sign 
     * 
     */
    double get_diag_tuple_energy_x();

    /**
     * @brief Returns the sum of all plaquette terms at imaginary time beta (or equivalently zero) WITHOUT an overall negative sign WITHOUT multiplying it by J.
     * 
     * @return plaquette energy with a minus sign 
     * 
     */
    double get_diag_tuple_energy_z();

    /**
     * @brief Returns the Fredenhagen-Marcu order parameter at equal imaginary time, see https://doi.org/10.1103/PhysRevLett.56.223.
     * 
     * @return complex number where the half Wilson/'t Hooft loop is the real part and the full Wilson/'t Hooft loop is the imaginary part
     * 
     */
    std::complex<double> fredenhagen_marcu();

    /**
     * @brief Returns the staggered imaginary plaquette flip imaginary time differences order parameter similar to https://doi.org/10.1103/PhysRevB.85.195104.
     * 
     * @return staggered imaginary times order parameter 
     * 
     */
    double get_staggered_imaginary_times_plaquette();

    /**
     * @brief Returns the staggered imaginary star flip imaginary time differences order parameter similar to https://doi.org/10.1103/PhysRevB.85.195104.
     * 
     * @return staggered imaginary times order parameter 
     * 
     */
    double get_staggered_imaginary_times_star();

    /**
     * @brief Searches for a winding percolating ("global") cluster of strings (i.e. spin = -1) in the lattice. Depends on the definition of coordinates in the lattice. 
     * 
     * @return bool which signals winding-percolation (1) or non-winding-percolation (0)
     * 
     */
    bool is_winding_percolating();

    /**
     * @brief Searches for a percolating ("global") cluster of strings (i.e. spin = -1) in the lattice. Depends on the definition of coordinates in the lattice. 
     * 
     * @return bool which signals percolation (1) or non-percolation (0)
     * 
     */
    bool is_percolating();

    /**
     * @brief Searches for a winding plaquette percolating ("global") cluster of plaquettes (i.e. connected by  spin = -1) in the lattice. Depends on the definition of coordinates in the lattice. 
     * 
     * @return bool which signals plaquette winding-percolation (1) or plaquette non-winding-percolation (0)
     * 
     */
    bool is_winding_plaquette_percolating();

    /**
     * @brief TODO.
     * 
     * @return TODO.
     * 
     */
    bool is_winding_cube_percolating();

    /**
     * @brief Returns the number of strings in the largest cluster of strings (i.e. spin = -1) in the lattice.
     * 
     * @return number of edges in the largest string-cluster.
     * 
     */
    int largest_cluster();

    /**
     * @brief Returns the number of plaquettes in the largest plaquette cluster connected strings (i.e. spin = -1) in the lattice.
     * 
     * @return number of plaquettes in the largest plaquette-cluster.
     * 
     */
    int largest_plaquette_cluster();

    /**
     * @brief Returns the percolation strength, i.e. the number of strings in the largest string-cluster divided by the total number of strings. If no cluster percolates it returns 0.
     * 
     * @return largest string cluster divided by number of strings if percolating, else 0.
     * 
     */
    double percolation_strength();

    /**
     * @brief Returns the percolation probability, i.e. whether there is at least one percolating cluster. The term probability does only make sense when averaging over many samples.
     * 
     * @return 1 if we have percolating cluster, else 0.
     * 
     */
    double percolation_probability();

    /**
     * @brief TODO.
     * 
     * @return TODO.
     * 
     */
    double plaquette_percolation_strength();

    /**
     * @brief Returns the plaquette percolation probability, i.e. whether there is at least one percolating plaquette cluster. The term probability does only make sense when averaging over many samples.
     * 
     * @return 1 if we have percolating plaquette cluster, else 0.
     * 
     */
    double plaquette_percolation_probability();

    /**
     * @brief TODO.
     * 
     * @return TODO.
     * 
     */
    double cube_percolation_strength();

    /**
     * @brief Returns the cube percolation probability, i.e. whether there is at least one percolating cube cluster. The term probability does only make sense when averaging over many samples.
     * 
     * @return 1 if we have percolating cube cluster, else 0.
     * 
     */
    double cube_percolation_probability();

    /**
     * @brief Perform a global rotation of the imaginary time by a random tau_0 with 0 <= tau_0 <= beta.
     * 
     */
    void rotate_imag_time();

    /**
     * @brief Appends the current spin to the spin_string for each edge.
     * 
     */
    void update_spin_string();

    /**
     * @brief Writes out a graphml XML file which contains the graph vertices (with occupation number and coordinates) and edges (with spins) 
     * 
     * @param file_name the filename of the XML file
     * @param output_directory the directory where the method will write the graphml file to
     * 
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

    // The object where all physical information is stored in
    LatticeGraph g;
    // This vector stores all elementary plaquettes in the lattice
    std::vector<std::vector<std::pair<int,int>>> plaquette_vector;
    // This vector stores the imaginary times of all plaquette spin flips
    std::vector<std::vector<double>> plaquette_flip_vector;
    // This vector stores the integrated potential energies of the plaquette spin flips
    std::vector<double> integrated_plaquette_energy_vector;
    //These vectors store the coordinates of the plaquettes
    std::vector<double> plaquette_x_vector;
    std::vector<double> plaquette_y_vector;
    std::vector<double> plaquette_z_vector;
    // This vector stores all elementary cubes in the lattice (3D)
    std::vector<std::vector<int>> cube_vector;
    //These vectors store the coordinates of the cubes
    std::vector<double> cube_x_vector;
    std::vector<double> cube_y_vector;
    std::vector<double> cube_z_vector;
    // In this vector we store all the cubes of which this plaquette is part of 
    std::vector<std::vector<int>> plaquette_part_of_cube_lookup;
    // In this vector we store all the plaquettes which are part of a given plaquette 
    std::vector<std::vector<int>> cube_has_plaquettes_lookup;

    // p -> edges of plaquette p (arbitrary length: 3/4/6/…)
    std::vector<std::vector<Edge>> plaquette_edges_cache_;
    // v -> incident edges (star at vertex v)
    std::vector<std::vector<Edge>> star_edges_cache_;
    // Edge descriptors
    std::vector<Edge> egde_cache_;

    // These two vectors are used to calculate the Fredenhagen-Marcu order parameter
    std::vector<VertexPair> half_path_vector;
    std::vector<VertexPair> full_path_vector;

    std::vector<double> MAX_COORDINATES;
    std::vector<double> MAX_PLAQUETTE_COORDINATES;

    /**
     * @brief Check the input validity before trying to create the LatticeGraph.
     * 
     * @details Will throw std::invalid_argument when input does not make sense.
     * 
     */
    void check_input_validity() const;
    void build_caches_();

    /**
     * @brief Constructs and returns LatticeGraph object which stores all physical information
     * 
     * @param basis the spin eigenbasis 
     * @param lattice_type the lattice type, e.g. "triangular"
     * @param L the system size of the lattice (in one dimension)
     * @param beta inverse temperature
     * @param boundaries the boundary condition of the lattice (periodic or open)
     * @param default_spin the default spin on the links (1 or -1)
     * 
     * @return LatticeGraph boost graph adjacency list
     * 
     */
    LatticeGraph init_lattice_graph(
        char basis,
        const std::string& lattice_type, 
        int L,
        double beta,
        const std::string& boundaries,
        int default_spin
        );
    std::vector<std::pair<Vertex, Vertex>> edge_vector;
    mutable std::uniform_int_distribution<int> edge_dist;
    mutable std::uniform_int_distribution<int> vertex_dist;
    mutable std::uniform_int_distribution<int> plaquette_dist;

    /**
     * @brief Helper function to construct 1) a half Wilson/'t Hooft loop and 2) a full Wilson/'t Hooft loop to calculate the Fredenhagen-Marcu order parameter.
     * 
     * @param start_y the upper bound of the full loop
     * @param end_y the lower bound of the half/full loop
     * @param middle_y the upper bound for the half loop
     * @param start_x the lower bound for the half/full loop
     * @param end_x the upper bound for the half/full loop
     * @return Half loop and full loop vertex pairs
     * 
     */
    std::pair<std::vector<VertexPair>, std::vector<VertexPair>> 
    construct_fredenhagen_marcu_loops(
        int start_y, int end_y, int middle_y, int start_x, int end_x, char basis
    );
    // RNG can be copied (with automatic reseed)
    std::shared_ptr<RNG> rng;
    std::uniform_real_distribution<double> uniform_dist{0., 1.};

    /**
     * @brief Return randomly selected element from iterable.
     * 
     * @tparam Iter the iterable from which elements are drawn
     * @tparam RandomGenerator the random number generator
     * @param start the lower index for the random element (included)
     * @param end the upper index for the random element (not included)
     * @param g the random number generator 
     * @return Iter advanced to a random element
     */
    template<typename Iter, typename RandomGenerator>
    Iter random_element(Iter start, Iter end, RandomGenerator& gen) {
        std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
        std::advance(start, dis(gen));
        return start;
    }

    /**
     * @brief Return randomly selected element from iterable.
     *
     * @tparam Iter the iterable from which elements are drawn
     * @param start the lower index for the random element (included)
     * @param end the upper index for the random element (not included)
     * @param g the random number generator
     * @return Iter advanced to a random element
     */
    template<typename Iter>
    Iter random_element(Iter start, Iter end) {
        static std::random_device rd;
        static RNG gen(rd());
        return random_element(start, end, gen);
    }

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

inline int Lattice::get_non_string_count() {
    int result = 0;
    for (const auto& edg : egde_cache_) {
        if (g[edg].spin == 1) result += 1;
    }
    return result;
}

inline int Lattice::get_string_count() {
    int result = 0;
    for (const auto& edg : egde_cache_) {
        if (g[edg].spin == -1) result += 1;
    }
    return result;
}

inline int Lattice::get_vertex_count() {
    return boost::num_vertices(g);
}

inline int Lattice::get_edge_count() {
    return boost::num_edges(g);
}

inline int Lattice::get_plaquette_count() {
    return plaquette_vector.size();
}

inline int Lattice::get_cube_count() {
    return cube_vector.size();
}
inline int Lattice::get_spin(const Edge& edg) {
    return g[edg].spin;
}

inline double Lattice::get_potential_edge_energy(const Edge& edg) {
    return g[edg].integrated_edge_energy;
}

inline void Lattice::set_potential_edge_energy(const Edge& edg, double potential_energy) {
    g[edg].integrated_edge_energy = potential_energy;
}

inline void Lattice::add_potential_edge_energy(const Edge& edg, double diff) {
    g[edg].integrated_edge_energy += diff;
}

inline double Lattice::get_potential_star_energy(int star_index) {
    return g[star_index].integrated_star_energy;
}
    
inline void Lattice::set_potential_star_energy(int star_index, double potential_energy) {
    g[star_index].integrated_star_energy = potential_energy;
}

inline void Lattice::add_potential_star_energy(int star_index, double diff) {
    g[star_index].integrated_star_energy += diff;
}

inline double Lattice::get_potential_plaquette_energy(int plaquette_index) {
    return integrated_plaquette_energy_vector[plaquette_index];
}
    
inline void Lattice::set_potential_plaquette_energy(int plaquette_index, double potential_energy) {
    integrated_plaquette_energy_vector[plaquette_index] = potential_energy;
}

inline void Lattice::add_potential_plaquette_energy(int plaquette_index, double diff) {
    integrated_plaquette_energy_vector[plaquette_index] += diff;
}

inline std::string Lattice::get_orientation(const Edge& edg) {
    return g[edg].orientation;
}

inline double Lattice::get_spin_flip_imag_time(const Edge& edg, int spin_flip_index) {
    return g[edg].spin_flips[spin_flip_index];
}

inline void Lattice::set_spin_flip_imag_time(const Edge& edg, int spin_flip_index, double imag_time) {
    g[edg].spin_flips[spin_flip_index] = imag_time;
}

inline void Lattice::set_single_spin_flip_imag_time(const Edge& edg, int spin_flip_index, double imag_time) {
    g[edg].single_spin_flips[spin_flip_index] = imag_time;
}

inline int Lattice::get_spin_flip_count(const Edge& edg) {
    return g[edg].spin_flips.size();
}

inline Lattice::Edge Lattice::edge_in_between(int v_1, int v_2) {
    const auto edg_full = boost::edge(v_1, v_2, g);
#ifndef NDEBUG
    if (!edg_full.second) {
        throw std::runtime_error(std::format("There is no edge between vertex {} and {}.", v_1, v_2));
    }
#endif
    return edg_full.first;
}

inline std::span<const Lattice::Edge> Lattice::get_plaquette_edges(int p_index) {
    const auto& v = plaquette_edges_cache_[static_cast<size_t>(p_index)];
    return {v.data(), v.size()};
}

inline std::span<const Lattice::Edge> Lattice::get_star_edges(int center_index) {
    const auto& s = star_edges_cache_[static_cast<size_t>(center_index)];
    return {s.data(), s.size()};
}

inline bool Lattice::exists_edge(int v_1, int v_2) {
    return boost::edge(v_1, v_2, g).second;
}

inline std::pair<int, int> Lattice::vertices_of_edge(const Edge& edg) {
    const Edge e = edg;
    const auto source_v = boost::source(e, g);
    const auto target_v = boost::target(e, g);
    return {source_v, target_v};
}

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

    // -- find initial spin at imag_time_1 ----------------------------------
    const auto& spin_flips = g[edg].spin_flips;   // sorted flip times
    auto it_edge = std::lower_bound(
        spin_flips.begin(), spin_flips.end(), imag_time_1);

    int base_spin  = get_spin(edg);
    int spin_prev  = ((it_edge - spin_flips.begin()) & 1) ? -base_spin : base_spin;
    int spin_after = spin_prev;

    // -- merge two sorted streams: edge flips + external flips ------------
    size_t i_delta = 0;
    const size_t N_delta = spin_flip_lookup.size();
    double t_prev = imag_time_1;
    double energy_before = 0.0, energy_after = 0.0;

    constexpr double no_imaginary_time_found = -1.0;

    auto next_edge_time = [&]() -> double {
        return (it_edge != spin_flips.end() && *it_edge < imag_time_2)
             ? *it_edge
             : no_imaginary_time_found;
    };
    auto next_delta_time = [&]() -> double {
        return (i_delta < N_delta && spin_flip_lookup[i_delta].first < imag_time_2)
             ? spin_flip_lookup[i_delta].first
             : no_imaginary_time_found;
    };

    while (true) {
        double t_edge  = next_edge_time();
        double t_delta = next_delta_time();
        double t_curr;

        if (t_edge < 0 && t_delta < 0) {
            break;
        } else if (t_edge < 0) {
            t_curr = t_delta;
        } else if (t_delta < 0) {
            t_curr = t_edge;
        } else {
            t_curr = std::min(t_edge, t_delta);
        }

        double dt = t_curr - t_prev;
        energy_before += dt * spin_prev;
        energy_after  += dt * spin_after;
        t_prev = t_curr;

        // edge flip toggles both histories
        const bool edge_event = (t_edge >= 0) && (t_delta < 0 || t_edge <= t_delta);
        if (edge_event) {
            spin_prev  = -spin_prev;
            spin_after = -spin_after;
            ++it_edge;
        }

        // consume all lookup flips at this same time
        while (i_delta < N_delta && spin_flip_lookup[i_delta].first == t_curr) {
            spin_after = -spin_after;
            ++i_delta;
        }
    }

    // -- final segment until imag_time_2 ----------------------------------
    double dt_tail = imag_time_2 - t_prev;
    energy_before += dt_tail * spin_prev;
    energy_after  += dt_tail * spin_after;

    return energy_after - energy_before;
}

[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_edge_energy(
    const Edge& edg, double imag_time_1, double imag_time_2
) {
    if (imag_time_1 == imag_time_2) {
        throw std::invalid_argument(
            "integrated_edge_energy: time interval must be non-zero");
    }

    const auto& spin_flips = g[edg].spin_flips;  
    auto lo = std::lower_bound(
        spin_flips.begin(), spin_flips.end(), imag_time_1);
    auto hi = std::upper_bound(
        lo, spin_flips.end(), imag_time_2);

    // determine spin just after imag_time_1
    int base_spin = get_spin(edg);
    int spin      = (((lo - spin_flips.begin()) & 1)
                     ? -base_spin
                     :  base_spin);

    double energy = 0.0;
    double t_prev = imag_time_1;

    // accumulate each flip interval
    for (auto it = lo; it != hi; ++it) {
        double t_curr = *it;
        energy += (t_curr - t_prev) * spin;
        spin   = -spin;
        t_prev = t_curr;
    }

    // tail interval until imag_time_2
    if (t_prev < imag_time_2) {
        energy += (imag_time_2 - t_prev) * spin;
    }

    return energy;
}

[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_tuple_energy_diff_combination(
    std::span<const Edge> tuple_edges, 
    double imag_time_1, 
    double imag_time_2, 
    const std::vector<std::pair<double, int>>& spin_flip_lookup
) {
    struct Node { double t; std::vector<double>::const_iterator it, end; };
    auto cmp = [](const Node& a, const Node& b){ return a.t > b.t; };

    thread_local std::vector<Node> heap_tls;
    heap_tls.clear();
    heap_tls.reserve(tuple_edges.size());

    int spin_prod = 1;

    for (Edge e : tuple_edges) {
        const auto& flips = g[e].spin_flips;
        auto lo = std::lower_bound(flips.begin(), flips.end(), imag_time_1);
        auto hi = std::upper_bound(lo,           flips.end(), imag_time_2);

        if (lo != hi) heap_tls.push_back({*lo, lo, hi});

        int s = get_spin(e);
        if (std::distance(flips.begin(), lo) & 1) s = -s;
        spin_prod *= s;
    }
    std::make_heap(heap_tls.begin(), heap_tls.end(), cmp);

    // VERY IMPORTANT: This tuple here does not refer to the input "tuple_edges" but instead to the tuple where the tuple flip takes place
    // In principle, we don't need it for 2D but for special 3D cases, like when looking at cubic 12-spin interactions which share 3 bonds with a neighboring star.
    // Remember that spin_flip_lookup contains the plaquette flip (once) that is marked with 1 as well as the single spin flips. These
    // single spin flips all take place on the tuple where the tuple flips take place, i.e. the number of these tells us the number links
    // we share with the tuple where the tuple_flip takes place. This is the number we need.
    const bool odd_tuple = ((spin_flip_lookup.size() & 1) == 0);
    double E_before = 0., E_after  = 0.;
    int    spin_before = spin_prod, spin_after = spin_prod;
    double t_prev = imag_time_1;
    size_t idx_delta = 0;

    while (!heap_tls.empty() || idx_delta < spin_flip_lookup.size()) {
        double t_next;
        bool edge_event = false;

        if (!heap_tls.empty() 
        && (idx_delta == spin_flip_lookup.size() || heap_tls.front().t <= spin_flip_lookup[idx_delta].first)) {
            edge_event = true;
            t_next     = heap_tls.front().t;
        } else {
            t_next = spin_flip_lookup[idx_delta].first;
        }

        double dt = t_next - t_prev;
        E_before += dt * spin_before;
        E_after  += dt * spin_after;
        t_prev    = t_next;

        if (edge_event) {
            /* flip both streams */
            spin_before = -spin_before;
            spin_after  = -spin_after;

            /* pop root */
            std::pop_heap(heap_tls.begin(), heap_tls.end(), cmp);
            Node n = heap_tls.back();
            heap_tls.pop_back();

            if (++n.it != n.end) {
                n.t = *n.it;
                heap_tls.push_back(n);
                std::push_heap(heap_tls.begin(), heap_tls.end(), cmp);
            }
        }

        while (idx_delta < spin_flip_lookup.size() 
        && spin_flip_lookup[idx_delta].first == t_next) {
            if (spin_flip_lookup[idx_delta].second != 1 || odd_tuple)
                spin_after = -spin_after;
            ++idx_delta;
        }
    }
    return E_after - E_before;
}

[[gnu::hot, gnu::always_inline]]
inline double Lattice::integrated_tuple_energy(
    std::span<const Edge> tuple_edges, double imag_time_1, double imag_time_2
) {
    struct Node {
        double t;
        std::vector<double>::const_iterator it, end;
    };
    auto cmp = [](const Node& a, const Node& b){
        return a.t > b.t;
    };

    // build a min‑heap over the first flip of each edge in the tuple
    thread_local std::vector<Node> heap_tls;
    heap_tls.clear();
    heap_tls.reserve(tuple_edges.size());

    int spin_prod = 1;
    for (auto const& edg : tuple_edges) {
        auto const& spin_flips = g[edg].spin_flips;
        auto lo = std::lower_bound(spin_flips.begin(),
                                   spin_flips.end(),
                                   imag_time_1);
        auto hi = std::upper_bound(lo,
                                   spin_flips.end(),
                                   imag_time_2);

        if (lo != hi) {
            heap_tls.push_back({ *lo, lo, hi });
        }

        int before_count = static_cast<int>(lo - spin_flips.begin());
        int s = get_spin(edg);
        if (before_count & 1) {
            s = -s;
        }
        spin_prod *= s;
    }
    std::make_heap(heap_tls.begin(), heap_tls.end(), cmp);

    // sweep through all flip times in ascending order
    double energy_before = 0.0;
    double t_prev        = imag_time_1;
    int    spin          = spin_prod;

    while (!heap_tls.empty()) {
        // pop the earliest event
        std::pop_heap(heap_tls.begin(), heap_tls.end(), cmp);
        Node node = heap_tls.back();
        heap_tls.pop_back();

        double t_curr = node.t;
        double dt     = t_curr - t_prev;
        energy_before += dt * spin;
        spin = -spin;                // flip spin for next interval

        // advance this edge’s iterator
        ++node.it;
        if (node.it != node.end) {
            node.t = *node.it;
            heap_tls.push_back(node);
            std::push_heap(heap_tls.begin(), heap_tls.end(), cmp);
        }

        t_prev = t_curr;
    }

    // final segment up to imag_time_2
    if (t_prev < imag_time_2) {
        energy_before += (imag_time_2 - t_prev) * spin;
    }

    return energy_before;
}

} // namespace paratoric
