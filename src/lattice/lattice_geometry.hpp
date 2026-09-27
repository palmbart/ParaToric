// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#pragma once

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace paratoric {

/**
 * Geometry only: no graph descriptors, spin state, RNG, or QMC dependencies.
 * Indices are positions in the corresponding vectors. Edge insertion order,
 * directed endpoint pairs, and plaquette/cube order preserve the lattice's
 * historical numbering, including parallel edges in small periodic systems.
 */
struct LatticeGeometry {
    using VertexPair = std::pair<int, int>;
    using Path = std::vector<VertexPair>;
    using Coordinates = std::array<double, 3>;

    struct Edge {
        int source;
        int target;
        std::string orientation = "x";
    };

    int dimensionality = 2;
    // Physical positions in vertex-index order; unused axes are zero.
    std::vector<Coordinates> vertices;
    // Insertion order, retaining multiplicity and endpoint direction.
    std::vector<Edge> edges;
    // Each boundary is an ordered, closed sequence of vertex pairs.
    std::vector<Path> plaquettes;
    // Historical percolation coordinates (cell labels, not face centroids).
    std::vector<Coordinates> plaquette_coordinates;
    std::vector<std::vector<int>> cubes;
    std::vector<Coordinates> cube_coordinates;
    // Ascending indices in each incidence list.
    std::vector<std::vector<int>> plaquette_cubes;
    std::vector<std::vector<int>> cube_plaquettes;
    std::vector<double> max_coordinates;
    std::vector<double> max_plaquette_coordinates;
};

/** Tile the supported unit cells, clipping open boundaries or wrapping periodic
 * ones. Throws std::invalid_argument for invalid geometry specifications.
 * Kagome appends hexagonal plaquettes after the triangles; both shapes retain
 * zero plaquette coordinates and can be distinguished by boundary size (3/6).
 * Open cubic cubes retain the historical wrap in x.
 */
LatticeGeometry make_lattice_geometry(
    const std::string& lattice_type, int system_size, const std::string& boundaries
);

/** Half/full Wilson (z) or dual 't Hooft (x) paths as vertex pairs.
 * The geometry must match the supplied type, size, and boundaries.
 * Kept separate so other lattice classes can use geometry without these
 * observables. Preserves the existing size/support checks in debug builds.
 */
std::pair<LatticeGeometry::Path, LatticeGeometry::Path> make_fredenhagen_marcu_paths(
    const LatticeGeometry& geometry, const std::string& lattice_type,
    int system_size, const std::string& boundaries, char basis
);

} // namespace paratoric
