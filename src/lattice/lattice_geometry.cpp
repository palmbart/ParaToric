// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#include "lattice/lattice_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <initializer_list>
#include <stdexcept>

namespace paratoric {
namespace {

using Coordinates = LatticeGeometry::Coordinates;
using VertexPair = LatticeGeometry::VertexPair;

int modulo(int a, int b) {
    const int result = a % b;
    return result < 0 ? result + b : result;
}

// A vertex of a translated unit cell: offsets from the cell origin and a
// site within its basis. Staggered triangular/kagome rows use two motifs,
// the two halves of their translational unit cell, to retain vertex numbering.
struct Site {
    int dx = 0, dy = 0, dz = 0, basis = 0;
};

void add_edge(LatticeGeometry& g, int source, int target, const char* orientation = "x") {
    if (source >= 0 && target >= 0) {
        g.edges.push_back({source, target, orientation});
    }
}

void add_plaquette(LatticeGeometry& g, const std::vector<int>& vertices, Coordinates coordinates) {
    if (std::ranges::find(vertices, -1) != vertices.end()) return;
    LatticeGeometry::Path boundary;
    boundary.reserve(vertices.size());
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        boundary.emplace_back(vertices[i], vertices[(i + 1) % vertices.size()]);
    }
    g.plaquettes.push_back(std::move(boundary));
    g.plaquette_coordinates.push_back(coordinates);
}

struct Cell {
    LatticeGeometry& geometry;
    int L, depth, sites;
    bool periodic;
    int x, y, z = 0;

    int vertex(Site s = {}) const {
        const int vx = x + s.dx, vy = y + s.dy, vz = z + s.dz;
        if (!periodic && (vx < 0 || vx >= L || vy < 0 || vy >= L || vz < 0 || vz >= depth)) {
            return -1;
        }
        return sites * (modulo(vz, depth) * L * L + modulo(vy, L) * L + modulo(vx, L)) + s.basis;
    }

    void edge(Site from, Site to, const char* orientation = "x") const {
        add_edge(geometry, vertex(from), vertex(to), orientation);
    }

    void face(std::initializer_list<Site> sites_in_face, Coordinates coordinates) const {
        std::vector<int> vertices;
        vertices.reserve(sites_in_face.size());
        for (const auto& site : sites_in_face) vertices.push_back(vertex(site));
        add_plaquette(geometry, vertices, coordinates);
    }
};

void make_cartesian(LatticeGeometry& g, int L, bool periodic, bool cubic) {
    const int depth = cubic ? L : 1;
    g.dimensionality = cubic ? 3 : 2;
    g.vertices.resize(L * L * depth);
    g.edges.reserve(g.dimensionality * g.vertices.size());
    for (int x = 0; x < L; ++x) {
        for (int y = 0; y < L; ++y) {
            for (int z = 0; z < depth; ++z) {
                const Cell cell{g, L, depth, 1, periodic, x, y, z};
                g.vertices[cell.vertex()] = {double(x), double(y), double(z)};
                cell.edge({}, {1, 0, 0}, "x");
                cell.edge({}, {0, 1, 0}, "y");
                if (cubic) cell.edge({}, {0, 0, 1}, "z");
            }
        }
    }

    const auto faces = [&](int x, int y, int z) {
        const Cell cell{g, L, depth, 1, periodic, x, y, z};
        const Coordinates coordinates{double(x), double(y), double(z)};
        cell.face({{}, {1, 0}, {1, 1}, {0, 1}}, coordinates);
        if (!cubic) return;
        cell.face({{}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1}}, coordinates);
        cell.face({{}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}}, coordinates);
        // Preserve the existing open-cube convention: clip y/z, wrap x.
        if (periodic || (y < L - 1 && z < L - 1)) {
            const Cell wrapped{g, L, depth, 1, true, x, y, z};
            g.cubes.push_back({wrapped.vertex(), wrapped.vertex({0, 1}),
                wrapped.vertex({0, 0, 1}), wrapped.vertex({1, 0}),
                wrapped.vertex({0, 1, 1}), wrapped.vertex({1, 0, 1}),
                wrapped.vertex({1, 1}), wrapped.vertex({1, 1, 1})});
            g.cube_coordinates.push_back(coordinates);
        }
    };
    // Square faces historically follow vertex order; cubic faces follow x/y/z.
    if (cubic) {
        for (int x = 0; x < L; ++x)
            for (int y = 0; y < L; ++y)
                for (int z = 0; z < L; ++z) faces(x, y, z);
    } else {
        for (int y = 0; y < L; ++y)
            for (int x = 0; x < L; ++x) faces(x, y, 0);
    }
}

void make_triangular_or_kagome(LatticeGeometry& g, int L, bool periodic, bool kagome) {
    const int sites = kagome ? 3 : 1;
    const double height = std::sqrt(3) / 2.;
    const std::array<Coordinates, 3> basis{{{0., 0., 0.}, {0.5, 0., 0.}, {0.25, 0.5 * height, 0.}}};
    g.vertices.resize(sites * L * L);
    g.edges.reserve((kagome ? 6 : 3) * L * L);
    for (int x = 0; x < L; ++x) {
        for (int y = 0; y < L; ++y) {
            const Cell cell{g, L, 1, sites, periodic, x, y};
            for (int s = 0; s < sites; ++s) {
                g.vertices[cell.vertex({0, 0, 0, s})] = {
                    (y % 2) / 2.0 + x + basis[s][0], height * y + basis[s][1], 0.};
            }
            if (kagome) {
                cell.edge({}, {0, 0, 0, 1});
                cell.edge({0, 0, 0, 1}, {0, 0, 0, 2});
                cell.edge({0, 0, 0, 2}, {});
                cell.edge({0, 0, 0, 1}, {1, 0});
                cell.edge({0, 0, 0, 2}, {0, 1, 0, y % 2});
                cell.edge({0, 0, 0, 2}, {y % 2 ? 1 : -1, 1, 0, 1 - y % 2});
            } else {
                cell.edge({}, {1, 0});
                cell.edge({}, {0, 1});
                cell.edge({}, {y % 2 ? 1 : -1, 1});
            }
        }
    }
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const Cell cell{g, L, 1, sites, periodic, x, y};
            if (kagome) {
                // KAGOME TRIANGULAR PLAQUETTES: keep this block separate from
                // the hexagons below for future independent plaquette terms.
                // Percolation coordinates remain unspecified (zero).
                cell.face({{}, {0, 0, 0, 1}, {0, 0, 0, 2}}, {});
                if (y % 2 == 0) {
                    cell.face({{0, 0, 0, 2}, {0, 1}, {-1, 1, 0, 1}}, {});
                } else {
                    cell.face({{0, 0, 0, 2}, {1, 1}, {0, 1, 0, 1}}, {});
                }
            } else {
                const Coordinates coordinates{double(x), double(y), 0.};
                cell.face({{}, {1, 0}, {y % 2, 1}}, coordinates);
                cell.face({{}, {1, 0}, {y % 2, -1}}, coordinates);
            }
        }
    }
    if (kagome) {
        // KAGOME HEXAGONAL PLAQUETTES: appended after all triangles so their
        // indices stay unchanged. Both shapes currently share the plaquette
        // registry and coupling J; distinguish them by boundary size (3/6)
        // when separating their terms in the future.
        for (int y = 0; y < L; ++y) {
            for (int x = 0; x < L; ++x) {
                const Cell cell{g, L, 1, 3, periodic, x, y};
                cell.face({{0, 0, 0, 1}, {1, 0}, {1, 0, 0, 2},
                    {y % 2, 1, 0, 1}, {y % 2, 1}, {0, 0, 0, 2}}, {});
            }
        }
    }
}

void make_honeycomb(LatticeGeometry& g, int L, bool periodic) {
    // Brick-wall numbering of two-site cells. Open boundaries add a closing
    // row/column; the first and last rows each omit one boundary site.
    const auto row_width = [&](int row) {
        return periodic ? 2 * L : 2 * L + ((row == 0 || row == L) ? 1 : 2);
    };
    const auto vertex = [&](int col, int row) {
        if (periodic) return modulo(row, L) * 2 * L + modulo(col, 2 * L);
        return row == 0 ? col : (2 * L + 1) + (row - 1) * (2 * L + 2) + col;
    };
    const int rows = periodic ? L : L + 1;
    const double height = std::sqrt(3) / 2.;
    g.vertices.resize(periodic ? 2 * L * L : (2 * L + 2) * (L + 1) - 2);
    g.edges.reserve(3 * L * L);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < row_width(row); ++col) {
            const int physical_col = col + (!periodic && row == L && L % 2 == 0 ? 1 : 0);
            // Integer row/2 is intentional: these are the original coordinates.
            g.vertices[vertex(col, row)] = {physical_col * height,
                0.5 + row + row / 2 + (physical_col % 2) * ((row % 2) - 0.5), 0.};
            if (periodic || col + 1 < row_width(row)) {
                add_edge(g, vertex(col, row), vertex(col + 1, row));
            }
        }
    }
    // Horizontal bonds precede all vertical bonds in the historical ordering.
    for (int row = 0; row < L; ++row) {
        for (int col = row % 2; col < row_width(row); col += 2) {
            const int next_col = col - (!periodic && row == L - 1 && L % 2 == 0 ? 1 : 0);
            add_edge(g, vertex(col, row), vertex(next_col, row + 1));
        }
    }
    // One hexagon per two-site cell; alternate the top vertex by row parity.
    for (int row = 0; row < L; ++row) {
        for (int cell_x = 0; cell_x < L; ++cell_x) {
            const int col = 2 * cell_x + (periodic ? (row + 1) % 2 : 1 + row % 2);
            const int lower_col = col - (!periodic && row == L - 1 && L % 2 == 0 ? 1 : 0);
            add_plaquette(g, {vertex(col, row), vertex(col + 1, row),
                vertex(lower_col + 1, row + 1), vertex(lower_col, row + 1),
                vertex(lower_col - 1, row + 1), vertex(col - 1, row)},
                {double(periodic ? cell_x : col), double(row), 0.});
        }
    }
}

void finish_geometry(LatticeGeometry& g) {
    g.max_coordinates.assign(g.dimensionality, 0.);
    g.max_plaquette_coordinates.assign(g.dimensionality, 0.);
    for (int axis = 0; axis < g.dimensionality; ++axis) {
        for (const auto& v : g.vertices) g.max_coordinates[axis] = std::max(g.max_coordinates[axis], v[axis]);
        for (const auto& p : g.plaquette_coordinates)
            g.max_plaquette_coordinates[axis] = std::max(g.max_plaquette_coordinates[axis], p[axis]);
    }

    // Only cubes incident to a face's first vertex can contain that face.
    // Preserve the old endpoint-containment rule, including small tori where
    // several cubes have identical vertex sets, and ascending incidence order.
    std::vector<std::vector<int>> vertex_cubes(g.vertices.size());
    for (std::size_t c = 0; c < g.cubes.size(); ++c) {
        auto vertices = g.cubes[c];
        std::ranges::sort(vertices);
        const auto duplicates = std::ranges::unique(vertices);
        vertices.erase(duplicates.begin(), duplicates.end());
        for (int v : vertices) vertex_cubes[v].push_back(static_cast<int>(c));
    }
    g.plaquette_cubes.resize(g.plaquettes.size());
    g.cube_plaquettes.resize(g.cubes.size());
    for (std::size_t p = 0; p < g.plaquettes.size(); ++p) {
        for (int c : vertex_cubes[g.plaquettes[p].front().first]) {
            const auto& vertices = g.cubes[c];
            if (std::ranges::all_of(g.plaquettes[p], [&](VertexPair edge) {
                return std::ranges::find(vertices, edge.first) != vertices.end()
                    && std::ranges::find(vertices, edge.second) != vertices.end();
            })) {
                g.plaquette_cubes[p].push_back(c);
                g.cube_plaquettes[c].push_back(static_cast<int>(p));
            }
        }
    }
}

} // namespace

LatticeGeometry make_lattice_geometry(const std::string& lattice_type, int L, const std::string& boundaries) {
    if (lattice_type != "square" && lattice_type != "cubic" && lattice_type != "triangular"
        && lattice_type != "honeycomb" && lattice_type != "kagome") {
        throw std::invalid_argument(std::format("Lattice type \"{}\" is not supported.", lattice_type));
    }
    if (L <= 0) throw std::invalid_argument("System size has to be strictly positive.");
    if (boundaries != "periodic" && boundaries != "open") {
        throw std::invalid_argument("Boundaries can be \"periodic\" or \"open\".");
    }
    const bool periodic = boundaries == "periodic";
    if (periodic && L % 2) {
        if (lattice_type == "triangular") {
            throw std::invalid_argument("Periodic triangular lattice has to have an even system size.");
        }
        if (lattice_type == "honeycomb") {
            throw std::invalid_argument("Periodic honeycomb lattice is only supported for even system size.");
        }
        if (lattice_type == "kagome") {
            throw std::invalid_argument("Periodic Kagome lattice has to have an even system size.");
        }
    }
    LatticeGeometry geometry;
    if (lattice_type == "square" || lattice_type == "cubic") {
        make_cartesian(geometry, L, periodic, lattice_type == "cubic");
    } else if (lattice_type == "honeycomb") {
        make_honeycomb(geometry, L, periodic);
    } else {
        make_triangular_or_kagome(geometry, L, periodic, lattice_type == "kagome");
    }
    finish_geometry(geometry);
    return geometry;
}

std::pair<LatticeGeometry::Path, LatticeGeometry::Path> make_fredenhagen_marcu_paths(
    const LatticeGeometry& geometry, const std::string& lattice_type,
    int L, const std::string& boundaries, char basis
) {
    const double max_x = geometry.max_coordinates[0];
    const double max_y = geometry.max_coordinates[1];
    const int start_y = static_cast<int>(max_y / 4.0);
    const int end_y = static_cast<int>(3 * max_y / 4.0);
    const int middle_y = static_cast<int>((start_y + end_y) / 2.0);
    const int start_x = static_cast<int>(max_x / 4.0);
    const int end_x = static_cast<int>(3 * max_x / 4.0);
    std::vector<VertexPair> half_loop;
    std::vector<VertexPair> full_loop;

    if (basis == 'z') {
        if (lattice_type == "square" || lattice_type == "cubic" || lattice_type == "triangular") {
            const int z_offset = geometry.dimensionality == 3
                ? static_cast<int>(geometry.max_coordinates[2] / 2.0) * L * L : 0;
            int previous = z_offset + middle_y * L + start_x;
            const auto append = [&](int x, int y) {
                const int next = z_offset + y * L + x;
                full_loop.emplace_back(previous, next);
                previous = next;
            };
            for (int y = middle_y + 1; y <= end_y; ++y) append(start_x, y);
            for (int x = start_x + 1; x <= end_x; ++x) append(x, end_y);
            for (int y = end_y - 1; y >= middle_y; --y) append(end_x, y);
            half_loop = full_loop;
            for (int y = middle_y - 1; y >= start_y; --y) append(end_x, y);
            for (int x = end_x - 1; x >= start_x; --x) append(x, start_y);
            for (int y = start_y + 1; y <= middle_y; ++y) append(start_x, y);
        } else if (lattice_type == "honeycomb") {
            if (boundaries == "periodic") {
                int prev_vertex, next_vertex;

                prev_vertex = (2 * L) * middle_y + start_x + 1;
                if (middle_y%2 == 0 && prev_vertex%2 == 1) {
                    prev_vertex = (2 * L) * middle_y + start_x;
                } else if (middle_y%2 == 1 && prev_vertex%2 == 0) {
                    prev_vertex = (2 * L) * middle_y + start_x;
                }

                for (int y = middle_y; y < end_y; y+=2) {
                    next_vertex = prev_vertex + 2 * L;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;

                    next_vertex = prev_vertex + 1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                next_vertex = prev_vertex + 1;
                full_loop.emplace_back(prev_vertex, next_vertex);
                prev_vertex = next_vertex;

                for (int x = start_x; x < end_x-1; x+=2) {
                    next_vertex = prev_vertex+1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;

                    next_vertex = prev_vertex+1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                for (int y = end_y; y > middle_y; y-=2) {
                    next_vertex = prev_vertex - 2 * L;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;

                    next_vertex = prev_vertex + 1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                // Store half_loop after half of the path
                half_loop = full_loop;

                for (int y = middle_y; y > start_y; y-=2) {
                    next_vertex = prev_vertex - 2 * L;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;

                    next_vertex = prev_vertex - 1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                next_vertex = prev_vertex - 1;
                full_loop.emplace_back(prev_vertex, next_vertex);
                prev_vertex = next_vertex;

                for (int x = end_x; x > start_x+1; x-=2) {
                    next_vertex = prev_vertex-1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;

                    next_vertex = prev_vertex-1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                for (int y = start_y; y < middle_y; y+=2) {
                    next_vertex = prev_vertex + 2 * L;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;

                    next_vertex = prev_vertex - 1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }
            } else {
                if (L < 3) {
#ifndef NDEBUG
                    throw std::runtime_error(
                        std::format(
                            "System size for open honeycomb Wilson loops has to be at least L=3 but is L={}.",
                            L
                        )
                    );
#endif
                }

                auto honeycomb_vertex = [L](int row, int col) {
                    if (row == 0) {
                        return col;
                    }
                    return (2 * L + 1)
                        + (row - 1) * (2 * L + 2)
                        + col;
                };

                auto vertical_down_col = [L](int row, int col) {
                    if (row == L - 1 && L % 2 == 0) {
                        return col - 1;
                    }
                    return col;
                };

                auto vertical_up_col = [L](int row, int col) {
                    if (row == L && L % 2 == 0) {
                        return col + 1;
                    }
                    return col;
                };

                int current_row = middle_y;
                int current_col = start_x + ((start_x & 1) == (middle_y & 1) ? 0 : 1);

                auto append_horizontal = [&](int step) {
                    const int next_col = current_col + step;
                    full_loop.emplace_back(
                        honeycomb_vertex(current_row, current_col),
                        honeycomb_vertex(current_row, next_col)
                    );
                    current_col = next_col;
                };

                auto append_vertical_down = [&]() {
                    const int next_col = vertical_down_col(current_row, current_col);
                    full_loop.emplace_back(
                        honeycomb_vertex(current_row, current_col),
                        honeycomb_vertex(current_row + 1, next_col)
                    );
                    ++current_row;
                    current_col = next_col;
                };

                auto append_vertical_up = [&]() {
                    const int next_col = vertical_up_col(current_row, current_col);
                    full_loop.emplace_back(
                        honeycomb_vertex(current_row, current_col),
                        honeycomb_vertex(current_row - 1, next_col)
                    );
                    --current_row;
                    current_col = next_col;
                };

                for (int y = middle_y; y < end_y; y += 2) {
                    append_vertical_down();
                    append_horizontal(1);
                }

                append_horizontal(1);

                for (int x = start_x; x < end_x - 1; x += 2) {
                    append_horizontal(1);
                    append_horizontal(1);
                }

                for (int y = end_y; y > middle_y; y -= 2) {
                    append_vertical_up();
                    append_horizontal(1);
                }

                half_loop = full_loop;

                for (int y = middle_y; y > start_y; y -= 2) {
                    append_vertical_up();
                    append_horizontal(-1);
                }

                append_horizontal(-1);

                for (int x = end_x; x > start_x + 1; x -= 2) {
                    append_horizontal(-1);
                    append_horizontal(-1);
                }

                for (int y = start_y; y < middle_y; y += 2) {
                    append_vertical_down();
                    append_horizontal(-1);
                }
            }

        } else if (lattice_type == "kagome") {
            // TODO: Construct Wilson paths for kagome geometry.
#ifndef NDEBUG
            throw std::invalid_argument("Kagome lattice is not supported for Wilson loops.");
#endif
            return std::make_pair(half_loop, full_loop);
        }
    } else if (basis == 'x') {
        if (L < 6) {
#ifndef NDEBUG
            throw std::runtime_error(
                std::format(
                    "System size for 't Hooft loops in x-basis has to be at least L=6 but is L={}.",
                    L
                )
            );
#endif
            return std::make_pair(half_loop, full_loop);
        }

        if (lattice_type == "square") {
            int start_vertex = 0;

            for (int y = middle_y; y < end_y+1; ++y) {
                start_vertex = y * L + start_x;
                full_loop.emplace_back(start_vertex, start_vertex-1);
            }

            for (int x = start_x; x < end_x+1; ++x) {
                start_vertex = end_y * L + x;
                full_loop.emplace_back(start_vertex, start_vertex+L);
            }

            for (int y = end_y; y > middle_y-1; --y) {
                start_vertex = y * L + end_x;
                full_loop.emplace_back(start_vertex, start_vertex+1);
            }

            // Store half_loop after half of the path
            half_loop = full_loop;

            for (int y = middle_y-1; y > start_y-1; --y) {
                start_vertex = y * L + end_x;
                full_loop.emplace_back(start_vertex, start_vertex+1);
            }

            for (int x = end_x; x > start_x-1; --x) {
                start_vertex = start_y * L + x;
                full_loop.emplace_back(start_vertex, start_vertex-L);
            }

            for (int y = start_y; y < middle_y+1; ++y) {
                start_vertex = y * L + start_x;
                full_loop.emplace_back(start_vertex, start_vertex-1);
            }
        } else if (lattice_type == "triangular") {
            int start_vertex = 0;
            int end_y_tri = end_y;
            int start_y_tri = start_y;

            if (end_y % 2 == 1) end_y_tri = end_y + 1;
            if (start_y % 2 == 0) start_y_tri = start_y + 1;

            for (int y = middle_y; y < end_y_tri+1; ++y) {
                if (y%2 == 0) {
                    start_vertex = y * L + start_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex+1, start_vertex+L);
                } else {
                    start_vertex = y * L + start_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex, start_vertex+L+1);
                }
            }

            for (int x = start_x+1; x < end_x; ++x) {
                start_vertex = end_y_tri * L + x;
                full_loop.emplace_back(start_vertex, start_vertex+L);

                full_loop.emplace_back(start_vertex+1, start_vertex+L);
            }

            // Last bond
            start_vertex = end_y_tri * L + end_x;
            full_loop.emplace_back(start_vertex, start_vertex+L);

            for (int y = end_y_tri; y > middle_y-1; --y) {
                if (y%2 == 0) {
                    start_vertex = y * L + end_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex+1, start_vertex-L);
                } else {
                    start_vertex = y * L + end_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex, start_vertex-L+1);
                }
            }

            // Store half_loop after half of the path
            half_loop = full_loop;

            for (int y = middle_y-1; y > start_y_tri-1; --y) {
                if (y%2 == 0) {
                    start_vertex = y * L + end_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex+1, start_vertex-L);
                } else {
                    start_vertex = y * L + end_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex, start_vertex-L+1);
                }
            }

            for (int x = end_x; x > start_x+1; --x) {
                start_vertex = start_y_tri * L + x;
                full_loop.emplace_back(start_vertex, start_vertex-L);

                full_loop.emplace_back(start_vertex-1, start_vertex-L);
            }

            // Last bond
            start_vertex = start_y_tri * L + start_x + 1;
            full_loop.emplace_back(start_vertex, start_vertex-L);

            for (int y = start_y_tri; y < middle_y; ++y) {
                if (y%2 == 0) {
                    start_vertex = y * L + start_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex+1, start_vertex+L);
                } else {
                    start_vertex = y * L + start_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex, start_vertex+L+1);
                }
            }

        } else if (lattice_type == "honeycomb") {
            const int plaquette_max = L - 1;
            const int start_y_hc = plaquette_max / 4;
            const int end_y_hc = (3 * plaquette_max) / 4;
            const int middle_y_hc = (start_y_hc + end_y_hc) / 2;
            const int start_x_hc = plaquette_max / 4;
            const int end_x_hc = (3 * plaquette_max) / 4;

            full_loop.reserve(
                2 * (end_y_hc - start_y_hc) + 2 * (end_x_hc - start_x_hc)
            );

            auto plaquette_index = [](int L, int y, int x) {
                return y * L + x;
            };

            auto append_common_edge = [&](int plaquette_1, int plaquette_2) {
                const auto& edges_1 = geometry.plaquettes[plaquette_1];
                const auto& edges_2 = geometry.plaquettes[plaquette_2];

                for (const auto& edge_1 : edges_1) {
                    const int source_1 = edge_1.first;
                    const int target_1 = edge_1.second;
                    for (const auto& edge_2 : edges_2) {
                        if ((source_1 == edge_2.first && target_1 == edge_2.second)
                            || (source_1 == edge_2.second && target_1 == edge_2.first)) {
                            full_loop.emplace_back(source_1, target_1);
                            return;
                        }
                    }
                }

#ifndef NDEBUG
                throw std::runtime_error("Honeycomb 't Hooft loop crosses non-adjacent plaquettes.");
#endif
            };

            int prev_plaquette = plaquette_index(L, middle_y_hc, start_x_hc);
            int next_plaquette;

            for (int y = middle_y_hc + 1; y < end_y_hc + 1; ++y) {
                next_plaquette = plaquette_index(L, y, start_x_hc);
                append_common_edge(prev_plaquette, next_plaquette);
                prev_plaquette = next_plaquette;
            }

            for (int x = start_x_hc + 1; x < end_x_hc + 1; ++x) {
                next_plaquette = plaquette_index(L, end_y_hc, x);
                append_common_edge(prev_plaquette, next_plaquette);
                prev_plaquette = next_plaquette;
            }

            for (int y = end_y_hc - 1; y > middle_y_hc - 1; --y) {
                next_plaquette = plaquette_index(L, y, end_x_hc);
                append_common_edge(prev_plaquette, next_plaquette);
                prev_plaquette = next_plaquette;
            }

            half_loop = full_loop;

            for (int y = middle_y_hc - 1; y > start_y_hc - 1; --y) {
                next_plaquette = plaquette_index(L, y, end_x_hc);
                append_common_edge(prev_plaquette, next_plaquette);
                prev_plaquette = next_plaquette;
            }

            for (int x = end_x_hc - 1; x > start_x_hc - 1; --x) {
                next_plaquette = plaquette_index(L, start_y_hc, x);
                append_common_edge(prev_plaquette, next_plaquette);
                prev_plaquette = next_plaquette;
            }

            for (int y = start_y_hc + 1; y < middle_y_hc + 1; ++y) {
                next_plaquette = plaquette_index(L, y, start_x_hc);
                append_common_edge(prev_plaquette, next_plaquette);
                prev_plaquette = next_plaquette;
            }
        } else if (lattice_type == "kagome") {
            // TODO: Construct dual 't Hooft paths for kagome geometry.
#ifndef NDEBUG
            throw std::invalid_argument("Kagome lattice is not supported for 't Hooft loops.");
#endif
            return std::make_pair(half_loop, full_loop);
        }
    }

    return std::make_pair(half_loop, full_loop);
}

} // namespace paratoric
