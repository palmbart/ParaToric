// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#include "lattice/lattice.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>

#include <algorithm> 
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <span>
#include <sstream>
#include <stack>
#include <tuple>
#include <utility>
#include <vector>

#define UNUSED(expr) do { (void)(expr); } while (0)

namespace paratoric {

void Lattice::check_input_validity() const {
    if (DEFAULT_SPIN != -1 && DEFAULT_SPIN != 1) {
        throw std::invalid_argument("Default spin has to be 1 or -1.");
    }

    if (BETA <= 0) {
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

Lattice::LatticeGraph Lattice::init_lattice_graph(
    char basis,
    const std::string& lattice_type, 
    int L,
    double beta,
    const std::string& boundaries,
    int default_spin
    ) {
    UNUSED(beta);

    if (lattice_type == "square" 
        || lattice_type == "triangular" 
        || lattice_type == "honeycomb" 
        || lattice_type == "kagome") {
        LATTICE_DIMENSIONALITY = 2;
    } else if (lattice_type == "cubic") {
        LATTICE_DIMENSIONALITY = 3;
    }

    if (lattice_type == "square") {
        std::vector<VertexPair> vertex_pair_vec;
        std::vector<std::string> edge_orientation_vec;
        vertex_pair_vec.reserve(2 * L * L);
        edge_orientation_vec.reserve(2 * L * L);

        g = LatticeGraph(L * L);
        if (BOUNDARIES == "periodic") {
            for (int x = 0; x < L; ++x) {
                for (int y = 0; y < L; ++y) {
                    for (int z = 0; z < 1; ++z) {
                        int vertex_index = z * L * L + y * L + x;
                        int vertex_index_right = z * L * L + y * L + (x + 1) % L;
                        int vertex_index_down = z * L * L + ((y + 1) % L) * L + x;
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_right );
                        edge_orientation_vec.emplace_back("x");
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_down );
                        edge_orientation_vec.emplace_back("y");
                    }
                }
            }
        } else {
            for (int x = 0; x < L; ++x) {
                for (int y = 0; y < L; ++y) {
                    for (int z = 0; z < 1; ++z) {
                        int vertex_index = z * L * L + y * L + x;
                        int vertex_index_right = z * L * L + y * L + (x + 1) % L;
                        int vertex_index_down = z * L * L + ((y + 1) % L) * L + x;

                        if (x < L - 1) {
                            vertex_pair_vec.emplace_back( vertex_index, vertex_index_right );
                            edge_orientation_vec.emplace_back("x");
                        }

                        if (y < L - 1) {
                            vertex_pair_vec.emplace_back( vertex_index, vertex_index_down );
                            edge_orientation_vec.emplace_back("y");
                        }
                    }
                }
            }
        }
        
        // Add edges to the graph
        for (size_t i = 0; i < vertex_pair_vec.size(); ++i) {
            boost::add_edge(
                vertex_pair_vec[i].first, vertex_pair_vec[i].second, 
                EdgeData(edge_orientation_vec[i]), 
                g
            );
        }

    } else if (lattice_type == "cubic") {
        std::vector<VertexPair> vertex_pair_vec;
        std::vector<std::string> edge_orientation_vec;
        vertex_pair_vec.reserve(3 * L * L * L);
        edge_orientation_vec.reserve(3 * L * L * L);

        g = LatticeGraph(L * L * L);
        if (BOUNDARIES == "periodic") {
            for (int x = 0; x < L; ++x) {
                for (int y = 0; y < L; ++y) {
                    for (int z = 0; z < L; ++z) {
                        int vertex_index = z * L * L + y * L + x;
                        int vertex_index_right = z * L * L + y * L + (x + 1) % L;
                        int vertex_index_down = z * L * L + ((y + 1) % L) * L + x;
                        int vertex_index_deep = ((z+1) % L) * L * L + y * L + x;
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_right );
                        edge_orientation_vec.emplace_back("x");
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_down );
                        edge_orientation_vec.emplace_back("y");
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_deep );
                        edge_orientation_vec.emplace_back("z");
                    }
                }
            }
        } else {
            for (int x = 0; x < L; ++x) {
                for (int y = 0; y < L; ++y) {
                    for (int z = 0; z < L; ++z) {
                        int vertex_index = z * L * L + y * L + x;
                        int vertex_index_right = z * L * L + y * L + (x + 1) % L;
                        int vertex_index_down = z * L * L + ((y + 1) % L) * L + x;
                        int vertex_index_deep = ((z+1) % L) * L * L + y * L + x;

                        if (x < L - 1) {
                            vertex_pair_vec.emplace_back( vertex_index, vertex_index_right );
                            edge_orientation_vec.emplace_back("x");
                        }

                        if (y < L - 1) {
                            vertex_pair_vec.emplace_back( vertex_index, vertex_index_down );
                            edge_orientation_vec.emplace_back("y");
                        }

                        if (z < L - 1) {
                            vertex_pair_vec.emplace_back( vertex_index, vertex_index_deep );
                            edge_orientation_vec.emplace_back("z");
                        }
                    }
                }
            }
        }

        // Add edges to the graph
        for (size_t i = 0; i < vertex_pair_vec.size(); ++i) {
            boost::add_edge(
                vertex_pair_vec[i].first, 
                vertex_pair_vec[i].second, 
                EdgeData(edge_orientation_vec[i]), 
                g
            );
        }
        
    } else if (lattice_type == "honeycomb") {
        std::vector<VertexPair> vertex_pair_vec;
        vertex_pair_vec.reserve(3 * L * L);

        if (BOUNDARIES == "periodic") {
            if (L % 2 == 1) {
                throw std::invalid_argument("Periodic honeycomb lattice is only supported for even system size.");
            }

            g = LatticeGraph(2 * L * L);

            //////////////////////
            // Horizontal edges //
            //////////////////////

            // First layer, it has one node less than the other layers!
            for (int i = 0; i < 2*L; ++i) {
                vertex_pair_vec.emplace_back( i, modulo(i+1, 2 * L) );
            }

            // Intermediate layers
            for (int y = 1; y < L; ++y) {
                for (int x = 0; x < 2*L; ++x) {
                    int vertex_index = y * 2 * L + x; // Note that first layer has one less node!
                    vertex_pair_vec.emplace_back( vertex_index, modulo(vertex_index+1, 2 * L) + y * (2 * L) );
                }
            }

            ////////////////////
            // Vertical edges //
            ////////////////////

            // First layer 
            for (int i = 0; i < 2*L; i+=2) {
                vertex_pair_vec.emplace_back( i, i+(2 * L) );
            }

            // Intermediate Layers (up until third (!!!) last layer)
            for (int y = 1; y < L-1; ++y) {
                for (int x = (y%2); x < 2*L; x+=2) {
                    int vertex_index = y * 2 * L + x;
                    vertex_pair_vec.emplace_back( vertex_index, vertex_index+(2 * L) );
                }
            }
            
            int count = 1;
            for (int x = (L - 1) % 2; x < 2*L; x+=2) {
                int y = (L - 1);
                int vertex_index_1 = (2 * L) + (y - 1) * (2 * L) + x;
                int vertex_index_2 = count;
                count += 2;
                
                vertex_pair_vec.emplace_back( vertex_index_1, vertex_index_2 );
            }
        } else {
            g = LatticeGraph((2 * L + 2) * (L + 1) - 2 );

            //////////////////////
            // Horizontal edges //
            //////////////////////

            // First layer, it has one node less than the other layers!
            for (int i = 0; i < 2*L; ++i) {
                vertex_pair_vec.emplace_back( i, i+1 );
            }

            // Intermediate layers
            for (int y = 1; y < L; ++y) {
                for (int x = 0; x < 2*L+1; ++x) {
                    // Note that first layer has one less node!
                    int vertex_index = (2 * L + 1) + (y - 1) * (2 * L + 2) + x; 
                    vertex_pair_vec.emplace_back( vertex_index, vertex_index+1 );
                }
            }

            // Last layer
            for (int x = 0; x < 2*L; ++x) {
                int y = L;
                int vertex_index = (2 * L + 1) + (y - 1) * (2 * L + 2) + x;
                vertex_pair_vec.emplace_back( vertex_index, vertex_index+1 );
            }

            ////////////////////
            // Vertical edges //
            ////////////////////

            // First layer 
            for (int i = 0; i < 2*L+1; i+=2) {
                vertex_pair_vec.emplace_back( i, i+(2 * L + 1) );
            }

            // Intermediate Layers (up until third (!!!) last layer)
            for (int y = 1; y < L-1; ++y) {
                for (int x = (y%2); x < 2*L+2; x+=2) {
                    int vertex_index = (2 * L + 1) + (y - 1) * (2 * L + 2) + x;
                    vertex_pair_vec.emplace_back( vertex_index, vertex_index+(2 * L + 2) );
                }
            }

            for (int x = (L - 1) % 2; x < 2*L+2; x+=2) {
                int y = (L - 1);
                int vertex_index_1 = (2 * L + 1) + (y - 1) * (2 * L + 2) + x;
                int vertex_index_2;
                if (L % 2 == 0) {
                    vertex_index_2 = vertex_index_1 + (2 * L + 1);
                } else {
                    vertex_index_2 = vertex_index_1 + (2 * L + 2);
                }
                vertex_pair_vec.emplace_back( vertex_index_1, vertex_index_2 );
            }
        }

        // Add edges to the graph
        for (const auto&[v1, v2] : vertex_pair_vec) {
            boost::add_edge(v1, v2, g);
        }

    } else if (lattice_type == "triangular") {
        std::vector<VertexPair> vertex_pair_vec;
        vertex_pair_vec.reserve(3 * L * L);

        if (BOUNDARIES == "periodic") {
            if (L % 2 == 1) {
                throw std::invalid_argument("Periodic triangular lattice has to have an even system size.");
            }
            g = LatticeGraph(L * L);

            for (int x = 0; x < L; ++x) {
                for (int y = 0; y < L; ++y) {
                    int vertex_index = y * L + x;
                    int vertex_index_right = y * L + (x + 1) % L;
                    int vertex_index_down_right = ((y+1) % L) * L + modulo(x + 1, L);
                    int vertex_index_down_mid = ((y+1) % L) * L + x;
                    int vertex_index_down_left = ((y+1) % L) * L + modulo(x - 1, L);
                    
                    vertex_pair_vec.emplace_back( vertex_index, vertex_index_right );
                    
                    if (y % 2 == 0) {
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_down_mid );
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_down_left );
                    } else {
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_down_mid );
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_down_right );
                    }  
                }
            }
        } else {
            g = LatticeGraph(L * L);

            for (int x = 0; x < L; ++x) {
                for (int y = 0; y < L; ++y) {
                    int vertex_index = y * L + x;
                    int vertex_index_right = y * L + (x + 1) % L;
                    int vertex_index_down_right = ((y+1) % L) * L + modulo(x + 1, L);
                    int vertex_index_down_mid = ((y+1) % L) * L + x;
                    int vertex_index_down_left = ((y+1) % L) * L + modulo(x - 1, L);
                    
                    if (x < L - 1) {
                        vertex_pair_vec.emplace_back( vertex_index, vertex_index_right );
                    }

                    if (y < L - 1) {
                        if (y % 2 == 0) {
                            vertex_pair_vec.emplace_back( vertex_index, vertex_index_down_mid );
                            if (x > 0) {
                                vertex_pair_vec.emplace_back( vertex_index, vertex_index_down_left );
                            }
                        } else {
                            vertex_pair_vec.emplace_back( vertex_index, vertex_index_down_mid );
                            if (x < L - 1) {
                                vertex_pair_vec.emplace_back( vertex_index, vertex_index_down_right );
                            }
                        }
                    }
                }
            }
        }

        // Add edges to the graph
        for (const auto&[v1, v2] : vertex_pair_vec) {
            boost::add_edge(v1, v2, g);
        }

    } else if (lattice_type == "kagome") {
        std::vector<VertexPair> vertex_pair_vec;
        vertex_pair_vec.reserve(6 * L * L);

        if (BOUNDARIES == "periodic") {
            if (L % 2 == 1) {
                throw std::invalid_argument("Periodic Kagome lattice has to have an even system size.");
            }

            g = LatticeGraph(3 * L * L);

            for (int x = 0; x < L; ++x) {
                for (int y = 0; y < L; ++y) {
                    int triangular_superlattice_index = y * L + x;

                    // Bonds in a triangle
                    int v_index_triangle_left = 3 * triangular_superlattice_index;
                    int v_index_triangle_right = v_index_triangle_left + 1;
                    int v_index_triangle_top = v_index_triangle_left + 2;
                    vertex_pair_vec.emplace_back( v_index_triangle_left, v_index_triangle_right );
                    vertex_pair_vec.emplace_back( v_index_triangle_right, v_index_triangle_top );
                    vertex_pair_vec.emplace_back( v_index_triangle_top, v_index_triangle_left );

                    //int vertex_index = y * L + x;
                    int triangular_superlattice_index_right = y * L + (x + 1) % L;
                    int v_index_right_triangle_left = 3 * triangular_superlattice_index_right;

                    vertex_pair_vec.emplace_back( v_index_triangle_right, v_index_right_triangle_left );

                    int triangular_superlattice_index_down_right = ((y+1) % L) * L + modulo(x + 1, L);
                    int triangular_superlattice_index_down_mid = ((y+1) % L) * L + x;
                    int triangular_superlattice_index_down_left = ((y+1) % L) * L + modulo(x - 1, L);

                    if (y % 2 == 0) {
                        int v_index_down_mid_triangle_left = 3 * triangular_superlattice_index_down_mid;
                        int v_index_down_left_triangle_right = 3 * triangular_superlattice_index_down_left + 1;
                        vertex_pair_vec.emplace_back( v_index_triangle_top, v_index_down_mid_triangle_left );
                        vertex_pair_vec.emplace_back( v_index_triangle_top, v_index_down_left_triangle_right );
                    } else {
                        int v_index_down_right_triangle_left = 3 * triangular_superlattice_index_down_right;
                        int v_index_down_mid_triangle_right = 3 * triangular_superlattice_index_down_mid + 1;
                        vertex_pair_vec.emplace_back( v_index_triangle_top, v_index_down_mid_triangle_right );
                        vertex_pair_vec.emplace_back( v_index_triangle_top, v_index_down_right_triangle_left );
                    }
                }
            }
        } else {
            g = LatticeGraph(3 * L * L);

            for (int x = 0; x < L; ++x) {
                for (int y = 0; y < L; ++y) {
                    int triangular_superlattice_index = y * L + x;

                    // Bonds in a triangle
                    int v_index_triangle_left = 3 * triangular_superlattice_index;
                    int v_index_triangle_right = v_index_triangle_left + 1;
                    int v_index_triangle_top = v_index_triangle_left + 2;
                    vertex_pair_vec.emplace_back( v_index_triangle_left, v_index_triangle_right );
                    vertex_pair_vec.emplace_back( v_index_triangle_right, v_index_triangle_top );
                    vertex_pair_vec.emplace_back( v_index_triangle_top, v_index_triangle_left );

                    //int vertex_index = y * L + x;
                    int triangular_superlattice_index_right = y * L + (x + 1) % L;
                    int v_index_right_triangle_left = 3 * triangular_superlattice_index_right;

                    if (x < L - 1) {
                        vertex_pair_vec.emplace_back( v_index_triangle_right, v_index_right_triangle_left );
                    }

                    int triangular_superlattice_index_down_right = ((y+1) % L) * L + modulo(x + 1, L);
                    int triangular_superlattice_index_down_mid = ((y+1) % L) * L + x;
                    int triangular_superlattice_index_down_left = ((y+1) % L) * L + modulo(x - 1, L);

                    if (y < L - 1) {
                        if (y % 2 == 0) {
                            int v_index_down_mid_triangle_left = 3 * triangular_superlattice_index_down_mid;
                            int v_index_down_left_triangle_right = 3 * triangular_superlattice_index_down_left + 1;
                            vertex_pair_vec.emplace_back( 
                                v_index_triangle_top, 
                                v_index_down_mid_triangle_left 
                            );
                            if (x > 0) {
                                vertex_pair_vec.emplace_back( 
                                    v_index_triangle_top, 
                                    v_index_down_left_triangle_right 
                                );
                            }
                        } else {
                            int v_index_down_right_triangle_left = 3 * triangular_superlattice_index_down_right;
                            int v_index_down_mid_triangle_right = 3 * triangular_superlattice_index_down_mid + 1;
                            vertex_pair_vec.emplace_back( 
                                v_index_triangle_top, 
                                v_index_down_mid_triangle_right 
                            );
                            if (x < L - 1) {
                                vertex_pair_vec.emplace_back( 
                                    v_index_triangle_top, 
                                    v_index_down_right_triangle_left 
                                );
                            }
                        }
                    }
                }
            }
        }

        // Add edges to the graph
        for (const auto&[v1, v2] : vertex_pair_vec) {
            boost::add_edge(v1, v2, g);
        }

    } else {
        throw std::runtime_error(std::format("Lattice type \"{}\" is not supported.", lattice_type));
    }

    // Assign coordinates
    if (lattice_type == "square") {
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                for (int z = 0; z < 1; ++z) {
                    int vertex_index = z * 1 * L + y * L + x;
                    g[vertex_index].x = x;
                    g[vertex_index].y = y;
                    g[vertex_index].z = z;
                }
            }
        }

        MAX_COORDINATES.emplace_back(L - 1 ); // x
        MAX_COORDINATES.emplace_back(L - 1 ); // y
        
    } else if (lattice_type == "cubic") {
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                for (int z = 0; z < L; ++z) {
                    int vertex_index = z * L * L + y * L + x;
                    g[vertex_index].x = x;
                    g[vertex_index].y = y;
                    g[vertex_index].z = z;
                }
            }
        }

        MAX_COORDINATES.emplace_back(L - 1 ); // x
        MAX_COORDINATES.emplace_back(L - 1 ); // y
        MAX_COORDINATES.emplace_back(L - 1 ); // z

    } else if (lattice_type == "honeycomb") {
        if (BOUNDARIES == "periodic") {
            const double h = std::sqrt(3) / 2.;
            std::vector<double> x_max_vec, y_max_vec;
            // First layer, it has one node less than the other layers!
            for (int i = 0; i < 2*L; ++i) {
                int x_old = i;
                int y_old = 0;
                double y_new = 0.5 + y_old + y_old / 2 + (x_old % 2) * ((y_old % 2) - 0.5);
                double x_new = x_old * h;
                g[i].x = x_new;
                g[i].y = y_new;
                x_max_vec.emplace_back(x_new);
                y_max_vec.emplace_back(y_new);
            }

            // Intermediate layers
            for (int y_old = 1; y_old < L; ++y_old) {
                for (int x_old = 0; x_old < 2*L; ++x_old) {
                    int vertex_index = y_old * 2 * L + x_old;
                    double y_new = 0.5 + y_old + y_old / 2 + (x_old % 2) * ((y_old % 2) - 0.5);
                    double x_new = x_old * h;
                    g[vertex_index].x = x_new;
                    g[vertex_index].y = y_new;
                    x_max_vec.emplace_back(x_new);
                    y_max_vec.emplace_back(y_new);
                }
            }

            auto x_max = std::max_element(x_max_vec.begin(), x_max_vec.end()); 
            auto y_max = std::max_element(y_max_vec.begin(), y_max_vec.end()); 
            MAX_COORDINATES.emplace_back( *x_max ); // x
            MAX_COORDINATES.emplace_back( *y_max ); // y
        } else {
            const double h = std::sqrt(3) / 2.;
            std::vector<double> x_max_vec, y_max_vec;
            // First layer, it has one node less than the other layers!
            for (int i = 0; i < 2*L+1; ++i) {
                int x_old = i;
                int y_old = 0;
                double y_new = 0.5 + y_old + y_old / 2 + (x_old % 2) * ((y_old % 2) - 0.5);
                double x_new = x_old * h;
                g[i].x = x_new;
                g[i].y = y_new;
                x_max_vec.emplace_back(x_new);
                y_max_vec.emplace_back(y_new);
            }

            // Intermediate layers
            for (int y_old = 1; y_old < L; ++y_old) {
                for (int x_old = 0; x_old < 2*L+2; ++x_old) {
                    int vertex_index = (2 * L + 1) + (y_old - 1) * (2 * L + 2) + x_old;
                    double y_new = 0.5 + y_old + y_old / 2 + (x_old % 2) * ((y_old % 2) - 0.5);
                    double x_new = x_old * h;
                    g[vertex_index].x = x_new;
                    g[vertex_index].y = y_new;
                    x_max_vec.emplace_back(x_new);
                    y_max_vec.emplace_back(y_new);
                }
            }

            // Last layer
            for (int x_old = 0; x_old < 2*L+1; ++x_old) {
                int y_old = L;
                int vertex_index = (2 * L + 1) + (y_old - 1) * (2 * L + 2) + x_old;
                double y_new, x_new;
                if (L % 2 == 0) {
                    x_new = (x_old+1) * h;
                    y_new = 0.5 + y_old + y_old / 2 + ((x_old+1) % 2) * ((y_old % 2) - 0.5);
                } else {
                    x_new = x_old * h;
                    y_new = 0.5 + y_old + y_old / 2 + (x_old % 2) * ((y_old % 2) - 0.5);
                }
                g[vertex_index].x = x_new;
                g[vertex_index].y = y_new;
                x_max_vec.emplace_back(x_new);
                y_max_vec.emplace_back(y_new);
            }

            auto x_max = std::max_element(x_max_vec.begin(), x_max_vec.end()); 
            auto y_max = std::max_element(y_max_vec.begin(), y_max_vec.end()); 
            MAX_COORDINATES.emplace_back( *x_max ); // x
            MAX_COORDINATES.emplace_back( *y_max ); // y
        } 
    } else if (lattice_type == "triangular") {
        const double h = std::sqrt(3) / 2.;
        std::vector<double> x_max_vec, y_max_vec;

        for (const auto& v : boost::make_iterator_range(boost::vertices(g))) {
            int x_old = v % L;
            int y_old = v / L;

            double x_new = (y_old % 2) / 2.0 + x_old;
            double y_new = h * y_old;

            g[v].x = x_new;
            g[v].y = y_new;
            x_max_vec.emplace_back(x_new);
            y_max_vec.emplace_back(y_new);
        }

        auto x_max = std::max_element(x_max_vec.begin(), x_max_vec.end()); 
        auto y_max = std::max_element(y_max_vec.begin(), y_max_vec.end()); 
        MAX_COORDINATES.emplace_back( *x_max ); // x
        MAX_COORDINATES.emplace_back( *y_max ); // y

    } else if (lattice_type == "kagome") {
        const double h = std::sqrt(3) / 2.;
        std::vector<double> x_max_vec, y_max_vec;

        for (const auto& v : boost::make_iterator_range(boost::vertices(g))) {
            int x_old = (v) / 3 % L;
            int y_old = (v) / 3 / L;

            double x_new = (y_old % 2) / 2.0 + x_old;
            double y_new = h * y_old;

            double x_sublattice = 0.;
            double y_sublattice = 0.;

            if ((v) % 3 == 0) {
                x_sublattice = 0.;
                y_sublattice = 0.;
            } else if ((v) % 3 == 1) {
                x_sublattice = 0.5;
                y_sublattice = 0.;
            } else {
                x_sublattice = 0.25;
                y_sublattice = 0.5*h;
            }

            g[v].x = x_new + x_sublattice;
            g[v].y = y_new + y_sublattice;
            x_max_vec.emplace_back(x_new + x_sublattice);
            y_max_vec.emplace_back(y_new + y_sublattice);
        }

        auto x_max = std::max_element(x_max_vec.begin(), x_max_vec.end()); 
        auto y_max = std::max_element(y_max_vec.begin(), y_max_vec.end()); 
        MAX_COORDINATES.emplace_back( *x_max ); // x
        MAX_COORDINATES.emplace_back( *y_max ); // y
    }

    // definition of plaquette_vector
    if (lattice_type == "square") {
        for (const auto& v : boost::make_iterator_range(boost::vertices(g))) {
            int x = v % L;
            int y = v / L;
            int pos1 = v;
            int pos2 = y * L + (x + 1) % L;
            int pos3 = ((y+1) % L) * L + (x + 1) % L;
            int pos4 = ((y+1) % L) * L + x;
            if (BOUNDARIES == "periodic") {
                plaquette_vector.emplace_back( 
                    std::vector<std::pair<int,int>> {
                        {pos1,pos2}, {pos2,pos3}, 
                        {pos3,pos4}, {pos4,pos1}
                    }
                );
                plaquette_flip_vector.emplace_back();
                integrated_plaquette_energy_vector.emplace_back(0.);
                plaquette_x_vector.emplace_back(x);
                plaquette_y_vector.emplace_back(y);
            } else {
                if (x < L - 1 && y < L - 1) {
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, 
                            {pos3,pos4}, {pos4,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x);
                    plaquette_y_vector.emplace_back(y);
                }
            }   
        }

    } else if (lattice_type == "cubic") {
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                for (int z = 0; z < L; ++z) {
                    int pos1_1 = z * L * L + y * L + x;
                    int pos1_2 = z * L * L + y * L + (x + 1) % L;
                    int pos1_3 = z * L * L + ((y + 1) % L) * L + (x + 1) % L;
                    int pos1_4 = z * L * L + ((y + 1) % L) * L + x;
                    if (BOUNDARIES == "periodic") {
                        plaquette_vector.emplace_back( 
                            std::vector<std::pair<int,int>> {
                                {pos1_1,pos1_2}, {pos1_2,pos1_3}, 
                                {pos1_3,pos1_4}, {pos1_4,pos1_1}
                            }
                        );
                        plaquette_flip_vector.emplace_back();
                        integrated_plaquette_energy_vector.emplace_back(0.);
                        plaquette_x_vector.emplace_back(x);
                        plaquette_y_vector.emplace_back(y);
                        plaquette_z_vector.emplace_back(z);
                    } else {
                        if (x < L - 1 && y < L - 1) {
                            plaquette_vector.emplace_back( 
                                std::vector<std::pair<int,int>> {
                                    {pos1_1,pos1_2}, {pos1_2,pos1_3}, 
                                    {pos1_3,pos1_4}, {pos1_4,pos1_1}
                                }
                            );
                            plaquette_flip_vector.emplace_back();
                            integrated_plaquette_energy_vector.emplace_back(0.);
                            plaquette_x_vector.emplace_back(x);
                            plaquette_y_vector.emplace_back(y);
                            plaquette_z_vector.emplace_back(z);
                        }
                    }

                    /////////////////////////////////////////////////////////////////////

                    int pos2_1 = z * L * L + y * L + x;
                    int pos2_2 = z * L * L + y * L + (x + 1) % L;
                    int pos2_3 = ((z+1) % L) * L * L + y * L + (x + 1) % L;
                    int pos2_4 = ((z+1) % L) * L * L + y * L + x;
                    if (BOUNDARIES == "periodic") {
                        plaquette_vector.emplace_back( 
                            std::vector<std::pair<int,int>> {
                                {pos2_1,pos2_2}, {pos2_2,pos2_3}, 
                                {pos2_3,pos2_4}, {pos2_4,pos2_1}
                            }
                        );
                        plaquette_flip_vector.emplace_back();
                        integrated_plaquette_energy_vector.emplace_back(0.);
                        plaquette_x_vector.emplace_back(x);
                        plaquette_y_vector.emplace_back(y);
                        plaquette_z_vector.emplace_back(z);
                    } else {
                        if (x < L - 1 && z < L - 1) {
                            plaquette_vector.emplace_back( 
                                std::vector<std::pair<int,int>> {
                                    {pos2_1,pos2_2}, {pos2_2,pos2_3}, 
                                    {pos2_3,pos2_4}, {pos2_4,pos2_1}
                                }
                            );
                            plaquette_flip_vector.emplace_back();
                            integrated_plaquette_energy_vector.emplace_back(0.);
                            plaquette_x_vector.emplace_back(x);
                            plaquette_y_vector.emplace_back(y);
                            plaquette_z_vector.emplace_back(z);
                        }
                    }

                    /////////////////////////////////////////////////////////////////////

                    int pos3_1 = z * L * L + y * L + x;
                    int pos3_2 = z * L * L + ((y + 1) % L) * L + x;
                    int pos3_3 = ((z+1) % L) * L * L + ((y + 1) % L) * L + x;
                    int pos3_4 = ((z+1) % L) * L * L + y * L + x;
                    if (BOUNDARIES == "periodic") {
                        plaquette_vector.emplace_back( 
                            std::vector<std::pair<int,int>> {
                                {pos3_1,pos3_2}, {pos3_2,pos3_3}, 
                                {pos3_3,pos3_4}, {pos3_4,pos3_1}
                            }
                        );
                        plaquette_flip_vector.emplace_back();
                        integrated_plaquette_energy_vector.emplace_back(0.);
                        plaquette_x_vector.emplace_back(x);
                        plaquette_y_vector.emplace_back(y);
                        plaquette_z_vector.emplace_back(z);
                    } else {
                        if (y < L - 1 && z < L - 1) {
                            plaquette_vector.emplace_back( 
                                std::vector<std::pair<int,int>> {
                                    {pos3_1,pos3_2}, {pos3_2,pos3_3}, 
                                    {pos3_3,pos3_4}, {pos3_4,pos3_1}
                                }
                            );
                            plaquette_flip_vector.emplace_back();
                            integrated_plaquette_energy_vector.emplace_back(0.);
                            plaquette_x_vector.emplace_back(x);
                            plaquette_y_vector.emplace_back(y);
                            plaquette_z_vector.emplace_back(z);
                        }
                    }

                    int pos_cube_1 = z * L * L + y * L + x;
                    int pos_cube_2 = z * L * L + ((y + 1) % L) * L + x;
                    int pos_cube_3 = ((z + 1) % L) * L * L + y * L + x;
                    int pos_cube_4 = z * L * L + y * L + ((x + 1) % L);
                    int pos_cube_5 = ((z + 1) % L) * L * L + ((y + 1) % L) * L + x;
                    int pos_cube_6 = ((z + 1) % L) * L * L + y * L + ((x + 1) % L);
                    int pos_cube_7 = z * L * L + ((y + 1) % L) * L + ((x + 1) % L);
                    int pos_cube_8 = ((z + 1) % L) * L * L + ((y + 1) % L) * L + ((x + 1) % L);


                    if (BOUNDARIES == "periodic") { 
                        cube_vector.emplace_back(
                            std::vector<int> {
                                pos_cube_1, pos_cube_2, pos_cube_3, pos_cube_4, 
                                pos_cube_5, pos_cube_6, pos_cube_7, pos_cube_8
                            }
                        );
                        cube_x_vector.emplace_back(x);
                        cube_y_vector.emplace_back(y);
                        cube_z_vector.emplace_back(z);
                    } else {
                        if (y < L - 1 && z < L - 1) {
                            cube_vector.emplace_back(
                                std::vector<int> {
                                    pos_cube_1, pos_cube_2, pos_cube_3, pos_cube_4, 
                                    pos_cube_5, pos_cube_6, pos_cube_7, pos_cube_8
                                }
                            );
                            cube_x_vector.emplace_back(x);
                            cube_y_vector.emplace_back(y);
                            cube_z_vector.emplace_back(z);
                        }
                    } 
                }
            }
        }
        
    } else if (lattice_type == "honeycomb") {
        if (BOUNDARIES == "periodic") {
            // First layer
            for (int i = 1; i < 2*L; i+=2) {
                int vertex_index = i;
                int pos1 = modulo(vertex_index, 2 * L);  // Always on top of plaquette
                int pos2 = modulo(vertex_index + 1, 2 * L);
                int pos3 = modulo(vertex_index + 1, 2 * L) + 2 * L;
                int pos4 = modulo(vertex_index, 2 * L) + 2 * L;
                int pos5 = modulo(vertex_index - 1, 2 * L) + 2 * L;
                int pos6 = modulo(vertex_index - 1, 2 * L);
                plaquette_vector.emplace_back( 
                    std::vector<std::pair<int,int>> {
                        {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                        {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                    }
                );
                plaquette_flip_vector.emplace_back();
                integrated_plaquette_energy_vector.emplace_back(0.);
                plaquette_x_vector.emplace_back(int(i/2.)); 
                plaquette_y_vector.emplace_back(0.);
            }
            
            // Last layers
            for (int y = 1; y < L; ++y) {
                for (int x = (y+1)%2; x < 2*L; x+=2) {
                    int pos1 = y * (2 * L) + modulo(x, 2 * L);  // Always on top of plaquette
                    int pos2 = y * (2 * L) + modulo(x + 1, 2 * L);
                    int pos3 = modulo(y+1, L) * (2 * L) + modulo(x + 1, 2 * L);
                    int pos4 = modulo(y+1, L) * (2 * L) + modulo(x, 2 * L);
                    int pos5 = modulo(y+1, L) * (2 * L) + modulo(x - 1, 2 * L);
                    int pos6 = y * (2 * L) + modulo(x - 1, 2 * L);
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                            {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(int(x/2.));
                    plaquette_y_vector.emplace_back(y);
                }
            }
        } else {
            // First layer
            for (int i = 1; i < 2*L; i+=2) {
                int vertex_index = i;
                int pos1 = vertex_index;  // Always on top of plaquette
                int pos2 = vertex_index + 1;
                int pos3 = vertex_index + 2 + 2 * L;
                int pos4 = vertex_index + 1 + 2 * L;
                int pos5 = vertex_index + 2 * L;
                int pos6 = vertex_index - 1;
                plaquette_vector.emplace_back( 
                    std::vector<std::pair<int,int>> {
                        {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                        {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                    }
                );
                plaquette_flip_vector.emplace_back();
                integrated_plaquette_energy_vector.emplace_back(0.);
                plaquette_x_vector.emplace_back(i);
                plaquette_y_vector.emplace_back(0.);
            }
            
            // Last layers
            for (int y = 1; y < L; ++y) {
                if (y%2 == 0) {
                    for (int x = 0; x < 2*L; ++x) {
                        int vertex_index = (2 * L + 1) + (y - 1) * (2 * L + 2) + x;
                        if (vertex_index % 2 == 1) {
                            continue;
                        } 

                        int pos1, pos2, pos3, pos4, pos5, pos6;

                        if (y == (L - 1)) {
                            pos1 = vertex_index;  // Always on top of plaquette
                            pos2 = vertex_index + 1;
                            if (L % 2 == 1) {
                                pos3 = vertex_index + 3 + 2 * L;
                                pos4 = vertex_index + 2 + 2 * L;
                                pos5 = vertex_index + 1 + 2 * L;
                            } else {
                                pos3 = vertex_index + 2 + 2 * L;
                                pos4 = vertex_index + 1 + 2 * L;
                                pos5 = vertex_index + 2 * L;
                            }
                            pos6 = vertex_index - 1;
                        } else {
                            pos1 = vertex_index;  // Always on top of plaquette
                            pos2 = vertex_index + 1;
                            pos3 = vertex_index + 3 + 2 * L;
                            pos4 = vertex_index + 2 + 2 * L;
                            pos5 = vertex_index + 1 + 2 * L;
                            pos6 = vertex_index - 1;
                        }

                        plaquette_vector.emplace_back( 
                            std::vector<std::pair<int,int>> {
                                {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                                {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                            }
                        );
                        plaquette_flip_vector.emplace_back();
                        integrated_plaquette_energy_vector.emplace_back(0.);
                        plaquette_x_vector.emplace_back(x);
                        plaquette_y_vector.emplace_back(y);
                    }

                } else {
                    for (int x = 1; x < 2*L+1; ++x) {
                        int vertex_index = (2 * L + 1) + (y - 1) * (2 * L + 2) + x;
                        if (vertex_index % 2 == 0) {
                            continue;
                        } 

                        int pos1, pos2, pos3, pos4, pos5, pos6;

                        if (y == (L - 1)) {
                            pos1 = vertex_index;  // Always on top of plaquette
                            pos2 = vertex_index + 1;
                            if (L % 2 == 1) {
                                pos3 = vertex_index + 3 + 2 * L;
                                pos4 = vertex_index + 2 + 2 * L;
                                pos5 = vertex_index + 1 + 2 * L;
                            } else {
                                pos3 = vertex_index + 2 + 2 * L;
                                pos4 = vertex_index + 1 + 2 * L;
                                pos5 = vertex_index + 2 * L;
                            }
                            pos6 = vertex_index - 1;
                        } else {
                            pos1 = vertex_index;  // Always on top of plaquette
                            pos2 = vertex_index + 1;
                            pos3 = vertex_index + 3 + 2 * L;
                            pos4 = vertex_index + 2 + 2 * L;
                            pos5 = vertex_index + 1 + 2 * L;
                            pos6 = vertex_index - 1;
                        }

                        plaquette_vector.emplace_back( 
                            std::vector<std::pair<int,int>> {
                                {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                                {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                            }
                        );
                        plaquette_flip_vector.emplace_back();
                        integrated_plaquette_energy_vector.emplace_back(0.);
                        plaquette_x_vector.emplace_back(x);
                        plaquette_y_vector.emplace_back(y);
                    }
                }
            }
        }
        
    } else if (lattice_type == "triangular") {
        for (const auto& v : boost::make_iterator_range(boost::vertices(g))) {
            int x_old = v % L;
            int y_old = v / L;
            int pos1, pos2, pos3;

            if (BOUNDARIES == "periodic") {
                if (y_old % 2 == 0) {
                    pos1 = y_old * L + x_old;
                    pos2 = y_old * L + modulo(x_old + 1, L);
                    pos3 = modulo(y_old+1, L) * L + x_old;
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x_old);
                    plaquette_y_vector.emplace_back(y_old);
                }

                if (y_old % 2 == 0) {
                    pos1 = y_old * L + x_old;
                    pos2 = y_old * L + modulo(x_old + 1, L);
                    pos3 = modulo(y_old-1, L) * L + x_old;
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x_old);
                    plaquette_y_vector.emplace_back(y_old);
                }

                if (y_old % 2 == 1) {
                    pos1 = y_old * L + x_old;
                    pos2 = y_old * L + modulo(x_old + 1, L);
                    pos3 = modulo(y_old+1, L) * L + modulo(x_old + 1, L);
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x_old);
                    plaquette_y_vector.emplace_back(y_old);
                }

                if (y_old % 2 == 1) {
                    pos1 = y_old * L + x_old;
                    pos2 = y_old * L + modulo(x_old + 1, L);
                    pos3 = modulo(y_old-1, L) * L + modulo(x_old + 1, L);
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x_old);
                    plaquette_y_vector.emplace_back(y_old);
                }
            } else { // open boundaries
                if (y_old % 2 == 0 && (y_old % L) < L - 1 && (x_old % L) < L - 1) {
                    pos1 = v;
                    pos2 = v + 1;
                    pos3 = v + SYSTEM_SIZE; 
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x_old);
                    plaquette_y_vector.emplace_back(y_old);
                }

                if (y_old % 2 == 1 && (y_old % L) < L - 1 && (x_old % L) < L - 1) {
                    pos1 = v;
                    pos2 = v + 1;
                    pos3 = v + SYSTEM_SIZE + 1; 
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x_old);
                    plaquette_y_vector.emplace_back(y_old);
                }

                if (y_old % 2 == 0 && y_old > 0 && (x_old % L) < L - 1) {
                    pos1 = v;
                    pos2 = v + 1;
                    pos3 = v - SYSTEM_SIZE; 
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x_old);
                    plaquette_y_vector.emplace_back(y_old);
                }

                if (y_old % 2 == 1 && y_old > 0 && (x_old % L) < L - 1) {
                    pos1 = v;
                    pos2 = v + 1;
                    pos3 = v - SYSTEM_SIZE + 1; 
                    plaquette_vector.emplace_back( 
                        std::vector<std::pair<int,int>> {
                            {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                        }
                    );
                    plaquette_flip_vector.emplace_back();
                    integrated_plaquette_energy_vector.emplace_back(0.);
                    plaquette_x_vector.emplace_back(x_old);
                    plaquette_y_vector.emplace_back(y_old);
                }   
            }         
        }
    } else if (lattice_type == "kagome") {
        // triangular plaquettes
        for (const auto& v : boost::make_iterator_range(boost::vertices(g))) {
            int pos1, pos2, pos3;
            if ((v) % 3 == 0) {
                pos1 = v;
                pos2 = v+1;
                pos3 = v+2;
                plaquette_vector.emplace_back( 
                    std::vector<std::pair<int,int>> {
                        {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                    }
                );
                plaquette_flip_vector.emplace_back();
                integrated_plaquette_energy_vector.emplace_back(0.);
                plaquette_x_vector.emplace_back(0.);
                plaquette_y_vector.emplace_back(0.);
                // TODO Add plaquette coordinates;
            } else if ((v) % 3 == 1) {
                int x_old = (v) / 3 % L;
                int y_old = (v) / 3 / L;

                int triangular_superlattice_index_down_right = ((y_old+1) % L) * L + modulo(x_old + 1, L);
                int triangular_superlattice_index_down_mid = ((y_old+1) % L) * L + x_old;
                int triangular_superlattice_index_down_left = ((y_old+1) % L) * L + modulo(x_old - 1, L);

                if (BOUNDARIES == "periodic") {
                    if (y_old % 2 == 0) {
                        int v_index_down_mid_triangle_left = 3 * triangular_superlattice_index_down_mid;
                        int v_index_down_left_triangle_right = 3 * triangular_superlattice_index_down_left + 1;
                            pos1 = v+1;
                            pos2 = v_index_down_mid_triangle_left;
                            pos3 = v_index_down_left_triangle_right;
                            plaquette_vector.emplace_back( 
                                std::vector<std::pair<int,int>> {
                                    {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                                }
                            );
                            plaquette_flip_vector.emplace_back();
                            integrated_plaquette_energy_vector.emplace_back(0.);
                            plaquette_x_vector.emplace_back(0.);
                            plaquette_y_vector.emplace_back(0.);
                    } else {
                        int v_index_down_right_triangle_left = 3 * triangular_superlattice_index_down_right;
                        int v_index_down_mid_triangle_right = 3 * triangular_superlattice_index_down_mid + 1;
                            pos1 = v+1;
                            pos2 = v_index_down_right_triangle_left;
                            pos3 = v_index_down_mid_triangle_right;
                            plaquette_vector.emplace_back( 
                                std::vector<std::pair<int,int>> {
                                    {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                                }
                            );
                            plaquette_flip_vector.emplace_back();
                            integrated_plaquette_energy_vector.emplace_back(0.);
                            plaquette_x_vector.emplace_back(0.);
                            plaquette_y_vector.emplace_back(0.);
                    }
                } else {
                    if (y_old < L - 1) {
                        if (y_old % 2 == 0) {
                            int v_index_down_mid_triangle_left = 3 * triangular_superlattice_index_down_mid;
                            int v_index_down_left_triangle_right = 3 * triangular_superlattice_index_down_left + 1;
                            if (x_old > 0) {
                                pos1 = v+1;
                                pos2 = v_index_down_mid_triangle_left;
                                pos3 = v_index_down_left_triangle_right;
                                plaquette_vector.emplace_back( 
                                    std::vector<std::pair<int,int>> {
                                        {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                                    }
                                );
                                plaquette_flip_vector.emplace_back();
                                integrated_plaquette_energy_vector.emplace_back(0.);
                                plaquette_x_vector.emplace_back(0.);
                                plaquette_y_vector.emplace_back(0.);
                            }
                        } else {
                            int v_index_down_right_triangle_left = 3 * triangular_superlattice_index_down_right;
                            int v_index_down_mid_triangle_right = 3 * triangular_superlattice_index_down_mid + 1;
                            if (x_old < L - 1) {
                                pos1 = v+1;
                                pos2 = v_index_down_right_triangle_left;
                                pos3 = v_index_down_mid_triangle_right;
                                plaquette_vector.emplace_back( 
                                    std::vector<std::pair<int,int>> {
                                        {pos1,pos2}, {pos2,pos3}, {pos3,pos1}
                                    }
                                );
                                plaquette_flip_vector.emplace_back();
                                integrated_plaquette_energy_vector.emplace_back(0.);
                                plaquette_x_vector.emplace_back(0.);
                                plaquette_y_vector.emplace_back(0.);
                            }
                        }
                    }
                }
            }
        }

        // honeycomb plaquettes
        /*
        for (const auto& v : boost::make_iterator_range(vs)) { 
            int pos1, pos2, pos3, pos4, pos5, pos6;
            if ((v) % 3 == 1) {
                int x_old = (v) / 3 % L;
                int y_old = (v) / 3 / L;

                int triangular_superlattice_index_down_right = ((y_old+1) % L) * L + modulo(x_old + 1, L);
                int triangular_superlattice_index_down_mid = ((y_old+1) % L) * L + x_old;
                int triangular_superlattice_index_level_right = ((y_old) % L) * L + modulo(x_old + 1, L);

                if (BOUNDARIES == "periodic") {
                    if (y_old % 2 == 0) { 
                        pos1 = v;
                        pos2 = 3 * triangular_superlattice_index_level_right;
                        pos3 = 3 * triangular_superlattice_index_level_right + 2;
                        pos4 = 3 * triangular_superlattice_index_down_mid + 1;
                        pos5 = 3 * triangular_superlattice_index_down_mid;
                        pos6 = v + 1;
                        plaquette_vector.emplace_back( 
                            std::vector<std::pair<int,int>> {
                                {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                                {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                            }
                        );
                        plaquette_flip_vector.emplace_back();
                        integrated_plaquette_energy_vector.emplace_back(0.);
                        plaquette_x_vector.emplace_back(0.);
                        plaquette_y_vector.emplace_back(0.);
                    } else {
                        pos1 = v;
                        pos2 = 3 * triangular_superlattice_index_level_right;
                        pos3 = 3 * triangular_superlattice_index_level_right + 2;
                        pos4 = 3 * triangular_superlattice_index_down_right + 1;
                        pos5 = 3 * triangular_superlattice_index_down_right;
                        pos6 = v + 1;
                        plaquette_vector.emplace_back( 
                            std::vector<std::pair<int,int>> {
                                {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                                {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                            }
                        );
                        plaquette_flip_vector.emplace_back();
                        integrated_plaquette_energy_vector.emplace_back(0.);
                        plaquette_x_vector.emplace_back(0.);
                        plaquette_y_vector.emplace_back(0.);
                    }
                } else {
                    if (y_old < L-1) {
                        if (y_old % 2 == 0) { 
                            if (x_old < L-1) {
                                pos1 = v;
                                pos2 = 3 * triangular_superlattice_index_level_right;
                                pos3 = 3 * triangular_superlattice_index_level_right + 2;
                                pos4 = 3 * triangular_superlattice_index_down_mid + 1;
                                pos5 = 3 * triangular_superlattice_index_down_mid;
                                pos6 = v + 1;
                                plaquette_vector.emplace_back( 
                                    std::vector<std::pair<int,int>> {
                                        {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                                        {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                                    }
                                );
                                plaquette_flip_vector.emplace_back();
                                integrated_plaquette_energy_vector.emplace_back(0.);
                                plaquette_x_vector.emplace_back(0.);
                                plaquette_y_vector.emplace_back(0.);
                            }
                        } else {
                            if (x_old < L-1) {
                                pos1 = v;
                                pos2 = 3 * triangular_superlattice_index_level_right;
                                pos3 = 3 * triangular_superlattice_index_level_right + 2;
                                pos4 = 3 * triangular_superlattice_index_down_right + 1;
                                pos5 = 3 * triangular_superlattice_index_down_right;
                                pos6 = v + 1;
                                plaquette_vector.emplace_back( 
                                    std::vector<std::pair<int,int>> {
                                        {pos1,pos2}, {pos2,pos3}, {pos3,pos4}, 
                                        {pos4,pos5}, {pos5,pos6}, {pos6,pos1}
                                    }
                                );
                                plaquette_flip_vector.emplace_back();
                                integrated_plaquette_energy_vector.emplace_back(0.);
                                plaquette_x_vector.emplace_back(0.);
                                plaquette_y_vector.emplace_back(0.);
                            }
                        }
                    }
                }
            }
        }
        */
    }

    auto plaquette_x_max = std::max_element(plaquette_x_vector.begin(), plaquette_x_vector.end()); 
    auto plaquette_y_max = std::max_element(plaquette_y_vector.begin(), plaquette_y_vector.end()); 
    MAX_PLAQUETTE_COORDINATES.emplace_back( *plaquette_x_max ); // x
    MAX_PLAQUETTE_COORDINATES.emplace_back( *plaquette_y_max ); // y
    if (LATTICE_DIMENSIONALITY == 3) {
        auto plaquette_z_max = std::max_element(plaquette_z_vector.begin(), plaquette_z_vector.end()); 
        MAX_PLAQUETTE_COORDINATES.emplace_back( *plaquette_z_max ); // z
    }

    double max_x = MAX_COORDINATES[0];
    double max_y = MAX_COORDINATES[1];
    int start_y = static_cast<int>(max_y/4.0);
    int end_y = static_cast<int>(3 * max_y/4.0);
    int middle_y = static_cast<int>((start_y + end_y)/2.0);
    int start_x = static_cast<int>(max_x/4.0);
    int end_x = static_cast<int>(3 * max_x/4.0);

    std::tie(half_path_vector, full_path_vector) 
    = construct_fredenhagen_marcu_loops(
        start_y, end_y, middle_y, start_x, end_x, basis
    );

    // This vector stores all elementary plaquettes edges in the lattice
    std::vector<std::vector<Edge>> plaquette_edge_lookup;
    for (size_t p_index = 0; p_index < plaquette_vector.size(); ++p_index) {
        std::vector<Edge> plaquette_edges;
        const auto plaquette_vertex_pairs = get_plaquette_vertex_pairs(p_index);
        for (size_t i = 0; i < plaquette_vertex_pairs.size(); ++i) {
            plaquette_edges.emplace_back(
                edge_in_between( plaquette_vertex_pairs[i].first, plaquette_vertex_pairs[i].second )
            );
        }
        plaquette_edge_lookup.emplace_back(plaquette_edges);
    }

    for (const auto& edg : boost::make_iterator_range(boost::edges(g))) {
        g[edg].spin = default_spin;
    }

    for (size_t p_index = 0; p_index < plaquette_vector.size(); ++p_index) {
        const auto plaquette_vertex_pairs = get_plaquette_vertex_pairs(p_index);
        std::vector<int> plaquette_part_of_cubes;
        for (size_t c_index = 0; c_index < cube_vector.size(); ++c_index) {
            const auto cube_vertices = get_cube_vertices(c_index);
            if (std::all_of(
                plaquette_vertex_pairs.begin(), 
                plaquette_vertex_pairs.end(), 
                [&cube_vertices](const std::pair<int, int> v_pair) {
                    return std::find(
                        cube_vertices.begin(), 
                        cube_vertices.end(), 
                        v_pair.first) != cube_vertices.end() 
                        && std::find(
                            cube_vertices.begin(), 
                            cube_vertices.end(), 
                            v_pair.second
                        ) != cube_vertices.end();
                }
            )) {
                plaquette_part_of_cubes.emplace_back(c_index);
            } 
        }
        plaquette_part_of_cube_lookup.emplace_back(plaquette_part_of_cubes); 
    }

    for (size_t c_index = 0; c_index < cube_vector.size(); ++c_index) {
        std::vector<int> cube_has_plaquettes;
        for (size_t p_index = 0; p_index < plaquette_vector.size(); ++p_index) {
            const auto plaquette_part_of_cubes = plaquette_part_of_cube_lookup[p_index];
            if (std::find(
                plaquette_part_of_cubes.begin(), 
                plaquette_part_of_cubes.end(), 
                c_index
                ) != plaquette_part_of_cubes.end()
            ) {
                cube_has_plaquettes.emplace_back(p_index);
            }
        }
        cube_has_plaquettes_lookup.emplace_back(cube_has_plaquettes);
    }

    edge_dist = std::uniform_int_distribution<int>(0, get_edge_count() - 1);
    vertex_dist = std::uniform_int_distribution<int>(0, get_vertex_count() - 1);
    plaquette_dist = std::uniform_int_distribution<int>(0, plaquette_vector.size() - 1);

    for (auto e : boost::make_iterator_range(boost::edges(g))) {
        edge_vector.emplace_back(boost::source(e, g), boost::target(e, g));
    }

    return g;
}

std::pair<std::vector<Lattice::VertexPair>, std::vector<Lattice::VertexPair>> 
Lattice::construct_fredenhagen_marcu_loops(
    int start_y, int end_y, int middle_y, int start_x, int end_x, char basis
) {
    std::vector<VertexPair> half_loop;
    std::vector<VertexPair> full_loop;

    if (basis == 'z') {
        if (LATTICE_TYPE == "square" || LATTICE_TYPE == "cubic" || LATTICE_TYPE == "triangular") {
            if (LATTICE_DIMENSIONALITY == 2) {
                int prev_vertex = middle_y * SYSTEM_SIZE + start_x;
                int next_vertex;

                for (int y = middle_y+1; y < end_y+1; ++y) {
                    next_vertex = y * SYSTEM_SIZE + start_x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 

                for (int x = start_x+1; x < end_x+1; ++x) {
                    next_vertex = end_y * SYSTEM_SIZE + x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                for (int y = end_y-1; y > middle_y-1; --y) {
                    next_vertex = y * SYSTEM_SIZE + end_x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 

                // Store half_loop after half of the path
                half_loop = full_loop;

                for (int y = middle_y-1; y > start_y-1; --y) {
                    next_vertex = y * SYSTEM_SIZE + end_x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 

                for (int x = end_x-1; x > start_x-1; --x) {
                    next_vertex = start_y * SYSTEM_SIZE + x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                for (int y = start_y+1; y < middle_y+1; ++y) {
                    next_vertex = y * SYSTEM_SIZE + start_x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 
            } else if (LATTICE_DIMENSIONALITY == 3) {
                double max_z = MAX_COORDINATES[2];
                int middle_z = static_cast<int>(max_z/2.0);

                int prev_vertex = middle_z * SYSTEM_SIZE * SYSTEM_SIZE + middle_y * SYSTEM_SIZE + start_x;
                int next_vertex;

                for (int y = middle_y+1; y < end_y+1; ++y) {
                    next_vertex = middle_z * SYSTEM_SIZE * SYSTEM_SIZE + y * SYSTEM_SIZE + start_x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 

                for (int x = start_x+1; x < end_x+1; ++x) {
                    next_vertex = middle_z * SYSTEM_SIZE * SYSTEM_SIZE + end_y * SYSTEM_SIZE + x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                for (int y = end_y-1; y > middle_y-1; --y) {
                    next_vertex = middle_z * SYSTEM_SIZE * SYSTEM_SIZE + y * SYSTEM_SIZE + end_x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 

                // Store half_loop after half of the path
                half_loop = full_loop;

                for (int y = middle_y-1; y > start_y-1; --y) {
                    next_vertex = middle_z * SYSTEM_SIZE * SYSTEM_SIZE + y * SYSTEM_SIZE + end_x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 

                for (int x = end_x-1; x > start_x-1; --x) {
                    next_vertex = middle_z * SYSTEM_SIZE * SYSTEM_SIZE + start_y * SYSTEM_SIZE + x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                }

                for (int y = start_y+1; y < middle_y+1; ++y) {
                    next_vertex = middle_z * SYSTEM_SIZE * SYSTEM_SIZE + y * SYSTEM_SIZE + start_x;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 
            }
        } else if (LATTICE_TYPE == "honeycomb") {
            if (BOUNDARIES == "periodic") {
                int prev_vertex, next_vertex;

                prev_vertex = (2 * SYSTEM_SIZE) * middle_y + start_x + 1;
                if (middle_y%2 == 0 && prev_vertex%2 == 1) {
                    prev_vertex = (2 * SYSTEM_SIZE) * middle_y + start_x;
                } else if (middle_y%2 == 1 && prev_vertex%2 == 0) {
                    prev_vertex = (2 * SYSTEM_SIZE) * middle_y + start_x;
                }

                for (int y = middle_y; y < end_y; y+=2) {
                    next_vertex = prev_vertex + 2 * SYSTEM_SIZE;
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
                    next_vertex = prev_vertex - 2 * SYSTEM_SIZE;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;

                    next_vertex = prev_vertex + 1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 

                // Store half_loop after half of the path
                half_loop = full_loop;

                for (int y = middle_y; y > start_y; y-=2) {
                    next_vertex = prev_vertex - 2 * SYSTEM_SIZE;
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
                    next_vertex = prev_vertex + 2 * SYSTEM_SIZE;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;

                    next_vertex = prev_vertex - 1;
                    full_loop.emplace_back(prev_vertex, next_vertex);
                    prev_vertex = next_vertex;
                } 
            } else {
#ifndef NDEBUG
                throw std::invalid_argument("Open honeycomb lattice is not supported for Wilson loops.");
#endif
                return std::make_pair(half_loop, full_loop);
            }

        } else if (LATTICE_TYPE == "kagome") {
            // TODO implement this
#ifndef NDEBUG
            throw std::invalid_argument("Kagome lattice is not supported for Wilson loops.");
#endif
            return std::make_pair(half_loop, full_loop);
        } 
    } else if (basis == 'x') {
        if (SYSTEM_SIZE < 6) {
#ifndef NDEBUG
            throw std::runtime_error(
                std::format(
                    "System size for 't Hooft loops in x-basis has to be at least L=6 but is L={}.", 
                    SYSTEM_SIZE
                )
            );
#endif
            return std::make_pair(half_loop, full_loop);
        }

        if (LATTICE_TYPE == "square") {
            int start_vertex = 0;

            for (int y = middle_y; y < end_y+1; ++y) {
                start_vertex = y * SYSTEM_SIZE + start_x;
                full_loop.emplace_back(start_vertex, start_vertex-1);
            } 

            for (int x = start_x; x < end_x+1; ++x) {
                start_vertex = end_y * SYSTEM_SIZE + x;
                full_loop.emplace_back(start_vertex, start_vertex+SYSTEM_SIZE);
            }

            for (int y = end_y; y > middle_y-1; --y) {
                start_vertex = y * SYSTEM_SIZE + end_x;
                full_loop.emplace_back(start_vertex, start_vertex+1);
            } 

            // Store half_loop after half of the path
            half_loop = full_loop;

            for (int y = middle_y-1; y > start_y-1; --y) {
                start_vertex = y * SYSTEM_SIZE + end_x;
                full_loop.emplace_back(start_vertex, start_vertex+1);
            } 

            for (int x = end_x; x > start_x-1; --x) {
                start_vertex = start_y * SYSTEM_SIZE + x;
                full_loop.emplace_back(start_vertex, start_vertex-SYSTEM_SIZE);
            }

            for (int y = start_y; y < middle_y+1; ++y) {
                start_vertex = y * SYSTEM_SIZE + start_x;
                full_loop.emplace_back(start_vertex, start_vertex-1);
            } 
        } else if (LATTICE_TYPE == "triangular") {
            int start_vertex = 0;
            int end_y_tri = end_y;
            int start_y_tri = start_y;

            if (end_y % 2 == 1) end_y_tri = end_y + 1;
            if (start_y % 2 == 0) start_y_tri = start_y + 1;

            for (int y = middle_y; y < end_y_tri+1; ++y) {
                if (y%2 == 0) {
                    start_vertex = y * SYSTEM_SIZE + start_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex+1, start_vertex+SYSTEM_SIZE);
                } else {
                    start_vertex = y * SYSTEM_SIZE + start_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex, start_vertex+SYSTEM_SIZE+1);
                }
            } 

            for (int x = start_x+1; x < end_x; ++x) {
                start_vertex = end_y_tri * SYSTEM_SIZE + x;
                full_loop.emplace_back(start_vertex, start_vertex+SYSTEM_SIZE);

                full_loop.emplace_back(start_vertex+1, start_vertex+SYSTEM_SIZE);
            }

            // Last bond
            start_vertex = end_y_tri * SYSTEM_SIZE + end_x;
            full_loop.emplace_back(start_vertex, start_vertex+SYSTEM_SIZE);

            for (int y = end_y_tri; y > middle_y-1; --y) {
                if (y%2 == 0) {
                    start_vertex = y * SYSTEM_SIZE + end_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex+1, start_vertex-SYSTEM_SIZE);
                } else {
                    start_vertex = y * SYSTEM_SIZE + end_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex, start_vertex-SYSTEM_SIZE+1);
                }
            } 

            // Store half_loop after half of the path
            half_loop = full_loop;

            for (int y = middle_y-1; y > start_y_tri-1; --y) {
                if (y%2 == 0) {
                    start_vertex = y * SYSTEM_SIZE + end_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex+1, start_vertex-SYSTEM_SIZE);
                } else {
                    start_vertex = y * SYSTEM_SIZE + end_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex, start_vertex-SYSTEM_SIZE+1);
                }
            } 

            for (int x = end_x; x > start_x+1; --x) {
                start_vertex = start_y_tri * SYSTEM_SIZE + x;
                full_loop.emplace_back(start_vertex, start_vertex-SYSTEM_SIZE);

                full_loop.emplace_back(start_vertex-1, start_vertex-SYSTEM_SIZE);
            }

            // Last bond
            start_vertex = start_y_tri * SYSTEM_SIZE + start_x + 1;
            full_loop.emplace_back(start_vertex, start_vertex-SYSTEM_SIZE);

            for (int y = start_y_tri; y < middle_y; ++y) {
                if (y%2 == 0) {
                    start_vertex = y * SYSTEM_SIZE + start_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex+1, start_vertex+SYSTEM_SIZE);
                } else {
                    start_vertex = y * SYSTEM_SIZE + start_x;
                    full_loop.emplace_back(start_vertex, start_vertex+1);

                    full_loop.emplace_back(start_vertex, start_vertex+SYSTEM_SIZE+1);
                }
            } 

        } else if (LATTICE_TYPE == "honeycomb") {
            // TODO implement this
#ifndef NDEBUG
            throw std::invalid_argument("Honeycomb lattice is not supported for 't Hooft loops.");
#endif
            return std::make_pair(half_loop, full_loop);
        } else if (LATTICE_TYPE == "kagome") {
            // TODO implement this
#ifndef NDEBUG
            throw std::invalid_argument("Kagome lattice is not supported for 't Hooft loops.");
#endif
            return std::make_pair(half_loop, full_loop);
        } 
    } 

    return std::make_pair(half_loop, full_loop);  
}

// TODO
void Lattice::build_caches_() {
    egde_cache_.clear();
    for (auto e : boost::make_iterator_range(boost::edges(g))) {
        g[e].part_of_plaquette_lookup.clear();
        egde_cache_.emplace_back(e);
    }

    // ----- 1) Plaquette -> edges (arbitrary length) -----
    plaquette_edges_cache_.clear();
    plaquette_edges_cache_.resize(plaquette_vector.size());

    for (size_t p = 0; p < plaquette_vector.size(); ++p) {
        const auto vpairs = get_plaquette_vertex_pairs(static_cast<int>(p));
        auto& pedges = plaquette_edges_cache_[p];
        pedges.clear();
        pedges.reserve(vpairs.size());

        for (const auto& pr : vpairs) {
            auto uv = boost::edge(pr.first, pr.second, g);
            if (!uv.second) {
                // try opposite orientation for directed graphs
                uv = boost::edge(pr.second, pr.first, g);
            }
            if (!uv.second) {
                // hard fail in debug; avoid silent UB on invalid descriptor
                throw std::runtime_error("build_caches_: missing edge between plaquette vertices");
            }
            pedges.emplace_back(uv.first);
            g[uv.first].part_of_plaquette_lookup.emplace_back(static_cast<int>(p));
        }
    }

    // ----- 2) Star (vertex) -> incident edges (use vertex_index map!) -----
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
    const auto& spin_flips = g[edg].spin_flips;
    const auto it = std::lower_bound(spin_flips.begin(), spin_flips.end(), tau);

    int index = 0;
    if (it != spin_flips.end() && *it == tau) [[likely]] {
        index = static_cast<int>(it - spin_flips.begin());
    } else [[unlikely]] {
        throw std::runtime_error(std::format("get_spin_flip_index: There is no spin flip at {}.", tau));
    }
    return index;
}

int Lattice::get_single_spin_flip_index(const Edge& edg, double tau) {
    const auto& single_spin_flips = g[edg].single_spin_flips;
    const auto it = std::lower_bound(single_spin_flips.begin(), single_spin_flips.end(), tau);

    int index = 0;
    if (it != single_spin_flips.end() && *it == tau) [[likely]] {
        index = static_cast<int>(it - single_spin_flips.begin());
    } else [[unlikely]] {
        throw std::runtime_error(std::format("get_single_spin_flip_index: There is no single spin flip at {}.", tau));
    }
    return index;
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
    if (imag_time_next_single_spin_flip < imag_time_single_spin_flip) [[unlikely]] {
        throw std::invalid_argument("delete_double_single_spin_flip: imag_time_next_tuple_flip has to be larger than imag_time_tuple_flip.");
    }

    auto& spin_flips = g[edg].spin_flips;
    auto it_2 = std::lower_bound(spin_flips.begin(), spin_flips.end(), imag_time_next_single_spin_flip);

    if (it_2 != spin_flips.end() && *it_2 == imag_time_next_single_spin_flip) [[likely]] {
        it_2 = spin_flips.erase(it_2);
    } else [[unlikely]] {
        throw std::runtime_error(std::format("delete_double_single_spin_flip: There is no spin flip at {}.", imag_time_next_single_spin_flip));
    }

    auto r_it = std::find(std::make_reverse_iterator(it_2), spin_flips.rend(), imag_time_single_spin_flip);
    if (r_it == spin_flips.rend()) [[unlikely]] {
        throw std::runtime_error(std::format("delete_double_single_spin_flip: There is no spin flip at {}.", imag_time_single_spin_flip));
    }
    auto it_1 = std::prev(r_it.base());
    spin_flips.erase(it_1);


    auto& single_spin_flips = g[edg].single_spin_flips;
    auto it_s_2 = std::lower_bound(single_spin_flips.begin(), single_spin_flips.end(), imag_time_next_single_spin_flip);

    if (it_s_2 != single_spin_flips.end() && *it_s_2 == imag_time_next_single_spin_flip) [[likely]] {
        it_s_2 = single_spin_flips.erase(it_s_2);
    } else [[unlikely]] {
        throw std::runtime_error(std::format("delete_double_single_spin_flip: There is no spin flip at {}.", imag_time_next_single_spin_flip));
    }

    auto r_s_it = std::find(std::make_reverse_iterator(it_s_2), single_spin_flips.rend(), imag_time_single_spin_flip);
    if (r_s_it == single_spin_flips.rend()) [[unlikely]] {
        throw std::runtime_error(std::format("delete_double_single_spin_flip: There is no spin flip at {}.", imag_time_single_spin_flip));
    }
    auto it_s_1 = std::prev(r_s_it.base());
    single_spin_flips.erase(it_s_1);
}

void Lattice::delete_single_spin_flip(const Edge& edg, int spin_flip_index) {
    auto& spin_flips = g[edg].spin_flips;
    const auto imag_time = spin_flips[spin_flip_index];
    spin_flips.erase(spin_flips.begin() + spin_flip_index);    

    auto& single_spin_flips = g[edg].single_spin_flips;
    auto it = std::lower_bound(single_spin_flips.begin(), single_spin_flips.end(), imag_time);

    if (it != single_spin_flips.end() && *it == imag_time) [[likely]] {
        it = single_spin_flips.erase(it);
    } else [[unlikely]] {
        throw std::runtime_error(std::format("delete_single_spin_flip: There is no single spin flip at {}.", imag_time));
    }
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
        auto& spin_flips = g[edg].spin_flips;
        auto it_2 = std::lower_bound(spin_flips.begin(), spin_flips.end(), imag_time_next_tuple_flip);

        if (it_2 != spin_flips.end() && *it_2 == imag_time_next_tuple_flip) [[likely]] {
            it_2 = spin_flips.erase(it_2);
        } else [[unlikely]] {
            throw std::runtime_error(std::format("delete_double_tuple_flip: There is no spin flip at {}.", imag_time_next_tuple_flip));
        }
        
        auto r_it = std::find(std::make_reverse_iterator(it_2), spin_flips.rend(), imag_time_tuple_flip);
        if (r_it == spin_flips.rend()) [[unlikely]] {
            throw std::runtime_error(std::format("delete_double_tuple_flip: There is no spin flip at {}.", imag_time_tuple_flip));
        }
        auto it_1 = std::prev(r_it.base());
        spin_flips.erase(it_1);
    }

    if (BASIS == 'x') {
        auto& tuple_spin_flips = plaquette_flip_vector[tuple_index];
        auto it_2 = std::lower_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), imag_time_next_tuple_flip);

        if (it_2 != tuple_spin_flips.end() && *it_2 == imag_time_next_tuple_flip) [[likely]] {
            it_2 = tuple_spin_flips.erase(it_2);
        } else [[unlikely]] {
            throw std::runtime_error(std::format("delete_double_tuple_flip: There is no spin flip at {}.", imag_time_next_tuple_flip));
        }
        
        auto r_it = std::find(std::make_reverse_iterator(it_2), tuple_spin_flips.rend(), imag_time_tuple_flip);
        if (r_it == tuple_spin_flips.rend()) [[unlikely]] {
            throw std::runtime_error(std::format("delete_double_tuple_flip: There is no spin flip at {}.", imag_time_tuple_flip));
        }
        auto it_1 = std::prev(r_it.base());
        tuple_spin_flips.erase(it_1);
    } else {
        auto& tuple_spin_flips = g[tuple_index].star_flips;
        auto it_2 = std::lower_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), imag_time_next_tuple_flip);

        if (it_2 != tuple_spin_flips.end() && *it_2 == imag_time_next_tuple_flip) [[likely]] {
            it_2 = tuple_spin_flips.erase(it_2);
        } else [[unlikely]] {
            throw std::runtime_error(std::format("delete_double_tuple_flip: There is no spin flip at {}.", imag_time_next_tuple_flip));
        }
        
        auto r_it = std::find(std::make_reverse_iterator(it_2), tuple_spin_flips.rend(), imag_time_tuple_flip);
        if (r_it == tuple_spin_flips.rend()) [[unlikely]] {
            throw std::runtime_error(std::format("delete_double_tuple_flip: There is no spin flip at {}.", imag_time_tuple_flip));
        }
        auto it_1 = std::prev(r_it.base());
        tuple_spin_flips.erase(it_1);
    }
}

void Lattice::delete_tuple_flip(int tuple_index, std::span<const Edge> tuple_edges, double imag_time_tuple_flip) {
    for (const Edge& edg : tuple_edges) {
        auto& spin_flips = g[edg].spin_flips;
        const auto it = std::lower_bound(spin_flips.begin(), spin_flips.end(), imag_time_tuple_flip);

        if (it != spin_flips.end() && *it == imag_time_tuple_flip) [[likely]] {
            spin_flips.erase(it);
        } else [[unlikely]] {
            throw std::runtime_error(std::format("delete_tuple_flip: There is no spin flip at {}.", imag_time_tuple_flip));
        }
    }

    if (BASIS == 'x') {
        auto& tuple_spin_flips = plaquette_flip_vector[tuple_index];
        const auto it = std::lower_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), imag_time_tuple_flip);

        if (it != tuple_spin_flips.end() && *it == imag_time_tuple_flip) [[likely]] {
            tuple_spin_flips.erase(it);
        } else [[unlikely]] {
            throw std::runtime_error(std::format("delete_tuple_flip: There is no spin flip at {}.", imag_time_tuple_flip));
        }
    } else {
        auto& tuple_spin_flips = g[tuple_index].star_flips;
        const auto it = std::lower_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), imag_time_tuple_flip);

        if (it != tuple_spin_flips.end() && *it == imag_time_tuple_flip) [[likely]] {
            tuple_spin_flips.erase(it);
        } else [[unlikely]] {
            throw std::runtime_error(std::format("delete_tuple_flip: There is no spin flip at {}.", imag_time_tuple_flip));
        }
    }
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

double Lattice::flip_next_imag_time(const Edge& edg, double tau) { // TODO replace by index based method
    auto& spin_flips = g[edg].spin_flips;
    if (spin_flips.empty()) [[unlikely]] {
        return tau;
    }

    auto it = std::upper_bound(spin_flips.begin(), spin_flips.end(), tau);

    if (it != spin_flips.end()) [[likely]] {
        return *it;
    } else [[unlikely]] {
        return spin_flips.front();
    }
}

std::vector<double> Lattice::flip_next_imag_times_tuple(std::span<const Edge> tuple_edges, double tau) {
    std::vector<double> imag_times;
    for (const Edge& edg : tuple_edges) {
        imag_times.emplace_back(flip_next_imag_time(edg, tau));
    }
    return imag_times;
}

double Lattice::flip_prev_imag_time(const Edge& edg, double tau) {
    auto& spin_flips = g[edg].spin_flips;
    if (spin_flips.empty()) [[unlikely]] {
        return tau;
    }

    auto it = std::lower_bound(spin_flips.begin(), spin_flips.end(), tau);

    if (it != spin_flips.begin()) [[likely]] {
        return *(std::prev(it));
    } else [[unlikely]] {
        return spin_flips.back();
    }
}

std::vector<double> Lattice::flip_prev_imag_times_tuple(std::span<const Edge> tuple_edges, double tau) {
    std::vector<double> imag_times;
    for (const Edge& edg : tuple_edges) {
        imag_times.emplace_back(flip_prev_imag_time(edg, tau));
    }
    return imag_times;
}

void Lattice::insert_double_single_spin_flip(const Edge& edg, double tau_left, double tau_right) {
    auto& spin_flips = g[edg].spin_flips;
    if (spin_flips.empty()) [[unlikely]] {
        spin_flips.push_back(tau_left);
        spin_flips.push_back(tau_right);
    } else [[likely]] {
        auto it_right = std::upper_bound(spin_flips.begin(), spin_flips.end(), tau_right);
        it_right = spin_flips.insert(it_right, tau_right);

        // Find insertion point for tau_left in [spin_flips.begin(), it_right)
        auto it_left = std::lower_bound(spin_flips.begin(), it_right, tau_left);
        spin_flips.insert(it_left, tau_left);
    }

    auto& single_spin_flips = g[edg].single_spin_flips;
    if (single_spin_flips.empty()) [[unlikely]] {
        single_spin_flips.push_back(tau_left);
        single_spin_flips.push_back(tau_right);
    } else [[likely]] {
        auto it_right = std::upper_bound(single_spin_flips.begin(), single_spin_flips.end(), tau_right);
        it_right = single_spin_flips.insert(it_right, tau_right);

        auto it_left = std::lower_bound(single_spin_flips.begin(), it_right, tau_left);
        single_spin_flips.insert(it_left, tau_left);
    }
}

void Lattice::insert_single_spin_flip(const Edge& edg, double tau) {
    auto& spin_flips = g[edg].spin_flips;
    if (spin_flips.empty()) [[unlikely]] {
        spin_flips.emplace_back(tau);
    } else [[likely]] {
        auto right = std::upper_bound(spin_flips.begin(), spin_flips.end(), tau);
        spin_flips.insert(right, tau);
    }

    auto& single_spin_flips = g[edg].single_spin_flips;
    if (single_spin_flips.empty()) [[unlikely]] {
        single_spin_flips.emplace_back(tau);
    } else [[likely]] {
        auto right = std::upper_bound(single_spin_flips.begin(), single_spin_flips.end(), tau);
        single_spin_flips.insert(right, tau);
    }
}

void Lattice::insert_double_tuple_flip(
    int tuple_index, 
    std::span<const Edge> tuple_edges, 
    double tau_left, 
    double tau_right) {
    // Is this really necessary? Should check...
    if (tau_right < tau_left) [[unlikely]] {
        throw std::invalid_argument("insert_double_tuple_flip: tau_right has to be larger than tau_left.");
    }

    for (const Edge& edg : tuple_edges) {
        auto& spin_flips = g[edg].spin_flips;
        if (spin_flips.empty()) [[unlikely]] {
            spin_flips.emplace_back(tau_left);
            spin_flips.emplace_back(tau_right);
        } else [[likely]] {
            auto it_right = std::upper_bound(spin_flips.begin(), spin_flips.end(), tau_right);
            it_right = spin_flips.insert(it_right, tau_right);

            auto it_left = std::lower_bound(spin_flips.begin(), it_right, tau_left);
            spin_flips.insert(it_left, tau_left);
        }
    }

    if (BASIS == 'x') {
        auto& tuple_spin_flips = plaquette_flip_vector[tuple_index];
        if (tuple_spin_flips.empty()) [[unlikely]] {
            tuple_spin_flips.emplace_back(tau_left);
            tuple_spin_flips.emplace_back(tau_right);
        } else [[likely]] {
            auto it_right = std::upper_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), tau_right);
            it_right = tuple_spin_flips.insert(it_right, tau_right);

            auto it_left = std::lower_bound(tuple_spin_flips.begin(), it_right, tau_left);
            tuple_spin_flips.insert(it_left, tau_left);
        }
    } else {
        auto& tuple_spin_flips = g[tuple_index].star_flips;
        if (tuple_spin_flips.empty()) [[unlikely]] {
            tuple_spin_flips.emplace_back(tau_left);
            tuple_spin_flips.emplace_back(tau_right);
        } else [[likely]] {
            auto it_right = std::upper_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), tau_right);
            it_right = tuple_spin_flips.insert(it_right, tau_right);

            auto it_left = std::lower_bound(tuple_spin_flips.begin(), it_right, tau_left);
            tuple_spin_flips.insert(it_left, tau_left);
        }
    }
}

void Lattice::insert_tuple_flip(int tuple_index, std::span<const Edge> tuple_edges, double tau) {
    for (const Edge& edg : tuple_edges) {
        auto& spin_flips = g[edg].spin_flips;
        if (spin_flips.empty()) [[unlikely]] {
            spin_flips.emplace_back(tau);
        } else [[likely]] {
            auto right = std::upper_bound(spin_flips.begin(), spin_flips.end(), tau);
            spin_flips.insert(right, tau);
        }
    }

    if (BASIS == 'x') {
        auto& tuple_spin_flips = plaquette_flip_vector[tuple_index];
        if (tuple_spin_flips.empty()) [[unlikely]] {
            tuple_spin_flips.emplace_back(tau);
        } else [[likely]] {
            auto right = std::upper_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), tau);
            tuple_spin_flips.insert(right, tau);
        }
    } else {
        auto& tuple_spin_flips = g[tuple_index].star_flips;
        if (tuple_spin_flips.empty()) [[unlikely]] {
            tuple_spin_flips.emplace_back(tau);
        } else [[likely]] {
            auto right = std::upper_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), tau);
            tuple_spin_flips.insert(right, tau);
        }
    }
}

void Lattice::move_spin_flip(
    const Edge& edg, 
    int spin_flip_index, 
    double tau_new, 
    bool no_move_over_beta
) {
    if (no_move_over_beta) [[likely]] {
        int single_spin_flip_index = get_single_spin_flip_index(edg, g[edg].spin_flips[spin_flip_index]); 
        set_single_spin_flip_imag_time(edg, single_spin_flip_index, tau_new);
        set_spin_flip_imag_time(edg, spin_flip_index, tau_new);
    } else [[unlikely]] {
        auto& spin_flips = g[edg].spin_flips;
        auto& single_spin_flips = g[edg].single_spin_flips;
        if (spin_flip_index == 0) { // rotate to the left 
            *spin_flips.begin() = tau_new;
            std::rotate(spin_flips.begin(), spin_flips.begin() + 1, spin_flips.end());

            *single_spin_flips.begin() = tau_new;
            std::rotate(single_spin_flips.begin(), single_spin_flips.begin() + 1, single_spin_flips.end());
        } else { // rotate to the right 
            std::rotate(spin_flips.rbegin(), spin_flips.rbegin() + 1, spin_flips.rend());
            *spin_flips.begin() = tau_new;

            std::rotate(single_spin_flips.rbegin(), single_spin_flips.rbegin() + 1, single_spin_flips.rend());
            *single_spin_flips.begin() = tau_new;
        }
    } 
}

void Lattice::move_tuple_flip(
    int tuple_index, 
    std::span<const Edge> tuple_edges, 
    double tau_old, 
    double tau_new, 
    bool no_move_over_beta
) {
    for (const Edge& edg : tuple_edges) {
        const auto& spin_flips = g[edg].spin_flips;
        const auto it = std::lower_bound(spin_flips.begin(), spin_flips.end(), tau_old);
        int index = -1;
        if (it != spin_flips.end()) [[likely]] {
            index = it - spin_flips.begin(); // TODO does not need the index, only whether its 0 or not
        } 
        if (no_move_over_beta) [[likely]] {
            set_spin_flip_imag_time(edg, index, tau_new);
        } else [[unlikely]] {
            auto& spin_flips = g[edg].spin_flips;
            if (index == 0) { // rotate to the left 
                *spin_flips.begin() = tau_new;
                std::rotate(spin_flips.begin(), spin_flips.begin() + 1, spin_flips.end());
            } else { // rotate to the right 
                std::rotate(spin_flips.rbegin(), spin_flips.rbegin() + 1, spin_flips.rend());
                *spin_flips.begin() = tau_new;
            }
        } 
    }

    if (BASIS == 'x') {
        auto& tuple_spin_flips = plaquette_flip_vector[tuple_index];
        const auto it_tuple = std::lower_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), tau_old);
        int index = -1;
        if (it_tuple != tuple_spin_flips.end()) [[likely]] {
            index = it_tuple - tuple_spin_flips.begin(); // TODO does not need the index, only whether its 0 or not
        } 
        if (no_move_over_beta) [[likely]] {
            tuple_spin_flips[index] = tau_new;
        } else [[unlikely]] {
            if (index == 0) { // rotate to the left 
                *tuple_spin_flips.begin() = tau_new;
                std::rotate(tuple_spin_flips.begin(), tuple_spin_flips.begin() + 1, tuple_spin_flips.end());
            } else { // rotate to the right 
                std::rotate(tuple_spin_flips.rbegin(), tuple_spin_flips.rbegin() + 1, tuple_spin_flips.rend());
                *tuple_spin_flips.begin() = tau_new;
            }
        } 
    } else {
        auto& tuple_spin_flips = g[tuple_index].star_flips;
        const auto it_tuple = std::lower_bound(tuple_spin_flips.begin(), tuple_spin_flips.end(), tau_old);
        int index = -1;
        if (it_tuple != tuple_spin_flips.end()) [[likely]] {
            index = it_tuple - tuple_spin_flips.begin(); // TODO does not need the index, only whether its 0 or not
        } 
        if (no_move_over_beta) [[likely]] {
            tuple_spin_flips[index] = tau_new;
        } else [[unlikely]] {
            if (index == 0) { // rotate to the left 
                *tuple_spin_flips.begin() = tau_new;
                std::rotate(tuple_spin_flips.begin(), tuple_spin_flips.begin() + 1, tuple_spin_flips.end());
            } else { // rotate to the right 
                std::rotate(tuple_spin_flips.rbegin(), tuple_spin_flips.rbegin() + 1, tuple_spin_flips.rend());
                *tuple_spin_flips.begin() = tau_new;
            }
        } 
    }
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
    size_t edg_index = edge_dist(*rng);
    auto [u, v] = edge_vector[edg_index];
    auto edg = egde_cache_[edg_index];
    return {edg, u, v};
}

int Lattice::get_random_vertex() {
    return vertex_dist(*rng);
}

int Lattice::get_random_plaquette_index() {
    return plaquette_dist(*rng);
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

std::tuple<double, std::vector<int>, std::vector<double>> 
Lattice::integrated_star_energy_diff(
    const Edge& edg, 
    double imag_time_1, 
    double imag_time_2, 
    bool total_cache
) {
    if (imag_time_1 == imag_time_2) [[unlikely]] {
        throw std::invalid_argument("integrated_star_energy_diff: Time interval has to be non-zero.");
    } 

    std::vector<int> star_centers;
    std::vector<double> star_potential_energy_diffs;

    auto [source_v, target_v] = vertices_of_edge(edg);
    double energy = 0.;
    double source_energy = 0.;
    double target_energy = 0.;
    if (total_cache) {
        source_energy = -2 * get_potential_star_energy(source_v);
        target_energy = -2 * get_potential_star_energy(target_v);
    } else {
        const auto& sedges_source = get_star_edges(source_v);
        source_energy = integrated_tuple_energy_diff(sedges_source, imag_time_1, imag_time_2);
        const auto& sedges_target = get_star_edges(target_v);
        target_energy = integrated_tuple_energy_diff(sedges_target, imag_time_1, imag_time_2);
    }
    star_centers.emplace_back(source_v);
    star_potential_energy_diffs.emplace_back(source_energy);
    energy += source_energy;
    star_centers.emplace_back(target_v);
    star_potential_energy_diffs.emplace_back(target_energy);
    energy += target_energy;

    return {energy, std::move(star_centers), std::move(star_potential_energy_diffs)};
}

std::tuple<double, std::vector<int>, std::vector<double>> 
Lattice::integrated_plaquette_energy_diff(
    const Edge& edg, 
    double imag_time_1, 
    double imag_time_2, 
    bool total_cache
) {
    if (imag_time_1 == imag_time_2) [[unlikely]] {
        throw std::invalid_argument("integrated_plaquette_energy_diff: Time interval has to be non-zero.");
    } 
    std::vector<int> plaquette_indices;
    std::vector<double> plaquette_potential_energy_diffs;

    double energy = 0.;
    for (int p_index : g[edg].part_of_plaquette_lookup ) {
        double energy_p = 0.;
        if (total_cache) {
            energy_p = -2*get_potential_plaquette_energy(p_index);
        } else {
            const auto& pedges = get_plaquette_edges(p_index);
            energy_p = integrated_tuple_energy_diff(pedges, imag_time_1, imag_time_2);
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
        total_integrated_star_energy += integrated_tuple_energy(sedges, 0, BETA);
    }
    return total_integrated_star_energy;
}

double Lattice::total_integrated_plaquette_energy() {
    double total_integrated_plaquette_energy = 0.;
    for (size_t plaquette_index = 0; plaquette_index < plaquette_vector.size(); ++plaquette_index) {
        const auto& pedges = get_plaquette_edges(plaquette_index);
        total_integrated_plaquette_energy += integrated_tuple_energy(pedges, 0, BETA);
    }
    return total_integrated_plaquette_energy;
}

namespace {
    // Linear search (K ≤ 6). Returns -1 if not found.
    inline int index_in_span(std::span<const paratoric::Lattice::Edge> s,
                             const paratoric::Lattice::Edge& e) noexcept {
        for (size_t i = 0; i < s.size(); ++i)
            if (s[i] == e) return static_cast<int>(i);
        return -1;
    }

    // Insert (t, tag) into sorted-by-time vector with linear insertion (m tiny)
    inline void insert_sorted_by_time(std::vector<std::pair<double,int>>& v,
                                      double t, int tag) {
        auto it = v.begin();
        // Keep ascending order; stable wrt equal times
        for (; it != v.end() && it->first <= t; ++it) {}
        v.insert(it, {t, tag});
    }
} // anonymous namespace

[[gnu::hot]]
std::tuple<double, std::vector<int>, std::vector<double>> 
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

    std::vector<int> unique_vertices;
    unique_vertices.reserve(plaquette_edges.size() * 2); // each edge contributes up to 2 vertices
    std::vector<double> star_potential_energy_diffs;

    for (const auto& edg : plaquette_edges) {
        auto [source_v, target_v] = vertices_of_edge(edg);
        unique_vertices.emplace_back(source_v);
        unique_vertices.emplace_back(target_v);
    }

    std::sort(unique_vertices.begin(), unique_vertices.end());
    unique_vertices.erase(std::unique(
        unique_vertices.begin(), 
        unique_vertices.end()), 
        unique_vertices.end()
    );

    double energy = 0.;
    std::vector<std::pair<double,int>> local_spin_flips;
    for (const auto& v : unique_vertices) {
        const auto& star_edges = get_star_edges(v);
        local_spin_flips.clear();
        insert_sorted_by_time(local_spin_flips, imag_time_tuple_flip, 1);

        // Add singles that lie on edges shared with the plaquette
        for (const auto& e : star_edges) {
            auto it = std::find(plaquette_edges.begin(), plaquette_edges.end(), e);
            if (it != plaquette_edges.end()) {
                const size_t idx = static_cast<size_t>(std::distance(plaquette_edges.begin(), it));
                insert_sorted_by_time(local_spin_flips, spin_flip_lookup[idx], 0);
            }
        }
        double energy_v = integrated_tuple_energy_diff_combination(
            star_edges, imag_time_1, imag_time_2, local_spin_flips
        );
        star_potential_energy_diffs.emplace_back(energy_v);
        energy += energy_v;
    }
    return {energy, std::move(unique_vertices), std::move(star_potential_energy_diffs)};
}

[[gnu::hot]]
std::tuple<double, std::vector<int>, std::vector<double>> 
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
    std::vector<int> unique_plaquettes;

    for (const auto& edg : star_edges) {
        const auto& part_of_plaquette_lookup = g[edg].part_of_plaquette_lookup;
        unique_plaquettes.insert(unique_plaquettes.end(), part_of_plaquette_lookup.begin(), part_of_plaquette_lookup.end());
    }

    std::sort(unique_plaquettes.begin(), unique_plaquettes.end());
    unique_plaquettes.erase(std::unique(unique_plaquettes.begin(), unique_plaquettes.end()), unique_plaquettes.end());

    thread_local std::vector<std::pair<double,int>> flips;
    flips.reserve(8);

    std::vector<double> plaquette_potential_energy_diffs;
    plaquette_potential_energy_diffs.reserve(unique_plaquettes.size());

    double energy_sum = 0.0;

    for (int p_index : unique_plaquettes) {
        const auto plaq = get_plaquette_edges(p_index); // span<const Edge>
        flips.clear();
        insert_sorted_by_time(flips, imag_time_tuple_flip, 1);

        // Add singles on edges that are shared with the star
        for (const auto& pe : plaq) {
            const int idx = index_in_span(star_edges, pe);
            if (idx >= 0) {
                insert_sorted_by_time(flips, spin_flip_lookup[static_cast<size_t>(idx)], 0);
            }
        }

        const double ep = integrated_tuple_energy_diff_combination(
            plaq, imag_time_1, imag_time_2, flips);
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
        // only integrate up to beta/2 for dynamical susceptibility
        total_integrated_edge_energy += integrated_edge_energy_weighted(edg, 0, 0.5 * BETA);
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
            set_potential_star_energy(star_center, integrated_tuple_energy(sedges, 0, BETA));
        }
    } else {
        for (size_t plaquette_index = 0; plaquette_index < plaquette_vector.size(); ++plaquette_index) {
            const auto& pedges = get_plaquette_edges(plaquette_index);
            set_potential_plaquette_energy(plaquette_index, integrated_tuple_energy(pedges, 0, BETA));
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
    return {static_cast<double>(integrated_magnetization / (double)(get_edge_count()) ), static_cast<double>(magnetization)}; 
}

std::complex<double> Lattice::get_diag_dynamical_M_M() {
    const int magnetization = std::accumulate(
        egde_cache_.begin(), 
        egde_cache_.end(), 
        0, 
        [&](int lhs, const Edge& rhs) {return lhs + get_spin(rhs);}
    );

    double integrated_magnetization = total_integrated_edge_energy_weighted();
    return {static_cast<double>(integrated_magnetization / (double)(get_edge_count()) ), static_cast<double>(magnetization)}; 
}

std::complex<double> Lattice::get_non_diag_M_M() {
    double k_total = 0.0;
    for (const auto& edg : egde_cache_) {
        k_total += g[edg].single_spin_flips.size();
    }
    // Return raw count in .real(); imag unused (set =0 or copy for compatibility)
    return {k_total, k_total};
}

// Antiderivative of w(tau) = tau (i.e. min(tau, beta - tau) on [0, beta/2])
[[gnu::always_inline]] inline double W_tri_half(double t) {
    // assumes 0 <= t <= beta/2
    return 0.5 * t * t;
}

// \int_a^b w(tau) dtau with w(tau)=tau on [0, beta/2]
[[gnu::always_inline]] inline double integral_w_tri_half(double a, double b) {
    return W_tri_half(b) - W_tri_half(a);
}

// ===== Dynamical / fidelity susceptibility kernel =====
// Integrate sigma_e(tau) * w(tau) over [imag_time_1, imag_time_2]
// using w(tau)=min(tau, beta-tau). We only integrate [0, beta/2] here
// (callers ensure that) to avoid double counting symmetric time ordering.
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

    const auto& spin_flips = g[edg].spin_flips;  // sorted in [0, BETA)
    // find flips in [imag_time_1, imag_time_2)
    auto lo = std::lower_bound(spin_flips.begin(), spin_flips.end(), imag_time_1);
    auto hi = std::upper_bound(lo,               spin_flips.end(),  imag_time_2);

    // spin just after imag_time_1 (same parity logic as your unweighted version)
    const int base_spin = get_spin(edg); // value just after tau=0
    int spin = (((lo - spin_flips.begin()) & 1) ? -base_spin : base_spin);

    double weighted = 0.0;
    double t_prev   = imag_time_1;

    // accumulate piecewise-constant segments with triangular weight on [0, beta/2]
    for (auto it = lo; it != hi; ++it) {
        const double t_curr = *it;
        // add spin * \int_{t_prev}^{t_curr} w(tau) dtau
        weighted += spin * integral_w_tri_half(t_prev, t_curr);
        spin = -spin;
        t_prev = t_curr;
    }

    // tail to imag_time_2
    if (t_prev < imag_time_2) {
        weighted += spin * integral_w_tri_half(t_prev, imag_time_2);
    }

    return weighted;
}

std::complex<double> Lattice::get_kL_kR_single() {
    // REQUIREMENT: call rotate_imag_time() just before measuring,
    // so the cut at tau=0, beta/2 is uniformly random each time.
    double kL = 0.0, kR = 0.0;
    const double half = 0.5 * BETA;

    for (const auto& edg : egde_cache_) {
        const auto& flips = g[edg].single_spin_flips; // sorted, in [0, beta)
        // Count into halves without branches in the inner loop
        // (linear scan is faster than two binary searches per edge here
        //  because we must visit all elements anyway).
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
    // energy_beta_lmbda is the gauge field energy multiplied by lmbda and beta. The returned value is the gauge field energy multiplied by lmbda
    return energy_beta_lmbda / BETA;
}

double Lattice::get_non_diag_single_energy_z() {
    double energy_beta_h = 0.;
    for (const auto& edg : egde_cache_) {
        energy_beta_h += g[edg].single_spin_flips.size();
    }
    // energy_beta_h is the electric field energy multiplied by h and beta. The returned value is the electric field energy multiplied by h
    return energy_beta_h / BETA;
}

double Lattice::get_non_diag_tuple_energy_x() {
    double energy_beta_J = 0.;
    for (size_t plaquette_index = 0; plaquette_index < plaquette_vector.size(); ++plaquette_index) {
        energy_beta_J += plaquette_flip_vector[plaquette_index].size();
    }
    // energy_beta_J is the plaquette energy term multiplied by J and beta. The returned value is the plaquette energy multiplied by J
    return energy_beta_J / BETA;
}

double Lattice::get_non_diag_tuple_energy_z() {
    double energy_beta_mu = 0.;
    for (size_t star_center = 0; star_center < (size_t)get_vertex_count(); ++star_center) {
        energy_beta_mu += g[star_center].star_flips.size();
    }
    // energy_beta_mu is the star energy term multiplied by mu and beta. The returned value is the star energy multiplied by mu
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
                    if (get_spin(edge_in_between(v, neighbor)) == -1 && !discovered[neighbor])
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
                    // Only traverse plaquettes with active status.
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
        //TODO
        percolating = 0;
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
        //TODO
        return -1.;
    }
}

double Lattice::cube_percolation_strength() {
    return -1.;
}

double Lattice::cube_percolation_probability() {
    if (BOUNDARIES == "periodic") {
        return is_winding_cube_percolating();
    } else {
        //TODO
        return -1.;
    }
}

void Lattice::rotate_imag_time() {
    std::uniform_real_distribution<double> new_times_dist(0, BETA);
    const double tau_0 = new_times_dist(*rng);
    for (const auto& edg : egde_cache_) {
        auto& single_spin_flips = g[edg].single_spin_flips;
        auto& spin_flips = g[edg].spin_flips;
        auto it_single = std::lower_bound(single_spin_flips.begin(), single_spin_flips.end(), tau_0);
        auto it = std::lower_bound(spin_flips.begin(), spin_flips.end(), tau_0);
        size_t pivot_index = std::distance(spin_flips.begin(), it);
        // flips crossing the cut = pivot_index
        if (pivot_index % 2 == 1) {
            g[edg].spin *= -1;
        }
        // rotate so that pivot becomes first element
        std::rotate(single_spin_flips.begin(), it_single, single_spin_flips.end());
        std::rotate(spin_flips.begin(), it, spin_flips.end());
        // now wrap times by subtracting tau_0 and bringing into [0, beta)
        for (double& t : single_spin_flips) {
            t = modulo(t - tau_0, BETA);
            if (t==0) t += std::numeric_limits<double>::epsilon();
        }
        for (double& t : spin_flips) {
            t = modulo(t - tau_0, BETA);
            if (t==0) t += std::numeric_limits<double>::epsilon();
        }
    }

    if (BASIS == 'x') {
        for (int p_index = 0; p_index < get_plaquette_count(); ++p_index) {
            auto& plaquette_spin_flips = plaquette_flip_vector[p_index];
            auto it = std::lower_bound(plaquette_spin_flips.begin(), plaquette_spin_flips.end(), tau_0);
            std::rotate(plaquette_spin_flips.begin(), it, plaquette_spin_flips.end());
            for (double& t : plaquette_spin_flips) {
                t = modulo(t - tau_0, BETA);
                if (t==0) t += std::numeric_limits<double>::epsilon();
            }
        }
    } else {
        for (int s_index = 0; s_index < get_vertex_count(); ++s_index) {
            auto& star_spin_flips = g[s_index].star_flips;
            auto it = std::lower_bound(star_spin_flips.begin(), star_spin_flips.end(), tau_0);
            std::rotate(star_spin_flips.begin(), it, star_spin_flips.end());
            for (double& t : star_spin_flips) {
                t = modulo(t - tau_0, BETA);
                if (t==0) t += std::numeric_limits<double>::epsilon();
            }
        } 
    }

}

void Lattice::update_spin_string() {
    for (const auto& edg : egde_cache_) {
        auto& g_edg = g[edg];
        if (g_edg.spin_string.empty()) {
            g_edg.spin_string = std::to_string(g_edg.spin);
        } else {
            g_edg.spin_string += " " + std::to_string(g_edg.spin);
        }
    }
}

void Lattice::write_graph(const std::string& file_name, const std::filesystem::path& output_directory) {
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
