// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#define BOOST_TEST_MODULE TestLatticeGeometry

#include "lattice/lattice_geometry.hpp"
#include "lattice/lattice.hpp"
#include "paratoric/mcmc/extended_toric_code.hpp"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>

namespace paratoric {
namespace {

// Stable FNV-1a over lengths, indices, labels and coordinates rounded to 1e-12.
// Rounding avoids depending on libm/compiler last-bit differences in sqrt(3).
struct GeometryHash {
    std::uint64_t value = 14695981039346656037ULL;
    void add(std::uint64_t n) {
        for (int i = 0; i < 8; ++i) {
            value = (value ^ (n & 255)) * 1099511628211ULL;
            n >>= 8;
        }
    }
    void add(int n) { add(static_cast<std::uint64_t>(n)); }
    void add(double n) { add(static_cast<std::uint64_t>(std::llround(n * 1e12))); }
    void add(const std::string& s) {
        add(static_cast<std::uint64_t>(s.size()));
        for (unsigned char c : s) add(static_cast<int>(c));
    }
    void add(const std::pair<int, int>& p) { add(p.first); add(p.second); }
    template<class T, std::size_t N> void add(const std::array<T, N>& a) {
        for (const auto& x : a) add(x);
    }
    template<class T> void add(const std::vector<T>& v) {
        add(static_cast<std::uint64_t>(v.size()));
        for (const auto& x : v) add(x);
    }
    void add(const paratoric::LatticeGeometry& g) {
        add(g.dimensionality); add(g.vertices); add(static_cast<std::uint64_t>(g.edges.size()));
        for (const auto& e : g.edges) { add(e.source); add(e.target); add(e.orientation); }
        add(g.plaquettes); add(g.plaquette_coordinates);
        add(g.cubes); add(g.cube_coordinates);
        add(g.plaquette_cubes); add(g.cube_plaquettes);
        add(g.max_coordinates); add(g.max_plaquette_coordinates);
    }
};

struct Fixture {
    const char* type;
    const char* boundaries;
    int size;
    std::uint64_t geometry, x_paths, z_paths;
};

// Frozen from Lattice before the unit-cell extraction. These cover the complete
// ordered mappings, coordinates, orientations, cube incidence, and both paths.
// In particular L=2 retains parallel edges; odd open patches exercise row cuts.
// Kagome's subsequently appended hexagons are tested separately below.
constexpr Fixture fixtures[] = {
    {"square", "open", 2, 200113429711695586ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"square", "open", 3, 10077877769477543934ULL, 9808874869469701221ULL, 10145272227952522563ULL},
    {"square", "open", 4, 17444478573154458582ULL, 9808874869469701221ULL, 7064727803264080779ULL},
    {"square", "open", 6, 5793185067320925546ULL, 1992873620315726922ULL, 5292724192319397515ULL},
    {"square", "open", 7, 4241993535162558462ULL, 1123047068572142259ULL, 372824341776044659ULL},
    {"square", "open", 8, 17164768143674437670ULL, 1257959852287882354ULL, 2997324090125297977ULL},
    {"square", "open", 10, 15980311002235452694ULL, 3772072762501218510ULL, 3824966090882466041ULL},
    {"square", "periodic", 2, 2317684135802314551ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"square", "periodic", 3, 7743146729057167928ULL, 9808874869469701221ULL, 10145272227952522563ULL},
    {"square", "periodic", 4, 3480313625180775991ULL, 9808874869469701221ULL, 7064727803264080779ULL},
    {"square", "periodic", 6, 16834230515723417631ULL, 1992873620315726922ULL, 5292724192319397515ULL},
    {"square", "periodic", 7, 10321605866418449616ULL, 1123047068572142259ULL, 372824341776044659ULL},
    {"square", "periodic", 8, 2680931725857083663ULL, 1257959852287882354ULL, 2997324090125297977ULL},
    {"square", "periodic", 10, 4372363485375877063ULL, 3772072762501218510ULL, 3824966090882466041ULL},
    {"cubic", "open", 2, 3731207052943570678ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"cubic", "open", 3, 10574663209668341391ULL, 9808874869469701221ULL, 6830983799841977441ULL},
    {"cubic", "open", 4, 5572079754351308500ULL, 9808874869469701221ULL, 5375012402013025163ULL},
    {"cubic", "open", 6, 84474618372386114ULL, 9808874869469701221ULL, 2276932759680995467ULL},
    {"cubic", "open", 7, 3620289494339327709ULL, 9808874869469701221ULL, 11871278925928115465ULL},
    {"cubic", "open", 8, 12084042511888095190ULL, 9808874869469701221ULL, 2984738049277565241ULL},
    {"cubic", "open", 10, 7930322342178127751ULL, 9808874869469701221ULL, 17869729717950989641ULL},
    {"cubic", "periodic", 2, 18100441861668158246ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"cubic", "periodic", 3, 801598179483732066ULL, 9808874869469701221ULL, 6830983799841977441ULL},
    {"cubic", "periodic", 4, 2211211732931204642ULL, 9808874869469701221ULL, 5375012402013025163ULL},
    {"cubic", "periodic", 6, 11921006434455048042ULL, 9808874869469701221ULL, 2276932759680995467ULL},
    {"cubic", "periodic", 7, 13812722845035732246ULL, 9808874869469701221ULL, 11871278925928115465ULL},
    {"cubic", "periodic", 8, 9441222207768944734ULL, 9808874869469701221ULL, 2984738049277565241ULL},
    {"cubic", "periodic", 10, 20251919558569774ULL, 9808874869469701221ULL, 17869729717950989641ULL},
    {"triangular", "open", 2, 3906376326282813282ULL, 9808874869469701221ULL, 6457745870321885959ULL},
    {"triangular", "open", 3, 16289935744546812982ULL, 9808874869469701221ULL, 10145272227952522563ULL},
    {"triangular", "open", 4, 4019634877140554371ULL, 9808874869469701221ULL, 16962088202193734405ULL},
    {"triangular", "open", 6, 8462180631885355084ULL, 12880706973509946286ULL, 2125895707832625207ULL},
    {"triangular", "open", 7, 16173532909036861991ULL, 17190742282341182807ULL, 8761791350414751095ULL},
    {"triangular", "open", 8, 15601456696237076223ULL, 2511599401459353751ULL, 5839415122659254823ULL},
    {"triangular", "open", 10, 17824903125779283153ULL, 16814467487480871960ULL, 2040421262301380123ULL},
    {"triangular", "periodic", 2, 6311343380598122347ULL, 9808874869469701221ULL, 6457745870321885959ULL},
    {"triangular", "periodic", 4, 14874946291079815893ULL, 9808874869469701221ULL, 16962088202193734405ULL},
    {"triangular", "periodic", 6, 3285829398290374365ULL, 12880706973509946286ULL, 2125895707832625207ULL},
    {"triangular", "periodic", 8, 11491728909115181918ULL, 2511599401459353751ULL, 5839415122659254823ULL},
    {"triangular", "periodic", 10, 16043840949376920772ULL, 16814467487480871960ULL, 2040421262301380123ULL},
    {"honeycomb", "open", 2, 10585379519528241336ULL, 9808874869469701221ULL, 18018062118808294625ULL},
    {"honeycomb", "open", 3, 18059023214442902580ULL, 9808874869469701221ULL, 17585130257612631915ULL},
    {"honeycomb", "open", 4, 7312734739200484171ULL, 9808874869469701221ULL, 12238530042779350871ULL},
    {"honeycomb", "open", 6, 8382934561351601419ULL, 13374754114434524173ULL, 11384278665230500623ULL},
    {"honeycomb", "open", 7, 2667298431016959897ULL, 537466728785440638ULL, 9273591727704858683ULL},
    {"honeycomb", "open", 8, 13515521134832330150ULL, 126387554707291381ULL, 8641323455277156101ULL},
    {"honeycomb", "open", 10, 1917242308757113835ULL, 1264457783743445365ULL, 15999468852864768835ULL},
    {"honeycomb", "periodic", 2, 5541556987098984060ULL, 9808874869469701221ULL, 6414116535217072903ULL},
    {"honeycomb", "periodic", 4, 4182356791143992845ULL, 9808874869469701221ULL, 6565743850930031945ULL},
    {"honeycomb", "periodic", 6, 5066089957510237966ULL, 5851487749558412817ULL, 11628831302145246261ULL},
    {"honeycomb", "periodic", 8, 1337723711847180655ULL, 10272009829432008445ULL, 13834735286393312089ULL},
    {"honeycomb", "periodic", 10, 1515734113536598989ULL, 7552562834634194221ULL, 11186151008720176835ULL},
    {"kagome", "open", 2, 16418360156981121458ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "open", 3, 3544469849414119040ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "open", 4, 14308442984121275134ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "open", 6, 13557412802273759672ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "open", 7, 7978587069643096611ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "open", 8, 15314494415902196165ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "open", 10, 11631914461368872581ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "periodic", 2, 1350727967018941178ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "periodic", 4, 12952448298251881827ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "periodic", 6, 4508083414803766684ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "periodic", 8, 4185152660220154320ULL, 9808874869469701221ULL, 9808874869469701221ULL},
    {"kagome", "periodic", 10, 3177135575880490573ULL, 9808874869469701221ULL, 9808874869469701221ULL},
};

bool paths_throw_in_debug(const Fixture& f, char basis) {
#ifdef NDEBUG
    (void)f;
    (void)basis;
    return false;
#else
    return (basis == 'x' && f.size < 6) || std::string(f.type) == "kagome"
        || (basis == 'z' && std::string(f.type) == "honeycomb"
            && std::string(f.boundaries) == "open" && f.size < 3);
#endif
}

} // namespace

BOOST_AUTO_TEST_CASE(geometry_and_paths_match_legacy_mappings) {
    for (const auto& f : fixtures) {
        BOOST_TEST_CONTEXT(f.type << " " << f.boundaries << " L=" << f.size) {
            const auto g = make_lattice_geometry(f.type, f.size, f.boundaries);
            GeometryHash hash;
            if (std::string(f.type) == "kagome") {
                // Keep checking the original triangles and their indices even
                // though Kagome now appends hexagons to the same registry.
                auto triangles = g;
                const int count = std::string(f.boundaries) == "periodic"
                    ? 2 * f.size * f.size : f.size * f.size + (f.size - 1) * (f.size - 1);
                triangles.plaquettes.resize(count);
                triangles.plaquette_coordinates.resize(count);
                triangles.plaquette_cubes.resize(count);
                hash.add(triangles);
            } else {
                hash.add(g);
            }
            BOOST_CHECK_EQUAL(hash.value, f.geometry);
            for (char basis : {'x', 'z'}) {
                if (paths_throw_in_debug(f, basis)) {
                    BOOST_CHECK_THROW(make_fredenhagen_marcu_paths(g, f.type, f.size, f.boundaries, basis), std::exception);
                } else {
                    const auto [half, full] = make_fredenhagen_marcu_paths(g, f.type, f.size, f.boundaries, basis);
                    GeometryHash paths;
                    paths.add(half);
                    paths.add(full);
                    BOOST_CHECK_EQUAL(paths.value, basis == 'x' ? f.x_paths : f.z_paths);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(lattice_consumes_mappings_and_initializes_state) {
    for (const auto& f : fixtures) {
        const auto g = make_lattice_geometry(f.type, f.size, f.boundaries);
        for (char basis : {'x', 'z'}) {
            if (paths_throw_in_debug(f, basis)) continue;
            for (int spin : {-1, 1}) {
                BOOST_TEST_CONTEXT(f.type << " " << f.boundaries << " L=" << f.size << " " << basis << " spin=" << spin) {
                    const double beta = 1.25;
                    Lattice lat(LatSpec{basis, f.type, f.size, beta, f.boundaries, spin},
                        std::make_shared<Lattice::RNG>(793));
                    BOOST_CHECK_EQUAL(lat.get_vertex_count(), g.vertices.size());
                    BOOST_CHECK_EQUAL(lat.get_edge_count(), g.edges.size());
                    BOOST_CHECK_EQUAL(lat.get_plaquette_count(), g.plaquettes.size());
                    BOOST_CHECK_EQUAL(lat.get_cube_count(), g.cubes.size());
                    for (int v = 0; v < lat.get_vertex_count(); ++v) {
                        // Incident edges preserve insertion order, including both
                        // occurrences of a self-loop and all parallel edges.
                        std::vector<const LatticeGeometry::Edge*> expected;
                        for (const auto& edge : g.edges) {
                            if (edge.source == v) expected.push_back(&edge);
                            if (edge.target == v) expected.push_back(&edge);
                        }
                        const auto star = lat.get_star_edges(v);
                        BOOST_REQUIRE_EQUAL(star.size(), expected.size());
                        for (std::size_t i = 0; i < star.size(); ++i) {
                            const auto& e = *expected[i];
                            // Boost orients an incident descriptor away from
                            // the queried star, even when inserted in reverse.
                            BOOST_CHECK(lat.vertices_of_edge(star[i]) == std::make_pair(
                                v, e.source == v ? e.target : e.source));
                            BOOST_CHECK_EQUAL(lat.get_orientation(star[i]), e.orientation);
                            BOOST_CHECK_EQUAL(lat.get_spin(star[i]), spin);
                            BOOST_CHECK_EQUAL(lat.get_potential_edge_energy(star[i]), beta * spin);
                            BOOST_CHECK_EQUAL(lat.get_spin_flip_count(star[i]), 0);
                            BOOST_CHECK(lat.get_single_spin_flips(star[i]).empty());
                        }
                        if (basis == 'x') {
                            BOOST_CHECK_EQUAL(lat.get_potential_star_energy(v), beta * (star.size() % 2 ? spin : 1));
                        } else {
                            BOOST_CHECK(lat.get_tuple_spin_flips(v).empty());
                        }
                    }
                    for (int p = 0; p < lat.get_plaquette_count(); ++p) {
                        BOOST_CHECK(lat.get_plaquette_vertex_pairs(p) == g.plaquettes[p]);
                        const auto edges = lat.get_plaquette_edges(p);
                        BOOST_REQUIRE_EQUAL(edges.size(), g.plaquettes[p].size());
                        for (std::size_t i = 0; i < edges.size(); ++i) {
                            const auto [u, v] = g.plaquettes[p][i];
                            BOOST_CHECK(edges[i] == lat.edge_in_between(u, v));
                            BOOST_CHECK_EQUAL(g.plaquettes[p][i].second,
                                g.plaquettes[p][(i + 1) % edges.size()].first);
                        }
                        if (basis == 'z') {
                            BOOST_CHECK_EQUAL(lat.get_potential_plaquette_energy(p), beta * (edges.size() % 2 ? spin : 1));
                        } else {
                            BOOST_CHECK(lat.get_tuple_spin_flips(p).empty());
                        }
                    }
                    for (int c = 0; c < lat.get_cube_count(); ++c)
                        BOOST_CHECK(lat.get_cube_vertices(c) == g.cubes[c]);

                    // A separate RNG gives expected insertion indices without
                    // relying on implementation-specific std distributions.
                    Lattice::RNG expected_rng(793);
                    for (int i = 0; i < 100; ++i) {
                        const auto index = rng::uniform_index(expected_rng, g.edges.size());
                        const auto [edge, u, v] = lat.get_random_edge();
                        BOOST_CHECK_EQUAL(u, g.edges[index].source);
                        BOOST_CHECK_EQUAL(v, g.edges[index].target);
                        BOOST_CHECK_EQUAL(lat.get_orientation(edge), g.edges[index].orientation);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(unit_cell_counts_and_boundary_clipping) {
    for (int L : {1, 2, 3, 6, 7}) {
        for (const std::string boundary : {"open", "periodic"}) {
            const bool periodic = boundary == "periodic";
            const auto square = make_lattice_geometry("square", L, boundary);
            BOOST_CHECK_EQUAL(square.vertices.size(), L * L);
            BOOST_CHECK_EQUAL(square.edges.size(), periodic ? 2 * L * L : 2 * L * (L - 1));
            BOOST_CHECK_EQUAL(square.plaquettes.size(), periodic ? L * L : (L - 1) * (L - 1));
            const auto cubic = make_lattice_geometry("cubic", L, boundary);
            BOOST_CHECK_EQUAL(cubic.vertices.size(), L * L * L);
            BOOST_CHECK_EQUAL(cubic.edges.size(), periodic ? 3 * L * L * L : 3 * L * L * (L - 1));
            BOOST_CHECK_EQUAL(cubic.plaquettes.size(), periodic ? 3 * L * L * L : 3 * L * (L - 1) * (L - 1));
            // Deliberate compatibility with the old x-wrapping open cubes.
            BOOST_CHECK_EQUAL(cubic.cubes.size(), periodic ? L * L * L : L * (L - 1) * (L - 1));
            if (periodic && L % 2) continue;
            const auto triangular = make_lattice_geometry("triangular", L, boundary);
            BOOST_CHECK_EQUAL(triangular.vertices.size(), L * L);
            BOOST_CHECK_EQUAL(triangular.edges.size(), periodic ? 3 * L * L : (L - 1) * (3 * L - 1));
            BOOST_CHECK_EQUAL(triangular.plaquettes.size(), periodic ? 2 * L * L : 2 * (L - 1) * (L - 1));
            const auto honeycomb = make_lattice_geometry("honeycomb", L, boundary);
            BOOST_CHECK_EQUAL(honeycomb.vertices.size(), periodic ? 2 * L * L : (2 * L + 2) * (L + 1) - 2);
            BOOST_CHECK_EQUAL(honeycomb.edges.size(), periodic ? 3 * L * L : 3 * L * L + 4 * L - 1);
            BOOST_CHECK_EQUAL(honeycomb.plaquettes.size(), L * L);
            const auto kagome = make_lattice_geometry("kagome", L, boundary);
            BOOST_CHECK_EQUAL(kagome.vertices.size(), 3 * L * L);
            BOOST_CHECK_EQUAL(kagome.edges.size(), periodic ? 6 * L * L : 3 * L * L + (L - 1) * (3 * L - 1));
            BOOST_CHECK_EQUAL(kagome.plaquettes.size(), periodic ? 3 * L * L : L * L + 2 * (L - 1) * (L - 1));
        }
    }
}

BOOST_AUTO_TEST_CASE(kagome_triangles_and_hexagons_tile_the_lattice) {
    for (int L : {1, 2, 3, 4, 6, 7}) {
        for (const std::string boundary : {"open", "periodic"}) {
            const bool periodic = boundary == "periodic";
            if (periodic && L % 2) continue;
            BOOST_TEST_CONTEXT(boundary << " L=" << L) {
                const auto g = make_lattice_geometry("kagome", L, boundary);
                const int triangles = periodic ? 2 * L * L : L * L + (L - 1) * (L - 1);
                const int hexagons = periodic ? L * L : (L - 1) * (L - 1);
                BOOST_REQUIRE_EQUAL(g.plaquettes.size(), triangles + hexagons);
                std::map<LatticeGeometry::VertexPair, std::array<int, 2>> incidence;
                for (const auto& edge : g.edges) {
                    BOOST_REQUIRE(incidence.emplace(std::minmax(edge.source, edge.target),
                        std::array<int, 2>{}).second);
                }
                std::set<std::set<int>> distinct_faces;
                for (std::size_t p = 0; p < g.plaquettes.size(); ++p) {
                    const auto& face = g.plaquettes[p];
                    BOOST_REQUIRE_EQUAL(face.size(), p < triangles ? 3 : 6);
                    std::set<int> vertices;
                    for (std::size_t i = 0; i < face.size(); ++i) {
                        const auto [u, v] = face[i];
                        BOOST_CHECK_EQUAL(v, face[(i + 1) % face.size()].first);
                        const auto edge = incidence.find(std::minmax(u, v));
                        BOOST_REQUIRE(edge != incidence.end());
                        ++edge->second[face.size() == 3 ? 0 : 1];
                        vertices.insert(u);
                    }
                    BOOST_CHECK_EQUAL(vertices.size(), face.size());
                    BOOST_CHECK(distinct_faces.insert(vertices).second);
                }
                for (const auto& [edge, counts] : incidence) {
                    if (periodic) {
                        BOOST_CHECK_EQUAL(counts[0], 1);
                        BOOST_CHECK_EQUAL(counts[1], 1);
                    } else {
                        BOOST_CHECK(counts[0] <= 1 && counts[1] <= 1);
                        BOOST_CHECK(counts[0] + counts[1] >= 1);
                    }
                }
                BOOST_CHECK_EQUAL(static_cast<int>(g.vertices.size())
                    - static_cast<int>(g.edges.size()) + static_cast<int>(g.plaquettes.size()),
                    periodic ? 0 : 1);
            }
        }
    }
}

#ifdef NDEBUG
// Kagome's unimplemented Fredenhagen-Marcu paths still prevent constructing
// Lattice in debug builds. Geometry itself is tested above in both builds.
BOOST_AUTO_TEST_CASE(kagome_hexagons_participate_in_qmc) {
    for (const std::string boundary : {"open", "periodic"}) {
        const int L = 6;
        const int triangles = boundary == "periodic" ? 2 * L * L : L * L + (L - 1) * (L - 1);
        const int hexagons = boundary == "periodic" ? L * L : (L - 1) * (L - 1);
        Lattice lat(LatSpec{'z', "kagome", L, 1.25, boundary, -1});
        BOOST_CHECK_EQUAL(lat.get_diag_tuple_energy_z(), hexagons - triangles);
        BOOST_CHECK_EQUAL(lat.total_integrated_plaquette_energy(), 1.25 * (hexagons - triangles));

        Lattice xlat(LatSpec{'x', "kagome", L, 1.25, boundary, 1});
        const auto hexagon = xlat.get_plaquette_edges(triangles);
        BOOST_REQUIRE_EQUAL(hexagon.size(), 6);
        xlat.insert_double_tuple_flip(triangles, hexagon, 0.25, 0.75);
        BOOST_CHECK_EQUAL(xlat.get_tuple_spin_flips(triangles).size(), 2);
        for (const auto& edge : hexagon) BOOST_CHECK_EQUAL(xlat.get_spin_flip_count(edge), 2);
        BOOST_CHECK_EQUAL(xlat.get_non_diag_tuple_energy_x(), 2. / 1.25);

        for (char basis : {'x', 'z'}) {
            Config cfg;
            cfg.lat_spec = {basis, "kagome", L, 1.25, boundary, 1};
            cfg.param_spec.h = 0.7;
            cfg.param_spec.lmbda = 0.6;
            cfg.sim_spec.seed = 19381;
            cfg.sim_spec.N_thermalization = 100;
            cfg.sim_spec.N_between_samples = 50;
            cfg.sim_spec.N_samples = 20;
            cfg.sim_spec.N_resamples = 10;
            cfg.sim_spec.observables = {"energy", "energy_J", "plaquette_z", "anyon_count"};
            const auto result = ExtendedToricCode::get_sample(cfg);
            BOOST_REQUIRE_EQUAL(result.series.size(), cfg.sim_spec.observables.size());
            for (const auto& series : result.series) {
                BOOST_REQUIRE_EQUAL(series.size(), cfg.sim_spec.N_samples);
                for (const auto& value : series) BOOST_CHECK(std::isfinite(std::get<double>(value)));
            }
        }
    }
}
#endif

BOOST_AUTO_TEST_CASE(fredenhagen_marcu_uses_generated_paths) {
    for (const auto& f : fixtures) {
        if (f.size != 6) continue;
        const auto geometry = make_lattice_geometry(f.type, f.size, f.boundaries);
        for (char basis : {'x', 'z'}) {
            if (paths_throw_in_debug(f, basis)) continue;
            BOOST_TEST_CONTEXT(f.type << " " << f.boundaries << " " << basis) {
                Lattice lat(LatSpec{basis, f.type, f.size, 1.25, f.boundaries, 1},
                    std::make_shared<Lattice::RNG>(983));
                const auto [half, full] = make_fredenhagen_marcu_paths(
                    geometry, f.type, f.size, f.boundaries, basis);
                const auto product = [&](const auto& path) {
                    int value = 1;
                    for (const auto& [u, v] : path) {
                        BOOST_REQUIRE(lat.exists_edge(u, v));
                        value *= lat.get_spin(lat.edge_in_between(u, v));
                    }
                    return value;
                };
                for (int sample = 0; sample < 20; ++sample) {
                    for (int i = 0; i < lat.get_edge_count() / 3; ++i)
                        lat.flip_spin(lat.get_random_edge_descriptor());
                    const auto fm = lat.fredenhagen_marcu();
                    BOOST_CHECK_EQUAL(fm.real(), product(half));
                    BOOST_CHECK_EQUAL(fm.imag(), product(full));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(geometry_input_validation) {
    BOOST_CHECK_THROW(make_lattice_geometry("unknown", 6, "open"), std::invalid_argument);
    BOOST_CHECK_THROW(make_lattice_geometry("square", 0, "open"), std::invalid_argument);
    BOOST_CHECK_THROW(make_lattice_geometry("square", -2, "open"), std::invalid_argument);
    BOOST_CHECK_THROW(make_lattice_geometry("square", 6, "unknown"), std::invalid_argument);
    for (const std::string type : {"triangular", "honeycomb", "kagome"})
        BOOST_CHECK_THROW(make_lattice_geometry(type, 3, "periodic"), std::invalid_argument);
}

} // namespace paratoric
