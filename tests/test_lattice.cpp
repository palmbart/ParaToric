// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#define BOOST_TEST_MODULE TestLattice

#include "lattice/lattice.hpp"
#include "paratoric/types/types.hpp"

#include <boost/test/unit_test.hpp>

#include <vector>
#include <string>

namespace paratoric {

BOOST_AUTO_TEST_CASE(get_non_string_count_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    }; 
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 4));
    lat.flip_spin(lat.edge_in_between(4, 5));
    BOOST_CHECK(67 == lat.get_non_string_count());
}

BOOST_AUTO_TEST_CASE(get_string_count_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK(0 == lat.get_string_count());
    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 4));
    lat.flip_spin(lat.edge_in_between(4, 5));
    BOOST_CHECK(5 == lat.get_string_count());
}

BOOST_AUTO_TEST_CASE(get_vertex_count_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK(36 == lat.get_vertex_count());
}

BOOST_AUTO_TEST_CASE(get_edge_count_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK(72 == lat.get_edge_count());
}

BOOST_AUTO_TEST_CASE(get_plaquette_count_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK(36 == lat.get_plaquette_count());
}

BOOST_AUTO_TEST_CASE(get_anyon_count_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK(0 == lat.get_anyon_count());
}

BOOST_AUTO_TEST_CASE(get_anyon_count_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    BOOST_CHECK(2 == lat.get_anyon_count());
}

BOOST_AUTO_TEST_CASE(get_spin_flip_imag_time_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    BOOST_CHECK(2.0 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 0));
}

BOOST_AUTO_TEST_CASE(get_spin_flip_index_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    BOOST_CHECK(0 == lat.get_spin_flip_index(lat.edge_in_between(0,1), 2.0));
    BOOST_CHECK(1 == lat.get_spin_flip_index(lat.edge_in_between(0,1), 3.0));
}

BOOST_AUTO_TEST_CASE(set_spin_flip_imag_time_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.set_spin_flip_imag_time(lat.edge_in_between(0,1), 1, 2.5);
    BOOST_CHECK(2.5 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 1));
}

BOOST_AUTO_TEST_CASE(delete_double_single_spin_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.delete_double_single_spin_flip(lat.edge_in_between(0,1), 2.0, 3.0);
    BOOST_CHECK(0 == lat.get_spin_flip_count(lat.edge_in_between(0,1)));
}

BOOST_AUTO_TEST_CASE(delete_double_single_spin_flip_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.delete_double_single_spin_flip(lat.edge_in_between(0,1), 2.5, 2.7);
    BOOST_CHECK(2 == lat.get_spin_flip_count(lat.edge_in_between(0,1)));
}

BOOST_AUTO_TEST_CASE(delete_single_spin_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.delete_single_spin_flip(lat.edge_in_between(0,1), 2);
    BOOST_CHECK(3.0 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 2));
}

BOOST_AUTO_TEST_CASE(delete_single_spin_flip_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.delete_single_spin_flip(lat.edge_in_between(0,1), 1);
    BOOST_CHECK(2.7 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 1));
}

BOOST_AUTO_TEST_CASE(delete_double_tuple_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 3.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.7);
    lat.delete_double_tuple_flip(0, lat.get_plaquette_edges(0), 2.5, 2.7);
    BOOST_CHECK(2 == lat.get_spin_flip_count(lat.edge_in_between(0,1)));
}

BOOST_AUTO_TEST_CASE(delete_tuple_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 3.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.7);
    lat.delete_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    BOOST_CHECK(3 == lat.get_spin_flip_count(lat.edge_in_between(0,1)));
}

BOOST_AUTO_TEST_CASE(check_tuple_flip_present_tuple_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 3.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.7);
    lat.delete_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    BOOST_CHECK(false == lat.check_tuple_flip_present_tuple(lat.get_plaquette_edges(0), 2.5));
    BOOST_CHECK(true == lat.check_tuple_flip_present_tuple(lat.get_plaquette_edges(0), 2.7));
    BOOST_CHECK(false == lat.check_tuple_flip_present_tuple(lat.get_plaquette_edges(0), 2.4));
}

BOOST_AUTO_TEST_CASE(check_plaquette_flip_at_edge_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 3.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.7);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(1), 2.5);
    BOOST_CHECK(true == lat.check_plaquette_flip_at_edge(lat.edge_in_between(1,7), 2.5));
}

BOOST_AUTO_TEST_CASE(check_spin_flips_present_tuple_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 3.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.7);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(1), 2.5);
    BOOST_CHECK(true == lat.check_spin_flips_present_tuple(lat.get_plaquette_edges(1), 2.5));
    BOOST_CHECK(false == lat.check_spin_flips_present_tuple(lat.get_plaquette_edges(0), 2.4));
}

BOOST_AUTO_TEST_CASE(flip_next_imag_time_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    BOOST_CHECK(3.0 == lat.flip_next_imag_time(lat.edge_in_between(0,1), 2.7));
}

BOOST_AUTO_TEST_CASE(flip_next_imag_time_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK(2.7 == lat.flip_next_imag_time(lat.edge_in_between(0,1), 2.7));
}

BOOST_AUTO_TEST_CASE(flip_next_imag_times_tuple_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    BOOST_CHECK( (std::vector<double>{3.0, 2.0, 2.0, 2.0}) 
    == lat.flip_next_imag_times_tuple(lat.get_plaquette_edges(0), 2.0));
    BOOST_CHECK( (std::vector<double>{2.0, 2.0, 2.0, 2.0}) 
    == lat.flip_next_imag_times_tuple(lat.get_plaquette_edges(0), 1.9));
    BOOST_CHECK( (std::vector<double>{2.0, 2.0, 2.0, 2.0}) 
    == lat.flip_next_imag_times_tuple(lat.get_plaquette_edges(0), 3.0));
}

BOOST_AUTO_TEST_CASE(flip_next_imag_times_tuple_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK( (std::vector<double>{2.0, 2.0, 2.0, 2.0}) 
    == lat.flip_next_imag_times_tuple(lat.get_plaquette_edges(0), 2.0));
}

BOOST_AUTO_TEST_CASE(flip_prev_imag_time_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    BOOST_CHECK(2.5 == lat.flip_prev_imag_time(lat.edge_in_between(0,1), 2.7));
}

BOOST_AUTO_TEST_CASE(flip_prev_imag_time_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK(2.7 == lat.flip_prev_imag_time(lat.edge_in_between(0,1), 2.7));
}

BOOST_AUTO_TEST_CASE(flip_prev_imag_times_tuple_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 1.0);
    BOOST_CHECK( (std::vector<double>{1.0, 2.0, 2.0, 2.0}) 
    == lat.flip_prev_imag_times_tuple(lat.get_plaquette_edges(0), 2.0));
    BOOST_CHECK( (std::vector<double>{2.0, 2.0, 2.0, 2.0}) 
    == lat.flip_prev_imag_times_tuple(lat.get_plaquette_edges(0), 2.1));
    BOOST_CHECK( (std::vector<double>{2.0, 2.0, 2.0, 2.0}) 
    == lat.flip_prev_imag_times_tuple(lat.get_plaquette_edges(0), 1.0));
}

BOOST_AUTO_TEST_CASE(flip_prev_imag_times_tuple_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK( (std::vector<double>{2.0, 2.0, 2.0, 2.0}) 
    == lat.flip_prev_imag_times_tuple(lat.get_plaquette_edges(0), 2.0));
}

BOOST_AUTO_TEST_CASE(insert_double_single_spin_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_double_single_spin_flip(lat.edge_in_between(0,1), 2.0, 3.0);
    BOOST_CHECK(2 == lat.get_spin_flip_count(lat.edge_in_between(0,1)));
}

BOOST_AUTO_TEST_CASE(insert_double_single_spin_flip_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_double_single_spin_flip(lat.edge_in_between(0,1), 2.0, 3.0);
    lat.insert_double_single_spin_flip(lat.edge_in_between(0,1), 2.5, 2.7);
    lat.delete_double_single_spin_flip(lat.edge_in_between(0,1), 2.5, 2.7);
    BOOST_CHECK(2 == lat.get_spin_flip_count(lat.edge_in_between(0,1)));
}

BOOST_AUTO_TEST_CASE(insert_single_spin_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.delete_single_spin_flip(lat.edge_in_between(0,1), 2);
    BOOST_CHECK(3.0 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 2));
}

BOOST_AUTO_TEST_CASE(insert_single_spin_flip_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.delete_single_spin_flip(lat.edge_in_between(0,1), 1);
    BOOST_CHECK(2.7 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 1));
}

BOOST_AUTO_TEST_CASE(insert_double_tuple_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_double_tuple_flip(0, lat.get_plaquette_edges(0), 2.0, 3.0);
    lat.insert_double_tuple_flip(0, lat.get_plaquette_edges(0), 2.5, 2.7);
    lat.delete_double_tuple_flip(0, lat.get_plaquette_edges(0), 2.5, 2.7);
    BOOST_CHECK(2 == lat.get_spin_flip_count(lat.edge_in_between(0,1)));
}

BOOST_AUTO_TEST_CASE(insert_tuple_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 3.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.7);
    lat.delete_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    BOOST_CHECK(3 == lat.get_spin_flip_count(lat.edge_in_between(0,1)));
}

BOOST_AUTO_TEST_CASE(move_spin_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.move_spin_flip(lat.edge_in_between(0,1), 3, 1.0, false);
    BOOST_CHECK(1.0 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 0));
}

BOOST_AUTO_TEST_CASE(move_spin_flip_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 3.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.move_spin_flip(lat.edge_in_between(0,1), 0, 4.0, false);
    BOOST_CHECK(4.0 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 3));
}

BOOST_AUTO_TEST_CASE(move_tuple_flip_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 3.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.7);
    lat.move_tuple_flip(0, lat.get_plaquette_edges(0), 3.0, 1.0, false);
    BOOST_CHECK(1.0 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 0));
}

BOOST_AUTO_TEST_CASE(move_tuple_flip_test_2) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 3.0);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.5);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.7);
    lat.move_tuple_flip(0, lat.get_plaquette_edges(0), 2.0, 4.0, false);
    BOOST_CHECK(4.0 == lat.get_spin_flip_imag_time(lat.edge_in_between(0,1), 3));
}

BOOST_AUTO_TEST_CASE(integrated_tuple_energy_diff_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 1.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,6), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.95);
    BOOST_CHECK_CLOSE(1.9, lat.integrated_tuple_energy_diff(lat.get_star_edges(0), 1.5, 2.95), 0.001);
}

BOOST_AUTO_TEST_CASE(integrated_star_energy_diff_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 1.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,6), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(6,7), 2.1);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.95);
    const auto& [bare_energy, star_centers, bare_star_potential_energy_diffs] 
    = lat.integrated_star_energy_diff(lat.edge_in_between(0,6), 1.5, 2.95, false);
    BOOST_CHECK_CLOSE(-0.6, bare_energy, 0.001);
}

BOOST_AUTO_TEST_CASE(integrated_tuple_energy_test_1) {
    LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 1.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,6), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(6,7), 2.1);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.95);
    BOOST_CHECK_CLOSE(-0.95, lat.integrated_tuple_energy(lat.get_star_edges(0), 1.5, 2.95), 0.001);
}

BOOST_AUTO_TEST_CASE(tuple_energy_single_combination_test_2) {
    LatSpec spec = {
        'x', // basis
        "triangular", // lattice type
        8, // system size
        8., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 2.1);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 1.2);
    lat.insert_single_spin_flip(lat.edge_in_between(1,8), 5.8);
    lat.insert_single_spin_flip(lat.edge_in_between(8,0), 7.2);
    BOOST_CHECK_CLOSE(3.4, lat.integrated_tuple_energy(lat.get_plaquette_edges(0), 0, 8.), 0.001);
}

BOOST_AUTO_TEST_CASE(integrated_star_energy_diff_combination_test_1) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    std::vector<double> single_spin_flip_lookup {0.1, 1.2, 2.5, 4.9};

    const auto& [bare_energy_1, star_centers_1, bare_star_potential_energy_diffs_1] 
    = lat.integrated_star_energy_diff_combination(0, 0.05, 4.95, single_spin_flip_lookup, 2.1);
    const auto& [bare_energy_2, star_centers_2, bare_star_potential_energy_diffs_2] 
    = lat.integrated_star_energy_diff_combination(0, 0.1, 4.9, single_spin_flip_lookup, 2.1);

    std::println("Test1");
    BOOST_CHECK_CLOSE(-19.2, bare_energy_1, 0.001);
    std::println("");
    BOOST_CHECK_CLOSE(-19.2, bare_energy_2, 0.001);
    std::println("");
}

BOOST_AUTO_TEST_CASE(integrated_star_energy_diff_combination_test_2) {
     LatSpec spec = {
        'x', // basis
        "triangular", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    std::vector<double> single_spin_flip_lookup {0.1, 1.2, 2.5};

    const auto& [bare_energy_1, star_centers_1, bare_star_potential_energy_diffs_1] 
    = lat.integrated_star_energy_diff_combination(0, 0.05, 4.95, single_spin_flip_lookup, 2.1);
    const auto& [bare_energy_2, star_centers_2, bare_star_potential_energy_diffs_2] 
    = lat.integrated_star_energy_diff_combination(0, 0.1, 4.9, single_spin_flip_lookup, 2.1);

    std::println("Test2");
    BOOST_CHECK_CLOSE(-9.6, bare_energy_1, 0.001);
    std::println("");
    BOOST_CHECK_CLOSE(-9.6, bare_energy_2, 0.001);
    std::println("");
}

BOOST_AUTO_TEST_CASE(integrated_star_energy_diff_combination_test_3) {
     LatSpec spec = {
        'x', // basis
        "honeycomb", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    std::vector<double> single_spin_flip_lookup {0.1, 1.2, 2.5, 4.9, 5.2, 5.4};

    const auto& [bare_energy_1, star_centers_1, bare_star_potential_energy_diffs_1] 
    = lat.integrated_star_energy_diff_combination(0, 0.05, 4.95, single_spin_flip_lookup, 2.1);
    const auto& [bare_energy_2, star_centers_2, bare_star_potential_energy_diffs_2] 
    = lat.integrated_star_energy_diff_combination(0, 0.1, 4.9, single_spin_flip_lookup, 2.1);

    std::println("Test3");
    BOOST_CHECK_CLOSE(-21.2, bare_energy_1, 0.001);
    std::println("");
    BOOST_CHECK_CLOSE(-21.2, bare_energy_2, 0.001);
    std::println("");
}

BOOST_AUTO_TEST_CASE(integrated_star_energy_diff_combination_test_4) {
     LatSpec spec = {
        'x', // basis
        "cubic", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    std::vector<double> single_spin_flip_lookup {0.1, 1.2, 2.5, 4.9};

    const auto& [bare_energy_1, star_centers_1, bare_star_potential_energy_diffs_1] 
    = lat.integrated_star_energy_diff_combination(0, 0.05, 4.95, single_spin_flip_lookup, 2.1);
    const auto& [bare_energy_2, star_centers_2, bare_star_potential_energy_diffs_2] 
    = lat.integrated_star_energy_diff_combination(0, 0.1, 4.9, single_spin_flip_lookup, 2.1);

    std::println("Test4");
    BOOST_CHECK_CLOSE(-19.2, bare_energy_1, 0.001);
    std::println("");
    BOOST_CHECK_CLOSE(-19.2, bare_energy_2, 0.001);
    std::println("");
}

BOOST_AUTO_TEST_CASE(integrated_edge_energy_diff_test_1) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 1.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,6), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(6,7), 2.1);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.95);
    BOOST_CHECK_CLOSE(0.1, lat.integrated_edge_energy_diff(lat.edge_in_between(0,1), 1.5, 2.95), 0.001);
}

BOOST_AUTO_TEST_CASE(integrated_edge_energy_test_1) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 1.5);
    lat.insert_single_spin_flip(lat.edge_in_between(0,6), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(6,7), 2.1);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.0);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.7);
    lat.insert_single_spin_flip(lat.edge_in_between(0,1), 2.95);
    BOOST_CHECK_CLOSE(-0.05, lat.integrated_edge_energy(lat.edge_in_between(0,1), 1.5, 2.95), 0.001);
}

BOOST_AUTO_TEST_CASE(integrated_edge_energy_diff_combination_test_1) {
     LatSpec spec = {
        'x', // basis
        "triangular", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    auto pot_edge_energy_before = - lat.total_integrated_edge_energy();
    auto pot_star_energy_before = - lat.total_integrated_star_energy();

    auto pedges = lat.get_plaquette_edges(0);
    lat.insert_single_spin_flip(pedges[0], 5.81);
    lat.insert_single_spin_flip(pedges[1], 5.75);
    lat.insert_single_spin_flip(pedges[2], 7.33);
    lat.insert_tuple_flip(0, lat.get_plaquette_edges(0), 6.23);

    auto pot_edge_energy_after = - lat.total_integrated_edge_energy();
    auto pot_star_energy_after = - lat.total_integrated_star_energy();

    BOOST_CHECK_CLOSE(4, pot_edge_energy_after-pot_edge_energy_before, 0.001);
    BOOST_CHECK_CLOSE(6.32, pot_star_energy_after-pot_star_energy_before, 0.001);
}

BOOST_AUTO_TEST_CASE(is_percolating_test_1) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 4));
    lat.flip_spin(lat.edge_in_between(4, 5));
    BOOST_CHECK(1 == lat.is_percolating());
}

BOOST_AUTO_TEST_CASE(is_percolating_test_2) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 9));
    lat.flip_spin(lat.edge_in_between(9, 8));
    lat.flip_spin(lat.edge_in_between(8, 14));
    lat.flip_spin(lat.edge_in_between(14, 15));
    lat.flip_spin(lat.edge_in_between(15, 16));
    lat.flip_spin(lat.edge_in_between(16, 17));
    BOOST_CHECK(1 == lat.is_percolating());
}

BOOST_AUTO_TEST_CASE(is_percolating_test_3) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    BOOST_CHECK(0 == lat.is_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_percolating_test_1) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        8, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 11));
    lat.flip_spin(lat.edge_in_between(11, 10));
    lat.flip_spin(lat.edge_in_between(10, 18));
    lat.flip_spin(lat.edge_in_between(18, 19));
    lat.flip_spin(lat.edge_in_between(19, 20));
    lat.flip_spin(lat.edge_in_between(20, 21));
    lat.flip_spin(lat.edge_in_between(21, 22));
    lat.flip_spin(lat.edge_in_between(22, 23));
    lat.flip_spin(lat.edge_in_between(23, 16));

    BOOST_CHECK(0 == lat.is_winding_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_percolating_test_2) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        8, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 11));
    lat.flip_spin(lat.edge_in_between(11, 10));
    lat.flip_spin(lat.edge_in_between(10, 18));
    lat.flip_spin(lat.edge_in_between(18, 19));
    lat.flip_spin(lat.edge_in_between(19, 20));
    lat.flip_spin(lat.edge_in_between(20, 21));
    lat.flip_spin(lat.edge_in_between(21, 22));
    lat.flip_spin(lat.edge_in_between(22, 23));
    lat.flip_spin(lat.edge_in_between(23, 16));
    lat.flip_spin(lat.edge_in_between(16, 8));
    lat.flip_spin(lat.edge_in_between(8, 0));

    BOOST_CHECK(1 == lat.is_winding_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_percolating_test_3) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 4));
    lat.flip_spin(lat.edge_in_between(4, 5));
    lat.flip_spin(lat.edge_in_between(5, 11));
    lat.flip_spin(lat.edge_in_between(11, 10));
    lat.flip_spin(lat.edge_in_between(10, 9));
    lat.flip_spin(lat.edge_in_between(9, 8));
    lat.flip_spin(lat.edge_in_between(8, 7));
    lat.flip_spin(lat.edge_in_between(7, 6));
    lat.flip_spin(lat.edge_in_between(6, 0));

    BOOST_CHECK(0 == lat.is_winding_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_percolating_test_4) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 4));
    lat.flip_spin(lat.edge_in_between(4, 5));
    lat.flip_spin(lat.edge_in_between(5, 11));
    lat.flip_spin(lat.edge_in_between(11, 10));
    lat.flip_spin(lat.edge_in_between(10, 9));
    lat.flip_spin(lat.edge_in_between(9, 8));
    lat.flip_spin(lat.edge_in_between(8, 7));
    lat.flip_spin(lat.edge_in_between(7, 6));
    lat.flip_spin(lat.edge_in_between(6, 0));
    lat.flip_spin(lat.edge_in_between(5, 0));

    BOOST_CHECK(1 == lat.is_winding_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_percolating_test_5) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 4));
    lat.flip_spin(lat.edge_in_between(4, 5));
    lat.flip_spin(lat.edge_in_between(5, 11));
    lat.flip_spin(lat.edge_in_between(11, 17));
    lat.flip_spin(lat.edge_in_between(17, 23));
    lat.flip_spin(lat.edge_in_between(23, 29));
    lat.flip_spin(lat.edge_in_between(29, 35));
    lat.flip_spin(lat.edge_in_between(35, 34));
    lat.flip_spin(lat.edge_in_between(34, 33));
    lat.flip_spin(lat.edge_in_between(33, 32));
    lat.flip_spin(lat.edge_in_between(32, 31));
    lat.flip_spin(lat.edge_in_between(31, 30));
    lat.flip_spin(lat.edge_in_between(30, 24));
    lat.flip_spin(lat.edge_in_between(24, 18));
    lat.flip_spin(lat.edge_in_between(18, 12));
    lat.flip_spin(lat.edge_in_between(12, 6));
    lat.flip_spin(lat.edge_in_between(6, 0));

    BOOST_CHECK(0 == lat.is_winding_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_percolating_test_6) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 6));
    lat.flip_spin(lat.edge_in_between(6, 12));
    lat.flip_spin(lat.edge_in_between(12, 18));
    lat.flip_spin(lat.edge_in_between(18, 24));
    lat.flip_spin(lat.edge_in_between(24, 30));
    lat.flip_spin(lat.edge_in_between(30, 0));

    BOOST_CHECK(1 == lat.is_winding_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_percolating_test_7) {
     LatSpec spec = {
        'x', // basis
        "cubic", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 36));
    lat.flip_spin(lat.edge_in_between(36, 72));
    lat.flip_spin(lat.edge_in_between(72, 108));
    lat.flip_spin(lat.edge_in_between(108, 144));
    lat.flip_spin(lat.edge_in_between(144, 180));
    lat.flip_spin(lat.edge_in_between(180, 0));

    BOOST_CHECK(1 == lat.is_winding_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_plaquette_percolating_test_1) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        8, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 8));
    lat.flip_spin(lat.edge_in_between(1, 9));
    lat.flip_spin(lat.edge_in_between(2, 10));
    lat.flip_spin(lat.edge_in_between(3, 11));
    lat.flip_spin(lat.edge_in_between(11, 12));
    lat.flip_spin(lat.edge_in_between(12, 20));
    lat.flip_spin(lat.edge_in_between(13, 21));
    lat.flip_spin(lat.edge_in_between(13, 14));
    lat.flip_spin(lat.edge_in_between(6, 14));
    lat.flip_spin(lat.edge_in_between(7, 15));

    BOOST_CHECK(1 == lat.is_winding_plaquette_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_plaquette_percolating_test_2) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        8, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 8));
    lat.flip_spin(lat.edge_in_between(1, 9));
    lat.flip_spin(lat.edge_in_between(3, 11));
    lat.flip_spin(lat.edge_in_between(11, 12));
    lat.flip_spin(lat.edge_in_between(12, 20));
    lat.flip_spin(lat.edge_in_between(13, 21));
    lat.flip_spin(lat.edge_in_between(13, 14));
    lat.flip_spin(lat.edge_in_between(6, 14));
    lat.flip_spin(lat.edge_in_between(7, 15));

    BOOST_CHECK(0 == lat.is_winding_plaquette_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_plaquette_percolating_test_3) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        8, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(10, 11));
    lat.flip_spin(lat.edge_in_between(18, 19));
    lat.flip_spin(lat.edge_in_between(26, 27));
    lat.flip_spin(lat.edge_in_between(34, 35));
    lat.flip_spin(lat.edge_in_between(42, 43));
    lat.flip_spin(lat.edge_in_between(50, 51));
    lat.flip_spin(lat.edge_in_between(58, 59));

    BOOST_CHECK(1 == lat.is_winding_plaquette_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_plaquette_percolating_test_4) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        8, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 11));
    lat.flip_spin(lat.edge_in_between(11, 12));
    lat.flip_spin(lat.edge_in_between(12, 20));
    lat.flip_spin(lat.edge_in_between(20, 21));
    lat.flip_spin(lat.edge_in_between(21, 29));
    lat.flip_spin(lat.edge_in_between(29, 30));
    lat.flip_spin(lat.edge_in_between(30, 38));
    lat.flip_spin(lat.edge_in_between(38, 39));
    lat.flip_spin(lat.edge_in_between(39, 47));
    lat.flip_spin(lat.edge_in_between(47, 40));
    lat.flip_spin(lat.edge_in_between(40, 48));
    lat.flip_spin(lat.edge_in_between(48, 49));
    lat.flip_spin(lat.edge_in_between(49, 57));
    lat.flip_spin(lat.edge_in_between(57, 58));
    lat.flip_spin(lat.edge_in_between(58, 2));

    BOOST_CHECK(1 == lat.is_winding_plaquette_percolating());
}

BOOST_AUTO_TEST_CASE(is_winding_plaquette_percolating_test_5) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        8, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 11));
    lat.flip_spin(lat.edge_in_between(11, 12));
    lat.flip_spin(lat.edge_in_between(12, 20));
    lat.flip_spin(lat.edge_in_between(20, 21));
    lat.flip_spin(lat.edge_in_between(21, 29));
    lat.flip_spin(lat.edge_in_between(29, 30));
    lat.flip_spin(lat.edge_in_between(30, 38));
    lat.flip_spin(lat.edge_in_between(38, 39));
    lat.flip_spin(lat.edge_in_between(39, 47));
    lat.flip_spin(lat.edge_in_between(47, 40));
    lat.flip_spin(lat.edge_in_between(48, 49));
    lat.flip_spin(lat.edge_in_between(49, 57));
    lat.flip_spin(lat.edge_in_between(57, 58));
    lat.flip_spin(lat.edge_in_between(58, 2));

    BOOST_CHECK(0 == lat.is_winding_plaquette_percolating());
} 

BOOST_AUTO_TEST_CASE(is_winding_cube_percolating_test_1) {
    LatSpec spec = {
        'x', // basis
        "cubic", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 6));
    lat.flip_spin(lat.edge_in_between(1, 7));
    lat.flip_spin(lat.edge_in_between(2, 8));
    lat.flip_spin(lat.edge_in_between(3, 9));
    lat.flip_spin(lat.edge_in_between(4, 10));
    lat.flip_spin(lat.edge_in_between(5, 11));

    BOOST_CHECK(1 == lat.is_winding_cube_percolating());
}

BOOST_AUTO_TEST_CASE(percolation_strength_test_1) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "open", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 1));
    lat.flip_spin(lat.edge_in_between(1, 2));
    lat.flip_spin(lat.edge_in_between(2, 3));
    lat.flip_spin(lat.edge_in_between(3, 9));
    lat.flip_spin(lat.edge_in_between(9, 8));
    lat.flip_spin(lat.edge_in_between(8, 14));
    lat.flip_spin(lat.edge_in_between(14, 15));
    lat.flip_spin(lat.edge_in_between(15, 16));
    lat.flip_spin(lat.edge_in_between(16, 17));

    lat.flip_spin(lat.edge_in_between(18, 19));
    lat.flip_spin(lat.edge_in_between(19, 25));

    BOOST_CHECK(9/(double)lat.get_edge_count() == lat.percolation_strength());
}

BOOST_AUTO_TEST_CASE(plaquette_percolation_strength_test_1) {
     LatSpec spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        6., // beta
        "periodic", // boundaries
        1 // default spin
    };
    auto lat = Lattice(spec);

    lat.flip_spin(lat.edge_in_between(0, 6));
    lat.flip_spin(lat.edge_in_between(1, 7));
    lat.flip_spin(lat.edge_in_between(2, 8));
    lat.flip_spin(lat.edge_in_between(3, 9));
    lat.flip_spin(lat.edge_in_between(4, 10));
    lat.flip_spin(lat.edge_in_between(5, 11));
    lat.flip_spin(lat.edge_in_between(10, 11));

    lat.flip_spin(lat.edge_in_between(13, 14));
    lat.flip_spin(lat.edge_in_between(14, 15));

    BOOST_CHECK(7/(double)lat.get_plaquette_count() == lat.plaquette_percolation_strength());
}

} // namespace paratoric