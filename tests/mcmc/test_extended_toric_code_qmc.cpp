// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#define BOOST_TEST_MODULE TestExtendedToricCodeQMC

#include "mcmc/extended_toric_code_qmc.hpp"
#include "paratoric/types/types.hpp"

#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <vector>
#include <string>

namespace paratoric {

// High field parameter regime outsie of the topological phase
double h = 1.1;  
double h_therm = 0.9;
std::vector<double> h_hys {0.9, 1., 1.1};
double lmbda = 0.8; 
double lmbda_therm = 1.0;
std::vector<double> lmbda_hys {1.0, 1.1, 1.2};
double mu = 0.9; 
double J = 1.1; 

// Low field parameter regime inside of the topological phase
double h_top = 0.1;  
double h_therm_top = 0.0;
std::vector<double> h_hys_top {0.9, 1., 1.1};
double lmbda_top = 0.0; 
double lmbda_therm_top = 0.2;
std::vector<double> lmbda_hys_top {0.05, 0.1, 0.15};
double mu_top = 0.9; 
double J_top = 1.1; 

// Negative J parameter set
double h_neg = 0.1;  
double h_therm_neg = 0.0;
std::vector<double> h_hys_neg {0.9, 1., 1.1};
double lmbda_neg = 0.0; 
double lmbda_therm_neg = 0.2;
std::vector<double> lmbda_hys_neg {0.05, 0.1, 0.15};
double mu_neg = 0.9; 
double J_neg = -1.1; 

BOOST_AUTO_TEST_CASE(get_thermalization_test_1) {
    auto mc = ExtendedToricCodeQMC<'x'>();;
    LatSpec lat_spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda
    };
    SimSpec sim_spec = {
        .N_thermalization = 50,
        .N_resamples = 100,
        .seed = 0,
        .observables = std::vector<std::string> { 
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_thermalization(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_thermalization_test_2) {
    auto mc = ExtendedToricCodeQMC<'x'>();;
    LatSpec lat_spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top
    };
    SimSpec sim_spec = {
        .N_thermalization = 50,
        .N_resamples = 100,
        .seed = 0,
        .observables = std::vector<std::string> { 
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_thermalization(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_thermalization_test_3) {
    auto mc = ExtendedToricCodeQMC<'z'>();;
    LatSpec lat_spec = {
        'z', // basis
        "square", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu_neg,
        .h = h_neg,
        .J = J_neg,
        .lmbda = lmbda_neg
    };
    SimSpec sim_spec = {
        .N_thermalization = 50,
        .N_resamples = 100,
        .seed = 0,
        .observables = std::vector<std::string> { 
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_thermalization(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_1) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        1., // beta --> Test global updates
        "periodic", // boundaries
        1 // default spin
    };
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          
        .lmbda_therm = lmbda_therm   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_2) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "square", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    };
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top,
        .h_therm = h_therm_top,          
        .lmbda_therm = lmbda_therm_top   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_3) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "square", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    };
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          
        .lmbda_therm = lmbda_therm   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = true,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_4) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "square", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    };
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top,
        .h_therm = h_therm_top,          
        .lmbda_therm = lmbda_therm_top   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = true,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_5) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "triangular", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    };
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          
        .lmbda_therm = lmbda_therm   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_6) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "triangular", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    };
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top,
        .h_therm = h_therm_top,          
        .lmbda_therm = lmbda_therm_top   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_7) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "triangular", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    };
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          
        .lmbda_therm = lmbda_therm   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_8) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "triangular", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    };
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top,
        .h_therm = h_therm_top,          
        .lmbda_therm = lmbda_therm_top   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_9) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "honeycomb", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          
        .lmbda_therm = lmbda_therm   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    }; 
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    // Honeycomb lattice does not yet support fredenhagen_marcu in x-basis and should throw std::invalid_argument in Debug mode
    BOOST_WARN_THROW( auto result = mc.get_sample(Config{sim_spec, param_spec, lat_spec, out_spec});,
        std::invalid_argument
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_10) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "honeycomb", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top,
        .h_therm = h_therm_top,          
        .lmbda_therm = lmbda_therm_top   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    }; 
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    // Honeycomb lattice does not yet support fredenhagen_marcu in x-basis and should throw std::invalid_argument in Debug mode
    BOOST_WARN_THROW( auto result = mc.get_sample(Config{sim_spec, param_spec, lat_spec, out_spec});,
        std::invalid_argument
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_11) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "honeycomb", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          
        .lmbda_therm = lmbda_therm   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_12) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "honeycomb", // lattice type
        6, // system size
        10., // beta
        "periodic", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top,
        .h_therm = h_therm_top,          
        .lmbda_therm = lmbda_therm_top   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_13) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "cubic", // lattice type
        6, // system size
        10., // beta
        "open", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          
        .lmbda_therm = lmbda_therm   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_14) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "cubic", // lattice type
        6, // system size
        10., // beta
        "open", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top,
        .h_therm = h_therm_top,          
        .lmbda_therm = lmbda_therm_top   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_15) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "cubic", // lattice type
        6, // system size
        10., // beta
        "open", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          
        .lmbda_therm = lmbda_therm   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_sample_test_16) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "cubic", // lattice type
        6, // system size
        10., // beta
        "open", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu_top,
        .h = h_top,
        .J = J_top,
        .lmbda = lmbda_top,
        .h_therm = h_therm_top,          
        .lmbda_therm = lmbda_therm_top   
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .path_out = std::filesystem::path{"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_sample(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_hysteresis_test_1) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "cubic", // lattice type
        6, // system size
        10., // beta
        "open", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu,
        .J = J,
        .h_hys = h_hys,
        .lmbda_hys = lmbda_hys
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .seed = 0,
        .observables = std::vector<std::string> { 
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .paths_out = std::vector<std::filesystem::path> {"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_hysteresis(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_hysteresis_test_2) {
    auto mc = ExtendedToricCodeQMC<'x'>();
    LatSpec lat_spec = {
        'x', // basis
        "cubic", // lattice type
        6, // system size
        10., // beta
        "open", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu_top,
        .J = J_top,
        .h_hys = h_hys_top,
        .lmbda_hys = lmbda_hys_top
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .seed = 0,
        .observables = std::vector<std::string> { 
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .paths_out = std::vector<std::filesystem::path> {"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_hysteresis(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_hysteresis_test_3) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "cubic", // lattice type
        6, // system size
        10., // beta
        "open", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu,
        .J = J,
        .h_hys = h_hys,
        .lmbda_hys = lmbda_hys
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .seed = 0,
        .observables = std::vector<std::string> { 
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .paths_out = std::vector<std::filesystem::path> {"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_hysteresis(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

BOOST_AUTO_TEST_CASE(get_hysteresis_test_4) {
    auto mc = ExtendedToricCodeQMC<'z'>();
    LatSpec lat_spec = {
        'z', // basis
        "cubic", // lattice type
        6, // system size
        10., // beta
        "open", // boundaries
        1 // default spin
    }; 
    ParamSpec param_spec = {
        .mu = mu_top,
        .J = J_top,
        .h_hys = h_hys_top,
        .lmbda_hys = lmbda_hys_top
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .seed = 0,
        .observables = std::vector<std::string> { 
            "percolation_strength", "percolation_probability", "plaquette_percolation_strength",
            "plaquette_percolation_probability", "largest_cluster", "largest_plaquette_cluster", 
            "string_number", "energy", "energy_h", "energy_mu", "energy_J", "energy_lmbda", 
            "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", 
            "sigma_x_static_susceptibility", "sigma_x_dynamical_susceptibility",
            "sigma_z_static_susceptibility", "sigma_z_dynamical_susceptibility"
        }
    };
    OutSpec out_spec = {
        .paths_out = std::vector<std::filesystem::path> {"./"},
        .save_snapshots = 0
    };

    auto result = mc.get_hysteresis(
        Config{sim_spec, param_spec, lat_spec, out_spec}
    );
}

} // namespace paratoric