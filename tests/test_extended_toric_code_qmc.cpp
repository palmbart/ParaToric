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

double h = 1.;  
double h_therm = 0.9;
std::vector<double> h_hys {0.9, 1., 1.1};
double lmbda = 1.; 
double lmbda_therm = 1.1;
std::vector<double> lmbda_hys {0.9, 1., 1.1};
double mu = 1.; 
double J = 1.; 

BOOST_AUTO_TEST_CASE(get_sample_test_1) {
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
        .mu = mu,
        .h = h,
        .J = J,
        .lmbda = lmbda,
        .h_therm = h_therm,          // only if you need it
        .lmbda_therm = lmbda_therm   // only if you need it
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility"
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
        .h_therm = h_therm,          // only if you need it
        .lmbda_therm = lmbda_therm   // only if you need it
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = true,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility"
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
        .h_therm = h_therm,          // only if you need it
        .lmbda_therm = lmbda_therm   // only if you need it
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility"
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
        .h_therm = h_therm,          // only if you need it
        .lmbda_therm = lmbda_therm   // only if you need it
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility"
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
        .h_therm = h_therm,          // only if you need it
        .lmbda_therm = lmbda_therm   // only if you need it
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "sigma_x_susceptibility",
            "sigma_z_susceptibility"
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

BOOST_AUTO_TEST_CASE(get_sample_test_6) {
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
        .h_therm = h_therm,          // only if you need it
        .lmbda_therm = lmbda_therm   // only if you need it
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility"
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
        .h_therm = h_therm,          // only if you need it
        .lmbda_therm = lmbda_therm   // only if you need it
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility"
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
        .h_therm = h_therm,          // only if you need it
        .lmbda_therm = lmbda_therm   // only if you need it
    };
    SimSpec sim_spec = {
        .N_samples = 50,
        .N_thermalization = 50,
        .N_between_samples = 50,
        .N_resamples = 100,
        .custom_therm = false,
        .seed = 0,
        .observables = std::vector<std::string>{
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility"
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
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility" 
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
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility" 
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
            "percolation_strength", "percolation_probability", "plaquette_percolation_probability", 
            "largest_cluster", "string_number", "energy", "energy_h", "energy_mu", "energy_J", 
            "energy_lmbda", "sigma_x", "sigma_z", "star_x", "plaquette_z", "staggered_imaginary_times", 
            "delta", "anyon_count", "anyon_density", "fredenhagen_marcu", "sigma_x_susceptibility", 
            "sigma_z_susceptibility"
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

} // namespace paratoric