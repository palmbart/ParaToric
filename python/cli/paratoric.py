# ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
# Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

import argparse as ap
import logging
from math import pi
import multiprocessing
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

if not logger.handlers:
    formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d - %(module)s - %(levelname)s - %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')

    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(logging.INFO)
    console.flush = sys.stdout.flush
    console.setFormatter(formatter)
    logger.addHandler(console)


def cap_processes(args, limit):
    args['processes'] = min(args['processes'], limit)


def main(args):
    from job_handler import JobHandler

    lattice_params_list = ('lattice_type',
                           'system_size',
                           'boundaries',
                           'default_spin')

    lattice_params = {k: args[k] for k in lattice_params_list}

    if args['simulation'] == 'etc_T_sweep':
        cap_processes(args, args['T_steps'])
        job_handler = JobHandler()
        job_handler.etc_T_sweep(N_samples=args['N_samples'],
                                   N_thermalization=args['N_thermalization'],
                                   N_between_samples=args['N_between_samples'],
                                   T_lower=args['T_lower'],
                                   T_upper=args['T_upper'],
                                   T_steps=args['T_steps'],
                                   mu=args['mu_constant'],
                                   h=args['h_constant'],
                                   h_therm=args['h_constant_therm'],
                                   J=args['J_constant'],
                                   lmbda=args['lmbda_constant'],
                                   lmbda_therm=args['lmbda_constant_therm'],
                                   N_resamples=args['N_resamples'],
                                   custom_therm=args['custom_therm'],
                                   observables=args['observables'],
                                   seed=args['seed'],
                                   basis=args['basis'],
                                   save_snapshots=args['snapshots'],
                                   full_time_series=args['full_time_series'],
                                   processes=args['processes'],
                                   output_dir=args['output_directory'],
                                   **lattice_params)
        
    elif args['simulation'] == 'etc_hysteresis':
        cap_processes(args, 2)
        job_handler = JobHandler()
        job_handler.etc_hysteresis(N_samples=args['N_samples'],
                                      N_thermalization=args['N_thermalization'],
                                      N_between_samples=args['N_between_samples'],
                                      temperature=args['temperature'],
                                      mu=args['mu_constant'],
                                      h_hys=args['h_hysteresis'],
                                      J=args['J_constant'],
                                      lmbda_hys=args['lmbda_hysteresis'],
                                      N_resamples=args['N_resamples'],
                                      observables=args['observables'],
                                      seed=args['seed'],
                                      basis=args['basis'],
                                      save_snapshots=args['snapshots'],
                                      full_time_series=args['full_time_series'],
                                      processes=args['processes'],
                                      output_dir=args['output_directory'],
                                      **lattice_params)

    elif args['simulation'] == 'etc_h_sweep':
        cap_processes(args, args['h_steps'])
        job_handler = JobHandler()
        job_handler.etc_h_sweep(N_samples=args['N_samples'],
                                   N_thermalization=args['N_thermalization'],
                                   N_between_samples=args['N_between_samples'],
                                   temperature=args['temperature'],
                                   mu=args['mu_constant'],
                                   h_lower=args['h_lower'],
                                   h_upper=args['h_upper'],
                                   h_steps=args['h_steps'],
                                   h_therm=args['h_constant_therm'],
                                   J=args['J_constant'],
                                   lmbda=args['lmbda_constant'],
                                   lmbda_therm=args['lmbda_constant_therm'],
                                   N_resamples=args['N_resamples'],
                                   custom_therm=args['custom_therm'],
                                   observables=args['observables'],
                                   seed=args['seed'],
                                   basis=args['basis'],
                                   save_snapshots=args['snapshots'],
                                   full_time_series=args['full_time_series'],
                                   processes=args['processes'],
                                   output_dir=args['output_directory'],
                                   **lattice_params)

    elif args['simulation'] == 'etc_lmbda_sweep':
        cap_processes(args, args['lmbda_steps'])
        job_handler = JobHandler()
        job_handler.etc_lmbda_sweep(N_samples=args['N_samples'],
                                       N_thermalization=args['N_thermalization'],
                                       N_between_samples=args['N_between_samples'],
                                       temperature=args['temperature'],
                                       mu=args['mu_constant'],
                                       h=args['h_constant'],
                                       h_therm=args['h_constant_therm'],
                                       J=args['J_constant'],
                                       lmbda_lower=args['lmbda_lower'],
                                       lmbda_upper=args['lmbda_upper'],
                                       lmbda_steps=args['lmbda_steps'],
                                       lmbda_therm=args['lmbda_constant_therm'],
                                       N_resamples=args['N_resamples'],
                                       custom_therm=args['custom_therm'],
                                       observables=args['observables'],
                                       seed=args['seed'],
                                       basis=args['basis'],
                                       save_snapshots=args['snapshots'],
                                       full_time_series=args['full_time_series'],
                                       processes=args['processes'],
                                       output_dir=args['output_directory'],
                                       **lattice_params)
        
    elif args['simulation'] == 'etc_circle_sweep':
        cap_processes(args, args['Theta_steps'])
        job_handler = JobHandler()
        job_handler.etc_circle_sweep(N_samples=args['N_samples'],
                                        N_thermalization=args['N_thermalization'],
                                        N_between_samples=args['N_between_samples'],
                                        temperature=args['temperature'],
                                        mu=args['mu_constant'],
                                        h=args['h_constant'],
                                        J=args['J_constant'],
                                        lmbda=args['lmbda_constant'],
                                        radius=args['radius'],
                                        Theta_lower=args['Theta_lower'],
                                        Theta_upper=args['Theta_upper'],
                                        Theta_steps=args['Theta_steps'],
                                        N_resamples=args['N_resamples'],
                                        observables=args['observables'],
                                        seed=args['seed'],
                                        basis=args['basis'],
                                        save_snapshots=args['snapshots'],
                                        full_time_series=args['full_time_series'],
                                        processes=args['processes'],
                                        output_dir=args['output_directory'],
                                        **lattice_params)

    elif args['simulation'] == 'etc_thermalization':
        cap_processes(args, args['repetitions'])
        job_handler = JobHandler()
        job_handler.etc_thermalization(N_thermalization=args['N_thermalization'],
                                          repetitions=args['repetitions'],
                                          temperature=args['temperature'],
                                          mu=args['mu_constant'],
                                          h=args['h_constant'],
                                          J=args['J_constant'],
                                          lmbda=args['lmbda_constant'],
                                          N_resamples=args['N_resamples'],
                                          observables=args['observables'],
                                          seed=args['seed'],
                                          basis=args['basis'],
                                          save_snapshots=args['snapshots'],
                                          processes=args['processes'],
                                          output_dir=args['output_directory'],
                                          **lattice_params)
    
    else:
        raise ValueError(f'The simulation "{args["simulation"]}" does not exist.')


if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Run extended toric code QMC simulation', formatter_class=ap.ArgumentDefaultsHelpFormatter)

    # Monte Carlo parameters
    mc_group = parser.add_argument_group('Simulation parameters')
    mc_group.add_argument('-sim', '--simulation',
                          help='Simulation you want to run.',
                          type=str, default='etc_T_sweep', 
                          choices=['etc_T_sweep', 
                                   'etc_hysteresis',
                                   'etc_h_sweep', 
                                   'etc_lmbda_sweep',
                                   'etc_circle_sweep',
                                   'etc_thermalization'])
    mc_group.add_argument('-Ns', '--N_samples', 
                          help='Number of snapshots of the simulation.', 
                          type=int, default=1000)
    mc_group.add_argument('-Nth', '--N_thermalization', 
                          help='Number of thermalization steps.', 
                          type=int, default=2000)
    mc_group.add_argument('-Nbs', '--N_between_samples', 
                          help='Number of steps between samples.', 
                          type=int, default=100)
    mc_group.add_argument('-reps', '--repetitions',
                          help='Number of MC simulations to average over.',
                          type=int, default=10)
    mc_group.add_argument('-T', '--temperature',
                          help='Temperature.',
                          type=float, default=1.)
    mc_group.add_argument('-Tl', '--T_lower',
                          help='Lower limit of temperature.',
                          type=float, default=1.)
    mc_group.add_argument('-Tu', '--T_upper',
                          help='Upper limit of temperature.',
                          type=float, default=2.)
    mc_group.add_argument('-Ts', '--T_steps',
                          help='Number of temperature steps.',
                          type=int, default=15)
    mc_group.add_argument('-muc', '--mu_constant',
                          help='Value of the constant mu (star term).',
                          type=float, default=1.)
    mc_group.add_argument('-Jc', '--J_constant', 
                          help='Value of the constant J (plaquette term).', 
                          type=float, default=1.)
    mc_group.add_argument('-hc', '--h_constant', 
                          help='Value of the constant h (electric field term).', 
                          type=float, default=0.)
    mc_group.add_argument('-hhys', '--h_hysteresis', nargs='+',
                          help='Hysteresis value of the constant h (electric field term).', 
                          type=float, default=[0.])
    mc_group.add_argument('-hct', '--h_constant_therm', 
                          help='Value of the constant h for thermalization (electric field term).', 
                          type=float, default=0.)
    mc_group.add_argument('-hl', '--h_lower',
                          help='Lower limit of h (electric field term).',
                          type=float, default=1.)
    mc_group.add_argument('-hu', '--h_upper',
                          help='Upper limit of h (electric field term).',
                          type=float, default=2.)
    mc_group.add_argument('-hs', '--h_steps',
                          help='Number of h steps (electric field term).',
                          type=int, default=15)
    mc_group.add_argument('-lmbdac', '--lmbda_constant', 
                          help='Value of the constant lmbda (gauge field term).', 
                          type=float, default=0.)
    mc_group.add_argument('-lmbdahys', '--lmbda_hysteresis', nargs='+',
                          help='Hysteresis value of the constant lmbda (gauge field term).', 
                          type=float, default=[0.])
    mc_group.add_argument('-lmbdact', '--lmbda_constant_therm',
                          help='Value of the constant lmbda for thermalization (gauge field term).',
                          type=float, default=0.)
    mc_group.add_argument('-lmbdal', '--lmbda_lower',
                          help='Lower limit of lmbda (gauge field term).',
                          type=float, default=1.)
    mc_group.add_argument('-lmbdau', '--lmbda_upper',
                          help='Upper limit of lmbda (gauge field term).',
                          type=float, default=2.)
    mc_group.add_argument('-lmbdas', '--lmbda_steps',
                          help='Number of lmbda steps (gauge field term).',
                          type=int, default=15)
    mc_group.add_argument('-rad', '--radius',
                          help='Radius of circle sweep in lmbda-h-space.',
                          type=float, default=0.5)
    mc_group.add_argument('-Thl', '--Theta_lower',
                          help='Lower limit of Theta.',
                          type=float, default=0.)
    mc_group.add_argument('-Thu', '--Theta_upper',
                          help='Upper limit of Theta.',
                          type=float, default=2*pi)
    mc_group.add_argument('-Ths', '--Theta_steps',
                          help='Number of Theta steps.',
                          type=int, default=15)
    mc_group.add_argument('-Nr', '--N_resamples', 
                          help='Number of bootstrap resamples.', 
                          type=int, default=1000)
    mc_group.add_argument('-cth', '--custom_therm', 
                          help='Whether custom thermalization is used (to probe hysteresis).', 
                          type=int, default=0, choices=[0, 1])
    mc_group.add_argument('-obs', '--observables', nargs='+',
                          help='Observables that are measured.',
                          type=str, default=['energy'])
    mc_group.add_argument('-s', '--seed', 
                          help='Seed for the pseudorandom number generator.', 
                          type=int, default=0)
    mc_group.add_argument('-bas', '--basis',
                          help='Spin basis (\"x\" or \"z\").',
                          type=str, default='x')
    mc_group.add_argument('-outdir', '--output_directory',
                          help='Directory where the output is stored.',
                          type=str, default=None)
    mc_group.add_argument('-snap', '--snapshots',
                          help='Whether snapshots should be saved.',
                          type=int, default=0, choices=[0, 1])
    mc_group.add_argument('-fts', '--full_time_series',
                          help='Whether full time series should be saved.',
                          type=int, default=0, choices=[0, 1])
    mc_group.add_argument('-proc', '--processes',
                          help='Number of processes. "0" will use all available cores. Negative numbers "-x" will use all available cores minus x.',
                          type=int, default=-4)

    # Lattice parameters
    lattice_group = parser.add_argument_group('Lattice parameters')
    lattice_group.add_argument('-lat', '--lattice_type', 
                               help='Type of lattice used.', 
                               type=str, default='square', choices=['square', 'cubic', 'honeycomb', 'triangular'])
    lattice_group.add_argument('-L', '--system_size',
                               help='System size of lattice (one coordinate).',
                               type=int, default=5)
    lattice_group.add_argument('-bound', '--boundaries', 
                               help='Boundary conditions.', 
                               type=str, default='periodic', choices=['periodic', 'open'])
    lattice_group.add_argument('-dsp', '--default_spin',
                               help='Default spin (electric field) for lattice initialization.',
                               type=int, default=1, choices=[-1, 1])

    args = vars(parser.parse_args())

    args['custom_therm'] = bool(args['custom_therm'])
    args['snapshots'] = bool(args['snapshots'])
    args['full_time_series'] = bool(args['full_time_series'])

    logical_cpu_core_count = multiprocessing.cpu_count()

    if args['processes'] == 0: 
        # Use all available cores when args['processes'] == 0
        args['processes'] = logical_cpu_core_count 
    elif logical_cpu_core_count < args['processes']:
        # Never use more processes than cores
        args['processes'] = logical_cpu_core_count
    elif args['processes'] < 0:
        # Use all available cores minus |args['processes']|
        if (logical_cpu_core_count + args['processes']) > 0:
            args['processes'] = logical_cpu_core_count + args['processes']
        else:
            args['processes'] = 1

    logger.info(f'Potentially use {args["processes"]} out of {logical_cpu_core_count} logical CPU cores.')

    main(args)
