# ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
# Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

from datetime import datetime
from datetime import timedelta
import h5py
import logging
from itertools import groupby
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
#import numpy.typing as npt
import os
import subprocess
import sys
import time
import traceback
import uuid


class JobHandler:
    def __init__(self):
        self.obs_func_list = None

        # Set up logging
        self.__set_up_logging()

        # Here, all the observables get the right output function and a name for the plots
        self.obs_dict = [{'name': 'percolation_strength',
                         'type': 'real',
                         'output_str': 'Percolation strength',
                         'output_func': self._real_output},
                        {'name': 'percolation_probability',
                         'type': 'real',
                         'output_str': 'Percolation probability',
                         'output_func': self._real_output},
                        {'name': 'plaquette_percolation_strength',
                         'type': 'real',
                         'output_str': 'Plaquette percolation strength',
                         'output_func': self._real_output},
                        {'name': 'plaquette_percolation_probability',
                         'type': 'real',
                         'output_str': 'Plaquette percolation probability',
                         'output_func': self._real_output}, 
                        {'name': 'cube_percolation_strength',
                         'type': 'real',
                         'output_str': 'Cube percolation strength',
                         'output_func': self._real_output},
                        {'name': 'cube_percolation_probability',
                         'type': 'real',
                         'output_str': 'Cube percolation probability',
                         'output_func': self._real_output}, 
                        {'name': 'largest_cluster',
                         'type': 'real',
                         'output_str': 'Largest cluster',
                         'output_func': self._real_output},
                        {'name': 'largest_plaquette_cluster',
                         'type': 'real',
                         'output_str': 'Largest plaquette cluster',
                         'output_func': self._real_output},
                        {'name': 'string_number',
                         'type': 'real',
                         'output_str': 'Number of strings',
                         'output_func': self._real_output},
                        {'name': 'anyon_count',
                         'type': 'real',
                         'output_str': 'Anyon count',
                         'output_func': self._real_output},
                        {'name': 'anyon_density',
                         'type': 'real',
                         'output_str': 'Anyon density',
                         'output_func': self._real_output},
                        {'name': 'fredenhagen_marcu',
                         'type': 'complex',
                         'output_str': 'Fredenhagen-Marcu',
                         'output_func': self._real_output},
                        {'name': 'staggered_imaginary_times',
                         'type': 'real',
                         'output_str': r'$\langle \mu_p \rangle$',
                         'output_func': self._real_output},
                        {'name': 'sigma_x',
                         'type': 'real',
                         'output_str': r'$\langle \sigma^x_l \rangle$',
                         'output_func': self._real_output},
                        {'name': 'sigma_x_susceptibility',
                         'type': 'real',
                         'output_str': r'$\langle \chi_{\sigma^x_l} \rangle$',
                         'output_func': self._real_output},
                        {'name': 'sigma_x_dynamical_susceptibility',
                         'type': 'real',
                         'output_str': r'$\langle \chi_{\sigma^x_l} \rangle$',
                         'output_func': self._real_output},
                        {'name': 'sigma_z',
                         'type': 'real',
                         'output_str': r'$\langle \sigma^z_l \rangle$',
                         'output_func': self._real_output},
                        {'name': 'sigma_z_susceptibility',
                         'type': 'real',
                         'output_str': r'$\langle \chi_{\sigma^z_l} \rangle$',
                         'output_func': self._real_output},
                        {'name': 'sigma_z_dynamical_susceptibility',
                         'type': 'real',
                         'output_str': r'$\langle \chi_{\sigma^z_l} \rangle$',
                         'output_func': self._real_output},
                        {'name': 'star_x',
                         'type': 'real',
                         'output_str': r'$\langle A_s \rangle$',
                         'output_func': self._real_output},
                        {'name': 'plaquette_z',
                         'type': 'real',
                         'output_str': r'$\langle B_p \rangle$',
                         'output_func': self._real_output},
                        {'name': 'delta',
                         'type': 'real',
                         'output_str': r'$\langle A_s \rangle - \langle B_p \rangle$',
                         'output_func': self._real_output},
                        {'name': 'energy',
                         'type': 'real',
                         'output_str': r'$\langle H \rangle$',
                         'output_func': self._real_output},
                        {'name': 'energy_h',
                         'type': 'real',
                         'output_str': r'$\langle -h \; \sum_l \; \sigma^x_l \rangle$',
                         'output_func': self._real_output},
                        {'name': 'energy_mu',
                         'type': 'real',
                         'output_str': r'$\langle -\mu \; \sum_s \; A_s \rangle$',
                         'output_func': self._real_output},
                        {'name': 'energy_lmbda',
                         'type': 'real',
                         'output_str': r'$\langle -\lambda \; \sum_l \; \sigma^z_l \rangle$',
                         'output_func': self._real_output},
                        {'name': 'energy_J',
                         'type': 'real',
                         'output_str': r'$\langle -J \; \sum_p \; B_p \rangle$',
                         'output_func': self._real_output}]
        
    def __set_up_logging(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d - %(module)s - %(levelname)s - %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S')

        console = logging.StreamHandler(stream=sys.stdout)
        console.setLevel(logging.INFO)
        console.flush = sys.stdout.flush
        console.setFormatter(formatter)
        logger.addHandler(console)

        self.log = logger

    def __set_lattice_params(self, **kwargs):
        self.lattice_params = kwargs

    def __get_observable_output_str(self, obs: str):
        for elem in self.obs_dict:
            if elem['name'] == obs:
                return elem['output_str']
        else:
            self.log.info(f"The observable \"{obs}\" has not been found, returning \"{obs}\" as output string!")
            return obs

    def __get_observable_output_func(self, obs: str):
        for elem in self.obs_dict:
            if elem['name'] == obs:
                return elem['output_func']
        else:
            self.log.info(f"The observable \"{obs}\" has not been found, returning self._real_output as output function!")
            return self._real_output
        
    def __get_observable_type(self, obs: str):
        for elem in self.obs_dict:
            if elem['name'] == obs:
                return elem['type']
        else:
            self.log.info(f"The observable \"{obs}\" has not been found, returning \"real\" as output function!")
            return 'real'

    def __dictionary_output(self, dictionary: dict):
        out = ''
        for key, value in dictionary.items():
            out += f"{key}: {value}\n"
        return out

    def __construct_output_directory(self, output_dir: str, name : str, begin_time: str, subpathname: str):
        if output_dir == '':
            if not os.path.exists('../out'):
                os.mkdir('../out')
            output_dir = './' + 'out'

        outname = f"{name}_" + begin_time + "_" + uuid.uuid4().hex
        if not os.path.exists(os.path.join(output_dir, subpathname)):
            os.mkdir(os.path.join(output_dir, subpathname)) # Needed for macOS
        os.mkdir(os.path.join(output_dir, subpathname, outname))
        return os.path.join(output_dir, subpathname, outname)

    def __get_datetime(self):
        return datetime.now().strftime('%d_%m_%Y-%H_%M_%S')

    def __write_hdf5_file(self, dict: dict, path: str, filename: str = 'simulation_data.h5'):
        hf = h5py.File(os.path.join(path, filename), 'w')
        for key, value in dict.items():
            hf.create_dataset(key, data=value, compression='gzip')
        hf.close()

    def __write_parameters_file(self, path: str, run_time: timedelta, begin_time: str, end_time: str, lattice_kwargs: dict, mc_kwargs):
        with open(os.path.join(path, 'parameters.txt'), 'w') as text_file:
            print(
                f"The simulation parameters can be found here."
                f"\n\nBegin: {begin_time}"
                f"\nEnd: {end_time}"
                f"\nSimulation runtime: {run_time}"
                f"\n\nLattice parameters:\n{self.__dictionary_output(lattice_kwargs)}"
                f"\nMonte Carlo parameters:\n{self.__dictionary_output(mc_kwargs)}",
                file=text_file)
            
    def _get_thermalization_cpp(self,
                                 verbose: int,
                                 N_thermalization: int,
                                 beta: float,
                                 mu: float,
                                 h: float,
                                 J: float,
                                 lmbda: float,
                                 N_resamples: int,
                                 output_dir : str,
                                 obs: list,
                                 seed: int,
                                 basis: str = 'x',
                                 save_snapshots: bool = False,
                                 process_index: int = 0):
        lattice_type = self.lattice_params['lattice_type']
        system_size = self.lattice_params['system_size']
        boundaries = self.lattice_params['boundaries']
        if save_snapshots:
            save_snapshots_cpp = 1
        else:
            save_snapshots_cpp = 0
        default_spin = self.lattice_params['default_spin']

        obs_string = ' '.join(obs)

        folder_name = f"lattice={lattice_type}_bounds={boundaries}_basis={basis}_L={system_size}_h={h}_lmbda={lmbda}_mu={mu}_J={J}_beta={beta}_{uuid.uuid4().hex}"
        script_dir = os.path.dirname(__file__)               
        project_root = os.path.dirname(os.path.dirname(script_dir))            
        exe = os.path.join(project_root, "bin", "paratoric")
        bashCommand = f"{exe} --simulation etc_thermalization --N_thermalization {N_thermalization} --beta {beta} --mu_constant {mu} --h_constant {h} --J_constant {J} --lmbda_constant {lmbda} --N_resamples {N_resamples} --observables {obs_string} --seed {seed} " \
                      f"--basis {basis} --lattice_type {lattice_type} --system_size {system_size} --boundaries {boundaries} --default_spin {default_spin} " \
                      f"--output_directory {output_dir} --folder_name {folder_name} --snapshots {save_snapshots_cpp} --process_index {process_index}"
        subprocess.run(bashCommand.split(), check=True)

        # Here we expect that each observable in the hdf5 has the structure 'r', 'i' (real and imaginary part)
        result = []
        with h5py.File(os.path.join(output_dir, folder_name, 'obs.h5'), "r") as f:
            # acc_ratio is just a double and not complex
            acc_ratio = np.array(list(f['simulation/results/acc_ratio']))
            
            for obs_name in obs:
                if self.__get_observable_type(obs_name) in ["real", "complex"]:
                    raw = f[f"simulation/results/{obs_name}/series"][()] # this is a plain float64 array
                    # reinterpret as 16‑byte records with two float64 fields
                    dt = np.dtype([('r','<f8'),('i','<f8')])
                    data_struct = raw.view('V16').view(dt)
                    r = data_struct['r']
                    i = data_struct['i']
                    result.append(r + 1j*i)
        
        result = np.array(result, dtype=np.complex128)

        sample_numbers = [i for i in range(acc_ratio.size)]
        
        return sample_numbers, result, acc_ratio    

    def _get_sample_cpp(self,
                        verbose: int,
                        N_samples: int,
                        N_thermalization: int,
                        N_between_samples: int,
                        beta: float,
                        mu: float,
                        h: float,
                        h_therm: float,
                        J: float,
                        lmbda: float,
                        lmbda_therm: float,
                        N_resamples: int,
                        custom_therm: bool,
                        output_dir : str,
                        obs: list,
                        seed: int,
                        basis: str = 'x',
                        save_snapshots: bool = False,
                        full_time_series: bool = False,
                        process_index: int = 0):
        lattice_type = self.lattice_params['lattice_type']
        system_size = self.lattice_params['system_size']
        boundaries = self.lattice_params['boundaries']
        if save_snapshots:
            save_snapshots_cpp = 1
        else:
            save_snapshots_cpp = 0
        if full_time_series:
            full_time_series_cpp = 1
        else:
            full_time_series_cpp = 0
        default_spin = self.lattice_params['default_spin']

        obs_string = ' '.join(obs)

        folder_name = f"lattice={lattice_type}_bounds={boundaries}_basis={basis}_L={system_size}_h={h}_lmbda={lmbda}_mu={mu}_J={J}_beta={beta}_{uuid.uuid4().hex}"
        script_dir = os.path.dirname(__file__)               
        project_root = os.path.dirname(os.path.dirname(script_dir))             
        exe = os.path.join(project_root, "bin", "paratoric")

        bashCommand = f"{exe} --simulation etc_sample --N_samples {N_samples} --N_thermalization {N_thermalization} --N_between_samples {N_between_samples} --beta {beta} " \
                      f"--mu_constant {mu} --h_constant {h} --h_constant_therm {h_therm} --J_constant {J} --lmbda_constant {lmbda} --lmbda_constant_therm {lmbda_therm} " \
                      f"--N_resamples {N_resamples} --custom_therm {custom_therm} --observables {obs_string} --seed {seed}  " \
                      f"--basis {basis} --lattice_type {lattice_type} --system_size {system_size} --boundaries {boundaries} --default_spin {default_spin} " \
                      f"--output_directory {output_dir} --folder_name {folder_name} --snapshots {save_snapshots_cpp} --full_time_series {full_time_series_cpp} --process_index {process_index}"
        subprocess.run(bashCommand.split(), check=True)

        mean_list = []
        mean_error_list = []
        binder_list = []
        binder_error_list = []
        autocorrelation_time_list = []
        with h5py.File(os.path.join(output_dir, folder_name, 'obs.h5'), "r") as f:
            for obs_name in obs:
                if self.__get_observable_type(obs_name) in ["real", "complex"]:
                    base = f[f"simulation/results/{obs_name}"]

                    mean = base['mean'][()]      
                    mean_error = base['mean_error'][()]
                    binder = base['binder'][()]
                    binder_error = base['binder_error'][()]
                    autocorrelation_time = base['autocorrelation_time'][()]

                    mean_list.append(mean)
                    mean_error_list.append(mean_error)
                    binder_list.append(binder)
                    binder_error_list.append(binder_error)
                    autocorrelation_time_list.append(autocorrelation_time)
        
        mean_list = np.array(mean_list, dtype=np.float64)
        mean_error_list = np.array(mean_error_list, dtype=np.float64)
        binder_list = np.array(binder_list, dtype=np.float64)
        binder_error_list = np.array(binder_error_list, dtype=np.float64)
        autocorrelation_time_list = np.array(autocorrelation_time_list, dtype=np.float64)
        
        return mean_list, mean_error_list, binder_list, binder_error_list, autocorrelation_time_list
    
    def _get_hysteresis_cpp(self,
                            verbose: int,
                            N_samples: int,
                            N_thermalization: int,
                            N_between_samples: int,
                            beta: float,
                            mu: float,
                            h_hys: float,
                            J: float,
                            lmbda_hys: float,
                            N_resamples: int,
                            output_dir : str,
                            obs: list,
                            seed: int,
                            basis: str = 'x',
                            save_snapshots: bool = False,
                            full_time_series: bool = False,
                            process_index: int = 0):
        lattice_type = self.lattice_params['lattice_type']
        system_size = self.lattice_params['system_size']
        boundaries = self.lattice_params['boundaries']
        if save_snapshots:
            save_snapshots_cpp = 1
        else:
            save_snapshots_cpp = 0
        if full_time_series:
            full_time_series_cpp = 1
        else:
            full_time_series_cpp = 0
        default_spin = self.lattice_params['default_spin']

        obs_string = ' '.join(obs)
        h_hys_string = ' '.join([str(i) for i in h_hys])
        lmbda_hys_string = ' '.join([str(i) for i in lmbda_hys])

        folder_names = [f"lattice={lattice_type}_bounds={boundaries}_basis={basis}_L={system_size}_h={h}_lmbda={lmbda}_mu={mu}_J={J}_beta={beta}_{uuid.uuid4().hex}" for (h, lmbda) in zip(h_hys, lmbda_hys)]
        folder_names_string = ' '.join(folder_names)
        script_dir = os.path.dirname(__file__)               
        project_root = os.path.dirname(os.path.dirname(script_dir))               
        exe = os.path.join(project_root, "bin", "paratoric")
        bashCommand = f"{exe} --simulation etc_hysteresis --N_samples {N_samples} --N_thermalization {N_thermalization} --N_between_samples {N_between_samples} " \
                      f"--beta {beta} --mu_constant {mu} --h_hysteresis {h_hys_string} --J_constant {J} --lmbda_hysteresis {lmbda_hys_string} --N_resamples {N_resamples} --observables {obs_string} --seed {seed} " \
                      f"--basis {basis} --lattice_type {lattice_type} --system_size {system_size} --boundaries {boundaries} --default_spin {default_spin} " \
                      f"--output_directory {output_dir} --folder_names {folder_names_string} --snapshots {save_snapshots_cpp} --full_time_series {full_time_series_cpp} --process_index {process_index}"
        subprocess.run(bashCommand.split(), check=True)

        # Here we expect that each number in the hdf5 has the structure 'r', 'i' (real and imaginary part)
        mean_result_array = np.empty(shape=[len(h_hys), len(obs)], dtype=np.float64)
        mean_error_result_array = np.empty(shape=[len(h_hys), len(obs)], dtype=np.float64)
        binder_result_array = np.empty(shape=[len(h_hys), len(obs)], dtype=np.float64)
        binder_error_result_array = np.empty(shape=[len(h_hys), len(obs)], dtype=np.float64)
        autocorrelation_time_result_array = np.empty(shape=[len(h_hys), len(obs)], dtype=np.float64)
        for idx, folder_name in enumerate(folder_names):
            mean_list = []
            mean_error_list = []
            binder_list = []
            binder_error_list = []
            autocorrelation_time_list = []
            
            with h5py.File(os.path.join(output_dir, folder_name, 'obs.h5'), "r") as f:
                for obs_name in obs:
                    if self.__get_observable_type(obs_name) in ["real", "complex"]:
                        base = f[f"simulation/results/{obs_name}"]

                        mean = base['mean'][()]      
                        mean_error = base['mean_error'][()]
                        binder = base['binder'][()]
                        binder_error = base['binder_error'][()]
                        autocorrelation_time = base['autocorrelation_time'][()]

                        mean_list.append(mean)
                        mean_error_list.append(mean_error)
                        binder_list.append(binder)
                        binder_error_list.append(binder_error)
                        autocorrelation_time_list.append(autocorrelation_time)
            
            mean_result_array[idx,:] = np.array(mean_list, dtype=np.float64)
            mean_error_result_array[idx,:] = np.array(mean_error_list, dtype=np.float64)
            binder_result_array[idx,:] = np.array(binder_list, dtype=np.float64)
            binder_error_result_array[idx,:] = np.array(binder_error_list, dtype=np.float64)
            autocorrelation_time_result_array[idx,:] = np.array(autocorrelation_time_list, dtype=np.float64)
        
        return mean_result_array, mean_error_result_array, binder_result_array, binder_error_result_array, autocorrelation_time_result_array

    def _real_output(self, 
                     obs: str, 
                     path: str, 
                     simulation: str, 
                     hdf5_dict: dict, 
                     variable, 
                     obs_array, 
                     obs_array_error, 
                     obs_binder, 
                     obs_binder_error, 
                     obs_tau_int,
                     temperature: float, 
                     h: float, 
                     mu: float, 
                     J: float, 
                     lmbda: float, 
                     radius: float, 
                     comment: str = ''):

        output_str = self.__get_observable_output_str(obs)

        fig, ax = plt.subplots()
        ax.xaxis.set_tick_params(direction='in', which='both')
        ax.yaxis.set_tick_params(direction='in', which='both')
        if simulation == 'etc_T_sweep':
            ax.set(title=f"$\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
            ax.set(xlabel=f"$T$")
            ax.plot(variable, obs_tau_int, 'o-', color='goldenrod')
        elif simulation == 'etc_h_hysteresis':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$")
            ax.set(xlabel=f"$h$")
            ax.plot(variable, obs_tau_int, 'o-', color='goldenrod')
        elif simulation == 'etc_lmbda_hysteresis':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$")
            ax.set(xlabel=f"$\\lambda$")
            ax.plot(variable, obs_tau_int, 'o-', color='goldenrod')
        elif simulation == 'etc_h_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $\\lambda = {lmbda}$")
            ax.set(xlabel=f"$h$")
            ax.plot(variable, obs_tau_int, 'o-', color='goldenrod')
        elif simulation == 'etc_lmbda_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$")
            ax.set(xlabel=f"$\\lambda$")
            ax.plot(variable, obs_tau_int, 'o-', color='goldenrod')
        elif simulation == 'etc_circle_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$, $r = {radius}$")
            ax.set(xlabel=f"$\\theta / \\pi$")
            ax.plot(variable/np.pi, obs_tau_int, 'o-', color='goldenrod')
        else:
            ax.plot(variable, obs_tau_int, 'o-', color='goldenrod')
       
        ax.grid(color='silver', linestyle='-', alpha=0.3)
        ax.set(ylabel=output_str + ' int. autocorr. time')
        fig.tight_layout()
        fig.savefig(os.path.join(path, f"{simulation}_{obs}_{comment}_int_ac_time.pdf"))
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.xaxis.set_tick_params(direction='in', which='both')
        ax.yaxis.set_tick_params(direction='in', which='both')
        if simulation == 'etc_T_sweep':
            ax.set(title=f"$\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
            ax.set(xlabel=f"$T$")
            ax.errorbar(variable, obs_array, yerr=obs_array_error, fmt='o-', color='firebrick', capsize=3)
        elif simulation == 'etc_h_hysteresis':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$")
            ax.set(xlabel=f"$h$")
            ax.errorbar(variable, obs_array, yerr=obs_array_error, fmt='o-', color='firebrick', capsize=3)
        elif simulation == 'etc_lmbda_hysteresis':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$")
            ax.set(xlabel=f"$\\lambda$")
            ax.errorbar(variable, obs_array, yerr=obs_array_error, fmt='o-', color='firebrick', capsize=3)
        elif simulation == 'etc_h_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $\\lambda = {lmbda}$")
            ax.set(xlabel=f"$h$")
            ax.errorbar(variable, obs_array, yerr=obs_array_error, fmt='o-', color='firebrick', capsize=3)
        elif simulation == 'etc_lmbda_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$")
            ax.set(xlabel=f"$\\lambda$")
            ax.errorbar(variable, obs_array, yerr=obs_array_error, fmt='o-', color='firebrick', capsize=3)
        elif simulation == 'etc_circle_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$, $r = {radius}$")
            ax.set(xlabel=f"$\\theta / \\pi$")
            ax.errorbar(variable/np.pi, obs_array, yerr=obs_array_error, fmt='o-', color='firebrick', capsize=3)
        else:
            ax.errorbar(variable, obs_array, yerr=obs_array_error, fmt='o-', color='firebrick', capsize=3)
       
        ax.grid(color='silver', linestyle='-', alpha=0.3)
        ax.set(ylabel=output_str)
        fig.tight_layout()
        fig.savefig(os.path.join(path, f"{simulation}_{obs}_{comment}.pdf"))
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.xaxis.set_tick_params(direction='in', which='both')
        ax.yaxis.set_tick_params(direction='in', which='both')
        if simulation == 'etc_T_sweep':
            ax.set(title=f"$\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
            ax.set(xlabel=f"$T$")
            ax.errorbar(variable, obs_binder, yerr=obs_binder_error, fmt='o-', color='seagreen', capsize=3) 
        elif simulation == 'etc_h_hysteresis':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$")
            ax.set(xlabel=f"$h$")
            ax.errorbar(variable, obs_binder, yerr=obs_binder_error, fmt='o-', color='seagreen', capsize=3) 
        elif simulation == 'etc_lmbda_hysteresis':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$")
            ax.set(xlabel=f"$\\lambda$")
            ax.errorbar(variable, obs_binder, yerr=obs_binder_error, fmt='o-', color='seagreen', capsize=3) 
        elif simulation == 'etc_h_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $\\lambda = {lmbda}$")
            ax.set(xlabel=f"$h$")
            ax.errorbar(variable, obs_binder, yerr=obs_binder_error, fmt='o-', color='seagreen', capsize=3) 
        elif simulation == 'etc_lmbda_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$")
            ax.set(xlabel=f"$\\lambda$")
            ax.errorbar(variable, obs_binder, yerr=obs_binder_error, fmt='o-', color='seagreen', capsize=3) 
        elif simulation == 'etc_circle_sweep':
            ax.set(title=f"$T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$, $r = {radius}$")
            ax.set(xlabel=f"$\\theta / \\pi$")
            ax.errorbar(variable/np.pi, obs_binder, yerr=obs_binder_error, fmt='o-', color='seagreen', capsize=3) 
        else:
            ax.errorbar(variable/np.pi, obs_binder, yerr=obs_binder_error, fmt='o-', color='seagreen', capsize=3) 

        ax.set(ylabel=output_str + ' Binder ratio')
        fig.tight_layout()
        fig.savefig(os.path.join(path, f"{simulation}_{obs}_{comment}_binder.pdf"))
        plt.close(fig)

        hdf5_dict[obs + f"_{comment}"] = obs_array.astype(np.float64)
        hdf5_dict[obs + f"_{comment}" + '_binder'] = obs_binder.astype(np.float64)
        hdf5_dict[obs + f"_{comment}" + '_error'] = obs_array_error.astype(np.float64)
        hdf5_dict[obs + f"_{comment}" + '_binder_error'] = obs_binder_error.astype(np.float64)
        hdf5_dict[obs + f"_{comment}" + '_int_ac_time'] = obs_tau_int.astype(np.float64)

    def etc_T_sweep(self,
                       N_samples: int,
                       N_thermalization: int,
                       N_between_samples: int,
                       T_low: float, T_high: float, T_steps: int,
                       mu: float,
                       h: float,
                       h_therm: float,
                       J: float,
                       lmbda: float,
                       lmbda_therm: float,
                       N_resamples: int,
                       custom_therm: bool,
                       observables: str = 'energy',
                       seed: int = 0,
                       basis: str = 'x',
                       save_snapshots: bool = False,
                       full_time_series: bool = False,
                       processes: int = 8,
                       output_dir: str = '',
                       **kwargs):

        begin_time = self.__get_datetime()

        self.__set_lattice_params(**kwargs)

        # only get name dont change...
        output_dir = self.__construct_output_directory(output_dir, f"lattice={self.lattice_params['lattice_type']}_bounds={self.lattice_params['boundaries']}_basis={basis}_L={self.lattice_params['system_size']}_h={h}_lmbda={lmbda}_mu={mu}_J={J}_T={T_low}_to_{T_high}", begin_time, 'etc_T_sweep')

        temperatures = np.linspace(T_low, T_high, T_steps)
        betas = 1 / temperatures

        start = time.perf_counter()

        self.log.info(f"Spawning {processes} workers for {len(betas)} subprocesses.")
        output_dir_t = os.path.join(output_dir, 'data')
        os.mkdir(output_dir_t)

        if T_steps > 1:
            with mp.Pool(processes=processes) as pool:
                arguments = [(0, N_samples, N_thermalization, N_between_samples, beta, mu, h, h_therm, J, lmbda, lmbda_therm, N_resamples, custom_therm, output_dir_t, observables, seed, basis, save_snapshots, full_time_series, i) for i, beta in enumerate(betas)]
                results = pool.starmap(self._get_sample_cpp, arguments, chunksize=1)

            results = np.array(results)
        else:
            results = self._get_sample_cpp(0, N_samples, N_thermalization, N_between_samples, 1/T_low, mu, h, h_therm, J, lmbda, lmbda_therm, N_resamples, custom_therm, output_dir_t, observables, seed, basis, save_snapshots, full_time_series, 0)

            results = np.array([results])

        end = time.perf_counter()
        run_time = timedelta(seconds=round(end - start, 0))
        self.log.info(f"Workers finished in {run_time} h.")

        end_time = self.__get_datetime()

        ###############
        #   Output    #
        ###############

        sweep_params = {'N_samples': N_samples,
                        'N_thermalization': N_thermalization,
                        'N_between_samples': N_between_samples,
                        'T_low': T_low, 'T_high': T_high, 'T_steps': T_steps,
                        'h': h, 'h_therm': h_therm,
                        'lmbda': lmbda, 'lmbda_therm': lmbda_therm,
                        'mu': mu,
                        'J': J,
                        'N_resamples': N_resamples,
                        'observables': observables,
                        'integrated_autocorrelation_times': results[:, 4, :].T,
                        'basis': basis}

        self.log.info('Writing parameters file...')

        self.__write_parameters_file(output_dir, run_time, begin_time, end_time, kwargs, sweep_params)

        hdf5_dict = {'temperatures': temperatures.astype(np.float64)}

        self.log.info('Plotting observables...')

        for i, obs in enumerate(observables):
            try:
                obs_func = self.__get_observable_output_func(obs)
                obs_func(obs, output_dir, 'etc_T_sweep', hdf5_dict, temperatures,
                        results[:, 0, i], results[:, 1, i], results[:, 2, i], results[:, 3, i], results[:, 4, i], 0, h, mu, J, lmbda, 0)
            except Exception as e:
                self.log.info(traceback.format_exc())
            
        self.log.info('Writing HDF5 file...')

        self.__write_hdf5_file(hdf5_dict, output_dir)

    def etc_hysteresis(self,
                          N_samples: int,
                          N_thermalization: int,
                          N_between_samples: int,
                          temperature: float,
                          mu: float,
                          h_hys: float,
                          J: float,
                          lmbda_hys: float,
                          N_resamples: int,
                          observables: str = 'energy',
                          seed: int = 0,
                          basis: str = 'x',
                          save_snapshots: bool = False,
                          full_time_series: bool = False,
                          processes: int = 8,
                          output_dir: str = '',
                          **kwargs):

        begin_time = self.__get_datetime()

        self.__set_lattice_params(**kwargs)

        # only get name dont change...
        output_dir = self.__construct_output_directory(output_dir, f"lattice={self.lattice_params['lattice_type']}_bounds={self.lattice_params['boundaries']}_basis={basis}_L={self.lattice_params['system_size']}_h={h_hys[0]}_to_{h_hys[-1]}_lmbda={lmbda_hys[0]}_to_{lmbda_hys[-1]}_mu={mu}_J={J}_T={temperature}", begin_time, 'etc_hysteresis')

        beta = 1 / temperature

        start = time.perf_counter()

        self.log.info(f"Spawning {processes} workers for 2 subprocesses.")
        output_dir_t = os.path.join(output_dir, 'data')
        os.mkdir(output_dir_t)
        with mp.Pool(processes=processes) as pool:
            arguments = [(0, N_samples, N_thermalization, N_between_samples, beta, mu, h_hys, J, lmbda_hys, N_resamples, output_dir_t, observables, seed, basis, save_snapshots, 0),
                         (0, N_samples, N_thermalization, N_between_samples, beta, mu, h_hys[::-1], J, lmbda_hys[::-1], N_resamples, output_dir_t, observables, seed, basis, save_snapshots, 1)]
            [results_forward, results_backward] = pool.starmap(self._get_hysteresis_cpp, arguments, chunksize=1)

        results_forward = np.array(results_forward)
        results_backward = np.array(results_backward)

        end = time.perf_counter()
        run_time = timedelta(seconds=round(end - start, 0))
        self.log.info(f"Workers finished in {run_time} h.")

        end_time = self.__get_datetime()

        ###############
        #   Output    #
        ###############

        def all_equal(iterable):
            g = groupby(iterable)
            return next(g, True) and not next(g, False)

        sweep_params = {'N_samples': N_samples,
                        'N_thermalization': N_thermalization,
                        'N_between_samples': N_between_samples,
                        'temperature': temperature,
                        'h_hys': h_hys,
                        'lmbda_hys': lmbda_hys,
                        'mu': mu,
                        'J': J,
                        'N_resamples': N_resamples,
                        'observables': observables,
                        'integrated_autocorrelation_times_forward': results_forward[4, :, :].T,
                        'integrated_autocorrelation_times_backward': results_backward[4, :, :].T,
                        'basis': basis}

        self.log.info('Writing parameters file...')

        self.__write_parameters_file(output_dir, run_time, begin_time, end_time, kwargs, sweep_params)

        hdf5_dict = {'h_hys': np.array(h_hys).astype(np.float64), 'lmbda_hys': np.array(lmbda_hys).astype(np.float64)}

        self.log.info('Plotting observables...')

        for i, obs in enumerate(observables):
            try:
                if all_equal(h_hys):
                    obs_func = self.__get_observable_output_func(obs)
                    obs_func(obs, output_dir, 'etc_lmbda_hysteresis', hdf5_dict, lmbda_hys,
                            results_forward[0, :, i], results_forward[1, :, i], results_forward[2, :, i], results_forward[3, :, i], results_forward[4, :, i], temperature, h_hys[0], mu, J, 0, 0, 'forward')
                    obs_func(obs, output_dir, 'etc_lmbda_hysteresis', hdf5_dict, lmbda_hys[::-1],
                            results_backward[0, :, i], results_backward[1, :, i], results_backward[2, :, i], results_backward[3, :, i], results_backward[4, :, i], temperature, h_hys[0], mu, J, 0, 0, 'backward')
                else:
                    obs_func = self.__get_observable_output_func(obs)
                    obs_func(obs, output_dir, 'etc_h_hysteresis', hdf5_dict, h_hys,
                            results_forward[0, :, i], results_forward[1, :, i], results_forward[2, :, i], results_forward[3, :, i], results_forward[4, :, i], temperature, 0, mu, J, lmbda_hys[0], 0, 'forward')
                    obs_func(obs, output_dir, 'etc_h_hysteresis', hdf5_dict, h_hys[::-1],
                            results_backward[0, :, i], results_backward[1, :, i], results_backward[2, :, i], results_backward[3, :, i], results_backward[4, :, i], temperature, 0, mu, J, lmbda_hys[0], 0, 'backward')
            except Exception as e:
                self.log.info(traceback.format_exc())
            
        self.log.info('Writing HDF5 file...')

        self.__write_hdf5_file(hdf5_dict, output_dir)

    def etc_h_sweep(self,
                       N_samples: int,
                       N_thermalization: int,
                       N_between_samples: int,
                       temperature: float,
                       mu: float,
                       h_low: float, h_high: float, h_steps: int,
                       h_therm: float,
                       J: float,
                       lmbda: float,
                       lmbda_therm: float,
                       N_resamples: int,
                       custom_therm: bool,
                       observables: str = 'energy',
                       seed: int = 0,
                       basis: str = 'x',
                       save_snapshots: bool = False,
                       full_time_series: bool = False,
                       processes: int = 8,
                       output_dir: str = '',
                       **kwargs):

        begin_time = self.__get_datetime()

        self.__set_lattice_params(**kwargs)

        output_dir = self.__construct_output_directory(output_dir, f"lattice={self.lattice_params['lattice_type']}_bounds={self.lattice_params['boundaries']}_basis={basis}_L={self.lattice_params['system_size']}_h={h_low}_to_{h_high}_lmbda={lmbda}_mu={mu}_J={J}_T={temperature}", begin_time, 'etc_h_sweep')

        hs = np.linspace(h_low, h_high, h_steps)
        beta = 1 / temperature

        start = time.perf_counter()

        self.log.info(f"Spawning {processes} workers for {len(hs)} subprocesses.")
        output_dir_t = os.path.join(output_dir, 'data')
        os.mkdir(output_dir_t)

        if h_steps > 1:
            with mp.Pool(processes=processes) as pool:
                arguments = [(0, N_samples, N_thermalization, N_between_samples, beta, mu, h, h_therm, J, lmbda, lmbda_therm, N_resamples, custom_therm, output_dir_t, observables, seed, basis, save_snapshots, full_time_series, i) for i, h in enumerate(hs)]
                results = pool.starmap(self._get_sample_cpp, arguments, chunksize=1)

            results = np.array(results)
        else:
            results = self._get_sample_cpp(0, N_samples, N_thermalization, N_between_samples, beta, mu, h_low, h_therm, J, lmbda, lmbda_therm, N_resamples, custom_therm, output_dir_t, observables, seed, basis, save_snapshots, full_time_series, 0)

            results = np.array([results])

        end = time.perf_counter()
        run_time = timedelta(seconds=round(end - start, 0))
        self.log.info(f"Workers finished in {run_time} h.")

        end_time = self.__get_datetime()

        ###############
        #   Output    #
        ###############

        sweep_params = {'N_samples': N_samples,
                        'N_thermalization': N_thermalization,
                        'N_between_samples': N_between_samples,
                        'temperature': temperature, 
                        'h_low': h_low, 'h_high': h_high, 'h_steps': h_steps, 'h_therm': h_therm,
                        'lmbda': lmbda, 'lmbda_therm': lmbda_therm,
                        'mu': mu,
                        'J': J,
                        'N_resamples': N_resamples,
                        'observables': observables,
                        'integrated_autocorrelation_times': results[:, 4, :].T,
                        'basis': basis}
        
        self.log.info('Writing parameters file...')

        self.__write_parameters_file(output_dir, run_time, begin_time, end_time, kwargs, sweep_params)

        hdf5_dict = {'hs': hs.astype(np.float64)}

        self.log.info('Plotting observables...')

        for i, obs in enumerate(observables):
            try:
                obs_func = self.__get_observable_output_func(obs)
                obs_func(obs, output_dir, 'etc_h_sweep', hdf5_dict, hs,
                        results[:, 0, i], results[:, 1, i], results[:, 2, i], results[:, 3, i], results[:, 4, i], temperature, 0, mu, J, lmbda, 0)
            except Exception as e:
                self.log.info(traceback.format_exc())
            
        self.log.info('Writing HDF5 file...')

        self.__write_hdf5_file(hdf5_dict, output_dir)

    def etc_lmbda_sweep(self,
                           N_samples: int,
                           N_thermalization: int,
                           N_between_samples: int,
                           temperature: float,
                           mu: float,
                           h: float,
                           h_therm: float,
                           J: float,
                           lmbda_low: float, lmbda_high: float, lmbda_steps: int,
                           lmbda_therm: float,
                           N_resamples: int,
                           custom_therm: bool,
                           observables: str = 'energy',
                           seed: int = 0,
                           basis: str = 'x',
                           save_snapshots: bool = False,
                           full_time_series: bool = False,
                           processes: int = 8,
                           output_dir: str = '',
                           **kwargs):

        begin_time = self.__get_datetime()

        self.__set_lattice_params(**kwargs)

        output_dir = self.__construct_output_directory(output_dir, f"lattice={self.lattice_params['lattice_type']}_bounds={self.lattice_params['boundaries']}_basis={basis}_L={self.lattice_params['system_size']}_h={h}_lmbda={lmbda_low}_to_{lmbda_high}_mu={mu}_J={J}_T={temperature}", begin_time, 'etc_lmbda_sweep')

        lmbdas = np.linspace(lmbda_low, lmbda_high, lmbda_steps)
        beta = 1 / temperature

        start = time.perf_counter()

        self.log.info(f"Spawning {processes} workers for {len(lmbdas)} subprocesses.")
        output_dir_t = os.path.join(output_dir, 'data')
        os.mkdir(output_dir_t)

        if lmbda_steps > 1:
            with mp.Pool(processes=processes) as pool:
                arguments = [(0, N_samples, N_thermalization, N_between_samples, beta, mu, h, h_therm, J, lmbda, lmbda_therm, N_resamples, custom_therm, output_dir_t, observables, seed, basis, save_snapshots, full_time_series, i) for i, lmbda in enumerate(lmbdas)]
                results = pool.starmap(self._get_sample_cpp, arguments, chunksize=1)

            results = np.array(results)
        else:
            results = self._get_sample_cpp(0, N_samples, N_thermalization, N_between_samples, beta, mu, h, h_therm, J, lmbda_low, lmbda_therm, N_resamples, custom_therm, output_dir_t, observables, seed, basis, save_snapshots, full_time_series, 0)

            results = np.array([results])

        end = time.perf_counter()
        run_time = timedelta(seconds=round(end - start, 0))
        self.log.info(f"Workers finished in {run_time} h.")

        end_time = self.__get_datetime()

        ###############
        #   Output    #
        ###############

        sweep_params = {'N_samples': N_samples,
                        'N_thermalization': N_thermalization,
                        'N_between_samples': N_between_samples,
                        'temperature': temperature, 
                        'h': h, 'h_therm': h_therm,
                        'lmbda_low': lmbda_low, 'lmbda_high': lmbda_high, 'lmbda_steps': lmbda_steps, 'lmbda_therm': lmbda_therm,
                        'mu': mu,
                        'J': J,
                        'N_resamples': N_resamples,
                        'observables': observables,
                        'integrated_autocorrelation_times': results[:, 4, :].T,
                        'basis': basis}
        
        self.log.info('Writing parameters file...')

        self.__write_parameters_file(output_dir, run_time, begin_time, end_time, kwargs, sweep_params)

        hdf5_dict = {'lmbdas': lmbdas.astype(np.float64)}

        self.log.info('Plotting observables...')

        for i, obs in enumerate(observables):
            try:
                obs_func = self.__get_observable_output_func(obs)
                obs_func(obs, output_dir, 'etc_lmbda_sweep', hdf5_dict, lmbdas,
                        results[:, 0, i], results[:, 1, i], results[:, 2, i], results[:, 3, i], results[:, 4, i], temperature, h, mu, J, 0, 0)
            except Exception as e:
                self.log.info(traceback.format_exc())
            
        self.log.info('Writing HDF5 file...')

        self.__write_hdf5_file(hdf5_dict, output_dir)

    def etc_circle_sweep(self,
                            N_samples: int,
                            N_thermalization: int,
                            N_between_samples: int,
                            temperature: float,
                            mu: float,
                            h: float,
                            J: float,
                            lmbda: float,
                            radius: float,
                            Theta_low: float, Theta_high: float, Theta_steps: int,
                            N_resamples: int,
                            observables: str = 'energy',
                            seed: int = 0,
                            basis: str = 'x',
                            save_snapshots: bool = False,
                            full_time_series: bool = False,
                            processes: int = 8,
                            output_dir: str = '',
                            **kwargs):

        begin_time = self.__get_datetime()

        self.__set_lattice_params(**kwargs)

        output_dir = self.__construct_output_directory(output_dir, f"lattice={self.lattice_params['lattice_type']}_bounds={self.lattice_params['boundaries']}_basis={basis}_L={self.lattice_params['system_size']}_radius={radius}_Theta={Theta_low}_to_{Theta_high}_h={h}_lmbda={lmbda}_mu={mu}_J={J}_T={temperature}", begin_time, 'etc_circle_sweep')

        angles_array = np.array([Theta_low + (Theta_high - Theta_low)/Theta_steps * x for x in range(0, Theta_steps)])

        lambda_h_pair_array = np.array([(lmbda + np.cos(angle)*radius, h + np.sin(angle)*radius) for angle in angles_array])

        beta = 1 / temperature

        start = time.perf_counter()

        self.log.info(f"Spawning {processes} workers for {len(lambda_h_pair_array)} subprocesses.")
        output_dir_t = os.path.join(output_dir, 'data')
        os.mkdir(output_dir_t)
        with mp.Pool(processes=processes) as pool:
            arguments = [(0, N_samples, N_thermalization, N_between_samples, beta, mu, loc_h, loc_h, J, loc_lmbda, loc_lmbda, N_resamples, False, output_dir_t, observables, seed, basis, save_snapshots, full_time_series, i) for i, (loc_lmbda, loc_h) in enumerate(lambda_h_pair_array)]
            results = pool.starmap(self._get_sample_cpp, arguments, chunksize=1)

        results = np.array(results)

        end = time.perf_counter()
        run_time = timedelta(seconds=round(end - start, 0))
        self.log.info(f"Workers finished in {run_time} h.")

        end_time = self.__get_datetime()

        ###############
        #   Output    #
        ###############

        sweep_params = {'N_samples': N_samples,
                        'N_thermalization': N_thermalization,
                        'N_between_samples': N_between_samples,
                        'temperature': temperature, 
                        'h': h,
                        'lmbda': lmbda,
                        'mu': mu,
                        'J': J, 
                        'radius': radius, 
                        'Theta_low': Theta_low, 'Theta_high': Theta_high, 'Theta_steps': Theta_steps,
                        'N_resamples': N_resamples,
                        'observables': observables,
                        'integrated_autocorrelation_times': results[:, 4, :].T,
                        'basis': basis}
        
        self.log.info('Writing parameters file...')

        self.__write_parameters_file(output_dir, run_time, begin_time, end_time, kwargs, sweep_params)

        hdf5_dict = {'angles_array': angles_array.astype(np.float64),
                     'radius': np.array([radius]).astype(np.float64),
                     'lambda_h_pair_array': lambda_h_pair_array.astype(np.float64)}

        self.log.info('Plotting observables...')

        for i, obs in enumerate(observables):
            try:
                obs_func = self.__get_observable_output_func(obs)
                obs_func(obs, output_dir, 'etc_circle_sweep', hdf5_dict, angles_array,
                        results[:, 0, i], results[:, 1, i], results[:, 2, i], results[:, 3, i], results[:, 4, i], temperature, h, mu, J, lmbda, radius)
            except Exception as e:
                self.log.info(traceback.format_exc())
            
        self.log.info('Writing HDF5 file...')

        self.__write_hdf5_file(hdf5_dict, output_dir)

    def etc_thermalization(self,
                              N_thermalization: int,
                              repetitions: int,
                              temperature: float,
                              mu: float,
                              h: float,
                              J: float,
                              lmbda: float,
                              N_resamples: int,
                              observables: str = 'energy',
                              seed: int = 0,
                              basis: str = 'x',
                              save_snapshots: bool = False,
                              processes: int = 8,
                              output_dir: str = '',
                              **kwargs):

        begin_time = self.__get_datetime()

        self.__set_lattice_params(**kwargs)

        output_dir = self.__construct_output_directory(output_dir, f"lattice={self.lattice_params['lattice_type']}_bounds={self.lattice_params['boundaries']}_basis={basis}_L={self.lattice_params['system_size']}_h={h}_lmbda={lmbda}_mu={mu}_J={J}_T={temperature}", begin_time, 'etc_thermalization')

        beta = 1 / temperature

        start = time.perf_counter()

        self.log.info(f"Spawning {processes} workers for {repetitions} subprocesses.")
        output_dir_t = os.path.join(output_dir, 'data')
        os.mkdir(output_dir_t)
        with mp.Pool(processes=processes) as pool:
            arguments = ((0, N_thermalization, beta, mu, h, J, lmbda, N_resamples, output_dir_t, observables, seed, basis, save_snapshots, i) for i in range(repetitions))
            results = pool.starmap(self._get_thermalization_cpp, arguments, chunksize=1)

        step_array = np.asarray(results[0][0], dtype=np.float64)          
        obs_stack = np.stack([r[1] for r in results])   
        acc_ratio_stack = np.stack([r[2] for r in results]) 

        obs_array = obs_stack.mean(axis=0) 
        acc_ratio_array = acc_ratio_stack.mean(axis=0)

        end = time.perf_counter()
        run_time = timedelta(seconds=round(end - start, 0))
        self.log.info(f"Workers finished in {run_time} h.")

        end_time = self.__get_datetime()

        ###############
        #   Output    #
        ###############

        therm_params = {'N_thermalization': N_thermalization,
                        'repetitions': repetitions,
                        'temperature': temperature,
                        'h': h,
                        'lmbda': lmbda,
                        'mu': mu,
                        'J': J,
                        'N_resamples': N_resamples,
                        'observables': observables,
                        'basis': basis}
        
        self.log.info('Writing parameters file...')

        self.__write_parameters_file(output_dir, run_time, begin_time, end_time, kwargs, therm_params)

        hdf5_dict = {'steps': step_array.astype(np.float64), 'acc_ratios': acc_ratio_array.astype(np.float64)}

        self.log.info('Plotting observables...')

        fig, ax = plt.subplots()
        ax.xaxis.set_tick_params(direction='in', which='both')
        ax.yaxis.set_tick_params(direction='in', which='both')
        ax.plot(step_array, acc_ratio_array, color='navy')
        ax.grid(color='silver', linestyle='-', alpha=0.3)
        ax.set(title=f"Thermalization for $T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
        ax.set(yscale='log')
        ax.set(xlabel='Monte Carlo step')
        ax.set(ylabel='Monte Carlo acceptance ratios')

        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, 'etc_thermalization_acc_ratio.pdf'))
        plt.close(fig)

        for i, obs in enumerate(observables):
            try:
                if self.__get_observable_type(obs) == 'real':
                    output_str = self.__get_observable_output_str(obs)

                    fig, ax = plt.subplots()
                    ax.xaxis.set_tick_params(direction='in', which='both')
                    ax.yaxis.set_tick_params(direction='in', which='both')
                    ax.plot(step_array, obs_array[i].real, color='firebrick')
                    ax.grid(color='silver', linestyle='-', alpha=0.3)
                    ax.set(title=f"Thermalization for $T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
                    ax.set(xlabel='Monte Carlo step')
                    ax.set(ylabel=output_str)

                    fig.tight_layout()
                    fig.savefig(os.path.join(output_dir, f"etc_thermalization_{obs}.pdf"))
                    plt.close(fig)

                    hdf5_dict[obs] = obs_array[i].real

                    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
                    ax1.xaxis.set_tick_params(direction='in', which='both')
                    ax1.yaxis.set_tick_params(direction='in', which='both')
                    ax2.xaxis.set_tick_params(direction='in', which='both')
                    ax2.yaxis.set_tick_params(direction='in', which='both')
                    fig.suptitle(f"Thermalization for $T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
                    ax1.set(yscale='log')
                    ax1.set(ylabel='Monte Carlo acceptance ratio')
                    ax1.grid(color='silver', linestyle='-', alpha=0.3)
                    ax1.plot(step_array, acc_ratio_array, color='navy')

                    ax2.set(xlabel='Monte Carlo step')
                    ax2.set(ylabel=output_str)
                    ax2.grid(color='silver', linestyle='-', alpha=0.3)
                    ax2.plot(step_array, obs_array[i].real, color='firebrick')

                    fig.tight_layout()
                    fig.savefig(os.path.join(output_dir, f"etc_thermalization_duplex_{obs}.pdf"))
                    plt.close(fig)
                else:
                    output_str = self.__get_observable_output_str(obs)

                    # Real part

                    fig, ax = plt.subplots()
                    ax.xaxis.set_tick_params(direction='in', which='both')
                    ax.yaxis.set_tick_params(direction='in', which='both')
                    ax.plot(step_array, obs_array[i].real, color='firebrick')
                    ax.grid(color='silver', linestyle='-', alpha=0.3)
                    ax.set(title=f"Thermalization for $T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
                    ax.set(xlabel='Monte Carlo step')
                    ax.set(ylabel=output_str+'_real')

                    fig.tight_layout()
                    fig.savefig(os.path.join(output_dir, f"etc_thermalization_{obs}.pdf"))
                    plt.close(fig)

                    hdf5_dict[obs+'_real'] = obs_array[i].real

                    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
                    ax1.xaxis.set_tick_params(direction='in', which='both')
                    ax1.yaxis.set_tick_params(direction='in', which='both')
                    ax2.xaxis.set_tick_params(direction='in', which='both')
                    ax2.yaxis.set_tick_params(direction='in', which='both')
                    fig.suptitle(f"Thermalization for $T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
                    ax1.set(yscale='log')
                    ax1.set(ylabel='Monte Carlo acceptance ratio')
                    ax1.grid(color='silver', linestyle='-', alpha=0.3)
                    ax1.plot(step_array, acc_ratio_array, color='navy')

                    ax2.set(xlabel='Monte Carlo step')
                    ax2.set(ylabel=output_str+'_real')
                    ax2.grid(color='silver', linestyle='-', alpha=0.3)
                    ax2.plot(step_array, obs_array[i].real, color='firebrick')

                    fig.tight_layout()
                    fig.savefig(os.path.join(output_dir, f"etc_thermalization_duplex_{obs}_real.pdf"))
                    plt.close(fig)

                    # Imaginary part

                    fig, ax = plt.subplots()
                    ax.xaxis.set_tick_params(direction='in', which='both')
                    ax.yaxis.set_tick_params(direction='in', which='both')
                    ax.plot(step_array, obs_array[i].imag, color='firebrick')
                    ax.grid(color='silver', linestyle='-', alpha=0.3)
                    ax.set(title=f"Thermalization for $T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
                    ax.set(xlabel='Monte Carlo step')
                    ax.set(ylabel=output_str+'_imag')

                    fig.tight_layout()
                    fig.savefig(os.path.join(output_dir, f"etc_thermalization_{obs}.pdf"))
                    plt.close(fig)

                    hdf5_dict[obs+'_imag'] = obs_array[i].imag

                    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
                    ax1.xaxis.set_tick_params(direction='in', which='both')
                    ax1.yaxis.set_tick_params(direction='in', which='both')
                    ax2.xaxis.set_tick_params(direction='in', which='both')
                    ax2.yaxis.set_tick_params(direction='in', which='both')
                    fig.suptitle(f"Thermalization for $T = {temperature}$, $\\mu = {mu}$, $J = {J}$, $h = {h}$, $\\lambda = {lmbda}$")
                    ax1.set(yscale='log')
                    ax1.set(ylabel='Monte Carlo acceptance ratio')
                    ax1.grid(color='silver', linestyle='-', alpha=0.3)
                    ax1.plot(step_array, acc_ratio_array, color='navy')

                    ax2.set(xlabel='Monte Carlo step')
                    ax2.set(ylabel=output_str+'_imag')
                    ax2.grid(color='silver', linestyle='-', alpha=0.3)
                    ax2.plot(step_array, obs_array[i].imag, color='firebrick')

                    fig.tight_layout()
                    fig.savefig(os.path.join(output_dir, f"etc_thermalization_duplex_{obs}_imag.pdf"))
                    plt.close(fig)

            except Exception as e:
                self.log.info(traceback.format_exc())
        
        self.log.info('Writing HDF5 file...')

        self.__write_hdf5_file(hdf5_dict, output_dir)
