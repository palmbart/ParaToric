// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#pragma once

#include "mcmc/extended_toric_code_qmc.hpp"

// Frozen pre-optimization transition for seeded differential tests.
// Retain generic gathering and sorting independently of the optimized pair
// fast path. Production remains private; only this test accessor is a friend.
namespace paratoric {
struct TupleParityTestAccess {
    template<char Basis>
    static void step(ExtendedToricCodeQMC<Basis>& mc, Lattice& lat, double& energy,
                     double& ratio, double beta, double h, double mu, double J, double lmbda) {
        mc.metropolis_step_tuple_parity(lat, energy, ratio, beta, h, mu, J, lmbda);
    }

    template<char Basis>
    static double energy(ExtendedToricCodeQMC<Basis>& mc, Lattice& lat,
                         double h, double mu, double J, double lmbda) {
        return mc.total_integrated_pot_energy(lat, h, mu, J, lmbda);
    }

    template<char Basis>
    static void reference(ExtendedToricCodeQMC<Basis>& mc, Lattice& lat,
                          double& integrated_pot_energy, double& acc_ratio, double beta,
                          double h, double mu, double J, double lmbda) {
        std::vector<unsigned char> create_events;
        std::vector<double> event_times, edge_diffs;
        acc_ratio = 0.;
        const double coupling = Basis == 'x' ? J : mu;
        if (coupling == 0.) return;
        const auto& relations = lat.get_tuple_parity_relations();
        if (relations.empty()) return;
        const auto relation_index = mc.random_index(static_cast<int>(relations.size()));
        const auto& relation = relations[relation_index];
        const auto tuple_edges = [&](int i) {
            if constexpr (Basis == 'x') return lat.get_plaquette_edges(i);
            else return lat.get_star_edges(i);
        };

        const detail::PoissonParityProposal proposal(beta, coupling);
        create_events.resize(relation.size());
        double log_ratio = 0.;
        for (std::size_t i = 0; i < relation.size(); ++i) {
            const auto count = lat.get_tuple_spin_flips(relation[i]).size();
            const bool create = count == 0 || mc.uniform_dist(*mc.rng) < proposal.birth_probability(count);
            create_events[i] = create;
            log_ratio += proposal.log_ratio(count, create);
        }

        // Without a diagonal edge field, acceptance needs no times, edge traversal,
        // cache reads or mutations. Generate event locations only if accepted.
        const double diagonal_field = Basis == 'x' ? h : lmbda;
        if (diagonal_field == 0.) {
            acc_ratio = std::exp(log_ratio);
            if (!mc.accept(acc_ratio)) return;
            for (std::size_t i = 0; i < relation.size(); ++i) {
                const int tuple = relation[i];
                const auto edges = tuple_edges(tuple);
                if (create_events[i]) lat.insert_tuple_flip(tuple, edges, mc.uniform_real(0., beta));
                else {
                    const auto times = lat.get_tuple_spin_flips(tuple);
                    const double time = times[mc.random_index(static_cast<int>(times.size()))];
                    lat.delete_tuple_flip(tuple, edges, time);
                }
            }
            return;
        }

        event_times.resize(relation.size());
        for (std::size_t i = 0; i < relation.size(); ++i) {
            if (create_events[i]) event_times[i] = mc.uniform_real(0., beta);
            else {
                const auto times = lat.get_tuple_spin_flips(relation[i]);
                event_times[i] = times[mc.random_index(static_cast<int>(times.size()))];
            }
        }

        const auto& edges = lat.get_tuple_parity_edges(relation_index);
        edge_diffs.resize(edges.size());
        double action_diff = 0.;
        Lattice::SmallEnergyVector times;
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto& affected = edges[i];
            times.clear();
            for (const int position : affected.tuple_positions) times.push_back(event_times[position]);
            std::sort(times.begin(), times.end());
            double diff = 0.;
            // Every edge occurs an even number of times in a kernel relation. The
            // proposed toggles reverse its old spin only between successive pairs.
            for (std::size_t j = 0; j < times.size(); j += 2) {
                if (times[j] != times[j + 1]) {
                    diff -= 2. * lat.integrated_edge_energy(affected.edge, times[j], times[j + 1]);
                }
            }
            edge_diffs[i] = diff;
            action_diff -= diagonal_field * diff;
        }
        acc_ratio = std::exp(log_ratio - action_diff);
        if (!mc.accept(acc_ratio)) return;
        for (std::size_t i = 0; i < relation.size(); ++i) {
            const auto affected = tuple_edges(relation[i]);
            if (create_events[i]) lat.insert_tuple_flip(relation[i], affected, event_times[i]);
            else lat.delete_tuple_flip(relation[i], affected, event_times[i]);
        }
        integrated_pot_energy += action_diff;
        for (std::size_t i = 0; i < edges.size(); ++i) {
            lat.add_potential_edge_energy(edges[i].edge, edge_diffs[i]);
        }
    }
};
} // namespace paratoric
