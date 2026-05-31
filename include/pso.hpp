#ifndef PSO_H
#define PSO_H
#include "particle.hpp"
#include <cmath>
#include <cstddef>
namespace pso {
struct Weight_range {
    double begin{0};
    double end{0};
};
constexpr Weight_range DEFAULT_WEIGHT_RANGE{.begin = 0.1, .end = 0.01};

template <size_t NUM_VARS, size_t SWARM_SIZE> struct Solution {
    Particle<NUM_VARS> gBest{};
    Swarm<SWARM_SIZE, NUM_VARS> swarm{};
};

template <size_t SWARM_SIZE, size_t NUM_VARS>
[[nodiscard]] constexpr Solution<NUM_VARS, SWARM_SIZE>
pso(const variables<NUM_VARS> &lower_bound, const variables<NUM_VARS> &upper_bound,
    const Problem &problem, const size_t max_iter = 1000,
    const Coefficient &coefficients = DEFAULT_COEFFICIENTS,
    const Weight_range &weight_range = DEFAULT_WEIGHT_RANGE, const double mu = 0.1) {
    auto calc_weight = [&](size_t iter) {
        return weight_range.begin + ((weight_range.end - weight_range.begin) *
                                     static_cast<double>(iter) / static_cast<double>(max_iter));
    };

    auto calc_mutation_probability = [&](size_t iter) {
        const double den = max_iter > 1 ? static_cast<double>(max_iter) - 1.0 : 1.0;
        return std::pow(1 - (static_cast<double>(iter) / den), 1.0 / mu);
    };

    auto swarm = Swarm<SWARM_SIZE, NUM_VARS>(lower_bound, upper_bound, problem);

    auto gBest = swarm.get_Best();

    for (size_t i = 0; i < max_iter; i++) {
        const auto current_weight = calc_weight(i);
        const auto current_mutation_probability = calc_mutation_probability(i);
        swarm.update_particles(gBest, problem, current_weight, coefficients,
                               current_mutation_probability);

        if (const auto current_best = swarm.get_Best(); current_best.dominates(gBest)) {
            gBest = current_best;
        }
    }
    return {gBest, swarm};
}
} // namespace pso
#endif // PSO_H
