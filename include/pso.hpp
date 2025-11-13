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
pso(const variables<NUM_VARS> &lower_bound,
    const variables<NUM_VARS> &upper_bound, const Problem &problem,
    const size_t max_iter = 1000,
    const Coefficient &coefficients = DEFAULT_COEFFICIENTS,
    const Weight_range &weight_range = DEFAULT_WEIGHT_RANGE,
    const double mu = 0.1) {
  auto calc_weight = [&](size_t iter) {
    return ((static_cast<double>(max_iter - iter) -
             (weight_range.begin - weight_range.end)) /
            static_cast<double>(max_iter)) +
           weight_range.end;
  };

  auto calc_mutation_propablity = [&](size_t iter) {
    const double den = max_iter > 1 ? static_cast<double>(max_iter) - 1.0 : 1.0;
    return std::pow(1 - (static_cast<double>(iter) / den), 1.0 / mu);
  };

  auto swarm = Swarm<SWARM_SIZE, NUM_VARS>(lower_bound, upper_bound, problem);

  auto gBest = swarm.particles[0];

  for (size_t i = 0; i < max_iter; i++) {
    if (const auto current_best = swarm.get_Best();
        current_best.dominates(gBest)) {
      gBest = current_best;
    }

    const auto current_weight = calc_weight(i);
    const auto current_mutation_propablity = calc_mutation_propablity(i);
    swarm.update_particles(gBest, problem, current_weight, coefficients,
                           current_mutation_propablity);
  }
  return {gBest, swarm};
}
} // namespace pso
#endif // PSO_H
