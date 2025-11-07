#ifndef PSO_H
#define PSO_H
#include "particle.hpp"
#include <array>
#include <cmath>
#include <cstddef>
namespace pso {
constexpr auto DEFAULT_SWARM_SIZE = 100;

struct Weight_range {
  double begin{0};
  double end{0};
};
constexpr Weight_range DEFAULT_WEIGHT_RANGE{.begin = 0.1, .end = 0.01};
template <size_t Num_Vars, size_t Swarm_Size>
using Swarm = std::array<Particle<Num_Vars>, Swarm_Size>;

template <size_t Num_Vars, size_t Swarm_Size> struct Solution {
  Particle<Num_Vars> gBest{};
  Swarm<Num_Vars, Swarm_Size> swarm{};
};

template <size_t Num_Vars, size_t Swarm_Size = DEFAULT_SWARM_SIZE>
[[nodiscard]] constexpr Solution<Num_Vars, Swarm_Size>
pso(const variables<Num_Vars> &lower_bound,
    const variables<Num_Vars> &upper_bound, const Problem &problem,
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

  Swarm<Num_Vars, Swarm_Size> swarm;
  for (auto &particle : swarm) {
    particle = Particle(lower_bound, upper_bound, problem);
  }
  auto gBest = swarm[0];
  for (size_t i = 0; i < max_iter; i++) {
    if (auto current_best = Particle<Num_Vars>::get_Best(swarm);
        current_best.dominates(gBest)) {
      gBest = current_best;
    }

    auto current_weight = calc_weight(i);
    auto current_mutation_propablity = calc_mutation_propablity(i);
    for (auto &particle : swarm) {
      particle.update(gBest, problem, current_weight, coefficients,
                      current_mutation_propablity);
    }
  }
  return {gBest, swarm};
}
} // namespace pso
#endif // PSO_H
