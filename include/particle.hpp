#ifndef PARTICLE_H
#define PARTICLE_H
#include "rand.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <iostream>
#include <span>
#include <utility>

namespace pso {
struct Coefficient {
  double personal;
  double global;
};

constexpr Coefficient DEFAULT_COEFFICIENTS{.personal = 0.2, .global = 0.2};
constexpr double DEFAULT_WEIGHT = 0.5;
constexpr double DEFAULT_MUTATION_PROBABLITY = 0.1;
constexpr auto THRESHOLD = 0.5;

struct Cost {
  double objective{0};
  double infeasiblity{0};
};

template <size_t NUM_VARS> using variables = std::array<double, NUM_VARS>;

using Problem = std::function<Cost(std::span<double>)>;

template <size_t NUM_VARS> class Particle {
public:
  constexpr Particle() = default;
  Particle(variables<NUM_VARS> lower, variables<NUM_VARS> upper,
           const Problem &problem)
      : lower_bound{std::move(lower)}, upper_bound{std::move(upper)} {
    for (size_t i = 0; i < NUM_VARS; i++) {
      position[i] = rnd::unifrnd(lower_bound[i], upper_bound[i]);
      velocity[i] = 0.0;
    }
    cost = problem(position);
    pBest_position = position;
    pBest_cost = cost;
  }

  constexpr void
  update(const Particle &gBest, const Problem &problem,
         const double weight = DEFAULT_WEIGHT,
         const Coefficient &coefficients = DEFAULT_COEFFICIENTS,
         const double mutation_probablity = DEFAULT_MUTATION_PROBABLITY) {
    updateV(gBest, weight, coefficients);
    updateX();
    cost = problem(position);
    Mutate(problem, mutation_probablity);
    updatePBest();
  }

  [[nodiscard]] constexpr bool dominates(const Particle &other) const & {
    return ((cost.infeasiblity <= other.cost.infeasiblity) &&
            (cost.objective < other.cost.objective));
  }

  [[nodiscard]] constexpr static Particle get_Best(std::span<Particle> swarm) {
    return *std::min_element(
        swarm.begin(), swarm.end(),
        [](const auto &particle_a, const auto &particle_b) {
          return particle_a.dominates(particle_b);
        });
  }

  constexpr void info(std::ostream &out = std::cout) const & {
    out << "particle info:\n";
    out << "\tcost = " << cost.objective << '\n';
    out << "\tinfeasiblity = " << cost.infeasiblity << '\n';
    out << "\tx=(";
    for (size_t i = 0; i < NUM_VARS - 1; i++) {
      out << position[i] << ", ";
    }
    out << position.back() << ")\n";
    out << "\tv=(";
    for (size_t i = 0; i < NUM_VARS - 1; i++) {
      out << velocity[i] << ", ";
    }
    out << velocity.back() << ")\n";
    out << "\tpBest:" << '\n';
    out << "\t\tcost = " << pBest_cost.objective << '\n';
    out << "\t\tinfeasiblity = " << pBest_cost.infeasiblity << '\n';
    out << "\t\tx=(";
    for (size_t i = 0; i < NUM_VARS - 1; i++) {
      out << pBest_position[i] << ", ";
    }
    out << pBest_position.back() << ")\n";
  }

  constexpr void csv_out(std::ostream &out) const & {
    out << '"';
    for (size_t i = 0; i < NUM_VARS - 1; i++) {
      out << position[i] << ',';
    }
    out << position.back() << "\"," << cost.objective << ','
        << cost.infeasiblity << ",\"";
    for (size_t i = 0; i < NUM_VARS - 1; i++) {
      out << pBest_position[i] << ',';
    }
    out << pBest_position.back() << "\"," << pBest_cost.objective << ","
        << pBest_cost.infeasiblity << '\n';
  }

  constexpr static void csv_out(std::ostream &out, std::span<Particle> swarm) {
    out << "x,cost,infeasiblity,pBest,pBest_cost,pBest_infeasiblity\n";
    for (const auto &particle : swarm) {
      particle.csv_out(out);
    }
  }

private:
  constexpr void
  updateV(const Particle &gBest, const double weight = DEFAULT_WEIGHT,
          const Coefficient &coefficients = DEFAULT_COEFFICIENTS) {
    for (size_t i = 0; i < NUM_VARS; i++) {
      velocity[i] = (weight * velocity[i]) +
                    (coefficients.personal * rnd::rand() *
                     (pBest_position[i] - position[i])) +
                    (coefficients.global * rnd::rand() *
                     (gBest.position[i] - position[i]));
    }
  }
  constexpr void updateX() {
    for (size_t i = 0; i < NUM_VARS; i++) {
      position[i] += velocity[i];
      if (position[i] > upper_bound[i] || position[i] < lower_bound[i]) {
        velocity[i] *= -1;
        position[i] += 2 * velocity[i];
        while (position[i] > upper_bound[i] || position[i] < lower_bound[i]) {
          position[i] -= velocity[i];
          velocity[i] *= -0.5; // NOLINT(readability-magic-numbers,
                               // cppcoreguidelines-avoid-magic-numbers)
          position[i] += velocity[i];
        }
      }
    }
  }
  constexpr void updatePBest() {
    if ((cost.infeasiblity <= pBest_cost.infeasiblity) &&
        (cost.objective < pBest_cost.objective)) {
      pBest_position = position;
      pBest_cost.objective = cost.objective;
      pBest_cost.infeasiblity = cost.infeasiblity;
    }
  }
  constexpr void
  Mutate(const Problem &problem,
         const double mutation_probablity = DEFAULT_MUTATION_PROBABLITY) {
    if (rnd::rand() > mutation_probablity) {
      return;
    }
    const auto candidate = rnd::unifrnd<size_t>(0, NUM_VARS - 1);

    const double delta_x =
        (upper_bound[candidate] - lower_bound[candidate]) * mutation_probablity;
    const double new_lower_bound =
        std::max(position[candidate] - delta_x, lower_bound[candidate]);
    const double new_upper_bound =
        std::min(position[candidate] + delta_x, upper_bound[candidate]);
    auto new_position = position;
    new_position[candidate] = rnd::unifrnd(new_lower_bound, new_upper_bound);
    auto new_cost = problem(new_position);
    if ((new_cost.infeasiblity < cost.infeasiblity &&
         new_cost.objective < cost.objective) ||
        (rnd::rand() < THRESHOLD)) {
      position[candidate] = new_position[candidate];
      cost.objective = new_cost.objective;
      cost.infeasiblity = new_cost.infeasiblity;
    }
  }

  variables<NUM_VARS> lower_bound{};
  variables<NUM_VARS> upper_bound{};
  variables<NUM_VARS> position{};
  variables<NUM_VARS> velocity{};
  variables<NUM_VARS> pBest_position{};
  Cost cost{};
  Cost pBest_cost{};
};
} // namespace pso
#endif // PARTICLE_H
