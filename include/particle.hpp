#ifndef PARTICLE_H
#define PARTICLE_H
#include <Random-Helper.hpp>
#include <array>
#include <cstddef>
#include <functional>
#include <iostream>
#include <span>
#include <utility>

namespace pso {
struct Coefficient {
    double personal{0};
    double global{0};
};

constexpr Coefficient DEFAULT_COEFFICIENTS{.personal = 0.2, .global = 0.2};
constexpr double DEFAULT_WEIGHT = 0.5;
constexpr double DEFAULT_MUTATION_PROBABLITY = 0.1;
constexpr double THRESHOLD = 0.05;

struct Cost {
    double objective{0};
    double infeasibility{0};
};

[[nodiscard]] constexpr bool operator<(const Cost &lhs, const Cost &rhs) {
    if (lhs.infeasibility == 0 && rhs.infeasibility > 0)
        return true;
    if (lhs.infeasibility > 0 && rhs.infeasibility == 0)
        return false;
    if (lhs.infeasibility > 0 && rhs.infeasibility > 0)
        return lhs.infeasibility < rhs.infeasibility;
    return lhs.objective < rhs.objective;
}

template <size_t NUM_VARS> using variables = std::array<double, NUM_VARS>;

using Problem = std::function<Cost(std::span<const double>)>;

template <size_t NUM_VARS> class Particle {
  public:
    constexpr Particle() = default;
    constexpr Particle(variables<NUM_VARS> lower, variables<NUM_VARS> upper, const Problem &problem)
        : lower_bound{std::move(lower)}, upper_bound{std::move(upper)} {
        for (size_t i = 0; i < NUM_VARS; i++) {
            position[i] = rnd::unifrnd(lower_bound[i], upper_bound[i]);
            velocity[i] = 0.0;
        }
        cost = problem(position);
        pBest_position = position;
        pBest_cost = cost;
    }

    constexpr void update(const Particle &gBest, const Problem &problem,
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
        return cost < other.cost;
    }

    constexpr void info(std::ostream &out = std::cout) const & {
        out << "particle info:\n";
        out << "\tcost = " << cost.objective << '\n';
        out << "\tinfeasibility = " << cost.infeasibility << '\n';
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
        out << "\t\tinfeasibility = " << pBest_cost.infeasibility << '\n';
        out << "\t\tx=(";
        for (size_t i = 0; i < NUM_VARS - 1; i++) {
            out << pBest_position[i] << ", ";
        }
        out << pBest_position.back() << ")\n";
    }

    constexpr void export_csv(std::ostream &out) const & {
        out << '"';
        for (size_t i = 0; i < NUM_VARS - 1; i++) {
            out << position[i] << ',';
        }
        out << position.back() << "\"," << cost.objective << ',' << cost.infeasibility << ",\"";
        for (size_t i = 0; i < NUM_VARS - 1; i++) {
            out << pBest_position[i] << ',';
        }
        out << pBest_position.back() << "\"," << pBest_cost.objective << ","
            << pBest_cost.infeasibility << '\n';
    }

    constexpr static void export_csv(std::ostream &out, std::span<const Particle> swarm) {
        out << "x,cost,infeasibility,pBest,pBest_cost,pBest_infeasibility\n";
        for (const auto &particle : swarm) {
            particle.export_csv(out);
        }
    }

  private:
    constexpr void updateV(const Particle &gBest, const double weight = DEFAULT_WEIGHT,
                           const Coefficient &coefficients = DEFAULT_COEFFICIENTS) {
        for (size_t i = 0; i < NUM_VARS; i++) {
            velocity[i] =
                (weight * velocity[i]) +
                (coefficients.personal * rnd::rand() * (pBest_position[i] - position[i])) +
                (coefficients.global * rnd::rand() * (gBest.position[i] - position[i]));
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
        if ((cost < pBest_cost) || !(pBest_cost < cost)) {
            pBest_position = position;
            pBest_cost = cost;
        }
    }
    constexpr void Mutate(const Problem &problem,
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
        if (auto new_cost = problem(new_position); (new_cost < cost) || (rnd::rand() < THRESHOLD)) {
            position[candidate] = new_position[candidate];
            cost = new_cost;
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

template <size_t SWARM_SIZE, size_t NUM_VARS> struct Swarm {
    using Particle = Particle<NUM_VARS>;
    std::array<Particle, SWARM_SIZE> particles;

    constexpr Swarm() = default;
    constexpr Swarm(const variables<NUM_VARS> &lower_bound, const variables<NUM_VARS> &upper_bound,
                    const Problem &problem) {
        for (auto &particle : particles) {
            particle = Particle(lower_bound, upper_bound, problem);
        }
    }

    constexpr void update_particles(const Particle &gBest, const Problem &problem,
                                    const double weight, const Coefficient &coefficients,
                                    const double mutation_probablity) {
        for (auto &particle : particles) {
            particle.update(gBest, problem, weight, coefficients, mutation_probablity);
        }
    }

    constexpr void export_csv(std::ostream &out) const { Particle::export_csv(out, particles); }

    [[nodiscard]] constexpr Particle get_Best() const {
        return *std::min_element(particles.begin(), particles.end(),
                                 [](const auto &particle_a, const auto &particle_b) {
                                     return particle_a.dominates(particle_b);
                                 });
    }

    constexpr explicit operator std::span<Particle>() { return particles; }
    constexpr explicit operator std::span<const Particle>() const { return particles; }
};
} // namespace pso
#endif // PARTICLE_H
