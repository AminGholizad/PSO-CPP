#include "particle.hpp"
#include "pso.hpp"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <numbers>
#include <span>
constexpr double pi = std::numbers::pi;
[[nodiscard]] constexpr pso::Cost cost_fcn(std::span<double> variables) {
  pso::Cost cost{.objective = 0.0, .infeasiblity = 0.0};
  for (const auto &var : variables) {
    cost.objective +=
        std::sin(var * 5.0) + // NOLINT(readability-magic-numbers,
                              // cppcoreguidelines-avoid-magic-numbers)
        std::sin(var * 7.0) + // NOLINT(readability-magic-numbers,
        // cppcoreguidelines-avoid-magic-numbers)
        std::sin(var * 11.0); // NOLINT(readability-magic-numbers,
                              // cppcoreguidelines-avoid-magic-numbers)
    cost.infeasiblity += std::sin(var);
  }
  cost.objective = std::abs(cost.objective);
  cost.infeasiblity =
      std::abs((cost.infeasiblity / static_cast<double>(variables.size())) -
               0.7); // NOLINT(readability-magic-numbers,
                     // cppcoreguidelines-avoid-magic-numbers)
  return cost;
}
int main() {
  const size_t swarm_size{200};
  const size_t max_iter{2000};
  const pso::variables lower_bound{0.0, 0.0, 0.0, 0.0};
  const pso::variables upper_bound{pi / 2, pi / 2, pi / 2, pi / 2};
  const auto solution =
      pso::pso<swarm_size>(lower_bound, upper_bound, cost_fcn, max_iter);
  solution.gBest.info();
  std::ofstream file("./resualt.csv");
  solution.swarm.export_csv(file);
  return 0;
}
