// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Fast tournament simulation in C++.
// Takes team ratings, seeds, region indices, statuses (0=pending, 1=advanced),
// bracket seed order, region order, and runs n_sims full tournament simulations.
// Returns a matrix: rows = teams * n_sims, cols = 6 round indicators + sim_id.

// [[Rcpp::export]]
IntegerMatrix simulate_tournament_cpp(
    NumericVector ratings_in,
    IntegerVector seeds_in,
    IntegerVector region_idx_in,  // 1-indexed region assignment
    IntegerVector status_in,      // 0 = pending, 1 = advanced
    IntegerVector seed_order,     // bracket pairing order: {1,16,8,9,5,12,4,13,6,11,3,14,7,10,2,15}
    IntegerVector region_order,   // order of regions for F4 pairing
    int n_sims,
    int n_rounds,                 // how many rounds to simulate
    double sd_margin = 11.2,
    double sd_rating_update = 0.17
) {
  int n_teams = ratings_in.size();

  // Output: n_teams rows * n_sims, columns = n_rounds
  // We'll flatten: result[sim * n_teams + team_idx, round] = 0 or 1
  IntegerMatrix result(n_teams * n_sims, n_rounds);

  for (int sim = 0; sim < n_sims; sim++) {
    // Copy working arrays for this sim
    std::vector<double> ratings(ratings_in.begin(), ratings_in.end());
    std::vector<int> seeds(seeds_in.begin(), seeds_in.end());
    std::vector<int> regions(region_idx_in.begin(), region_idx_in.end());
    std::vector<int> statuses(status_in.begin(), status_in.end());
    std::vector<int> team_indices(n_teams);
    for (int i = 0; i < n_teams; i++) team_indices[i] = i;

    int n = n_teams;

    for (int round = 0; round < n_rounds && n > 1; round++) {
      // Order teams by bracket position
      std::vector<int> ord(n);
      for (int i = 0; i < n; i++) ord[i] = i;

      if (n == 4) {
        // F4: pair by region order (1v3, 2v4)
        // Sort by region_order index
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
          return regions[a] < regions[b];
        });
        // Cross-bracket: swap positions 1 and 2 (0-indexed: swap [1] and [2])
        if (ord.size() == 4) {
          std::swap(ord[1], ord[2]);
        }
      } else if (n == 2) {
        // Championship: keep as-is
      } else {
        // Regional rounds: sort by (region_order_idx, seed_bracket_position)
        // seed_bracket_position = position in seed_order array
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
          int reg_a = regions[a], reg_b = regions[b];
          if (reg_a != reg_b) return reg_a < reg_b;
          // Find position of seed in seed_order
          int pos_a = 99, pos_b = 99;
          for (int j = 0; j < seed_order.size(); j++) {
            if (seeds[a] == seed_order[j]) pos_a = j;
            if (seeds[b] == seed_order[j]) pos_b = j;
          }
          return pos_a < pos_b;
        });
      }

      // Reorder all arrays by ord
      std::vector<double> r2(n); std::vector<int> s2(n), rg2(n), st2(n), ti2(n);
      for (int i = 0; i < n; i++) {
        r2[i] = ratings[ord[i]];
        s2[i] = seeds[ord[i]];
        rg2[i] = regions[ord[i]];
        st2[i] = statuses[ord[i]];
        ti2[i] = team_indices[ord[i]];
      }
      ratings = r2; seeds = s2; regions = rg2; statuses = st2; team_indices = ti2;

      // Simulate games: odd vs even
      int n_games = n / 2;
      std::vector<int> winners(n_games);

      for (int g = 0; g < n_games; g++) {
        int i1 = 2 * g, i2 = 2 * g + 1;

        // Auto-advance logic
        if (statuses[i1] == 1 && statuses[i2] == 0) {
          winners[g] = i1;
        } else if (statuses[i2] == 1 && statuses[i1] == 0) {
          winners[g] = i2;
        } else {
          // Simulate
          double diff = ratings[i1] - ratings[i2];
          double margin = R::rnorm(diff, sd_margin);
          winners[g] = (margin > 0) ? i1 : i2;
        }
      }

      // Compact winners
      std::vector<double> new_ratings(n_games);
      std::vector<int> new_seeds(n_games), new_regions(n_games), new_statuses(n_games), new_ti(n_games);
      for (int g = 0; g < n_games; g++) {
        int w = winners[g];
        new_ratings[g] = ratings[w] + R::rnorm(0, sd_rating_update);
        new_seeds[g] = seeds[w];
        new_regions[g] = regions[w];
        new_statuses[g] = 0;  // all pending after first round
        new_ti[g] = team_indices[w];

        // Record progress
        int base_row = sim * n_teams + new_ti[g];
        if (round < n_rounds) {
          result(base_row, round) = 1;
        }
      }

      ratings = new_ratings;
      seeds = new_seeds;
      regions = new_regions;
      statuses = new_statuses;
      team_indices = new_ti;
      n = n_games;
    }
  }

  return result;
}
