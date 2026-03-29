#include <Rcpp.h>
using namespace Rcpp;

// Fast tournament simulation in C++.
// rounds_won_in: how many rounds each team has already won.
//   -1 = eliminated, 0 = pending (hasn't played yet),
//   1 = won R64, 2 = won R64+R32, etc.
// Teams auto-advance through rounds they've already won.
// Eliminated teams always lose.

// [[Rcpp::export]]
IntegerMatrix simulate_tournament_cpp(
    NumericVector ratings_in,
    IntegerVector seeds_in,
    IntegerVector region_idx_in,
    IntegerVector rounds_won_in,
    IntegerVector seed_order,
    IntegerVector region_order,
    int n_sims,
    int n_rounds,
    double sd_margin = 11.2,
    double sd_rating_update = 0.17
) {
  int n_teams = ratings_in.size();
  IntegerMatrix result(n_teams * n_sims, n_rounds);

  for (int sim = 0; sim < n_sims; sim++) {
    std::vector<double> ratings(ratings_in.begin(), ratings_in.end());
    std::vector<int> seeds(seeds_in.begin(), seeds_in.end());
    std::vector<int> regions(region_idx_in.begin(), region_idx_in.end());
    std::vector<int> rounds_won(rounds_won_in.begin(), rounds_won_in.end());
    std::vector<int> team_indices(n_teams);
    for (int i = 0; i < n_teams; i++) team_indices[i] = i;

    int n = n_teams;

    for (int round = 0; round < n_rounds && n > 1; round++) {
      // Order teams by bracket position
      std::vector<int> ord(n);
      for (int i = 0; i < n; i++) ord[i] = i;

      if (n == 4) {
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
          return regions[a] < regions[b];
        });
        if (ord.size() == 4) std::swap(ord[1], ord[2]);
      } else if (n == 2) {
        // Championship: keep as-is
      } else {
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
          int reg_a = regions[a], reg_b = regions[b];
          if (reg_a != reg_b) return reg_a < reg_b;
          int pos_a = 99, pos_b = 99;
          for (int j = 0; j < seed_order.size(); j++) {
            if (seeds[a] == seed_order[j]) pos_a = j;
            if (seeds[b] == seed_order[j]) pos_b = j;
          }
          return pos_a < pos_b;
        });
      }

      // Reorder
      std::vector<double> r2(n);
      std::vector<int> s2(n), rg2(n), rw2(n), ti2(n);
      for (int i = 0; i < n; i++) {
        r2[i] = ratings[ord[i]];
        s2[i] = seeds[ord[i]];
        rg2[i] = regions[ord[i]];
        rw2[i] = rounds_won[ord[i]];
        ti2[i] = team_indices[ord[i]];
      }
      ratings = r2; seeds = s2; regions = rg2; rounds_won = rw2; team_indices = ti2;

      // Simulate games
      int n_games = n / 2;
      std::vector<int> winners(n_games);

      for (int g = 0; g < n_games; g++) {
        int i1 = 2 * g, i2 = 2 * g + 1;
        int rw1 = rounds_won[i1], rw2v = rounds_won[i2];

        if (rw1 == -1 && rw2v == -1) {
          // Both eliminated — shouldn't happen, pick i1
          winners[g] = i1;
        } else if (rw1 == -1) {
          // Team 1 eliminated
          winners[g] = i2;
        } else if (rw2v == -1) {
          // Team 2 eliminated
          winners[g] = i1;
        } else if (rw1 > round && rw2v <= round) {
          // Team 1 already won this round, team 2 hasn't
          winners[g] = i1;
        } else if (rw2v > round && rw1 <= round) {
          // Team 2 already won this round, team 1 hasn't
          winners[g] = i2;
        } else if (rw1 > round && rw2v > round) {
          // Both already won this round — both auto-advance, but only one can.
          // This means one of them lost to someone else in reality (shouldn't happen
          // in correct bracket). Pick the one with more rounds won, else simulate.
          if (rw1 > rw2v) winners[g] = i1;
          else if (rw2v > rw1) winners[g] = i2;
          else {
            double diff = ratings[i1] - ratings[i2];
            double margin = R::rnorm(diff, sd_margin);
            winners[g] = (margin > 0) ? i1 : i2;
          }
        } else {
          // Both pending for this round — simulate
          double diff = ratings[i1] - ratings[i2];
          double margin = R::rnorm(diff, sd_margin);
          winners[g] = (margin > 0) ? i1 : i2;
        }
      }

      // Compact winners
      std::vector<double> new_ratings(n_games);
      std::vector<int> new_seeds(n_games), new_regions(n_games), new_rw(n_games), new_ti(n_games);
      for (int g = 0; g < n_games; g++) {
        int w = winners[g];
        new_ratings[g] = ratings[w] + R::rnorm(0, sd_rating_update);
        new_seeds[g] = seeds[w];
        new_regions[g] = regions[w];
        new_rw[g] = rounds_won[w];  // preserve rounds_won (don't reset)
        new_ti[g] = team_indices[w];

        int base_row = sim * n_teams + new_ti[g];
        if (round < n_rounds) {
          result(base_row, round) = 1;
        }
      }

      ratings = new_ratings;
      seeds = new_seeds;
      regions = new_regions;
      rounds_won = new_rw;
      team_indices = new_ti;
      n = n_games;
    }
  }

  return result;
}
