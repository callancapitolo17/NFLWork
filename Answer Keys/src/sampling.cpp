// =============================================================================
// sampling.cpp - Rcpp implementations of Answer Key sampling hot path
// =============================================================================
//
// Replaces R implementations of mean_match() and balance_sample() in Tools.R.
//
// Key optimizations:
// 1. mean_match: std::nth_element (O(M)) replaces setorder full sort (O(M log M))
// 2. balance_sample: incremental sum tracking (O(1)) replaces full rescan (O(M))
//
// =============================================================================

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>

// [[Rcpp::export]]
Rcpp::List mean_match_cpp(
    Rcpp::NumericVector spread_vals,
    Rcpp::NumericVector total_vals,
    Rcpp::NumericVector game_date_num,
    int N,
    double parent_spread,
    double parent_total,
    double ss,
    double st,
    int max_iter_mean,
    double tol_mean
) {
    int M = spread_vals.size();
    if (N > M) N = M;
    if (N <= 0) Rcpp::stop("N must be positive");

    // Work with raw pointers for speed (avoid Rcpp proxy overhead in loops)
    const double* sv = REAL(spread_vals);
    const double* tv = REAL(total_vals);
    const double* gd = REAL(game_date_num);

    // Index array: 0-based internally
    std::vector<int> idx(M);
    std::iota(idx.begin(), idx.end(), 0);

    // Pre-allocate distance array
    std::vector<double> dist(M);

    double adj_spread = parent_spread;
    double adj_total = parent_total;
    double inv_ss = 1.0 / ss;
    double inv_st = 1.0 / st;

    bool mm_converged = false;

    for (int iter = 0; iter < max_iter_mean; ++iter) {
        // Step 1: Compute squared Euclidean distances (vectorized, O(M))
        for (int i = 0; i < M; ++i) {
            double ds = (sv[i] - adj_spread) * inv_ss;
            double dt_val = (tv[i] - adj_total) * inv_st;
            dist[i] = ds * ds + dt_val * dt_val;
        }

        // Step 2: Partial sort — move N smallest-distance elements to front (O(M))
        // Comparator: distance ascending, game_date descending (more recent first as tiebreaker)
        auto cmp = [&](int a, int b) {
            if (dist[a] != dist[b]) return dist[a] < dist[b];
            return gd[a] > gd[b];
        };
        std::nth_element(idx.begin(), idx.begin() + N, idx.end(), cmp);

        // Step 3: Compute means of top-N elements (O(N))
        double sum_s = 0.0, sum_t = 0.0;
        for (int i = 0; i < N; ++i) {
            sum_s += sv[idx[i]];
            sum_t += tv[idx[i]];
        }
        double mean_s = sum_s / N;
        double mean_t = sum_t / N;

        // Step 4: Check convergence
        double err_s = mean_s - parent_spread;
        double err_t = mean_t - parent_total;
        if (std::abs(err_s) < tol_mean && std::abs(err_t) < tol_mean) {
            mm_converged = true;
            break;
        }

        // Step 5: Adjust targets
        adj_spread -= err_s;
        adj_total -= err_t;
    }

    // Final sort: sort top-N by distance for balance_sample ordering (O(N log N))
    auto cmp_final = [&](int a, int b) {
        if (dist[a] != dist[b]) return dist[a] < dist[b];
        return gd[a] > gd[b];
    };
    std::sort(idx.begin(), idx.begin() + N, cmp_final);
    // Sort remaining M-N by distance too (O((M-N) log (M-N)))
    std::sort(idx.begin() + N, idx.end(), cmp_final);

    // Convert to 1-based R indices
    Rcpp::IntegerVector order(M);
    for (int i = 0; i < M; ++i) {
        order[i] = idx[i] + 1;
    }

    return Rcpp::List::create(
        Rcpp::Named("order") = order,
        Rcpp::Named("N") = N,
        Rcpp::Named("converged") = mm_converged
    );
}


// [[Rcpp::export]]
Rcpp::List balance_sample_cpp(
    Rcpp::IntegerVector order_indices,
    Rcpp::IntegerVector actual_cover,
    Rcpp::IntegerVector actual_over,
    int N,
    double target_cover,
    double target_over,
    int tol_error
) {
    int M = order_indices.size();

    // included[i] tracks whether position i in the distance-sorted order is in the sample
    std::vector<bool> included(M, false);
    int sum_cover = 0, sum_over = 0;

    // Initialize: first N positions are included
    for (int i = 0; i < N && i < M; ++i) {
        included[i] = true;
        int row = order_indices[i] - 1;  // convert to 0-based
        sum_cover += actual_cover[row];
        sum_over += actual_over[row];
    }

    int target_cover_count = (int)std::round(target_cover * N);
    int target_over_count = (int)std::round(target_over * N);
    int cover_error = sum_cover - target_cover_count;
    int over_error = sum_over - target_over_count;

    bool converged = false;

    while (true) {
        if (std::abs(cover_error) <= tol_error && std::abs(over_error) <= tol_error) {
            converged = true;
            break;
        }

        bool removal_failed = true;
        bool addition_failed = true;

        // ---- REMOVE: scan top-N positions, worst distance first (i = N-1 downto 0) ----
        // Per Feustel: "Start the search at the 500th game and work upward"
        int remove_both = -1, remove_either = -1;
        for (int i = N - 1; i >= 0; --i) {
            if (!included[i]) continue;
            int row = order_indices[i] - 1;
            int cov = actual_cover[row];
            int ov = actual_over[row];
            bool cov_helps = (cover_error > 0 && cov == 1) || (cover_error < 0 && cov == 0);
            bool ov_helps = (over_error > 0 && ov == 1) || (over_error < 0 && ov == 0);

            if (cov_helps && ov_helps) { remove_both = i; break; }
            if ((cov_helps || ov_helps) && remove_either == -1) remove_either = i;
        }

        int rem = (remove_both >= 0) ? remove_both : remove_either;
        if (rem >= 0) {
            included[rem] = false;
            int row = order_indices[rem] - 1;
            // Incremental update: O(1) instead of O(M) rescan
            sum_cover -= actual_cover[row];
            sum_over -= actual_over[row];
            cover_error = sum_cover - target_cover_count;
            over_error = sum_over - target_over_count;
            removal_failed = false;
        }

        // ---- ADD: scan positions N..M-1, closest first ----
        // Per Feustel: "start searching at game 501"
        int add_both = -1, add_either = -1;
        for (int i = N; i < M; ++i) {
            if (included[i]) continue;
            int row = order_indices[i] - 1;
            int cov = actual_cover[row];
            int ov = actual_over[row];
            bool cov_helps = (cover_error > 0 && cov == 0) || (cover_error < 0 && cov == 1);
            bool ov_helps = (over_error > 0 && ov == 0) || (over_error < 0 && ov == 1);

            if (cov_helps && ov_helps) { add_both = i; break; }
            if ((cov_helps || ov_helps) && add_either == -1) add_either = i;
        }

        int add = (add_both >= 0) ? add_both : add_either;
        if (add >= 0) {
            included[add] = true;
            int row = order_indices[add] - 1;
            // Incremental update: O(1) instead of O(M) rescan
            sum_cover += actual_cover[row];
            sum_over += actual_over[row];
            cover_error = sum_cover - target_cover_count;
            over_error = sum_over - target_over_count;
            addition_failed = false;
        }

        // Per Feustel: if both add and remove fail, cannot improve this sample
        if (removal_failed && addition_failed) {
            converged = false;
            break;
        }
    }

    // Collect included indices (1-based positions in the reordered data.table)
    std::vector<int> result;
    result.reserve(N);
    for (int i = 0; i < M; ++i) {
        if (included[i]) result.push_back(i + 1);  // 1-based
    }

    return Rcpp::List::create(
        Rcpp::Named("included_indices") = Rcpp::wrap(result),
        Rcpp::Named("converged") = converged,
        Rcpp::Named("cover_error") = cover_error,
        Rcpp::Named("over_error") = over_error
    );
}
