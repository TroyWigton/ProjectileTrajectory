#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include "math.h"
#include <functional>
#include <string>
#include "constants.hpp"
#include "types.hpp"
#include "derivative_functions.hpp"

//#define DEBUG // define before local includes to enable debug mode in those files
#include "golden_section.hpp"
#include "integrators.hpp"
#include "simulation.hpp"

struct OptimalAngleSample {
    double k_over_m;
    double optimal_angle_deg;
    double max_distance_m;
};

struct AnglePrecisionResult {
    double refined_optimal_angle;
    double max_distance;
    double required_angle_precision;
    double dist_minus;
    double dist_optimal;
    double dist_plus;
};

AnglePrecisionResult refine_angle_precision(
    const std::function<double(double)>& distance_func,
    double initial_optimal_angle,
    double initial_max_distance,
    double distance_tolerance,
    double bracket_a,
    double bracket_b,
    bool verbose) {

    if (verbose) {
        std::cout << "\nPhase 2: Finding angle precision for distance tolerance "
                  << distance_tolerance << " m:\n";
    }

    double angle_offset = 0.001;
    double low = 0.0;
    double high = 5.0;

    while (high - low > 1e-9) {
        angle_offset = (low + high) / 2.0;

        double dist_plus = distance_func(initial_optimal_angle + angle_offset);
        double dist_minus = distance_func(initial_optimal_angle - angle_offset);

        double max_deviation = std::max(
            std::abs(initial_max_distance - dist_plus),
            std::abs(initial_max_distance - dist_minus));

        if (max_deviation > distance_tolerance) {
            high = angle_offset;
        } else {
            low = angle_offset;
        }
    }

    if (verbose) {
        std::cout << "Angle precision threshold: ±" << angle_offset << " degrees\n";
        std::cout << "\nPhase 3: Refining optimal angle to required precision:\n";
    }

    double refined_a = initial_optimal_angle - 2.0 * angle_offset;
    double refined_b = initial_optimal_angle + 2.0 * angle_offset;
    double refined_optimal_angle = golden_section_search_max(distance_func, refined_a, refined_b, angle_offset / 10.0);
    double refined_max_distance = distance_func(refined_optimal_angle);

    double dist_at_optimal = distance_func(refined_optimal_angle);
    double dist_plus = distance_func(refined_optimal_angle + angle_offset);
    double dist_minus = distance_func(refined_optimal_angle - angle_offset);

    if (verbose) {
        int angle_precision_digits = 1;
        if (angle_offset > 0.0) {
            angle_precision_digits = std::max(1, static_cast<int>(-std::log10(angle_offset)) + 1);
        }

        std::cout << "\n=== RESULTS ===\n";
        std::cout << "Optimal launch angle: " << refined_optimal_angle << " degrees\n";
        std::cout << "Maximum distance: " << refined_max_distance << " m\n";
        std::cout << std::setprecision(angle_precision_digits);
        std::cout << "Required angle precision: ±" << angle_offset << " degrees\n";
        std::cout << std::setprecision(6);
        std::cout << "  (to maintain distance within ±" << distance_tolerance << " m)\n";

        std::cout << "\nVerification at angle tolerance boundaries:\n";
        std::cout << "  Angle: " << (refined_optimal_angle - angle_offset)
                  << "° → Distance: " << dist_minus
                  << " m (Δ = " << (refined_max_distance - dist_minus) << " m)\n";
        std::cout << "  Angle: " << refined_optimal_angle
                  << "° → Distance: " << dist_at_optimal << " m (optimal)\n";
        std::cout << "  Angle: " << (refined_optimal_angle + angle_offset)
                  << "° → Distance: " << dist_plus
                  << " m (Δ = " << (refined_max_distance - dist_plus) << " m)\n";
    }

    return {
        refined_optimal_angle,
        refined_max_distance,
        angle_offset,
        dist_minus,
        dist_at_optimal,
        dist_plus
    };
}

OptimalAngleSample evaluate_optimum_for_drag(
    double k_over_m,
    double v0,
    double h0,
    double g,
    double deltaT,
    double coarse_angle_tol,
    double distance_tolerance)
{
    Simulation simulation(v0, drag_deriv, rk4_step_with_params, g, k_over_m, h0);

    auto distance_func = [&](double angle_deg) {
        return simulation.run(angle_deg, deltaT).distance;
    };

    const double bracket_a = 10.0;
    const double bracket_b = 46.0;

    double optimal_angle = golden_section_search_max(distance_func, bracket_a, bracket_b, coarse_angle_tol);
    double max_distance = distance_func(optimal_angle);

    // tighten the result (binary search + refinement)
    // ... reuse Phase 2 / Phase 3 helpers here 

    return {k_over_m, optimal_angle, max_distance};
}

std::vector<OptimalAngleSample> sweep_drag_values(
    const std::vector<double>& drag_values,
    double v0,
    double h0,
    double g,
    double deltaT,
    double coarse_angle_tol,
    double distance_tolerance)
{
    std::vector<OptimalAngleSample> samples;
    samples.reserve(drag_values.size());

    for (double k_over_m : drag_values) {
        samples.push_back(evaluate_optimum_for_drag(
            k_over_m, v0, h0, g, deltaT, coarse_angle_tol, distance_tolerance));
    }
    return samples;
}

int main(int argc, char* argv[]) {
    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        }
    }
    const double g = 9.81;
    const double deltaT = 0.001;
    const double v0 = 100;
    const double h0 = 0.0;
    const double distance_tolerance = 0.01; // Maximum distance loss (m) tolerable at angle precision boundary

    double k_over_m; // 
    // Vacuum 0.0, GolfBall .004, PingPongBall 0.1, beachball 0.3
    // k_over_m = 0.0;   // No Drag scenario (validation testing)
    k_over_m = 0.0025;
    if (verbose) {
        std::cout << "Using drag coefficient k/m = " << k_over_m << "\n";
        std::cout << "Target distance precision: " << distance_tolerance << " m\n";
    }

    Simulation simulation(v0, drag_deriv, rk4_step_with_params, g, k_over_m, h0);

    auto distance_func = [&](double angle_deg) {
        return simulation.run(angle_deg, deltaT).distance;
    };

    // Create a wide initial bracket that definitely contains the maximum
    double a = 10.0;
    double b = 46.0;
    std::cout << std::fixed << std::setprecision(6);

    double bracket_angle_tol = 0.001;
    double optimal_angle = golden_section_search_max(distance_func, a, b, bracket_angle_tol);
    double max_distance = distance_func(optimal_angle);

    if (verbose) {
        std::cout << "Approximate optimal angle: " << optimal_angle << " degrees\n";
        std::cout << "Maximum distance: " << max_distance << " m\n";
    }

    auto precision_result = refine_angle_precision(
        distance_func,
        optimal_angle,
        max_distance,
        distance_tolerance,
        a,
        b,
        verbose);

    optimal_angle = precision_result.refined_optimal_angle;
    max_distance = precision_result.max_distance;
    double required_angle_precision = precision_result.required_angle_precision;

    std::vector<double> drag_grid = {0.0, 0.001, 0.0025, 0.005, 0.01, 0.1, 0.3};
    auto samples = sweep_drag_values(drag_grid, v0, h0, g, deltaT, bracket_angle_tol, distance_tolerance);
    std::cout << std::endl;
    for (const auto& sample : samples) {
        std::cout << "k/m=" << sample.k_over_m
                  << " → optimal angle " << sample.optimal_angle_deg
                  << "°, distance " << sample.max_distance_m << " m\n";
    }

    return 0;
}