// Modified golden-section search for maximization
template<typename Func>
double golden_section_search_max(Func f, double a, double b, double tol, bool verbose = false) {
    
    #ifdef DEBUG
    verbose = true;  // Force verbose in debug mode
    #endif

    if (!(a < b)) {
        throw std::invalid_argument("golden_section_search_max: require a < b");
    }
    if (!(tol > 0.0)) {
        throw std::invalid_argument("golden_section_search_max: require tol > 0");
    }

    // Golden ratio phi = (1 + sqrt(5)) / 2
    constexpr double golden_ratio = 1.61803398874989484820;
    constexpr double inv_golden_ratio = 1.0 / golden_ratio;
    
    double c = b - (b - a) * inv_golden_ratio;
    double d = a + (b - a) * inv_golden_ratio;
    double fc = f(c);
    double fd = f(d);
    
    if (verbose) {
        std::cout << "Starting golden-section maximization in [" << a << ", " << b << "]\n";
        std::cout << "Initial points: " << c << "° (" << fc << "m), " 
                  << d << "° (" << fd << "m)\n";
    }
    
    unsigned int iteration = 0;
    // Use interval width as the stopping criterion
    while ((b - a) > tol) {
        ++iteration;
        if (fc > fd) {  // keep the higher value for maximization
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) * inv_golden_ratio;
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) * inv_golden_ratio;
            fd = f(d);
        }
    
        if (verbose) {
            std::cout << "Iteration " << iteration << ": New bracket: [" << a << ", " << b << "], "
                      << "points: " << c << "° (" << fc << "m), "
                      << d << "° (" << fd << "m)\n";
        }
    }
    
    // Prefer the best evaluated interior point; fallback to midpoint if needed
    double optimal;
    if (std::isfinite(fc) && std::isfinite(fd)) {
        optimal = (fc > fd) ? c : d;
    } else {
        optimal = (a + b) / 2.0;
    }

    if (verbose) {
        std::cout << "Converged to " << optimal << " degrees with distance " 
                  << f(optimal) << " m\n";
    }
    return optimal;
}
