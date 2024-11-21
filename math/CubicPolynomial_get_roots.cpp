// This file contains the body of CubicPolynomial::get_roots
// in order to avoid code duplication in tests that need(ed)
// to add test code to it.

// We classify cubic into four class in order to determine which algorithm to use.
//
// Let P(x) be the monic third degree polynomial:
//
//   P(x) = c0 + c1 x + c2 x² + x³
//
// In that case x_max < x_min (the x coordinates of respectively
// the local maximum and local minimum of the cubic).
//
// P'(x) = c1 + 2 c2 x + 3 x², set to zero to find
//      x_max = (-c2 - sqrt(c2^2 - 3 c1)) / 3
//      x_min = (-c2 + sqrt(c2^2 - 3 c1)) / 3
//
// Assume this cubic has three real roots, r₀, r₁ and r₂ where
// |r₀| ⩽ |r₁| ⩽ |r₂|.
//
// Assume for now that r₂ > 0.
//
// Class A:
//      x_max > 0 and 2 P(x_max) > 3 |P(0)| + P0        (r₀ < r₁ < r₂).
//
//   |r₀| ≪ |r₁| ⩽ |r₂|
//
//   0 <= |r0|/|r1| <= 0.171468
//   0 <= |r0|/|r2| <= 0.105892
//
// Class B:
//      x_max > 0 and 2 P(x_max) <= 3 |P(0)| + P0       (r₀ < r₁ < r₂).
//
//   |r₀| ≲ |r₁| ⩽ |r₂|
//
//   0.0980763 <= |r0|/|r1| <= 1
//   0 <= |r0|/|r2| <= 1
//
// Class C:
//      x_max < 0 and 2 P(x_max) <= 3 |P(0)| + P0       (r₁ < r₀ < r₂).
//
//   |r₀| ≲ |r₁| ⩽ |r₂|
//
//   0.171662 <= |r0|/|r1| <= 1
//   0 <= |r0|/|r2| <= 1
//
// Class D:
//      x_max < 0 and 2 P(x_max) > 3 |P(0)| + P0        (r₁ < r₀ < r₂).
//
//   |r₀| ≪ |r₁| ⩽ |r₂|
//
//   0 <= |r0|/|r1| <= 0.300282
//   0 <= |r0|/|r2| <= 0.300282
//

#if defined(CWDEBUG) || defined(RANDOM_CUBICS_TEST)
#define GETROOTS_ASSIGN_INITIAL_GUESS
#endif

  static constexpr double sqrt3 = 1.7320508075688773;    // √3

#if defined(CWDEBUG) && defined(RESTART)
  RESTART
#endif

#ifndef HALLEY_ITERATIONS_TEST
  if (coefficients_[3] == 0.0)
  {
    // The cubic is actually a quadratic.
    QuadraticPolynomial qp(coefficients_[0], coefficients_[1], coefficients_[2]);
    return qp.get_roots(*reinterpret_cast<std::array<double, 2>*>(&roots_out[0]));
  }

  // Step one: divide all coefficients by coefficients_[3]. This does not change the roots.
  double const c0 = coefficients_[0] / coefficients_[3];
  double const c1 = coefficients_[1] / coefficients_[3];
  double const c2 = coefficients_[2] / coefficients_[3];
  double const d = utils::square(c2) - 3.0 * c1;

  // The cubic is now monic:
  //
  //   p(x) = c0 + c1 x + c2 x² + x³
  //
  // The first derivative is,
  //
  //   p'(x) = c1 + 2 c2 x + 3 x²
  //
  // Setting this to zero gives:
  //                                  -c2 +/- sqrt(c2^2 - 3c1)   -c2 +/- sqrt(d)
  //   0 = c1 + 2 c2 x + 3 x² --> x = ------------------------ = ---------------
  //                                             3                     3
  // Remember if we have local extrema or not.
  bool const cubic_has_local_extrema = d > 0.0;

  // Calculate the inflection point (where the second derivative is zero).
  double const Ix = -c2 / 3.0;
  // Magic (took me several days to find this).
  double const M = 27.0 * c0 + (2.0 * d - 3.0 * c1) * c2;

  double C0, C1;
  double scale;

  // Avoid the local extrema and the inflection point because the derivative might be zero there too.
  // We go for the root that is the furthest away from the inflection point.

  // Note: the Halley_iterations test defines u itself.
  double u;

  // If cubic_has_local_extrema then this value is the root when S(C0) = 0,
  // otherwise it is an approximation of the root for small values of C0.
  double root0;

  if (cubic_has_local_extrema)
  {
    // Transform the cubic into
    //
    //   Q(u) = C0 - 3u + u³
    //
    // After transforming the cubic we have the following four possibilities:
    //
    //          starting point
    //                ↓
    //  --------------+--O--> u-axis (for example)
    //      .-.         /
    //     /   \       /
    //    /     \←----/----- Inflection point
    //   /       \   /
    //            `-´←------ The extreme with the largest absolute y value.
    //       ↑  ↑  ↑  ↑
    //      -1  0  1  2
    //
    // The starting point is set on the positive slope on the side of the local extreme that is the furthest away from the x-axis.
    //
    double const sqrt_d = std::sqrt(d);
    C0 = M / (d * sqrt_d);
    C1 = -3.0;

    // The applied transform means that any root found must be scaled back by multiplying with
    scale = sqrt_d / 3.0;
    // and then adding Ix back.

    // Define the value of the root in the case that C0 is zero.
    root0 = std::copysign(sqrt3, -C0);

    // A special case of the initial (and final) guess of the root in the case that C0 is zero.
    u = root0;
  }
  else
  {
    // Transform the cubic into
    //
    //   Q(u) = C0 + 3u + u³
    //
    //          starting point
    //                ↓
    //                /
    //  -------------O+-> u-axis (for example)
    //              /
    //             /
    //            /
    //         .⋅´←--------- Inflection point
    //        /
    //       /
    //      /
    //       ↑  ↑  ↑  ↑
    //      -1  0  1  2

    double const sqrt_md = std::sqrt(-d);
    C0 = M / (-d * sqrt_md);
    C1 = 3.0;

    // The applied transform means that any root found must be scaled back by multiplying with
    scale = sqrt_md / 3.0;
    // and then adding Ix back.

    // A special case of the initial guess of the root in the case that C0 is zero.
    u = 0.0;
  }
#endif  // no HALLEY_ITERATIONS_TEST
  Dout(dc::notice, "C0 = " << std::setprecision(18) << C0 << ", C1 = " << C1);

  // Determine if the inflection point is above or below the x-axis.
  bool const inflection_point_y_larger_than_zero = C0 > 0.0;

#ifdef GETROOTS_ASSIGN_INITIAL_GUESS
#ifndef RANDOM_CUBICS_TEST
  double
#endif
    initial_guess = u;
#endif

#ifdef RANDOM_CUBICS_TEST
  if (plot_fitter)
  {
    // Plot the cubic.
    plot_fitter->solve(
        [=](double u) -> cairowindow::Point { return {u, C0 + (C1 + utils::square(u)) * u}; },
        plot.viewport());
    plot.add_bezier_fitter(layer, line_style({.line_color = cairowindow::color::green}), *plot_fitter);
  }
#endif

#ifdef HALLEY_ITERATIONS_TEST
  int max_limit = 100;
#else
  if (AI_LIKELY(C0 != 0.0))       // Is 0 not a root of the cubic?
  {
    constexpr int max_limit = 10;
#endif
    double const abs_C0 = std::abs(C0);
    if (!cubic_has_local_extrema && abs_C0 < 2.60987)
      u = (utils::square(C0) / 81.0 - 1.0 / 3.0) * C0;
    else
    {
      // We need the cube root of C0.
      double const cbrtC0 = std::cbrt(C0);

      if (!cubic_has_local_extrema)
        u = -cbrtC0 + 1.0 / cbrtC0;
      else
      {
        // This value is the root when S(C0) = 1.
        double const root1 = -(cbrtC0 + 1.0 / cbrtC0);

#ifndef HALLEY_ITERATIONS_TEST
        double SC0_approximation;
        if (abs_C0 < 6.82)
        {
          if (abs_C0 <= 0.90375741846)
          {
            // If cubic_has_local_extrema and |C0| is less than 0.90375741845959156233304814223072905692
            // then the following second degree polynomial approximates S(C0) such that the maximum
            // relative error in the guessed root is 0.005.
            //
            //    0.2307668111362540090428220 * |C0|   +
            //    0.3991077472117580786261304 * |C0|^2
            //
            // where S(C0) is defined as that `root0 + S(C0)⋅(root1 - root0)` is the exact root;
            // aka S(C0) = (root - root0) / (root1 - root0).

            SC0_approximation = (0.230766811136254009 + 0.399107747211758079 * abs_C0) * abs_C0;

            // This polynomial has been determined with `minus_three_case.cxx` with
            // TILL_C0_half set to 1 (and polynomial_degree = 2).
          }
          else
          {
            // If cubic_has_local_extrema and |C0| is larger than 0.90 but less than 6.82
            // then the following third degree polynomial approximates S(C0) such that the maximum
            // relative error in the guessed root is 0.0042.
            //
            //    0.1336242865900584967672058          +
            //    0.5309536153375157312519206 * |C0|   +
            //   -0.1103045467148938868482953 * |C0|^2 +
            //    0.0074894907819310989348705 * |C0|^3

            SC0_approximation =
              0.133624286590058497 + (0.530953615337515731 + (-0.110304546714893887 + 0.0074894907819310989 * abs_C0) * abs_C0) * abs_C0;

            // This polynomial has been determined with `minus_three_case.cxx` with
            // TILL_C0_half set to 0 (and polynomial_degree = 3).
          }
          u = (1.0 - SC0_approximation) * root0 + SC0_approximation * root1;
        }
        else
          u = root1;
      }
    }
#endif

#ifdef GETROOTS_ASSIGN_INITIAL_GUESS
    initial_guess = u;
    Dout(dc::notice, "Initial guess: " << std::setprecision(18) << initial_guess);
#endif
    int limit = max_limit;

    double step = std::numeric_limits<double>::infinity();
    double prev_u;
    double prev_step;
    // If the relative error is down to this or better than only one more iteration is needed.
    constexpr double max_relative_error_before_last_iteration = 3.5e-6;
    do
    {
      prev_u = u;
      prev_step = step;
      // Calculate Q(u) = C0 + C1 * u + u^3.
      double Q_u = C0 + u * (utils::square(u) + C1);
      // Calculate Q''(u) = 6 * u;
      double half_Qpp_u = 3.0 * u;
      // Calculate Q'(u) = C1 + 3 * u^2.
      double Qp_u = half_Qpp_u * u + C1;
      // Apply Halley's method.
      step = -Q_u * Qp_u / (utils::square(Qp_u) - Q_u * half_Qpp_u);
      u += step;                                                                // uₙ₊₁ = uₙ - Q(u)Q'(u) / (Q'(u)² - ½Q(u)Q"(u))
#ifdef CWDEBUG
      Dout(dc::notice, "Halley: u = " << std::setprecision(18) << u << " (" <<
          std::nextafter(u, std::numeric_limits<double>::infinity()) << "); step = " << step << "; Δu = " << (u - prev_u));
      if (std::abs(prev_u - u) <= std::abs(max_relative_error_before_last_iteration * u))
      {
        // This was the last iteration. Check that one more iteration would not have gotten a better answer.

        // Detect what the CORRECT answer is.
        double Qu = C0 + (C1 + u * u) * u;
        double correct_root = u;
        // The derivative of Q(u) around the root is positive: dQ/du = C1 + 3 u^2, where u^2 >= 3.
        double direction = Qu < 0.0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
        while (Qu != 0.0)
        {
          double next_root = std::nextafter(correct_root, direction);
          double next_Qu = C0 + (C1 + next_root * next_root) * next_root;
          double next_direction = next_Qu < 0.0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
          if (next_direction != direction)
          {
            if (std::abs(next_Qu) < 0.99 * std::abs(Qu))
              correct_root = next_root;
            break;
          }
          Qu = next_Qu;
          correct_root = next_root;
        }
        if (C0 != 0)
        {
          bool close_enough = u == correct_root;
          if (!close_enough)
          {
            double nextafter = std::nextafter(u, std::numeric_limits<double>::infinity());
            double nextbefore = std::nextafter(u, -std::numeric_limits<double>::infinity());
            close_enough = nextafter == correct_root || nextbefore == correct_root;
            if (!close_enough)
            {
              nextafter = std::nextafter(nextafter, std::numeric_limits<double>::infinity());
              nextbefore = std::nextafter(nextbefore, -std::numeric_limits<double>::infinity());
              close_enough = nextafter == correct_root || nextbefore == correct_root;
            }
          }
          if (!close_enough)
          {
            // Do one more Halley iteration.
            double Q_u = C0 + u * (utils::square(u) + C1);
            double half_Qpp_u = 3.0 * u;
            double Qp_u = half_Qpp_u * u + C1;
            u += -Q_u * Qp_u / (utils::square(Qp_u) - Q_u * half_Qpp_u);

            // If now u is within one or two resolution steps from the correct root then we DID stop too soon.
            double nextafter = std::nextafter(u, std::numeric_limits<double>::infinity());
            double nextbefore = std::nextafter(u, -std::numeric_limits<double>::infinity());
            ASSERT(!(u == correct_root || nextafter == correct_root || nextbefore == correct_root ||
              std::nextafter(nextafter, std::numeric_limits<double>::infinity()) == correct_root ||
              std::nextafter(nextbefore, -std::numeric_limits<double>::infinity()) == correct_root));
          }
        }
        else
          u = correct_root;
      }
#endif
    }
    while (--limit && std::abs(prev_u - u) > std::abs(max_relative_error_before_last_iteration * u));
    ASSERT(limit > 0);
#if defined(RANDOM_CUBICS_TEST) || defined(GETROOTS_ASSIGN_ITERATIONS)
    iterations = max_limit - limit;
#ifndef HALLEY_ITERATIONS_TEST
    // Although very little, it happens that we get three iterations.
    ASSERT(max_limit != 10 || iterations <= 3);
#endif
#endif

    Dout(dc::notice|continued_cf, "Root found: " << std::setprecision(18) << u);
#if defined(RANDOM_CUBICS_TEST) || defined(GETROOTS_ASSIGN_ITERATIONS)
    Dout(dc::continued, "; iterations: " << iterations);
#endif
#ifdef GETROOTS_ASSIGN_INITIAL_GUESS
    Dout(dc::finish, "; guess: " << initial_guess << " with relative error: " << ((initial_guess - u) / std::abs(u)));
#else
    Dout(dc::finish, "");
#endif
#ifndef HALLEY_ITERATIONS_TEST
  }

#ifdef RANDOM_CUBICS_TEST
  stop = std::abs(u) > 15.0;
#endif

  roots_out[0] = u * scale + Ix;
  int number_of_roots = 1;

  if (cubic_has_local_extrema)
  {
    // Find the other two roots, if any.
    [[maybe_unused]] double remainder;
    QuadraticPolynomial qp =
#if defined(RANDOM_CUBICS_TEST) || defined(GETROOTS_ASSIGN_ITERATIONS)
      cubic.
#endif
      long_division(roots_out[0], remainder);
    number_of_roots += qp.get_roots(*reinterpret_cast<std::array<double, 2>*>(&roots_out[1]));
  }

  return number_of_roots;
#endif  // HALLEY_ITERATIONS_TEST

#undef GETROOTS_ASSIGN_INITIAL_GUESS
