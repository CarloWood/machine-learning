// This file contains the body of CubicPolynomial::get_roots
// in order to avoid code duplication in tests that need(ed)
// to add test code to it.

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
  }
  Dout(dc::notice, "C0 = " << C0 << ", C1 = " << C1);

  // Determine if the inflection point is above or below the x-axis.
  bool const inflection_point_y_larger_than_zero = C0 > 0.0;

  static constexpr double sqrt3 = 1.7320508075688773;    // -√3
  // Define the value of the root in the case that C0 is zero.
  double root0 = std::copysign(sqrt3, -C0);

  // Avoid the local extrema and the inflection point because the derivative might be zero there too.
  // We go for the root that is the furthest away from the inflection point.

  // Special case for if zero is a root (i.e. Q(u) = u⋅(u² ± 3)).
  double u = cubic_has_local_extrema ? root0 : 0.0;
#ifdef CWDEBUG
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

  if (AI_LIKELY(C0 != 0.0))       // Is 0 not a root of the cubic?
  {
    // We need the cube root of C0.
    double const cbrtC0 = std::cbrt(C0);

    // Define the value of the root in the case that C0 is large.
    double const root1 = -(cbrtC0 + 1.0 / cbrtC0);

    int max_limit = 10;
    if (std::abs(C0) <= 0.90375741846)
    {
      // If cubic_has_local_extrema and |C0| is less than 0.90375741845959156233304814223072905692
      // then the following third degree polynomial approximates S(C0) such that the maximum
      // relative error in the guessed root is 0.00059.
      //
      //    0.0736968737245806627769712 * |C0|   +
      //    1.1816417812075004469387217 * |C0|^2 +
      //   -0.7259036707909135067571286 * |C0|^3
      //
      // where S(C0) is defined as that `(1 - S(C0))⋅root0 + S(C0)⋅root1` is the exact root.

      double abs_C0 = std::abs(C0);
      double SC0_approximation = (0.073696873724580669 + (1.1816417812075004 - 0.72590367079091356 * abs_C0) * abs_C0) * abs_C0;
      u = (1.0 - SC0_approximation) * root0 + SC0_approximation * root1;
    }
    else
    {
      max_limit = 100;
      u = root1;
    }

#ifdef CWDEBUG
    initial_guess = u;
    Dout(dc::notice, "Initial guess: " << initial_guess);
#endif
    int limit = max_limit;

    // Since the initial guess is already very accurate, it is more than sufficient to
    // determine the resolution of a double around the value of the root (nextafter(u) - u).
    // If we add less than 1.5 times that to u (but more than 0.5 times that) then u will
    // be incremented with more than this resolution delta, which is when another iteration
    // can still improve the result.
    double const epsilon = 6.0 * (std::nextafter(u, std::numeric_limits<double>::infinity()) - u);
    ASSERT(epsilon > 0.0);
    double prev_u;
    double step = std::numeric_limits<double>::infinity();
    double prev_step;
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
      // Make sure that comparing with epsilon doesn't do worse than detecting that u only changed by a single resolution delta.
      double near = step > 0.0 ? std::nextafter(prev_u, std::numeric_limits<double>::infinity())
                               : std::nextafter(prev_u, -std::numeric_limits<double>::infinity());
      // If u only changed a single bit, then step should NOT be larger than epsilon!
      // Unless the difference between the initial guess and the current value of u caused the resolution delta to have become larger.
      if (u == near)
        ASSERT(!(std::abs(step) > epsilon));
#endif
    }
    while (--limit && std::abs(step) > epsilon);
    ASSERT(limit > 0);
    iterations = max_limit - limit;
    ASSERT(max_limit != 10 || iterations <= 3);

    Dout(dc::notice, "Root found: " << u << "; guess: " << initial_guess << " with relative error: " << ((initial_guess - u) / std::abs(u)));
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
#ifdef RANDOM_CUBICS_TEST
      cubic.
#endif
      long_division(roots_out[0], remainder);
    number_of_roots += qp.get_roots(*reinterpret_cast<std::array<double, 2>*>(&roots_out[1]));
  }

  return number_of_roots;
