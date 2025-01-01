#include <cassert>
#include <math.h>

#include "gcem.hpp"

namespace Lookups
{
    // static constexpr double min_speed = 1.0 / 100000; // in steps per microsecond / steptime_high = 100 milliseconds
    // static constexpr double M = 1.0 / 500;            // Assuming steptime_low = 100 microseconds

    static constexpr double max_speed = 5000; // in steps/s

    constexpr double f(double x)
    {
        if (x > 0.0)
        {
            return gcem::pow(M_E, -1.0 / x);
        }
        else
        {
            return 0.0;
        }
    }

    constexpr double f_d(double x)
    {
        if (x > 0.0)
        {
            return f(x) / (x * x);
        }
        else
        {
            return 0.0;
        }
    }

    // y = g(x) describes the position of the blind at time x (for x,y in [0,1])
    constexpr double g(double x)
    {
        return f(x) / (f(x) + f(1.0 - x));
    }

    // y = g'(x) describes the velocity of the blind at time x (for x in [0,1])
    // g'(x) is the derivate of g(x)
    constexpr double g_d(double x)
    {
        if (x > 0 && x < 1)
        {
            double dividend = f_d(x) * (f(x) + f(1.0 - x)) - f(x) * (f_d(x) - f_d(1.0 - x));
            double divisor = f(x) * f(x);

            return dividend / divisor;
        }
        else
        {
            return 0.0;
        }
    }

    constexpr double approx_inv_g(double y, double max_error)
    {
        assert(y >= 0.0 && y <= 1.0);

        // initially approximate g(x) with g(x) ~= x
        double x = y;

        while (gcem::abs(g(x) - y) > max_error)
        {
            // refine approximation with newton's method
            x = x - g(x) / g_d(x);
        }

        return x;
    }

    template <int N>
    struct smoothLookup
    {
        double arr[N + 1];
        int length = N + 1;

        // when regarding the step count as a continuous number (with the highest step count = N)
        // this accuracy number describes the maximum error (by how much we may deviate from the exact
        // integer step value) while approximating the absolute (normalized) time when this step occurs
        static constexpr double relative_step_accuracy = 1.0 / 500.0;

        constexpr smoothLookup() : arr(), length(N)
        {
            for (int i = 0; i <= N; i++)
            {
                double normalized_y = (double)i / N;
                arr[i] = approx_inv_g(normalized_y, 1.0 / N * relative_step_accuracy);
            }
        }
    }

    static constexpr int max_steps = 25500;
    static constexpr Lookups::smoothLookup<max_steps> smooth_lookup = Lookups::smoothLookup<max_steps>();

    constexpr double absolute_normalized_time_at_step(int step, int total_steps)
    {
        assert(step >= 0 && total_steps >= 0 && step <= total_steps);

        if (step == total_steps)
        {
            return smooth_lookup[max_steps];
        }

        // use linear interpolation
        double scale = (double)max_steps / total_steps;
        double scaled_step = step * scale;
        int scaled_step_floored = (int)gcem::floor(scaled_step);

        double x1 = smooth_lookup[scaled_step_floored];
        double x2 = smooth_lookup[scaled_step_floored + 1];

        double gradient = x2 - x1;
        double partial_step = scaled_step - (double)scaled_step_floored;

        return x1 + partial_step * gradient;
    }

    // calculate minimum total time for procedure to adhere to `max_speed` when traveling the given `total_steps`
    constexpr double minimum_total_time_for_steps(int total_steps)
    {
        // determine the maximum speed the normalized function reaches (in normalized steps/s)
        constexpr double normalized_top_speed = g_d(0.5); /* maximum value of g'(x) is at x=0.5 */

        // To determine the total time to achieve a given `max_speed` we must solve for total_time:
        // max_speed = [total_steps * g(x/total_time)]'(total_time/2) = total_steps/total_time * g'(0.5)
        // => total_time = total_steps * g'(0.5) / max_speed

        double total_time = normalized_top_speed / max_speed * total_steps;
        return total_time;
    }
}