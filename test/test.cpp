#include <sgolay/sgolay.hpp>
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <array>
#include <cmath>


TEST_CASE("Weight generation vs. Reference from paper")
{
    SECTION("5pt, quadratic, smooth")
    {  
        constexpr auto m = sgf::M(2);
        constexpr auto n = sgf::N(2);
        constexpr auto s = sgf::S(0);
        constexpr auto i = sgf::I(-2);
        /* Reference weights pulled from paper (Table I) with the
         * parameters specified above */
        constexpr auto ref = std::array<float, sgf::getWindowSize(m)>{
        /*            T            */
        /*  -2 | -1 | 0  | 1  | 2  */
            31 , 9  , -3 , -5 , 3
        };
        constexpr auto norm = std::array<float, sgf::getWindowSize(m)>{
            35 , 35 , 35 , 35 , 35
        };

        for (int idx = 0; idx < ref.size(); idx++) {
            float w = sgf::detail::weight(i, sgf::T(idx - m), m, n, s);
            REQUIRE_THAT(w1 * norm[idx], WithinAbs(ref[idx], 0.0001));
        }
    }

    SECTION("5pt, quadratic, 1st deriv")
    {  
        constexpr auto m = sgf::M(2);
        constexpr auto n = sgf::N(2);
        constexpr auto s = sgf::S(1);
        constexpr auto i = sgf::I(-2);
        /* Reference weights pulled from paper (Table I) with the
         * parameters specified above */
        constexpr auto ref = std::array<float, sgf::getWindowSize(m)>{
        /*            T             */
        /*  -2 |  -1 | 0  |  1 |  2 */
           -54 , -34 , -2 ,  6 , 26
        };
        constexpr auto norm = std::array<float, sgf::getWindowSize(m)>{
            70 ,  70 , 10 , 70 , 70
        };

        for (int idx = 0; idx < ref.size(); idx++) {
            float w = sgf::detail::weight(i, sgf::T(idx - m), m, n, s);
            REQUIRE_THAT(w1 * norm[idx], WithinAbs(ref[idx], 0.0001));
        }
    }
}

//==============================================================================

int main(int argc, char* argv[])
{
  Catch::Session session; // There must be exactly one instance

  // writing to session.configData() here sets defaults
  // this is the preferred way to set them

  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0) // Indicates a command line error
        return returnCode;

  // writing to session.configData() or session.Config() here
  // overrides command line args
  // only do this if you know you need to

  int numFailed = session.run();

  // numFailed is clamped to 255 as some unices only use the lower 8 bits.
  // This clamping has already been applied, so just return it here
  // You can also do any post run clean-up here
  return numFailed;
}
