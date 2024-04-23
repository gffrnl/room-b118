/*   libb118
 *
 *   modules/frlap/b118/frlap/gdm/coefficients/spectral.hpp
 *
 *   Huang & Oberman coefficients - quadratic interpolation
 *
 *   Copyright (C) 2024   Guilherme F. Fornel        <gffrnl@gmail.com>
 *                        Fabio Souto de Azevedo     <fazedo@gmail.com>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <b118/almost_equal.hpp>
#include <b118/numbers.hpp>
#include "../../../frlap.hpp"
#include "./method.hpp"
#include "./generator.hpp"

namespace b118 { namespace frlap { namespace gdm { namespace coefficients {

template<typename Real> struct spectral;

template<>
class spectral<double> final : public method<double>,
                               public generator<double, spectral> {
    using generator<double, spectral>::ealpha;
    using generator<double, spectral>::deltax;

 public:
    spectral(double ealpha, double deltax)
        : generator<double, spectral>(ealpha, deltax) {
        if (ealpha == 0.0) {
            sin_v = 0.0;
            cos_v = 1.0;
            sinc_v = b118::numbers::pi / 2.0;
            sinm1_v = -1.0;
            tgamma_v = 1.0;  // gamma(alpha + 1) = 0!
            // gamma(0) is not defined
            pi_alpha = 1.0;  // pi^0 = 1
        } else if (ealpha == 1.0) {
            sin_v = 1.0;
            cos_v = 0.0;
            sinc_v = 1.0;
            sinm1_v = 0.0;
            tgamma_v = 1.0;  // gamma(alpha + 1) = 1!
            tgammam1_v = 0.0;  // gamma(alpha) - 1 = 0! - 1
            pi_alpha = b118::numbers::pi;
        } else if (ealpha == 2.0) {
            sin_v = 0.0;
            cos_v = -1.0;
            sinc_v = 0.0;
            sinm1_v = -1.0;
            tgamma_v = 2.0;  // gamma(alpha + 1) = 2!
            tgammam1_v = 0.0;  // gamma(alpha) - 1 = 0! - 1
            pi_alpha = b118::numbers::pi * b118::numbers::pi;
        } else {
            sin_v = std::sin(b118::numbers::pi * ealpha / 2);
            cos_v = std::cos(b118::numbers::pi * ealpha / 2);
            sinc_v = sin_v / ealpha;
            sinm1_v = -cos_v*cos_v/(1.0 + sin_v);
            tgamma_v = std::tgamma(ealpha + 1);
            tgammam1_v = std::expm1(std::lgamma(ealpha));  // gamma(alpha) - 1
            pi_alpha = std::pow(b118::numbers::pi, ealpha);
        }
        const_mult = std::pow(deltax, -ealpha) * ealpha / M_PI;
    }

    double operator()(std::size_t const & k) const override {
        if (k == 0) {
            return std::pow(b118::numbers::pi/deltax, ealpha) / (ealpha + 1.0);
        } else if (k <= k_small) {
            return mu_from_fit(k);
        } else {
            return mu_from_assymptotics(k) * const_mult;
        }
    }

 private:
    // Constant object members
    double sin_v;        // std::sin(M_PI * alpha / 2)
    double cos_v;        // std::cos(M_PI * alpha / 2)
    double sinc_v;       // std::sinc(M_PI * alpha / 2)/alpha

    double sinm1_v;      // -sq(std::cos(M_PI * alpha / 2.0))
                         //          /(1.0 + std::sin(M_PI * alpha / 2))

    double tgamma_v;     // std::tgamma(alpha + 1)
    double tgammam1_v;   // std::lgamma(alpha) - 1
    double pi_alpha;     // std::pow(M_PI, alpha)

    double const_mult;   // std::pow(deltax, -alpha) * alpha / M_PI;

    // Constant class members
    static constexpr std::size_t k_small = 10;

    // Mathemathical constants, TODO(gffrnl): PUT in b118::numbers;
    static constexpr double gamma_v = 0.577215664901532861;

    // Coefficients up to 10
    static double coef_1(double ealpha) {
        static constexpr double P[] = {
            -6.366197723675813440769e-01,
            -2.059297839733384303484e-01,
            -2.292961545633005969750e-02,
            -8.681444663645901367220e-04,
            -1.582316408178944320388e-09,
            +3.918977599710780955764e-11};
        static constexpr double Q[] = {
            1.000000000000000000000e+00,
            1.147612567595240351553e+00,
            4.749649121779834721790e-01,
            9.072899514665789767032e-02,
            8.127542932033228903687e-03,
            2.763512269646208214064e-04};

        double x = ealpha - 1;
        double numer{P[5]};
        double denom{Q[5]};
        for (int i = 4; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_2(double ealpha) {
        static constexpr double P[] = {
            -7.759291740995763658089e-01,
            -2.829217703086196658121e-01,
            -5.637457780907778450443e-02,
            -7.007412223701844205798e-03,
            -6.193816376638685183221e-04,
            -3.553519937251785659484e-05,
            -1.367979103048156252422e-06};
        static constexpr double Q[] = {
            +1.000000000000000000000e+00,
            +7.075928393432178267543e-01,
            +7.411465822590414278618e-02,
            -1.733490210731306556679e-02,
            -7.574067579592451475499e-04,
            +2.249821654221352245262e-04,
            -9.909862606418638557347e-06};

        double x = ealpha - 1;
        double numer{P[6]};
        double denom{Q[6]};
        for (int i = 5; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_3(double ealpha) {
        double P[] = {
            -6.366197723675813431208e-01,
            -1.181856265161497932998e-01,
            -3.983716909337714229956e-01,
            -1.563482404044122387383e-01,
            -3.419095330183582834564e-02,
            -4.674327221765267860909e-03,
            -4.426622442597902427640e-04,
            -2.677498711812969262154e-05,
            -9.775269296293205906909e-07};
        double Q[] = {
            +1.000000000000000000000e+00,
            +4.920020775703081555415e-01,
            -6.167273791302243451617e-02,
            -2.177844016897405044047e-02,
            +3.972871056744789555343e-03,
            +6.067168191454400087860e-05,
            -6.673454339958298168532e-05,
            +6.476978402921136039895e-06,
            -2.218498388536791430567e-07};

        double x = ealpha - 1;
        double numer{P[8]};
        double denom{Q[8]};
        for (int i = 7; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_4(double ealpha) {
        static constexpr double P[] = {
            -9.913304792854260644113e-01L,
            -3.604465765936957834959e-01,
            -1.926408341999085247102e-01,
            -6.676082389958020871953e-02,
            -1.544232597594299740411e-02,
            -2.380213741653500629811e-03,
            -2.565051478564905242352e-04,
            -1.792922349756083858332e-05,
            -7.290004803716587430137e-07};
        static constexpr double Q[] = {
            +1.000000000000000000000e+00L,
            +3.982369752877753372544e-01,
            -1.081539473708352530033e-01,
            -1.619145220444328044554e-02,
            +6.086186982755490323047e-03,
            -2.790053161910219408667e-04,
            -7.944363051877226811732e-05,
            +1.192960335632360789574e-05,
            -5.377198840244050502600e-07};

        double x = ealpha - 1;
        double numer{P[8]};
        double denom{Q[8]};
        for (int i = 7; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_5(double ealpha) {
        static constexpr double P[] = {
            -6.366197723675813431116e-01,
            -1.132929113788398219819e-01,
            -6.386856002702607426713e-01,
            -2.319230172068818600686e-01,
            -7.107160881892270445913e-02,
            -1.683779180893674303425e-02,
            -2.999191267001301728827e-03,
            -3.661212171405373457705e-04,
            -2.891219623822112043671e-05,
            -1.196294595409810400987e-06};
        static constexpr double Q[] = {
            +1.000000000000000000000e+00,
            +2.322333051171824434270e-01,
            -1.686247406724822322708e-01,
            +5.034737292322581875708e-03,
            +8.433450741235804247939e-03,
            -1.533264413828341522702e-03,
            +2.413414174372103081538e-06,
            +2.988167850255876128913e-05,
            -3.886356922263679390340e-06,
            +1.775387177125740664526e-07};

        double x = ealpha - 1;
        double numer{P[9]};
        double denom{Q[9]};
        for (int i = 8; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_6(double ealpha) {
        static constexpr double P[] = {
            -1.119328559217149816456e+00,
            -3.774027075730150242996e-01,
            -3.036889963674641126427e-01,
            -1.128760354600728122843e-01,
            -2.851553552537801321030e-02,
            -5.334526442559530124315e-03,
            -7.815085890729471000657e-04,
            -8.652175625909772587552e-05,
            -6.617987967027415450828e-06,
            -2.872208941784029013527e-07};
        static constexpr double Q[] = {
            +1.000000000000000000000e+00,
            +1.847105866560383586572e-01,
            -1.839046179052694885812e-01,
            +1.210987308706122274892e-02,
            +8.950647649050410760246e-03,
            -1.951263799975775258156e-03,
            +3.535091652686901164834e-05,
            +3.579811918968853779267e-05,
            -5.045309440502321649547e-06,
            +2.360850467196300692369e-07};

        double x = ealpha - 1;
        double numer{P[9]};
        double denom{Q[9]};
        for (int i = 8; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_7(double ealpha) {
        static constexpr double P[] = {
            -6.366197723675813433915e-01,
            -1.313982442799666563455e-01,
            -8.192210713297660515892e-01,
            -3.093250531212002864565e-01,
            -1.161930635137453542677e-01,
            -3.267468318586470084177e-02,
            -6.417603679860046273259e-03,
            -8.619419248285984557777e-04,
            -7.474804962830469146716e-05,
            -3.718248126938508891003e-06};
        static constexpr double Q[] = {
            +1.000000000000000000000e+00,
            +9.339602775251610107010e-02,
            -2.069999684787325933570e-01,
            +2.852588834259701921674e-02,
            +9.135575279360643552115e-03,
            -2.974015018870659194391e-03,
            +1.660624005887600592604e-04,
            +4.920957378633904081284e-05,
            -9.358564500451040259930e-06,
            +5.313222025292598314613e-07};

        double x = ealpha - 1;
        double numer{P[9]};
        double denom{Q[9]};
        for (int i = 8; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_8(double ealpha) {
        static constexpr double P[] = {
            -1.210518378665849886399e+00,
            -4.705979536760370495798e-01,
            -4.318731991300256296135e-01,
            -1.795611320887010742946e-01,
            -5.165456879556205502871e-02,
            -1.159728053603521662630e-02,
            -2.002804431097856561209e-03,
            -2.602712585795104360641e-04,
            -2.262287494211915644076e-05,
            -1.165203046489235364453e-06};
        static constexpr double Q[] = {
            +1.000000000000000000000e+00,
            +1.014995029616300718625e-01,
            -2.094113759867870555031e-01,
            +2.681273948915745606830e-02,
            +9.972781039797130578565e-03,
            -3.036908541967370299212e-03,
            +1.348674477897301571139e-04,
            +5.749085315677258350305e-05,
            -1.011569660649543222883e-05,
            +5.512116802699264269481e-07};

        double x = ealpha - 1;
        double numer{P[9]};
        double denom{Q[9]};
        for (int i = 8; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_9(double ealpha) {
        static constexpr double P[] = {
            -6.366197723675813442203e-01,
            -1.570064393190070275588e-01,
            -9.657821889034147624603e-01,
            -3.943420774033056905967e-01,
            -1.656581796412698690422e-01,
            -5.251279224546032568510e-02,
            -1.098591490399947662539e-02,
            -1.640883122409694358546e-03,
            -1.584710416747311185655e-04,
            -9.808245277913587130547e-06};
        static constexpr double Q[] = {
            +1.000000000000000000000e+00,
            +8.364767804663631321282e-03,
            -2.282992755732345746752e-01,
            +4.688523114789032207492e-02,
            +9.550033789972392915305e-03,
            -4.511312915271598335823e-03,
            +3.713301944162669106982e-04,
            +8.044591510960301968304e-05,
            -1.939912103338426163548e-05,
            +1.281282482064112927509e-06};

        double x = ealpha - 1;
        double numer{P[9]};
        double denom{Q[9]};
        for (int i = 8; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }

    static double coef_10(double ealpha) {
        static constexpr double P[] = {
            -1.281368483991673320070e+00,
            -5.806619792384043326929e-01,
            -5.545544036957355009114e-01,
            -2.581853378548760646070e-01,
            -8.045577459119633669278e-02,
            -2.044372857644134313770e-02,
            -3.882942557215586938672e-03,
            -5.598000213437109457659e-04,
            -5.290328388126457794829e-05,
            -3.101987499652022207225e-06};
        static constexpr double Q[] = {
            +1.000000000000000000000e+00,
            +6.042724705199052308765e-02,
            -2.314978046932805685846e-01,
            +3.603456480649801602474e-02,
            +1.254058490890153216257e-02,
            -4.381848103520755694525e-03,
            +2.075652415696184525883e-04,
            +1.041056590982561133958e-04,
            -1.991846677612194718337e-05,
            +1.183243311078138411900e-06};

        double x = ealpha - 1;
        double numer{P[9]};
        double denom{Q[9]};
        for (int i = 8; i >= 0; --i) {
            numer = (numer * x) + P[i];
            denom = (denom * x) + Q[i];
        }
        return numer/denom;
    }


    // a partir de 10.
    static double Y(std::size_t k, double beta) {
        double const f = 1.0 / ((k * b118::numbers::pi)
                             *  (k * b118::numbers::pi));
        std::size_t const kpi2 = std::ceil(1.0/f);
        double prod = 1.0;
        double sum = 0.0;
        std::size_t j;
        // critério com alpha = 1
        for (j = 0; (j+2)*(j+1) <= kpi2; j += 2) {
            prod *= -(beta - j) * (beta - j - 1) * f;
            const double old_sum = sum;
            sum += prod;
            if (sum == old_sum) {
                return sum;
            }
        }

        prod *= -(beta - j) * (beta - j - 1) * f;
        return sum + prod;  //  * approx_J(k, beta - j - 2);
        return sum;
    }

    double mu_from_fit(std::size_t k) const {
        static double (*coef[])(double)  = {
            coef_1, coef_2, coef_3, coef_4, coef_5,
            coef_6, coef_7, coef_8, coef_9, coef_10};

        double I_c = coef[k-1](ealpha)
                * pi_alpha
                * std::pow(k*deltax, -ealpha)
                / (b118::numbers::pi * k);

        if (k % 2 == 0) {
            return ealpha * (1.0 - ealpha) * I_c;
        } else {
            return ealpha * I_c;
        }
    }

    double mu_from_assymptotics(std::size_t k) const {
        double sgn = (k % 2 == 0) ? -1. : 1.;

        if (ealpha < 1e-18) {
            return -(b118::numbers::pi/2
                        + sgn
                        * (1 + Y(k, -1))/(b118::numbers::pi*k)
                    )
                * (k % 2 == 0 ? (1.0 - ealpha) : 1.0)
                / k;

        } else if (ealpha == 1 && k % 2 == 0) {
            return 0.0;

        } else if (ealpha == 1) {
            return -2.0/(k*k);

        } else if (k % 2 != 0 || ealpha < 0.5) {
            double E_1 = std::pow(k, 1.0 - ealpha)
                    * tgamma_v
                    * sinc_v;
            double E_2 = (1.0 + Y(k, ealpha - 1.0))
                    * (pi_alpha/b118::numbers::pi);

            return -(E_1 + sgn *  E_2)/(k*k);

        } else {
            double E_1 = tgammam1_v * sin_v;
            double E_2 = sinm1_v;
            double E_3 = -std::expm1((ealpha - 1.0) * std::log(k * M_PI));
            double E_4 = -Y(k, ealpha - 1.0) * pi_alpha /(b118::numbers::pi*k);
            double E_5 = (E_1 + E_2 + E_3) * std::pow(k, -ealpha);

            return -(E_4 + E_5) / k;
        }
    }
};

}}}}  // end namespace b118::frlap::gdm::coefficients
