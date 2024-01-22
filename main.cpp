//
//  main.cpp
//  SIC_integral
//
//  Created by イグノタスペベレル on 2024/01/11.
//

#include <bits/stdc++.h>
//#include "MT.h"

using namespace std;
float inf = numeric_limits<float>::infinity();
using ll = long long;
const double radius = 100;
const double GAMMA = pow(10, 6.0 / 10.0);
const double end_time = 10000;
constexpr double PI = 3.14159265358979323846264338;
const double theta = pow(10, 0.4);
//const double delta = 2.0 / pass_loss_exponent;
const double lambda_BS = 1.;
const double POWER = 1.0;
const double Euler_constant = 0.5772156649;
const double alpha = 4.5;
const double delta = 2. / alpha;
double n = 10000.0;


double urand(){
    double m, a;
    m = RAND_MAX + 1.0;
    a = (rand() + 0.5)/m;
    a = (rand() + a)/m;
    return (rand() + a)/m;
}

double exp_dist(double lambda) {
    //double g = genrand_real3();
    double g = urand();
    double tau = - log(1 - g) / lambda;
    return tau;
}


double gx(double x) {
    double dphi = atan(x) / n;
    double sum = 0;
    for (double i = 1; i < n - 1; i++) {
        double d1 = i * dphi;
        double d2 = (i+1) * dphi;
        double sec1 = 1. / cos(d1) / cos(d1);
        double a1 = 1. + pow(tan(d1), alpha / 2.);
        double sec2 = 1. / cos(d2) / cos(d2);
        double a2 = 1. + pow(tan(d2), alpha / 2.);
        sum += (sec2 / a2 + sec1 / a1) * dphi / 2.;
    }
    return sum;
}

double L_u(double lambda, double u2) {
    double a = - lambda * PI * pow(theta, delta) * u2;
    double b = PI * delta / sin(PI * delta) - gx(pow(theta, -delta));
    return exp(a * b);
}

double L_r(double lambda, double r2, double u_r) {
    double a = - lambda * PI * pow(theta, delta) * r2;
    double b = PI * delta / sin(PI * delta) - gx(pow(theta, -delta) * u_r);
    return exp(a * b);
}


int main() {
    double ds = 1. / n;
    double dt = 1. / n;
    
    for (double lambda = 0.5; lambda <= 5.0; lambda += 0.5) {
        double sum1 = 0;
        for (double s = 1; s < n - 1; s++) {
            double sum2 = 0;
            for (double t = 1; t < n - 1; t++) {
                double u_r1 = log(dt * t) / lambda / log(ds * s);
                double u_r2 = log(dt * (t+1)) / lambda / log(ds * s);
                double a1 = 1. / (1. + theta * pow(u_r1, 1. / delta));
                double a2 = 1. / (1. + theta * pow(u_r2, 1. / delta));
                double joutei = a1 * L_u(lambda, -log(dt * t) / lambda / PI) * L_r(lambda, -log(ds * s) / PI, u_r1);
                double katei = a2 * L_u(lambda, -log(dt * (t+1)) / lambda / PI) * L_r(lambda, -log(ds * s) / PI, u_r2);
                sum2 += (joutei + katei) * dt / 2.;
            }
            sum1 += sum2 * ds;
        }
        cout << lambda << " " << sum1 << " " << sum1 * lambda << endl;
    }
    cout << endl;
    
//    double dr = 1. / n, dth = 2 * PI / n;
//    for (double i = 0; i < n; i++) {
//        for (double j = 0; j < n; j++) {
//            
//        }
//    }
    
}
