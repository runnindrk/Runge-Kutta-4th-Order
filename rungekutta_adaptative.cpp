#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <omp.h>
#include "/home/darkness/Desktop/ClusterLibraries/eigen-3.4.0/Eigen/Dense"

// ------------------------------------------------------------------------------------
int number_of_differential_equations = 3;
Eigen::VectorXd f(number_of_differential_equations);
Eigen::MatrixXd k(number_of_differential_equations, 4);
double h;
double t;
Eigen::VectorXd r;
double t_h_h1;
Eigen::VectorXd r_h_h1;
double t_h_h2;
Eigen::VectorXd r_h_h2;
double t_2h;
Eigen::VectorXd r_2h;

Eigen::VectorXd differential_function(Eigen::VectorXd r, double t) {

    f[0] = 10*(r[1] - r[0]);
    f[1] = -r[0]*r[2] + 28*r[0] - r[1];
    f[2] = r[0]*r[1] - 8.0/3.0*r[2];
    
    return f;
}

Eigen::VectorXd time_step(Eigen::VectorXd r, double t, double h) {
    
    k.col(0) = h*differential_function(r, t);
    k.col(1) = h*differential_function(r + 0.5*k.col(0), t + 0.5);
    k.col(2) = h*differential_function(r + 0.5*k.col(1), t + 0.5);
    k.col(3) = h*differential_function(r + k.col(2), t );
    
    return r + 1.0/6.0*(k.col(0) + 2*k.col(1) + 2*k.col(2) + k.col(3));
}

void solver(Eigen::VectorXd r_init, double t_init, double t_fint, double h_init, double delta) {
    
    std::ostringstream fn;
    std::ofstream myfile;
	fn << "solution.txt";
	myfile.open(fn.str());

    h = h_init;
    t = t_init;
    r = r_init;

    t_h_h1 = t_init;
    r_h_h1 = r_init;
    t_h_h2 = t_init;
    r_h_h2 = r_init;

    t_2h = t_init;
    r_2h = r_init;

    while (t < t_fint) {
        for (int j = 0; j < 2; j++) {
            if (j == 0) {
                r_h_h1 = time_step(r, t, h);
                t_h_h1 += h;

                r_h_h2 = time_step(r_h_h1, t_h_h1, h);
                t_h_h2 = t_h_h1 + h;
            }

            if (j == 1) {
                r_2h = time_step(r, t, 2*h);
                t_2h += 2*h;
            }
        }

        double rho = 1;
        double error = sqrt((r_h_h2 - r_2h).cwiseAbs2().sum()); // this criteria might change
        if (error > 0) {
            rho = 30*h*delta/error;
        }
        else {
            rho = 10000;
        }
        
        if (rho >= 1) {
            myfile << t_h_h1 << " ";
            for (int j = 0; j < number_of_differential_equations; j++) {
                myfile << r_h_h1(j) << "  ";
            }
            myfile << "\n";
            myfile << t_h_h2 << " ";
            for (int j = 0; j < number_of_differential_equations; j++) {
                myfile << r_h_h2(j) << "  ";
            }
            myfile << "\n";

            t = t_h_h2;
            r = r_h_h2;
            double h_temp = h*(pow(rho, 1.0/4.0));

            if (h_temp > 2) {
                h = 2*h;
            }
            else
                h = h_temp;
        }
        else {
            h = h*(pow(rho, 1.0/4.0));
        }
    }
}

int main() {

    Eigen::VectorXd r_init(3);
    double t_init = 0;
    double t_fint = 100;
    double h_init = 0.001;
    double delta  = 0.000001;
    r_init << 0, 1, 1;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    solver(r_init, t_init, t_fint, h_init, delta);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
    std::cout << "Solved in " << (static_cast<float>(elapsed_time)*0.000001) << " ms" <<"\n";

}
