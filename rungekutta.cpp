#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <omp.h>
#include "/home/darkness/Desktop/ClusterLibraries/eigen-3.4.0/Eigen/Dense"

// ------------------------------------------------------------------------------------
int number_of_differential_equations = 3;

Eigen::VectorXd differential_function(Eigen::VectorXd r, double t) {
    Eigen::VectorXd f(number_of_differential_equations);

    f[0] = 10*(r[1] - r[0]);
    f[1] = -r[0]*r[2] + 28*r[0] - r[1];
    f[2] = r[0]*r[1] - 8.0/3.0*r[2];
    
    return f;
}

Eigen::VectorXd time_step(Eigen::VectorXd r, double t, double h) {
    Eigen::MatrixXd k(number_of_differential_equations, 4);
    
    k.col(0) = h*differential_function(r, t);
    k.col(1) = h*differential_function(r + 0.5*k.col(0), t + 0.5);
    k.col(2) = h*differential_function(r + 0.5*k.col(1), t + 0.5);
    k.col(3) = h*differential_function(r + k.col(2), t );
    
    return r + 1.0/6.0*(k.col(0) + 2*k.col(1) + 2*k.col(2) + k.col(3));
}

Eigen::MatrixXd solver(Eigen::VectorXd r_init, double t_init, double t_fint, double h) {
    
    unsigned int size_loop = static_cast<unsigned int>((t_fint-t_init)/h);
    double t = t_init;
    Eigen::VectorXd r = r_init;
    Eigen::MatrixXd solution(size_loop, 1 + number_of_differential_equations);

    for (int i = 0; i < size_loop; i++) {
        solution(i, 0) = t;
        solution(i, Eigen::seq(1, number_of_differential_equations)) = r;
        
        r = time_step(r, t, h);
        t = t + h;
    }
    return solution;
}

int main() {

    Eigen::VectorXd r_init(3);
    double t_init = 0;
    double t_fint = 100;
    double h = 0.001;
    r_init << 0, 1, 1;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    Eigen::MatrixXd sol = solver(r_init, t_init, t_fint, h);
    
    std::ostringstream fn;
    std::ofstream myfile;
	fn << "solution.txt";
	myfile.open(fn.str());

    unsigned int size_loop = static_cast<unsigned int>((t_fint-t_init)/h);
    for (int i = 0; i < size_loop; i++) {
        for (int j = 0; j < 1 + number_of_differential_equations; j++) {
            myfile << sol(i, j) << "    ";
        }
        myfile << "\n";
    }

    myfile.close();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
    std::cout << "Solved in " << (static_cast<float>(elapsed_time)*0.000001) << " ms" <<"\n";


}
