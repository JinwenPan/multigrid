#include "iostream"
#include "cmath"
#include "cstring"
#include "Grid.h"

double get_time() {	// Gets clock time with milliseconds accuracy and returns it as a double 
    struct timespec a;
    clock_gettime(CLOCK_MONOTONIC, &a);
    double t = (double) a.tv_nsec * 1e-6 + (double) a.tv_sec*1e3;
    return t;
}

void v_cycle(Grid *levels, int L, int iter){ // First paramter is a pointer to 'levels[]' which holds Grid objects
	// 'L' is number of levels and 'iter' is number of v-cycle iterations
	
	double resnorm_new=0., resnorm_old, q, error_norm; // resnorm (residual norm), q (convergence factor)

	levels[L].set_BC_for_u();	// Set Dirichlet boundaries for the finest grid
	levels[L].initialise_b();	// Set the right-hand-side vector for the finest grid

	int N = pow(2,L)+1;
	Grid analytic_sol(N);	// Create a seperate Grid object just for holding the analytical solution (for finding error)
	analytic_sol.set_solution();	// Set the analytic solution

	for (int it=1; it<=iter; ++it){		// Do 'iter' number of v-cycle iterations

		levels[L].relax(2);		// Relax twice on the finest grid

		for (int i=L; i>1 ; --i){	// Go down the v-cycle until the first level (coarsest grid) is relaxed
			levels[i].compute_r();	// Compute r on the current grid
			levels[i].restrict_r(levels[i-1]); 	// Restrict r to the coarser grid
			levels[i-1].set_u_to_zero();	// Set u (i.e. error on the coarser grid) to zero
			levels[i-1].relax(2);	// Relax the coarser grid twice
		}

		for (int i=1; i<L; ++i){	// Go up the v-cycle until the L-th level (finest grid) is relaxed
			levels[i].prolong_u(levels[i+1]); // Prolong u (i.e. error) from the current grid to the finer grid
			levels[i+1].correct_u();	// Correct u on the finer grid
			levels[i+1].relax(1);		// Relax the finer grid once
		}

		levels[L].compute_r();	// Compute residual on the finest grid
		resnorm_old = resnorm_new;	// Save old residual from previous v-cycle
		resnorm_new = levels[L].L2_norm(levels[L].r); // Find the new discrete L2 norm of residual for the finest grid
		q = resnorm_new / resnorm_old; 	// Compute convergence factor
		
		levels[L].compute_e(analytic_sol); 	// Find the error of solution on finest grid w.r.t the analytic solution
		error_norm = levels[L].L2_norm(levels[L].e); // Compute L2_norm of the error on the finest grid

		std::cout << it << " - Res-norm = " << resnorm_new << 
			"\t Conv-factor = " << q << "\t Error = "<< error_norm << std::endl;
	}
}

void v_cycle_Nm(Grid *levels, int L, int iter){	// V-cycle for the Neumann case (every step is the same except boundary conditions and relaxation)

	double resnorm_new=0., resnorm_old, q, error_norm;

	levels[L].set_BC_for_u_Nm();
	levels[L].initialise_b();

	int N = pow(2,L)+1;
	Grid analytic_sol(N);
	analytic_sol.set_solution();

	for (int it=1; it<=iter; ++it){

		levels[L].relax_Nm(2);

		for (int i=L; i!=1 ; --i){
			levels[i].compute_r();
			levels[i].restrict_r(levels[i-1]);
			levels[i-1].set_u_to_zero();
			levels[i-1].relax_Nm(2);
		}

		for (int i=1; i!=L; ++i){
			levels[i].prolong_u(levels[i+1]);
			levels[i+1].correct_u();
			levels[i+1].relax_Nm(1);
		}

		levels[L].compute_r();
		resnorm_old = resnorm_new;
		resnorm_new = levels[L].L2_norm(levels[L].r);
		q = resnorm_new / resnorm_old;

		levels[L].compute_e(analytic_sol);
		error_norm = levels[L].L2_norm(levels[L].e);

		std::cout << it << " - Res-norm = " << resnorm_new << 
			"\t Conv-factor = " << q << "\t Error = "<< error_norm << std::endl;
	}
}

int main(int argc, char *argv[]){

	double runtime;
    double start = get_time();	// Get clock time and save it in 'start'

	int L = std::stoi(argv[1]);		// First command line argument (number of levels) assigned to L
	int iter = std::stoi(argv[2]);	// Second command line argument (number of v-cycles) assigned to iter

	Grid levels[L+1];	// levels[] holds Grid objects for every level
	for (int i=1; i<=L ; ++i){
		int n = (pow(2,i)+1); // Number of grid points in x and y direction is calculated form the level number
		levels[i] = Grid(n); // The i-th level is stored in levels[i]. Note that levels[0] is unused so that the numbering of levels[] is convenient 
	}

	if (argc>3){	// If more than 2 arguments are passed in command line, do the bonous task
		std::cout<< "------ Solving solution for the case of Neumann plus Dirichlet boundary conditions ------" << std::endl;
		v_cycle_Nm(levels, L, iter);}	// Defined in main.cpp
	else {		// else do the base case
		std::cout<< "------ Solving solution for the case of Dirichlet boundary conditions ------" << std::endl;
		v_cycle(levels, L, iter);}	// // Defined in main.cpp 

	double end = get_time(); // Get clock time and save it in 'end'
    runtime = end - start;	// Find elapsed time
    std::cout << "Total computation time: " << runtime << " ms" <<std::endl;

    levels[L].write_output();	// Write solution to file

	return 0;
}