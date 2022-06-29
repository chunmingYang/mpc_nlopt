#include <stdio.h>
#include <math.h>
#include "nlopt.h"

/*
    Code Description:
        finished the single step optimization, then verified with matlab.
    Problem Statement:
        we have cost J(u) from MPC concept, according to (https://github.com/chunmingYang/diffCar_controls/blob/main/diffCar_MPC_costFunction.m)
        then we have to min(J(u)), having more details inside the "nodes.txt" on github.
*/

// MPC parameters
#define ts 0.01                         // step time
#define T 5                             // receding horizon
#define iter 200                        // iteration timess
double Q[9] = {1,0,0, 0,1,0, 0,0,1};    // quadrutic function parameters can be tuned
double R[9] = {0,0,0, 0,0,0, 0,0,0};    // quadrutic function parameters can be tuned

// initialization
double X[3] = {0, 0, 0};                // states initialization
double Xref[3] = {0.9, 0, 0};           // reference initialization
double U[15] = {};                      // we have the special index format here (counts by column)

// cost function setup for nlopt optimizer
double mycost(unsigned n, const double *U)
{
    // manual iteration initialization --- the total iteration time decided by the receding horizon
    double XX[3] = {X[0], X[1], X[2]}; // using local variables to avoid error
    double X_Xref[3];                  // X-Xref
    double mul[3];                     // transpose(X-Xref)*Q
    double cost;                       // cost initialization
    double cost_state;                 // transpose(X-Xref)*Q*(X-Xref)
    double A[9];                       // this diff car model doesn't have A matrix
    double B[9];                       // B matrix of the state space model of diff car
    double B_U[3];                     // B*U
    double Xdot[3];                    // state velocity
    
    // iteration 01
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[0]+B[1]*U[1]+B[2]*U[2];
    B_U[1] = B[3]*U[0]+B[4]*U[1]+B[5]*U[2];
    B_U[2] = B[6]*U[0]+B[7]*U[1]+B[8]*U[2];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    // iteration 02
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[3]+B[1]*U[4]+B[2]*U[5];
    B_U[1] = B[3]*U[3]+B[4]*U[4]+B[5]*U[5];
    B_U[2] = B[6]*U[3]+B[7]*U[4]+B[8]*U[5];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    // iteration 03
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[6]+B[1]*U[7]+B[2]*U[8];
    B_U[1] = B[3]*U[6]+B[4]*U[7]+B[5]*U[8];
    B_U[2] = B[6]*U[6]+B[7]*U[7]+B[8]*U[8];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    // iteration 04
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[9]+B[1]*U[10]+B[2]*U[11];
    B_U[1] = B[3]*U[9]+B[4]*U[10]+B[5]*U[11];
    B_U[2] = B[6]*U[9]+B[7]*U[10]+B[8]*U[11];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    // iteration 05
    X_Xref[0] = XX[0]-Xref[0];
    X_Xref[1] = XX[1]-Xref[1];
    X_Xref[2] = XX[2]-Xref[2];
    mul[0] = X_Xref[0]*Q[0]+X_Xref[1]*Q[3]+X_Xref[2]*Q[6];
    mul[1] = X_Xref[0]*Q[1]+X_Xref[1]*Q[4]+X_Xref[2]*Q[7];
    mul[2] = X_Xref[0]*Q[2]+X_Xref[1]*Q[5]+X_Xref[2]*Q[8];
    cost_state = X_Xref[0]*mul[0]+mul[1]*X_Xref[1]+mul[2]*X_Xref[2];
    cost = cost + cost_state;
    B[0] = cos(XX[2]);
    B[1] = -sin(XX[2]);
    B[2] = 0;
    B[3] = sin(XX[2]);
    B[4] = cos(XX[2]);
    B[5] = 0;
    B[6] = 0;
    B[7] = 0;
    B[8] = 1;
    B_U[0] = B[0]*U[12]+B[1]*U[13]+B[2]*U[14];
    B_U[1] = B[3]*U[12]+B[4]*U[13]+B[5]*U[14];
    B_U[2] = B[6]*U[12]+B[7]*U[13]+B[8]*U[14];
    Xdot[0] = B_U[0];
    Xdot[1] = B_U[1];
    Xdot[2] = B_U[2];
    XX[0] = XX[0]+ts*Xdot[0];
    XX[1] = XX[1]+ts*Xdot[1];
    XX[2] = XX[2]+ts*Xdot[2];

    return cost;
}

// optimization main loop
void main()
{
    // optimization main loop iteration 01
    unsigned n = 15;                                                                         // number of decision variables
    double lb[] = {-50,-50,-50,   -50,-50,-50,   -50,-50,-50,   -50,-50,-50,   -50,-50,-50}; // low bound for decision variables
    double ub[] = {50,50,50,   50,50,50,   50,50,50,   50,50,50,   50,50,50};                // up bound for decision variables
    nlopt_opt opt;                                                                           // generate the optimization problem using nlopt api
    opt = nlopt_create(NLOPT_LN_COBYLA, n);                                                  // generate the optimization problem using nlopt api
    nlopt_set_lower_bounds(opt, lb);                                                         // set the low bound to optimizer
    nlopt_set_upper_bounds(opt, ub);                                                         // set the up bound to optimizer
    nlopt_set_min_objective(opt, mycost, NULL);                                              // set the cost function
    nlopt_set_xtol_rel(opt, 1e-4);                                                           // set the optimizer valve to converge
    double minf;                                                                             // final optimized J result
    double U[15] = {0,0,0,   0,0,0,   0,0,0,   0,0,0,   0,0,0};
    if (nlopt_optimize(opt, U, &minf) < 0) 
    {
        printf("nlopt failed!\n");
    }
    else 
    {
        printf("we have the optimized min(J) is %f\n", minf);
        printf("in this case we have corresponding decision variables are:\n");
        printf("*****************\n");
        printf("%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n", U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11],U[12],U[13],U[14]);
        printf("*****************\n");
    }
    nlopt_destroy(opt);
}
