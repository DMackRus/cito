// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_PARAMS class defines parameters and types that are customized w.r.t.
// the robot and the environment.

// ***** CLASS TYPE ************************************************************
// Robot and environment specific

// This file is SPECIFIC to the model flymanoid.xml

#include "mujoco.h"
#include <Eigen/Dense>

#ifndef CITO_PARAMS_H
#define CITO_PARAMS_H

// USER-SPECIFIC PATHS *********************************************************
namespace paths {
    // model file
    const char *const modelFile = "/home/aykut/Development/cito/src/cito/model/flymanoid.xml";
    // log files
    const char *const logFile = "/home/aykut/Development/cito/mjLogs/mjLog_flymanoid";
    const char *const trajFile = "/home/aykut/Development/cito/mjLogs/traj_flymanoid.txt";
}
// TASK PARAMETERS *************************************************************
namespace task {
    // properties of the joint to be controlled
    const double desiredPoseInput[6] = {0.0, 1.25, 0.0, 0.0, 0.0, 0.0};
    const double desiredVeloInput[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const Eigen::Matrix<double, 6, 1> desiredPose(desiredPoseInput);
    const Eigen::Matrix<double, 6, 1> desiredVelo(desiredVeloInput);
    const int controlJointPos0 = 0;   // index of the first element of the joint position
    // cost function weights
    const double w1 = 1e2;    // weight on deviations in x and y directions
    const double w2 = 1e0;    // weight on deviations in z and orientation
    const double w3 = 5e-3;   // weight on virtual stiffness
    const double w4 = 1e2;    // weight on final velocities
}
// SIMULATION AND MODEL PARAMETERS *********************************************
namespace params {
    // simulation ==============================================================
    const double tf = 2.00;           // [s] final time
    const double tc = 1e-1;           // [s] control sampling period
    const double dt = 2e-3;           // [s] dynamic sampling period
    const int ncts = (int) floor(tf/tc);    // number of control time steps
    const int ndpc = (int) floor(tc/dt);    // number of dynamic time steps per a control step
    // model ===================================================================
    const int nact  = 8;              // number of actuated dof
    const int nfree = 1;              // number of free joints
    // contact related
    const int ncrbt = 4;              // number of contact candidates on the robot (end effectors)
    const int ncenv = 8;              // number of contact candidates in the environment
    const int nsite = ncrbt + ncenv;  // number of sites
    const int npair = 16;             // number of contact pairs
    // specific indices
    const int jact[nact]    = {6,7,8,9,10,11,12,13};  // actuated dof indices
    const int jfree[nfree]  = {0};                    // index of the free joint
    const int bfree[nfree]  = {5};                    // index of the free body
    const int spair1[npair] = {8,8,8,8,               // indices of the sites on the robot
                               9,9,9,9,
                               10,10,10,10,
                               11,11,11,11};
    const int spair2[npair] = {2,3,6,7,               // corresponding indices of the sites in the environment
                               0,1,4,5,
                               2,3,6,7,
                               0,1,4,5};
    // initial pose of the robot and controls
    const double kCon0[npair] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10};
    const double acon[npair]  = {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10};
    const double phi_r = 200;     // parameter for the radius of the distance sphere (200 ~ 1 cm)
    // contact surface normals for each pair
    const double csn[npair*3] = {0,1,0, 0,-1,0, 0,1,0, 0,-1,0,
                                 0,1,0, 0,-1,0, 0,1,0, 0,-1,0,
                                 0,1,0, 0,-1,0, 0,1,0, 0,-1,0,
                                 0,1,0, 0,-1,0, 0,1,0, 0,-1,0};
}
// constant variables for types
const int NU    = params::nact;           // number of actuated joints
const int NPAIR = params::npair;          // number of contact pairs
const int NV    = NU + 6*params::nfree;   // degrees of freedom
const int N     = 2*NV;                   // dimensionality of states
const int M     = NU + NPAIR;             // dimensionality of controls
const int NTS   = params::ncts;           // number of control time
const int NTRAJ = (NTS+1)*N + NTS*M;      // number of trajectory variables
// ***** TYPES *****************************************************************
// eigen+mujoco types for a time instant
typedef Eigen::Matrix<mjtNum, N, 1>                  stateVec_t;
typedef Eigen::Matrix<mjtNum, N, N, Eigen::ColMajor> stateDer_t;
typedef Eigen::Matrix<mjtNum, M, 1>                  ctrlVec_t;
typedef Eigen::Matrix<mjtNum, N, M, Eigen::ColMajor> ctrlDer_t;
typedef Eigen::Matrix<mjtNum, NPAIR, 1>              kConVec_t;
// threaded types for multiple time steps
typedef std::vector<stateVec_t, Eigen::aligned_allocator<stateVec_t>> stateVecThread;
typedef std::vector<stateDer_t, Eigen::aligned_allocator<stateDer_t>> stateDerThread;
typedef std::vector<ctrlVec_t,  Eigen::aligned_allocator<ctrlVec_t>>  ctrlVecThread;
typedef std::vector<ctrlDer_t,  Eigen::aligned_allocator<ctrlDer_t>>  ctrlDerThread;
typedef std::vector<kConVec_t,  Eigen::aligned_allocator<kConVec_t>>  kConVecThread;
// ***** STRUCTURES ************************************************************
struct trajectory
{
    stateVecThread X;
    ctrlVecThread  U;
    stateDerThread Fx;
    ctrlDerThread  Fu;
};

#endif //CITO_PARAMS_H
