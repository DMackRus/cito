/*! Parameters */
/**
 *  \brief CitoParams contains user-specific and general definitions
 *
 *  This class defines global variables that are specific to simulation,
 *  robot, and environment as well as general types and structures.
 *
 *  \author Aykut Onol
 */

#include <stdio.h>
#include <iostream>
#include <string>

#include "mujoco.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <yaml-cpp/yaml.h>

#ifndef CITO_PARAMS_H
#define CITO_PARAMS_H

typedef Eigen::MatrixXd eigDbl;
typedef Eigen::Matrix<mjtNum, Eigen::Dynamic, Eigen::Dynamic> eigMjc;

class CitoParams{
public:
    /// Constructor
    CitoParams(const mjModel *model);
    /// Destructor
    ~CitoParams() {}
    /// Parameters
    const mjModel *model;
    double tf, tc, dt;
    int N, ndpc, nu, nv, n, m, nTraj, nPair,
        nFree, *pFree, *bFree, *dAct,
        *quatAdr, *dofAdr;
    Eigen::VectorXi sPair1, sPair2;
    eigMjc nCS;
};

/// User-specific paths
namespace paths {
    const std::string workspaceDir = std::getenv("CITO_WS");
    const std::string modelFile    = "sawyer_push.xml";
}

/// Specify the parameters related to simulation, model, and contact candidates
/// Parameters
namespace params {
    /// Simulation
    const double tf = 1.00;                 // [s] final time
    const double tc = 5e-2;                 // [s] control sampling period
    const double dt = 2e-3;                 // [s] dynamic sampling period
    const int ncts = (int) floor(tf/tc);    // number of control time steps
    const int ndpc = (int) floor(tc/dt);    // number of dynamic time steps per a control step
    /// Model
    const int nact  = 7;                    // number of actuated joints
    const int nfree = 1;                    // number of free joints
    // specific joint and body indices
    const int jact[nact] = {6, 7, 8, 9, 10, 11, 12};    // indices of actuated DOF
    const int jfree[nfree] = {0};                       // indices of free joints
    const int bfree[nfree] = {1};                       // indices of free bodies
    /// Contact
    const int npair = 1;                // number of contact pairs
    // site indices for each contact pair
    const int spair1[npair] = {4};      // indices of the sites on the robot
    const int spair2[npair] = {0};      // indices of the corresponding sites in the environment
    // contact surface normals for each contact pair
    const double csn[npair*3] = {1, 0, 0};
}
// The following constants and types are not changed */
/// Constants
const int NU    = params::nact;         // number of actuated joints
const int NPAIR = params::npair;        // number of contact pairs
const int NV    = NU + 6*params::nfree; // degrees of freedom
const int N     = 2*NV;                 // dimensionality of states
const int M     = NU + NPAIR;           // dimensionality of controls
const int NTS   = params::ncts;         // number of control time
const int NTRAJ = (NTS+1)*N + NTS*M;    // number of trajectory variables
/// Types
// eigen+mujoco types for a time instant
typedef Eigen::Matrix<mjtNum, N, 1>                  stateVec;
typedef Eigen::Matrix<mjtNum, N, N, Eigen::ColMajor> stateDer;
typedef Eigen::Matrix<mjtNum, M, 1>                  ctrlVec;
typedef Eigen::Matrix<mjtNum, N, M, Eigen::ColMajor> ctrlDer;
typedef Eigen::Matrix<mjtNum, NPAIR, 1>              kConVec;
// threaded types for trajectories
typedef std::vector<stateVec, Eigen::aligned_allocator<stateVec>> stateTraj;
typedef std::vector<stateDer, Eigen::aligned_allocator<stateDer>> stateDerTraj;
typedef std::vector<ctrlVec,  Eigen::aligned_allocator<ctrlVec>>  ctrlTraj;
typedef std::vector<ctrlDer,  Eigen::aligned_allocator<ctrlDer>>  ctrlDerTraj;
typedef std::vector<kConVec,  Eigen::aligned_allocator<kConVec>>  KConTraj;
/// Structs
struct trajectory
{
    stateTraj X;
    ctrlTraj  U;
    stateDerTraj Fx;
    ctrlDerTraj  Fu;
};

#endif //CITO_PARAMS_H
