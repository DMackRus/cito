/*! Parameters */
/**
 *  \brief CitoParams contains definitions that are used across classes.
 *
 *  This header contains global types, structs, and paths as well as the
 *  CitoParams class that parses the model and config files and defines
 *  parameters that are used across classes.
 *
 *  \author Aykut Onol
 */

#ifndef CITO_PARAMS_H
#define CITO_PARAMS_H

#include <iostream>
#include <string>
#include <chrono>

#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>
#include "mujoco.h"

/// Types
typedef Eigen::VectorXd eigVd;
typedef Eigen::MatrixXd eigMd;
typedef std::vector<eigMd> eigTd;

/// Structs
struct trajectory
{
    eigMd X;
    eigMd U;
    eigTd Fx;
    eigTd Fu;
};

/// Paths
namespace paths
{
    const std::string workspaceDir = std::getenv("CITO_WS");
}

/// CitoParams class
class CitoParams{
public:
    /// Constructor
    CitoParams(const mjModel *model);
    /// Destructor
    ~CitoParams();
    /// Simulation and model parameters
    const mjModel *model;
    double tf, tc, dt;
    int N, ndpc, nu, nv, n, m, nTraj, nPair,
        nFree, *pFree, *bFree, *dAct,
        *quatAdr, *dofAdr;
    Eigen::VectorXi sPair1, sPair2;
    eigMd nCS;
    /// Task parameters
    eigVd desiredPos;
    int controlJointDOF0;
    double weight[4];
};

#endif //CITO_PARAMS_H
