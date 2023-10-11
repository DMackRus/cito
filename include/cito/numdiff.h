/*! Numerical Differentiation */
/**
 *  \brief NumDiff consists of methods for numerical differentiation.
 *
 *  This class defines functions for numerical differentiation of the MuJoCo
 *  dynamics including the forces imposed by the contact model.
 *
 *  \author Aykut Onol
 */

#ifndef NUMDIFF_H
#define NUMDIFF_H

#include "cito/control.h"

struct derivative_interpolator{
    std::string keyPoint_method;
    int min_n;
    int max_n;
    std::vector<double> jerk_thresholds;
    std::vector<double> velChange_thresholds;
    double error_threshold;
};

struct indexTuple{
    int startIndex;
    int endIndex;
};

class NumDiff
{
public:
    /// Constructor
    NumDiff(const mjModel *m_, Params *cp_, Control *cc_);
    /// Destructor
    ~NumDiff() {}
    /// This function calculates derivatives of the state and control trajectories
    void linDyn(const mjData *dMain, const eigVd &uMain, double *Fxd, double *Fud, double compensateBias, std::vector<int> cols);

    void saveLinearisation(const std::string file_name, eigTd Fxd, eigTd Fud, eigMd X, eigMd U, int horizon);

    // Generate keypoints functions
    std::vector<std::vector<int>> generateKeypoints(derivative_interpolator di, const eigMd X, int horizon);
    std::vector<std::vector<int>> generateKeypointsAdaptiveJerk(derivative_interpolator di, const eigMd jerk_profile, int horizon);
    std::vector<std::vector<int>> generateKeypointsMagnitudeVelChange(derivative_interpolator di, const eigMd X, int horizon);
    std::vector<std::vector<int>> generateKeypointsIterativeError(derivative_interpolator di, int horizon);

    // Supporting functions for generating keypoints
    eigMd getJerkProfile(const eigMd X, int horizon);
    bool checkDoFColumnError(derivative_interpolator di, indexTuple indices, int dof);

    void interpolateDerivs(std::vector<std::vector<int>> keypoints, eigTd &Fxd, eigTd &Fud, int horizon);

private:
    /// This function sets xNew to the integration of data given a control input
    void copyTakeStep(const mjData *dMain, const eigVd &u, double *xNew, double compensateBias);
    /// This function calculates central differences by taking full steps
    void hardWorker(const mjData *dMain, const eigVd &uMain, double *deriv, double compensateBias, std::vector<int> cols);
    /// MuJoCo model
    const mjModel *m;
    /// Perturbation for central differences
    double eps = 1e-6;
    /// Variables for differentiation
    eigVd xNewTemp, xNewP, xNewN;
    eigVd uTemp;
    /// Objects
    Params *cp;
    Control *cc;

    // Keypoint variables
    std::vector<std::vector<int>> computedKeyPoints;        // Stores the keypoints computed via iterative error so re-computation isn't performed
    std::vector<mjData*> rollout_data;
};

#endif //NUMDIFF_H
