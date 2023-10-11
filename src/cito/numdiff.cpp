// ***** DESCRIPTION ***********************************************************
// NumDiff defines methods for numerical differentiation of MuJoCo
// dynamics including the forces imposed by the contact model.

#include "cito/numdiff.h"

// ***** CONSTRUCTOR ***********************************************************
NumDiff::NumDiff(const mjModel *m_, Params *cp_, Control *cc_) : m(m_), cp(cp_), cc(cc_)
{
    // initialize Eigen variables
    xNewTemp.resize(cp->n);
    xNewP.resize(cp->n);
    xNewN.resize(cp->n);
    uTemp.resize(cp->m);

    for(int i = 0; i < cp->N + 1; i++)
    {
        rollout_data.push_back(mj_makeData(m));
    }

    std::cout << "rollout data made from model \n";
}

// ***** FUNCTIONS *************************************************************
// copyTakeStep: sets xNew to the integration of data given a control input
void NumDiff::copyTakeStep(const mjData *dMain, const eigVd &u, double *xNew, double compensateBias)
{
    // create new data
    mjData *d;
    d = mj_makeData(m);
    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, m->nv);
    mju_copy(d->qacc, dMain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6 * m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);
    // run full computation at center point (usually faster than copying dMain)
    mj_forward(m, d);
    cc->setControl(d, u, compensateBias);
    // take a full control step (i.e., tc/dt steps)
    cc->takeStep(d, u, 0, compensateBias);
    // get new state
    xNewTemp.setZero();
    xNewTemp = cp->getState(d);
    mju_copy(xNew, xNewTemp.data(), cp->n);
    // delete data
    mj_deleteData(d);
}

// hardWorker: for full, slow finite-difference computation
void NumDiff::hardWorker(const mjData *dMain, const eigVd &uMain, double *deriv, double compensateBias, std::vector<int> cols)
{
    // create data
    mjData *d;
    d = mj_makeData(m);
    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, m->nv);
    mju_copy(d->qacc, dMain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6 * m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);
    // finite-difference over positions

    for (int i = 0; i < m->nv; i++)
    {
        bool computeCol = false;
        for(int j = 0; j < cols.size(); j++){
            if(i == cols[j]){
                computeCol = true;
            }
        }

        if(!computeCol){
            continue;
        }

        // get joint id for this dof
        int jID = m->dof_jntid[i];
        // apply quaternion or simple perturbation
        if (cp->quatAdr[i] >= 0)
        {
            mjtNum angvel[3] = {0, 0, 0};
            angvel[cp->dofAdr[i]] = eps;
            mju_quatIntegrate(d->qpos + cp->quatAdr[i], angvel, 1);
        }
        else
        {
            d->qpos[m->jnt_qposadr[jID] + i - m->jnt_dofadr[jID]] += eps;
        }
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uMain, xNewP.data(), compensateBias);
        // undo perturbation
        mju_copy(d->qpos, dMain->qpos, m->nq);
        // apply quaternion or simple perturbation
        if (cp->quatAdr[i] >= 0)
        {
            mjtNum angvel[3] = {0, 0, 0};
            angvel[cp->dofAdr[i]] = -eps;
            mju_quatIntegrate(d->qpos + cp->quatAdr[i], angvel, 1);
        }
        else
        {
            d->qpos[m->jnt_qposadr[jID] + i - m->jnt_dofadr[jID]] -= eps;
        }
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uMain, xNewN.data(), compensateBias);
        // undo perturbation
        mju_copy(d->qpos, dMain->qpos, m->nq);
        // compute column i of dx/dqpos
        for (int j = 0; j < cp->n; j++)
        {
            deriv[i * cp->n + j] = (xNewP(j) - xNewN(j)) / (2 * eps);
        }
    }
    // finite-difference over velocities
    for (int i = 0; i < m->nv; i++)
    {
        bool computeCol = false;
        for(int j = 0; j < cols.size(); j++){
            if(i == cols[j]){
                computeCol = true;
            }
        }

        if(!computeCol){
            continue;
        }

        // perturb velocity
        d->qvel[i] += eps;
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uMain, xNewP.data(), compensateBias);
        // perturb velocity
        d->qvel[i] = dMain->qvel[i] - eps;
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uMain, xNewN.data(), compensateBias);
        // undo perturbation
        d->qvel[i] = dMain->qvel[i];
        // compute column i of dx/dqvel
        for (int j = 0; j < cp->n; j++)
        {
            deriv[cp->n * m->nv + i * cp->n + j] = (xNewP(j) - xNewN(j)) / (2 * eps);
        }
    }
    // finite-difference over control variables
    // copy uMain to uTemp for perturbations
    uTemp = uMain;
    for (int i = 0; i < cp->m; i++)
    {
        bool computeCol = false;
        for(int j = 0; j < cols.size(); j++){
            if(i == cols[j]){
                computeCol = true;
            }
        }

        if(!computeCol){
            continue;
        }

        // perturbation in the positive direction
        uTemp(i) += eps;
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uTemp, xNewP.data(), compensateBias);
        // perturbation in the negative direction
        uTemp(i) -= 2 * eps;
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uTemp, xNewN.data(), compensateBias);
        // compute column i of dx/du
        for (int j = 0; j < cp->n; j++)
        {
            deriv[cp->n * cp->n + i * cp->n + j] = (xNewP(j) - xNewN(j)) / (2 * eps);
        }
    }
    // delete data
    mj_deleteData(d);
}

void NumDiff::save_data_to_rollout_data(mjData *dmain, int index){

    rollout_data[index]->time = dmain->time;
    mju_copy(rollout_data[index]->qpos, dmain->qpos, m->nq);
    mju_copy(rollout_data[index]->qvel, dmain->qvel, m->nv);
    mju_copy(rollout_data[index]->qacc, dmain->qacc, m->nv);
    mju_copy(rollout_data[index]->qacc_warmstart, dmain->qacc_warmstart, m->nv);
    mju_copy(rollout_data[index]->qfrc_applied, dmain->qfrc_applied, m->nv);
    mju_copy(rollout_data[index]->xfrc_applied, dmain->xfrc_applied, 6 * m->nbody);
    mju_copy(rollout_data[index]->ctrl, dmain->ctrl, m->nu);
}

// linDyn: calculates derivatives of the state and control trajectories
void NumDiff::linDyn(const mjData *dMain, const eigVd &uMain, double *Fxd, double *Fud, double compensateBias, std::vector<int> cols)
{
    // TODO: consider doing the memory allocation/freeing in the constructor/destructor
    double *deriv = (double *)mju_malloc(sizeof(double) * cp->n * (cp->n + cp->m));
    this->hardWorker(dMain, uMain, deriv, compensateBias, cols);
    mju_copy(Fxd, deriv, cp->n * cp->n);
    mju_copy(Fud, deriv + cp->n * cp->n, cp->n * cp->m);
    mju_free(deriv);
}

std::vector<std::vector<int>> NumDiff::generateKeypoints(derivative_interpolator di, const eigMd X, int horizon){
    std::vector<std::vector<int>> keypoints;
    std::vector<int> fullRow;
    std::vector<int> emptyRow;

    int dof = X.rows() / 2;

    if(di.keyPoint_method == "set_interval"){
        for(int i = 0; i < dof; i++){
            fullRow.push_back(i);
        }
        keypoints.push_back(fullRow);


        for(int i = 1; i < horizon; i++){
            if(i % di.min_n == 0){
                keypoints.push_back(fullRow);
            }
            else{
                keypoints.push_back(emptyRow);
            }
        }

    }
    else if(di.keyPoint_method == "adaptive_jerk"){
        eigMd jerk_profile = getJerkProfile(X, horizon);
        keypoints = generateKeypointsAdaptiveJerk(di, jerk_profile, horizon);

    }
    else if(di.keyPoint_method == "magvel_change"){
        keypoints = generateKeypointsMagnitudeVelChange(di, X, horizon);

    }
    else if(di.keyPoint_method == "iterative_error"){
//        keypoints = generateKeypointsIterativeError(di, horizon);
        std::cout << "method iterative error recognized, but it shouldn't have reached this function" << std::endl;
    }
    else{
        std::cout << "keyPoint_method not recognized" << std::endl;
    }

    return keypoints;
}

std::vector<std::vector<int>> NumDiff::generateKeypointsAdaptiveJerk(derivative_interpolator di, const eigMd jerk_profile, int horizon){
    std::vector<std::vector<int>> keypoints;
    int dof = jerk_profile.rows();

    for(int t = 0; t < horizon; t++){
        keypoints.push_back(std::vector<int>());
    }

    // Enforce that first keypoint, all dofs are computed
    for(int i = 0; i < dof; i++){
        keypoints[0].push_back(i);
    }

    int lastIndices[dof];
    for(int i = 0; i < dof; i++){
        lastIndices[i] = 0;
    }

    // Loop over the trajectory
    for(int j = 0; j < dof; j++){
        for(int i = 0; i < jerk_profile.cols(); i++){

            if((i - lastIndices[j]) >= di.min_n) {
                // Check if the jerk is above the threshold
                if (jerk_profile(i, j) > di.jerk_thresholds[j]) {
                    keypoints[i].push_back(j);
                    lastIndices[j] = i;
                }
            }

            if((i - lastIndices[j]) >= di.max_n){
                keypoints[i].push_back(j);
                lastIndices[j] = i;
            }
        }
    }

    // Enforce that the last time-step, all dofs are computed
    for(int i = 0; i < dof; i++){
        keypoints.back().push_back(i);
    }

    return keypoints;
}

eigMd NumDiff::getJerkProfile(const eigMd X, int horizon){
    int dof = cp->n / 2;
    eigMd jerk_profile = eigMd::Zero(dof, horizon - 2);

    for(int i = 0; i < horizon - 2; i++){
        for(int j = 0; j < dof; j++){
            jerk_profile(j, i) = X(j + dof, i + 2) - 2 * X(j + dof, i + 1) + X(j + dof, i);
        }
    }

    return jerk_profile;
}

std::vector<std::vector<int>> NumDiff::generateKeypointsMagnitudeVelChange(derivative_interpolator di, const eigMd X, int horizon){
    std::vector<std::vector<int>> keypoints;

    int dof = cp->n / 2;

    for(int t = 0; t < horizon; t++){
        keypoints.push_back(std::vector<int>());
    }

    for(int i = 0; i < dof; i++){
        keypoints[0].push_back(i);
    }

    // Keeps track of interval from last keypoint for this dof
    std::vector<int> lastKeypointCounter = std::vector<int>(dof, 0);
    std::vector<double> lastVelValue = std::vector<double>(dof, 0);
    std::vector<double> lastVelDirection = std::vector<double>(dof, 0);

    for(int i = 0; i < dof; i++){
        lastVelValue[i] = X(i + dof, 0);
    }

    // Loop over the velocity dofs
    for(int i = 0; i < dof; i++){
        // Loop over the horizon
        for(int t = 1; t < horizon; t++){

            lastKeypointCounter[i]++;
            double currentVelDirection = X(i + dof, t) - X(i + dof, t - 1);
            double currentVelChangeSinceKeypoint = X(i + dof, t) - lastVelValue[i];

            // If the vel change is above the required threshold
            if(lastKeypointCounter[i] >= di.min_n){
                if(abs(currentVelChangeSinceKeypoint) > di.velChange_thresholds[i]){
                    keypoints[t].push_back(i);
                    lastVelValue[i] = X(i + dof, t);
                    lastKeypointCounter[i] = 0;
                    continue;
                }
            }

            // If the interval is greater than minN
            if(lastKeypointCounter[i] >= di.min_n){
                // If the direction of the velocity has changed
                if(currentVelDirection * lastVelDirection[i] < 0){
                    keypoints[t].push_back(i);
                    lastVelValue[i] = X(i + dof, t);
                    lastKeypointCounter[i] = 0;
                    continue;
                }
            }
            else{
                lastVelDirection[i] = currentVelDirection;
            }

            // If interval is greater than maxN
            if(lastKeypointCounter[i] >= di.max_n){
                keypoints[t].push_back(i);
                lastVelValue[i] = X(i + dof, t);
                lastKeypointCounter[i] = 0;
                continue;
            }
        }
    }

    // Enforce last keypoint for all dofs at horizonLength - 1
    for(int i = 0; i < dof; i++){
        keypoints.back().push_back(i);
    }

    return keypoints;
}

void NumDiff::interpolateDerivs(std::vector<std::vector<int>> keypoints, eigTd &Fxd, eigTd &Fud, int horizon){
    int dof = Fxd[0].rows() / 2;
    int num_ctrl = Fud[0].cols();

    // Declare memory of interpolation variables
    eigMd Fxd_start_col1;
    eigMd Fxd_end_col1;
    eigMd Fxd_add_col1;

    eigMd Fxd_start_col2;
    eigMd Fxd_end_col2;
    eigMd Fxd_add_col2;

    eigMd Fud_start_col;
    eigMd Fud_end_col;
    eigMd Fud_add_col;


    int start_indices[dof];
    for(int i = 0; i < dof; i++){
        start_indices[i] = 0;
    }

    for(int t = 1; t < horizon; t++) {
        for(int i = 0; i < dof; i++){
            std::vector<int> columns = keypoints[t];

            if(columns.size() == 0){
                continue;
            }

            for(int j = 0; j < columns.size(); j++){

                if(i == columns[j]){
                    Fxd_start_col1 = Fxd[start_indices[i]].block(0, i, (2*dof), 1);
                    Fxd_end_col1 = Fxd[t].block(0, i, (2*dof), 1);
                    Fxd_add_col1 = (Fxd_end_col1 - Fxd_start_col1) / (t - start_indices[i]);

                    // Same again for column 2 which is dof + i
                    Fxd_start_col2 = Fxd[start_indices[i]].block(0, dof + i, (2*dof), 1);
                    Fxd_end_col2 = Fxd[t].block(0, dof + i, (2*dof), 1);
                    Fxd_add_col2 = (Fxd_end_col2 - Fxd_start_col2) / (t - start_indices[i]);

                    if(i < num_ctrl){
                        Fud_start_col = Fud[start_indices[i]].block(0, i, (2*dof), 1);
                        Fud_end_col = Fud[t].block(0, i, (2*dof), 1);
                        Fud_add_col = (Fud_end_col - Fud_start_col) / (t - start_indices[i]);
                    }

                    // Loop over the indices between the start and end
                    for(int k = start_indices[i] + 1; k < t; k++){
                        Fxd[k].block(0, i, (2*dof), 1) = Fxd_start_col1 + Fxd_add_col1 * (k - start_indices[i]);
                        Fxd[k].block(0, dof + i, (2*dof), 1) = Fxd_start_col2 + Fxd_add_col2 * (k - start_indices[i]);
                        if(i < num_ctrl){
                            Fud[k].block(0, i, (2*dof), 1) = Fud_start_col + Fud_add_col * (k - start_indices[i]);
                        }
                    }

                    start_indices[i] = t;
                }
            }
        }
    }

    // Loop over the keypoints
//    for(int i = 0; i < num_keypoints - 1; i++){
//        int start_index = keypoints[i];
//        int end_index = keypoints[i+1];
//        int interval = end_index - start_index;
//
//        eigMd Fxd_start = Fxd[start_index];
//        eigMd Fxd_end = Fxd[end_index];
//        eigMd Fxd_add = (Fxd_end - Fxd_start) / interval;
//
//        eigMd Fud_start = Fud[start_index];
//        eigMd Fud_end = Fud[end_index];
//        eigMd Fud_add = (Fud_end - Fud_start) / interval;
//
//        for(int j = 1; j < interval; j++){
//            Fxd[start_index + j] = Fxd_start + Fxd_add * j;
//            Fud[start_index + j] = Fud_start + Fud_add * j;
//        }
//    }

}

void NumDiff::saveLinearisation(const std::string file_prefix, eigTd Fxd, eigTd Fud,  eigMd X, eigMd U, int horizon){
    std::string projectParentPath = __FILE__;
    projectParentPath = projectParentPath.substr(0, projectParentPath.find_last_of("/\\"));
    projectParentPath = projectParentPath.substr(0, projectParentPath.find_last_of("/\\"));
    projectParentPath = projectParentPath.substr(0, projectParentPath.find_last_of("/\\"));
    std::string rootPath = projectParentPath + "/savedTrajecInfo/" + file_prefix;

    std::string filename = rootPath + "/A_matrices.csv";
    std::cout << "filename: " << filename << std::endl;
    std::ofstream fileOutput;
    fileOutput.open(filename);

    int dof = Fxd[0].rows() / 2;
    int num_ctrl = Fud[0].cols();

    // Save the fx matrices
    for(int i = 0; i < horizon - 1; i++){
        // Row
        for(int j = 0; j < (dof); j++){
            // Column
            for(int k = 0; k < (2 * dof); k++){
                fileOutput << Fxd[i](j + dof, k) << ",";
            }
        }
        fileOutput << std::endl;
    }

    fileOutput.close();

    filename = rootPath + "/B_matrices.csv";
    std::cout << "filename: " << filename << std::endl;
    fileOutput.open(filename);

    // Save the fu matrices
    for(int i = 0; i < horizon - 1; i++){
        // Row
        for(int j = 0; j < (dof); j++){
            // Column
            for(int k = 0; k < num_ctrl; k++){
                fileOutput << Fud[i](j + dof, k) << ",";
            }
        }
        fileOutput << std::endl;
    }

    fileOutput.close();

    std::cout << "X: " << X << std::endl;

    filename = rootPath + "/states.csv";
    std::cout << "filename: " << filename << std::endl;
    fileOutput.open(filename);

    // Save the fu matrices
    for(int i = 0; i < horizon - 1; i++){
        // Row
        for(int j = 0; j < 2 * dof; j++){
            // Column
            fileOutput << X(j, i) << ",";
        }
        fileOutput << std::endl;
    }

    fileOutput.close();

    std::cout << "U: " << U << std::endl;

    filename = rootPath + "/controls.csv";
    std::cout << "filename: " << filename << std::endl;
    fileOutput.open(filename);

    // Save the fu matrices
    for(int i = 0; i < horizon - 1; i++){
        // Row
        for(int j = 0; j < num_ctrl; j++){
            // Column
            fileOutput << U(j, i) << ",";
        }
        fileOutput << std::endl;
    }

    fileOutput.close();
}
