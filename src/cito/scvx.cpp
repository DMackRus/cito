// ***** DESCRIPTION ***********************************************************
// SCVX class defines functions that are used to roll-out the dynamics,
// evaluate the cost, and execute the SCVX algorithm.

#include "cito/scvx.h"

// ***** CONSTRUCTOR ***********************************************************
SCVX::SCVX(const mjModel *m_, Params *cp_, Control *cc_) : m(m_), cp(cp_), cc(cc_), nd(m_, cp_, cc_), sq(m_, cp_)
{
    // initialize Eigen variables
    finalPos.resize(6);
    // get SCVX parameters
    YAML::Node paramSCVX = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/scvx.yaml");
    beta_expand = paramSCVX["beta_expand"].as<double>();
    beta_shrink = paramSCVX["beta_shrink"].as<double>();
    maxIter = paramSCVX["maxIter"].as<int>();
    dLTol = paramSCVX["dLTol"].as<double>();
    rho0 = paramSCVX["rho0"].as<double>();
    rho1 = paramSCVX["rho1"].as<double>();
    rho2 = paramSCVX["rho2"].as<double>();
    rMin = paramSCVX["rMin"].as<double>();
    rMax = paramSCVX["rMax"].as<double>();
    r0 = paramSCVX["r0"].as<double>();
    J = new double[maxIter + 1];
    JTemp = new double[maxIter + 1];
    JTilde = new double[maxIter + 1];
    dJ = new double[maxIter + 1];
    dL = new double[maxIter + 1];
    rho = new double[maxIter + 1];
    r = new double[maxIter + 1];
    accept = new bool[maxIter];
    time_derivs = new double[maxIter];
    time_qp = new double[maxIter];
    // set initial trust-region radius
    r[0] = r0;
    // set bounds
    cc->getBounds();
    // resize trajectories
    XSucc.resize(cp->n, cp->N + 1);
    dX.resize(cp->n, cp->N + 1);
    XTilde.resize(cp->n, cp->N + 1);
    USucc.resize(cp->m, cp->N);
    UTemp.resize(cp->m, cp->N);
    dU.resize(cp->m, cp->N);
    Fx.resize(cp->N);
    Fu.resize(cp->N);
    for (int i = 0; i < cp->N; i++)
    {
        Fx[i].resize(cp->n, cp->n);
        Fu[i].resize(cp->n, cp->m);
    }
}
// ***** DESTRUCTOR ************************************************************
SCVX::~SCVX()
{
    delete[] J;
    delete[] JTemp;
    delete[] JTilde;
    delete[] dJ;
    delete[] dL;
    delete[] rho;
    delete[] r;
    delete[] accept;
}

// ***** FUNCTIONS *************************************************************
// getCost: returns the nonlinear cost given control trajectory and final state
double SCVX::getCost(const eigMd &X, const eigMd &U)
{
    // final cost
    finalPos = X.col(cp->N).segment(cp->controlJointDOF0, 6);
    Jf = 0.5 * (cp->weight[0] * (cp->desiredPos.head(2) - finalPos.head(2)).squaredNorm() +
                cp->weight[1] * (cp->desiredPos.tail(4) - finalPos.tail(4)).squaredNorm());
    // integrated cost
    Ji = cp->weight[2] * X.leftCols(cp->N).bottomRows(m->nv).squaredNorm() +
         cp->weight[3] * U.bottomRows(cp->nPair).sum();
    // total cost
    Jt = Jf + Ji;
    return Jt;
}

// runSimulation: rolls-out and linearizes the dynamics given control trajectory
trajectory SCVX::runSimulation(const eigMd &U, bool linearize, int save, double compensateBias)
{
    // make mjData
    mjData *d = NULL;
    d = mj_makeData(m);
    // initialize d
    mju_copy(d->qpos, m->key_qpos, m->nq);
    mj_forward(m, d);
    cc->setControl(d, U.col(0), compensateBias);

//    std::cout << "before rollout occurs \n";

    // Rollout the dynamics and save data to buffer
    nd.save_data_to_rollout_data(d, 0);
    for (int i = 0; i < cp->N; i++)
    {
        mj_forward(m, d);
        // get the current state values
        XSucc.col(i).setZero();
        XSucc.col(i) = cp->getState(d);
        // take tc/dt steps
        cc->takeStep(d, U.col(i), save, compensateBias);

        // save data to buffer
        nd.save_data_to_rollout_data(d, i + 1);
    }

//    std::cout << "rollout and save occured \n";


    std::vector<std::vector<int>> keypoints;
    std::vector<double> jerkThresholds = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    std::vector<double> velChange_thresholds = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    derivative_interpolator interpolator = {"iterative_error", 1, 5, jerkThresholds, velChange_thresholds, 0.1};
    if(linearize){

        if(interpolator.keyPoint_method == "iterative_error"){
            computedKeyPoints.clear();
            keypoints = generateKeypointsIterativeError(interpolator, cp->N, U, compensateBias);
        }
        else{
            keypoints = nd.generateKeypoints(interpolator, XSucc, cp->N);
        }

        for(int i = 0; i < keypoints.size(); i++){
            std::cout << i << " : ";
            for(int j = 0; j < keypoints[i].size(); j++){
                std::cout << keypoints[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Compute the derivatives, but not for iterative error keypoint method as that has already computed
    // the derivatives
    if(linearize && interpolator.keyPoint_method != "iterative_error"){
        for (int i = 0; i < cp->N; i++){
            Fx[i].setZero();
            Fu[i].setZero();
            // Only linearize the dynamics when a keypoint is reached
            std::vector<int> cols = keypoints[i];
            if(cols.size() > 0){
                nd.linDyn(nd.rollout_data[i], U.col(i), Fx[i].data(), Fu[i].data(), compensateBias, cols);
            }
        }
    }

    // rollout (and linearize) the dynamics
//    for (int i = 0; i < cp->N; i++)
//    {
//        mj_forward(m, d);
//        // get the current state values
//        XSucc.col(i).setZero();
//        XSucc.col(i) = cp->getState(d);
//        // linearization
//        if (linearize)
//        {
//            Fx[i].setZero();
//            Fu[i].setZero();
//            // Only linearize the dynamics when a keypoint is reached
//            std::vector<int> cols = keypoints[i];
//            if(cols.size() > 0){
//                nd.linDyn(d, U.col(i), Fx[i].data(), Fu[i].data(), compensateBias, cols);
//
//            }
//        }
//        // take tc/dt steps
//        cc->takeStep(d, U.col(i), save, compensateBias);
//    }

//    std::cout << Fx[0] << std::endl;
//    std::cout << Fu[0] << std::endl;

    if(linearize){
        if(!(interpolator.keyPoint_method == "set_interval" && interpolator.min_n == 1)){
            nd.interpolateDerivs(keypoints, Fx, Fu, cp->N);
        }
    }

    XSucc.col(cp->N).setZero();
    XSucc.col(cp->N) = cp->getState(d);
    // delete data
    mj_deleteData(d);
    // build trajectory
    traj.X = XSucc;
    traj.U = U;
    if (linearize)
    {
        traj.Fx = Fx;
        traj.Fu = Fu;
    }
    return traj;
}

// solveSCVX: executes the successive convexification algorithm
eigMd SCVX::solveSCVX(const eigMd &U0)
{
    auto optStart = std::chrono::system_clock::now();
    // initialize USucc for the first succession
    USucc = U0;
    // start the SCVX algorithm
    int iter = 0;
    while (!stop)
    {
        std::cout << "Iteration " << iter + 1 << ":" << '\n';
        // simulation and convexification ======================================
        if (iter == 0 || accept[iter - 1])
        {
            std::cout << "INFO: convexification in progress\n";
            auto tDiffStart = std::chrono::system_clock::now();
            trajS = {};
            trajS = this->runSimulation(USucc, true, 0, 1);
            auto tDiffEnd = std::chrono::system_clock::now();
            auto tDiff = std::chrono::duration<double>(tDiffEnd - tDiffStart).count();
            time_derivs[iter] = tDiff;
            std::cout << "INFO: convexification took " << tDiff << " s \n";
        }
        else{
            time_derivs[iter] = 0.0f;
        }
        // get the nonlinear cost if the first iteration
        if (iter == 0)
        {
            J[iter] = this->getCost(trajS.X, USucc);
        }
        // convex optimization =================================================
        double *dTraj = new double[cp->nTraj];
        std::cout << "INFO: QP solver in progress\n\n";
        auto tQPStart = std::chrono::system_clock::now();
        sq.solveCvx(dTraj, r[iter], trajS.X, USucc, trajS.Fx, trajS.Fu, cc->isJFree, cc->isAFree,
                    cc->qposLB, cc->qposUB, cc->tauLB, cc->tauUB);
        auto tQPEnd = std::chrono::system_clock::now();
        auto tQP = std::chrono::duration<double>(tQPEnd - tQPStart).count();
        std::cout << "\nINFO: QP solver took " << tQP << " s \n\n";
        time_qp[iter] = tQP;

        // apply the change
        for (int i = 0; i < cp->N + 1; i++)
        {
            // states
            dX.col(i).setZero();
            XTilde.col(i).setZero();
            for (int j = 0; j < cp->n; j++)
            {
                dX.col(i)[j] = dTraj[i * cp->n + j];
            }
            XTilde.col(i) = trajS.X.col(i) + dX.col(i);
            // controls
            if (i < cp->N)
            {
                dU.col(i).setZero();
                UTemp.col(i).setZero();
                for (int j = 0; j < cp->m; j++)
                {
                    dU.col(i)[j] = dTraj[(cp->N + 1) * cp->n + i * cp->m + j];
                }
                UTemp.col(i) = USucc.col(i) + dU.col(i);
            }
        }
        // evaluate the dynamics for the change and get the cost values ========
        trajTemp = {};
        trajTemp = this->runSimulation(UTemp, false, 0, 1);
        // get the linear and nonlinear costs
        JTilde[iter] = this->getCost(XTilde, UTemp);
        JTemp[iter] = this->getCost(trajTemp.X, UTemp);
        // similarity measure ==================================================
        dJ[iter] = J[iter] - JTemp[iter];
        dL[iter] = J[iter] - JTilde[iter];
        rho[iter] = dJ[iter] / dL[iter];
        if (fabs(dL[iter]) < dLTol)
        {
            dLTolMet = 1;
        }
        // accept or reject the solution =======================================
        // reject
        if (rho[iter] <= rho0 || (dL[iter] < 0 && dJ[iter] < 0))
        {
            accept[iter] = false;
            r[iter + 1] = r[iter] / beta_shrink;
            J[iter + 1] = J[iter];
        }
        else
        {
            accept[iter] = true;
        }
        // accept
        if (accept[iter])
        {
            J[iter + 1] = JTemp[iter];
            USucc = UTemp;
            if (rho[iter] < rho1)
            {
                r[iter + 1] = r[iter] / beta_shrink;
            }
            else if (rho[iter] >= rho1 && rho[iter] < rho2)
            {
                r[iter + 1] = r[iter];
            }
            else if (rho[iter] >= rho2)
            {
                r[iter + 1] = r[iter] * beta_expand;
            }
        }
        // bound the trust region radius r
        r[iter + 1] = std::max(r[iter + 1], rMin);
        r[iter + 1] = std::min(r[iter + 1], rMax);
        // stopping criteria check =============================================
        if (iter + 1 == maxIter)
        {
            stop = true;
            std::cout << "\n\n\tINFO: Maximum number of iterations reached.\n\n";
        }
        if (dLTolMet)
        {
            stop = true;
            std::cout << "\n\n\tINFO: |dL| = |" << dL[iter] << "| < dLTol = " << dLTol << "\n\n";
        }
        // screen output for the iteration =====================================
        std::cout << "Actual:\nFinal pos: " << trajTemp.X.col(cp->N).head(m->nv).transpose() << "\n";
        std::cout << "Final vel: " << trajTemp.X.col(cp->N).tail(m->nv).transpose() << "\n";
        std::cout << "Predicted:\nFinal pos: " << XTilde.col(cp->N).head(m->nv).transpose() << "\n";
        std::cout << "Final vel: " << XTilde.col(cp->N).tail(m->nv).transpose() << "\n";
        std::cout << "L = " << JTilde[iter] << ", J = " << JTemp[iter] << ", kmax = " << trajTemp.U.bottomRows(cp->nPair).maxCoeff() << ", kavg = " << trajTemp.U.bottomRows(cp->nPair).sum() / (cp->nPair * cp->N) << "\n\n\n";
        // next iteration ======================================================
        iter++;
        delete[] dTraj;
    }
    // summary screen output ===============================================
    std::cout << "\n\nSCVX Summary\nJ0=" << J[0] << "\n\n";
    for (int i = 0; i < iter; i++)
    {
        if (i % 10 == 0)
        {
            printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n",
                   "Iteration", "L", "J", "dL", "dJ", "rho", "r", "accept", "derivs time", "qp time");
        }
        printf("%-12d%-12.6g%-12.6g%-12.3g%-12.3g%-12.3g%-12.3g%-12d%-12.6g%-12.6g\n",
               i + 1, JTilde[i], JTemp[i], dL[i], dJ[i], rho[i], r[i], accept[i], time_derivs[i], time_qp[i]);
    }

    optTime = 0.0f;
    derivsTime = 0.0f;
    qpTime = 0.0f;
    costReduction = 1 - (J[iter - 1]) / J[0];

    for(int i = 0; i < iter; i++){
        derivsTime += time_derivs[i];
        qpTime += time_qp[i];
        optTime += time_derivs[i] + time_qp[i];
    }
    std::cout << "\n\nOptimization time: " << optTime << " seconds\n\n";
    std::cout << "Derivatives time: " << derivsTime << " seconds\n\n";
    std::cout << "QP time: " << qpTime << " seconds\n\n";
    return USucc;
}

std::vector<std::vector<int>> SCVX::generateKeypointsIterativeError(derivative_interpolator di, int horizon, const eigMd U, double compensateBias){
    int dof = cp->n / 2;
    std::vector<std::vector<int>> keypoints;
    bool binsComplete[dof];
    std::vector<indexTuple> indexTuples;
    int startIndex = 0;
    int endIndex = horizon - 1;

    // Initialise variables
    for(int i = 0; i < dof; i++){
        binsComplete[i] = false;
        computedKeyPoints.push_back(std::vector<int>());
    }

    for(int i = 0; i < horizon; i++){
        keypoints.push_back(std::vector<int>());
    }

    // Loop through all dofs in the system
//#pragma omp parallel for
    for(int i = 0; i < dof; i++){
        std::vector<indexTuple> listOfIndicesCheck;
        indexTuple initialTuple;
        initialTuple.startIndex = startIndex;
        initialTuple.endIndex = endIndex;
        listOfIndicesCheck.push_back(initialTuple);

        std::vector<indexTuple> subListIndices;
        std::vector<int> subListWithMidpoints;

        while(!binsComplete[i]){
            bool allChecksComplete = true;

            for(int j = 0; j < listOfIndicesCheck.size(); j++) {

                int midIndex = (listOfIndicesCheck[j].startIndex + listOfIndicesCheck[j].endIndex) / 2;
//                cout <<"dof: " << i <<  ": index tuple: " << listOfIndicesCheck[j].startIndex << " " << listOfIndicesCheck[j].endIndex << endl;
                bool approximationGood = checkDoFColumnError(di, listOfIndicesCheck[j], i, U, compensateBias);

                if (!approximationGood) {
                    allChecksComplete = false;
                    indexTuple tuple1;
                    tuple1.startIndex = listOfIndicesCheck[j].startIndex;
                    tuple1.endIndex = midIndex;
                    indexTuple tuple2;
                    tuple2.startIndex = midIndex;
                    tuple2.endIndex = listOfIndicesCheck[j].endIndex;
                    subListIndices.push_back(tuple1);
                    subListIndices.push_back(tuple2);
                }
                else{
                    subListWithMidpoints.push_back(listOfIndicesCheck[j].startIndex);
                    subListWithMidpoints.push_back(midIndex);
                    subListWithMidpoints.push_back(listOfIndicesCheck[j].endIndex);
                }
            }

            if(allChecksComplete){
                binsComplete[i] = true;
                subListWithMidpoints.clear();
            }

            listOfIndicesCheck = subListIndices;
            subListIndices.clear();
        }
    }

    // Loop over the horizon
    for(int i = 0; i < horizon; i++){
        // Loop over the dofs
        for(int j = 0; j < dof; j++){
            // Loop over the computed key points per dof
            for(int k = 0; k < computedKeyPoints[j].size(); k++){
                // If the current index is a computed key point
                if(i == computedKeyPoints[j][k]){
                    keypoints[i].push_back(j);
                }
            }
        }
    }

    // Sort list into order
    for(int i = 0; i < horizon; i++){
        std::sort(keypoints[i].begin(), keypoints[i].end());
    }

    // Remove duplicates
    for(int i = 0; i < horizon; i++){
        keypoints[i].erase(std::unique(keypoints[i].begin(), keypoints[i].end()), keypoints[i].end());
    }

    return keypoints;
}

bool SCVX::checkDoFColumnError(derivative_interpolator di, indexTuple indices, int dofIndex, const eigMd U, double compensateBias){
    int dof = cp->n / 2;
    bool approximationGood = false;
    double errorSum = 0.0f;
    int counter = 0;

    Eigen::MatrixXd midColumnsApprox[2];
    for(int i = 0; i < 2; i++){
        midColumnsApprox[i] = Eigen::MatrixXd::Zero(cp->n, 1);
    }

    int midIndex = (indices.startIndex + indices.endIndex) / 2;
    if((indices.endIndex - indices.startIndex) <=  di.min_n){
        return true;
    }

    bool startIndexExists = false;
    bool midIndexExists = false;
    bool endIndexExists = false;

    for(int i = 0; i < computedKeyPoints[dofIndex].size(); i++){
        if(computedKeyPoints[dofIndex][i] == indices.startIndex){
            startIndexExists = true;
        }

        if(computedKeyPoints[dofIndex][i] == midIndex){
            midIndexExists = true;
        }

        if(computedKeyPoints[dofIndex][i] == indices.endIndex){
            endIndexExists = true;
        }
    }

    std::vector<int> cols;
    cols.push_back(dofIndex);

    if(!startIndexExists){
        nd.linDyn(nd.rollout_data[indices.startIndex], U.col(indices.startIndex), Fx[indices.startIndex].data(), Fu[indices.startIndex].data(), compensateBias, cols);
        computedKeyPoints[dofIndex].push_back(indices.startIndex);
    }

    if(!midIndexExists){
        nd.linDyn(nd.rollout_data[midIndex], U.col(midIndex), Fx[midIndex].data(), Fu[midIndex].data(), compensateBias, cols);
        computedKeyPoints[dofIndex].push_back(midIndex);
    }

    if(!endIndexExists){
        nd.linDyn(nd.rollout_data[indices.endIndex], U.col(indices.endIndex), Fx[indices.endIndex].data(), Fu[indices.endIndex].data(), compensateBias, cols);
        computedKeyPoints[dofIndex].push_back(indices.endIndex);
    }

    midColumnsApprox[0] = (Fx[indices.startIndex].block(0, dofIndex, dof*2, 1) + Fx[indices.endIndex].block(0, dofIndex, dof*2, 1)) / 2;
    midColumnsApprox[1] = (Fx[indices.startIndex].block(0, dofIndex + dof, dof*2, 1) + Fx[indices.endIndex].block(0, dofIndex + dof, dof*2, 1)) / 2;

    for(int i = 0; i < 2; i++){
        int A_col_indices[2] = {dofIndex, dofIndex + dof};
        for(int j = 0; j < cp->n; j++){
            double sqDiff = pow((Fx[midIndex](j, A_col_indices[i]) - midColumnsApprox[i](j, 0)),2);

            counter++;
            errorSum += sqDiff;
        }
    }

    double averageError;
    if(counter > 0){
        averageError = errorSum / counter;
    }
    else{
        averageError = 0.0f;
    }

    // 0.00005
    if(averageError <  di.error_threshold){ //0.00001
        approximationGood = true;
    }

    return approximationGood;
}

// refresh: refreshes SCVX variables for a new run
void SCVX::refresh()
{
    // delete old variables
    delete[] J;
    delete[] JTemp;
    delete[] JTilde;
    delete[] dJ;
    delete[] dL;
    delete[] rho;
    delete[] r;
    delete[] accept;
    delete[] time_derivs;
    delete[] time_qp;
    // create new variables
    J = new double[maxIter + 1];
    JTemp = new double[maxIter + 1];
    JTilde = new double[maxIter + 1];
    dJ = new double[maxIter + 1];
    dL = new double[maxIter + 1];
    rho = new double[maxIter + 1];
    r = new double[maxIter + 1];
    accept = new bool[maxIter];
    time_derivs = new double[maxIter];
    time_qp = new double[maxIter];
    // set initial trust-region radius
    r[0] = r0;
    // reset flags
    stop = false;
    dLTolMet = false;
}