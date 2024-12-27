/*
 * Volcanic Skies
 * Copyright (C) 2024 P. Cilliers Pretorius, University of Cape Town
 *
 * This file is part of the Volcanic Skies project.
 *
 * Volcanic Skies is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License (GPL) as published by the Free Software 
 * Foundation, either version 2 of the License, or (at your discretion) any later version.
 *
 * Volcanic Skies is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with 
 * Volcanic Skies in the LICENSE file. If not, see <https://www.gnu.org/licenses/>.
 *
 * Additional information and disclaimers regarding liability and third-party 
 * components can be found in the NOTICE file included with this project.
 *
 */
#ifndef PWM_PRESSURE_ENGINE_H
#define PWM_PRESSURE_ENGINE_H

#define peEpsilon
#define solverTolerance 0.1
#define solverIterations 4000

#include "abstractEngine.h"
#include "airLayer.h"
#include "flatStaggeredGrid.h"
#include "square2DArray.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <omp.h>
#include <memory>

namespace PWM{
    namespace Engine{
        template<typename T, typename V> class pressureEngine : public AbstractEngine{
            private:
                //Data structure to hold reference to each of the layers in the simulation
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;

                //Function that calculates the divergence of the layer
                void calcDiv(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);

                void multiGridSolve(T& pressure, const T& divergence, int iteration = 0);

                //Function that calculates the pressures needed to correct for nabla u = 0.
                void solvePressureMatrix(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);

                //Solve the spherical Poisson equation using the Intel Math Kernel Library
                void sphericalPoissonSolver(PWM::Model::airLayer<T, V>& l);

                //Calculate and apply the divergence correction of velocities
                void solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<T, V>>& l, V scaleFactor = 1);

                std::shared_ptr<T> bufferVelTheta;
                std::shared_ptr<T> bufferVelPhi;
                std::shared_ptr<T> divergence;
                Eigen::SparseMatrix<double> A;

                Eigen::VectorXd div;
                Eigen::VectorXd pressure;
                Eigen::VectorXd PrevPressure;
            protected:
                void step_internal() override;
            public:
                pressureEngine(float dt, bool active);
                pressureEngine(const size_t width, const V WrldSizeSize, const float dt, bool active);
                pressureEngine(const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, const float dt, bool active);

                //Calculate divergence and pressure (without solving) to simulate uplift/downdrafts to some extent
                void calcPressures();

                //Calculate divergence and apply a coarse smoothing to the pressure field
                void smoothPressure();

                //Calculate the velocity adjustments resulting from the pressure field
                void solvePressureProjection();

                //Function to add an air layer to the engine
                void addLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);
                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;
        };

        template<typename T, typename V>
        inline pressureEngine<T, V>::pressureEngine(float dt, bool active) : AbstractEngine(dt, active){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>();
            bufferVelPhi = std::make_shared<T>();
            divergence = std::make_shared<T>();
            div = Eigen::VectorXd();
            pressure = Eigen::VectorXd();
            PrevPressure = Eigen::VectorXd();
        }

        template<typename T, typename V>
        inline pressureEngine<T, V>::pressureEngine(const size_t width, const V WrldSize, const float dt, bool active) : AbstractEngine(dt, active){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            bufferVelTheta = std::make_shared<T>(width, WrldSize);
            bufferVelPhi = std::make_shared<T>(width, WrldSize);
            divergence = std::make_shared<T>(width, WrldSize);
        }

        template<typename T, typename V>
        inline pressureEngine<T, V>::pressureEngine(const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, const float dt, bool active) : AbstractEngine(dt, active){
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();
            bufferVelTheta = std::make_shared<T>(width, height, 0.0, 0.5, xWrldSize, yWrldSize);
            bufferVelPhi = std::make_shared<T>(width, height, 0.5, 0.0, xWrldSize, yWrldSize);
            divergence = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);

            div = Eigen::VectorXd(width * height);
            pressure = Eigen::VectorXd(width * height);
            PrevPressure = Eigen::VectorXd(width * height);
            this->A = Eigen::SparseMatrix<double>(width * height, width * height);
            this->A.reserve(Eigen::VectorXd::Constant(width * height, 5));
            for (int i = 0; i < height; ++i){
                for (int j = 0; j < width; ++j){
                    size_t cIdx, aIdx, bIdx, lIdx, rIdx;
                    cIdx = PWM::Utils::convert2Dto1DUtil(height, width, i, j);
                    aIdx = PWM::Utils::convert2Dto1DUtil(height, width, i - 1, j);
                    bIdx = PWM::Utils::convert2Dto1DUtil(height, width, i + 1, j);
                    lIdx = PWM::Utils::convert2Dto1DUtil(height, width, i, j - 1);
                    rIdx = PWM::Utils::convert2Dto1DUtil(height, width, i, j + 1);
                    A.insert(cIdx, cIdx) = 0;
                    A.insert(cIdx, aIdx) = 0;
                    A.insert(cIdx, bIdx) = 0;
                    A.insert(cIdx, lIdx) = 0;
                    A.insert(cIdx, rIdx) = 0;
                }
            }
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::addLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            airLayers.push_back(l);
        }

        template<typename T, typename V>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& pressureEngine<T, V>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::calcDiv(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::calcDiv(PWM::Model::airLayer<T, V>& l) for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<PWM::PWMDataStructure::square2DArray<double>, double>::calcDiv(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>& l){
            int cellNum = l->getObstacles().getWidth();
            double cellSize = 1.0f / cellNum;
            double velSolid = 0;
            #pragma omp parallel for
                for (int i = 0; i < cellNum; ++i){
                    for (int j = 0; j < cellNum; ++j){
                        if (l->getObstacles(i, j)){
                            divergence->setData(i, j, 0);
                            continue;
                        }
                        double velBot, velTop, velRight, velLeft;
                        if (l->getObstacles(i + 1, j))
                            velBot = velSolid;
                        else
                            velBot = l->getVelocityTheta(i + 1, j);
                        if (l->getObstacles(i - 1, j))
                            velTop = velSolid;
                        else
                            velTop = l->getVelocityTheta(i - 1, j);
                        if (l->getObstacles(i, j - 1))
                            velLeft = velSolid;
                        else
                            velLeft = l->getVelocityPhi(i, j - 1);
                        if (l->getObstacles(i, j + 1))
                            velRight = velSolid;
                        else
                            velRight = l->getVelocityPhi(i, j + 1);
                        double div = ((velBot - velTop) + (velRight - velLeft)) / 2.f / cellSize;
                        divergence->setData(i, j, div);
                    }
                }
            #pragma omp barrier
        }

        template<>
        inline void pressureEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::calcDiv(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>& l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            double velSolid = 0;
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        if (l->getObstacles(i, j)){
                            divergence->setData(i, j, 0);
                            continue;
                        }
                        double velBot, velTop, velRight, velLeft;
                        velTop = (l->getObstacles(i - 1, j)) ? velSolid : l->getVelocityTheta(i, j);
                        velBot = (l->getObstacles(i + 1, j)) ? velSolid : l->getVelocityTheta(i + 1, j);
                        velLeft = (l->getObstacles(i, j - 1)) ? velSolid : l->getVelocityPhi(i, j);
                        velRight = (l->getObstacles(i, j + 1)) ? velSolid : l->getVelocityPhi(i, j + 1);
                        double div = (velBot - velTop) / cellSizeY + (velRight - velLeft) / cellSizeX;
                        divergence->setData(i, j, div);
                    }
                }
            #pragma omp barrier
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::multiGridSolve(T& pressure, const T& divergence, int iteration){
            static T* fine_r[MAX_LAYER];
            static T* fine_e[MAX_LAYER];
            static T* coarse_r[MAX_LAYER];
            static T* coarse_e[MAX_LAYER];

            static float initialized = false;
            if (!initialized){
                for (int i = 0; i < MAX_LAYER; ++i){
                    fine_r[i] = nullptr;
                    fine_e[i] = nullptr;
                    coarse_r[i] = nullptr;
                    coarse_e[i] = nullptr;
                }
                initialized = true;
            }
            auto w = pressure.getWidth();
            auto s = pressure.getSize();
            if (! fine_r[iteration])
                fine_r [iteration] = new T(w, s);
            if (! fine_e[iteration])
                fine_e [iteration] = new T(w, s);
            if (! coarse_r[iteration])
                coarse_r [iteration] = new T(w / 2, s);
            if (! coarse_e[iteration])
                coarse_e [iteration] = new T(w / 2, s);

            if (! pressure.compareSizes(*fine_r[iteration]))
                fine_r[iteration]->resize(pressure);

            *fine_r[iteration] = 0;
            *fine_e[iteration] = 0;
            *coarse_r[iteration] = 0;
            *coarse_e[iteration] = 0;

            pressure.gaussSeidelSmoothBSI(divergence, 4);

            //Compute residual
            pressure.laplacian(*fine_r[iteration]);
            fine_r[iteration]->mul(-1);
            fine_r[iteration]->add(divergence);

            if (pressure.getWidth() <= 2){
                coarse_e[iteration]->gaussSeidelSmoothBSI(*coarse_r[iteration], 10, divergence.getWidth());
            }
            else{
                multiGridSolve(*coarse_e[iteration], *coarse_r[iteration], iteration + 1);
            }

            // Interpolate?
            coarse_e[iteration]->expand(*fine_e[iteration]);

            //Apply correction (p = p + e)
            pressure += *fine_e[iteration];

            pressure.gaussSeidelSmoothBSI(divergence, 4);
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::solvePressureMatrix(std::shared_ptr<PWM::Model::airLayer<T, V> > &l){
           std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::solvePressureMatrix() for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::solvePressureMatrix(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double> > &l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            int nCell = nX * nY;
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            double zero = 0, negOne = -1, velSolid = 0;
            auto dt = getDt();
            double rho = PWM::Utils::altitudeAdjustedDensity(l->getHeight(), l->getPlanet());
            double dtRho = dt / rho;;

//            ATest = Eigen::MatrixXd::Zero(nCell, nCell);

            div.setZero(nCell);
            pressure.setZero(nCell);

            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        size_t cIdx, aIdx, bIdx, lIdx, rIdx;
                        cIdx = PWM::Utils::convert2Dto1DUtil(nY, nX, i, j);
                        aIdx = PWM::Utils::convert2Dto1DUtil(nY, nX, i - 1, j);
                        bIdx = PWM::Utils::convert2Dto1DUtil(nY, nX, i + 1, j);
                        lIdx = PWM::Utils::convert2Dto1DUtil(nY, nX, i, j - 1);
                        rIdx = PWM::Utils::convert2Dto1DUtil(nY, nX, i, j + 1);
                        if (l->getObstacles(i, j)){
                            div[cIdx] =  0;

                            //Because Eigen Sparse Matrices assume any values not set explicitly to be 0, we don't have to specify it. Left here for completeness.
                            //size_t aIDX, bIDX, lIDX, rIDX;
                            A.coeffRef(cIdx, aIdx) = 0;
                            A.coeffRef(cIdx, bIdx) = 0;
                            A.coeffRef(cIdx, lIdx) = 0;
                            A.coeffRef(cIdx, rIdx) = 0;
                            A.coeffRef(cIdx, cIdx) = 0;

                            continue;
                        }
                        int pijCoeff = 4;
                        double velBot = 0, velTop = 0, velRight = 0, velLeft = 0;
                        double coeffAbove = 0, coeffBelow = 0, coeffLeft = 0, coeffRight = 0;

                        //Check above the cell for obstacles and adjust formula as needed
                        if (l->getObstacles(i - 1, j)){
                            velTop = velSolid;
                            --pijCoeff;
                            coeffAbove = zero;
                        }
                        else{
                            velTop = l->getVelocityTheta(i, j);
                            coeffAbove = negOne;
                        }

                        //Check below the cell for obstacles and adjust formula as needed
                        if (l->getObstacles(i + 1, j)){
                            velBot = velSolid;
                            --pijCoeff;
                            coeffBelow = zero;
                        }
                        else{
                            velBot = l->getVelocityTheta(i + 1, j);
                            coeffBelow = negOne;
                        }

                        //Check left of the cell for obstacles and adjust formula as needed
                        if (l->getObstacles(i, j - 1)){
                            velLeft = velSolid;
                            --pijCoeff;
                            coeffLeft = zero;
                        }
                        else{
                            velLeft = l->getVelocityPhi(i, j);
                            coeffLeft = negOne;
                        }

                        //Check right of the cell for obstacles and adjust formula as needed
                        if (l->getObstacles(i, j + 1)){
                            velRight = velSolid;
                            --pijCoeff;
                            coeffRight = zero;
                        }
                        else{
                            velRight = l->getVelocityPhi(i, j + 1);
                            coeffRight = negOne;
                        }

                        double dive = (velBot - velTop) / cellSizeY + (velRight - velLeft) / cellSizeX;//Calculate divergence
                        double adjustedDive = (dive / dtRho) * (cellSizeX * cellSizeY);//Isolate coefficients to be either 4/3/2/1 for ij or -1 (or 0) for the neighbours.


                        //Set b value
                        div[cIdx] = -adjustedDive;
                        //Set values in A (fill in relevant values in row ij)

                        A.coeffRef(cIdx, aIdx) = coeffAbove;
                        A.coeffRef(cIdx, bIdx) = coeffBelow;
                        A.coeffRef(cIdx, lIdx) = coeffLeft;
                        A.coeffRef(cIdx, rIdx) = coeffRight;
                        A.coeffRef(cIdx, cIdx) = pijCoeff;

//                        trips.emplace_back(Eigen::Triplet<double>(cidx, aIdx, coeffAbove));
//                        trips.emplace_back(Eigen::Triplet<double>(cidx, bIdx, coeffBelow));
//                        trips.emplace_back(Eigen::Triplet<double>(cidx, lIdx, coeffLeft));
//                        trips.emplace_back(Eigen::Triplet<double>(cidx, rIdx, coeffRight));
//                        trips.emplace_back(Eigen::Triplet<double>(cidx, cidx, pijCoeff));

//                        ATest(cidx, aIdx) = coeffAbove;
//                        ATest(cidx, bIdx) = coeffBelow;
//                        ATest(cidx, lIdx) = coeffLeft;
//                        ATest(cidx, rIdx) = coeffRight;
//                        ATest(cidx, cidx) = pijCoeff;
                    }
                }
             #pragma omp barrier

            A.makeCompressed();
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper/*, Eigen::IncompleteCholesky<double, Eigen::Lower|Eigen::Upper>*/> solver;
            solver.setTolerance(solverTolerance);
            solver.setMaxIterations(solverIterations);
            solver.compute(A);
            if (solver.info() != Eigen::Success)
                std::cerr << "Eigen solver decomposition failed! Solver instead has " << solver.info() << std::endl;

            #pragma omp parallel for
                for (int i = 0; i < nCell; ++i){
                    PrevPressure[i] = l->getPressure(i);
                }
            #pragma omp barrier

            pressure = solver.solveWithGuess(div, PrevPressure);
            if (solver.info() != Eigen::Success){
                std::cerr << "Eigen solver solve failed in layer at height " << l->getHeight() << " metre!" << std::endl;
            }
            if (!(solver.error() < solverTolerance))
                std::cerr << "Solver error: " << solver.error() << std::endl;
//            if (x.isApprox(pressure)){
//                std::cerr << "ColPivHouseholder and CG give same results! Hooray!" << std::endl;
//            }
//            if (div.isApprox(A * pressure)){
//                std::cerr << "ConjugateGradient correctly solves Ap = div in " << solver.iterations() << " iterations." << std::endl;
//            }

            #pragma omp parallel for
                for (int i = 0; i < nCell; ++i){
                    l->setPressure(i, pressure[i]);
                }
            #pragma omp barrier
            #ifdef DEBUG
                std::cerr << "Solver used " << solver.iterations() << " iterations to reach an error of " << solver.error() << "." << std::endl;
                std::cerr << "Pressure has a range of " << l->getPressure().range() << ", a mean of " << l->getPressure().mean() << ", and a standard deviation of " << l->getPressure().stdDev() << std::endl;
            #endif
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::smoothPressure(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::smoothPressure() for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<PWM::PWMDataStructure::square2DArray<double>, double>::smoothPressure(){
            startComputation();
            for (auto l : airLayers){
                calcDiv(l);
                l->getPressure().gaussSeidelSmoothBSI(*divergence, 4);
            }
            endComputation();
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<T, V>>& l, V scaleFactor){
            std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<T, V>>& l, V scaleFactor) for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<PWM::PWMDataStructure::square2DArray<double>, double>::solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>& l, double scaleFactor){
            int cellNum = l->getObstacles().getWidth();
            auto dt = getDt();
            auto rho = PWM::Utils::altitudeAdjustedDensity(l->getHeight(), l->getPlanet());
            double cellSize = 1.0f / cellNum;
            #pragma omp parallel for
                for (int i = 0; i < cellNum; ++i){
                    for (int j = 0; j < cellNum; ++j){
                        double gradX = scaleFactor * (dt / rho) * (l->getPressure(i + 1, j) - l->getPressure(i - 1, j)) / cellSize;
                        double gradY = scaleFactor * (dt / rho) * (l->getPressure(i, j + 1) - l->getPressure(i, j - 1)) / cellSize;
                        bufferVelTheta->setData(i, j, (l->getVelocityTheta(i, j) - gradX));
                        bufferVelPhi->setData(i, j, (l->getVelocityPhi(i, j) - gradY));
                    }
                }
            #pragma omp barrier

            //swap
            l->swapVels(bufferVelTheta, bufferVelPhi);
        }

        template<>
        inline void pressureEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>& l, double scaleFactor){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            double velSolid = 0;
            auto dt = getDt();
            auto rho = PWM::Utils::altitudeAdjustedDensity(l->getHeight(), l->getPlanet());
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
//                        if (i == 256 && j == 256){
//                            std::cerr << "Debug point reached" << std::endl;
//                        }
                        //Check if the central cell is solid
                        if (l->getObstacles(i, j)){
                            bufferVelTheta->setData(i, j, velSolid);
                            bufferVelPhi->setData(i, j, velSolid);
                            continue;
                        }
                        //If not, check if cell above is solid
                        if (l->getObstacles(i - 1, j)){
                            bufferVelTheta->setData(i, j, velSolid);
                        }
                        else{//Pressure update as per normal formula
                            double gradX = scaleFactor * (dt / rho) * (l->getPressure(i, j) - l->getPressure(i - 1, j)) / cellSizeX;
                            double newX = l->getVelocityTheta(i, j) - gradX;
                            bufferVelTheta->setData(i, j, newX);
                        }
                        //Also check if cell left is solid
                        if (l->getObstacles(i, j - 1)){
                            bufferVelPhi->setData(i, j, velSolid);
                        }
                        else{//Pressure update as per normal formula
                            double gradY = scaleFactor * (dt / rho) * (l->getPressure(i, j) - l->getPressure(i, j - 1)) / cellSizeY;
                            double newY = l->getVelocityPhi(i, j) - gradY;
                            bufferVelPhi->setData(i, j, newY);
                        }
                    }
                }
            #pragma omp barrier

            //swap
            l->swapVels(bufferVelTheta, bufferVelPhi);
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::solvePressureProjection(){
            startComputation();
            for (auto l : airLayers){
                solvePressureProjection(l);
            }
            endComputation();
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::calcPressures(){
            startComputation();
            for (auto l : airLayers){
                solvePressureMatrix(l);
                // std::cout << "X velocities after advection:\n" << l->getVelocityTheta().print() << std::endl;
            }
            endComputation();
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::step_internal(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::step_internal() for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<PWM::PWMDataStructure::square2DArray<double>, double>::step_internal(){
            for (auto l : airLayers){
                calcDiv(l);
                // std::cout << "Divergence:\n" << divergence->print() << std::endl;
                l->getPressure().mul(0);
                multiGridSolve(l->getPressure(), *divergence);
                //std::cout << "Pressure:\n" << l->getPressure().print() << std::endl;
                solvePressureProjection(l);
            }
        }

        template<>
        inline void pressureEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::step_internal(){
            for (auto l : airLayers){
                solvePressureMatrix(l);
                solvePressureProjection(l);
                // std::cout << "X velocities after advection:\n" << l->getVelocityTheta().print() << std::endl;
            }
        }
    }
}
#endif // PRESSUREENGINE_H
