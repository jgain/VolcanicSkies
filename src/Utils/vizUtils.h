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
#ifndef VIZUTILS_H
#define VIZUTILS_H

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include "flatStaggeredGrid.h"
#include <opencv2/imgcodecs.hpp>
#include "planet.h"
#include "skirt.hpp"
#include "square2DArray.h"
#include "terrain_structure.hpp"
#include <tuple>
#include <vector>
#include "world.h"

namespace PWM {
    namespace Utils {
        inline std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> cloudColourMap = {
            std::make_tuple(255, 255, 255), //EMPTY
            std::make_tuple(0, 255, 0),//STRATUS,
            std::make_tuple(0, 153, 51),//ALTOSTRATUS,
            std::make_tuple(0, 204, 255),//CIRROSTRATUS,
            std::make_tuple(153, 255, 51),//NIMBOSTRATUS,
            std::make_tuple(0, 0, 255),//CIRRUS,
            std::make_tuple(255, 0, 0),//CUMULUS,
            std::make_tuple(255, 255, 0),//STRATOCUMULUS,
            std::make_tuple(204, 0, 102),//ALTOCUMULUS,
            std::make_tuple(255, 0, 255),//CIRROCUMULUS,
            std::make_tuple(255, 153, 0),//CUMULONIMBUS,
            std::make_tuple(0, 0, 0)//CTYPEEND
        };

        template<typename T>
        inline int writeMoisImage(std::string file, std::shared_ptr<PWM::Model::planet> planetos, const PWM::PWMDataStructure::square2DArray<T>& data, std::pair<T, T>& minMax){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = minMax;
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.size(); ++i){
                    float moisNormed = (std::abs(spread) < 1 ? 0 : ((data.getData(i) - mM.first) / spread) * 255);
                    int r = 255 - (255 - moisNormed);
                    int g = 0;
                    int b = 255 - moisNormed;
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeMoisImage(std::string file, std::shared_ptr<PWM::Model::planet> planetos, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(0.2, 1.0);
                T max = 2.0;
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;
                    float moisNormed, mois = data.getData(i);
                    if (mois <= 0){
                        r = g = b = 0;
                    }
                    else if (mois < mM.first){
                        moisNormed = ((mois) / mM.first) * 255;
                        r = 0;
                        g = 0;
                        b = 255 - (255 - moisNormed);
                    }
                    else if (mois < mM.second){
                        moisNormed = ((mois - mM.first) / spread) * 255;
                        r = 0;
                        g = 255 - (255 - moisNormed);
                        b = 255 - moisNormed;
                    }
                    else if (mois < max){
                        moisNormed = (((mois - mM.second) / (max - mM.second)) * 255);
                        r = 255 - (255 - moisNormed);
                        g = 255 - moisNormed;
                        b = 0;
                    }
                    else{
                        r = g = b = 255;
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeMoisImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            //Image data
            auto mM = std::make_pair<T, T>(0.1, 1.0);
            T max = 2.0;
            T spread = mM.second - mM.first;
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    int r, g, b;
                    float moisNormed, mois = data.getData(i, j);
                    if (mois <= 0){
                        r = g = b = 0;
                    }
                    else if (mois < mM.first){
                        moisNormed = ((mois) / mM.first) * 255;
                        r = 0;
                        g = 0;
                        b = 255 - (255 - moisNormed);
                    }
                    else if (mois < mM.second){
                        moisNormed = ((mois - mM.first) / spread) * 255;
                        r = 0;
                        g = 255 - (255 - moisNormed);
                        b = 255 - moisNormed;
                    }
                    else if (mois < max){
                        moisNormed = (((mois - mM.second) / (max - mM.second)) * 255);
                        r = 255 - (255 - moisNormed);
                        g = 255 - moisNormed;
                        b = 0;
                    }
                    else{
                        r = g = b = 255;
                    }

                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T1, typename T2>
        inline int writeMoisImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T1>& data, const PWM::PWMDataStructure::flatStaggeredGrid<T2>& obs){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            //Image data
            auto mM = std::make_pair<T1, T1>(0.1, 1.0);
            T1 max = 2.0;
            T1 spread = mM.second - mM.first;
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    int r, g, b;
                    if (obs.getData(i, j)){
                        r = g = b = 0;
                    }
                    else{
                        float moisNormed, mois = data.getData(i, j);
                        if (mois <= 0){
                            r = g = b = 0;
                        }
                        else if (mois < mM.first){
                            moisNormed = ((mois) / mM.first) * 255;
                            r = 0;
                            g = 0;
                            b = 255 - (255 - moisNormed);
                        }
                        else if (mois < mM.second){
                            moisNormed = ((mois - mM.first) / spread) * 255;
                            r = 0;
                            g = 255 - (255 - moisNormed);
                            b = 255 - moisNormed;
                        }
                        else if (mois < max){
                            moisNormed = (((mois - mM.second) / (max - mM.second)) * 255);
                            r = 255 - (255 - moisNormed);
                            g = 255 - moisNormed;
                            b = 0;
                        }
                        else{
                            r = g = b = 255;
                        }
                    }
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeTempImage(std::string file, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(240, 350);
                T spread = mM.second - mM.first;
                T max = 1200;

                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;

                    T t = data.getData(i);

                    // JG - change of viz to help build intuition
                    // colour band approach
                    if(t < 213.0f) // < -60C
                    {
                        r = 0; g = 0; b = 0; // black
                    }
                    else if(t < 263.15f) // [-60C, -10C]
                    {
                        r = 0; g = 0; b = 255; // blue
                    }
                    else if(t < 283.15f) // [-10C,10C]
                    {
                        r = 255; g = 0; b = 255; // purple
                    }
                    else if(t < 373.15)
                    {
                        r = 255; g = (t-283.15) / 90.0f; b = (t-283.15) / 90.0f; // red shading to white
                    }
                    else
                    {
                        r = 255; g = 255; b = 255; // white
                    }

                    /*
                            float tNormed, t = data.getData(i);
                            if (t < mM.first){
                                tNormed = (std::abs(spread) < 1 ? 0 : ((t) / mM.first) * 255);
                                r = 0;
                                g = 0;
                                b = 255 - (255 - tNormed);
                            }
                            else if (t < mM.second){
                                tNormed = (std::abs(spread) < 1 ? 0 : ((t - mM.first) / spread) * 255);
                                r = 255 - (255 - tNormed);
                                g = 0;
                                b = 255 - tNormed;
                            }
                            else if (t < max){
                                tNormed = (((t - mM.second) / (max - mM.second)) * 255);
                                r = 255;
                                g = 255 - (255 - tNormed);
                                b = 255 - (255 - tNormed);
                            }
                            else{
                                r = g = b = 255;
                            } */
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeTempImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            //Image data
            auto mM = std::make_pair<T, T>(240, 350);
            T spread = mM.second - mM.first;
            T max = 1200;
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    int r, g, b;
                    int tNormed;
                    T t = data.getData(i, j);
                    if (t < mM.first){
                        tNormed = ((t) / mM.first) * 255;
                        r = 0;
                        g = 0;
                        b = 255 - (255 - tNormed);
                    }
                    else if (t < mM.second){
                        tNormed = ((t - mM.first) / spread) * 255;
                        r = 255 - (255 - tNormed);
                        g = 0;
                        b = 255 - tNormed;
                    }
                    else if (t < max){
                        tNormed = ((t - mM.second) / (max - mM.second)) * 255;
                        r = 255;
                        g = 255 - (255 - tNormed);
                        b = 255 - (255 - tNormed);
                    }
                    else{
                        r = g = b = 255;
                    }

                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T1, typename T2>
        inline int writeTempImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T1>& data, const PWM::PWMDataStructure::flatStaggeredGrid<T2>& obs){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            //Image data
            auto mM = std::make_pair<T1, T1>(240, 350);
            T1 spread = mM.second - mM.first;
            T1 max = 1200;
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    int r, g, b;
                    if (obs.getData(i, j)){
                        r = 100;
                        g = 182;
                        b = 106;
                    }
                    else{
                        int tNormed;
                        T1 t = data.getData(i, j);
                        if (t < mM.first){
                            tNormed = ((t) / mM.first) * 255;
                            r = 0;
                            g = 0;
                            b = 255 - (255 - tNormed);
                        }
                        else if (t < mM.second){
                            tNormed = ((t - mM.first) / spread) * 255;
                            r = 255 - (255 - tNormed);
                            g = 0;
                            b = 255 - tNormed;
                        }
                        else if (t < max){
                            tNormed = ((t - mM.second) / (max - mM.second)) * 255;
                            r = 255;
                            g = 255 - (255 - tNormed);
                            b = 255 - (255 - tNormed);
                        }
                        else{
                            r = g = b = 255;
                        }
                    }
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeAshImage(std::string file, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(240, 350);
                T spread = mM.second - mM.first;
                T max = 1;
                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;
                    T rho = data.getData(i);
                    if (rho > max)
                    {
                        r,g=0;
                        b = 255;
                    }
                    else
                    {
                        r=255 - 255*rho;
                        g=255 - 255*rho;
                        b=255;
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeAshImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            T max = 10;
            //Image data
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    int r, g, b;
                    T rho = data.getData(i, j) * 3;
                    if (rho > max)
                    {
                        g = b = 0;
                        r = 255;
                    }
                    else
                    {
                        r = 255;
                        g = 255 - 255 * rho;
                        b = 255 - 255 * rho;
                    }

                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }

            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeTerrAshImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            T max = 10;
            //Image data
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    int r, g, b;
                    T rho = data.getData(i, j) * 3;
                    if (rho > max)
                    {
                        r = g = 0;
                        b = 255;
                    }
                    else
                    {
                        r = 255 - 255 * rho;
                        g = 255 - 255 * rho;
                        b = 255;
                    }

                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }

            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeTerrAshFile(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data, std::pair<double, double> plumeLoc){
            size_t nX = data.getX(), nY = data.getY();
            std::stringstream out;
            std::ofstream fil;
            try {
                fil.open(file);
                fil << "# Elev map of ash thickness" << std::endl;
                fil << nX << "\t" << nY << std::endl;
                fil << data.cellSizeX() << "\t" << data.cellSizeY() << std::endl;
                fil << plumeLoc.first << "\t" << plumeLoc.second << std::endl;
                fil << std::setprecision(5);
                for (int i = 0; i < nX; ++i){
                    for (int j = 0; j < nY; ++j){
                        fil << data.getData(i, j) << "\t";
                    }
                    fil << std::endl;
                }
                fil.close();
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write terrain ash values to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
            return 0;
        }

        template<typename T1, typename T2>
        inline int writeAshImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T1>& data, const PWM::PWMDataStructure::flatStaggeredGrid<T2>& obs){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
                //Image data
            auto mM = std::make_pair<T1, T1>(240, 350);
            T1 spread = mM.second - mM.first;
            T1 max = 1;
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    int r, g, b;
                    if (obs.getData(i, j)){
                        r = g = b = 0;
                    }
                    T1 rho = data.getData(i, j);
                    if (rho > max)
                    {
                        r = g = 0;
                        b = 255;
                    }
                    else
                    {
                        r = 255 - 255 * rho;
                        g = 255 - 255 * rho;
                        b = 255;
                    }

                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }

            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
            return -1;

        }

        template<typename T>
        inline int writeVelImage(std::string file, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                T max = 40;
                int r, g, b;
                for (int i = 0; i < data.size(); ++i){
                    T t = data.getData(i);
                    int tNormed = ((std::abs(t) / max) * 255);
                    r = 0;
                    if (PWM::Utils::getSign(t) > 0){
                        if (t < max){
                            g = 255 - std::min(255 - tNormed, 255);
                            b = 0;
                        }
                        else{
                            g = 255;
                            b = 0;
                        }
                    }
                    else{
                        if (t > -max){
                            g = 0;
                            b = 255 - std::min(255 - tNormed, 255);

                        }
                        else{
                            g = 0;
                            b = 255;
                        }
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeVelImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            //Image data
            T max = 140;
            int r, g, b;
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    T t = data.getData(i, j);
                    int tNormed = ((std::abs(t) / max) * 255);
                    r = 0;
                    if (PWM::Utils::getSign(t) > 0){
                        g = 255 - std::min(255 - tNormed, 255);
                        b = 0;
                    }
                    else{
                        g = 0;
                        b = 255 - std::min(255 - tNormed, 255);
                    }

                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeVelImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data, const PWM::PWMDataStructure::flatStaggeredGrid<T>& obs){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            //Image data
            T max = 140;
            int r, g, b;
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    if (obs.getData(i, j)){
                        r = 83;
                        g = b = 19;
                    }
                    else{
                        T t = data.getData(i, j);
                        int tNormed = ((std::abs(t) / max) * 255);
                        r = 0;
                        if (PWM::Utils::getSign(t) > 0){
                            g = 255 - std::min(255 - tNormed, 255);
                            b = 0;
                        }
                        else{
                            g = 0;
                            b = 255 - std::min(255 - tNormed, 255);
                        }
                    }
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeCloudImages(std::string outputDirectory, int step, int lay, float height, const PWM::PWMDataStructure::square2DArray<T>& mask){
            std::stringstream fC;
            fC << outputDirectory << "/Layer_" << lay << "_(" << height << "m)_Clouds_Step_" << step << ".ppm";
            std::ofstream fil;
            try {
                //                std::cerr << "Writing layer " << j << " to " << fC.str() << "." << std::endl;
                fil.open(fC.str());
                //Header
                fil << "P3\n" << mask.getWidth() << " " << mask.getWidth() << "\n255" << std::endl;

                //Image data
                for (int x = 0; x < mask.getWidth(); ++x){
                    for (int y = 0; y < mask.getWidth(); ++y){
                        auto cloudType = static_cast<int>(mask.getData(x, y));
                        auto colour = cloudColourMap[cloudType];
                        fil << std::get<0>(colour) << " " << std::get<1>(colour) << " " << std::get<2>(colour) << " ";
                    }
                }
                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << fC.str() << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeCloudImages(std::string outputDirectory, int step, int lay, float height, const PWM::PWMDataStructure::flatStaggeredGrid<T>& mask){
            std::stringstream fC;
            fC << outputDirectory << "/Layer_" << lay << "_(" << height << "m)_Clouds_Step_" << step << ".png";
            size_t nX = mask.getX(), nY = mask.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    auto cloudType = static_cast<int>(mask.getData(i, j));
                    auto colour = cloudColourMap[cloudType];
                    int r = std::get<0>(colour);
                    int g = std::get<1>(colour);
                    int b = std::get<2>(colour);
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(fC.str(), img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << fC.str() << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeCloudImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& mask){
            size_t nX = mask.getX(), nY = mask.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    auto cloudType = static_cast<int>(mask.getData(i, j));
                    auto colour = cloudColourMap[cloudType];
                    int r = std::get<0>(colour);
                    int g = std::get<1>(colour);
                    int b = std::get<2>(colour);
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeTerrElevImage(std::string file, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = data.minmax();
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.size(); ++i){
                    float tNormed = ((data.getData(i) - mM.first) / spread) * 255;
                    int r = 255 - (255 - tNormed);
                    int g = 255 - (255 - tNormed);
                    int b = 255 - (255 - tNormed);
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeTerrElevImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << nX << " " << nY << "\n255" << std::endl;

                //Image data
                auto mM = data.minmax();
                T spread = mM.second - mM.first;
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        float tNormed = ((data.getData(i, j) - mM.first) / spread) * 255;
                        int r = 255 - (255 - tNormed);
                        int g = 255 - (255 - tNormed);
                        int b = 255 - (255 - tNormed);
                        fil << r << " " << g << " " << b << " ";
                    }
                    fil << '\n';
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeTerrElevImage(std::string file, const terrain_structure& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.field_size << " " << data.field_size << "\n255" << std::endl;

                //Image data
                T min = std::numeric_limits<T>::max(), max = std::numeric_limits<T>::min();
                for (int i = 0; i < data.field_size; ++i)
                    for (int j = 0; j < data.field_size; ++j){
                        double val = data.height_field[i][j];
                        if (val < min)
                            min = val;
                        if (val > max)
                            max = val;
                    }
                auto mM = std::make_pair(min, max);
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.field_size; ++i)
                    for (int j = 0; j < data.field_size; ++j){
                        float tNormed = (std::abs(spread) < 1 ? 0 : ((data.height_field[i][j] - mM.first) / spread) * 255);
                        int r = 255 - (255 - tNormed);
                        int g = 255 - (255 - tNormed);
                        int b = 255 - (255 - tNormed);
                        fil << r << " " << g << " " << b << " ";
                    }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T, typename T2>
        inline int writePresImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data, const PWM::PWMDataStructure::flatStaggeredGrid<T2>& obs){
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX, CV_8UC4);
            //Image data
            auto mM = data.minmax();
            T spread = mM.second - mM.first;
            for (int i = 0; i < nY; ++i){
                for (int j = 0; j < nX; ++j){
                    int r, g, b;
                    if (obs.getData(i, j)){
                        r = g = 206;
                        b = 25;
                    }
                    else{
                        float tNormed = ((data.getData(i, j) - mM.first) / spread) * 255;
                        r = 255 - (255 - tNormed);
                        g = 255 - (255 - tNormed);
                        b = 255 - (255 - tNormed);
                    }
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int temperatureStatistics(const PWM::PWMDataStructure::flatStaggeredGrid<T>& data)
        {
            float meant = 0.0f; float maxt = std::numeric_limits<float>::min();

            for (int i = 0; i < data.size(); ++i)
            {
                float t = (float) data.getData(i);
                meant += t;
                if(t > maxt)
                    maxt = t;
            }
            meant /= (float) data.size();
            std::cerr << "Temperature: mean = " << meant << " max = " << maxt << std::endl;
            return 0;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline int writeAirImages(std::string outputDirectory, int i, const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld){
            #pragma omp parallel for
            for (int j = 0; j < wld->getAirLayers().size(); ++j){
                int hght = (int) wld->getAirLayer(j)->getHeight();
                std::stringstream fT, fM, fCW, fV1, fV2, fCV, fCV1, fCV2, fP, fA;
                fT << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Temperature_Step_" << i << ".png";
                fM << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Moisture_Step_" << i << ".png";
                fCW << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_CondensedWater_Step_" << i << ".png";
                fV1 << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_VelocityY_Step_" << i << ".png";
                fV2 << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_VelocityX_Step_" << i << ".png";
                //                fP << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Pressure_Step_" << i << ".png";
                fA << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Particulates_Step_" << i << ".png";
                fCV << outputDirectory << "/ConvectionLayer_" << j << "_Vertical_Velocity_Step_" << i << ".png";
                fCV1 << outputDirectory << "/ConvectionLayer_" << j << "_Rainfall_Step_" << i << ".png";
                fCV2 << outputDirectory << "/ConvectionLayer_" << j << "_Ashfall_Step_" << i << ".png";
//                writeTempImage<V>(fT.str(), wld->getAirLayer(j)->getTemperature());
                if (j < wld->getAirLayers().size() - 1){
//                    writeVelImage(fCV.str(), wld->getConvectionLayer(j)->getVerticalVelocities());
//                    writeMoisImage(fCV1.str(), wld->getConvectionLayer(j)->getRainfall());
//                    writeAshImage(fCV2.str(), wld->getConvectionLayer(j)->getAshfall());
                }

//                writeMoisImage<V>(fCW.str(), wld->getAirLayer(j)->getCondensedWater());
//                writeMoisImage<V>(fM.str(), wld->getAirLayer(j)->getMoisture());
//                writeVelImage<V>(fV1.str(), wld->getAirLayer(j)->getVelocityTheta());
//                writeVelImage<V>(fV2.str(), wld->getAirLayer(j)->getVelocityPhi());
//                writeAshImage(fA.str(), wld->getAirLayer(j)->getParticulates());
                writeCloudImages<V>(outputDirectory, i, j, wld->getAirLayer(j)->getHeight(), wld->getAirLayer(j)->getClouds());
            }
            #pragma omp barrier
            return 0;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline int writeAirImagesSelect(std::string outputDirectory, int i, const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld){
            #pragma omp parallel for
            for (int j = 0; j < wld->getAirLayers().size(); ++j){ // ++j
                // JG change
                int hght = (int) wld->getAirLayer(j)->getHeight();
                std::stringstream fT, fM, fCW, fV1, fV2, fCV, fCV1, fCV2, fP, fA;
                fT << outputDirectory << "/Layer_" << j << "_(" << hght  << "m)_Temperature_Step_" << i << ".png";
                fM << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_Moisture_Step_" << i << ".png";
                // fCW << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_CondensedWater_Step_" << i << ".png";
                fV1 << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_VelocityX_Step_" << i << ".png";
                fV2 << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_VelocityY_Step_" << i << ".png";
//                fP << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Pressure_Step_" << i << ".png";
                fA << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_Particulates_Step_" << i << ".png";
                fCV << outputDirectory << "/ConvectionLayer_" << j << "_Vertical_Velocity_Step_" << i << ".png";
                fCV1 << outputDirectory << "/ConvectionLayer_" << j << "_Rainfall_Step_" << i << ".png";
                fCV2 << outputDirectory << "/ConvectionLayer_" << j << "_Ashfall_Step_" << i << ".png";
//                writeTempImage<V>(fT.str(), wld->getAirLayer(j)->getTemperature());
                if (j < wld->getAirLayers().size() - 1){
//                    writeVelImage(fCV.str(), wld->getConvectionLayer(j)->getVerticalVelocities());
//                    writeMoisImage(fCV1.str(), wld->getConvectionLayer(j)->getRainfall());
//                    writeAshImage(fCV2.str(), wld->getConvectionLayer(j)->getAshfall());
                }
               // writeMoisImage<V>(fCW.str(), wld->getAirLayer(j)->getCondensedWater());
               // writeMoisImage<V>(fM.str(), wld->getAirLayer(j)->getMoisture());
                // writeVelImage<V>(fV1.str(), wld->getAirLayer(j)->getVelocityTheta());
                // writeVelImage<V>(fV2.str(), wld->getAirLayer(j)->getVelocityPhi());
               writeAshImage(fA.str(), wld->getAirLayer(j)->getParticulates());
                // writeCloudImages<V>(outputDirectory, i, j, wld->getAirLayer(j)->getHeight(), wld->getAirLayer(j)->getClouds());
                // JG change to handle memory
            }
            #pragma omp barrier
            return 0;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline int writeSkirts(std::string file, int step, const skirt& skir, const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld){
            size_t nX = wld->getAirLayer(0)->getObstacles().getX(), nY = wld->getAirLayer(0)->getObstacles().getX();
            return 1;
        }
    
        template<typename T, typename TT, typename V, typename VV>
        inline void airImagesAnalysis(int i, int cx, int cy, float r, const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld){
            // report proportion of clouds within a specified radius of the vent across all layers, and ashfall over the entire domain
            
            int cpixelcount = 0, cloudcount = 0, ashfallcount = 0, apixelcount = 0;
            for (int j = 0; j < wld->getAirLayers().size(); j++){
                auto & cmask = wld->getAirLayer(j)->getClouds();
                size_t nx = cmask.getX(), ny = cmask.getY();
               
                // look for clouds within a pixel radius of the vent
                for (int x = 0; x < nx; ++x){
                    for (int y = 0; y < ny; ++y){
                        float d = (cx - x) * (cx - x) + (cy - y) * (cy - y);
                        if(d < r*r)
                        {
                            auto cloudType = static_cast<int>(cmask.getData(x, y));
                            if(cloudType > 0)
                                cloudcount++;
                            cpixelcount++;
                        }
                    }
                }
                
                if(j < wld->getAirLayers().size()-1)
                {
                    // PWM::PWMDataStructure::square2DArray<T>&
                    auto & amask = wld->getConvectionLayer(j)->getAshfall();
                    nx = amask.getX(), ny = amask.getY();
                    
                    // search for ashfall over the entire domain
                    for (int x = 0; x < nx; ++x){
                        for (int y = 0; y < ny; ++y){
                            auto ashval = static_cast<int>(amask.getData(x, y));
                            if(ashval > 0.0000001f)
                                ashfallcount++;
                            apixelcount++;
                        }
                    }
                }
             
            }
            std::cerr << "TIMESTEP " << i << std::endl;
            std::cerr << "Cloud proportion = " << (float) cloudcount / (float) cpixelcount << std::endl;
            std::cerr << "Ashfall proportion = " << (float) ashfallcount / (float) apixelcount << std::endl;
            if(ashfallcount > 0 || cloudcount > 0)
                std::cerr << "*********************" << std::endl;
        }
    }
}
#endif // VIZUTILS_H
