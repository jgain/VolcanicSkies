/**
 * @file
 *
 * Startup code that must be called by @c main.
 */

/*
#include <thread>
#ifdef _WIN32
#include <Eigen/Core>
#elif
#include <eigen3/Eigen/Core>
#endif
#include <ImfThreading.h>
*/
#include "initialize.h"

void utsInitialize()
{
    // Eigen::initParallel();
    // Imf::setGlobalThreadCount(std::thread::hardware_concurrency());
}
