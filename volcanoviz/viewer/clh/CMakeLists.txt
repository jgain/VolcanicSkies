if (BUILD_CLH)

    add_library(clh clh.cpp timer.cpp)
    target_link_libraries(clh common ${Boost_PROGRAM_OPTIONS_LIBRARY} ${OPENCL_LIBRARIES})

    add_executable(climgformats climgformats.cpp)
    target_link_libraries(climgformats ${OPENCL_LIBRARIES})

endif()
