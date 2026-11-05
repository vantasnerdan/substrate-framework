# P254 A40: unmodified pinned libcntr plus the actual real-time IVP consumer.
cmake_minimum_required(VERSION 3.16)
project(P254_NESSi_Gaussian_IVP LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(NESSI_SOURCE "/home/administrator/.cache/dft/p254-nessi-v2-source" CACHE PATH "Pinned NESSi 6bb4c832 source")
set(omp OFF CACHE BOOL "No OpenMP" FORCE)
set(mpi OFF CACHE BOOL "No MPI" FORCE)
set(hdf5 OFF CACHE BOOL "No HDF5" FORCE)
set(ness OFF CACHE BOOL "No steady-state/FFTW mode" FORCE)
set(BUILD_TEST OFF CACHE BOOL "No generic upstream suite" FORCE)
set(BUILD_DOC OFF CACHE BOOL "No upstream documentation generation" FORCE)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
# Exclusion suppresses upstream install rules (one has a nonexistent cntr/utils
# relative directory in this pin). We install the same headers/library correctly
# below, without modifying or copying the third-party source tree.
add_subdirectory("${NESSI_SOURCE}/libcntr" nessi-libcntr EXCLUDE_FROM_ALL)
add_executable(nessi_gaussian_embedding nessi_gaussian_embedding.cpp)
target_include_directories(nessi_gaussian_embedding PRIVATE "${NESSI_SOURCE}/libcntr")
target_link_libraries(nessi_gaussian_embedding PRIVATE cntr Eigen3::Eigen)
set_target_properties(nessi_gaussian_embedding PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
install(TARGETS cntr LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/cntr)
install(DIRECTORY "${NESSI_SOURCE}/libcntr/cntr/utils/" DESTINATION include/cntr/utils
        FILES_MATCHING PATTERN "*.hpp")
install(TARGETS nessi_gaussian_embedding RUNTIME DESTINATION bin)
install(FILES nessi_gaussian_ivp.hpp DESTINATION include/p254)
install(FILES "${NESSI_SOURCE}/LICENSE.md" DESTINATION share/licenses/nessi RENAME LICENSE.md)
