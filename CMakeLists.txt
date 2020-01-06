project(h-tautau)
cmake_minimum_required(VERSION 3.6)

get_filename_component(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools" ABSOLUTE)
include("${AnalysisTools_DIR}/cmake/include/common.cmake")

file(GLOB_RECURSE ANA_TOOLS_SRC "${AnalysisTools_DIR}/*.cpp")
add_library(AnalysisTools STATIC ${ANA_TOOLS_SRC} "${RootDict}")
add_dependencies(AnalysisTools GenerateRootDict)

file(GLOB_RECURSE HTT_UTILITIES_SRC "${CMSSW_BASE_SRC}/HTT-utilities/*.cc")
add_library(HTT-utilities STATIC ${HTT_UTILITIES_SRC})
set_source_files_properties(${HTT_UTILITIES_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE HH_KINFIT2_SRC "${CMSSW_BASE_SRC}/HHKinFit2/*.cpp")
add_library(HHKinFit2 STATIC ${HH_KINFIT2_SRC})
target_include_directories(HHKinFit2 PRIVATE "${CMSSW_BASE_SRC}/HHKinFit2/HHKinFit2/interface")
set_source_files_properties(${HH_KINFIT2_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE SVFIT_SRC "${CMSSW_BASE_SRC}/TauAnalysis/ClassicSVfit/*.cc")
list(FILTER SVFIT_SRC EXCLUDE REGEX ".*test.*")
add_library(SVfit STATIC ${SVFIT_SRC})
set_source_files_properties(${SVFIT_SRC} PROPERTIES COMPILE_FLAGS "-w")

execute_process(COMMAND scram tool info Eigen OUTPUT_VARIABLE SCRAM_E_INFO)
string(REGEX MATCH "\nINCLUDE=([^\n]*)" E_STR ${SCRAM_E_INFO})
set(Eigen_INCLUDE_DIRS "${CMAKE_MATCH_1}")

execute_process(COMMAND scram tool info protobuf OUTPUT_VARIABLE SCRAM_PB_INFO)
string(REGEX MATCH "\nINCLUDE=([^\n]*)" PB_STR ${SCRAM_PB_INFO})
set(Protobuf_INCLUDE_DIRS "${CMAKE_MATCH_1}")
string(REGEX MATCH "\nLIBDIR=([^\n]*)" PB_STR ${SCRAM_PB_INFO})
set(Protobuf_LIBRARY_DIRS "${CMAKE_MATCH_1}")
set(Protobuf_LIBRARIES "${CMAKE_MATCH_1}/libprotobuf.so")

execute_process(COMMAND scram tool info tensorflow-c OUTPUT_VARIABLE SCRAM_TFC_INFO)
string(REGEX MATCH "\nINCLUDE=([^\n]*)" TFC_STR ${SCRAM_TFC_INFO})
set(APPEND Tensorflow-c_INCLUDE_DIRS "${CMAKE_MATCH_1}")
string(REGEX MATCH "\nLIBDIR=([^\n]*)" TFC_STR ${SCRAM_TFC_INFO})
set(Tensorflow-c_LIBRARY_DIRS "${CMAKE_MATCH_1}")
set(Tensorflow-c_LIBRARIES "${CMAKE_MATCH_1}/libtensorflow.so")

execute_process(COMMAND scram tool info tensorflow-cc OUTPUT_VARIABLE SCRAM_TFCC_INFO)
string(REGEX MATCH "\nINCLUDE=([^\n]*)" TFCC_STR ${SCRAM_TFCC_INFO})
set(Tensorflow-cc_INCLUDE_DIRS "${CMAKE_MATCH_1}")
string(REGEX MATCH "\nLIBDIR=([^\n]*)" TFCC_STR ${SCRAM_TFCC_INFO})
set(Tensorflow-cc_LIBRARY_DIRS "${CMAKE_MATCH_1}")
set(Tensorflow-cc_LIBRARIES "${Tensorflow-cc_LIBRARY_DIRS}/libtensorflow_cc.so")

include_directories(SYSTEM ${Tensorflow-cc_INCLUDE_DIRS} ${Tensorflow-c_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS} ${Protobuf_INCLUDE_DIRS})

target_link_libraries("nutella" ${Tensorflow-cc_LIBRARIES} ${Tensorflow-c_LIBRARIES} ${Protobuf_LIBRARIES} pthread $ENV{CMSSW_RELEASE_BASE}/lib/slc6_amd64_gcc630/libPhysicsToolsTensorFlow.so)

foreach(exe_name ${EXE_LIST})
    target_link_libraries("${exe_name}" h-tautau AnalysisTools HTT-utilities SVfit HHKinFit2)
endforeach()

add_library(h-tautau STATIC ${SOURCE_LIST})

set_target_properties(Print_Graph Print_Graph_2 Print_Graph_3 Print_TreeBranch2D PROPERTIES EXCLUDE_FROM_ALL 1)
