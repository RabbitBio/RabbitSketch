cmake_minimum_required(VERSION 3.16)

project(rabbitsketch VERSION 2.0.0 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# The codebase uses C++17 facilities; spelling -std=c++20 directly broke the
# supported GCC 8 toolchain (which only recognizes the draft c++2a spelling).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(CXXAPI "build pure c++ library instead of python bind interface" OFF)
# FastKMV key pipeline (only one of these applies to fastkmv.cpp):
#   default     = single-round murmur3-style fmix after canonical ntHash
#   DOUBLE_FMUX = two-round fmix (-DFASTKMV_DOUBLE_FMUX)
#   NO_FMUX     = raw (canonical ntHash >> 11), ablation (-DFASTKMV_NO_FMUX)
option(FASTKMV_NO_FMUX "FastKMV: skip fmix (raw ntHash keys)" OFF)
option(FASTKMV_DOUBLE_FMUX "FastKMV: two-round fmix instead of one" OFF)
option(RABBITSKETCH_BUILD_TESTS "Build RabbitSketch tests" OFF)
option(RABBITSKETCH_NATIVE_ARCH
       "Opt in to -march=native instead of the portable runtime-dispatched build" OFF)
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
find_package(OpenMP QUIET)
#SET(CMAKE_CXX_COMPILER /usr/bin/g++)
#SET(CMAKE_C_COMPILER /usr/bin/gcc)

include_directories(src)

aux_source_directory(src DIR_SRCS)
aux_source_directory(src/api API_SRCS)
aux_source_directory(src/estimators ESTIMATOR_SRCS)
aux_source_directory(src/io IO_SRCS)
list(APPEND DIR_SRCS ${API_SRCS} ${ESTIMATOR_SRCS} ${IO_SRCS})
list(REMOVE_ITEM DIR_SRCS "src/treeminhash.cpp")

# No -flto by default: LTO object code in .a is tied to the exact GCC that built it;
# linking that archive with another toolchain (e.g. conda g++) fails at lto-wrapper.
# Optional LTO (same compiler for lib + app): cmake -DCMAKE_CXX_FLAGS="-flto ..." ...
#linux only 
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -std=c++14 -static-libstdc++ -march=native")
if(RABBITSKETCH_NATIVE_ARCH)
	if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
		add_compile_options(-march=native)
		add_compile_definitions(RABBITSKETCH_NATIVE_ARCH=1)
	else()
		message(WARNING "RABBITSKETCH_NATIVE_ARCH is not implemented for this compiler")
	endif()
endif()

#set flags for specific files
#set_source_files_properties(src/MinHash.cpp PROPERTIES COMPILE_FLAGS -mavx2)
if(FASTKMV_NO_FMUX AND FASTKMV_DOUBLE_FMUX)
	message(FATAL_ERROR "FASTKMV_NO_FMUX and FASTKMV_DOUBLE_FMUX are mutually exclusive")
endif()
# Always refresh so toggling options does not leave stale -D flags on fastkmv.cpp
if(FASTKMV_NO_FMUX)
	set_source_files_properties(src/fastkmv.cpp PROPERTIES COMPILE_DEFINITIONS FASTKMV_NO_FMUX)
elseif(FASTKMV_DOUBLE_FMUX)
	set_source_files_properties(src/fastkmv.cpp PROPERTIES COMPILE_DEFINITIONS FASTKMV_DOUBLE_FMUX)
else()
	# Clear per-file defs (no value = remove property; portable vs. unset_property)
	set_property(SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/fastkmv.cpp" PROPERTY COMPILE_DEFINITIONS)
endif()


if (CXXAPI)
	add_compile_options(-DNOPYTHON)
	add_library("${PROJECT_NAME}" SHARED ${DIR_SRCS})
	add_library("${PROJECT_NAME}_static" STATIC ${DIR_SRCS})
	add_library(RabbitSketch::rabbitsketch ALIAS "${PROJECT_NAME}")
	add_library(RabbitSketch::rabbitsketch_static ALIAS "${PROJECT_NAME}_static")
	target_include_directories("${PROJECT_NAME}" PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
	target_include_directories("${PROJECT_NAME}_static" PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
	target_link_libraries("${PROJECT_NAME}" PRIVATE ZLIB::ZLIB)
	target_link_libraries("${PROJECT_NAME}_static" PUBLIC ZLIB::ZLIB)
	if(OpenMP_CXX_FOUND)
		target_link_libraries("${PROJECT_NAME}" PRIVATE OpenMP::OpenMP_CXX)
		target_link_libraries("${PROJECT_NAME}_static" PUBLIC OpenMP::OpenMP_CXX)
	endif()
	file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
	add_custom_command(TARGET "${PROJECT_NAME}_static" POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			"$<TARGET_FILE:${PROJECT_NAME}_static>"
			"${CMAKE_BINARY_DIR}/lib/")
	
  
	#set_target_properties("${PROJECT_NAME}_dynamic" PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
	#set_target_properties("${PROJECT_NAME}_static" PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

	install(TARGETS "${PROJECT_NAME}" "${PROJECT_NAME}_static"
		   EXPORT RabbitSketchTargets
		   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
		   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
		   INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)
	# Public headers include several header-only support libraries (phmap,
	# robin_hood, xxhash, ...). Install the complete header closure.
	install(DIRECTORY "src/" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
		FILES_MATCHING
		PATTERN "*.h"
		PATTERN "pybind.h" EXCLUDE)
	# Preserve the legacy source-header installation contract for consumers
	# that historically compiled these implementation units directly.
	install(FILES "src/common.cpp" "src/Kssd.cpp" "src/shuffle.cpp"
		DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

	configure_package_config_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake/RabbitSketchConfig.cmake.in"
		"${CMAKE_CURRENT_BINARY_DIR}/RabbitSketchConfig.cmake"
		INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/RabbitSketch")
	write_basic_package_version_file(
		"${CMAKE_CURRENT_BINARY_DIR}/RabbitSketchConfigVersion.cmake"
		VERSION ${PROJECT_VERSION}
		COMPATIBILITY SameMajorVersion)
	install(EXPORT RabbitSketchTargets
		FILE RabbitSketchTargets.cmake
		NAMESPACE RabbitSketch::
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/RabbitSketch")
	install(FILES
		"${CMAKE_CURRENT_BINARY_DIR}/RabbitSketchConfig.cmake"
		"${CMAKE_CURRENT_BINARY_DIR}/RabbitSketchConfigVersion.cmake"
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/RabbitSketch")
	configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/rabbitsketch.pc.in"
		"${CMAKE_CURRENT_BINARY_DIR}/rabbitsketch.pc" @ONLY)
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/rabbitsketch.pc"
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

	if(RABBITSKETCH_BUILD_TESTS)
		enable_testing()
		add_executable(unified_query_test tests/unified_query_test.cpp)
		target_link_libraries(unified_query_test rabbitsketch_static)
		add_test(NAME unified_query COMMAND unified_query_test)
		add_executable(fastx_builder_test tests/fastx_builder_test.cpp)
		target_link_libraries(fastx_builder_test
			rabbitsketch_static ${ZLIB_LIBRARIES})
		add_test(NAME fastx_builder COMMAND fastx_builder_test)
		add_executable(orderminhash_query_test
			tests/orderminhash_query_test.cpp)
		target_link_libraries(orderminhash_query_test rabbitsketch_static)
		add_test(NAME orderminhash_query COMMAND orderminhash_query_test)
		add_executable(runtime_dispatch_test
			tests/runtime_dispatch_test.cpp)
		target_link_libraries(runtime_dispatch_test rabbitsketch_static)
		add_test(NAME runtime_dispatch COMMAND runtime_dispatch_test)
		add_test(NAME runtime_dispatch_scalar COMMAND runtime_dispatch_test)
		set_tests_properties(runtime_dispatch_scalar PROPERTIES
			ENVIRONMENT "RABBITSKETCH_SIMD=scalar")
		find_package(Python3 COMPONENTS Interpreter QUIET)
		if(Python3_Interpreter_FOUND)
			add_test(NAME rankstream_python_vectors
				COMMAND ${Python3_EXECUTABLE}
					"${CMAKE_CURRENT_SOURCE_DIR}/tests/test_rankstream_vectors.py")
		endif()
	endif()

	#install(FILES "src/hll/hyperloglog.h"
	#	   DESTINATION include/hll
	#	   )
else ()
  if(RABBITSKETCH_BUILD_TESTS)
    message(FATAL_ERROR "RABBITSKETCH_BUILD_TESTS requires -DCXXAPI=ON")
  endif()
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(rabbitsketch ${DIR_SRCS})
  target_link_libraries(rabbitsketch PRIVATE ZLIB::ZLIB)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(rabbitsketch PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif ()
