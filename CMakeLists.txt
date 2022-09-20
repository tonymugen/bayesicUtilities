cmake_minimum_required(VERSION 3.11)
project(brUtilities VERSION 1.0 LANGUAGES CXX)

get_property(isMultiConfig GLOBAL
	PROPERTY GENERATOR_IS_MULTICONFIG
)
if(NOT isMultiConfig)
	set(allowedBuildTypes Debug Release)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
		STRINGS "${allowedBuildTypes}"
	)
	if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE) # default build type is Release
	elseif(NOT CMAKE_BUILD_TYPE IN_LIST allowedBuildTypes)
		message(FATAL_ERROR
			"Unknown build type ${CMAKE_BUILD_TYPE}"
		)
	endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

add_library(brUtilities
	src/random.cpp
	src/index.cpp
	src/utilities.cpp
)

set_target_properties(brUtilities PROPERTIES
	CXX_STANDARD         14
	CXX_STANDARD_REQURED ON
	CXX_EXTENSIONS       OFF
)

target_include_directories(brUtilities
	PRIVATE include
)
