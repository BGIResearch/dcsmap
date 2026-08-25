cmake_minimum_required(VERSION 3.10)
project(dcsmap LANGUAGES CXX VERSION 1.0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Pure static linking: statically link libstdc++, libgcc, libpthread and libc
# into a single self-contained binary with no runtime shared-library deps.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

find_package(Threads REQUIRED)

add_executable(dcsmap
    src/main.cpp
    src/Task.cpp
    src/Scripts.cpp
    src/simp_logger.cpp
)

target_include_directories(dcsmap PRIVATE src)
target_compile_definitions(dcsmap PRIVATE DCSMAP_VERSION="${PROJECT_VERSION}")
target_link_libraries(dcsmap PRIVATE Threads::Threads)

if(NOT MSVC)
    target_compile_options(dcsmap PRIVATE -O2 -Wall -Wextra)
endif()
