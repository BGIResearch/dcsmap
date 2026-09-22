cmake_minimum_required(VERSION 3.10)
project(dcsmap LANGUAGES CXX VERSION 1.0.2)

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

# Git commit appended to the version string (e.g. "1.0.1-ac999ff").
set(DCSMAP_GIT_COMMIT "unknown")
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE _git_commit
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(_git_commit)
        set(DCSMAP_GIT_COMMIT "${_git_commit}")
        # "-dirty" when tracked files differ from HEAD.
        execute_process(
            COMMAND ${GIT_EXECUTABLE} diff --quiet HEAD --
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE _git_dirty
            ERROR_QUIET)
        if(NOT _git_dirty EQUAL 0)
            set(DCSMAP_GIT_COMMIT "${DCSMAP_GIT_COMMIT}-dirty")
        endif()
    endif()
    # Re-run cmake when HEAD moves or the index changes, so the embedded
    # commit and dirty flag stay in sync.
    foreach(_git_file HEAD index)
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git/${_git_file})
            set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                ${CMAKE_CURRENT_SOURCE_DIR}/.git/${_git_file})
        endif()
    endforeach()
endif()

add_executable(dcsmap
    src/main.cpp
    src/Task.cpp
    src/Scripts.cpp
    src/simp_logger.cpp
)

target_include_directories(dcsmap PRIVATE src)
target_compile_definitions(dcsmap PRIVATE
    DCSMAP_VERSION="${PROJECT_VERSION}"
    DCSMAP_GIT_COMMIT="${DCSMAP_GIT_COMMIT}")
target_link_libraries(dcsmap PRIVATE Threads::Threads)

if(NOT MSVC)
    target_compile_options(dcsmap PRIVATE -O2 -Wall -Wextra)
endif()
