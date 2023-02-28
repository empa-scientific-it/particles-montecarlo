# Tests
# Add/enable GoogleTest library
include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG b796f7d44681514f58a683a3a71ff17c94edb0c1 # v1.13.0
)

FetchContent_MakeAvailable(googletest)

# need to link the test_XXX.cc to "libpart" AND to GoogleTest
# FIXME: create and adjust the real test src file(s)!
add_executable(gtest_run test_evaporation.cc)
target_link_libraries(gtest_run gtest gtest_main part)

# Enable tests discovery
include(GoogleTest)
gtest_discover_tests(gtest_run) # the target must be an executable known to CMake