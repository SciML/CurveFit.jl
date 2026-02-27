using TestItemRunner

# First try to get the group from the argument list, otherwise fall back to the
# environment variable.
# Example command for running the linfit tests: Pkg.test(; test_args=["nonlinfit"])
const GROUP = Symbol(get(ARGS, 1, get(ENV, "GROUP", "All")))

if GROUP == :All
    # Run all tests except those tagged with :nopre
    @run_package_tests filter = ti -> !(:nopre in ti.tags)
else
    @run_package_tests verbose=true filter = ti -> (GROUP in ti.tags)
end
