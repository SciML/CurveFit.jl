using TestItemRunner

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All"
    # Run all tests except those tagged with :nopre
    @run_package_tests filter = ti -> !(:nopre in ti.tags)
elseif GROUP == "nopre"
    # Run only tests tagged with :nopre
    @run_package_tests filter = ti -> (:nopre in ti.tags)
else
    error("Unknown test group: $GROUP")
end
