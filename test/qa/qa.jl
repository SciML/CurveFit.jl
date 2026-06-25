using SciMLTesting, CurveFit, Test

run_qa(
    CurveFit;
    explicit_imports = true,
    ei_kwargs = (;
        # Non-public names of upstream packages that CurveFit qualifies; ignored here
        # and expected to become public as those base libraries declare their APIs.
        all_explicit_imports_are_public = (;
            ignore = (
                :AbstractLinearAlgorithm, :AbstractNonlinearAlgorithm,  # SciMLBase
                :init, :solve, :solve!,                                  # CommonSolve
                :coef, :confint, :dof, :dof_residual, :fitted, :nobs,    # StatsAPI
                :predict, :residuals, :rss, :stderror, :vcov,            # StatsAPI
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                :AutoSpecializeCallable, :NonlinearSolvePolyAlgorithmCache,  # NonlinearSolveBase
                :Utils, :get_fu, :clean_sprint_struct,                       # NonlinearSolveBase(.Utils)
                :Success, :T, :successful_retcode,                           # SciMLBase(.ReturnCode)
                :init, :solve!,                                              # CommonSolve
                :coef, :confint, :dof, :dof_residual, :fitted, :nobs,        # StatsAPI
                :predict, :residuals, :rss, :stderror, :vcov,                # StatsAPI
                :rtoldefault,                                               # Base
            ),
        ),
    ),
)
