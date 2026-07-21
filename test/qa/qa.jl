using SciMLTesting, CurveFit, Test

run_qa(
    CurveFit;
    reexports_allow = (
        :coef, :confint, :dof, :dof_residual, :fitted, :init, :nobs, :predict,
        :residuals, :rss, :solve, :solve!, :stderror, :vcov,
    ),
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                :AutoSpecializeCallable, :BoundedWrapper,                    # NonlinearSolveBase
                :NonlinearSolvePolyAlgorithmCache,                           # NonlinearSolveBase
                :Utils, :get_fu, :clean_sprint_struct,                       # NonlinearSolveBase(.Utils)
                :rtoldefault,                                                # Base
            ),
        ),
    ),
)
