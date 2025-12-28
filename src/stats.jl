#StatsAPI Interface
""" Returns the coefficients of CurveFitSolution"""
StatsAPI.coef(sol::CurveFitSolution) = sol.u

""" Returns the residuals of CurveFitSolution"""
StatsAPI.residuals(sol::CurveFitSolution) = sol.resid

""" Returns the number of independent observations on which the model was fitted(length of the residual)."""
StatsAPI.nobs(sol::CurveFitSolution) = length(sol.resid)

""" Returns the number of degrees of freedom present in CurveFitSolution."""
StatsAPI.dof(sol::CurveFitSolution) = length(coef(sol))

""" Returns the residual degrees of freedom present in CurveFitSolution."""
StatsAPI.dof_residual(sol::CurveFitSolution) = nobs(sol) - dof(sol)

""" Returns the residual sum of squares of CurveFitSolution"""
StatsAPI.rss(sol::CurveFitSolution) = sum(abs2, sol.resid)

""" Returns the mean squared error of CurveFitSolution"""
mse(sol::CurveFitSolution) = rss(sol) / dof(sol)

""" Returns the covariance matrix of the coefficients of CurveFitSolution"""
function StatsAPI.vcov(sol::CurveFitSolution)
    is_nonlinear_problem(sol.prob) || return nothing
    if (SciMLBase.isinplace(sol.prob))
        y = sol.prob.y
        f!(y, p) = sol.prob.nlfunc(y, p, sol.prob.x)
        J = ForwardDiff.jacobian(f!, y, sol.u)
    else
        f(p) = sol.prob.nlfunc(p, sol.prob.x)
        J = ForwardDiff.jacobian(f, sol.u)  
    end
    F = qr!(J)
    R = F.R
    covr = inv(R' * R) * mse(sol)
    return covr       
end    

""" Returns the standard errors for coefficients of CurveFitSolution"""
function StatsAPI.stderror(sol::CurveFitSolution; rtol::Real=NaN, atol::Real=0)
    covar = vcov(sol)
    vars = diag(covar)
    vratio = minimum(vars) / maximum(vars)
    if !isapprox(
        vratio,
        0.0,
        atol=atol,
        rtol=isnan(rtol) ? Base.rtoldefault(vratio, 0.0, 0) : rtol
    ) && vratio < 0.0
        error("Covariance matrix is negative for atol=$atol and rtol=$rtol")
    end
    return sqrt.(abs.(vars))
end

""" Returns the margin of error at alpha significance level of CurveFitSolution """
function margin_of_error(sol::CurveFitSolution, alpha=0.05; rtol::Real=NaN, atol::Real=0)
    std_errors = stderror(sol, rtol=rtol, atol=atol)
    dist = Distributions.TDist(dof(sol))
    crit_val = eltype(coef(sol))(quantile(dist, Float64(1 - alpha/2)))
    return std_errors * crit_val
end    

""" Returns the confidence intervals for the coefficients of CurveFitSolution. Confidence level by default is 95% """
function StatsAPI.confint(sol::CurveFitSolution; level=0.95, rtol::Real=NaN, atol::Real=0)
    margin_err = margin_of_error(sol, 1-level, rtol=rtol, atol=atol)
    return collect(zip(coef(sol) - margin_err, coef(sol) + margin_err))
end 
