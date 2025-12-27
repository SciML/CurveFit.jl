function coef(sol::CurveFitSolution)
    return sol.u
end

function residuals(sol::CurveFitSolution)
    return sol.resid
end

function predict(sol::CurveFitSolution, x = sol.prob.x)
    return sol(x)
end

function fitted(sol::CurveFitSolution)
    return sol(sol.prob.x)
end

function nobs(sol::CurveFitSolution)
    return length(sol.prob.y)
end

function dof(sol::CurveFitSolution)
    return length(sol.u)
end

function dof_residual(sol::CurveFitSolution)
    return nobs(sol) - dof(sol)
end

function rss(sol::CurveFitSolution)
    return sum(abs2, residuals(sol))
end

function mse(sol::CurveFitSolution)
    return rss(sol) / dof_residual(sol)
end

function vcov(sol::CurveFitSolution)
    error("`vcov` is not yet implemented for `CurveFitSolution`.")
end

function stderror(sol::CurveFitSolution)
    error("`stderror` is not yet implemented for `CurveFitSolution`.")
end
