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

function jacobian(sol::CurveFitSolution{<:LinearCurveFitAlgorithm})
    x = sol.prob.x
    xfun = sol.alg.xfun
    J = Matrix{eltype(x)}(undef, length(x), 2)
    J[:, 1] .= xfun.(x) # Slope
    J[:, 2] .= 1        # Intercept
    return J
end

function jacobian(sol::CurveFitSolution{<:PolynomialFitAlgorithm})
    x = sol.prob.x
    n = sol.alg.degree
    J = Matrix{eltype(x)}(undef, length(x), n + 1)
    J[:, 1] .= 1
    for i in 1:n
        @. J[:, i + 1] = x^i 
    end
    return J
end

function jacobian(sol::CurveFitSolution)
    # Fallback for nonlinear
    # The residuals are r_i = model(u, x_i) - y_i
    # We need J_ij = dr_i/du_j
    # This is equivalent to d(model)/du since y is constant
    u = sol.u
    x = sol.prob.x
    
    # We need a function f(u) -> residuals
    # CurveFitProblem has nlfunc which is f(u, x) or f(resid, u, x)
    # We construct a wrapper for ForwardDiff
    
    if SciMLBase.isinplace(sol.prob)
        # In-place: f(resid, u, x)
        f_resid! = (resid, u_curr) -> sol.prob.nlfunc(resid, u_curr, x)
        resid_proto = similar(sol.resid)
        # Use ForwardDiff.jacobian(f!, y, x) -> J
        return ForwardDiff.jacobian(f_resid!, resid_proto, u)
    else
        # Out-of-place: f(u, x) -> resid (or predictions)
        # Note: nlfunc usually returns predictions. sol.resid = pred - y.
        # So d(resid)/du = d(pred)/du.
        f_pred = u_curr -> sol.prob.nlfunc(u_curr, x)
        return ForwardDiff.jacobian(f_pred, u)
    end
end


function isconverged(sol::CurveFitSolution)
    return sol.retcode == ReturnCode.Success
end



function vcov(sol::CurveFitSolution)
    J = jacobian(sol)
    
    # Compute the covariance matrix from the QR decomposition
    # This is numerically more stable than inv(J'J)
    Q, R = LinearAlgebra.qr(J)
    
    # Check for rank deficiency or other issues?
    # LinearAlgebra.qr usually handles full rank. 
    # R is upper triangular. Rinv = inv(R)
    
    # Ideally checking rank(R) would be good, but assuming J is full rank for now.
    
    Rinv = inv(R)
    covar = Rinv * Rinv' * mse(sol)
    
    return covar
end



function stderror(sol::CurveFitSolution; rtol::Real=NaN, atol::Real=0)
    covar = vcov(sol)
    vars = LinearAlgebra.diag(covar)
    
    # Safety check from LsqFit.jl
    vratio = minimum(vars) / maximum(vars)
    if !isapprox(
        vratio,
        0.0,
        atol=atol,
        rtol=isnan(rtol) ? Base.rtoldefault(vratio, 0.0, 0) : rtol,
    ) && vratio < 0.0
        error("Covariance matrix is negative for atol=$atol and rtol=$rtol")
    end
    
    return sqrt.(abs.(vars))
end


function margin_error(sol::CurveFitSolution, alpha=0.05; rtol::Real=NaN, atol::Real=0)
    std_errors = stderror(sol; rtol=rtol, atol=atol)
    dist = TDist(dof(sol))
    critical_values = quantile(dist, 1 - alpha / 2)
    return std_errors * critical_values
end

function confint(sol::CurveFitSolution; level=0.95, rtol::Real=NaN, atol::Real=0)
    margin_of_errors = margin_error(sol, 1 - level; rtol=rtol, atol=atol)
    return collect(zip(coef(sol) .- margin_of_errors, coef(sol) .+ margin_of_errors))
end
