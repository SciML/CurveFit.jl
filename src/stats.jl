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

function jacobian(sol::CurveFitSolution{<:RationalPolynomialFitAlgorithm})
    # Rational fit might solve a linear problem y*q(x) = p(x), but vcov
    # requires Jacobian of the actual rational model y = p(x)/q(x).
    # Since prob.nlfunc might be nothing (for linear rational fit),
    # we define the model explicitly and assume u contains [num_coeffs; den_coeffs].
    
    u = sol.u
    x = sol.prob.x
    num_deg = sol.alg.num_degree
    # den_deg = sol.alg.den_degree
    
    # Model function f(u) -> predictions
    function model_rational(u_curr, x_val)
        # u is [num_coeffs..., den_coeffs...]
        # num has num_deg + 1 coeffs
        # den has den_deg coeffs (implicit 1.0 constant term is handled in call, 
        # but sol.u layout depends on implementation. 
        # RationalPolynomialFitAlgorithm (linear) usually assumes:
        # constant term of q(x) is 1. 
        # Let's check call:
        # view(sol.u, 1:(sol.alg.num_degree + 1)) -> numerator
        # vcat(one, view(sol.u, (sol.alg.num_degree + 2):end)) -> denominator
        
        num_c = view(u_curr, 1:(num_deg + 1))
        den_c_params = view(u_curr, (num_deg + 2):length(u_curr))
        
        # We need to construct den with 1.0 at start, but 1.0 is constant.
        # evalpoly requires a vector.
        # Constructing [1.0, den_c_params...] with Duals might trigger allocation/conversion issues.
        # Better to eval explicitly: 1.0 + evalpoly(x*x, den_c_params)*x ? 
        # No, evalpoly(x, [1, c...]) = 1 + c1*x + c2*x^2 ... = 1 + x * evalpoly(x, c)
        
        val_num = evalpoly(x_val, num_c)
        val_den = one(eltype(u_curr)) + x_val * evalpoly(x_val, den_c_params)
        return val_num / val_den
    end
    
    range = 1:length(x)
    f_pred = u_curr -> map(i -> model_rational(u_curr, x[i]), range)
    
    return DifferentiationInterface.jacobian(f_pred, AutoForwardDiff(), u)
end

function jacobian(sol::CurveFitSolution{<:ExpSumFitAlgorithm})
    # ExpSumFitAlgorithm solves y = k + sum(p_i * exp(lambda_i * x))
    # u is a ComponentArray/NamedArrayPartition with (k, p, λ)
    # prob.nlfunc is likely nothing.
    
    u = sol.u
    x = sol.prob.x
    
    # We must access u fields. Since u might be generic array or ComponentArray,
    # we need to handle access carefully matching the call implementation.
    # sol.u has fields :k, :p, :λ if NamedArrayPartition.
    # If using ForwardDiff, u_curr will be a Vector{Dual}.
    # We need to reshape/interpret u_curr based on sol.u structure.
    # But NamedArrayPartition structure isn't preserved in AD usually if passed as vector.
    # We know the sizes from sol.alg (n, m is irrelevant here).
    
    n = sol.alg.n
    withconst = sol.alg.withconst
    
    function model_expsum(u_curr, x_val)
        # Extract parameters from flat vector u_curr
        # Layout: k (if withconst), p (n), λ (n)
        # Check src/expsumfit.jl backing: (; k, p, λ)
        # NamedArrayPartition stores them sequentially.
        
        idx = 1
        if withconst
            k = u_curr[idx]
            idx += 1
        else
            k = zero(eltype(u_curr))
            # k doesn't advance idx
        end
        
        # p is next n
        p = view(u_curr, idx:(idx + n - 1))
        idx += n
        
        # λ is next n
        λ = view(u_curr, idx:(idx + n - 1))
        
        # Computation: k + sum(p .* exp.(λ .* x))
        # Use sum generator to avoid allocation
        return k + sum(p[i] * exp(λ[i] * x_val) for i in 1:n)
    end
    
    range = 1:length(x)
    f_pred = u_curr -> map(i -> model_expsum(u_curr, x[i]), range)
    
    return DifferentiationInterface.jacobian(f_pred, AutoForwardDiff(), u)
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
    # We construct a wrapper for DifferentiationInterface
    
    if SciMLBase.isinplace(sol.prob)
        # In-place: f(resid, u, x)
        f_resid! = (resid, u_curr) -> sol.prob.nlfunc(resid, u_curr, x)
        resid_proto = similar(sol.resid)
        return DifferentiationInterface.jacobian(f_resid!, resid_proto, AutoForwardDiff(), u)
    else
        # Out-of-place: f(u, x) -> resid (or predictions)
        # Note: nlfunc usually returns predictions. sol.resid = pred - y.
        # So d(resid)/du = d(pred)/du.
        f_pred = u_curr -> sol.prob.nlfunc(u_curr, x)
        return DifferentiationInterface.jacobian(f_pred, AutoForwardDiff(), u)
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
