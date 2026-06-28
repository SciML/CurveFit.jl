"""
    coef(sol::CurveFitSolution)

Return the fitted coefficients.

The ordering of coefficients depends on the fitting algorithm used.
"""
function StatsAPI.coef(sol::CurveFitSolution)
    return sol.u
end

"""
    residuals(sol::CurveFitSolution; weighted::Bool = true)

Return the residuals of the fitted model.

When the problem was constructed with a `sigma`, residuals are returned as the
weighted form `(y - ŷ) / σ` (i.e. the quantity actually minimized by the
solver). Pass `weighted = false` to instead get the raw `y - ŷ`. Without a
`sigma`, both options return the same thing.
"""
function StatsAPI.residuals(sol::CurveFitSolution; weighted::Bool = true)
    if weighted || isnothing(sol.prob.sigma)
        return sol.resid
    end
    return sol.resid .* sol.prob.sigma
end

"""
    predict(sol::CurveFitSolution, x = sol.prob.x)

Evaluate the fitted model with new data.

If `x` is not provided, predictions are returned at the original data points
used during fitting.
"""
function StatsAPI.predict(sol::CurveFitSolution, x = sol.prob.x)
    return sol(x)
end

"""
    fitted(sol::CurveFitSolution)

Return the fitted values at the original data points.
"""
function StatsAPI.fitted(sol::CurveFitSolution)
    return sol(sol.prob.x)
end

"""
    nobs(sol::CurveFitSolution)

Return the number of observations used in the fit.
"""
function StatsAPI.nobs(sol::CurveFitSolution)
    return length(sol.prob.y)
end

# Return the positions in `sol.u` that are held fixed. These coefficients
# are still reported by `coef` for a consistent parameter layout, but they don't
# have a degree of freedom and contribute a zero row/column to the covariance.
fixed_param_indices(::CurveFitSolution) = ()
function fixed_param_indices(sol::CurveFitSolution{<:ExpSumFitAlgorithm})
    # `k` is always stored as the first coefficient, but it's only fitted when
    # `withconst` is true; otherwise it's held at 0.
    return sol.alg.withconst ? () : (1,)
end

"""
    dof(sol::CurveFitSolution)

Return the number of degrees of freedom of the model.

Note that this counts the fitted coefficients. Coefficients that are held fixed
(e.g. the constant of [`ExpSumFitAlgorithm`](@ref) with `withconst = false`) are
excluded.
"""
function StatsAPI.dof(sol::CurveFitSolution)
    return length(sol.u) - length(fixed_param_indices(sol))
end

"""
    dof_residual(sol::CurveFitSolution)

Return the residual degrees of freedom.

This is defined as `nobs(sol) - dof(sol)`.
"""
function StatsAPI.dof_residual(sol::CurveFitSolution)
    return nobs(sol) - dof(sol)
end

"""
    rss(sol::CurveFitSolution; weighted::Bool = true)

Return the residual sum of squares (RSS), defined as `sum(abs2, residuals(sol; weighted))`.

When `sigma` was provided to the problem, the default weighted RSS is χ². Pass
`weighted = false` to get the unweighted `sum(abs2, y - ŷ)` instead.
"""
function StatsAPI.rss(sol::CurveFitSolution; weighted::Bool = true)
    return sum(abs2, residuals(sol; weighted))
end

"""
    mse(sol::CurveFitSolution; weighted::Bool = true)

Return the [mean squared error](https://en.wikipedia.org/wiki/Mean_squared_error) of the fit.

Computed as `rss(sol; weighted) / dof_residual(sol)`. With a `sigma` on the
problem this is the reduced χ²; pass `weighted = false` for the plain MSE.
"""
function mse(sol::CurveFitSolution; weighted::Bool = true)
    return rss(sol; weighted) / dof_residual(sol)
end

function jacobian(sol::CurveFitSolution{<:LinearCurveFitAlgorithm})
    x = sol.prob.x
    alg = sol.alg
    if alg.yfun === identity
        # Linear model, compute Jacobian analytically
        J = Matrix{eltype(x)}(undef, length(x), 2)
        J[:, 1] .= alg.xfun.(x)
        J[:, 2] .= 1
        return J
    else
        # When there's a y-transform applied (e.g. power/exp fits) the residuals
        # are in the original y-space so it would be incorrect to use the linear
        # analytical Jacobian, instead we use DI.
        a, b_stored = sol.u
        u = [a, b_stored]
        f_pred = u_curr -> begin
            alg.yfun_inverse.(alg.yfun(u_curr[2]) .+ u_curr[1] .* alg.xfun.(x))
        end
        return DifferentiationInterface.jacobian(f_pred, AutoForwardDiff(), u)
    end
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

    function model_expsum(u_curr, x_val)
        # `sol.u` always stores the constant first, regardless of `withconst`
        # (when `withconst` is false it's just held at 0). Layout: k, p (n), λ (n).
        k = u_curr[1]
        p = view(u_curr, 2:(n + 1))
        λ = view(u_curr, (n + 2):(2n + 1))

        # Computation: k + sum(p .* exp.(λ .* x))
        # Use sum generator to avoid allocation
        return k + sum(p[i] * exp(λ[i] * x_val) for i in 1:n)
    end

    range = 1:length(x)
    f_pred = u_curr -> map(i -> model_expsum(u_curr, x[i]), range)

    return DifferentiationInterface.jacobian(f_pred, AutoForwardDiff(), u)
end

function jacobian(sol::CurveFitSolution{<:KingCurveFitAlgorithm})
    # King's law: U = t² where t = (E²−A)/B, so ∂U/∂A = −2t/B and ∂U/∂B = −2t²/B
    A, B = sol.u
    x = sol.prob.x
    J = Matrix{eltype(x)}(undef, length(x), 2)
    @inbounds for i in eachindex(x)
        t = (x[i]^2 - A) / B
        J[i, 1] = -2t / B    # ∂U/∂A
        J[i, 2] = -2t^2 / B  # ∂U/∂B
    end
    return J
end

function jacobian(sol::CurveFitSolution{<:ModifiedKingCurveFitAlgorithm})
    # Modified King: E^2 = A + B * U^n
    # x corresponds to E (Voltage) - Model predicts x^2
    # y corresponds to U (Velocity)
    # Model: x^2 = A + B * y^n
    # This is an implicit model `f(x, y, p) = 0` or `g(y) = h(x)`.
    # Our fitting minimizes (A + B*y^n - x^2)^2.
    # The "residual" is (Prediction - Observation).
    # Observation is x^2. Prediction is A + B*y^n.
    # We need Jacobian of Prediction w.r.t parameters.

    u = sol.u
    x = sol.prob.x # E
    y = sol.prob.y # U

    # Check if y is valid (it should be)
    @assert y !== nothing "Modified King fit requires valid `y` data (Velocity)"

    function model_mod_king(u_curr, y_val)
        A = u_curr[1]
        B = u_curr[2]
        n = u_curr[3]
        return A + B * y_val^n
    end

    range = 1:length(y)
    f_pred = u_curr -> map(i -> model_mod_king(u_curr, y[i]), range)

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

"""
    isconverged(sol::CurveFitSolution)

Return `true` if the underlying solver successfully converged.

This is determined from the solver return code.
"""
function isconverged(sol::CurveFitSolution)
    return SciMLBase.successful_retcode(sol)
end


"""
    vcov(sol::CurveFitSolution; absolute_sigma::Bool = false)

Return the variance–covariance matrix of the fitted coefficients.

The covariance matrix is computed via QR on the (weighted) Jacobian. When
`sigma` was provided to the problem, the Jacobian is scaled by `1 / σ` so the
result is `(JᵀWJ)⁻¹` with `W = diag(1 / σ²)`.

The covariance is then rescaled by reduced χ² (i.e. `mse(sol)`) so that `sigma`
acts as a relative weight. Pass `absolute_sigma = true` to skip this rescaling
when `sigma` carries absolute physical uncertainties (analogous to scipy's
`curve_fit(absolute_sigma=true)`).

!!! note "Fits estimated in a transformed space"
    Some fits estimate their coefficients by ordinary least squares (OLS) in a
    *transformed* space rather than in original `y`-space:

    - [`PowerCurveFitAlgorithm`](@ref) and [`ExpCurveFitAlgorithm`](@ref) fit a
      straight line in log space, minimizing relative rather than absolute error.
    - [`KingCurveFitAlgorithm`](@ref) fits `sqrt(U) = b + a·E²` by OLS and
      recovers `(A, B)` from `(b, a)`. [`ModifiedKingCurveFitAlgorithm`](@ref)
      is not affected since it's a regular nonlinear fit.
    - Any [`LinearCurveFitAlgorithm`](@ref) with `yfun ≠ identity`.

    To keep the reported uncertainty consistent with how the coefficients were
    actually estimated, the covariance is computed in that transformed space and
    then mapped back to the original parameters via the delta method, rather than
    from the original-space Jacobian used above.

    Note that [`residuals`](@ref), [`rss`](@ref), and [`mse`](@ref) remain in
    original `y`-space as fit diagnostics, so for these fits `mse(sol)` is *not*
    the variance scaling behind `vcov`.
"""
# (JᵀJ)⁻¹ via QR, which is more numerically stable than inv(J'J).
function _vcov_from_jacobian(J)
    R = LinearAlgebra.qr(J).R
    Rinv = inv(R)
    return Rinv * Rinv'
end

function StatsAPI.vcov(sol::CurveFitSolution; absolute_sigma::Bool = false)
    J = jacobian(sol)

    if !isnothing(sol.prob.sigma)
        J ./= sol.prob.sigma
    end

    fixed = fixed_param_indices(sol)
    if isempty(fixed)
        covar = _vcov_from_jacobian(J)
    else
        # Fixed coefficients are known exactly: drop them from the inversion
        # (their Jacobian columns would make JᵀJ singular) and scatter the
        # free-parameter covariance back, leaving zero rows/columns for them.
        free = setdiff(axes(J, 2), fixed)
        covar_free = _vcov_from_jacobian(J[:, free])
        covar = zeros(eltype(covar_free), size(J, 2), size(J, 2))
        covar[free, free] .= covar_free
    end

    if !absolute_sigma
        covar .*= mse(sol)
    end

    return covar
end

# Covariance for fits whose coefficients are estimated by ordinary least squares
# in a *transformed* space but whose parameters are reported in original
# space. Computing `vcov` from the original-space Jacobian (as the generic
# method does) gives a Gauss–Newton covariance of a different, nonlinear fit,
# not the uncertainty of the estimator that actually produced `sol.u`. Instead
# we form the transformed-space least-squares covariance `σ²·(XᵀX)⁻¹` for the
# straight line `Y = β0 + β1·t` and delta-method it to the stored
# parameterization via `param_map(β0, β1) -> [p1, p2]`.
function _transformed_ols_vcov(t, Y, β0, β1, param_map, absolute_sigma::Bool)
    n = length(t)
    T = float(promote_type(eltype(t), eltype(Y)))
    st = sum(t)
    stt = sum(abs2, t)
    det = n * stt - st^2

    # Residual variance in transformed space (sigma is unsupported for these
    # fits, so this is the unweighted estimate). absolute_sigma assumes unit
    # variance instead, mirroring the generic method.
    σ2 = if absolute_sigma
        one(T)
    else
        rss_t = sum(i -> abs2(Y[i] - (β0 + β1 * t[i])), eachindex(t))
        rss_t / (n - 2)
    end

    # σ²·(XᵀX)⁻¹ for (β0 intercept, β1 slope)
    Σβ = Matrix{T}(undef, 2, 2)
    Σβ[1, 1] = σ2 * stt / det          # Var(β0)
    Σβ[2, 2] = σ2 * n / det            # Var(β1)
    Σβ[1, 2] = Σβ[2, 1] = -σ2 * st / det  # Cov(β0, β1)

    # Delta method: G = ∂(p1, p2)/∂(β0, β1); Cov(p) = G·Σβ·Gᵀ.
    g = β -> param_map(β[1], β[2])
    G = DifferentiationInterface.jacobian(g, AutoForwardDiff(), [β0, β1])
    return G * Σβ * G'
end

# vcov() for linear fits of a transformed model, falls back to the default
# implementation if there's no y-transform.
function StatsAPI.vcov(sol::CurveFitSolution{<:LinearCurveFitAlgorithm}; absolute_sigma::Bool = false)
    alg = sol.alg
    if alg.yfun === identity
        return @invoke StatsAPI.vcov(sol::CurveFitSolution; absolute_sigma)
    end

    a, b_stored = sol.u
    t = alg.xfun.(sol.prob.x)
    Y = alg.yfun.(sol.prob.y)
    β0 = alg.yfun(b_stored)  # intercept
    β1 = a                   # slope
    # Generic in yfun_inverse so custom inverses work.
    param_map = (β0, β1) -> [β1, alg.yfun_inverse(β0)]
    return _transformed_ols_vcov(t, Y, β0, β1, param_map, absolute_sigma)
end

# King's law: OLS of `sqrt(U) = b + a·E²` (intercept b, slope a), with stored
# params (A, B) = (-b/a, 1/a). Same transformed-space estimator as power/exp.
function StatsAPI.vcov(sol::CurveFitSolution{<:KingCurveFitAlgorithm}; absolute_sigma::Bool = false)
    A, B = sol.u
    t = abs2.(sol.prob.x)   # E²
    Y = sqrt.(sol.prob.y)   # sqrt(U)
    a = 1 / B               # slope
    b = -A / B              # intercept (A = -b·B ⇒ b = -A/B)
    param_map = (β0, β1) -> [-β0 / β1, 1 / β1]   # (A, B)
    return _transformed_ols_vcov(t, Y, b, a, param_map, absolute_sigma)
end

"""
    stderror(sol::CurveFitSolution; absolute_sigma = false, rtol = NaN, atol = 0)

Return the standard errors of the fitted coefficients.

Standard errors are computed as the square roots of the diagonal elements of the
variance–covariance matrix. See [`vcov`](@ref) for the meaning of `absolute_sigma`.
"""
function StatsAPI.stderror(sol::CurveFitSolution; absolute_sigma::Bool = false, rtol::Real = NaN, atol::Real = 0)
    covar = vcov(sol; absolute_sigma)
    vars = LinearAlgebra.diag(covar)

    # Safety check from LsqFit.jl
    vratio = minimum(vars) / maximum(vars)
    if !isapprox(
            vratio,
            0.0,
            atol = atol,
            rtol = isnan(rtol) ? Base.rtoldefault(vratio, 0.0, 0) : rtol,
        ) && vratio < 0.0
        error("Covariance matrix is negative for atol=$atol and rtol=$rtol")
    end

    return sqrt.(abs.(vars))
end

"""
    margin_error(sol::CurveFitSolution, alpha = 0.05; absolute_sigma = false, rtol = NaN, atol = 0)

Returns the margin of error of the fitted coefficients, computed as
`stderror(sol) * t` where `t` is the critical value of the t-distribution for `1 - alpha / 2`.
See [`vcov`](@ref) for the meaning of `absolute_sigma`.
"""
function margin_error(sol::CurveFitSolution, alpha = 0.05; absolute_sigma::Bool = false, rtol::Real = NaN, atol::Real = 0)
    std_errors = stderror(sol; absolute_sigma, rtol, atol)
    dist = TDist(dof_residual(sol))
    critical_values = quantile(dist, 1 - alpha / 2)
    return std_errors * critical_values
end

"""
    confint(sol::CurveFitSolution; level = 0.95, absolute_sigma = false, rtol = NaN, atol = 0)

Return confidence intervals for the fitted parameters.

The confidence intervals are returned as a vector of `(lower, upper)` tuples,
computed as `coef(sol) ± margin_error(sol)`. See [`vcov`](@ref) for the meaning
of `absolute_sigma`.
"""
function StatsAPI.confint(sol::CurveFitSolution; level = 0.95, absolute_sigma::Bool = false, rtol::Real = NaN, atol::Real = 0)
    margin_of_errors = margin_error(sol, 1 - level; absolute_sigma, rtol, atol)
    return collect(zip(coef(sol) .- margin_of_errors, coef(sol) .+ margin_of_errors))
end
