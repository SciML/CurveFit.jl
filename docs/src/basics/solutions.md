# Solutions
 
Curve fitting results are returned as `CurveFitSolution` objects. These objects store the
fitted coefficients, residuals, the problem solved and the algorithm used to solve it. It also
stores a return code which holds information about the success of the solver. For more information
on the `ReturnCode.T` type, see [SciMLBase.jl](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/).
`CurveFitSolution` can be evaluated as a callable function on new input data. 

`CurveFitSolution` objects can be treated as statistical models. CurveFit defines methods for 
the `CurveFitSolution` type. To see all of the StatsAPI.jl functions implemented by CurveFit,
and how they are used, see [Tutorial](@ref).
