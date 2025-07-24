# References

## Bibliography

```@bibliography
```

## Key References for Curve Fitting

### Levenberg-Marquardt Algorithm
The Levenberg-Marquardt algorithm [more1978levenberg](@cite) is the workhorse of nonlinear least squares optimization, combining the advantages of gradient descent and Gauss-Newton methods.

### Numerical Methods
For comprehensive coverage of numerical methods in curve fitting, see [press2007numerical](@cite) and [bjorck1996numerical](@cite).

### Statistical Aspects
The statistical theory of nonlinear regression is thoroughly covered in [seber2003nonlinear](@cite).

## Algorithm Implementations

CurveFit.jl implements several key algorithms:

1. **Linear Least Squares**: Using QR decomposition for numerical stability
2. **Levenberg-Marquardt**: Following [gavin2019levenberg](@cite) with improvements from [transtrum2011improvements](@cite)
3. **Trust Region**: Based on the methods described in [madsen2004methods](@cite)

## Related Julia Packages

- [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl): General non-linear least squares fitting
- [GLM.jl](https://github.com/JuliaStats/GLM.jl): Generalized linear models
- [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl): Advanced nonlinear solvers
- [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl): General optimization

## Additional Resources

### Online Resources
- [Julia Documentation](https://docs.julialang.org/)
- [SciML Organization](https://sciml.ai/)
- [JuliaNLSolvers](https://julianlsolvers.github.io/)

### Tutorials and Courses
- [MIT 18.335J: Introduction to Numerical Methods](https://ocw.mit.edu/courses/mathematics/)
- [Curve Fitting in Julia](https://julialang.org/blog/)

### Community
- [Julia Discourse](https://discourse.julialang.org/)
- [Julia Slack](https://julialang.org/slack/)
- [Stack Overflow Julia Tag](https://stackoverflow.com/questions/tagged/julia-lang)