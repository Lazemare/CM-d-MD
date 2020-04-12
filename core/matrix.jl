function sort_eigens_v2(eigenval::Array, eigenvec::Array)

    # This function sorts eigenvectors by their corresponding eigenvalues, which are generated 
    # by the function 'LinearAlgebra.eigen'.
    order = sortperm(eigenval,rev=true)
    eigenval = eigenval[order]
    eigenvec = eigenvec[:,order]

return eigenval, eigenvec
end

function sort_eigens_abs(eigenval::Array, eigenvec::Array)

    # This function sorts eigenvectors by the absolute values of their corresponding eigenvalues, 
    # which are generated by the function 'LinearAlgebra.eigen'. 
    order = sortperm(broadcast(abs, eigenval),rev=true)
    eigenval = eigenval[order]
    eigenvec = eigenvec[:,order]

return eigenval, eigenvec
end

function spectral_decomposition(T::Array)

    # This model is used to apply the spectral decomposition of transition matrix.
    if size(T,1) != size(T,2)
        error("ERROR! Only be used for transition matrix!")
    end
    eigenval,eigenvec = eigen(T)
    eigenval,eigenvec = sort_eigens_v2(eigenval,eigenvec)
    D_lambda =  Diagonal(eigenval)
    T_target = eigenvec * D_lambda * inv(eigenvec)
    # Test the result.
    for i = 1:size(T,1)^2
        if abs(T_target[i] - T[i]) > 1^-7
            error("ERROR! Bad decomposition results!")
        end
    end

return D_lambda, eigenvec
end

function GHEP_PD(A::Array{Float64,2}, B::Array{Float64,2})

    C = inv(A) * B
    eigenval, eigenvec = eigen(C)
    eigenval = broadcast(real, eigenval)
    eigenvec = broadcast(real, eigenvec)
    eigenval, eigenvec = sort_eigens_abs(eigenval,eigenvec)
    # Normalize the eigenvectors by \mathbf{r}_i \mathbf{A} \mathbf{r}_i^T = 1
    for i = 1:size(A,2)
        c = eigenvec[:,i]' * A * eigenvec[:,i]
        eigenvec[:,i] = eigenvec[:,i] ./ (c ^ 0.5)
    end

return eigenval, eigenvec
end

function GaussianDistritution(dim1::Int=1,dim2::Int=1)

    function G(a::Float64)

        rsq = 0.0; v1 = 0; v2 = 0.0
        while rsq >= 1.0 || rsq == 0.0
            v1 = 2.0 * rand() - 1.0
            v2 = 2.0 * rand() - 1.0
            rsq = v1 ^ 2 + v2 ^ 2
        end
        fac = (-2.0 * log10(rsq) / rsq) ^ 0.5
        gset = v1 * fac
        iset = 1
        output = v2 * fac

    return output
    end

    out = zeros(Float64,dim1,dim2)
    out = broadcast(G,out)

return out
end

function svd_randn(inp::Array{Float64,2};cutoff::Union{Float64,Int})
    
    if typeof(cutoff) == Float64
        if cutoff <= 1.0
            cutoff = (minimum(size(inp)) * cutoff) |> floor |> Int
        else 
            @warn "Your cutoff is larger than 1.0, thus has been set to 1.0!"
            @warn "If you want to assign the rank of the output matrix, use a Int cutoff value."
            cutoff = minimum(size(inp))
        end
    elseif cutoff == 1
        cutoff = minimum(size(inp))
    end
    omega = GaussianDistritution(size(inp,2),cutoff)
    Y = inp * omega
    F = qr(Y)
    B = (F.Q)' * inp
    F_ = svd(B)
    U = F.Q * F_.U

return U, F_.S, F_.Vt
end

