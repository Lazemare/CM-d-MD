"""
    function cmap(param::Array, tau::Int; reversible::Bool=true,
    reweight::String="default", commute_map::Bool=true)

Calculate the commute map.

# Arguments:
- `param::Array`: a n×m array contains m order parameters with n frames.
- `tau::Int`: the lag time.
- `reversible::Bool=true`: if use the reversible moment matrices.
- `reweight::String="default"`: reweight scheme.
- `commute_map::Bool=true`: use commute map algorithm or TICA algorithm.
"""
function cmap(param::Array,tau::Int;reversible::Bool=true,
              reweight::String="default",commute_map::Bool=true)

    if reversible
        if reweight == "default"
            # Decorrelation of basis functions using PCA
            X = param[1:(size(param,1)-tau),:]
            pi = X' * fill(1.0,size(X,1),1) ./ size(X,1)
            COV = (X' * X) ./ size(X,1) - (pi * pi')
            S_d, Q_d = eigen(COV)
            param_new = (param .- repeat(pi',size(param,1),1)) * Q_d * diagm(abs.(S_d).^0.5)'
            # Update
            X_0 = param_new[1:(size(param_new,1)-tau),:]
            X_tau = param_new[(tau+1):size(param_new,1),:]
            # purge the mean from the coordinates.
            mu_ = 0.5.*(sum(X_0,dims=1) / size(X_0,1) + sum(X_tau,dims=1) / size(X_tau,1))
            mu = repeat(mu_,size(X_0,1),1)
            mu_param = repeat(mu_,size(param_new,1),1)
            X_0 = X_0 .- mu
            X_tau = X_tau .- mu
            param_new = param_new .- mu_param
            # Compute the moment matrices
            C_00 = (X_0' * X_0 + X_tau' * X_tau) ./ (2*size(X_0,1))
            C_0tau = (X_0' * X_tau + X_tau' * X_0) ./ (2*size(X_0,1))
            # Solve the GHEP
            eigenval, eigenvec = GHEP_PD(C_00,C_0tau)
            # Sort the eigenval and eigenvec
            eigenval, eigenvec = sort_eigens_abs(eigenval,eigenvec)
            if commute_map
                # Relaxation Time Scales
                rts = - tau ./ broadcast(log,broadcast(abs,eigenval))
                # Regularization of Time Scales
                factor = broadcast(tanh, (MathConstants.pi .* (rts .- tau) ./ tau .+ 1)) ./ 2
                factor = (broadcast(abs, factor) .+ factor) ./ 2
                rts = rts .* factor
                # The commute map
                CMAP = param_new * eigenvec * Diagonal((rts./2).^0.5)
                return rts, eigenvec * Diagonal((rts./2).^0.5), CMAP
            else 
                # TICA
                CMAP = param_new * eigenvec
                return rts, eigenvec, CMAP
            end
        elseif reweight == "koopman"
            # Decorrelation of basis functions using PCA
            X = param[1:(size(param,1)-tau),:]
            pi = X' * fill(1.0,size(X,1),1) ./ size(X,1)
            COV = (X' * X) ./ size(X,1) - (pi * pi')
            S_d, Q_d = eigen(COV)
            param_new = (param .- repeat(pi',size(param,1),1)) * Q_d * diagm(abs.(S_d).^0.5)'
            param_new = [param_new fill(1.0,size(param_new,1),1)]
            # Update X and Y
            X = param_new[1:(size(param_new,1)-tau),:]
            Y = param_new[(tau+1):size(param_new,1),:]
            # Get rewrighting coefficients
            K = (X' * Y) ./ size(X,1)
            vals, vecs = eigen(K')
            vecs = real.(vecs)
            vals = real.(vals)
            flag = 0
            for i = 1:size(vals,1)
                if vals[i] ≈ 1
                    flag = i
                end
            end
            u = vecs[:,flag]

            # Weight matrix
            W = zeros(size(X,1))
            for i = 1:size(X,1)
                W[i] = X[i,:]'*u
            end
            # Rewright
            pi_eq = 1/2 * (X.+Y)' * diagm(W) * fill(1.0,size(W,1),1)
            COV_eq = 1/2 * (X'*diagm(W)*X + Y'*diagm(W)*Y) - pi_eq*pi_eq'
            S_d_eq, Q_d_eq = eigen(COV_eq)
            param_eq = ((param_new .- repeat(pi_eq',size(param_new,1),1)) * 
                        Q_d_eq * diagm(abs.(S_d_eq).^0.5))
            param_eq = [param_eq fill(1.0,size(param_eq,1),1)]
            # Update X and Y
            X_eq = param_eq[1:(size(param_eq,1)-tau),:]
            Y_eq = param_eq[(tau+1):size(param_eq,1),:]
            C_00   = (X_eq' * X_eq + Y_eq' * Y_eq) ./ (2*size(X_eq,1))
            C_0tau = (X_eq' * Y_eq + Y_eq' * X_eq) ./ (2*size(X_eq,1))
            eigenval, eigenvec = CMdMD.GHEP_PD(C_00,C_0tau)
            # Sort the eigenval and eigenvec
            eigenval, eigenvec = CMdMD.sort_eigens_abs(eigenval,eigenvec)
            # Relaxation Time Scales
            rts = - tau ./ broadcast(log,broadcast(abs,eigenval))
            # Regularization of Time Scales
            factor = broadcast(tanh, (MathConstants.pi .* (rts .- tau) ./ tau .+ 1)) ./ 2
            factor = (broadcast(abs, factor) .+ factor) ./ 2
            rts = rts .* factor
            # The commute map
            CMAP = param_eq * eigenvec * Diagonal((rts./2).^0.5)
            return (rts[2:(size(rts,1)-1)], 
                   (eigenvec * Diagonal((rts./2).^0.5))[:,2:(size(rts,1)-1)], 
                    CMAP[:,2:(size(rts,1)-1)])
        else
            error("Unknown rewrighting method \"$reweight\"!")
        end
    else
        X_0 = param[1:(size(param,1)-tau),:]
        X_tau = param[(tau+1):size(param,1),:]
        # purge the mean from the coordinates.
        mu = sum(X_0,dims=1) / size(X_0,1)
        X_0 = X_0 .- mu
        X_tau = X_tau .- mu
        param = param .- mu
        # Compute the moment matrices
        C_00 = cov(X_0,corrected=false)
        C_0tau = cov(X_0,X_tau,corrected=false)
        # Solve the GHEP
        eigenval, eigenvec = GHEP_PD(C_00,C_0tau)
        # Sort the eigenval and eigenvec
        eigenval, eigenvec = sort_eigens_abs(eigenval,eigenvec)
        if commute_map
            # Relaxation Time Scales
            rts = - tau ./ broadcast(log,broadcast(abs,eigenval))
            # Regularization of Time Scales
            factor = broadcast(tanh, (MathConstants.pi .* (rts .- tau) ./ tau .+ 1)) ./ 2
            factor = (broadcast(abs, factor) .+ factor) ./ 2
            rts = rts .* factor
            # The commute map
            CMAP = param * eigenvec * Diagonal((rts./2).^0.5)
            return rts, eigenvec * Diagonal((rts./2).^0.5), CMAP
        else
            # TICA
            CMAP = param * eigenvec
            return rts, eigenvec, CMAP
        end
    end
end
