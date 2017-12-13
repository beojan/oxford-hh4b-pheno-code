-- Code from 'unsup' package: github.com/koraykv/unsup
-- ZCA-Whitening
--
-- Input: 
--  - data tensor M x N1 [x N2 x ...] (required); at least 2D.
--  - means: 1D tensor of size N = N1 x N2 x ... (flattned).
--  - P: ZCA-transfor matrix of size N x N.
--
-- Behavior: 
--  - if both means and P are provided, the ZCA-transformed data is returned, alongside means and P (unchanged). 
--  - otherwise, means and P are computed and returned, preceded by the transformed data. 
--
-- Input arguments are never changed.
--
local unsup = {}
-- PCA using covariance matrix
-- x is supposed to be MxN matrix, where M is the number of samples (trials) and each sample (trial) is N dim
-- returns the eigen values and vectors of the covariance matrix in *INCREASING* order
function unsup.pcacov(x)
   local mean = torch.mean(x,1)
   local xm = x - mean:expandAs(x)
   local c = torch.mm(xm:t(),xm)
   c:div(x:size(1)-1)
   local ce,cv = torch.symeig(c,'V')
   return ce,cv
end
function unsup.zca_whiten(data, means, P, invP, epsilon)
    local epsilon = epsilon or 1e-5
    local auxdata = data:clone()
    local dims = data:size()
    local nsamples = dims[1]
    local n_dimensions = data:nElement() / nsamples
    if data:dim() >= 3 then
       auxdata = auxdata:view(nsamples, n_dimensions)
    end
    if not means or not P or not invP then 
        -- compute mean vector if not provided 
        means = torch.mean(auxdata, 1):squeeze()
        -- compute transformation matrix P if not provided
        local ce, cv = unsup.pcacov(auxdata)
        ce:add(epsilon):sqrt()
        local invce = ce:clone():pow(-1)
        local invdiag = torch.diag(invce)
        P = torch.mm(cv, invdiag)
        P = torch.mm(P, cv:t())

        -- compute inverse of the transformation
        local diag = torch.diag(ce)
        invP = torch.mm(cv, diag)
        invP = torch.mm(invP, cv:t())
    end
    -- remove the means
    local xmeans = means:new():view(1,n_dimensions):expand(nsamples,n_dimensions)
    auxdata:add(-1, xmeans)
    -- transform in ZCA space
    auxdata = torch.mm(auxdata, P)

    auxdata:resizeAs(data)
    return auxdata, means, P, invP
end
return unsup
