G(r,x) = r == 0 ? x : log(1+r*x)/r
invG(r,x) = r == 0 ? x : (exp(r*x)-1)/r
fun(t) = 0.2*t[1]^2
jud0(x) = x == 0
judno0(x) = x != 0
jud1(x) = x == 1
jud2(x) = x == 2
judall(x,a) = x == a
function Data_generation(seed,N_S,exact_rate,true_beta,fun1,r,W,tau)
  Random.seed!(seed)
  N = size(W)[1]
  area = repeat([i for i=1:N],N_S)
  location = unique(area)
  n = N_S * N
  x = hcat(rand(Binomial(1,0.5),n),rand(Normal(0,1),n))
  obv_point = 1 .+ rand(Poisson(2),n)
  TransMatrix = hcat([(i == j ? 1 : 0) for i = 1:N-1,j = 1:N-1],repeat([-1.0],N-1))
  A = [(i == j ? sum(W[:,i]) : 0) for i=1:N,j=1:N]
  Sigmamatrix = Matrix(Hermitian(inv(TransMatrix * Matrix(A.-W) * TransMatrix')./tau))
  phicof_true = rand(MvNormal(Sigmamatrix),1)
  phicof_true = vcat(phicof_true,-sum(phicof_true))
  phicof = repeat([0.0],n)
  for i = 1:n
    for j = 1:N
      if area[i] == location[j]
        phicof[i] = phicof_true[j]
      end 
    end
  end
  linPred = exp.(x * true_beta .+ phicof)
  u1 = rand(Uniform(0,1),n)
  tt = invG.(r, -log.(1 .- u1)) ./ linPred
  function f!(F, t)
    F[1] = fun1(t) - tt[1]
  end
  trueTimes = nlsolve(f!, [1.0]).zero
  global  k = 2
  for i = 2:n
    global num = tt[k]
    function g!(F, x)
      F[1] = fun1(x) - num
    end
    trueTimes = vcat(trueTimes,nlsolve(g!, [1.0]).zero)
    k = k + 1
  end
  l = repeat([0.0],n)
  u = repeat([0.0],n)
  delta = repeat([2],n)
  ind = unique(sort(rand(DiscreteUniform(1, n),Int(round(n*exact_rate,digits = 0)))))
  for i = 1:n
    if (i in ind) + (trueTimes[i]<10) == 2
      l[i] = trueTimes[i]
      u[i] = trueTimes[i]
      delta[i] = 0
    else
      time_lag = rand(Exponential(1),obv_point[i])
      time_seq = vcat(0,cumsum(time_lag),Inf)
      if trueTimes[i] == 0
        l[i] = time_seq[1]
        u[i] = time_seq[2]
      end
      for j in 2:length(time_seq)
        if time_seq[j-1] < trueTimes[i] <= time_seq[j]
          l[i] = time_seq[j-1]
          u[i] = time_seq[j]
        end
      end
    end
    if u[i] == Inf 
      delta[i] = 3
    end
    if l[i] == 0
      delta[i] = 1
    end
  end
  data = hcat(x,l,u,delta,area)
  return(data)
end
function Ispline(x, order, knots)
  k = order + 1  
  m = length(knots)  
  n = m - 2 + k    
  t = vcat(repeat([1],k) * knots[1], knots[2:(m - 1)], repeat([1],k) * knots[m]) 
  yy_1 = begin
    local yy1 = zeros(n + k - 1, length(x)) 
    for l = k:n 
      yy1[l,:] = (t[l] .<=  x .< t[l + 1])/(t[l + 1] - t[l])
    end
    yy1
  end
  yytem_1 = begin
    local yytem1 = yy_1
    for ii = 1:order
      yytem1 = begin
        local yytem2 = zeros(n + k - 1 - ii, length(x))
        for i = (k - ii):n
          yytem2[i,:] = (ii + 1) .* ((x .- t[i]) .* yytem1[i, :] + (t[i + ii + 1] .- x) .* yytem1[i + 1,: ])./(t[i + ii + 1] - t[i])./ii
        end
        yytem2
      end
    end
    yytem1
  end
  index = begin
    local index_1 = zeros(length(x))
    for i = 1:length(x) 
      index_1[i] = sum(t .<= x[i])
    end
    index_1
  end 
  if order == 1 
    yy2 = begin 
      local yy = zeros(n - 1, length(x))
      for i = 2:n
        yy[i - 1, :] = (i .< index .- order .+ 1) .+ (i == index) .* (t[i + order + 1] - t[i]) .* yytem_1[i, :]./(order + 1)
      end
      yy
    end
  else
    yy2 = begin
      local yy = zeros(n - 1, length(x))
      for j = 1:length(x)
        for i = 2:n
          if i < (index[j] - order + 1)
            yy[i - 1, j] = 1
          elseif index[j] >= i >= (index[j] - order + 1)
            yy[i - 1, j] = (t[(i + order + 1):Int(index[j] + order + 1)] - t[i:Int(index[j])])' * yytem_1[i:Int(index[j]), j]/(order + 1)
          else
            yy[i - 1, j] = 0
          end
        end
      end
      yy
    end
  end
  return(yy2)
end
function Mspline(x, order, knots)
  k = order 
  m = length(knots)  
  n = m - 2 + k    
  t = vcat(repeat([1],k) * knots[1], knots[2:(m - 1)], repeat([1],k) * knots[m]) 
  yy_1 = begin
    local yy1 = zeros(n + k - 1, length(x)) 
    for l = k:n 
      yy1[l,:] = (t[l] .<=  x .< t[l + 1])/(t[l + 1] - t[l])
    end
    yy1
  end
  if order == 1
    yytem_1 = yy1
  else
    yytem_1 = begin
      local yytem1 = yy_1
      for ii = 1:(order-1)
        yytem1 = begin
          local yytem2 = zeros(n + k - 1 - ii, length(x))
          for i = (k - ii):n
            yytem2[i,:] = (ii + 1) .* ((x .- t[i]) .* yytem1[i, :] + (t[i + ii + 1] .- x) .* yytem1[i + 1,: ])./(t[i + ii + 1] - t[i])./ii
          end
          yytem2
        end
      end
      yytem1
    end
  end
  return yytem_1
end
function poissrnpositive(lambda)
  q = 200
  t = [i for i = 0:(q+1)]
  p = [pdf(Poisson(lambda), t[i]) for i = 1:length(t)]
  global pp = cumsum(p[2:(q + 1)])/(1 - p[1])
  u = rand(1)
  while u[1] > pp[q] 
    q = q + 1
    pp = vcat(pp, pp[q - 1] + pdf(Poisson(lambda), q)/(1 - p[1]))
  end
  ll = sum(u[1] .> pp) + 1
  return(ll)
end
function log_Posterior(beta, xcov, z, te2, spatialcof,frailtycof)
  sum1 = sum(z .* xcov * beta)
  sum2 = sum(te2 .* exp.(xcov * beta .+ spatialcof) .* frailtycof)
  sum3 = sum(beta.^2/(2*sigma^2))
  result = sum1-sum2-sum3
  return result
end
function Hessian_beta(beta, xcov, z, te2, spatialcof,frailtycof) 
  hessian1 = xcov' * (xcov .* (te2 .* exp.(xcov * beta + spatialcof) .* frailtycof))
  hessian2 = [(i == j ? sigma^(-2) : 0) for i = 1:length(beta), j = 1:length(beta)]
  result = hessian1 + hessian2
  return(result)
end
function beta_sampler(beta, xcov, z, te2, spatialcof, frailtycof, parbeta, iter, accept_rate)
  log_ini = log_Posterior(beta, xcov, z, te2, spatialcof, frailtycof)
  H_beta = iter < MH_iter ? Hessian_beta(beta, xcov, z, te2, spatialcof, frailtycof) : cov(parbeta[1000:(iter - 1),:])
  H_beta_temp = Matrix(Hermitian(c0^2*inv(H_beta)))
  beta_prop = rand(MvNormal(beta,H_beta_temp),1)
  log_prop = log_Posterior(beta_prop, xcov, z, te2, spatialcof, frailtycof)
  H_beta_prop = iter < MH_iter ? Hessian_beta(beta_prop, xcov, z, te2, spatialcof, frailtycof) : H_beta
  H_beta_prop_temp = Matrix(Hermitian(c0^2*inv(H_beta_prop)))
  logprop_ini = logpdf(MvNormal(beta,H_beta_temp), beta_prop[:,1])
  logini_prop = logpdf(MvNormal(beta_prop[:,1],H_beta_prop_temp),beta)
  
  logR = log_prop .- log_ini .- logprop_ini .+ logini_prop 

  if rand(Uniform(0,1),1)[1] < exp(logR)[1]
    beta = beta_prop[:,1]
    accept_rate = accept_rate + 1
  end
  res = Dict("beta" => beta, "accept_rate" => accept_rate)
  return(res)
end
function spatialcof_sampler(cof, j, z, w, te, tau, parcof,accept_rate,iter)
  w[j] = -sum(w)
  log_ini = z * cof[j] - te * exp(cof[j]) - 0.5 * w[j] * tau * ((w' * cof)[1]/w[j])^2
  H_cof = iter < MH_iter_spatial ? 1/sqrt(w[j] * tau + te * exp(cof[j])) : sqrt(cov(parcof[1000:(iter - 1)]))
  cof_prop = copy(cof)
  cof_prop[j] = rand(Normal(cof[j],c1 ^2 * H_cof),1)[1]
  log_prop = z * cof_prop[j] - te * exp(cof_prop[j]) - 0.5 * w[j] * tau * ((w' * cof_prop)[1]/w[j])^2
  H_cof_prop = iter < MH_iter_spatial ? 1/sqrt(w[j] * tau + te * exp(cof_prop[j])) : H_cof
  logprop_ini = logpdf(Normal(cof[j],c1 ^2 * H_cof), cof_prop[j])
  logini_prop = logpdf(Normal(cof_prop[j],c1 ^2 * H_cof_prop),cof[j])
  
  logR = log_prop .- log_ini .- logprop_ini .+ logini_prop 
  
  if rand(Uniform(0,1),1)[1] < exp(logR)[1]
    accept_rate = accept_rate + 1
    res = Dict("cof" => cof_prop[j], "accept_rate" => accept_rate)
  else
    res = Dict("cof" => cof[j], "accept_rate" => accept_rate)
  end
  return(res)
end
function loglikelihood(beta,gamma,spatial,r,bisL,bisR,xcov,n,status)
  LambdaL = bisL' * gamma
  LambdaR = bisR' * gamma
  temp = exp.(xcov * beta .+ spatial)
  tempL = LambdaL .* temp
  tempR = LambdaR .* temp
  constantL = G.(r,  tempL)
  constantR = G.(r,  tempR)
  FL = 1 .- exp.(- constantL)
  FR = 1 .- exp.(- constantR)
  ft = tempR .* exp.(-(r+1) .* constantL)
  F_1 = zeros(n)
  for i in 1:n
    if status[i] == 1
        F_1[i] = FR[i]
    end
    if status[i] == 2
      F_1[i] = FR[i]-FL[i]
    end
    if status[i] == 3
      F_1[i] = 1-FL[i]
    end
    if status[i] == 0
      F_1[i] = ft[i]
    end
  end
  return(F_1)
end
function est(r, order_1, n_knots, niter, a_eta, b_eta, a_lambda, b_lambda, data)
  xcov = Matrix(data[:,1:2])
  L = data[:,3]
  R = data[:,4]
  p = size(xcov)[2]
  status = data[:,5]
  area = data[:,6]
  location = unique(area)
  n = length(L)
  N = length(location)
  ind = [(R[i] != Inf ? R[i] : min(R...)) for i = 1:n]
  len = max(vcat(L,ind)...) + 0.1 - min(vcat(L,ind)...)
  knots = [min(vcat(L,ind)...) + i * len/(n_knots-1) for i = 0:(n_knots-1)]
  k = length(knots) - 2 + order_1
  t = R[findall(jud0,status)]
  R_true = R[findall(judno0,status)]
  bmst = Mspline(t, order_1, knots)
  bisL = Ispline(L, order_1, knots)
  bisR_true = Ispline(R_true, order_1, knots)
  bisR = zeros(k,n)
  bisR[:,findall(jud0,status)] = bmst
  bisR[:,findall(judno0,status)] = bisR_true
  te3 = copy(bisL)
  ind1 = findall(jud1,status) 
  ind2 = findall(jud2,status) 
  te3[:,ind1] = bisR[:,ind1]
  te3[:,ind2] = bisR[:,ind2]
  eta = rand(Gamma(a_eta, 1/b_eta),1)[1]
  tau = rand(Gamma(a_tau, 1/b_tau),1)[1]
  gamcoef = rand(Gamma(1,1),k) 
  LambdaL = bisL' * gamcoef
  LambdaR = bisR' * gamcoef 
  te2 = te3' * gamcoef
  beta = ones(size(xcov)[2]) .* 0.1#rand(Uniform(-1,1),p)#rand(MvNormal([(i==j ? 1 : 0) for i = 1:p,j = 1:p]))[:,1]#ones(size(xcov)[2]) .* 0.1
  frailtycof = r == 0 ? repeat([1.0],n) : rand(Gamma(1/r,r),n)
  #Sigmamatrix = [(i == j ? 1 : 0) for i = 1:(N-1), j = 1:(N-1)]
  spatialcof_true = repeat([0.0],N)#rand(Normal(0,1),N-1)
  #spatialcof_true = vcat(spatialcof_true,-sum(spatialcof_true))
  spatialcof = repeat([0.0],n)
  for i = 1:n
    for j = 1:N
      if area[i] == location[j]
        spatialcof[i] = spatialcof_true[j]
      end 
    end
  end
  A = [(i == j ? sum(W[:,i]) : 0) for i=1:N,j=1:N]
  parbeta = zeros(niter, p)
  pargam = zeros(niter, k)
  partau = zeros(niter, 1)
  parspatial = zeros(niter, N)
  parspatial_all = zeros(niter, n)
  iter = 1
  accept_rate = 0.0
  accept_rate_spatial = zeros(N)
  while iter < niter + 1
    #iter % 100 == 0 ? println(iter) : iter
    z = zeros(n, 1)
    zz = zeros(n, k)
    te1 = exp.(xcov * beta .+ spatialcof) .* frailtycof
    for i = 1:n
      if status[i] == 1.0
        templam1 = LambdaR[i] * te1[i]
        z[i] = poissrnpositive(templam1)[1]
        zz[i,:] = rand(Multinomial(Int(z[i]), gamcoef .* bisR[:,i]/sum(gamcoef .* bisR[:,i])),1)
      elseif status[i] == 2.0
        templam1 = (LambdaR[i]-LambdaL[i]) * te1[i]
        z[i] = poissrnpositive(templam1)[1]
        zz[i,:] = rand(Multinomial(Int(z[i]), gamcoef .* (bisR[:,i]-bisL[:,i])/sum(gamcoef .* (bisR[:,i]-bisL[:,i]))),1)
      elseif status[i] == 0.0
        z[i] = 1
        zz[i,:] = rand(Multinomial(Int(z[i]), gamcoef .* bisR[:,i]/sum(gamcoef .* bisR[:,i])),1)
      end
    end
    beta_res = beta_sampler(beta, xcov, z, te2, spatialcof, frailtycof, parbeta, iter, accept_rate)
    beta = beta_res["beta"]
    accept_rate = beta_res["accept_rate"]
    te1 = exp.(xcov * beta .+ spatialcof) .* frailtycof
    for l = 1:k
      tempa = 1 + sum(zz[:,l])
      tempb = eta + sum(te3[l,:] .* te1)
      gamcoef[l] = rand(Gamma(tempa,1/tempb),1)[1]
    end
    LambdaL = bisL' * gamcoef
    LambdaR = bisR' * gamcoef
	  te2 = te3' * gamcoef
	  te4 = te2 .* exp.(xcov * beta .+ spatialcof) 
	  if r != 0
	    frailtycof = [(rand(Gamma(z[i]+1/r,1/(te4[i]+1/r)),1))[1] for i = 1:n]
	  end
	  te5 = te2 .* exp.(xcov * beta) .* frailtycof
	  for i = 1:N
      jud(x) = judall(x,location[i])
      ind = findall(jud,area)
	    z_i = sum(z[ind])
	    te5_i = sum(te5[ind])
	    spatialcof_res = spatialcof_sampler(spatialcof_true, i, z_i, -W[:,i], te5_i, tau, parspatial[:,i],accept_rate_spatial[i],iter)
	    spatialcof_true[i] = spatialcof_res["cof"]
	    spatialcof_true = spatialcof_true .- mean(spatialcof_true)
	    accept_rate_spatial[i] = spatialcof_res["accept_rate"]
	    spatialcof[ind] .= spatialcof_true[i]
	  end
    eta = rand(Gamma(a_eta + k, 1/(b_eta + sum(gamcoef))),1)[1]
	  tau = rand(Gamma(0.5 * (N - g) + a_tau,1/(0.5 * (spatialcof_true' * (A .- W) * spatialcof_true)[1] + b_tau)),1)[1]
    pargam[iter,:] = gamcoef
    parbeta[iter,:] = beta
    partau[iter] = tau
    parspatial[iter,:] = spatialcof_true
    parspatial_all[iter,:] = spatialcof
    iter = iter + 1
  end
  #scatter((parbeta[1000:6000,1],parbeta[1000:6000,2]),mc=:red, ma=0.3)
  #histogram(parbeta[5000:10000,1])
  #histogram(parbeta[5000:10000,2])
  #marginalkde(parbeta[1000:6000,1], parbeta[1000:6000,2], clip=((-1.6, 3), (-3, 3)), xlabel=L"X", ylabel=L"Y")
  result = hcat(parbeta[[(1000+2*i) for i = 1:5000],:],partau[[(1000+2*i) for i = 1:5000]])
  chain_trans = Chains(result, [:X1, :X2, :tau])
  result_summary = summarystats(chain_trans)
  result_quantile = quantile(chain_trans)
  lower = result_quantile[:,2]
  upper = result_quantile[:,6]
  estimator = result_summary[:,2]
  
  gamcoef_bar = [(mean(pargam[1000:6000,i])) for i = 1:k]
  spatialcof_bar = [(mean(parspatial_all[1000:6000,i])) for i = 1:n]
  D_thetabar = -2 * sum(log.(loglikelihood(estimator[1:2],gamcoef_bar,spatialcof_bar,r,bisL,bisR,xcov,n,status)))
  likelihood_iter = zeros(5000,n)
  likelihood = zeros(5000)
  for i = 1001:6000
    likelihood_iter[i-1000,:] = loglikelihood(parbeta[i,1:2],pargam[i,:],parspatial_all[i,:],r,bisL,bisR,xcov,n,status)
    likelihood[i-1000] = sum(log.(likelihood_iter[i-1000,:]))
  end
  fin = zeros(n)
  for j = 1:n
    fin[j] = mean(1 ./ likelihood_iter[:,j])
  end
  D_theta = -2 .* mean(likelihood)
  DIC = 2 * D_theta - D_thetabar

  LPML = sum(log.(fin))
  res = Dict("mean" => estimator, "lower" => lower, "upper" => upper, "parvar" => result_summary[:,3]
    , "mcse" => result_summary[:,5], "DIC" => DIC, "LPML" => LPML)
  return res
end