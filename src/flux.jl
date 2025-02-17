using Statistics
using StatsBase
using PhysicalConstants
using Unitful
import Loess

const λ = 40660 / 1000 # "J.mmol^-1" latent heat of evaporation of water
const C_p = 29.07 # Molar Heat Capacity at constant pressure J.mol^-1.K^-1
const R = ustrip(PhysicalConstants.CODATA2018.R)

"Convert vector `F` to flux of type `fluxtype` given the density `density`"
function tofluxunits(F, density, fluxtype)
  if fluxtype == :SensibleHeat
    F = F .* density * C_p
    units = u"J/m^2/s"
  elseif fluxtype == :LatentHeat
    F = F .* density * λ
    units = u"J/m^2/s"
  elseif fluxtype == :CO2
    F = F .* density
    units = u"μmol" / u"m^2" / u"s"
  elseif fluxtype == :kinematic
    units = u"m^2" / u"s^2"
  else
    throw(error("wrong flux type"))
  end
  return (F, units)
end

function map_idx(CI, mapping)
  X = Array{Float64}(undef, length(CI))
  Y = Array{Float64}(undef, length(CI))
  Z = Array{Float64}(undef, length(CI))
  for (i, c) in enumerate(CI)
    x, y, z = mapping(c)
    X[i] = x
    Y[i] = y
    Z[i] = z
  end
  return (X, Y, Z)
end

"""
  time_integrate_flux(decomp,mask)

Integrate along the scale domain the flux `decomp` decomposed in the time-scale domain given the integration mask `mask`.

Optionnaly a symbol `fluxtype`can be given to indicate the type of flux such that it is converted to common flux units, see `tofluxunits`.
"""
function time_integrate_flux(decomp, mask)
  size(mask) == size(decomp) || throw(error("Wrong size between mask and decomp"))
  F = Vector{Float64}(undef, size(decomp, 1))
  for j in axes(F, 1)
    values = decomp[j, mask[j, :]]
    if isempty(values)
      F[j] = NaN
    else
      F[j] = sum(decomp[j, mask[j, :]])
    end
  end
  return F
end


"""
    compute_wind_amplitude(wind_speeds,time_params)

Compute the wind amplitude given the three wind speed components stored in `wind_speeds`, with the averaging convolutional kernel parameters `time_params`.
"""
function compute_wind_amplitude(wind_speeds::AbstractArray{<:Real,2}, time_params)
  work_dim = size(wind_speeds)
  (time_sampling, _), kernel_params, KernelC =
    init_averaging_conv_kernel(work_dim; time_params...)
  wind_amplitude = KernelC(wind_speeds, kernel_params)[time_sampling, 1, :]
  wind_amplitude = mapslices(x -> sqrt(sum(abs2, x)), wind_amplitude, dims=2)
  wind_amplitude = dropdims(wind_amplitude, dims=2)
  return wind_amplitude
end


"""
    compute_density(P,T,time_params)

Compute the density given the pressure `P` (Pa) and the temperature `T` (K) with the averaging convolutional kernel parameters `time_params`.
"""
function compute_density(P, T, time_params)
  length(P) == length(T) || throw(error("Signals of different size"))
  work_dim = length(P)
  (time_sampling, _), kernel_params, KernelC =
    init_averaging_conv_kernel(work_dim; time_params...)
  density = P ./ (R * T)
  density = KernelC(density, kernel_params)[time_sampling]
  return density
end


"""
    timescale_flux_decomp(w,θ,time_params,scale_params;with_info=false)

Compute the time-scale decomposition of the flux `wθ` given averaging and wavelet kernel parameters `time_params` and `scale_params`.

# Arguments

  - `w::Vector`: First signal (e.g. the vertical wind speed)
  - `θ::Vector`: Second signal (e.g. temperature)
  - `time_params::NamedTuple`: Named Tuple of the parameters for initializing the averaging convolutional kernel `\\phi`, see `init_averaging_conv_kernel`
  - `scale_params::NamedTuple`: Named Tuple of parameters for initializing the wavelet convolutional kernel `\\psi`, see `init_wave_conv_kernel`
  - `with_info::Bool=false`: Output informations about the decomposition
    - `func=nothing`: Apply function `func` before averaging step
"""
function timescale_flux_decomp(
  w,
  θ,
  time_params,
  scale_params;
  with_info=false,
  func=nothing,
)
  length(w) == length(θ) || throw(error("Signals must be of the same size."))
  work_dim = length(w)
  (freq_peak, σ_waves), wave_params, WaveC =
    init_wave_conv_kernel(work_dim; scale_params..., with_sigma=with_info)
  (time_sampling, σ_averaging), kernel_params, KernelC = init_averaging_conv_kernel(
    (work_dim, size(wave_params, 2) + 1);
    time_params...,
    with_sigma=with_info,
  )
  σ_t = (σ_waves, σ_averaging)
  w_ξ = WaveC(w, wave_params)
  θ_ξ = WaveC(θ, wave_params)
  flux = dropdims(w_ξ .* θ_ξ, dims=3)
  if !(isnothing(func))
    flux = func.(flux)
  end
  flux = KernelC(flux, kernel_params)
  flux = dropdims(flux, dims=2)
  flux = flux[time_sampling, :]
  if with_info
    return (time_sampling, (freq_peak, σ_t), flux)
  else
    return (time_sampling, nothing, flux)
  end
end


"""
  get_timescale_mask(work_dim,σ_waves,σ_averaging,factor=(3,3),max_sigma=false)

Construct a time-scale mask given the time deviation of the wavelet filters and the averaging kernel.

# Arguments

 - `work_dim`: dimension of analysis.
 - `σ_waves`: wavelet filters time-deviation
 - `σ_averaging`: averaging kernel time-deviation
 - `factor=(3,3)`: amounts by which time-deviations of wavelets and filters are multiplied
 - `max_sigma=false`: return a border error mask with the maximum time deviation
"""
function get_timescale_mask(
  work_dim,
  σ_waves,
  σ_averaging,
  factor=(3, 3),
  max_sigma=false,
)
  mask = falses(work_dim, length(σ_waves))
  max_sigma_val = ceil(Int, maximum(σ_waves) * factor[1] + σ_averaging * factor[2])
  for (i, s) in enumerate(σ_waves)
    if max_sigma
      s = max_sigma_val
    else
      s = ceil(Int, s * factor[1] + σ_averaging * factor[2])
    end
    mask[1:s, i] .= true
    mask[(end-s+1):end, i] .= true
  end
  return mask
end


"""
    amplitude_reynolds_w(u,v,w,time_params,scale_params)

Compute the amplitude of the vertical components of the Reynold's tensor using the three wind speed components `u`,`v` and `w` and the time-scale decomposition parameters `time_params` and `scale_params`.

# Arguments

  - `u,v,w::Vector`: wind speed components signal (e.g. the vertical wind speed)
  - `time_params::NamedTuple`: Named Tuple of the parameters for initializing the averaging convolutional kernel `\\phi`, see`init_averaging_conv_kernel`
  - `scale_params::NamedTuple`: Named Tuple of parameters for initializing the wavelet convolutional kernel `\\psi`, see `init_wave_conv_kernel`
"""
function amplitude_reynolds_w(u, v, w, time_params, scale_params)
  _, _, uv = timescale_flux_decomp(u, w, time_params, scale_params)
  _, _, vw = timescale_flux_decomp(v, w, time_params, scale_params)
  time_sampling, (freq_peak, σ_t), ww =
    timescale_flux_decomp(w, w, time_params, scale_params; with_info=true)
  τ_rey_w = sqrt.(uv .^ 2 .+ vw .^ 2 + ww .^ 2)
  return (time_sampling, (freq_peak, σ_t), τ_rey_w)
end

"""
    turbulence_mask_extraction(u,v,w,time_params,scale_params,method,method_params...)

Extract a time-scale mask of the vertical turbulent transport using the three wind speed components `u`,`v` and `w` using the time-scale decomposition parameters `time_params` and `scale_params` (see `init_averaging_conv_kernel` and `init_wave_conv_kernel`) and the turbulence extraction methods `method` with parameters `method_params`.

# Arguments

  - `u,v,w::Vector`: wind speed components signal (e.g. the vertical wind speed)
  - `time_params::NamedTuple`: Named Tuple of the parameters for initializing the averaging convolutional kernel `\\phi`, see`init_averaging_conv_kernel`
  - `scale_params::NamedTuple`: Named Tuple of parameters for initializing the wavelet convolutional kernel `\\psi`, see `init_wave_conv_kernel`
  - `method::Function`: Method used to extract the turbulent transport signal, see `turbu_extract_threshold`,`turbu_extract_laplacian` and `turbu_extract_diffusion`.
"""
function turbulence_mask_extraction(
  u,
  v,
  w,
  time_params,
  scale_params;
  method::Function,
  method_params...,
)
  _, (freq_peak, σ_t), τ_rey_w = amplitude_reynolds_w(u, v, w, time_params, scale_params)
  out = (τ_rey_w, σ_t)
  return method(τ_rey_w; method_params...)
end
function turbulence_mask_extraction(τ_rey_w, ; method::Function, method_params...)
  return method(τ_rey_w; method_params...)
end

function turbu_extract_threshold(τ_rey_w; threshold)
  mask = τ_rey_w .> threshold
  return mask
end

function _locally_weighted_regression(t, eta, span=0.25)
  model = Loess.loess(t, eta, span=span)
  tmin, tmax = extrema(t)
  function g(t)
    if t < tmin
      return Loess.predict(model, tmin)
    elseif t > tmax
      return Loess.predict(model, tmax)
    else
      return Loess.predict(model, t)
    end
  end
  return g
end

function turbu_extract_laplacian(
  t,
  eta,
  log10τ_w;
  δ_Δτ=1,
  δ_τ=1e-3,
  span=0.25,
  mask_error=falses(size(log10τ_w)),
)

  S = size(log10τ_w)
  τ_mapped = view(log10τ_w, :)

  # Reject low-pass filter i.e. at freq_p[end]
  mask = trues(S)
  mask[:, end] .= false

  adj_mat = grid_adj_mat(S, mask) # 9-point grid adjacency matrix with removed vertices from mask
  weights_mat = adj_mat # Using the adjacency matrix as the weight matrix amounts to compute a normal laplacian

  g = MyGraph(adj_mat, weights_mat)
  L = laplacian_matrix(g)
  Δτ = reshape(L * τ_mapped, S) # the laplacian is zero where mask is false
  τ_mapped = reshape(τ_mapped, S)

  detected = reshape(δ_Δτ .< Δτ, S) # Look at (t,eta) points with important minimas
  detected[mask_error] .= false # remove points with convolution errors
  detected[(eta.>0)] .= false # remove points above eta = 0

  #itp = _interpolate_eta(t[detected[:]], eta[detected[:]], λ,d) # Old way: Bspline interpolation + smoothness regularization to get interpolated value at each time t, extrapolate with constant values on the borders
  itp = _locally_weighted_regression(t[detected], eta[detected], span)

  mask_advec = (itp.(t) .< eta) .&& mask# Get the mask removing the advection + removing the mean value
  mask_lowcoeff = (log10(δ_τ) .< τ_mapped) .&& mask_advec
  masks = (; minimas=detected, advection=mask_advec, turbulence=mask_lowcoeff)

  return (masks, Δτ, itp)
end

function turbu_extract_diffusion(
  τ_w;
  time_sampling,
  freq_peak,
  ref_dist=1,
  mean_wind=nothing,
)
  S = size(τ_w)
  CI = CartesianIndices(S)
  # Reject low-pass filter i.e. at freq_p[end]
  mask = trues(S)
  mask[:, end] .= false
  if isnothing(mean_wind)
    vertex_mapping =
      (c::CartesianIndex) -> Float64[
        time_sampling[c[1]],
        log.(freq_peak[c[2]] * ref_dist),
        log(τ_w[c[1], c[2]]),
      ]
  else
    vertex_mapping =
      (c::CartesianIndex) -> Float64[
        time_sampling[c[1]],
        log.(freq_peak[c[2]] * ref_dist / mean_wind[c[1]]),
        log(τ_w[c[1], c[2]]),
      ]
  end
  t, eta, τ_mapped = map_idx(CI, vertex_mapping)
  adj_mat = grid_adj_mat(S, mask)
  σ_τ = std(τ_mapped[mask[:]])
  function weight_func(i::Int, j::Int)
    c_i = CI[i]
    c_j = CI[j]
    v_i = vertex_mapping(c_i)[3]
    v_j = vertex_mapping(c_j)[3]
    v = exp(-(v_i - v_j) / σ_τ) # Asymetric Potential
    return v
  end
  weights_mat = generate_weight_mat(adj_mat, weight_func; normalize=true)
  g = MyGraph(adj_mat, weights_mat)
  tau_rey = (reshape(t, S), reshape(eta, S), reshape(τ_mapped, S))
  return (tau_rey, g)
  #
  #  s=-0.1 .< Y .< 0.1
  #  s=sparse(vec(s))
  #  M=sum(s)
  #  func_acc(s,i)=begin
  #      x=s .+ droptol!(g.weights*s,1e-6)
  #      x=x*(M/sum(x))
  #      return x
  #  end
  #  #all_s=accumulate(func_acc,1:10,init=s);
  #  @warn  println("Not Fully implemented yet")
  #  return (g,Δv,Σ,(X,Y,τ_mapped))
end


"""
    flux_estimation(data,z_d,fs,time_params,scale_params)

Wavelet based estimation of the flux given averaging parameters `time_params` and wavelet parameters `scale_params`.

# Arguments
 - `data::DataFrame`: wind speed, pressure and gas concentrations measurements
 - `z_d::Real`: measurement height above the zeros place displacement height
 - `fs::Integer`: sampling frequency
 - `time_params::NamedTuple`: averaging parameters
 - `scale_params::NamedTuple`: wavelet decomposition parameters
"""
function flux_estimation(
  data;
  z_d,
  fs,
  time_params,
  scale_params,
  time_params_mean_wind=time_params,
  time_params_density=time_params,
  time_params_turbu=time_params,
  freq_tl=(0.1, 1),
  dates=(nothing, nothing),
  with_decomp=false,
  analysis_range=Colon(),
  kwargs...,
)
  sdate, edate = dates
  t0 = time()
  (; U, V, W, T, CO2, H2O, P) = data
  work_dim = length(U)

  contains_nan = any(i -> !isnothing(findfirst(isnan, data[:, i])), 2:size(data, 2))
  if contains_nan
    throw(error("Data contains NaN values. Aborting."))
  end
  T = T .+ 274.15 # °C TO K
  P = 1000 * P # kPa to Pa

  time_h = (0:(work_dim-1))
  time_h = time_h ./ (60 * 60 * fs)
  wind_speeds = hcat(U, V)

  # Wind Amplitude signal and density signal
  mean_wind = compute_wind_amplitude(wind_speeds, time_params_mean_wind)
  density = compute_density(P, T, time_params_density)
  # Time-Lag optimisation
  τ_max = 60 * fs # 1min max timelag search
  τ_arr, corr_H2O = optim_timelag(W, H2O, scale_params, freq_tl, τ_max)
  _, corr_CO2 = optim_timelag(W, CO2, scale_params, freq_tl, τ_max)

  tl_max = -2 * fs # 2s maximum timelag 
  m_tl = τ_arr .<= 0 # search only maximum in negative timelag, gas analyser is always late
  tl_H2O = τ_arr[m_tl][argmax(abs.(corr_H2O[m_tl]))]
  tl_H2O = tl_H2O < tl_max ? 0 : tl_H2O # If we reach a maximum, better not to take it
  tl_CO2 = τ_arr[m_tl][argmax(abs.(corr_CO2[m_tl]))]
  tl_CO2 = tl_CO2 < tl_max ? 0 : tl_CO2
  timelags = Dict(
    :H2O => (tl_H2O, τ_arr, corr_H2O),
    :CO2 => (tl_CO2, τ_arr, corr_CO2))
  circshift!(H2O, timelags[:H2O][1])
  circshift!(CO2, timelags[:CO2][1])

  max_tl = max(abs(timelags[:H2O][1]), abs(timelags[:CO2][1]))

  # Time-Scale Analyses
  time_sampling, (freq_peak, σ_t), decomp_CO2 =
    timescale_flux_decomp(W, CO2, time_params, scale_params; with_info=true)
  _, _, decomp_T = timescale_flux_decomp(W, T, time_params, scale_params)
  _, _, decomp_H2O = timescale_flux_decomp(W, H2O, time_params, scale_params)
  _, _, τ_w = amplitude_reynolds_w(U, V, W, time_params_turbu, scale_params)
  Z = log10.(τ_w)

  σ_t = (σ_t[1] .+ max_tl, σ_t[2]) # add maximum timelag estimated as border error 
  mask_σ_t = get_timescale_mask(work_dim, σ_t..., (4, 4), false)[time_sampling, :] # convolution border errors mask

  to_eta(i_t, j_ξ) = log10(((z_d * freq_peak[j_ξ]) / mean_wind[i_t]))
  S = size(decomp_T) # Dimension
  CI = CartesianIndices(S)

  t = map(c -> time_h[time_sampling[c[1]]], CI) # get the time values
  eta = map(c -> to_eta(c[1], c[2]), CI) # get the normalized freq

  (masks, Δτ, itp) =
    turbu_extract_laplacian(t, eta, Z, δ_Δτ=1, δ_τ=1e-3, mask_error=mask_σ_t)
  mask_minima, mask_noadvection, mask_turbulence = masks
  advec_line = itp.(time_h[time_sampling])
  advec_line = advec_line[analysis_range]

  mask_analysis = falses(size(decomp_T)) # Restriction to period of analysis
  mask_analysis[analysis_range, 1:(end-1)] .= true # We take everything during the analysis range period and without the first frequency peak at 0
  mask_nomean = copy(mask_analysis) .&& .!(mask_σ_t) # Remove border errros
  mask_turbulence = mask_turbulence .&& .!(mask_σ_t) # Remove border errros

  decomp_T,units_T = tofluxunits(decomp_T,density,:SensibleHeat)
  decomp_CO2,units_CO2 = tofluxunits(decomp_CO2,density,:CO2)
  decomp_H2O,units_H2O = tofluxunits(decomp_H2O,density,:H2O)

  F_T_nomean, units_T = time_integrate_flux(decomp_T, mask_nomean)
  F_H2O_nomean, units_H2O =
    time_integrate_flux(decomp_H2O, mask_nomean)
  F_CO2_nomean, units_CO2 = time_integrate_flux(decomp_CO2, mask_nomean)


  F_T_noadvection, units_T =
    time_integrate_flux(decomp_T, mask_noadvection)
  F_H2O_noadvection, units_H2O =
    time_integrate_flux(decomp_H2O, mask_noadvection)
  F_CO2_noadvection, units_CO2 =
    time_integrate_flux(decomp_CO2, mask_noadvection)

  F_T_turbulence, units_T =
    time_integrate_flux(decomp_T, mask_turbulence)
  F_H2O_turbulence, units_H2O =
    time_integrate_flux(decomp_H2O, mask_turbulence)
  F_CO2_turbulence, units_CO2 =
    time_integrate_flux(decomp_CO2, mask_turbulence)

  if with_decomp
    decomp = Dict(
      pairs((;
        t,
        eta, #in log10
        H=decomp_T,
        FC=decomp_CO2,
        LE=decomp_H2O,
        TAU_W=Z, #in log10
        DELTA_TAU=Δτ,
        mask_analysis,
        mask_minima,
        mask_noadvection,
        mask_turbulence,
        mask_σ_t,
      )),
    )
  else
    decomp = nothing
  end

  fluxes = Dict(
    pairs((;
      F_T_nomean=F_T_nomean[analysis_range],
      F_H2O_nomean=F_H2O_nomean[analysis_range],
      F_CO2_nomean=F_CO2_nomean[analysis_range],
      F_T_noadvection=F_T_noadvection[analysis_range],
      F_H2O_noadvection=F_H2O_noadvection[analysis_range],
      F_CO2_noadvection=F_CO2_noadvection[analysis_range],
      F_T_turbulence=F_T_turbulence[analysis_range],
      F_H2O_turbulence=F_H2O_turbulence[analysis_range],
      F_CO2_turbulence=F_CO2_turbulence[analysis_range],
    )),
  )

  time_analysis_h = time_h[time_sampling][analysis_range]

  t1 = time()
  results = Dict(
    pairs((;
      dates=(sdate, edate),
      time_analysis_h,
      time_execution=t1 - t0,
      advec_line,
      fluxes,
      decomp,
      timelags,
    )),
  )
  return results
end

#TODO: Optimise this, freq_tl is usually fixed WaveC can be used for multiple timelag optimisations
function optim_timelag(w, θ, scale_params, freq_tl, τ_max)
  length(w) == length(θ) || throw(error("Signals must be of the same size."))
  work_dim = length(w)
  (freq_peak, σ_waves), wave_params, _ =
    init_wave_conv_kernel(work_dim; scale_params..., with_sigma=true)
  mask_waves = freq_tl[1] .<= freq_peak .<= freq_tl[2]
  WaveC = WaveletConv((work_dim, 1), (work_dim, size(wave_params, 2)))

  τ = 0:(τ_max-1)
  τ_arr = vcat(reverse(-τ .- 1), τ)

  w_ξ = WaveC(w, wave_params)[:, mask_waves]
  θ_ξ = WaveC(θ, wave_params)[:, mask_waves]
  out = irfft(sum(rfft(w_ξ, 1) .* conj(rfft(θ_ξ, 1)), dims=2)[:], work_dim)

  out = vcat(out[end-length(τ)+1:end], out[1:length(τ)])
  return (τ_arr, out)
end

## UTILS, TODO: put in utils.jl
function find_nan_regions(F)
  s = Int64[]
  e = Int64[]
  i = 1
  L = length(F)
  while !isnothing(i)
    i = findnext(isnan, F, i)
    if isnothing(i)
      break
    else
      si = i > 1 ? i - 1 : i
      push!(s, si)
    end

    j = findnext(!isnan, F, i)
    sj = isnothing(j) ? L : j
    push!(e, sj)
    i = j # End loop if j==nothing
  end
  return (s, e)
end
