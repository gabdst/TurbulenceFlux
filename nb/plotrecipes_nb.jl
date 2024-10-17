using CairoMakie
using DelaunayTriangulation
using Colors
using ColorSchemes
using Measurements

using Unitful
using LaTeXStrings
using Printf
using FFTW
using KernelDensity
using StatsBase
using Statistics
using LinearAlgebra

cmap_flux = Makie.Reverse(:bam)
cmap_tau = Makie.Reverse(:roma)
_center_edges(e) = (e[2:end] .+ e[1:(end-1)]) / 2

cmap_Δhist = :cork
cmap_hist = Makie.Reverse(:lapaz)

mask_color = :black
mask_alpha = 0.25
_colors_ = Makie.wong_colors()
color_flux = _colors_[1]
color_ref = :black

label_flux_turb = L"\text{HR\!-\!TM}"
label_flux_noturb = L"\text{HR}"
label_flux_thresh = L"\text{HR\!-\!thresh}"
label_eddycov30_nou = L"\text{EC30}"

set_theme!(theme_latexfonts())
update_theme!(fontsize=12, size=(15, 7) .* 72)

nounits = (L"\left[\text{-}\right]",)
getunits(s::Symbol) = begin
  s in [:SensibleHeat, :LatentHeat, :SW_IN] && return (
    L"\left[\mathrm{W}\, \mathrm{m}^{-2}\right]",
  )
  s == :FC && return (
    L"\left[\mathrm{{\mu}mol}\, \mathrm{m}^{-2}\, \mathrm{s}^{-1}\right]",
  )
  s == :Ta && return (
    L"\left[\mathrm{K}\right]",)
  s == :H2O && return (
    L"\left[\mathrm{mmol}\, \mathrm{mol}^{-1}\right]",
  )
  s == :CO2 && return (
    L"\left[\mathrm{{\mu}mol}\, \mathrm{mol}^{-1}\right]",
  )
  s == :time && return (L"\left[\mathrm{h}\right]",)
  s == :freq && return (L"\left[\mathrm{Hz}\right]",)
  s == :tau && return (L"\left[\mathrm{m}^2\, \mathrm{s}^{-2}\right]",)
  s in [:w, :u, :v, :ustar] && return (L"\left[\mathrm{m}\, \mathrm{s}^{-1}\right]",)
  s in [:proba, :nfreq, :deltatau] && return nounits
  throw(error("Wrong"))
end
getname(s::Symbol) = begin
  s == :SensibleHeat && return (L"F_H",)
  s == :LatentHeat && return (L"F_\mathrm{LE}",)
  s == :Ta && return (L"T_a",)
  s == :FC && return (L"F_\mathrm{C}",)
  s == :CO2 && return (L"\mathrm{CO_2}",)
  s == :H2O && return (L"\mathrm{H_2O}",)
  s == :SW_IN && return (L"\mathrm{SW_{IN}}",)
  s == :w && return (L"w",)
  s == :time && return (L"\mathrm{time}",)
  s == :ustar && return (L"u^*",)
  s == :nfreq && return (L"\eta",)
  s == :freq && return (L"\mathrm{frequency}",)
  s == :tau && return (L"\tau_w(t,\,\eta)",)
  s == :deltatau && return (L"\Delta \log \tau_w(t,\,\eta)",)
  s == :proba && return (L"\mathrm{Probability}",)
  s == :deltaproba && return (L"\mathrm{\Delta Probability}",)
  throw(error("Wrong"))
end
getlabel(s::Symbol; with_name=true, with_units=true) = begin
  names = getname(s)
  units = getunits(s)
  if with_name && with_units
    out = ([L"%$n %$u" for (n, u) in zip(names, units)]...,)
  elseif with_name
    out = ([L"%$n" for n in names]...,)
  elseif with_units
    out = ([L"%$u" for u in units]...,)
  else
    out = ("",)
  end
  return length(out) == 1 ? out[1] : out
end

log10_format = values -> [L"10^{%$(Int(value))}" for value in values]
to_h_format(x) = begin
  x = x * 1u"hr"
  h = floor(u"hr", x)
  h = Int(h.val)
  L"%$(h)\mathrm{h}"
end
time_format = values -> map(to_h_format, values)
hour_format = time_format
h_format = time_format
to_hmin_format(x) = begin
  x = x * 1u"hr"
  h = floor(u"hr", x)
  min = floor(u"minute", x - h)
  h = Int(h.val)
  min = Int(min.val)
  L"%$(h)\mathrm{h}%$(min)\mathrm{min}"
end
hmin_format = values -> map(to_hmin_format, values)
to_s_format(x) = begin
  x = x * 1u"s"
  h = floor(u"s", x)
  h = Int(h.val)
  L"%$(h)\mathrm{s}"
end
s_format = values -> map(to_s_format, values)

function add_labels!(gls, padding=(0, 5, 5, 0))
  L = length(gls)
  labs = range('A', length=L, step=1)
  for (l, gl) in zip(labs, gls)
    Label(gl[1, 1, TopLeft()], string(l), halign=:left, font=:bold, padding=padding)
  end
end

function plot_speeds(Signals; xticks=Makie.automatic, xtickformat=Makie.automatic)
  f = Figure()
  ax = Axis(f[1, 1]; xticks, xtickformat, ylabel=getlabel(:w,with_name=false))
  t = time_h[sampling_1min]
  lines!(ax, t, view(Signals.U, sampling_1min), label="U")
  lines!(ax, t, view(Signals.V, sampling_1min), label="V")
  lines!(ax, t, view(Signals.W, sampling_1min), label="W")
  axislegend()
  return f
end

function plot_scalars(Signals; xticks=Makie.automatic, xtickformat=Makie.automatic)

  t = time_h[sampling_1min]
  f = Figure()
  color = :green
  ax1 = Axis(
    f[1, 1];
    ylabel=getlabel(:CO2; with_name=false),
    ylabelcolor=color,
    yticklabelcolor=color,
  )
  lines!(ax1, t, Signals.CO2[sampling_1min], color=color)
  color = :blue
  ax2 = Axis(
    f[1, 1];
    ylabel=getlabel(:H2O; with_name=false),
    ylabelcolor=color,
    yticklabelcolor=:blue,
    yaxisposition=:right,
  )
  lines!(ax2, t, Signals.H2O[sampling_1min], color=color)
  color = :red
  ax3 = Axis(
    f[1, 1];
    ylabel=getlabel(:Ta; with_name=false),
    ylabelcolor=color,
    yticklabelcolor=:red,
    yaxisposition=:right,
    ylabelpadding=15,
    xticks,
    xtickformat,
  )
  lines!(ax3, t, Signals.T[sampling_1min] .+ 273.15, color=color)
  hidexdecorations!(ax2)
  hidexdecorations!(ax1)
  return f
end

function plot_wave(
  σ_t,
  wave;
  f=nothing,
  factor=3,
  colors=Iterators.cycle(
    distinguishable_colors(5, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true),
  ),
)
  s_t = ceil(Int, σ_t * factor)
  x = vcat(wave[(end-s_t+2):end], wave[1:s_t])
  x = x / maximum(x)
  t = (((-s_t+2):s_t) .- 1) / (fs * 60 * 60)
  lab = Printf.@sprintf("σ_t: %4.2e h", σ_t / (fs * 60 * 60))
  s_f = σ_t / (fs * 60 * 60)
  if isnothing(f)
    f = Figure(fontsize=18)
    color = iterate(colors, 1)[1]
    position = :bottom
    xlabel = getlabel(:time)
    hide_y = false
  else
    color = iterate(colors, length(f.content) + 1)[1]
    position = :top
    xlabel = ""
    hide_y = true
  end

  ax_out = Axis(
    f[1, 1],
    xlabel=xlabel,
    ylabel=latexstring("\\frac{\\psi_\\xi(t)}{\\psi_\\xi(0)}"),
    ylabelrotation=0,
    ylabelsize=18.0,
    xaxisposition=position,
    xtickcolor=color,
    xticklabelcolor=color,
  )
  ax_in = Axis(
    f[1, 1],
    xaxisposition=position,
    xtickcolor=color,
    xticklabelcolor=color,
    xticks=([-s_f, s_f], [L"-\sigma_t", L"\sigma_t"]),
    xtickalign=1,
  )
  linkxaxes!(ax_in, ax_out)
  linkyaxes!(ax_in, ax_out)

  if hide_y
    hideydecorations!(ax_in)
    hideydecorations!(ax_out)
  end
  lines!(t, x, color=color, label=lab)
  return f
end

function plot_waves_fft(waves)
  waves_fft = rfft(waves, 1)
  freq_sampling = 2:100:size(waves_fft, 1)
  freqs = range(0, fs / 2, length=size(waves_fft, 1))[freq_sampling]
  f = Figure(fontsize=18)
  ax = Axis(
    f[1, 1],
    xscale=log10,
    xlabel=getlabel(:freq),
    ylabel=latexstring("\\hat{\\psi}_\\xi(f)"),
    ylabelrotation=0,
  )
  for i in axes(waves_fft, 2)
    wave = abs2.(view(waves_fft, freq_sampling, i))
    lines!(ax, freqs, wave[:])
  end
  return f
end

function plot_sum_waves_fft(waves)
  waves_fft = rfft(waves, 1)
  combined_spec = sum(abs2, waves_fft, dims=2)
  freq_sampling = 2:5:size(waves_fft, 1)
  freqs = range(0, fs / 2, length=size(waves_fft, 1))[freq_sampling]
  f = Figure(fontsize=18)
  ax = Axis(
    f[1, 1],
    xscale=log10,
    xlabel=getlabel(:freq),
    ylabel=latexstring("\\sum_\\xi\\left|\\hat{\\psi}_\\xi(f)\\right|^2"),
    ylabelrotation=0,
    limits=(nothing, nothing, 0.75, 1.25),
  )
  lines!(ax, freqs, combined_spec[freq_sampling])
  return f
end

function plot_averaging_kernel(KernelC, kernel_params, factor=3)
  work_dim = KernelC.input_dim[1]
  u = div(work_dim, 2)
  δ_dirac = circshift(vcat(1, zeros(work_dim - 1)), u) # need to add a phase since the output is cropped on the border
  kernel = circshift(dropdims(KernelC(δ_dirac, kernel_params), dims=2), -u)
  σ_t = KernelC.kernel_dim[1] * kernel_params[1]
  s_t = ceil(Int, σ_t * factor)
  x = vcat(kernel[(end-s_t+2):end], kernel[1:s_t])
  x /= maximum(x)
  t = (((-s_t+2):s_t) .- 1) / (fs * 60 * 60)
  f = Figure(fontsize=18)
  ax = Axis(
    f[1, 1],
    xlabel=getlabel(:time),
    ylabel=L"\frac{\phi(t)}{\phi(0)}",
    ylabelrotation=0,
  )
  s_f = σ_t / (fs * 60 * 60)
  ax_in = Axis(f[1, 1], xticks=([-s_f, s_f], [L"-\sigma_t", L"\sigma_t"]))
  linkxaxes!(ax_in, ax)
  linkyaxes!(ax_in, ax)
  lines!(ax, t, x)
  return f
end

function add_colorbar!(
  figpos,
  co;
  vmin,
  vmax,
  label_z=label_momentum,
  ticks=LogTicks(round(Int, vmin):round(Int, vmax)),
  tickformat=log10_format,
  vertical=false,
)
  Colorbar(
    figpos,
    co,
    label=label_z,
    ticks=ticks,
    tickformat=tickformat,
    vertical=vertical,
    size=8,
  )
  return nothing
end

function add_mask!(
  ax_co,
  t,
  eta,
  mask;
  mask_color=:black,
  mask_alpha=0.25,
  tri=nothing,
)
  if !isnothing(tri)
    tricontourf!(ax_co, tri, mask, colormap=[(:white, 0.0), (mask_color, mask_alpha)])
  else
    tricontourf!(
      ax_co,
      t,
      eta,
      mask,
      colormap=[(:white, 0.0), (mask_color, mask_alpha)],
    )
  end
  return nothing
end

function plot_time_eddyscale!(
  fig,
  t,
  eta,
  Z;
  vmin,
  vmax,
  ylabel=getlabel(:nfreq),
  xlabel=getlabel(:time),
  label_z="",
  title="",
  mask=nothing,
  colormap=cmap_flux,
  nlevels=20,
  levels=range(vmin, vmax, length=nlevels),
  ticks=Makie.automatic,
  tickformat=Makie.automatic,
  xticks=Makie.automatic,
  xtickformat=Makie.automatic,
  yticks=LogTicks(-3:2),
  ytickformat=log10_format,
  vertical=false,
  tri=nothing,
  mask_alpha=0.25,
)

  ax_decomp =
    Axis(fig[1, 1]; xlabel, ylabel, title, yticks, ytickformat, xticks, xtickformat)
  if !isnothing(tri)
    co = tricontourf!(
      ax_decomp,
      tri,
      Z;
      levels,
      colormap,
      extendhigh=:auto,
      extendlow=:auto,
    )
  else
    co = tricontourf!(
      ax_decomp,
      t,
      eta,
      Z;
      levels,
      colormap,
      extendhigh=:auto,
      extendlow=:auto,
    )
  end
  fig_cb_pos = vertical ? fig[1, 2] : fig[0, 1]
  add_colorbar!(
    fig_cb_pos,
    co;
    vmin,
    vmax,
    label_z=label_z,
    ticks,
    tickformat,
    vertical=vertical,
  )
  if !isnothing(mask)
    add_mask!(ax_decomp, t, eta, mask; tri, mask_alpha)
  end
  return ax_decomp
end

function make_constrained_triangulation(x, y, S)
  LI = LinearIndices(S)
  pts = [x'; y']
  boundary_nodes = reverse(vcat(LI[:, 1], LI[end, :], LI[end:-1:1, end], LI[1, end:-1:1]))
  tri = triangulate(pts; boundary_nodes)
  return tri
end
