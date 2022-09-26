
include("LoadPackages.jl")

###################
#### FUNCTIONS ####
###################

#### PARAMETER POSTERIORS ####

function plot_inf_posterior_compare(other_res, from, to, plottitle, legendpos, alphaval)

  p1 = plot(other_res[from:to, :post_inf], title = plottitle, label = "post", alpha = alphaval, lw = 3)
        plot!(other_res[from:to, :post_prime_inf], label = "post_prime", alpha = alphaval, lw = 3, legend = legendpos, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

  return(p1)
end

#### CUMULATIVE ACCEPTANCE RATES ####

function calc_acc_inf(other_res)

  acc_p = cumsum(other_res[:,:is_accepted_inf])./other_res[:,:sample]

  return(acc_p)
end

function plot_acc(acc, plotttitle)
  plot(acc, title = plotttitle, lw = 3, label = "Params", size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)
  println(" ", " ")
  println("Cumulative Acceptance Rate")
  println("Base acc = ", acc[end])

  return current()
end

#### BATCH ACCEPTANCE RATES ####

function calc_acc_batch_inf(other_res, N_its)

    blocks = convert(Int64, N_its/25)

    acc_p = fill(-Inf, blocks)

    for n_tune in 1:blocks

      batch_start = ((n_tune - 1)*25) + 1
      batch_end = (n_tune*25) - 1
      bot = batch_end - batch_start

      acc_p[n_tune] = sum(other_res[batch_start:batch_end,:is_accepted_inf])/bot

    end

    acc_p = repeat(acc_p, inner = 25)

    return(acc_p)
end

function plot_acc_batch_inf(other_res, from, to, N_its, plottitle, legendpos, printtrue)

  acc = calc_acc_batch_inf(other_res, N_its)

  plot(acc[from:to], title = plottitle, label = "Params", lw = 3, legend = legendpos, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)
  if printtrue == true
    println(" ", " ")
    println("Final Batch Acceptance Rate")
    println("Base acc = ", acc[1][end])
  end

  return current()
end

function plot_acc_rates_against_λ_and_m_inf_scaled(other_res, tuning_res, from, to, N_its, plottitle, legendpos, printtrue)

  plot_acc_batch_inf(other_res, from, to, N_its, plottitle, legendpos, printtrue)
    plot!(tuning_res[from:to, :λ_inf], label = "λ_inf", lw = 3)
    plot!([0.33], seriestype = :hline, label = :none, lw = 3)
    plot!( ( tuning_res[from:to,:m_inf]./mean(tuning_res[from:to,:m_inf]) ), label = "m_inf", lw = 3)
end

function plot_acc_rates_against_λ_and_m_inf(other_res, tuning_res, from, to, N_its, plottitle, legendpos, printtrue)

  plot_acc_batch_inf(other_res, from, to, N_its, plottitle, legendpos, printtrue)
    plot!(tuning_res[from:to, :λ_inf], label = "λ_inf", lw = 3)
    plot!([0.33], seriestype = :hline, label = :none, lw = 3)
    plot!(tuning_res[from:to,:m_inf], label = "m_inf", lw = 3)
end

#### REJECTION REASONS ####

function plot_reason(aug_stats, plottitle)

  gdf = combine(groupby(aug_stats, :reason), nrow)
  gdf = sort(gdf, 1)

  labels = ["MH rejected" "Accepted" "Out of Bounds" "Invalid" "No Events to Move"]

  pos = convert.(Int64, gdf[:,1]) .+1

  labels_gdf = labels[:, pos]

  gdf.name = labels_gdf[1, :]
  gdf.percent = gdf[:, 2]/size(aug_stats, 1)*100

  bar(gdf[:, 3],
      gdf[:, 4],
      group = gdf[:, 1],
      labels = labels_gdf,
      title = plottitle)
end

#### PARAMETERS ####

function plot_trace(res, from, to)

  l = @layout [a ; b c]
  p1 = plot(res[from:to,1], title = "β", label = "Trace", lw = 0.5)
          plot!([β_tr], seriestype = :hline, lw = 3, label = "True")
  p2 = plot(res[from:to,2], title = "γ", label = "Trace", lw = 0.5)
          plot!([γ_tr], seriestype = :hline, lw = 3, label = "True")
  p3 = plot(res[from:to,3], title = "δ", label = "Trace", lw = 0.5)
          plot!([δ_tr], seriestype = :hline, lw = 3, label = "True")

  plot(p1, p2, p3, layout = l)
end

function plot_hist(res, from, to, num_bins)

  l = @layout [a ; b c]
  hist1 = fit(Histogram, res[from:to,:β], nbins = num_bins)
  p1 = plot(hist1, title = "β", label = "β")
          plot!([β_tr], seriestype = :vline, lw = 3, label = "True")
  hist2 = fit(Histogram, res[from:to,2], nbins = num_bins)
  p2 = plot(hist2, title = "γ", label = "γ")
          plot!([γ_tr], seriestype = :vline, lw = 3, label = "True")
  hist3 = fit(Histogram, res[from:to,3], nbins = num_bins)
  p3 = plot(hist3, title = "δ", label = "δ")
          plot!([δ_tr], seriestype = :vline, lw = 3, label = "True")

  i1 = argmax(hist1.weights)
  p1max = hist1.edges[1][i1:i1+1]

  i2 = argmax(hist2.weights)
  p2max = hist2.edges[1][i2:i2+1]

  i3 = argmax(hist3.weights)
  p3max = hist3.edges[1][i3:i3+1]

  println("    ")
  println("β Max Posterior Weight = ", p1max)
  println("γ Max Posterior Weight = ", p2max)
  println("δ Max Posterior Weight = ", p3max)

  plot(p1, p2, p3, layout = l)
end


#### HACKED MULTI PLOT ####

@userplot CornerPlot
@recipe function f(cp::CornerPlot)
    m = cp.args[1]

    nl = get(plotattributes, :levels, 10)
    N = size(m, 1)

    labs = pop!(plotattributes, :label, ["x$i" for i=1:N])
    if labs!=[""] && length(labs)!=N
        error("Number of labels not identical to number of datasets")
    end

    legend := false
    layout := (N,N)

    for i in 1:N
        # Do the diagonals
        @series begin
            subplot := i + (i-1)*N
            seriestype := :density
            xlims := (minimum(m[i]), maximum(m[i]))
            ylims := (0, Inf)
            xguide := labs[i]
            x := m[i]
        end
    end

    for i in 1:N
        for j in 1:(i-1)
            # Do the kdeplots
            k = kde((m[j], m[i]))
            dv = vec(k.density)
            inds = reverse(sortperm(dv))
            cd = cumsum(dv[inds])
            C = cd[end]

            levels = []
            for i in 1:nl
                f = i/(nl+1)
                cf = f*C
                ind = searchsortedfirst(cd, cf)
                push!(levels, dv[inds[ind]])
            end
            levels = reverse(levels)

            @series begin
                seriestype := :contour
                subplot := (i-1)*N + j
                seriescolor --> :viridis
                x := k.x
                y := k.x
                z := permutedims(k.density)
                levels := levels
                #xlims := (minimum(m[j]), maximum(m[j]))
                #ylims := (minimum(m[i]), maximum(m[i]))
                xlims := (minimum(m[j]), quantile(m[j], [0.95])[1])
                ylims := (minimum(m[i]), quantile(m[i], [0.95])[1])
                xguide := labs[j]
                yguide := labs[i]
                k.x, k.y, permutedims(k.density)
            end
        end
    end

    for i in 1:N
        for j in (i+1):N
            # Do the scatterplots
            @series begin
                seriestype := scatter
                subplot := (i-1)*N + j
                x := m[j]
                y := m[i]
                markersize --> 0.1
                xlims := (minimum(m[j]), maximum(m[j]))
                ylims := (minimum(m[i]), maximum(m[i]))
                xguide := labs[j]
                yguide := labs[i]
                m[j], m[i]
            end
        end
    end
end
