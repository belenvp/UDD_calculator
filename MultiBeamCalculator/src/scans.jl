function align_profiles(sim_results::Vector)
    xs = Float64[]
    ys = Float64[]
    x_axes = Vector{Float64}[]
    ref_i = ceil(Int, length(sim_results) / 2)

    for (i, result) in enumerate(sim_results)
        y_profile = result.y_profile

        # Find the approximate middle position by finding
        # the midpoint between the two edges of the signal.
        threshold = mean(y_profile[1:10])
        start_idx = findfirst(>(threshold), y_profile)
        end_idx = length(y_profile) - findfirst(>(threshold), reverse(y_profile))
        middle_idx = round(Int, (end_idx - start_idx) / 2)

        # Commented out for now because the cross-correlation alignment seems fishy
        # if i != ref_i
        #     # Divide by 2 because we're taking the middle 50% of the signal, and then by 2 again
        #     # because this variable is half the window length.
        #     central_window_half_length = round(Int, (end_idx - start_idx) / 2 / 6)
        #     # @show central_window_half_length
        #     center_idxs = max(1, middle_idx - central_window_half_length):min(length(y_profile), middle_idx + central_window_half_length)
        #     reference_profile_center = @view sim_results[ref_i].y_profile[center_idxs]
        #     profile_center = @view y_profile[center_idxs]
        #     lags = collect(-(length(profile_center)-1):(length(profile_center) - 1))
        #     correlation = DSP.xcorr(reference_profile_center, profile_center)

        #     middle_idx = first(center_idxs) + lags[argmax(correlation)]
        # end

        x_axis = Vector{Float64}(1:length(y_profile))
        x_axis .-= middle_idx
        x_axis .+= 240 # -minimum(x_axis) + 1
        x_axis ./= 2
        push!(x_axes, x_axis)
        append!(xs, x_axis)
        append!(ys, y_profile)
    end

    y_signal_counts = fit(Histogram, xs, uweights(length(xs)), -99:900)
    y_signal_hist = fit(Histogram, xs, Weights(ys), -99:900)
    mean_signal = y_signal_hist.weights ./ y_signal_counts.weights
    filter!(isfinite, mean_signal)

    return (; x_axes, mean_signal)
end
