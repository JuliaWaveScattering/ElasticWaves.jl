@recipe function plot(bearing::RollerBearing, t = 0.0)

    grid --> false
    xguide --> "x"
    yguide --> "y"
    aspect_ratio := 1.0
    label --> ""
    linecolor --> :grey

    if bearing.number_of_rollers > 0
        θs = LinRange(0,2pi,bearing.number_of_rollers + 1)[1:end-1]

        θs = θs .+ (bearing.angular_speed * t)

        r = bearing.roller_radius
        R = bearing.rollers_inside ? (bearing.inner_radius - r) : (bearing.outer_radius + r)

        xs = R .* [ [cos(θ),sin(θ)] for θ in θs]

        circles = Circle.(xs,r)
        for circle in circles
            @series begin
                circle
            end
        end
    end

    circle_outer = Circle(bearing.inner_radius)
    @series begin
        circle_outer
    end

    circle_outer = Circle(bearing.outer_radius)
    @series begin
        circle_outer
    end

end