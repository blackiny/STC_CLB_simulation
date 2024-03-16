import MeshCat as mc

using LinearAlgebra, Plots
using GeometryBasics: HyperRectangle, Vec, Point, Mesh, Rect
using FileIO

include(joinpath(@__DIR__,"./params.jl"))

function skew(v)
    [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end
function dcm_from_phi(ϕ)
    theta = norm(ϕ)
    r = (abs(theta)>1e-12) ? ϕ/theta : zeros(3)
    Q = (I + sin(theta) * skew(r) + (1.0 - cos(theta)) *
    skew(r) * skew(r))
    return Q 
end

function build_CAV_string(vis, num_HV::Int64, initial_pos::Vector{Float64})
    """
    string structure be like:
    HV HV HV ... HV CAV HV
    ---------------
          num_HV
    initial_pos direction: from leading HV to tail HV
    vis: MeshCat visualizer
    """
    @assert length(initial_pos) == num_HV + 2
    @assert num_HV > 0
    @assert maximum(diff(initial_pos)) < 0.0 # initial pos match string structure
    new_car_obj = mc.MeshFileGeometry(joinpath(@__DIR__,"Car.obj"))
    mc.setobject!(vis["Car"]["CAV"], new_car_obj, mc.MeshPhongMaterial(color=mc.RGBA(0.6, 0.6, 1.0, 1.0)))
    mc.setobject!(vis["Car"]["LHV"], new_car_obj, mc.MeshPhongMaterial(color=mc.RGBA(1.0, 0.3, 0.3, 1.0)))
    for i=1:num_HV
        mc.setobject!(vis["Car"]["HV_$i"], new_car_obj, mc.MeshPhongMaterial(color=mc.RGBA(0.3, 0.3, 0.3, 1.0)))
    end
    string_len = initial_pos[1] - initial_pos[end]
    norminal_pos = zeros(length(initial_pos))
    car_width = car_params.w
    # initial middle car at origin
    mid_car_ind = (num_HV + 3) ÷ 2
    for i=1:num_HV
        norminal_pos[i+2] = initial_pos[i+2] - initial_pos[mid_car_ind]
        mc.settransform!(vis["Car"]["HV_$i"],
            mc.Translation([norminal_pos[i+2], -(0.4+0.5*car_width), 0])∘ mc.LinearMap(dcm_from_phi(pi/2*[0,0,1]))∘ mc.LinearMap(dcm_from_phi(pi/2*[1,0,0])))
    end
    norminal_pos[2] = initial_pos[2] - initial_pos[mid_car_ind]
    norminal_pos[1] = initial_pos[1] - initial_pos[mid_car_ind]
    mc.settransform!(vis["Car"]["CAV"],
        mc.Translation([norminal_pos[2], -(0.4+0.5*car_width), 0])∘ mc.LinearMap(dcm_from_phi(pi/2*[0,0,1]))∘ mc.LinearMap(dcm_from_phi(pi/2*[1,0,0])))
    mc.settransform!(vis["Car"]["LHV"],
        mc.Translation([norminal_pos[1], -(0.4+0.5*car_width), 0])∘ mc.LinearMap(dcm_from_phi(pi/2*[0,0,1]))∘ mc.LinearMap(dcm_from_phi(pi/2*[1,0,0])))
    return string_len, norminal_pos
end

function build_road(vis, road_length::Float64, road_width::Float64)
    road_obj = HyperRectangle(Vec(0., 0, 0), Vec(road_length, road_width+0.2, 0.01))
    mc.setobject!(vis["Road"]["base"], road_obj, mc.MeshBasicMaterial(color=mc.RGBA(102/255, 102/255, 102/255, 1.0)))
    mc.settransform!(vis["Road"]["base"], mc.Translation([-road_length/2, -0.5*road_width-0.1, -0.01]))
    road_curb_obj = HyperRectangle(Vec(0., 0, 0), Vec(road_length, 0.05, 0.05))
    mc.setobject!(vis["Road"]["curb_down"], road_curb_obj, mc.MeshBasicMaterial(color=mc.RGBA(255/255, 255/255, 255/255, 1.0)))
    mc.settransform!(vis["Road"]["curb_down"], mc.Translation([-road_length/2, -0.5*road_width-0.05, 0]))
    mc.setobject!(vis["Road"]["curb_up"], road_curb_obj, mc.MeshBasicMaterial(color=mc.RGBA(255/255, 255/255, 255/255, 1.0)))
    mc.settransform!(vis["Road"]["curb_up"], mc.Translation([-road_length/2, 0.5*road_width, 0]))
end

function animate_traffic_flow(num_HV::Int64, X::Vector{Float64}, X_HV::Vector{Vector{Float64}}, ref_X::Vector{Float64}, dt::Float64)
    @assert length(X_HV) == num_HV && num_HV > 0
    @assert length(ref_X) == length(X) && length(ref_X) == length(X_HV[1])
    
    vis = mc.Visualizer()
    mc.setprop!(vis["/Background"], "top_color", mc.RGBA(0.9, 0.9, 0.9, 0.7))
    mc.setprop!(vis["/Background"], "bottom_color", mc.RGBA(0.9, 0.9, 0.9, 0.7))
    # set camera
    mc.settransform!(vis["/Cameras/default"], mc.Translation([0, 0, 0]))
    mc.setprop!(vis["/Cameras/default/rotated/<object>"], "position", [0, 30, 0])
    mc.setvisible!(vis["/Grid"], false)
    mc.setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", 0.85)
    mc.setprop!(vis["/Lights/PointLightPositiveX/<object>"], "intensity", 0.0)

    initial_pos = [ref_X[1];X[1];[X_HV[i][1] for i=1:num_HV]]
    string_len, nominal_pos = build_CAV_string(vis, num_HV, initial_pos)
    #build road
    road_length = ceil(string_len*2.5)
    road_width = lane_params.w
    build_road(vis, road_length, 2.0*road_width)

    car_width = car_params.w
    line_length = lines_params.l
    line_width = lines_params.w
    spacing = lines_params.space
    line_buffer = maximum(diff(ref_X))
    num_lines = floor(Int, 1.1 * (road_length + line_buffer + ref_X[end]) / (line_length + spacing))
    lines_pos = zeros(num_lines)
    lines_od_ind = 0
    lines_bak_pos = 0.0
    lines_bak_ind = 1
    for i=1:num_lines
        lines_pos[i] = (i-1) * (line_length + spacing) - road_length / 2
        if lines_pos[i] ≥ road_length / 2 + line_buffer
            lines_bak_ind = i
            lines_bak_pos = lines_pos[i]
            break
        end
    end
    for i=1:num_lines
        line_obj = HyperRectangle(Vec(0., 0, 0), Vec(line_length, line_width, 0.01))
        mc.setobject!(vis["lines"]["line_$i"], line_obj, mc.MeshBasicMaterial(color=mc.RGBA(250/255, 200/255, 55/255, 1.0)))
        if i < lines_bak_ind
            mc.settransform!(vis["lines"]["line_$i"], mc.Translation([lines_pos[i], -line_width/2, 0]))
        else
            mc.settransform!(vis["lines"]["line_$i"], mc.Translation([lines_bak_pos, -line_width/2, 0]))
        end
    end
    
    anim = mc.Animation(floor(Int,1/dt))
    ref_lines = ref_X .- ref_X[1] 
    for k=1:length(ref_lines)
        mc.atframe(anim, k) do
            #move lines backwards
            lines_od_ind_new = lines_od_ind
            for i=(lines_od_ind+1):(lines_bak_ind-1)
                if lines_pos[i]-ref_lines[k] ≤ -road_length / 2 - line_length
                    lines_od_ind_new = i
                end
                mc.settransform!(vis["lines"]["line_$i"], mc.Translation([lines_pos[i]-ref_lines[k], -line_width/2, 0]))
            end
            lines_od_ind = lines_od_ind_new
            
            tail_pos = lines_pos[lines_bak_ind-1]-ref_lines[k] + (line_length + spacing)
            lines_bak_ind_new = lines_bak_ind
            for j=lines_bak_ind:num_lines
                if lines_bak_pos > tail_pos
                    lines_pos[j] = tail_pos + ref_lines[k]
                    mc.settransform!(vis["lines"]["line_$j"], mc.Translation([lines_pos[j]-ref_lines[k], -line_width/2, 0]))
                    tail_pos += (line_length + spacing)
                    lines_bak_ind_new = j+1
                else
                    mc.settransform!(vis["lines"]["line_$j"], mc.Translation([lines_bak_pos, -line_width/2, 0]))
                end
            end
            lines_bak_ind = lines_bak_ind_new
            
            #move cars relative to LHV
            CAV_pos = nominal_pos[1] + X[k] - ref_X[k]
            mc.settransform!(vis["Car"]["CAV"],
                mc.Translation([CAV_pos, -(0.4+0.5*car_width), 0])∘ mc.LinearMap(dcm_from_phi(pi/2*[0,0,1]))∘ mc.LinearMap(dcm_from_phi(pi/2*[1,0,0])))
            for i=1:num_HV
                HV_pos = nominal_pos[1] + X_HV[i][k] - ref_X[k]
                mc.settransform!(vis["Car"]["HV_$i"],
                    mc.Translation([HV_pos, -(0.4+0.5*car_width), 0])∘ mc.LinearMap(dcm_from_phi(pi/2*[0,0,1]))∘ mc.LinearMap(dcm_from_phi(pi/2*[1,0,0])))
            end
            
        end
    end
    
    mc.setanimation!(vis, anim)
    return mc.render(vis)
end



