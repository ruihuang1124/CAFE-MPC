function gait_res = flip_left_and_right(gait)
gait_res = gait;
tau_sz = size(gait.body_states, 1);
for k = 1:tau_sz
    gait_res.qJs(k, :) = flip_vector(gait.qJs(k, :));
    gait_res.qJds(k, :) = flip_vector(gait.qJds(k, :));
    gait_res.foot_placements(k, :) = flip_vector(gait.foot_placements(k, :));
    gait_res.foot_velocities(k, :) = flip_vector(gait.foot_velocities(k, :));
    gait_res.GRFs(k, :) = flip_vector(gait.GRFs(k, :));
    gait_res.torques(k, :) = flip_vector(gait.torques(k, :));
    gait_res.contacts(k, :) = flip_vector(gait.contacts(k, :));
    gait_res.status_durations(k, :) = flip_vector(gait.status_durations(k,:));
end
end

function v_res = flip_vector(v_in)
if length(v_in) == 12
    v_res = v_in([4:6, 1:3, 10:12, 7:9]);
end
if length(v_in) == 4
    v_res = v_in([2,1,4,3]);
end
end

