function gait = read_gait_from_file(gait)
    body_states = readmatrix(gait + "body_state.csv");
    contacts = readmatrix(gait + "contact.csv");
    foot_placements = readmatrix(gait + "ee_pos.csv");
    foot_velocities = readmatrix(gait + "ee_vel.csv");
    qJs = readmatrix(gait + "jnt.csv");
    t = readmatrix(gait + "time.csv", "Delimiter",",");
    if isfile(gait + "jntvel.csv")
        qJds = readmatrix(gait + "jntvel.csv");
    else
        qJds = zeros(size(qJs, 1), 12);
    end
    gait = struct("body_states", body_states, ...
                  "contacts", contacts, ...
                  "foot_placements", foot_placements, ...
                  "foot_velocities", foot_velocities, ...
                  "qJs", qJs, ...
                  "qJds", qJds, ...
                  "t", t);    
end
