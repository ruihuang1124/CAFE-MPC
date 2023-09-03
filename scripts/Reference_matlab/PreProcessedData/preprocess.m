clear all;

% Gait ="RunJump/";
% Gait = "MixedHopping/";
% Gait = "RunJump_ICRA23/";
% Gait = "Trot/";
Gait = "Bound/";
% Gait = "Pace/";
% Gait = "FlyPace/";
% Gait = "Pronk/";
% Gait = "Stand/";
% Gait = "FlyTrot/";
body_states = readmatrix(Gait + "body_state.csv");
contacts = readmatrix(Gait + "contact.csv");
foot_placements = readmatrix(Gait + "ee_pos.csv");
foot_velocities = readmatrix(Gait + "ee_vel.csv");
qJs = readmatrix(Gait + "jnt.csv");
% foot_heights = readmatrix(Gait + "height.csv");
t = readmatrix(Gait + "time.csv", "Delimiter",",");