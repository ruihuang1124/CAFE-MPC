addpath(genpath("PreProcessedData/"));
addpath("PostProcessedData/");
%% Lists of gaits
clear
Gaits{1} ="RunJump/";
Gaits{2} = "MixedHopping/";
Gaits{3} = "RunJump_ICRA23/";
Gaits{4} = "Trot/";
Gaits{5} = "Bound/";
Gaits{6} = "Pace/";
Gaits{7} = "FlyPace/";
Gaits{8} = "Pronk/";
Gaits{9} = "Stand/";
Gaits{10} = "FlyTrot/";
Gaits{11} = "Barrel/";
Gaits{12} = "PreBarrelTrot/";
Gaits{13} = "PreBarrelPace/";
Gaits{14} = "PreBarrelPronk/";
gait_prepross_path = "PreProcessedData/";

% First running stage of running barrel roll
gait0_num = 7;
gait0 = read_gait_from_file(gait_prepross_path + Gaits{gait0_num});
gait0 = truncate_gait(gait0, 1, 160);

% Barrel roll stage of running barrel roll 
gait1_num = 11;
gait1 = read_gait_from_file(gait_prepross_path + Gaits{gait1_num});
gait1 = truncate_gait(gait1, 14, 130);

fields = fieldnames(gait1);
for i = 1:length(fields)
    field = fields{i};
    val = getfield(gait1, field);
    val_new = [val(1:67, :); repmat(val(67,:),[15,1]); val(68:end,:)];
    gait1 = setfield(gait1, field, val_new);
end

% Second running state of running barrel roll 
gait2_num = 7;
gait2 = read_gait_from_file(gait_prepross_path + Gaits{gait2_num});

%% Regular locomotion gait
% gait = gait0;

%% Barrel roll -> loco
% Offset the x,y positions in gait 2 to start at the end of gait 1
% gait = gait1;
% gait2.body_states(:,[4,5]) = gait2.body_states(:,[4,5]) + gait.body_states(end,[4,5]);
% gait2.foot_placements(:,[1,4,7,10]) = gait2.foot_placements(:,[1,4,7,10]) ... 
%                                     + gait.foot_placements(end, [1,4,7,10]);
% gait2.foot_placements(:,[2,5,8,11]) = gait2.foot_placements(:,[2,5,8,11]) ... 
%                                     + gait.foot_placements(end, [2,5,8,11]);
% gait2.body_states(:,3) = 2*pi;
% gait = combine_two_gaits(gait, gait2);

%% Loco -> barrel roll -> loco
% Offset the x,y positions in gait 1 to start at the end of gait 0
gait = gait0;
gait1.body_states(:,[4,5]) = gait1.body_states(:,[4,5]) + gait.body_states(end,[4,5]);
gait1.foot_placements(:,[1,4,7,10]) = gait1.foot_placements(:,[1,4,7,10]) ...
                                    + gait.foot_placements(end, [1,4,7,10]);
gait1.foot_placements(:,[2,5,8,11]) = gait1.foot_placements(:,[2,5,8,11]) ...
                                    + gait.foot_placements(end, [2,5,8,11]);
gait = combine_two_gaits(gait, gait1);

gait2.body_states(:,[4,5]) = gait2.body_states(:,[4,5]) + gait.body_states(end,[4,5]);
gait2.foot_placements(:,[1,4,7,10]) = gait2.foot_placements(:,[1,4,7,10]) ... 
                                    + gait.foot_placements(end, [1,4,7,10]);
gait2.foot_placements(:,[2,5,8,11]) = gait2.foot_placements(:,[2,5,8,11]) ... 
                                    + gait.foot_placements(end, [2,5,8,11]);
gait2.body_states(:,3) = 2*pi;
gait = combine_two_gaits(gait, gait2);

%% write to file
write_gait_to_file(gait);

%%
function write_gait_to_file(gait)
    tau_sz = size(gait.body_states, 1);
    gait.status_durations = zeros(tau_sz, 4);
    
    if ~isfield(gait,"torques")
        gait.torques = zeros(tau_sz, 12);
    end
    
    if ~isfield(gait,"GRFs")
        gait.GRFs= zeros(tau_sz, 12);
        % Calculate the GRF reference
        mass = 9; g = 10;
        for k = 1:tau_sz
            F = [0, 0, mass * g / sum(gait.contacts(k, :))];
            for leg = 1: 4
                if gait.contacts(k, leg)
                    gait.GRFs(k, 3*(leg-1)+1:3*leg) = F;
                end
            end
        end
    end
    
    % Calculate contact status durations for each leg
    dt = gait.t(2) - gait.t(1);
    for leg = 1:4
        gait.status_durations(:, leg) = Induce_status_duration_per_leg(gait.contacts(:,leg), dt);
    end
    
    %% Used for HKD
    % gait = flip_left_and_right(gait);
    % for k = 1:tau_sz
    %     eul = gait.body_states(k, 1:3);
    %     euld = gait.body_states(k, 7:9);
    %     gait.body_states(k, 7:9) = euld2omegab(eul', euld');
    % end
    
    %% Write to file
    % write contact information to csv file
    fname = "PostProcessedData/"+"quad_reference.csv";
    fid = fopen(fname, 'w');
    fprintf(fid, 'dt\n');
    fprintf(fid, '%4.3f\n', dt);
    for i = 1:tau_sz
        fprintf(fid, 'body_state \n');
        fprintf_array(fid, gait.body_states(i, :), '%6.3f ');
    
        fprintf(fid, 'jnt_angle\n');
        fprintf_array(fid, gait.qJs(i, :), '%6.3f ');
    
        if isfield(gait, "qJds")
            fprintf(fid, 'jnt_vel\n');
            fprintf_array(fid, gait.qJds(i, :), '%6.3f ');
        end
    
        fprintf(fid, 'foot_placements\n');
        fprintf_array(fid, gait.foot_placements(i, :), '%6.3f ');
    
        fprintf(fid, 'foot_velocities\n');
        fprintf_array(fid, gait.foot_velocities(i, :), '%6.3f ');
    
        fprintf(fid, 'grf\n');
        fprintf_array(fid, gait.GRFs(i, :), '%6.3f ');
    
        fprintf(fid, 'torque\n');
        fprintf_array(fid, gait.torques(i, :), '%6.3f ');
    
        fprintf(fid, 'contact\n');
        fprintf_array(fid, gait.contacts(i, :), '%d ');
    
        fprintf(fid, 'status_dur\n');
        fprintf_array(fid, gait.status_durations(i, :), '%6.3f ');
    end
    fclose(fid);
end

function gait = truncate_gait(gait_in, kstart, kend)
gait = gait_in;
fields = fieldnames(gait_in);
for j=1:numel(fields)
    fieldname = fields{j};
    gait.(fieldname) = gait_in.(fieldname)(kstart:kend, :);
end
end

%% Help functions
function status_durs = Induce_status_duration_per_leg(contacts_leg, dt)
tau_sz = length(contacts_leg);
status_durs = zeros(tau_sz, 1);
status_dur = 0;
status_start = 1;
contact_prev = contacts_leg(1);

for k = 2:tau_sz
    status_dur = status_dur + dt;
    contact_cur = contacts_leg(k);    
    if (contact_cur ~= contact_prev) 
        status_end = k - 1;
        status_durs(status_start:status_end) = status_dur;
        status_start = k;
        status_dur = 0;
        contact_prev = contact_cur;
    end
    
    if k == tau_sz
        status_end = k;
        status_durs(status_start:status_end) = status_dur;
    end        
end
end

function fprintf_array(fid, a, format)
for j = 1:length(a)
    fprintf(fid, format, a(j));
end
fprintf(fid, '\n');
end
