function gait_comb = combine_two_gaits(gait1, gait2)
gait_comb.body_states = [gait1.body_states; gait2.body_states];
gait_comb.contacts = [gait1.contacts; gait2.contacts];
gait_comb.foot_placements = [gait1.foot_placements; gait2.foot_placements];
gait_comb.foot_velocities = [gait1.foot_velocities; gait2.foot_velocities];
gait_comb.qJs = [gait1.qJs; gait2.qJs];
dt = gait1.t(2) - gait1.t(1);
gait_comb.t = [gait1.t; gait2.t + gait1.t(end) + dt];
end

