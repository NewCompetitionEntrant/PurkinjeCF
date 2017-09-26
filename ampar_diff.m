function Synapse=ampar_diff(sim_time, impaired)

% Rate constants
K_exo=(0.1/60);                     % per sec
K_endo=0.02;                        % per sec
K_on=1;                             % per sec
K_off=0.02;                         % per sec
D_psd=0.0002;                       % um^2 per sec
D_in=0.03;                          % um^2 per sec
D_out=0.05;                         % um^2 per sec

% Neuron parameters
NeuronParams.soma_dia=22;           % um
NeuronParams.psd_width=0.25;        % um
NeuronParams.synapse_width=0.5;     % um
NeuronParams.AMPARs_per_vesicle=50;

% Set the location of the synapse
radius=NeuronParams.soma_dia/2;
s=NeuronParams.synapse_width;
NeuronParams.synapse_loc=[-s/2, radius-s, s/2, radius];

% Set the location of the PSD
s=NeuronParams.psd_width;
NeuronParams.psd_loc=[-s/2, radius-s, s/2, radius];


% Simulation parameters
dt=1; % sec
N_syn_tot=4;

% Create the synapses with different rates
Synapse(1).K_exo=K_exo;
Synapse(1).K_endo=K_endo;
Synapse(1).K_on=K_on;
Synapse(1).K_off=K_off;
Synapse(1).D_psd=D_psd;
Synapse(1).D_in=D_in;
Synapse(1).D_out=D_out;
Synapse(1).AMPAR={};
Synapse(1).N_AMPAR_vesc=0;

if (impaired == 0)
	Synapse(2).K_exo=K_exo*2;
	Synapse(2).K_endo=K_endo;
	Synapse(2).K_on=K_on*2;
	Synapse(2).K_off=K_off;
	Synapse(2).D_psd=D_psd;
	Synapse(2).D_in=D_in*2;
	Synapse(2).D_out=D_out*2;
	Synapse(2).AMPAR={};
	Synapse(2).N_AMPAR_vesc=0;
	
	Synapse(3).K_exo=K_exo*5;
	Synapse(3).K_endo=K_endo;
	Synapse(3).K_on=K_on*5;
	Synapse(3).K_off=K_off;
	Synapse(3).D_psd=D_psd;
	Synapse(3).D_in=D_in*5;
	Synapse(3).D_out=D_out*5;
	Synapse(3).AMPAR={};
	Synapse(3).N_AMPAR_vesc=0;
	
	Synapse(4).K_exo=K_exo*10;
	Synapse(4).K_endo=K_endo;
	Synapse(4).K_on=K_on*10;
	Synapse(4).K_off=K_off;
	Synapse(4).D_psd=D_psd;
	Synapse(4).D_in=D_in*10;
	Synapse(4).D_out=D_out*10;
	Synapse(4).AMPAR={};
	Synapse(4).N_AMPAR_vesc=0;
else
	Synapse(2).K_exo=K_exo;
	Synapse(2).K_endo=K_endo;
	Synapse(2).K_on=K_on;
	Synapse(2).K_off=K_off;
	Synapse(2).D_psd=D_psd;
	Synapse(2).D_in=D_in;
	Synapse(2).D_out=D_out;
	Synapse(2).AMPAR={};
	Synapse(2).N_AMPAR_vesc=0;
	
	Synapse(3).K_exo=K_exo;
	Synapse(3).K_endo=K_endo;
	Synapse(3).K_on=K_on;
	Synapse(3).K_off=K_off;
	Synapse(3).D_psd=D_psd;
	Synapse(3).D_in=D_in;
	Synapse(3).D_out=D_out;
	Synapse(3).AMPAR={};
	Synapse(3).N_AMPAR_vesc=0;
	
	Synapse(4).K_exo=K_exo;
	Synapse(4).K_endo=K_endo;
	Synapse(4).K_on=K_on;
	Synapse(4).K_off=K_off;
	Synapse(4).D_psd=D_psd;
	Synapse(4).D_in=D_in;
	Synapse(4).D_out=D_out;
	Synapse(4).AMPAR={};
	Synapse(4).N_AMPAR_vesc=0;
end


% Set up the exocytosis and endocytosis timers for each synapse
for (N_syn=1:N_syn_tot)
	if (Synapse(N_syn).K_exo < 1.0)
		init_exo_timer=round(1/Synapse(N_syn).K_exo);
	else
		init_exo_timer=1;
	end

	if (Synapse(N_syn).K_endo < 1.0)
		init_endo_timer=round(1/Synapse(N_syn).K_endo);
	else
		init_endo_timer=1;
	end

	Synapse(N_syn).exo_timer=init_exo_timer;
	Synapse(N_syn).init_exo_timer=init_exo_timer;
	Synapse(N_syn).endo_timer=init_endo_timer;
	Synapse(N_syn).init_endo_timer=init_endo_timer;
end


% Main simulation loop
h=waitbar(0, 'Running the AMPAR diffusion simulation...');
for (sim_time_completed=0:dt:sim_time)
	if (mod(sim_time_completed, 10) == 0)
		waitbar(sim_time_completed/sim_time);
	end

	for (N_syn=1:N_syn_tot)
		if (Synapse(N_syn).N_AMPAR_vesc > 0)
			Synapse(N_syn) = DiffuseAMPARs(Synapse(N_syn), NeuronParams);
		end

		Synapse(N_syn).exo_timer=Synapse(N_syn).exo_timer-1;
		if (Synapse(N_syn).exo_timer == 0)
			Synapse(N_syn) = AMPAR_exocytosis(Synapse(N_syn), NeuronParams);
			Synapse(N_syn).exo_timer=Synapse(N_syn).init_exo_timer;
		end

		Synapse(N_syn).endo_timer=Synapse(N_syn).endo_timer-1;
		if (Synapse(N_syn).endo_timer == 0)
			Synapse(N_syn) = AMPAR_endocytosis(Synapse(N_syn), NeuronParams);
			Synapse(N_syn).endo_timer=Synapse(N_syn).init_endo_timer;
		end
	end
end

close(h);

[X, Y]=meshgrid(-radius:0.5:radius, -radius:0.5:radius);
Z=zeros(size(X, 1), size(Y, 1));
Z_bound=Z;
grid_res_X=size(X, 1)/NeuronParams.soma_dia;
grid_res_Y=size(Y, 1)/NeuronParams.soma_dia;
grid_center_X=size(Z, 1)/2;
grid_center_Y=size(Z, 2)/2;

for (N_syn=1:N_syn_tot)
	for (i=1:Synapse(N_syn).N_AMPAR_vesc)
		switch (N_syn)
			case 1
				% No rotation
				AMPAR_x=Synapse(N_syn).AMPAR{i}.x;
				AMPAR_y=Synapse(N_syn).AMPAR{i}.y;
			case 2
				% rotate the co-ordinates by 90 degrees clockwise
				AMPAR_x=Synapse(N_syn).AMPAR{i}.y;
				AMPAR_y=-1*Synapse(N_syn).AMPAR{i}.x;
			case 3
				% rotate the co-ordinates by 180 degrees clockwise
				AMPAR_x=-1*Synapse(N_syn).AMPAR{i}.x;
				AMPAR_y=-1*Synapse(N_syn).AMPAR{i}.y;
			case 4
				% rotate the co-ordinates by 270 degrees clockwise
				AMPAR_x=-1*Synapse(N_syn).AMPAR{i}.y;
				AMPAR_y=Synapse(N_syn).AMPAR{i}.x;
			otherwise
				% impossible to come here
		end

		grid_x=round((AMPAR_x * grid_res_X) + grid_center_X);
		grid_y=round((AMPAR_y * grid_res_Y) + grid_center_Y); 

		if (grid_x < 1)
			grid_x = 1;
		elseif (grid_x > size(Z, 1))
			grid_x = size(Z, 1);
		end

		if (grid_y < 1)
			grid_y = 1;
		elseif (grid_y > size(Z, 1))
			grid_y = size(Z, 1);
		end

		Z(grid_x, grid_y)=Z(grid_x, grid_y)+Synapse(N_syn).AMPAR{i}.N_AMPARs;
		Z_bound(grid_x, grid_y)=Z_bound(grid_x, grid_y)+sum(Synapse(N_syn).AMPAR{i}.bound_vec);
	end
end

figure(1);
surf(X, Y, Z);
shading interp;
colormap(jet(20));
xlabel('X (\muM)');
ylabel('Y (\muM)');
colorbar;
title(sprintf('Spatial distribution of AMPARs for 4 Purkinje cell-Climbing fiber synapses after %d seconds', sim_time))

figure(2);
surf(X, Y, Z_bound);
shading interp;
colormap(jet(20));
xlabel('X (\muM)');
ylabel('Y (\muM)');
colorbar;
title(sprintf('Bounded AMPARs for 4 Purkinje cell-Climbing fiber synapses after %d seconds', sim_time))


% Release an AMPAR vesicle from the Endoplasmic Reticulum
function synapse_ret=AMPAR_exocytosis(synapse, NeuronParams)

% Generate a random AMPAR exocytosis location
AMPAR_vesc.x=2*randn();
if (AMPAR_vesc.x > 4)
	AMPAR_vesc.x = 4;
elseif (AMPAR_vesc.x < -4)
	AMPAR_vesc.x = -4;
end

if (AMPAR_vesc.x >= 0)
	AMPAR_vesc.y=AMPAR_vesc.x;
else
	AMPAR_vesc.y=-AMPAR_vesc.x;
end

AMPAR_vesc.N_AMPARs=NeuronParams.AMPARs_per_vesicle;
AMPAR_vesc.y_cntr=0;
AMPAR_vesc.in_synapse=0;
AMPAR_vesc.in_PSD=0;
AMPAR_vesc.bound_time=zeros(NeuronParams.AMPARs_per_vesicle, 1);
AMPAR_vesc.unbound_time=zeros(NeuronParams.AMPARs_per_vesicle, 1);
AMPAR_vesc.bound_vec=zeros(NeuronParams.AMPARs_per_vesicle, 1);
AMPAR_vesc.available=ones(NeuronParams.AMPARs_per_vesicle, 1);

synapse_ret=synapse;
synapse_ret.N_AMPAR_vesc=synapse.N_AMPAR_vesc+1;
synapse_ret.AMPAR{synapse_ret.N_AMPAR_vesc}=AMPAR_vesc;


% Absorb an AMPAR into the Endoplasmic Reticulum
function synapse_ret=AMPAR_endocytosis(synapse, NeuronParams)

synapse_ret=synapse;

for (i=1:synapse_ret.N_AMPAR_vesc)
	if ((synapse_ret.AMPAR{i}.in_synapse == 0) && (synapse_ret.AMPAR{i}.in_PSD == 0))
		if (synapse_ret.AMPAR{i}.N_AMPARs > 0)
			% Select the first available AMPAR in this vesicle for endocytosis
			idx = find(synapse_ret.AMPAR{i}.available == 1);
			synapse_ret.AMPAR{i}.available(idx(1))=0;

			% Decrement the number of AMPARs in this vesicle
			synapse_ret.AMPAR{i}.N_AMPARs=synapse_ret.AMPAR{i}.N_AMPARs-1;

			% disp('Endocytosis of AMPAR');
			break;
		end
	end
end

% Diffuse AMPAR vesicles towards the synapse
function synapse_ret=DiffuseAMPARs(synapse, NeuronParams)

synapse_ret=synapse;

bind_count=0;
unbind_count=0;

% Determine the number of AMPARs to bind to and unbind from the scaffolding proteins in the PSD
if (synapse.K_on < 1)
	bind_time_thresh=(1/synapse.K_on);
	tot_bind_count = 1;
else
	bind_time_thresh=1;
	tot_bind_count=synapse.K_on;
end

if (synapse.K_off < 1)
	unbind_time_thresh=(1/synapse.K_off);
	tot_unbind_count = 1;
else
	unbind_time_thresh=1;
	tot_unbind_count = synapse.K_off;
end

for (i=1:synapse_ret.N_AMPAR_vesc)
	% Skip empty vesicles
	if (synapse_ret.AMPAR{i}.N_AMPARs == 0)
		continue;
	end

	if (synapse_ret.AMPAR{i}.in_PSD == 1)
		% AMPARs in the PSD

		% Vesicle will move in the PSD if no AMPAR is bounded
		if (sum(synapse_ret.AMPAR{i}.bound_vec) == 0)
			curr_x=synapse_ret.AMPAR{i}.x;
			curr_y=synapse_ret.AMPAR{i}.y;
			y_cntr=synapse_ret.AMPAR{i}.y_cntr;
			[new_x, new_y, new_y_cntr]=move_AMPAR(synapse_ret.D_psd, curr_x, curr_y, y_cntr, NeuronParams);

			synapse_ret.AMPAR{i}.x=new_x;
			synapse_ret.AMPAR{i}.y=new_y;
			synapse_ret.AMPAR{i}.y_cntr=new_y_cntr;

			% Check if the AMPAR is still in the PSD	
			if ((new_x >= NeuronParams.psd_loc(1)) && ...
		    	(new_y >= NeuronParams.psd_loc(2)) && ...
		    	(new_x <= NeuronParams.psd_loc(3)) && ...
		    	(new_y <= NeuronParams.psd_loc(4)))
				synapse_ret.AMPAR{i}.in_PSD = 1; % redundant
			else
				synapse_ret.AMPAR{i}.in_PSD = 0;

				% Check if the AMPAR has moved into the synapse	
				if ((new_x >= NeuronParams.synapse_loc(1)) && ...
		    		(new_y >= NeuronParams.synapse_loc(2)) && ...
		    		(new_x <= NeuronParams.synapse_loc(3)) && ...
		    		(new_y <= NeuronParams.synapse_loc(4)))
					synapse_ret.AMPAR{i}.in_synapse = 1;
				else
					synapse_ret.AMPAR{i}.in_synapse = 0;
				end
			end
		end

		AMPAR_bound_times=synapse_ret.AMPAR{i}.bound_time;
		AMPAR_unbound_times=synapse_ret.AMPAR{i}.unbound_time;
		bounded_AMPARs=synapse_ret.AMPAR{i}.bound_vec;
		available_AMPARs=synapse_ret.AMPAR{i}.available;

		% Find all AMPAR indices in this vesicle that have not undergone endocytosis
		avail_idx = find(available_AMPARs == 1);

		% Update bound and unbound times of the AMPARs in the PSD
		[AMPAR_bound_times, AMPAR_unbound_times]=update_bound_unbound_times(AMPAR_bound_times, AMPAR_unbound_times, bounded_AMPARs, avail_idx);

		% Bind the AMPARs that have waited long enough in the PSD
		if ((tot_bind_count-bind_count) > 0)
			[AMPAR_bound_times, AMPAR_unbound_times, bounded_AMPARs, N_bounded]=...
				bind_AMPARs(AMPAR_bound_times, AMPAR_unbound_times, bounded_AMPARs, avail_idx, bind_time_thresh, tot_bind_count-bind_count);
			bind_count=bind_count+N_bounded;

			synapse_ret.AMPAR{i}.bound_time=AMPAR_bound_times;
			synapse_ret.AMPAR{i}.unbound_time=AMPAR_unbound_times;
			synapse_ret.AMPAR{i}.bound_vec=bounded_AMPARs;
		end

		% Unbind the AMPARs that have stayed bounded long enough in the PSD
		if ((tot_unbind_count-unbind_count) > 0)
			[AMPAR_bound_times, AMPAR_unbound_times, bounded_AMPARs, N_unbounded]=...
				unbind_AMPARs(AMPAR_bound_times, AMPAR_unbound_times, bounded_AMPARs, avail_idx, unbind_time_thresh, tot_unbind_count-unbind_count);
			unbind_count=unbind_count+N_unbounded;

			synapse_ret.AMPAR{i}.bound_time=AMPAR_bound_times;
			synapse_ret.AMPAR{i}.unbound_time=AMPAR_unbound_times;
			synapse_ret.AMPAR{i}.bound_vec=bounded_AMPARs;
		end
	elseif (synapse_ret.AMPAR{i}.in_synapse == 1) 
		% AMPARs in the synapse
		curr_x=synapse_ret.AMPAR{i}.x;
		curr_y=synapse_ret.AMPAR{i}.y;
		y_cntr=synapse_ret.AMPAR{i}.y_cntr;
		[new_x, new_y, new_y_cntr]=move_AMPAR(synapse_ret.D_in, curr_x, curr_y, y_cntr, NeuronParams);

		synapse_ret.AMPAR{i}.x=new_x;
		synapse_ret.AMPAR{i}.y=new_y;
		synapse_ret.AMPAR{i}.y_cntr=new_y_cntr;

		% Check if the AMPAR vesicle has moved into the PSD	
		if ((new_x >= NeuronParams.psd_loc(1)) && ...
		    (new_y >= NeuronParams.psd_loc(2)) && ...
		    (new_x <= NeuronParams.psd_loc(3)) && ...
		    (new_y <= NeuronParams.psd_loc(4)))
			synapse_ret.AMPAR{i}.in_PSD = 1;
			synapse_ret.AMPAR{i}.in_synapse = 0;
		elseif ((new_x >= NeuronParams.synapse_loc(1)) && ...
		       (new_y >= NeuronParams.synapse_loc(2)) && ...
		       (new_x <= NeuronParams.synapse_loc(3)) && ...
		       (new_y <= NeuronParams.synapse_loc(4)))
			synapse_ret.AMPAR{i}.in_synapse = 1; % redundant
			synapse_ret.AMPAR{i}.in_PSD = 0;
	  	else
			% AMPAR vescicle has moved out of the synapse
			synapse_ret.AMPAR{i}.in_synapse = 0; 
			synapse_ret.AMPAR{i}.in_PSD = 0;
		end
	else
		% AMPARs out of the synapse
		curr_x=synapse_ret.AMPAR{i}.x;
		curr_y=synapse_ret.AMPAR{i}.y;
		y_cntr=synapse_ret.AMPAR{i}.y_cntr;
		[new_x, new_y, new_y_cntr]=move_AMPAR(synapse_ret.D_out, curr_x, curr_y, y_cntr, NeuronParams);

		synapse_ret.AMPAR{i}.x=new_x;
		synapse_ret.AMPAR{i}.y=new_y;
		synapse_ret.AMPAR{i}.y_cntr=new_y_cntr;

		% Check if the AMPAR has moved into the synapse	
		if ((new_x >= NeuronParams.synapse_loc(1)) && ...
		    (new_y >= NeuronParams.synapse_loc(2)) && ...
		    (new_x <= NeuronParams.synapse_loc(3)) && ...
		    (new_y <= NeuronParams.synapse_loc(4)))
			synapse_ret.AMPAR{i}.in_synapse = 1;
		end
	end
end


% Function to compute the new x and y coordinates of the AMPAR
function [new_x, new_y, new_y_cntr]=move_AMPAR(D, curr_x, curr_y, y_cntr, NeuronParams)

% Radius of the soma
radius=NeuronParams.soma_dia/2;

% The AMPAR will laterally diffuse along the surface of the soma
% if it has reached the surface
if (abs(((curr_x^2)+(curr_y^2)) - (radius^2)) < 1e-8)
	% AMPAR is at the surface of the soma
	if (curr_x < 0)
		del_x=D/sqrt(2);
	else
		del_x=-D/sqrt(2);
	end

	new_x=curr_x+del_x;
	new_y=sqrt((radius^2)-(new_x^2));
	new_y_cntr=0;

	return;
end

% Check if AMPAR is at the 45 degree boundary
if (curr_x == curr_y)
	% Move to the left of the boundary line
	del_x=-D/sqrt(2);
elseif (-curr_x == curr_y)
	% Move to the right of the boundary line
	del_x=D/sqrt(2);
else
	% toss a coin to determine if the AMPAR moves along the +ve or -ve
	% x-axis
	rand_val=rand();
	if (rand_val > 0.5)
		% Move along the +ve x-axis
		del_x=D/sqrt(2);
	else
		% Move along the -ve x-axis
		del_x=-D/sqrt(2);
	end
end

if (curr_y == 0)
	% If the AMPAR vesicle is at the origin then move along the +ve y-axis
	del_y=D/sqrt(2);

	if (y_cntr < 2)
		new_y_cntr=y_cntr+1;
	else
		new_y_cntr=0;
	end
else
	% Move more often along the positive y-direction
	if (y_cntr < 2)
		del_y=D/sqrt(2);
		new_y_cntr=y_cntr+1;
	else
		new_y_cntr=0;

		% toss a coin to determine if the AMPAR moves along the +ve or -ve
		% y-axis
		rand_val=rand();
		if (rand_val > 0.5)
			% Move along the +ve y-axis
			del_y=D/sqrt(2);
		else
			% Move along the -ve y-axis
			del_y=-D/sqrt(2);
		end
	end
end

% Compute new location
new_x=curr_x+del_x;
new_y=curr_y+del_y;

% Check if we are crossing the boundaries
max_x=radius*cosd(45);
min_x=-1*max_x;

if (new_x > max_x)
	new_x=max_x;
elseif (new_x < min_x)
	new_x=min_x;
end

if (new_y > radius)
	new_y=radius;
	new_x=0;
elseif (new_y < 0)
	new_y=0;
end

% Check if the AMPAR is crossing the 45 degree lines
if (new_x >= 0)
	if (new_y < new_x)
		new_y = new_x;
	end
else
	if (new_y < -new_x)
		new_y = -new_x;
	end
end

% Check that the AMPAR does not go beyond the soma surface
if (((new_x^2)+(new_y^2)) > (radius^2))
	new_y=sqrt((radius^2)-(new_x^2));
end


% Update the bound and unbound times of the AMPARs in the PSD
function [bound_times_ret, unbound_times_ret]=update_bound_unbound_times(bound_times, unbound_times, bounded, avail_idx)

bound_idx=find(bounded == 1);
bound_avail_idx=intersect(avail_idx, bound_idx); 

unbound_idx=find(bounded == 0);
unbound_avail_idx=intersect(avail_idx, unbound_idx); 

bound_times_ret=bound_times;
unbound_times_ret=unbound_times;

for (i=1:length(bound_avail_idx))
	idx=bound_avail_idx(i);
	bound_times_ret(i)=bound_times_ret(i)+1;
end

for (i=1:length(unbound_avail_idx))
	idx=unbound_avail_idx(i);
	unbound_times_ret(i)=unbound_times_ret(i)+1;
end


% Bind available AMPARs in the PSD that have stayed unbounded long enough
function [bound_times_ret, unbound_times_ret, bounded_ret, N_bounded]=...
	bind_AMPARs(bound_times, unbound_times, bounded, avail_idx, bind_time_thresh, N_bindings)

unbound_idx=find(bounded == 0);
idx_avail_for_binding=intersect(avail_idx, unbound_idx); 
max_bindings=length(idx_avail_for_binding);

bound_times_ret=bound_times;
unbound_times_ret=unbound_times;
bounded_ret=bounded;
N_bounded=0;

for (i=1:max_bindings)
	idx=idx_avail_for_binding(i);

	if (unbound_times_ret(idx) >= bind_time_thresh)
		bounded_ret(idx)=1;
		bound_times_ret(idx)=0;
		unbound_times_ret(idx)=0;
		N_bounded=N_bounded+1;

		if (N_bounded >= N_bindings)
			return;
		end
	end
end



% Unbind available AMPARs in the PSD that have stayed bounded long enough
function [bound_times_ret, unbound_times_ret, bounded_ret, N_unbounded]=...
	unbind_AMPARs(bound_times, unbound_times, bounded, avail_idx, unbind_time_thresh, N_unbindings)

bound_idx=find(bounded == 1);
idx_avail_for_unbinding=intersect(avail_idx, bound_idx); 
max_unbindings=length(idx_avail_for_unbinding);

bound_times_ret=bound_times;
unbound_times_ret=unbound_times;
bounded_ret=bounded;
N_unbounded=0;

for (i=1:max_unbindings)
	idx=idx_avail_for_unbinding(i);

	if (bound_times_ret(idx) >= unbind_time_thresh)
		bounded_ret(idx)=0;
		bound_times_ret(idx)=0;
		unbound_times_ret(idx)=0;
		N_unbounded=N_unbounded+1;

		% fprintf(1, 'Unbinding AMPAR\n');

		if (N_unbounded >= N_unbindings)
			return;
		end
	end
end







