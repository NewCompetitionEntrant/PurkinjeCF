close all; clear all; clc;

g_init=8.5;
g_delta=4;
g_track=zeros(1000/0.25, 1);
t_track=zeros(1000/0.25, 1);
g=g_init;
t_fire=0;
g_fire=0;
tau=150;
iter=0;
for t=0:0.25:1000-0.25
	if (g > g_init)
		g=g_fire*exp(-(t-t_fire)/tau);
	end

	if (g < g_init)
		g = g_init;
	end

	if ((t ~= 0) && (t < 500) && (mod(t, 20) == 0))
		g=g+g_delta;
		g_fire=g;
		t_fire=t;
	end

	iter=iter+1;
	g_track(iter)=g;
	t_track(iter)=t;
end

figure();
plot(t_track, g_track, 'b-');
grid on;
