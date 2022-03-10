 
###############################
# parameters and sets
###############################

param bar_t_max > 0, integer;   # number of time periods
param bar_t integer <= bar_t_max;
set T := 1..bar_t;				#set of time periods in the considered time horizon
set T_max := 1..bar_t_max;
param bar_n > 0, integer;     	# number of turbine/pump units
set J := 1..bar_n;            	# set of generation units

param c{J} > 0;					#cost for generating electricity by unit j
param c_u{J} > 0;				#startup cost of unit j
param D{T_max};						#electricity demand to be satisfied at time period t
param R{T_max}; 					#amount of spinning reserves at time period t
param v0{J} default 0; 			#on/off status for unit j at time period 0
param R_U{J}; 					#maximum ramp-up rate of unit j
param R_D{J}; 					#maximum ramp-down rate of unit j
param S_U{J}; 					#maximum startup rate of unit j
param S_D{J}; 					#maximum shutdown rate of unit j
param p0{J}; 					#production of unit j before the planning horizon
param T_U{J}; 					#minimum time the unit j should be on
param T_D{J}; 					#minimum time the unit j should be off
param U_t{J}; 					#number of time periods that j is required to be on at the start of the planning horizon
param D_t{J}; 					#number of time periods that j is required to remain off at the start of the planning horizon
param p_max_bound{J}; 			#maximum power output by unit j
param p_min_bound{J}; 			#minimum power output by unit j

###############################
# variables and simple bounds
###############################

var p{j in J, t in T} >= 0;		#electricity produced by unit j at time period t
var p_max{j in J, t in T} >= 0;	#maximum electricity generated at each time period by unit j
var v{j in J, t in T} binary; 	#indicates if unit j is on in time period t
var y{j in J, t in T} binary; 	#indicates if unit j is switching on in time period t
var z{j in J, t in T} binary; 	#indicates if unit j is shutting down in time period t


###############################
# Objective function
###############################

minimize cost_function:
	sum{j in J, t in T} (c[j]*p[j,t] + c_u[j]*y[j,t]);

################################
# Constraints
################################

subject to demand{t in T}:
	sum{j in J}p[j, t] = D[t];

subject to spinning_reserve{t in T}:
	sum{j in J} p_max[j, t] >= D[t] + R[t];
	
	
####################### Start up constrainst #############################
subject to status_update{j in J, t in T}:
	if t==1 then v0[j] - v[j,t] + y[j,t] - z[j,t]
          else v[j,t-1] - v[j,t] + y[j,t] - z[j,t]
	= 0;

#Redundant with uptime/downtime constraint ?
subject to status_disjoint{j in J, t in T}:
	y[j,t] + z[j,t] <=1;
	
	
	
#========================================================================#


####################### Ramping constraint ##############################

subject to ramping_up{j in J, t in T}:
	if t==1 then p[j,t] - p0[j] - (R_U[j]*v0[j] + S_U[j] * y[j, t])
          else p[j,t] - p[j, t-1] - (R_U[j]*v[j, t-1] + S_U[j] * y[j, t])
	<= 0;
	
subject to ramping_down{j in J, t in T}:
	if t==1 then p0[j] - p[j,t]  - (R_D[j]*v[j, t] + S_D[j] * z[j, t])
          else p[j, t-1] - p[j,t]  - (R_D[j]*v[j, t] + S_D[j] * z[j, t])
	<= 0;
#=======================================================================#


subject to uptime{j in J, t in T : t>U_t[j]}:
	sum{k in {(t-T_U[j]+1)..t} : k>=1} y[j, k] <= v[j, t];

subject to downtime{j in J, t in T : t>D_t[j]}:
	sum{k in {(t-T_D[j]+1)..t} : k>=1} z[j, k] <= (1 - v[j, t]);
	
########################### generation limit ################################

subject to lower_p{j in J, t in T}:
	p_min_bound[j] * v[j, t] <= p[j, t];
	
subject to upper_p{j in J, t in T}:
	p[j, t] <= p_max_bound[j] * v[j, t];


subject to lower_p_max{j in J, t in T}:
	p_min_bound[j] * v[j, t] <= p_max[j, t];
	
subject to upper_p_max{j in J, t in T}:
	p_max[j, t] <= p_max_bound[j] * v[j, t];
	
subject to p_max_ramping_up{j in J, t in T}:
	if t==1 then p_max[j,t] - p0[j] - (R_U[j]*v0[j] + S_U[j] * y[j, t])
          else p_max[j,t] - p[j, t-1] - (R_U[j]*v[j, t-1] + S_U[j] * y[j, t])
	<= 0;

subject to p_max_ramping_down{j in J, t in 1..bar_t-1}:
	p_max[j, t] - (p_max_bound[j] * (v[j, t] - z[j, t+1]) + z[j, t+1] * S_D[j])<=0;


#===========================================================================#







