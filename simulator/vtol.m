classdef vtol < multirotor
    properties
        WingDirections = {[0; 1; 0], [0; -1; 0]};
        WingLengths = [1; 1]; % in meters
        WingSurfaceArea = 0.44; % in m^2
        AngleOfIncident = 0; % in degrees
        Wingspan = 2; % in meters
        MeanChord = 0.22; % in meters
%         C_A = 0.03; %Aileron constant
%         C_E = 0.03; %Elevator constant
%         C_R = 0.03; %Rudder constant
        C_A = 0.11730; %Aileron constant
        C_E = 0.55604; %Elevator constant
        C_R = 0.08810; %Rudder constant
        r_w = 0.0625; %wing center of lift wrt, center of gravity in meters
        r_t = 0.6385; %tail center of lift wrt, center of gravity in meters
        S_A = 0.0720; %aileron surface area in m^2
        S_E = 0.03; %elevator surface area in m^2
        S_R = 0.008; %rudder surface area in m^2
        l3 = 0.5189; %front props delta x from center of mass in meters
        l4 = 0.4574; %back props delta x from center of mass in meters
        L0 = 0.3642; %tilting rotors y offset
        fw_p_lim_min = -45;
        fw_p_lim_max = 45;
        fw_t_climb_r_sp = 3;
        fw_t_sink_r_sp = 2;
        is_transition = false;
    end

    properties (SetAccess=protected, GetAccess=public)
        % Add your properties that can only be set internally
    end
    
    properties (SetAccess=protected, GetAccess=protected)
        % Add your properties that can only be seen internally
    end
    
    %% Public methods
    methods
        function obj = vtol(ArmAngles, RotationDirections)
            obj = obj@multirotor(ArmAngles, RotationDirections);
        end
        
        function [wrench, aeromoments,rotor_speed_squared] = CalcGeneratedWrench(obj, plantinput)
            wrench = CalcGeneratedWrench@multirotor(obj, plantinput);
            
            %wrench(6) = wrench(6) - 0.4 * physics.AirDensity * 0.52 / 2 * norm(obj.State.Velocity(1:2))^2; % Add vertical lift
%             obj.CalcAerodynamicMoment(obj.State.AirVelocity)
%             wrench(1:3) = wrench(1:3)/(1 + norm(obj.State.AirVelocity)) + obj.CalcDeflectionMoment(obj.State.AirVelocity, plantinput);
            aeromoments = obj.CalcAerodynamicMoment(obj.State.AirVelocity);
            %aeromoments = obj.CalcAerodynamicMoment(obj.State.AirVelocity,lin_accel, euler_accel);
            wrench(1:3) = wrench(1:3) + obj.CalcDeflectionMoment(obj.State.AirVelocity, plantinput);
            [force,rotor_speed_squared] = obj.CalcAerodynamicForce(obj.State.AirVelocity);
            wrench(4:6) = wrench(4:6) + force;
            
            %% Checking force in Body Frame X-dir
            %% Changing properties for checking dynamics
%             wrench(1:3) = [0;0;0];
%             wrench(4:5) = [-sind(obj.State.RPY(3))*29.01265,cosd(obj.State.RPY(3))*29.01625];
%             wrench(6) = 0; %balance in z
        end
        
        function new_state = CalcNextState(obj, wrench, Thrust_Body,tf_sensor_wrench, ...
                wind_force, plantinput, dt, is_collision, collision_normal, ...
                air_velocity)
            
            new_state = CalcNextState@multirotor(obj, wrench, Thrust_Body,tf_sensor_wrench, ...
                wind_force, plantinput, dt, is_collision, collision_normal, air_velocity);
            [~, alpha, beta] = CalcWindToBodyRotation(obj, air_velocity);
            new_state.AngleOfAttack = alpha;
            new_state.SideSlipAngle = beta;
            da = min(0.1, abs(plantinput.AileronLeftRate));
            if(plantinput.AileronLeftRate < 0)
                da = -da;
            end
            new_state.AileronLeftPos = obj.State.AileronLeftPos + da;
            new_state.AileronRightPos = -new_state.AileronLeftPos;
            new_state.ElevatorPos = plantinput.ElevatorRate;
            new_state.RudderPos = plantinput.RudderRate;
            new_state.AeroMoments = plantinput.Aeromoment;
        end
    
    end
    %% Private methods
    methods (Access = private)
        function [force,rotor_speed_squared] = CalcAerodynamicForce(obj,Va_i)
            [rwb, alpha, ~] = obj.CalcWindToBodyRotation(Va_i);
            alpha = deg2rad(alpha);
            RIB = obj.GetRotationMatrix(); % inertial to body
            RBI = RIB'; % body to inertial
            q_bar = (Va_i' * Va_i) * physics.AirDensity / 2;
            
            c_y = 0; % TODO: Later
            %c_x = get_cd(alpha); % TODO: Add beta-dependent part later
            %c_z = get_cl(alpha);
            c_z = 0.35 + 0.11 * alpha; %C
            c_d = 0.01 + 0.2 * alpha * alpha; %C_D,a from paper 
            
            
            drag = q_bar * obj.WingSurfaceArea * c_d;
            lateral = q_bar * obj.WingSurfaceArea * c_y;
            lift = q_bar * obj.WingSurfaceArea * c_z;

            force = rwb * [-drag; lateral;-lift];
            
            weight = obj.Mass*physics.Gravity;
            
            BF_weight = RBI * weight;
            
            % Update rotor_speed_squared based on dynamics
            rotor_speed_squared = zeros(obj.NumOfRotors,1);
            for i = 1 : obj.NumOfRotors
                rotor_speed_squared(i) = (-force(1) + BF_weight(1)) / (obj.NumOfRotors * obj.Rotors{i}.ThrustConstant);
            end
            
            force = RBI * force;
        end
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function moment = CalcAerodynamicMoment(obj, Va_i) %TODO
            [rwb, ~, ~] = obj.CalcWindToBodyRotation(Va_i);
            
            Va_b = rwb*Va_i;
    
            q_bar = (Va_b' * Va_b) * physics.AirDensity / 2;
            
            p = obj.State.Omega(1);
            q = obj.State.Omega(2);
            r = obj.State.Omega(3);
            
            %dimensionless angular rate
            p_tilde = (obj.Wingspan * p) / (2 * norm(Va_b)^2); % why there's a square for Va_b
            q_tilde = (obj.MeanChord * q) / (2 * norm(Va_b)^2);
            r_tilde = (obj.Wingspan * r) / (2 * norm(Va_b)^2);
                      
            roll_moment = q_bar * obj.WingSurfaceArea * obj.Wingspan * p_tilde;
            pitch_moment = q_bar * obj.WingSurfaceArea * obj.Wingspan * q_tilde;
            yaw_moment = q_bar * obj.WingSurfaceArea * obj.Wingspan * r_tilde;
%             c_i = get_ci(alpha);
%             c_m = get_cm(alpha);
%             c_n = get_cn(alpha);
%             
%             roll_moment = q_bar * obj.WingSurfaceArea * obj.Wingspan * c_i;
%             pitch_moment = q_bar * obj.WingSurfaceArea * obj.MeanChord * c_m;
%             yaw_moment = q_bar * obj.WingSurfaceArea * obj.Wingspan * c_n;
            
            moment = [roll_moment; pitch_moment; yaw_moment];
        end
        
        function moment = CalcDeflectionMoment(obj, Va_i, plantinput)
            [rwb, ~, ~] = obj.CalcWindToBodyRotation(Va_i);
            Va_b = rwb*Va_i;
            
            q_bar = (Va_b' * Va_b) * physics.AirDensity / 2;
            roll_moment = q_bar * obj.S_A * obj.C_A * 2 * plantinput.AileronLeftRate;
            pitch_moment = q_bar * obj.S_E * obj.C_E * plantinput.ElevatorRate;
            yaw_moment = q_bar * obj.S_R * obj.C_R * plantinput.RudderRate;
            
            moment = [roll_moment; pitch_moment; yaw_moment];
        end

        function [R_WB, alpha, beta] = CalcWindToBodyRotation(obj, Va_i)
            rbi = obj.GetRotationMatrix();
            Va_b = rbi * Va_i;
            %Calculate angle of attack a and sideslip angle b
            a = CalcAngleOfAttack(Va_b) + obj.AngleOfIncident;
            b = CalcSideSlipAngle(Va_b);
%             a = 10;
            b = 0;
           
            %Rotation matrix from body to wind

            R_BW = [cosd(a) * cosd(b),    sind(b),     sind(a)*cosd(b);
                    -sind(b) * cosd(a),   cosd(b),     -sind(a)*sind(b);
                    -sind(a)         ,   0  ,        cosd(a)];
%             R_BW = [cosd(b) * cosd(a),  -sind(b) * cosd(a),     -sind(a)
%                     sind(b),            cosd(b),                0
%                     cosd(b) * sind(a),  - sind(b) * sind(a),    cosd(a) ];
            R_WB = R_BW.';
         
            if nargout > 1
                alpha = a;
                beta = b;
            end
        end
    end
    
end

%% Other function

function alpha = CalcAngleOfAttack(Va_b)
%     alpha = atan2d(Va_b(3), Va_b(1));
    alpha = atand(Va_b(3)/Va_b(1));
    if isnan(alpha)
        alpha = 0;
    end
end

function beta = CalcSideSlipAngle(Va_b)
    beta = asind(Va_b(2) / norm(Va_b));
    if isnan(beta)
        beta = 0;
    end
end

function cl = get_cl(alpha)
    data = load('C_l');
    cl = fixpt_interp1([-8.5:.25:13.75, 14.5,14.75, 15],data.C_l,alpha,sfix(8),2^-3,sfix(16), 2^-14,'Nearest');
end
function cd = get_cd(alpha)
    data = load('C_d');
    cd = fixpt_interp1([-8.5:.25:13.75, 14.5,14.75, 15],data.C_d,alpha,sfix(8),2^-3,sfix(16), 2^-14,'Nearest');
end
function cm = get_cm(alpha)
data = load('C_m');
cm = fixpt_interp1([-8.5:.25:13.75, 14.5,14.75, 15],data.C_m,alpha,sfix(8),2^-3,sfix(16), 2^-14,'Nearest');
end