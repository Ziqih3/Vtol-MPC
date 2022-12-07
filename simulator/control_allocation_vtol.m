classdef control_allocation_vtol < handle
    properties
    end
    
    properties(SetAccess=protected, GetAccess=public)
        Method control_allocation_types % The control allocation method
    end
    
    properties(SetAccess=protected, GetAccess=protected)
        NDI_L                 % L matrix (related to body-fixed thrust forces)
        NDI_M                 % M matrix (related to body-fixed thrust and reaction moments)
    end

    methods
        function obj = control_allocation_vtol(multirotor)
            obj.SetMethod(multirotor, control_allocation_types.Unified_PseudoInv);
        end
        
        function SetMethod(obj, multirotor, method)
            obj.Method = method;
            if method == control_allocation_types.NDI
                obj.InitializeNDIMethod(multirotor);
            end
        end
        
        function [rotor_speeds_squared, deflections , tilt,saturated] = CalcActuators(obj, mult, lin_accel, ang_accel)
        % Calculate the rotor speeds from the desired linear and angular accelerations
            
            persistent lin_accel_last
            persistent ang_accel_last
            persistent rotor_speeds_squared_last
            persistent deflections_last
            persistent saturated_last
            if isempty(lin_accel_last)
                lin_accel_last = zeros(3, 1);
            end
            if isempty(deflections_last)
                deflections_last = zeros(3, 1);
            end
            if isempty(ang_accel_last)
                ang_accel_last = zeros(3, 1);
            end
            if isempty(rotor_speeds_squared_last)
                rotor_speeds_squared_last = zeros(mult.NumOfRotors, 1);
            end
            if isempty(saturated_last)
                saturated_last = false;
            end
             
%             if isequal(lin_accel, lin_accel_last) && isequal(ang_accel, ang_accel_last)
%                 rotor_speeds_squared = rotor_speeds_squared_last;
%                 deflections = deflections_last;
%                 saturated = saturated_last;
%                 return;
%             end
        
            if obj.Method == control_allocation_types.NDI
                [rotor_speeds_squared, deflections] = obj.ActuatorCommands(mult, lin_accel, ang_accel);
            elseif obj.Method == control_allocation_types.Unified_PseudoInv
                [rotor_speeds_squared, deflections,tilt] = obj.ActuatorCommands_PseudoInv(mult, lin_accel, ang_accel);
            end

            saturation_flag = false;
            max_rotor_speeds = cell2mat(cellfun(@(s)s.MaxSpeedSquared, mult.Rotors, 'uni', 0));
            if any(rotor_speeds_squared > max_rotor_speeds)
                %mx = max((rotor_speeds_squared - max_rotor_speeds) ./ max_rotor_speeds);
                %rotor_speeds_squared = rotor_speeds_squared - mx * max_rotor_speeds - 1e-5;
                ind = rotor_speeds_squared > max_rotor_speeds;
                rotor_speeds_squared(ind) = max_rotor_speeds(ind);
                saturation_flag = true;
            end
            min_rotor_speeds = cell2mat(cellfun(@(s)s.MinSpeedSquared, mult.Rotors, 'uni', 0));
            if any(rotor_speeds_squared < min_rotor_speeds)
                ind = rotor_speeds_squared < min_rotor_speeds;
                rotor_speeds_squared(ind) = min_rotor_speeds(ind);
                saturation_flag = true;
            end
            
            if nargin > 1
                saturated = saturation_flag;
            end
            lin_accel_last = lin_accel;
            ang_accel_last = ang_accel;
            saturated_last = saturated;
            rotor_speeds_squared_last = rotor_speeds_squared;
        end
    end
    
    %% Private Methods
    methods(Access=protected)
        function InitializeNDIMethod(obj, multirotor)
        % Initialize the NDI method
        
            % Calculate L matrix (related to body thrust forces)
            obj.NDI_L = zeros(3, multirotor.NumOfRotors);
            for i = 1 : multirotor.NumOfRotors
               obj.NDI_L(:, i) = multirotor.Rotors{i}.GetThrustForcePerUnitInput();
            end

            % Calculate G matrix (related to body reaction moments)
            NDI_G = zeros(3, multirotor.NumOfRotors);
            for i = 1 : multirotor.NumOfRotors
               NDI_G(:, i) = multirotor.Rotors{i}.GetReactionMomentPerUnitInput();
            end
            
            % Calculate F matrix (related to body thrust moments)
            NDI_F = zeros(3, multirotor.NumOfRotors);
            for i = 1 : multirotor.NumOfRotors
                r = multirotor.Rotors{i}.Position;
                F = multirotor.Rotors{i}.GetThrustForcePerUnitInput();
                NDI_F(:, i) = cross(r, F);
            end
            
            obj.NDI_M = NDI_F + NDI_G;
        end
        function rotor_speeds_squared = NDIRotorSpeeds(obj, multirotor, lin_accel, ang_accel)
        % Calculate the rotor speeds from the desired linear and angular accelerations
        % using NDI method
            
            % Create the desired output matrix y
            y = [lin_accel; 
                ang_accel];
        
            % Get the rotation matrix
            RBI = multirotor.GetRotationMatrix();
            
            % Calculate eta_dot
            phi = deg2rad(multirotor.State.RPY(1));
            theta = deg2rad(multirotor.State.RPY(2));
            
            phi_dot = deg2rad(multirotor.State.EulerRate(1));
            theta_dot = deg2rad(multirotor.State.EulerRate(2));
            eta_dot = calc_eta_dot(phi, theta, phi_dot, theta_dot);
            
            % Calculate eta
            eta = [1,   sin(phi)*tan(theta), cos(phi)*tan(theta);
                   0, cos(phi), -sin(phi);
                   0, sin(phi) / cos(theta), cos(phi) / cos(theta)];

            % Calculate the A matrix in y = A + Bu
            NDI_M_Grav = zeros(3, 1);
            for i = 1 : multirotor.NumOfRotors
                r = multirotor.Rotors{i}.Position;
                G_motor = multirotor.Rotors{i}.MotorMass * physics.Gravity;
                G_motorB = RBI * G_motor;
                G_arm = multirotor.Rotors{i}.ArmMass * physics.Gravity;
                G_armB = RBI * G_arm;
                NDI_M_Grav = NDI_M_Grav + cross(r, G_motorB) + cross(r/2, G_armB);
            end
            if multirotor.HasEndEffector()
                r = multirotor.EndEffector.EndEffectorPosition;
                G_eeI = multirotor.EndEffector.EndEffectorMass * physics.Gravity;
                G_eeB = RBI * G_eeI;
                G_armI = multirotor.EndEffector.ArmMass * physics.Gravity;
                G_armB = RBI * G_armI;
                NDI_M_Grav = NDI_M_Grav + cross(r, G_eeB) + cross(r/2, G_armB);
            end
            
            A_moment = eta_dot * multirotor.State.Omega + eta * multirotor.I_inv * ...
                (NDI_M_Grav - cross(multirotor.State.Omega, multirotor.I * multirotor.State.Omega));
            A = [physics.Gravity; A_moment];
            
            % Calculate the B matrix
            B_force = RBI' * obj.NDI_L / multirotor.TotalMass;
            B_moment = eta * multirotor.I_inv * obj.NDI_M;
            B = [B_force; B_moment];
            
            % Calculate the rotor speeds
            rotor_speeds_squared = pinv(B) * (y - A); 
        end
        
        function [rotor_speeds_squared, deflections] = ActuatorCommands(obj, multirotor, lin_accel, ang_accel)
        % Calculate the rotor speeds from the desired linear and angular accelerations
        % using NDI method
            
            % Create the desired output matrix y
            y = [lin_accel; ang_accel];
            M_des = (ang_accel' * multirotor.I)';
            F_des = lin_accel * multirotor.TotalMass;
%             [rbw, ~, ~] = multirotor.CalcWindToBodyRotation(multirotor.State.AirVelocity);
            Va_i = multirotor.State.AirVelocity;
            q_bar = (Va_i' * Va_i) * physics.AirDensity / 2;

            del_a = M_des(1) / (multirotor.C_A * multirotor.WingSurfaceArea * multirotor.Wingspan * q_bar);
            del_e = ( M_des(2) - (multirotor.l3 - multirotor.l4)/2*F_des(3) ) / (multirotor.C_E * multirotor.WingSurfaceArea * multirotor.MeanChord * q_bar);
            del_r = M_des(3) / (multirotor.C_R * multirotor.WingSurfaceArea * multirotor.Wingspan * q_bar);
            
            del_a = min(max(del_a, -1), 1);
            del_e = min(max(del_e, -1), 1);
            del_r = min(max(del_r, -1), 1);
            
            a = 0.0185; %ramp function 1 slope parameter
            b = 35.217; %ramp function 1 position parameter
            
            f1 = min(max(0, a*(q_bar - b) + 0.5), 1);
            
            del_a = del_a * f1;
            del_e = del_e * f1;
            del_r = del_r * f1;
            
            deflections = [del_a, del_e, del_r];
%             disp(deflections)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            M_def = [multirotor.C_A * multirotor.WingSurfaceArea * multirotor.Wingspan *  q_bar * del_a;
                     multirotor.C_E * multirotor.WingSurfaceArea * multirotor.MeanChord * q_bar * del_e;
                     multirotor.C_R * multirotor.WingSurfaceArea * multirotor.Wingspan *  q_bar * del_r];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            M_r = M_des - M_def;
            
            angular_accel_r = (M_r' * multirotor.I_inv)';
            rotor_speeds_squared = NDIRotorSpeeds(obj, multirotor, lin_accel, angular_accel_r);
            
            control_sp = [M_des; F_des];
            
        end
        
        function [rotor_speeds_squared, deflections,tilt] = ActuatorCommands_PseudoInv(obj, multirotor, lin_accel, ang_accel)
            % Calculate the rotor speeds from the desired linear and angular accelerations
            % using NDI method
                
                % Create the desired output matrix y
                y = [lin_accel; ang_accel];
                M_des = (ang_accel' * multirotor.I)';
                F_des = lin_accel * multirotor.TotalMass;
                control_sp = [M_des; F_des];

                Va_i = multirotor.State.AirVelocity;
                q_bar = (Va_i' * Va_i) * physics.AirDensity / 2;
                
                actuator_trim = zeros(12, 1); % Change trim values if neccessary
                if any(multirotor.State.ServoAngles) ~= 0
                    actuator_trim(1:4) = 1;
                else
                    actuator_trim(1:4) = 0;
                end
                actuator_trim(5:6) = multirotor.State.ServoAngles(1);
                actuator_trim(7:8) = multirotor.State.ServoAngles(2);
                
                effectiveness_matrix = calc_eff_mat(q_bar, actuator_trim);
                control_trim = effectiveness_matrix * actuator_trim;
                actuator_sp = actuator_trim + pinv(effectiveness_matrix) * (control_sp - control_trim);
                
                rotor_speeds_squared = actuator_sp(1:4);
                tilt = actuator_sp(5:8);
                deflections = [abs(actuator_sp(9)+actuator_sp(10))/2  actuator_sp(11) actuator_sp(12)];
                
        end
    end
end

%% Other functions
function eta_dot = calc_eta_dot(phi, theta, phi_dot, theta_dot)
    eta_dot_11 = 0;
    eta_dot_12 = sin(phi)*(tan(theta)^2 + 1)*theta_dot + cos(phi)*tan(theta)*phi_dot;
    eta_dot_13 = cos(phi)*(tan(theta)^2 + 1)*theta_dot - sin(phi)*tan(theta)*phi_dot;

    eta_dot_21 = 0;
    eta_dot_22 = -phi_dot*sin(phi);
    eta_dot_23 = -phi_dot*cos(phi);

    eta_dot_31 = 0;
    eta_dot_32 = (cos(phi)*phi_dot)/cos(theta) + (sin(phi)*sin(theta)*theta_dot)/cos(theta)^2;
    eta_dot_33 = (cos(phi)*sin(theta)*theta_dot)/cos(theta)^2 - (sin(phi)*phi_dot)/cos(theta);

    eta_dot = [eta_dot_11 eta_dot_12 eta_dot_13;
               eta_dot_21 eta_dot_22 eta_dot_23;
               eta_dot_31 eta_dot_32 eta_dot_33];
end

function effectiveness_matrix = calc_eff_mat(q_bar, trim)
    
    
    %           1         2        3         4
    
    % Px = [ 0.1515   -0.1515   0.1515   -0.1515];
    % Py = [-0.245     0.245    0.245    -0.245];
    Px = [ 0.1515    0.1515  -0.1515   -0.1515];
    Py = [ 0.245    -0.245   -0.245     0.245];
    Pz = zeros(1, 4);
    Ct = 5.0;
    Km = 0.07;
    ro = 1.225;
    S = 0.4266;
    b = 2.0;
    c_bar = 0.2;
    Cla = 0.1;
    Cme = 0.5;
    Cnr = 0.5;
    effectiveness_matrix = [-Py(1) * Ct*cos(trim(5)) - Ct * Km * sin(trim(5)),				-Py(2) * Ct*cos(trim(6)) - Ct * Km * sin(trim(6)),			 -Py(3) * Ct*cos(trim(7)) + Ct * Km * sin(trim(7)),			    -Py(4) * Ct*cos(trim(8)) + Ct * Km * sin(trim(8)),			 Py(1) * Ct*trim(1) *sin(trim(5)) - Ct * Km * trim(1) *cos(trim(5)),    	     Py(2) * Ct*trim(2) *sin(trim(6)) - Ct * Km * trim(2) *cos(trim(6)),		 		Py(3) * Ct*trim(3) *sin(trim(7)) + Ct * Km * trim(3) *cos(trim(7)),		            Py(4) * Ct*trim(4) *sin(trim(8)) + Ct * Km * trim(4) *cos(trim(8)),  		-q_bar*S*b*Cla,	    q_bar*S*b*Cla,	    0.0, 		 		    0.0;
                             Ct*(Px(1) * cos(trim(5)) + Pz(1) * sin(trim(5))),  			 Ct*(Px(2) * cos(trim(6)) + Pz(2) * sin(trim(6))),			  Ct*(Px(3) * cos(trim(7)) + Pz(3) * sin(trim(7))),			 	 Ct*(Px(4) * cos(trim(8)) + Pz(4) * sin(trim(8))),			 Ct*trim(1) *(-Px(1) * sin(trim(5)) + Pz(1) * cos(trim(5))), 	     		     Ct*trim(2) *(-Px(2) *sin(trim(6)) + Pz(2) * cos(trim(6))),							Ct*trim(3) *(-Px(3) *sin(trim(7)) + Pz(3) * cos(trim(7))),							Ct*trim(4) *(-Px(4) *sin(trim(8)) + Pz(4) *cos(trim(8))),  					 0.0, 		 	    0.0, 		 	    q_bar*S*c_bar*Cme,	 	0.0; 			
                            -Py(1) * Ct*sin(trim(5)) + Ct * Km * cos(trim(5)),				-Py(2) * Ct*sin(trim(6)) + Ct * Km * cos(trim(6)),			 -Py(3) * Ct*sin(trim(7)) - Ct * Km * cos(trim(7)),			    -Py(4) * Ct*sin(trim(8)) - Ct * Km * cos(trim(8)),			-Py(1) * Ct*trim(1) *cos(trim(5)) - Ct * Km * trim(1) *sin(trim(5)),   	        -Py(2) * Ct*trim(2) *cos(trim(6)) - Ct * Km * trim(2) *sin(trim(6)),		  	   -Py(3) * Ct*trim(3) *cos(trim(7)) + Ct * Km * trim(3) *sin(trim(7)), 	 		   -Py(4) * Ct*trim(4) *cos(trim(8)) + Ct * Km * trim(4) *sin(trim(8)),  		 0.0, 		 	    0.0, 		 	    0.0, 		 		 	q_bar*S*b*Cnr; 	
                             Ct * sin(trim(5)),	 										     Ct * sin(trim(6)),			   							      Ct * sin(trim(7)),											 Ct * sin(trim(8)),										     Ct * trim(1) *cos(trim(5)), 									     		     Ct * trim(2) *cos(trim(6)),											 			Ct * trim(3) *cos(trim(7)), 											 			Ct * trim(4) *cos(trim(8)),   												 0.0, 		 	    0.0, 		 	    0.0, 		 		 	0.0; 			
                             0.0,  			 												 0.0,  							   			 				  0.0,  													 	 0.0,	 													 0.0,		 		  					     					     			 0.0, 			      							 									0.0, 			    	  											   	 		   	0.0, 			           								 					 0.0, 		 	    0.0, 		 	    0.0, 		 		 	0.0; 			
                            -Ct * cos(trim(5)),	 						    			    -Ct * cos(trim(6)),			   							     -Ct * cos(trim(7)),											-Ct * cos(trim(8)),										     Ct * trim(1) *sin(trim(5)), 									     		     Ct * trim(2) *sin(trim(6)),											 			Ct * trim(3) *sin(trim(7)), 											 			Ct * trim(4) *sin(trim(8)),   								 				 0.0, 		 	    0.0, 		 	    0.0, 		 		 	0.0]; 			
    
end