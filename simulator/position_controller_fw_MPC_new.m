classdef position_controller_fw_MPC_new < pid_controller
    properties
        RateLimits = [7; 7; 9]; % in m/s
        OutputMax = [5; 5; 10]; % in m/s^2
        VelMax = [40;40;20];
        Q = zeros(6);
        q_x = 1;
        c = [1 1 1];
        k_x = 1;
        WingSurfaceArea = 0.44;
        Mass;
%%%%%%%%
        q_ref = [5, 5, 10]';
        alpha = 0.1; % for cost function calculation\

        Q_rpy = diag([20 20 0]);
        Q_rpy_dot = diag([5 5 0]);
        Q_u = diag([0.0025 1 100 200 50]);
        Q_t = 40;
        
        I_inv = diag([10.685,5.7465,4.6678]);
        R = [2  0   0;
             0  1   0;
             0  0   2];
        aoi = deg2rad(10);

        acc_Max = [1; 1; 6]; % in m/s^2

    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    methods

        function [thrust_des,rpy_des] = CalculateControlCommand(obj, mult, pos_des, vel_des, yaw_des, acc_des, time)
            import casadi.*
                    if isempty(vel_des)
                            vel_des = zeros(3, 1);

                    end
        if isempty(acc_des)
            acc_des = zeros(3, 1);
        end


        obj.Mass = mult.Mass;

        Init_Vel = mult.State.AirVelocity;
        Init_rpy = mult.State.RPY;
        Init_rpy_dot = mult.State.EulerRate;
        Init_Tilt = mult.State.ServoAngles(1);
        Init_last_rpy_dot = mult.State.LastEulerRate;
        Init_last_thrust = norm(mult.State.LastThrust); % Assuming tilting is the same in all servos

        dt = time - obj.LastTime;


        % Horizon
        N = 20;
        
        % Velocity
        V_x = SX.sym('V_x');
        V_y = SX.sym('V_y');
        V_z = SX.sym('V_z');
        V = [V_x; V_y; V_z];
        
        % Attitude 
        Roll = SX.sym('Roll');
        Pitch = SX.sym('Pitch');
        Yaw = SX.sym('Yaw');
        rpy = [Roll; Pitch; Yaw];
        
        % Body rate setpoint parametrization
        Roll_dot = SX.sym('Roll_dot');
        Pitch_dot = SX.sym('Pitch_dot');
        Yaw_dot = SX.sym('Yaw_dot');
        rpy_dot = [Roll_dot; Pitch_dot; Yaw_dot];
        
        % Tilt angle parametrization 
        tilt_angle = SX.sym('tilt_angle');
        
        % Last body rate parametrization
        Roll_dot_last = SX.sym('Roll_dot_last');
        Pitch_dot_last = SX.sym('Pitch_dot_last');
        Yaw_dot_last = SX.sym('Yaw_dot_last');
        last_rpy_dot = [Roll_dot_last; Pitch_dot_last; Yaw_dot_last];
        
        % 1x1

        last_thrust = SX.sym('last_thrust');

        % 14x1
        state = [V; rpy; rpy_dot; tilt_angle; last_rpy_dot; last_thrust];
        n_states = length(state);
        % 3x1
        Roll_MPC = SX.sym('Roll_MPC');
        Pitch_MPC = SX.sym('Pitch_MPC');
        Yaw_MPC = SX.sym('Yaw_MPC');
        rpy_MPC = [Roll_MPC; Pitch_MPC; Yaw_MPC];
        
        % Tilt speed parametrization
        tilt_speed = SX.sym('tilt_speed'); 
        % Thrust
        Thrust = SX.sym('Thrust');
        % Control input
        controls= [Thrust; tilt_speed; rpy_MPC];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J cost
           
%         % reference tracking
%         
%          x = norm(obj.q_ref' * ( - V_sp)); % before approximation
%         
%         J_ref = 2*obj.alpha*log(1 + exp(x/obj.alpha)) - x - 2*obj.alpha*log(2);
% 
%         % state and input cost
%         J_state_input = rpy'*obj.Q_rpy*rpy+ rpy_dot'*obj.Q_rpy_dot*rpy_dot + u'*obj.Q_u*u +obj.q_t*(last_thrust - Thrust)^2;    
%         % soft constraints
%         %a = −0.332, b = 13.35, c = −0.477, d = −2.303.
%         % coeff q_x = 1
%         J_soft_tilt= obj.q_x*exp(-0.332*V_x*tilt_angle + 13.35*tilt_angle + (-0.477*V_x) + (-2.303));
%         % coeff k_x ,offset c[x,y,z] = 1 
%         J_soft_vel = obj.q_ref(1)*(exp(-3*(V_x + 1)) + obj.k_x*V_x - obj.c(1))...
%             + obj.q_ref(2)*(exp((-V_y - 2)) + exp(V_y - 2) - obj.c(2))...
%             + obj.q_ref(3)*(exp((-V_z - 2)) + exp(V_y - 2) - obj.c(3));
%         



        n_controls = length(controls);

        F_aero = obj.CalcAeroForce(V);
%       disp(tilt_angle)
        accel = obj.CalcAccel(F_aero,Thrust,rpy,tilt_angle);
        accel = accel + physics.Gravity;
       
        rpy_sp = rpy_MPC + [0;0;Yaw];
        body_rate_sp = diag([2,2,2]) * (rpy_sp - rpy);
        torque = diag([2,2,2]) * (body_rate_sp - rpy_dot) + diag([1,1,1]) * (last_rpy_dot - rpy_dot);
        
        %state update EQ
        rhs = [accel;
               rpy_dot;
               obj.I_inv*(torque);
               tilt_speed;
               last_rpy_dot;
               last_thrust];
        
        f = Function('f',{state,controls},{rhs});  %nonlinear mapping function f(x,u)
        U = SX.sym('U',n_controls,N);
        P = SX.sym('P',n_states + n_states);
        X = SX.sym('X',n_states,(N+1));


        % A vector that represents the states over the optimization problem
%        objective_function = J_ref + J_state_input + J_soft;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % eq contraint 
        %init var!!!!!!!!!!!!
        st = X(:,1);
        g = []; % constraints vector
        g = [g;st-P(1:14)];  %init constraint
            
        objective_function = 0;

        h = dt;

        for k = 1:N
            st = X(:,k);
            con = U(:,k);
            st_next = X(:,k+1);
            obj_soft_tilt= obj.q_x*exp(-0.332*P(n_states +1)*P(n_states+10) + 13.35*(n_states+10) + (-0.477*P(n_states+1)) + (-2.303));
            obj_soft_vel = obj.q_ref(1)*(exp(-3*(P(n_states +1) + 1)) + obj.k_x*P(n_states +1) - obj.c(1))...
                        + obj.q_ref(2)*(exp((-P(n_states +2)- 2)) + exp(P(n_states +2) - 2) - obj.c(2))...
                        + obj.q_ref(3)*(exp((-P(n_states +3) - 2)) + exp(P(n_states +2) - 2) - obj.c(3));
            obj_stateinput = (st(4:6) - P((n_states+4):(n_states+6)))' * obj.Q_rpy * (st(4:6) - P((n_states+4):(n_states+6)))...
                            +(st(7:9) - P((n_states+7):(n_states+9)))' * obj.Q_rpy_dot * ...
                            (st(7:9) - P((n_states+7):(n_states+9)))...
                            + con' * obj.Q_u * con + obj.Q_t * (P(n_states + 14) - con(1)) ^ 2;
             
            x = obj.q_ref'*(st(1:3)-P((n_states):(n_states+2)));

            obj_reference = 2*obj.alpha*log(1 + exp(x/obj.alpha)) - x - 2*obj.alpha*log(2);



            objective_function = objective_function +  obj_stateinput  + obj_soft_tilt + obj_soft_vel+obj_reference;

            k1 = f(st, con);
            k2 = f(st + h/2*k1,con);
            k3 = f(st + h/2*k2,con);
            k4 = f(st + h*k3, con);
            st_RK4_next = st + h/6* (k1 +2*k2 + 2*k3 +k4);
            % g = [f(N-1) - x(N)]
            g = [g;st_next - st_RK4_next];
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

        
        % make the decision variables one column vector
        OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_controls*N,1)];
        nlp_prob = struct('f',objective_function,'x',OPT_variables,'g',g,'p',P);
        opts = struct;
        opts.ipopt.max_iter = 100;
        opts.ipopt.print_level = 0;
        opts.print_time = 0;
        opts.ipopt.acceptable_tol = 1e-8;
        opts.ipopt.acceptable_obj_change_tol = 1e-6;
        solver = nlpsol('solver','ipopt',nlp_prob,opts);
        %------------------------------------------------------
        args = struct;
        % equality constraints

        args.lbg(1:14*(N+1)) = -1e-20;
        args.ubg(1:14*(N+1)) = 1e-20;

        % input and states constraints
        args.lbx(1:14:14*(N+1),1) = -40; args.ubx(1:14:14*(N+1),1) = 40; %V in m/s
        args.lbx(2:14:14*(N+1),1) = -40; args.ubx(2:14:14*(N+1),1) = 40;
        args.lbx(3:14:14*(N+1),1) = -10; args.ubx(3:14:14*(N+1),1) = 10;
        args.lbx(4:14:14*(N+1),1) = -45; args.ubx(4:14:14*(N+1),1) = 45; %rpy in degree
        args.lbx(5:14:14*(N+1),1) = -45; args.ubx(5:14:14*(N+1),1) = 45;
        args.lbx(6:14:14*(N+1),1) = -inf; args.ubx(6:14:14*(N+1),1) = inf;
        args.lbx(7:14:14*(N+1),1) = -180; args.ubx(7:14:14*(N+1),1) = 180; %rpy_dot in degree
        args.lbx(8:14:14*(N+1),1) = -180; args.ubx(6:14:14*(N+1),1) = 180;
        args.lbx(9:14:14*(N+1),1) = -180; args.ubx(9:14:14*(N+1),1) = 180;
        args.lbx(10:14:14*(N+1),1) = -7; args.ubx(10:14:14*(N+1),1) = 90; %servo angle in degree
        args.lbx(11:14:14*(N+1),1) = -180; args.ubx(11:14:14*(N+1),1) = 180; %rpy_dot last set points in degree
        args.lbx(12:14:14*(N+1),1) = -180; args.ubx(12:14:14*(N+1),1) = 180;
        args.lbx(13:14:14*(N+1),1) = -180; args.ubx(13:14:14*(N+1),1) = 180;
        args.lbx(14:14:14*(N+1),1) = 0; args.ubx(14:14:14*(N+1),1) = 105; %thrust in N calculated from max pitch angle and weight
        args.lbx(14*(N+1)+1:5:14*(N+1)+5*N,1) = 0; args.ubx(14*(N+1)+1:5:14*(N+1)+5*N,1) = 105; %thrust in N calculated from max pitch angle and weight
        args.lbx(14*(N+1)+2:5:14*(N+1)+5*N,1) = -45; args.ubx(14*(N+1)+2:5:14*(N+1)+5*N,1) = 45; %servo angle speed in degree/s
        args.lbx(14*(N+1)+3:5:14*(N+1)+5*N,1) = -60; args.ubx(14*(N+1)+3:5:14*(N+1)+5*N,1) = 60; %rpy_MPC in degree
        args.lbx(14*(N+1)+4:5:14*(N+1)+5*N,1) = -60; args.ubx(14*(N+1)+4:5:14*(N+1)+5*N,1) = 60;
        args.lbx(14*(N+1)+5:5:14*(N+1)+5*N,1) = -90; args.ubx(14*(N+1)+5:5:14*(N+1)+5*N,1) = 90;



        %run mpc



        x0 = [Init_Vel ; Init_rpy ; Init_rpy_dot;Init_Tilt;Init_last_rpy_dot;Init_last_thrust]; % initial condition.
        x_ref= [vel_des; [0;0;0] ; [0;0;0];0;[0;0;0];70];
        u0 = zeros(N,5); %  control inputs for each robot

        X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables


        args.p = [x0;x_ref]; % set the values of the parameters vector
            % initial value of the optimization variables
        args.x0 = [reshape(X0',14*(N+1),1);reshape(u0',5*N,1)];
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

        u = reshape(full(sol.x(14*(N+1)+1:end))',5,N)'; % get controls only from the solution
        
        R_BR =[cos(Init_Tilt)  0   sin(Init_Tilt);
            0           1   0         ;
            sin(Init_Tilt)   0   -cos(Init_Tilt)];

        angles = [Init_rpy(3) Init_rpy(2) Init_rpy(1)];

        RNI = zeros(3, 3);
        cang = cos(angles);
        sang = sin(angles);


        RNI = [cang(2).*cang(1)                            , cang(2).*sang(1)                            , -sang(2);
            sang(3).*sang(2).*cang(1) - cang(3).*sang(1), sang(3).*sang(2).*sang(1) + cang(3).*cang(1), sang(3).*cang(2);
            cang(3).*sang(2).*cang(1) + sang(3).*sang(1), cang(3).*sang(2).*sang(1) - sang(3).*cang(1), cang(3).*cang(2)];

        rbi = RNI'; % Body to inertial matrix

        Body = R_BR * [0;0;u(1,1)'];
        Iner = rbi*Body;


        thrust_des = obj.LimitOutput(Iner);


        rpy_des = u(1,3:5)';


        end

        function des_vel = SetVelDes(obj,mult,waypoint_des,time)
            obj.SetPID(10,0,3);
      

            dt = time - obj.LastTime;
            %disp(vel_des)
            % Calculate the error
            pos_err = waypoint_des - mult.State.Position;


            vel_err = 0 - mult.State.Velocity;

            % Update the error integral
            obj.ErrorIntegral = obj.ErrorIntegral + pos_err * dt;

            % Calculate the PID result
            des_vel =  obj.P * pos_err + ...
                obj.D * vel_err + obj.I * obj.ErrorIntegral;

            des_vel = obj.LimitVel(des_vel);
        end





        function force = CalcAeroForce(obj,V)
            %         V = V';
            %         a = atand(V(3)/V(1)) + obj.aoi;
        a = i/2*log((i+V(3)/V(1)/(i-V(3)/V(1)))) + obj.aoi;
        b = 0;
%         R_BW = [cosd(a) * cosd(b),    sind(b),     sind(a)*cosd(b);
%                     -sind(b) * cosd(a),   cosd(b),     -sind(a)*sind(b);
%                     -sind(a)         ,   0  ,        cosd(a)];
        R_BW = [cos(a) * cos(b),    sin(b),     sin(a)*cos(b);
                    -sin(b) * cos(a),   cos(b),     -sin(a)*sin(b);
                    -sin(a)         ,   0  ,        cos(a)];
        
        R_WB = R_BW.';
%       disp((V' * V))
        q_bar = (V' * V) * physics.AirDensity / 2;
        c_y = 0; 
        c_z = 0.35 + 0.11 * a;
        c_d = 0.01 + 0.2 * a * a;
        drag = q_bar * obj.WingSurfaceArea * c_d;
        lateral = q_bar * obj.WingSurfaceArea * c_y;
        lift = q_bar * obj.WingSurfaceArea * c_z;
%       disp(size(q_bar))
        force = R_WB * [-drag; lateral;-lift];
    end

    function accel = CalcAccel(obj,F_a,T,rpy,tilt)
        
       
        T = [0;0;T];
        
        R_BR = [cos(tilt)  0   sin(tilt);
                0           1   0         ;
                sin(tilt)   0   -cos(tilt)];
        T = R_BR * T;

%       disp(F_a)
        F_tot = F_a + T;
        angles = [rpy(3) rpy(2) rpy(1)];
        
        RNI = zeros(3, 3);
        cang = cos(angles);
        sang = sin(angles);

        
        RNI = [cang(2).*cang(1)                            , cang(2).*sang(1)                            , -sang(2);
               sang(3).*sang(2).*cang(1) - cang(3).*sang(1), sang(3).*sang(2).*sang(1) + cang(3).*cang(1), sang(3).*cang(2);
               cang(3).*sang(2).*cang(1) + sang(3).*sang(1), cang(3).*sang(2).*sang(1) - sang(3).*cang(1), cang(3).*cang(2)];

        rbi = RNI'; % Body to inertial matrix
        F_Inertial = rbi * F_tot;

        accel = F_Inertial / obj.Mass;
    end
    end
end


%% Helper functions


