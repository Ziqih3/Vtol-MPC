classdef attitude_controller_vtol < pid_controller

    properties

        MC_attitude_controller attitude_controller
        FW_attitude_controller attitude_controller_vtol
        RateLimits = [70; 70; 30]; % in deg/s
        OutputMax = [8; 8; 8]; % in rad/s^2

        %%%%%%%%%%%%%%%%%%%%%
        blending_speed = 16;%m/s
        Transition_air_speed = 20;
        weight_MC = 1;
        weight_FW = 0;
        vtol_mode = 'MC mode';


    end

    methods
        function euler_accel = CalculateControlCommand(obj,mult, rpy_des, rpy_dot_des, eul_acc_des, time)



            obj.MC_attitude_controller = attitude_controller;
            obj.FW_attitude_controller = attitude_controller_fw;


            obj.MC_attitude_controller.SetPID(obj.P,obj.I,obj.D);
            obj.FW_attitude_controller.SetPID(obj.P,obj.I,obj.D);

            
%             disp('p')
%             disp(obj.P)
%             disp('I')
%             disp(obj.I)
%             disp('D')
%             disp(obj.D)



            if isempty(rpy_dot_des)
                rpy_dot_des = zeros(3, 1);
            end
            if isempty(eul_acc_des)
                eul_acc_des = zeros(3, 1);
            end
            dt = time - obj.LastTime;

            %%% current transition airspeed is 20 m/s
            obj.weight_MC = 1;


            %if is vtol_front_transition
            if sqrt(mult.State.AirVelocity(1)^2 + mult.State.AirVelocity(2)^2 + mult.State.AirVelocity(3)^2) >= obj.blending_speed;

                if obj.vtol_mode == 'Phase_1'
                    obj.weight_MC = min(((sqrt(mult.State.AirVelocity(1)^2 + mult.State.AirVelocity(2)^2 + mult.State.AirVelocity(3)^2)-obj.blending_speed)/(obj.Transition_air_speed-obj.blending_speed)),1);
                    obj.weight_FW = 1 - obj.weight_MC;

                elseif obj.vtol_mode == 'MC_mode'
                    obj.weight_MC = 1;
                    obj.weight_FW = 0;
                else
                    obj.weight_MC = 0;
                    obj.weight_FW = 1;
                end

            end
             %disp(obj.weight_MC)
             %disp(obj.weight_FW)
             %disp(sqrt(mult.State.AirVelocity(1)^2 + mult.State.AirVelocity(2)^2 + mult.State.AirVelocity(3)^2))
             %disp(obj.vtol_mode)

            % MC
            %%%%%%%%%%% keep thrust as the MC current thrust, but cannot be lower
            % than 0.25 maximum (position controller all off, give a hard
            % number)

            if mult.Servos{1}.CurrentAngle == 0
                obj.vtol_mode = 'MC_mode';
                %disp(obj.vtol_mode)


                %euler_accel = obj.MC_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des, dt);

                euler_accel = obj.weight_MC * obj.MC_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des, dt) + obj.weight_FW * obj.FW_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des,dt);

                % Phase_1_needs_blending
            elseif mult.Servos{1}.CurrentAngle > 0 & mult.Servos{1}.CurrentAngle <=27
                obj.vtol_mode = 'Phase_1';
                %disp(obj.vtol_mode)



                %euler_accel = obj.MC_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des, dt);

                euler_accel = obj.weight_MC * obj.MC_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des, dt) + obj.weight_FW * obj.FW_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des,dt);

                % Phase_2
            elseif mult.Servos{1}.CurrentAngle > 27 & mult.Servos{1}.CurrentAngle < 90

                obj.vtol_mode = 'Phase_2';
                %disp(obj.vtol_mode)

                %euler_accel = obj.FW_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des,dt);

                euler_accel = obj.weight_MC * obj.MC_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des, dt) + obj.weight_FW * obj.FW_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des,dt);


            elseif mult.Servos{1}.CurrentAngle == 90


                % Phase_FW
                %%%%%%%%%%% keep thrust as 100% FW mode (position controller
                %%%%%%%%%%% on)

                obj.vtol_mode = 'FW_mode';
                %disp(obj.vtol_mode)

                %euler_accel = obj.FW_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des,dt);

                euler_accel = obj.weight_MC * obj.MC_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des, dt) + obj.weight_FW * obj.FW_attitude_controller.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des,dt);




            end

            obj.LastTime = time;
        end

    end
end
