classdef controller < handle
    properties
        ControlAllocation control_allocation_vtol
        AttitudeController attitude_controller_vtol
        PositionController position_controller_fw_MPC_new
        HMFController hmf_controller
    end
    
    methods
        function obj = controller(mult)
            obj.ControlAllocation = control_allocation_vtol(mult);
            obj.AttitudeController = attitude_controller_vtol;
            obj.PositionController = position_controller_fw_MPC_new;
            obj.HMFController = hmf_controller;
        end
        
        function [rotor_speeds_squared, deflections,tilt,saturated] = ControlAcceleration(obj, mult, lin_acc_des, euler_acc_des)
            [rotor_speeds_squared, deflections,tilt,saturated] = obj.ControlAllocation.CalcActuators(mult, lin_acc_des, euler_acc_des);
        end
        
        function euler_accel = ControlAttitude(obj, mult, rpy_des, rpy_dot_des, eul_acc_des, dt)
            euler_accel = obj.AttitudeController.CalculateControlCommand(mult, rpy_des, rpy_dot_des, eul_acc_des, dt);
        end

        function [lin_accel, rpy_des] = ControlPosition(obj, mult, pos_des, yaw_des, vel_des, acc_des, dt)
            vel_des = obj.PositionController.SetVelDes(mult,pos_des,dt);
            disp('vel_des')
            disp(vel_des)
            [lin_accel, rpy_des] = obj.PositionController.CalculateControlCommand(mult, pos_des, vel_des, yaw_des,acc_des, dt);
        end
        
        function [lin_accel, rpy_des] = ControlMotionAndForce(obj, mult, force_des, pos_des, yaw_des, vel_des, acc_des, ...
                contact_normal, vel_mat, force_constraint, dt)
            [lin_accel, rpy_des] = obj.HMFController.ControlMotionAndForce(mult, ...
                force_des, pos_des, yaw_des, vel_des, acc_des, contact_normal, vel_mat, force_constraint, dt);
        end
        
        function Reset(obj)
            obj.AttitudeController.Reset();
            obj.PositionController.Reset();
        end
        
        function SetAttitudeStrategy(obj, attitude_strategy)
            obj.PositionController.SetAttitudeStrategy(attitude_strategy);
            obj.HMFController.SetAttitudeStrategy(attitude_strategy);
        end
    end
end
