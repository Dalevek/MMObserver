classdef identificationlibrary < handle
    %IDENTIFICATIONLIBRARY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %%chyba nie potrzebne...
        n,m,forgetting, theta;
        K,P,In,Im,setInitalTheta;
        p0;
        
        %setInitialTheta;
        error;
    end
    
    methods
        function obj = identificationlibrary(n,m,p0,forgetting,theta0)
            %default value
            %obj.p0 = 1000.0;
            %obj.forgetting = 1.0;
            obj.p0 = p0;
            obj.forgetting = forgetting;
            
            obj.n = n;
            obj.m = m;
            obj.forgetting = forgetting;
            obj.theta = zeros(obj.n,1);
            obj.K = zeros(obj.n,obj.m);
            obj.P = p0*eye(obj.n,obj.n);
            obj.In = eye(obj.n,obj.n);
            obj.Im = eye(obj.m,obj.m);
            obj.theta = theta0;
            %obj.setInitialTheta(theta0);
            
        end
        
        function setInitialTheta(obj,theta0)
            %theta0 = squeeze(struct2table(theta0,'AsArray',true));
            theta0 = squeeze(theta0);
            if obj.n == 1
                obj.theta(0) = theta0;
            elseif length(theta0) == obj.n
                for i = 1:obj.n
                    obj.theta(i) = theta0(i);
                end
            end   
        end
        
        function update = update(obj,y,phi)
            obj.error = y - phi*obj.theta;
            %update gain
            tmp = obj.forgetting*obj.Im + phi*obj.P*phi.';
            obj.K = obj.P*phi.'/tmp; %inverse - znacznik "/"
            %update P
            obj.P = ((obj.In - obj.K*phi)*obj.P)/obj.forgetting;
            %update parameters
            obj.theta = obj.theta + obj.K*obj.error;
            %%%%%%%%%%%%%%%% TESTOWO - return
            update = obj.theta;
            return
        end
        
        function output = output(obj,phi)
            output = phi*obj.theta;
            return
        end
        
        function restart(obj,p0)
            obj.P = p0*eye(obj.n,obj.n);
        end
    end
    
end

