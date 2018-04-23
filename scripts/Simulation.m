
classdef Simulation    < handle
		properties
           name,Tp,Np,N,M,P;
           R,L,Kt,fd,J,mi;
           A,B,C;
           Ao,Bo,Co,Lo;
           %yrd
           x0_1, x0_2, x0_3;
        end
        
        methods
            
        
    %logger = logging.getLogger('Simulation') 
    
        function obj = Simulation(name, x0_1, x0_2, x0_3)
            %default parameters
            obj.x0_1 = x0_1;
            obj.x0_2 = x0_2;
            obj.x0_3 = x0_3;
            
            
            obj.name = name;
            obj.Tp = 0.0001;
            obj.Np = 6*450;
            obj.N = 3; %# state size
            obj.M = 1;%# output number
            obj.P = 1; %# input number
            %# motor parameters
            obj.R  = 3.2         ;%# Ohm
            obj.L  = 0.0086      ;%# mH
            obj.Kt = 0.0319      ;%# Nm/A       
            obj.fd = 0.00012    ; %# Nms/rad
            obj.J  = 30*10.0^-6; %# kgm2
            obj.mi = -0.06;     %# Nm/As
            %# DC motor 
            obj.A = zeros(obj.N-1,obj.N-1);
            obj.A(1,1) = - obj.fd/obj.J   ;     
            obj.A(1,2) =   obj.Kt/obj.J    ;    
            obj.A(2,1) = - obj.Kt/obj.L    ;    
            obj.A(2,2) = -  obj.R/obj.L    ;    
            obj.B = zeros(obj.N-1,obj.N-1) ;     
            obj.B(1,1) = - 1.0/obj.J;
            obj.B(2,2) =   1.0/obj.L;
            obj.C = zeros(1,obj.N-1);
            obj.C(1,1) = 0.0;
            obj.C(1,2) = 1.0; 
            %# observer definition
            obj.Ao = zeros(obj.N,obj.N);
            obj.Ao(1,1) = - obj.fd/obj.J;        
            obj.Ao(1,2) =   obj.Kt/obj.J;
            obj.Ao(1,3) = -     1.0/obj.J;
            obj.Ao(2,1) = - obj.Kt/obj.L;        
            obj.Ao(2,2) = -  obj.R/obj.L;        
            obj.Bo = zeros(obj.N,1);
            obj.Bo(1,1) = 0.0;
            obj.Bo(2,1) = 1.0/obj.L;
            obj.Bo(3,1) = 0.0;
            obj.Co = zeros(1,obj.N);
            obj.Co(1,1) = 0.0;
            obj.Co(1,2) = 1.0;
            obj.Co(1,3) = 0.0;   
            %# observer gain
            obj.Lo = zeros(obj.N,1);
            obj.Lo(3,1) = -obj.mi; 
            %# check closed loop
            %self.logger.debug('Ao-LoCo:')
            %self.logger.debug(self.Ao-self.Lo*self.Co)
            %w, _ = la.eig(self.Ao-self.Lo*self.Co)
            %self.logger.debug('Poles of matrix Ao-LoCo:')
            %self.logger.debug(w)
            return
        end 
    end
end