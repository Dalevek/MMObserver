classdef MMObserver < handle
 %    '''
 %     Multi Observer for linear time invariant system
 %    '''
    
    properties
        Tp,N,M,MM,A,B,C,L;
        mapping,observerType;
        xemCurr = {}, xemNext = {}, yem, ksim = {};
        memory = {},finiteMemoryDelay;
        xe0,alfa;
        %integralChainInitialization() param.
        Az,Bz;
        izmCurr , izmNext ;
                   
           
        %createMatrixT() param.
        T;
        
        
        idyRLS;
    end
    
    methods
        function obj = MMObserver(N, M, Tp, A, B, C, L, mapping, observerType, RLSGain, RLSForgetting, xe0)
%         mapping = 'integral-finite'
%         observerType = 'RLS'
%         RLSGain = 1000.0
%         RLSForgetting = 1.0
%         xe0 = []
        
            obj.Tp = Tp;% time sample
            obj.N = N;% number of estimated state variables
            obj.M = M; % number of estimated output variables
            obj.MM = obj.N + 1; % number of multi observers
            obj.A = A;% observer matrixes
            obj.B = B;
            obj.C = C;
            obj.L = L;
            obj.mapping = mapping;% type of mapping between output and xi
            obj.observerType = observerType; % method of update of alpha weights
            if obj.mapping == 'integral-finite' %size of memory
                obj.finiteMemoryDelay = 10;
            end

            obj.xemCurr = []; %estimated state for m-observer
            obj.xemNext = [];
            obj.yem = zeros(obj.M,1);%estimated output for m-observer
            obj.ksim = [];%xi signal for m-observer
            
            if obj.mapping == 'integral-finite'
               obj.memory = []; %memory for xi signal for m-observer   
            end

            for i = 1:obj.MM %[Poprawiona pêtla] initialized signals with zeros for m-observers
                %obj.xemCurr.append(zeros(obj.N,1));
                %obj.xemNext.append(zeros(obj.N,1));
                %obj.ksim.append(zeros(obj.N,1)); <------ z Pythona (brak
                %metody append)
                obj.xemCurr(:,:,i) = zeros(obj.N,1);
                obj.xemNext(:,:,i) = zeros(obj.N,1);
                obj.ksim(:,:,i) = zeros(obj.N,1);
                if mapping == 'integral-finite'                        
                    %obj.memory.append(zeros(obj.N,obj.finiteMemoryDelay));
                    obj.memory(:,:,i) = zeros(obj.N,obj.finiteMemoryDelay);
                end
            end

            obj.integralChainInitialization();

            obj.xe0 = xe0; %initial points for m-observers
            obj.setXemPoints(obj.xe0)  ;                                                               

            obj.alfa = ((1.0/obj.MM)*ones(obj.MM,1));%weights for m-observers
            if obj.observerType == 'RLS'     
                obj.idyRLS = identificationlibrary(obj.MM-1,obj.MM-1, RLSGain, RLSForgetting, obj.alfa(1:obj.MM-1,1));
            end
        end
      
        function [izmCurr, izmNext] = integralChainInitialization(obj)
            
            Nz = obj.M*(obj.N-1); %integral of output error for m-observer
            Nones = obj.M*(obj.N-2);
            obj.Az = diag(ones(Nones), obj.M);
            obj.Bz = zeros(Nz,obj.M);
            obj.Bz(Nones:Nz,:) = diag(obj.M); %????
            obj.createMatrixT();
            %obj.izmCurr = [];
            %obj.izmNext = [];   
            for i = 1:obj.MM
                %obj.izmCurr.append(zeros(Nz,1));
                %obj.izmNext.append(zeros(Nz,1));
                obj.izmCurr(:,:,i) = zeros(Nz,1);
                obj.izmNext(:,:,i) = zeros(Nz,1);
            end
            return;
        end
            
        function createMatrixT(obj)
        
            NM = obj.N*obj.M;
            obj.T = zeros(NM, NM);
            for j = 1:NM
                for k = 1:NM
                    js = (j-1)*obj.M+1;
                    je = j*obj.M;
                    ks = (k-1)*obj.M+1;
                    ke = k*obj.M;               
                    if j == k                    
                        obj.T(js:je,ks:ke) = eye(obj.M);
                    elseif j > k
                        p = j - k;
                        obj.T(js:je,ks:ke) = obj.C*(obj.A^(p-1))*obj.L; 
                    end
                end
            end
         end
        
        function setXemPoints(obj,xe0)
            if isempty(xe0)%initial points for m-observers
                %default vectors for 3D: [1 0 0]^T, [0 1 0]^T, [0 0 1]^T  [0 0 0]^T
                for m = 1:obj.MM-1
                    obj.xemNext(m,1) = 1.0;
                end
            else         
                for m = 1:obj.MM-1%user defined (initial condition multiplied by 1.0/alfa_k)
                    obj.xemNext(m,1) = obj.MM*xe0(m);
                end
            end
        end
        
        function singleObserver(obj, m, u, etam)                     
%             dxem = obj.A*obj.xemCurr(m,:,0) + obj.B*u - obj.L*(etam(:,0));%derivative of estimated state        
%             obj.xemNext(m,:,0) = obj.xemCurr(m,:,0) + obj.Tp*dxem; %estimated state
            dxem = obj.A*obj.xemCurr(:,1,m) + obj.B*u - obj.L*(etam(:,1));%derivative of estimated state        
            obj.xemNext(:,1,m) = obj.xemCurr(:,1,m) + obj.Tp*dxem; %estimated state
        end
        
        function calculateMapping(obj, m, etam)
%             tmp = zeros(obj.N*obj.M,1); %mapping between output and xi for m-observer
%             tmp(0:(obj.N-1)*obj.M,0) = obj.izmCurr(m,:,0);
%             tmp((obj.N-1)*obj.M:obj.N*obj.M) = etam(:,0);
%             obj.ksim(m,:,0) = obj.T*tmp     ;           
%             if obj.mapping == 'integral-finite'
%                 ksimd = ctls.delay(obj.ksim(m,:,0), obj.memory(m), obj.finiteMemoryDelay); %calculate integral of xi between from t-td to t
%                 obj.ksim(m,:,0) = obj.ksim(m,:,0) - ksimd;
%             end
            tmp = zeros(obj.N*obj.M,1); %mapping between output and xi for m-observer
            tmp(1:(obj.N-1)*(obj.M),1) = obj.izmCurr(:,:,m);
            tmp((obj.N-1)*obj.M:obj.N*obj.M) = etam(:,1);
            obj.ksim(:,:,m) = obj.T*tmp     ;           
            if obj.mapping == 'integral-finite'
                ksimd = delay(obj.ksim(:,1,m), obj.memory(:,:,m), obj.finiteMemoryDelay); %calculate integral of xi between from t-td to t
                obj.ksim(:,1,m) = obj.ksim(:,1,m) - ksimd;
            end
        end
        
        function weightsRLS(obj)
            idyY = -obj.ksim(:,1,obj.MM);%left side of regressor equation

            M = zeros(obj.MM-1,obj.MM-1);%right side of regressor equation 
            for m = 1:obj.MM-1
                M(:,m) = obj.ksim(:,1,m) - obj.ksim(:,1,obj.MM);
            end

            tmp = obj.idyRLS.update(idyY, M);%weights update
            for m = 1:(obj.MM-1) 
                obj.alfa(m,1) = tmp(m);
            end
            obj.alfa(obj.MM,1) = 1.0 - sum(obj.alfa(1:obj.MM,1));%.sum();   
        end
        
        function xe = observe(obj, u, y)
            xe = zeros(obj.N,1);  %estimated state       

            for m = 1:(obj.MM) %calculate for m-observer

                %obj.xemCurr(m,:,1) = obj.xemNext(m,:,1);%change
                %obj.izmCurr(m,:,1) = obj.izmNext(m,:,1);
                obj.xemCurr(:,:,m) = obj.xemNext(:,:,m);%change
                obj.izmCurr(:,:,m) = obj.izmNext(:,:,m);
                %obj.yem(:,1) = obj.C*obj.xemCurr(m,:,1);
                obj.yem(:,1) = obj.C*obj.xemCurr(:,:,m);%observer estimated output

                etam = obj.yem(:,1) - y; %observer estimated output error

%                 obj.izmNext(m,:,1)  = obj.izmCurr(m,:,1)  + obj.Tp*(obj.Az*obj.izmCurr(m,:,1) + obj.Bz*etam(:,1));%integral of error 
%                 obj.singleObserver(m, u, etam);            
%                 obj.calculateMapping(m, etam);
% 
%                 xe(:,1) = xe(:,1) + np.squeeze(np.asarray(obj.alfa(m,1)))*obj.xemCurr(m,:,1);%virtual single model
                obj.izmNext(:,:,m)  = obj.izmCurr(:,1,m)  + obj.Tp*(obj.Az*obj.izmCurr(:,1,m) + obj.Bz*etam(:,1));%integral of error 
                obj.singleObserver(m, u, etam);            
                obj.calculateMapping(m, etam);

                xe(:,1) = xe(:,1) + squeeze(obj.alfa(m,1))*obj.xemCurr(:,:,m);%virtual single model
            end
        
        % if obj.observerType == 'none':
        % no weights estimation (equal to single mode)        
            if obj.observerType == 'RLS'

                obj.weightsRLS();%weights estimation based on RLS

            end
            return; %xe(:,0);%return estimated state
        end
        
        function restart(obj, p0)
            obj.setXemPoints(obj.xe0);
            if obj.observerType == 'RLS'   
                obj.idyRLS.restart(p0);
            end
        end
    
 end

end      