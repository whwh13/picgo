classdef runs
    properties
        number
        s=snapshot
        time_vector;
        chem;
    end
    methods
        function obj=runs(N,varargin)
            if(nargin==2)
            obj(N)=run;
            M=varargin{2};
            for idx=1:N
                obj(idx).s(M)=snapshot;
            end
            elseif(nargin==1)
            for idx=1:N
                obj(idx).number=idx;
            end 
            end
        end
        function obj=appendsnapshot(obj,s1)
            if(nargin>0)
                if(~isempty(obj.s(numel(obj.s)).serial))
                obj.s(numel(obj.s)+1)=s1;
                else
                obj.s(numel(obj.s))=s1;
                end
            end
        end
        function s1=getsnapshot(obj,serial)
            s1=obj.s(find([obj.s(:).serial]==serial));
        end
    end
end
