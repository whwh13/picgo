classdef snapshot
    properties
        serial
        f=filaments;
       l=elements
        m=elements
    end
    methods
        function obj=snapshot(M)
            if(nargin>0)
                obj(M)=snapshot;
                for idx=1:M
                    obj(idx).id=idx;
                end
            end
        end
    end
end
