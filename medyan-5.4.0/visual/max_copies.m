classdef max_copies
    properties
        f_max=0;
        l_max=0;
        m_max=0;
    end
    methods
         function obj=max_copies(nf,nl,nm) %number of fil, linker, motor types
             obj.f_max=nf;
             obj.l_max=nl;
             obj.m_max=nm;
         end
    end
end