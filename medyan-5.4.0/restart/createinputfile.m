function createinputfile(Ifile,Ofile,varargin)  
%% function takes in input file, processes it to Medyan readable format.
% createinputfile(Ifile,Ofile,L,Lf,M,Mf,B,Bf)
% Ifile - input file eg. 'restartinput.txt'
% Ofile - output file eg. 'restartoutput.txt'
% L - cell list of SPECIESLINKER names. eg. {'LA1','LA2'}
% Lf -list of filament IDs corresponding to L. eg. [0,1]
% M - cell list of SPECIESMOTOR names. eg. {'MOA1','MOA2'}
% Mf -list of filament IDs corresponding to M. eg. [0,0]
% B - cell list of SPECIESBRANCHER names. eg. {'BA1,'BA2'}
% Bf -list of filament IDs corresponding to B eg. [1,1]
if(length(varargin)>=2)
    L=varargin{1};
    Lf=cell2mat(varargin(2));
end
if(length(varargin)>=4)
    M=varargin{3};
    Mf=cell2mat(varargin(4));
end
if(length(varargin)==6)
    B=varargin{5};
    Bf=cell2mat(varargin(6));
end
    
f1=fopen(Ifile,'r');
f2=fopen(Ofile,'w');
while(~feof(f1))
    line=fgetl(f1);
    if(numel(line)>5)
        if(strcmp(line(1:8),'FILAMENT'))
            ftype=str2num(line(9:end));
            ftype=ftype(2);
            fprintf(f2,'%s\n',['FILAMENT ',num2str(ftype),' ',fgetl(f1)]);
        elseif(strcmp(line(1:6),'LINKER'))
            ltype=str2num(line(7:end));
            ltype=ltype(2)+1;
            fprintf(f2,'%s\n',[L{ltype},' ',num2str(Lf(ltype)),' ',fgetl(f1)]);
        elseif(strcmp(line(1:5),'MOTOR'))
            mtype=str2num(line(6:end));
            mtype=mtype(2)+1;
            fprintf(f2,'%s\n',[M{mtype},' ',num2str(Mf(mtype)),' ',fgetl(f1)]);
        elseif(strcmp(line(1:5),'BRANCHER'))
            btype=str2num(line(6:end));
            btype=btype(2);
            fprintf(f2,'%s\n',[B{btype},' ',num2str(Bf(btype)),' ',fgetl(f1)]);
        end
    end
end
fclose('all');
end
