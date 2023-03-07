function catsnaptraj(file1,file2,catfile,frame,time_limit)
file1=fopen([file1,'/snapshot.traj'],'r');
file2=fopen([file2,'/snapshot.traj'],'r');
cfile=fopen([catfile,'/snapshot.traj'],'a');
line=fgetl(file1);
curr_time=str2num(line);
curr_time=curr_time(2);
fprintf(cfile,'%s\n',line);
stop=0;
while(~feof(file1)&& ~stop)
    line=fgetl(file1);
    if(isempty(line))
        line=fgetl(file1);
        curr_time=str2num(line);
        curr_time=curr_time(2);
        if(round(curr_time,2)<=frame)
        fprintf(cfile,'\n');
        end
    end
    if(round(curr_time,2)<=frame)
    fprintf(cfile,'%s\n',line);
    else
        stop=1;
    end
end
fclose(file1);
f=0;
t=time_limit*2;
stop=0;
while(~feof(file2) && ~stop)
    line=fgetl(file2);
    if(isempty(line));
        if(f)
            fprintf(cfile,'%s',line);
        end
        f=1;
        linex=fgetl(file2);
        if(isstr(linex))
            linex=str2num(linex);
            t=linex(2)+frame;
            if(t<=time_limit+5)
            fprintf(cfile,'\n%i %7.4f',linex(1),linex(2)+frame); 
            else
                stop=1;
            end
	    else
	      stop=1;
              break;
            end
    end
    if(f && ~stop)
    fprintf(cfile,'%s\n',line);    
    end                    
end    
fclose('all');
end
