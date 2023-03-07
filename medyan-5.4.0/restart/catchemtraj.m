function catchemtraj(file1,file2,catfile,frame,time_limit)
file1=fopen([file1,'/chemistry.traj'],'r');
file2=fopen([file2,'/chemistry.traj'],'r');
cfile=fopen([catfile,'/chemistry.traj'],'a');
line=fgetl(file1);
curr_time=str2num(line);
curr_time=curr_time(2);
fprintf(cfile,'%s\n',line);
stop=0;
while(~feof(file1)&& ~stop)
    line=fgetl(file1);
    if(isempty(line)&&~feof(file1))
      line=fgetl(file1);
        curr_time=str2num(line);
        curr_time=curr_time(2);
        if(curr_time<=frame)
        fprintf(cfile,'\n');
        end
    end
    if(curr_time<=frame)
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
     if(isempty(line) && isstr(line));
        if(f)
            fprintf(cfile,'%s',line);
        end
        f=1;
            linex=fgetl(file2);
	      if(isstr(linex))
	    linex=str2num(linex);
            t=linex(2)+frame;
            if(t<=time_limit+5)
            fprintf(cfile,'\n%i %7.4f\n',linex(1),linex(2)+frame);
            else
                stop=1;
                break;
            end
	    else
	      stop=1;
            break;
	    end
    elseif(f && ~stop)
    fprintf(cfile,'%s\n',line);
    elseif(stop && ~f) break;
    end                   
end    
fclose('all');
end
