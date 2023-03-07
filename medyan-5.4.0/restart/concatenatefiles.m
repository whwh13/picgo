function concatenatefiles(frame1,time_limit,file1dir,file2dir,outputdir)
%% Concatenates snapshot, chemistry trajectory files of restarted simulations.
%% concatenatefiles(frame1,time_limit,file1dir,file2dir,outputdir)
%frame1 - is the last time step to consider in file1. Any frame with time stamp lesser
% or equal will be written. eg. 1502.37
%time_limit - desired time duration of simulation eg. 2000.
%file1dir - full path leading to directory of original simulation.
%file2dir - full path leading to directory of restarted simulation
%outputdir - output directory (add full path). eg. 'restartdir'
mkdir(outputdir);
catfile=outputdir;
list=dir(file1dir);
for i=1:size(list,1)
    if(numel(list(i).name)>5)
    if(strcmp('traj',list(i).name(end-3:end)))
        if(strcmp(list(i).name(1:end-5),'chemistry'))
            catchemtraj(file1dir,file2dir,catfile,frame1,time_limit);
        elseif(strcmp(list(i).name(1:end-5),'snapshot'))
            catsnaptraj(file1dir,file2dir,catfile,frame1,time_limit);
        elseif(strcmp(list(i).name(1:end-5),'walltensions'))
            catwalltraj(file1dir,file2dir,catfile,frame1,time_limit);
        end
    end
    end
end
end
