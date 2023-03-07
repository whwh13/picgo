function Medyan2VMD(path,outputfile)
%% Converts snapshot.traj to VMD readable pdb file.
% Medyan2VMD(path,outputfile)
% path - path to snapshot.traj.
% outuputfile - desired outputfile name.

RoC_parted10={};
RoC10={};
Length10={};
time_vector10={};
polarity_profile={};
orientation_profile={};
N = 1;
r=runs(N);

for i=1:N
	i
    r(i)=readsnapshot([path,'/snapshot.traj'],r(i));
    time_vector10=[time_vector10;{r(i).time_vector}];
end
VMDstylesnap_polym(r,outputfile);
save([outputfile,'.mat']);
end
