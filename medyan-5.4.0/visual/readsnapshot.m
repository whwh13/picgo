function r=readsnapshot(filename,r)
f1=fopen(filename,'r');
 time_step=str2num(fgetl(f1));
       r.time_vector=time_step(2);
       line=fgetl(f1);
       Coord_cell={};
       Coord_myosin_cell=[];
       myosin_id_cell={};
       Coord_linker_cell=[];
       linker_id_cell={};
       myosin_id=[];
       linker_id=[];
       s_no=1;
       Coord=[];
       filament_numbeads=[];
while(~feof(f1))
%   while(~feof(f1))
    if(length(line)==0)
       time_step=str2num(fgetl(f1));
       r.time_vector=[r.time_vector,time_step(2)];
       s_=snapshot;
       s_.serial=s_no; % frame number;
       s_no=s_no+1;
%       [coord_cell1,coord_cell2]=cluster_filaments(Coord_cell,Coord,filament_numbeads);
       s_.f.coord_cell1=Coord_cell;
       s_.f.coord_cell2={};
       s_.m.coord_cell=Coord_myosin_cell;
       s_.m.id=myosin_id;
       s_.l.id=linker_id;
       s_.l.coord_cell=Coord_linker_cell;
       r=appendsnapshot(r,s_);
       clear s_;
       Coord=[];
       Coord_cell={};
       Coord_myosin_cell=[];
       Coord_linker_cell=[];
       myosin_id=[];
       linker_id=[];
       filament_numbeads=[];
    elseif(strcmp(line(1),'F')==1)
        dummy=str2num(fgetl(f1));
       Coord_cell=[Coord_cell;{dummy}];
       Coord=[Coord;reshape(dummy,3,[])'];
       filament_numbeads=[filament_numbeads,size(Coord,1)];
       clear dummy;
    elseif(strcmp(line(1),'M')==1)
        Coord_myosin_cell=[Coord_myosin_cell;str2num(fgetl(f1))];
       dummy=strsplit(line,' ');
       myosin_id=[myosin_id,str2double(dummy(2))];       
       elseif(strcmp(line(1),'L')==1)
        Coord_linker_cell=[Coord_linker_cell;str2num(fgetl(f1))];
       dummy=strsplit(line,' ');
     linker_id=[linker_id,str2double(dummy(2))];
   end
   line=fgetl(f1);
end
       if(feof(f1))
       s_=snapshot;
       s_.serial=s_no; % frame number;
       s_no=s_no+1;
%       [coord_cell1,coord_cell2]=cluster_filaments(Coord_cell,Coord,filament_numbeads);
       s_.f.coord_cell1=Coord_cell;
       s_.f.coord_cell2={};
       s_.m.coord_cell=Coord_myosin_cell;
       s_.m.id=myosin_id;
       s_.l.id=linker_id;
       s_.l.coord_cell=Coord_linker_cell;
       
       r=appendsnapshot(r,s_);
       size(r.s)
       else
           r.time_vector=[r.time_vector(1:end-1)];
       end
       fclose(f1);
       clear s_ Coord_cell Coord_myosin_cell myosin_id filament_numbeads;
end
