function VMDstylesnap(r,outputfile)
%% Eventhough the MATLAB code that gives the class r cannot take more than one type of filament,...
%% linker or motor, we create a class to hold more than one type of those for future expansion.
f1=fopen([outputfile,'.pdb'],'w');
count_model=1;
speciesvec=max_copies(1,1,1);
for runs=1
    for snap=1:size(r(runs).s,2)
        dummy=r(runs).s(snap).f.coord_cell1;
%         dummy2=size(r(runs).s(snap).f.coord_cell1,1);
        dummy2=[];
        for it=1:size(dummy,1)
            dummy2=[dummy2,numel(dummy{it})/3];
        end
        if(size(speciesvec.f_max,2)<size(dummy2,2))
            speciesvec.f_max=padarray(speciesvec.f_max,[0 abs(numel(speciesvec.f_max)-numel(dummy2))],'post');
        elseif(size(speciesvec.f_max,2)>size(dummy2,2))
            dummy2=padarray(dummy2,[0 abs(numel(speciesvec.f_max)-numel(dummy2))],'post');
        end
        speciesvec.f_max=max(dummy2,speciesvec.f_max);
        dummy2=size(r(runs).s(snap).l.coord_cell,1);
        speciesvec.l_max(1)=max(dummy2,speciesvec.l_max(1));
        dummy=r(runs).s(snap).m.coord_cell;
        dummy2=size(r(runs).s(snap).m.coord_cell,1);
        speciesvec.m_max(1)=max(dummy2,speciesvec.m_max(1));
    end
    clear dummy dummy2 it;
    for snap=1:size(r(runs).s,2)
       fprintf(f1,'MODEL     %4i\n',count_model);%can only pring 10k frames.
       count_f=0;
       fils=r(runs).s(snap).f.coord_cell1;
       chain=['F'];
       chnum=1;
       for f=1:size(fils,1)
           A=reshape(fils{f},3,[])'./10;% A model that is 10 times smaller.
           b=0.0;
           for i=1:size(A,1)
             fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,A(i,1),A(i,2),A(i,3),b);   
             count_f=count_f+1;
           end
           for i=size(A,1)+1:speciesvec.f_max(f)
              fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,A(end,1),A(end,2),A(end,3),b);   
              count_f=count_f+1;
           end
           count_f=count_f+1;
       end
       for f=size(fils,1)+1:size(speciesvec.f_max,2)
           b=0.0;
           for i=1:speciesvec.f_max(f)
              fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,0.0,0.0,0.0,b);   
              count_f=count_f+1;
           end
              count_f=count_f+1;
       end
       fprintf(f1,'TER   %5i      ARG %s%4i\n',count_f,chain(chnum),count_f);
       %% LINKER
       count_f=0;
       link=r(runs).s(snap).l.coord_cell./10;
       chain=['l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'];
       chnum=1;
       for i=1:speciesvec.l_max(1)
           if(count_f>=9998)
               fprintf(f1,'TER   %5i      ARG %s%4i\n',count_f,chain(chnum),count_f);
               count_f=0;
               chnum=chnum+1;
           end
           if(i<=size(link,1))
               fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,link(i,1),link(i,2),link(i,3),b);   
               count_f=count_f+1;
               fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,link(i,4),link(i,5),link(i,6),b);   
               count_f=count_f+2;
           else
               fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,0.0,0.0,0.0,b);   
               count_f=count_f+1;
               fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,0.0,0.0,0.0,b);   
               count_f=count_f+2;
           end
       end
       if(count_f<9998)
       fprintf(f1,'TER   %5i      ARG %s%4i\n',count_f,chain(chnum),count_f);
       end
       clear link;
       %% MOTOR
       count_f=0;
       motor=r(runs).s(snap).m.coord_cell./10;
       chain=['L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'];
       chnum=1;
       for i=1:speciesvec.m_max(1)
            if(count_f>=9998)
               fprintf(f1,'TER   %5i      ARG %s%4i\n',count_f,count_f);
               count_f=0;
               chnum=chnum+1;
           end
           if(i<=size(motor,1))
               fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,motor(i,1),motor(i,2),motor(i,3),b);   
               count_f=count_f+1;
               fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,motor(i,4),motor(i,5),motor(i,6),b);   
               count_f=count_f+2;
           else
               fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,0.0,0.0,0.0,b);   
               count_f=count_f+1;
               fprintf(f1,'ATOM  %5i  CA  ARG %s%4i    %8.3f%8.3f%8.3f  1.00%6.2f\n',count_f,chain(chnum),count_f,0.0,0.0,0.0,b);   
               count_f=count_f+2;
           end
       end
       if(count_f<9998)
       fprintf(f1,'TER   %5i      ARG %s%4i\n',count_f,chain(chnum),count_f);
       end
       fprintf(f1,'ENDMDL\n');
       count_model=count_model+1;
end
end
fclose(f1);
end