classdef filaments
    properties
        coord_cell1
        coord_cell2
        coord_cellinterp
        forces
    end
    methods
        function coord_matrix=get_matrix1(obj)
            coord_matrix=[];
                   for fils=1:length(obj.coord_cell1)
                       coord_matrix=[coord_matrix;reshape(cell2mat(obj.coord_cell1(fils)),3,[])'];
                   end
        end
         function coord_matrix=get_matrix2(obj)
            coord_matrix=[];
                   for fils=1:length(obj.coord_cell2)
                       coord_matrix=[coord_matrix;reshape(cell2mat(obj.coord_cell2(fils)),3,[])'];
                   end
        end
    end
end
