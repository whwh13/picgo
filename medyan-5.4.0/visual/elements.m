classdef elements
    properties
        coord_cell
        id
    end
    methods
        function coord_matrix=get_matrix(obj)
            coord_matrix=[];
        for fils=obj.coord_cell
            coord_matrix=[coord_matrix;reshape(cell2mat(fils),3,[])'];
        end
        end
    end
end
