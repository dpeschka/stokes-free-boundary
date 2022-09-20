% output matlab tetraeder mesh with scalar values assigned to vertices
% to a VTK mesh for visualization, e.g. using paraview

function [ id ] = outvtk_tri_scalar(fname,x,y,e2p,u,sname)

npoint = max(max(e2p(:,1:3)));
ntet   = size(e2p,1);

fprintf('npoint = %i   ntet = %i\n',npoint,ntet);

[fid,id] = fopen(fname,'w+');

% output header
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'tet mesh created from MATLAB\n');
fprintf(fid,'ASCII\n');

% create unstructured grid with points
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i float\n',npoint);
for i=1:npoint
    fprintf(fid,'%f %f 0\n',x(i),y(i));
end

% output cell connectivity
fprintf(fid,'CELLS %i %i\n',ntet,4*ntet);
for i=1:ntet
    fprintf(fid,'3 %i %i %i\n',e2p(i,1:3)-1);
end

% output cell type 10 = tetraeder
fprintf(fid,'CELL_TYPES %i\n',ntet);
for i=1:ntet
    fprintf(fid,'5\n');
end

fprintf(fid,'POINT_DATA %i\n',npoint);

% output scalar data
if ~iscellstr(sname)
    fprintf(fid,'SCALARS %s float 1\n',sname);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:npoint
        fprintf(fid,'%f\n',u(i));
    end
else
    for k=1:length(sname)
        fprintf(fid,'SCALARS %s float 1\n',sname{k});
        fprintf(fid,'LOOKUP_TABLE default\n');
        for i=1:npoint
            fprintf(fid,'%f\n',u{k}(i));
        end
    end
end
%end

%if ~isempty(vx)
%    fprintf(fid,'VECTORS vectors float\n');
%    for i=1:npoint
%        fprintf(fid,'%f %f %f\n',vx(i),vy(i),vz(i));
%    end
%end

fclose(fid);