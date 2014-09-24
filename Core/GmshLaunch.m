function GmshLaunch(gmeshFile,ParameterName,ParameterValue,NrOfParameters,NrOfdimensions,outfile) 
% launches gmesh meshing
% Gmesh should be recognized as an internal command -> put the executable in the System32 folder in Windows 
% 1. Modifies the parameters in the .geo file
% 2. Launches gmesh mesh with NrOfdimensions
% 3. Convert the .msh to the CodeFem format 

% 1. Modify the entries of the .geo file
fid = fopen([gmeshFile '.geo'],'r+');
frewind(fid);
for iParameter=1:NrOfParameters
    fprintf(fid,[ParameterName{iParameter} ' = %10.5f;\n'],ParameterValue(iParameter));
end
fprintf(fid,'\n\n\n\n ///////////////////');
%
fclose(fid);
% 2. Launch gmsh
dos(['gmsh ' gmeshFile ['.geo -' num2str(NrOfdimensions) ' -v 0']]);
% 3. Convert the .msh to the CodeFem format 
Gmsh2CodeFem([gmeshFile '.msh'],outfile); 