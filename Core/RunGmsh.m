function RunGmsh(GeoFile, Parameters, NrOfdimensions, Output)

% Open the output file
outfile = fopen([Output '.geo'], 'w+');

% Write the definitions of the parameters
names = fieldnames(Parameters);
for n=1:length(names)
    if isfloat(Parameters.(names{n}))
        fprintf(outfile, [names{n} ' = %.16f;\n'], Parameters.(names{n}));
    end
end
fprintf(outfile, '\n');

% Open the original file
infile = fopen(GeoFile, 'r+');
% Copy its content into the output file
while ~feof(infile)
    x = fgets(infile);
    fprintf(outfile, '%s', x);
end

fclose(infile);
fclose(outfile);

% Run gmsh to produce the .msh file
system(['gmsh ' Output ['.geo -' num2str(NrOfdimensions) ' -v 0']]);

% Parse the .msh file to create the elements
Gmsh2CodeFem([Output '.msh'], Output);
