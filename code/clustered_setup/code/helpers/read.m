function output = read(filename)

input=fopen(filename);
output=fscanf(input,'%g');
fclose(input);

end