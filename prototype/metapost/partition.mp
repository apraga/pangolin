prologues := 2; % Necessary for good BB !
filenametemplate "%j-%c.mps"; % for v0.9 or inferior
%outputtemplate := "%j-%c.mps"; % for v1.0 or superior
warningcheck := 0; % No warning about numbers > 4096
u:=1.5cm;

% Partitions
beginfig(1);
  input lib_plot
  nb_lat2 = 90;
  nb_lat = 2*nb_lat;
  draw_partition_fill(nb_lat2);
endfig;

end
