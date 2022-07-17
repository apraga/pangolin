%
% calcul des coordonnees des noeuds
% middle contient le milieu de 2 noeuds
% step contient la coordonnees d'un noeuds en fonction du precedent
% mesh_prev contient le numero de la maille superieure a laquelle le noeud appartient
% mesh_next contient le numero de la maille inferieure a laquelle le noeud appartient
function [middle,step,nb_nodes,mesh_prev,mesh_next] = sub_grille(nb_lat,...
    nb_lat2)
  nb_nodes = get_max_nodes(nb_lat,nb_lat2);

  mesh_prev = zeros(nb_lat,nb_nodes);
  mesh_next = zeros(nb_lat,nb_nodes);
  mesh_next = zeros(nb_lat,nb_nodes);
  step = zeros(nb_lat,nb_nodes);
  middle = zeros(nb_lat,nb_nodes);

  % Algorithm for computing the middle, previous et next cell at a given
  % latitude

  for i = 1:nb_lat2-1
    nb_cells = get_nb_cells(i,nb_lat,nb_lat2);
    nb_cells_next = get_nb_cells(i+1,nb_lat,nb_lat2);
    step_cur = 1./nb_cells;
    step_next = 1./nb_cells_next;
    i_cur = 0;
    i_next = 0;
    k = 1;
    cursor = 0; % Used for computing the steps and middle
    cursor_prev = 0;
    while (cursor < 1-eps)
      cursor_prev = cursor;
      tmp = (i_cur+1)*step_cur;
      tmp2 = (i_next+1)*step_next;
      mesh_prev(i,k) = i_cur+1;
      mesh_next(i,k) = i_next+1;

      if (tmp < tmp2)
        i_cur = i_cur + 1;
        cursor = tmp;
      else
        i_next = i_next + 1;
        cursor = tmp2;
      end
      step_tmp = cursor - cursor_prev;
      step(i,k) = step_tmp;
      middle(i,k) = cursor_prev + 0.5*step_tmp;
      k = k + 1;
    end
  end

  % Special case at the equator
  nb_cells = get_nb_cells(nb_lat2,nb_lat,nb_lat2);
  step_cur = 1./nb_cells;
  step(nb_lat2,1:nb_cells) = step_cur;
  middle(nb_lat2,1:nb_cells) = 0.5*step_cur:step_cur:1;
  mesh_next(nb_lat2,1:nb_cells) = 1:nb_cells;
  mesh_prev(nb_lat2,1:nb_cells) = 1:nb_cells;

  % South hemisphere
  i_prev = nb_lat2-1;
  for i = nb_lat2+1:nb_lat-1
    mesh_prev(i,:) = mesh_next(i_prev,:);
    mesh_next(i,:) = mesh_prev(i_prev,:);
    step(i,:) = step(i_prev,:);
    middle(i,:) = middle(i_prev,:);
    i_prev = i_prev - 1;
  end

  %% Checking
  %name = strcat('grille2.dat');
  %format = strcat('%f %f %f %f \n');
  %fid= fopen (name,'w');
  %if (fid == -1)
  %  disp('error opening file ');
  %end
  %count= fprintf (fid,'prev next step middle \n');


  %for i = 1:nb_lat-1
  %  nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
  %  for k = 1:nb_nodes
  %    count= fprintf (fid,format,mesh_prev(i,k),mesh_next(i,k),step(i,k),middle(i,k));
  %  end
  %end
  %fclose(fid);

end
