!*******************************************************************************
! Contains subroutines and function for output (writing files mostly : graph,
! geometry, postscript files
!*******************************************************************************
module output
  use grid
  implicit none

  contains
  !****************************************************************************
  ! Writing graph 
  ! Format : a graph at the Scotch format :
  !  nb_version
  !  n_vertices n_arcs
  !  basevalue flag
  !  [node i] n_arcs_i end_arc_i
  !****************************************************************************
  subroutine write_graph(id_file,filename,n_max,nb_lat,nb_lat2)
    integer :: n_vertices,n_arcs,n_max
    integer :: i,j,k,n_cells,n_cells_cur,n_zone
    integer :: size_max,sum_neighbours, id_file,nb_lat,nb_lat2
    character(len=*) :: filename
    integer,allocatable :: neighbours(:)

    ! At most 3 + 2 + 2 neighbours but we allocate 8 for symmetry
    size_max = 8
    allocate (neighbours(size_max))
    n_vertices = 3*nb_cells(n_max,nb_lat,nb_lat2)
    n_arcs = nb_arcs(n_max,nb_lat,nb_lat2)

    open(id_file,file=filename//".grf")
    n_cells = nb_cells(n_max,nb_lat,nb_lat2)
    sum_neighbours = 0
    ! Grite headers
    write (id_file,'(i4)') 0
    write (id_file,'(i7,i7)') n_vertices,n_arcs
    write (id_file,'(i4,a4)') 1,'000'
    do i = 1,n_max
      n_cells_cur = nb_cells_lat(i,nb_lat,nb_lat2)
      do n_zone = 0,2
        do j = 1,n_cells_cur
          call find_neighbours(j,i,n_zone,n_cells,n_max,nb_lat,nb_lat2,size_max,neighbours)
          call get_neighbours(id_file,size_max,neighbours,n_cells,sum_neighbours)
        end do
      end do
    end do
    if (sum_neighbours /= n_arcs) then
      write (*,*) "Wrong number of arcs : ",sum_neighbours,n_arcs
    else
      write (*,*) "number of arcs : ",sum_neighbours,n_arcs
    end if
    close(id_file)
    deallocate (neighbours)
  end subroutine

  !****************************************************************************
  ! Writing geometry file (scotch format)
  !****************************************************************************
  subroutine write_geom(id_file,filename,n_max,nb_lat,nb_lat2)

    integer :: n_vertices,n_max, id_file,nb_lat,nb_lat2
    integer :: i,j,k,n_cells,size_max
    real :: lat, lon, dlat, dlon
    character(len=*),parameter  :: format_i='(i5,f10.5,f10.5)'
    character(*) :: filename

    n_vertices = 3*nb_cells(n_max,nb_lat,nb_lat2)
    dlat = 90./nb_lat2

    open(id_file,file=filename//".xyz")
    ! Write headers
    write (id_file,'(i7)') 2
    write (id_file,'(i7)') n_vertices
    k = 0
    do i = 1,n_max
      n_cells = 3*nb_cells_lat(i,nb_lat,nb_lat2)
      dlon = 360./n_cells
      lat = 90.5-i*dlat
      do j = 1,n_cells
        lon = (j-0.5)*dlon
        ! k is the label
        write (id_file,format_i) k,lon,lat
        k = k + 1
      end do
    end do
    close(id_file)
  end subroutine

  !****************************************************************************
  ! Writing macros for drawing
  !****************************************************************************
  subroutine write_macros(id_file,r)
    integer,intent(in) :: id_file,r
    integer :: i,nb_colors
    ! Color list for partitioning
    ! Beware, we need to transpose as the array is initialized lines by lines
    real,dimension(16,3) :: colors_map
    colors_map= transpose(reshape((/&
    &1.00, 0.00, 0.00,& ! Red
    &0.00, 1.00, 0.00,& ! Green        
    &1.00, 1.00, 0.00,& ! Yellow       
    &0.00, 0.00, 1.00,& ! Blue         
    &1.00, 0.00, 1.00,& ! Magenta      
    &0.00, 1.00, 1.00,& ! Cyan         
    &1.00, 0.50, 0.20,& ! Orange       
    &0.30, 0.55, 0.00,& ! Olive        
    &0.72, 0.47, 0.47,& ! Dark pink    
    &0.33, 0.33, 0.81,& ! Sea blue     
    &1.00, 0.63, 0.63,& ! Pink         
    &0.62, 0.44, 0.65,& ! Violet       
    &0.60, 0.80, 0.70,& ! Pale green   
    &0.47, 0.20, 0.00,& ! Brown        
    &0.00, 0.68, 0.68,& !Turquoise    
    &0.81, 0.00, 0.40&  ! Purple     
    &/), (/size(colors_map,2),size(colors_map,1)/))) 

    ! Function for drawing circles
    write (id_file,'(a)') "% Draw circle of radius 10"
    write (id_file,'(a)') "/circle {"
    write (id_file,'(i3,a)') r,"  0 360 arc closepath"
    write (id_file,'(a)') "  fill"
    write (id_file,'(a)') "  stroke"
    write (id_file,'(a)') "} def"
    write (id_file,'(a)') 

    ! Write color functions
    nb_colors = size(colors_map,1)
    do i = 1,nb_colors
      write (id_file,'(a,a,a)',advance='no') "/",color_name(i-1),"{ "
      write (id_file,'(f5.2,f5.2,f5.2)',advance='no') colors_map(i,1),colors_map(i,2),colors_map(i,3)
      write (id_file,'(a)') " setrgbcolor } def"
    end do
    write (id_file,*)
  end subroutine

  !****************************************************************************
  ! Compute the name of the color macro number i (with i >= 0)
  !****************************************************************************
  function color_name(i) 
    integer, intent(in) :: i
    integer i2
    character(7) color_name
    i2 = modulo(i,16)
    color_name = char(iachar('a')+i2)
    color_name = "color_" // color_name
  end function
  
  
  !****************************************************************************
  ! Get dimensions of bounding box, offset and scaling
  !****************************************************************************
  subroutine get_dimensions(width,height,offset_y,scaling_y,offset_x,scaling_x, &
    r,projection,nb_lat,nb_lat2)
    integer, intent(out) :: offset_y, offset_x, width,height
    real, intent(out) :: scaling_y, scaling_x
    integer :: nb_lat,nb_lat2,r,projection
    real :: min_dlon, min_dlat

    if (projection == 0) then 
      offset_y = 91
      offset_x = 1
      min_dlon = 120./nb_cells_lat(nb_lat2,nb_lat,nb_lat2)
      min_dlat = 1.
      ! We want dlon*scaling_x >= 2*r
      scaling_x = 2*r/min_dlon
      width = 360*scaling_x
      ! Idem for y
      scaling_y = 2*r/min_dlat
      height = 180*scaling_y
    else
      offset_y = 180
      offset_x = 180
      scaling_x = 90.
      scaling_y = 90.
      width = 400
      height = 400
    end if

  end subroutine

  !****************************************************************************
  ! Writing mesh on postscript format
  !****************************************************************************
  subroutine write_ps(filename,parttab,n_max,nb_lat,nb_lat2)
    ! 0 for lat-lon, 1 for orthographic
    integer :: projection = 0
    integer :: n_max, nb_lat,nb_lat2,n_cells
    integer :: id_file,i,j,i_color,k
    character(len=*),parameter  :: format_i='(f10.2,f10.2,a)'
    character(*) :: filename
    integer :: width,height,offset_x,offset_y
    real :: scaling_x,scaling_y
    integer,allocatable :: parttab(:)
    integer :: r = 1

    ! Set variables for display
    call get_dimensions(width,height,offset_y,scaling_y,offset_x,scaling_x,r,&
    projection,nb_lat,nb_lat2)

    id_file = 3
    open(id_file,file=filename//".ps")
    
    ! Bounding box
    write (id_file,'(a)') "%!PS-Adobe-3.0 EPSF-3.0"
    write (id_file,'(a)') "%%Title: Partition of the reduced grid"
    write (id_file,'(a,i5,i5,i5,i5)') "%%BoundingBox: ",0,0,width,height
    write (id_file,'(a)') "%%Pages: 1"
    write (id_file,'(a)') "%%Page: 1 1"
    write (id_file,'(a)') "%%EOF"
    write (id_file,*) 

    ! Write macros
    call write_macros(id_file,r)

    ! Write mesh 
    if (projection == 0) then
      call write_mesh_latlon(parttab,id_file,offset_y,scaling_y,offset_x,&
      scaling_x,n_max,nb_lat,nb_lat2)
    else
      call write_mesh_ortho(parttab,id_file,offset_y,scaling_y,offset_x,&
      scaling_x,n_max,nb_lat,nb_lat2)
    end if
    !call write_mesh2(parttab,id_file,offset_y,scaling_y,offset_x,scaling_x,n_max,nb_lat,nb_lat2)
    close(id_file)
  end subroutine

  
  !****************************************************************************
  ! Compute orthographic tttHection from point in spherical coordinates
  ! (lat,lon)
  ! With clipping if the point should not be displayed
  !****************************************************************************
  subroutine orthographic(x,y,clipping,lat,lon)
   real :: x,y, coef,r,lon0,lat0
   real :: lat2,lon2, cos_c
   real, intent(in) :: lat,lon 
   integer :: clipping 
   real :: pi = 3.14159265
   r = 2.
   lat0 = 0.
   lon0 = 0.
   coef = pi/180
   lat2 = lat*coef
   lon2 = lon*coef
   cos_c = sin(lat0)*sin(lat2) + cos(lat0)*cos(lat2)*cos(lon2-lon0)
   clipping = 0
   if (cos_c < 0.) then 
     clipping = -1
     return
   end if
   x = r*cos(lat2)*sin(lon2 - lon0)
   y = r*(cos(lat0)*sin(lat2)-sin(lat0)*cos(lat2)*cos(lon2 - lon0))
  end subroutine

  !****************************************************************************
  ! Write in orthographic coordinates
  !****************************************************************************
  subroutine write_mesh_ortho(parttab,id_file,offset_y,scaling_y,offset_x,&
    scaling_x,n_max,nb_lat,nb_lat2)
    integer :: n_max, nb_lat,nb_lat2,n_cells
    integer,intent(in) :: id_file
    integer :: i,j,k,i_color
    character(len=*),parameter :: format_i='(f10.2,f10.2,a)'
    real :: dlat,dlon,lat,lon
    integer :: offset_x,offset_y,clipping
    real :: scaling_x,scaling_y,pos_x,pos_y
    integer,allocatable :: parttab(:)
 
    dlat = 90./nb_lat2
    k = 0
    do i = 1,n_max
      n_cells = 3*nb_cells_lat(i,nb_lat,nb_lat2)
      dlon = 360./n_cells
      lat = 90.5-i*dlat
      do j = 1,n_cells
        lon = (j-0.5)*dlon
        call orthographic(pos_x,pos_y,clipping,lat,lon)
        if (clipping > -1) then
          pos_x = offset_x + pos_x*scaling_x
          pos_y = offset_y + pos_y*scaling_y
          i_color = parttab(k)
          write (id_file,'(a)') color_name(i_color)
          write (id_file,format_i) pos_x,pos_y," circle"
        end if
        k = k + 1
      end do
    end do
  end subroutine

  !****************************************************************************
  ! Write in lat-lon coordinates
  !****************************************************************************
  subroutine write_mesh_latlon(parttab,id_file,offset_y,scaling_y,offset_x,&
    scaling_x,n_max,nb_lat,nb_lat2)
    integer :: n_max, nb_lat,nb_lat2,n_cells
    integer,intent(in) :: id_file
    integer :: i,j,k,i_color
    character(len=*),parameter :: format_i='(f10.2,f10.2,a)'
    real :: dlat,dlon,lat,lon
    integer :: offset_x,offset_y,clipping
    real :: scaling_x,scaling_y,pos_x,pos_y
    integer,allocatable :: parttab(:)
 
    dlat = 90./nb_lat2
    k = 0
    do i = 1,n_max
      n_cells = 3*nb_cells_lat(i,nb_lat,nb_lat2)
      dlon = 360./n_cells
      lat = (offset_y + 90.5-i*dlat)*scaling_y
      do j = 1,n_cells
        lon = offset_x + (j-0.5)*dlon*scaling_x
        i_color = parttab(k)
        write (id_file,'(a)') color_name(i_color)
        write (id_file,format_i) lon,lat," circle"
        k = k + 1
      end do
    end do
  end subroutine
end module
