/* Converts a 2D cgns point based 2.4 to ascii format.
 * Author: Leonardo Motta */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"


/* Prototypes */
void  cgns2ascii(const char * filename);


/* Main program */
int main(int argc, char * argv[]){

    if (argc < 2) {    
        printf("ERROR: One argument expected.\n");
        exit(1);
    }

    char * filename = argv[1];

    printf("----------------------------------------------------------\n");
    printf(" ++ Reading mesh file: %s\n", filename);
    printf("----------------------------------------------------------\n\n");

    cgns2ascii(filename);

}


/* Function that converts cgns to sparta code */
void cgns2ascii(const char * filename){

    int cg_file, nbases, nzones, nsections, ndataset, normallist, spaTyp, npe;
    int zoneid, baseid, index_sect, nbndry, iparent_flag, normalindex[3];

    char zone_name[25];
    char sectionname[30];

    /* CGNS datatypes */
    ElementType_t itype, iend;
    GridLocation_t igr;
    cgsize_t size[3],istart,iparentdata;
    PointSetType_t iptset;
    cgsize_t npts,normallistflag,end,ElementDataSize;
    DataType_t normaldatatype;
    
    FILE * f_points;
    FILE * f_connec;
    FILE * f_boundc;


    /* Openning our mesh files ...*/
    f_points = fopen("pointCord.dat","w+");
    f_connec = fopen("elemnConn.dat","w+");
    f_boundc = fopen("elemnBonc.dat","w+");


    printf("----------------------------------------------------------\n");
    printf(" ++ Reading basic mesh data...\n ");
    printf("----------------------------------------------------------\n\n");


    /* Openning the cgnsfile. In the mid-level library, all native cgns functio
     * ns are responsable for navigate through the data tree and get everything
     * I need. */

    if ( cg_open( filename, CG_MODE_READ, &cg_file ) ) cg_error_exit();


    /* Now lets get the number of bases in the mesh. Until now I dont know exac
     * ly what these bases actually are but I know that for my simple code
     * I dont need or want more than one of these, so lets check if these are 
     * one. */

    cg_nbases( cg_file, &nbases );

    if ( nbases > 1 ) {
        printf("CGNS file nbases greater than one is not supported.\n");
        exit(1);
    } else {
        printf("  + The number of bases in the file is: %d\n", nbases);
    }

    /* Now lets get the number of zones in the data-tree. */
    for ( baseid = 1; baseid <= nbases; baseid++ ){
        cg_nzones( cg_file, baseid, &nzones );
        if ( nzones > 1 ) {
          printf("CGNS file nzones greater than one is not supported.\n");
          exit(1);
        }
    }


    /* Baseid is being incremented at the end of the loop since it is a very 
     * important variable, we need it to stay as the value given by cg_zones 
     * function. */

    baseid = baseid - 1;


    /* Now its big time for us ! We will get some very nice grid information !
     * This basic information will be usefull for us later when building the 
     * arrays that will support the mesh */

    for ( zoneid = 1; zoneid <= nzones; zoneid++ ){

        cg_zone_read( cg_file, baseid, zoneid, zone_name, size );
        cg_nsections( cg_file, baseid, zoneid, &nsections );

        printf("   + CGNS Zone name  : %s\n", zone_name);
        printf("   + Number of nodes : %d\n", size[0]);
        printf("   + Number of elemns: %d\n", size[1]);
        printf("   + Number of parts : %d\n", nsections);
    }

    zoneid = zoneid - 1;


    /* Now we will get the node coordinates of the mesh so lets allocate some arrays
     * to hold useful information such as connectivity and node coordinates. Lets st
     * art with the node coordinates. */ 

    double * x_coord = (double*)malloc(size[0]*sizeof(double));
    double * y_coord = (double*)malloc(size[0]*sizeof(double));

    cgsize_t start = 1;

    printf("----------------------------------------------------------\n");
    printf(" ++ Reading node coordnates...\n ");
    printf("----------------------------------------------------------\n\n");


    /* Lets get the node coordinates */
    cg_coord_read( cg_file, baseid, zoneid, "CoordinateX", RealDouble, &start, &size[0], x_coord);
    cg_coord_read( cg_file, baseid, zoneid, "CoordinateY", RealDouble, &start, &size[0], y_coord);


    /* Writing point data */
    int coord;
    fprintf(f_points,"%d\n",size[0]);
    for ( coord = 0; coord<size[0]; coord++)
        fprintf(f_points,"%d %.17e %.17e\n",coord+1,x_coord[coord],y_coord[coord]);

    printf("  + Nodes were red succesfully ! \n\n ");

    /* If structured grids were disered the problem was solved but there is no
     * free lunch for us so lets get the boundary conditions */ 
    printf("----------------------------------------------------------\n");
    printf(" ++ Reading element connectivity...\n ");
    printf("----------------------------------------------------------\n\n");

    cg_nsections( cg_file, baseid, zoneid, &nsections);

    printf(" + Number of sections found in the mesh: %i\n",nsections);

    /* Lets loop through sections and find out our element connectivity */

    int ele;
    fprintf(f_connec,"%d\n",size[1]);

    for (index_sect=1; index_sect <= nsections; index_sect++){

        cg_section_read( cg_file, baseid, zoneid, index_sect, sectionname, 
                         &itype, &istart, &end, &nbndry, &iparent_flag);

        printf("  + Reading section data...\n");
        printf("  + section name=%s\n",sectionname);
        printf("  + section type=%s\n",ElementTypeName[itype]);
        printf("  + istart,end=%i, %i\n\n",(int)istart,(int)end);

        /* Here we will write the connectivity information. It's important to point out
         * that this convertion code works for 2D grids only. This means that, any BAR_2
         * element is boundary and any MIXED, TRI_3 or QUAD_4 is interior element */

        if (itype == QUAD_4){

            /* Declare our elemn connectivity vector */
            cgsize_t ielem[size[1]][4];
            int n_vol = (int)size[1];

            /* Put connectivity information in ielem */
            cg_elements_read( cg_file, baseid, zoneid, index_sect, ielem[0], &iparentdata);

            /* Output the connectivity to the file */
            for (ele = 0; ele < n_vol; ele++){
                fprintf(f_connec,"%d %d %d %d %d %d\n",ele+1,4,
                        (int)ielem[ele][0],(int)ielem[ele][1],(int)ielem[ele][2],(int)ielem[ele][3]);
            }

        }else if (itype == TRI_3){

            /* Declare our elemn connectivity vector */
            cgsize_t ielem[size[1]][3];
            int n_vol = (int)size[1];

            /* Put connectivity information in ielem */
            cg_elements_read( cg_file, baseid, zoneid, index_sect, ielem[0], &iparentdata);

            /* Output the connectivity to the file */
            for (ele = 0; ele < n_vol; ele++){
                fprintf(f_connec,"%d %d %d %d %d\n",ele+1,3,
                        (int)ielem[ele][0],(int)ielem[ele][1],(int)ielem[ele][2]);
            }

        }else if (itype == MIXED){

            /* Note: The way cgns deals with MIXED elements is, until now, very
             * obscure for me. In the manual they tell us that one extra integer
             * is used up front the ielem array in order to tell me the element 
             * type. It worked perfectly fine for TRI_3 elements but, for QUAD_4
             * one strange 7 appears in the diagonal. I fiexd it in a very un-el
             * elegant way so, any suggestions are welcome ! */

            /* Declare our elemn connectivity vector */
            cgsize_t ielem[size[1]][4];
            int n_vol = (int)size[1];

            /* Put connectivity information in ielem */
            cg_elements_read( cg_file, baseid, zoneid, index_sect, ielem[0], &iparentdata);

            int j = 0; 

            for (ele = 0; ele < n_vol; ele++){

                /* If the element is mixed, cgns will give 5 as value for this type of element */
                if ((int)ielem[ele][0] == 5){
                    fprintf(f_connec,"%d %d %d %d %d\n",ele+1,3,
                            (int)ielem[ele][1],(int)ielem[ele][2],(int)ielem[ele][3]);
                } else {
                    fprintf(f_connec,"%d %d %d %d %d %d\n",ele+1,4,
                            (int)ielem[ele][j+1],(int)ielem[ele][j+2],(int)ielem[ele][j+3],(int)ielem[ele][j+4]);
                    j = j + 1;
                }
            }
        }else if (itype == BAR_2){

            int spaTyp;

            /* Get boundary condition type from the user */
            printf("\n===============================================\n");
            printf(" Please, select the type of B.C for: %s \n",sectionname);
            printf("   - Farfield          : 1\n");
            printf("   - Outlet            : 2\n");
            printf("   - Inviscid Wall     : 3\n");
            printf("   - Symmetry          : 4\n");
            printf("===============================================\n");
            printf("Please, input boundary condition type: ");
            scanf("%d",&spaTyp);

            /* Discover the number of B.Cs */
            int bc_size = (int)end - (int)istart;

            /* Declare our elemn connectivity vector */
            cgsize_t ielem[bc_size][2];

            /* Put connectivity information in ielem */
            cg_elements_read( cg_file, baseid, zoneid, index_sect, ielem[0], &iparentdata);

            /* Output the connectivity to the file */
            for (ele = 1; ele <=  bc_size; ele++){
                fprintf(f_boundc,"%d %d %d\n",spaTyp,
                        (int)ielem[ele][0],(int)ielem[ele][1]);
            }
        }
    }
    
    cg_close (cg_file);

    free(x_coord);
    free(y_coord);

    fclose(f_points);
    fclose(f_connec);
    fclose(f_boundc);

}
