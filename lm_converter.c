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
    int zoneid, baseid, index_sect, nbndry, iparent_flag, nbocos, normalindex[3];

    char zone_name[25];
    char sectionname[30];
    char boconame[33];


    /* CGNS datatypes */

    ElementType_t itype, iend;
    GridLocation_t igr;
    cgsize_t size[3],istart,iparentdata;
    BCType_t ibocotype;
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


    /* Now we will get the node coordnates of the mesh so lets allocate some arrays
     * to hold useful information such as connectivity and node coordnates. Lets st
     * art with the node coordnates. */ 

    double * x_coord = (double*)malloc(size[0]*sizeof(double));
    double * y_coord = (double*)malloc(size[0]*sizeof(double));

    cgsize_t start = 1;

    printf("----------------------------------------------------------\n");
    printf(" ++ Reading node coordnates...\n ");
    printf("----------------------------------------------------------\n\n");


    /* Lets get the node coordnates */

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

    int kindof;

    int ele;
    fprintf(f_connec,"%d\n",size[1]);

    for (index_sect=1; index_sect <= nsections; index_sect++){

        cg_section_read( cg_file, baseid, zoneid, index_sect, sectionname, 
                         &itype, &istart, &end, &nbndry, &iparent_flag);

        printf("  + Reading section data...\n");
        printf("  + section name=%s\n",sectionname);
        printf("  + section type=%s\n",ElementTypeName[itype]);
        printf("  + istart,end=%i, %i\n\n",(int)istart,(int)end);


        /* Check if the elemnt type is two dimensional, here I am considering 
         * the mesh will have triangles and/or quad elements. Cheer up it 
         * works fine. */

        if (itype == QUAD_4){


            /* Declare our elemn connectivity vector */

            cgsize_t ielem[size[1]][4];


            /* Put connectivity information in ielem */

            cg_elements_read( cg_file, baseid, zoneid, index_sect, ielem[0], &iparentdata);


            /* Output the connectivity to the file */

            for (ele = 0; ele < end; ele++){

                kindof = 4;

                fprintf(f_connec,"%d %d %d %d %d %d\n",ele+1,kindof,(int)ielem[ele][0],(int)ielem[ele][1],(int)ielem[ele][2],(int)ielem[ele][3]);
            }

        }else if (itype == TRI_3){


            /* Declare our elemn connectivity vector */

            cgsize_t ielem[size[1]][3];


            /* Put connectivity information in ielem */

            cg_elements_read( cg_file, baseid, zoneid, index_sect, ielem[0], &iparentdata);


            /* Output the connectivity to the file */

            for (ele = 0; ele < end; ele++)

                kindof = 3;

                fprintf(f_connec,"%d %d %d %d %d\n",ele+1,kindof,(int)ielem[ele][0],(int)ielem[ele][1],(int)ielem[ele][2]);

        }else if (itype == MIXED){

            cgsize_t ielem[size[1]][4];
        }
    } 

    
    /* Now that we have (almost) our connectivity, lets get the last information
     * that is useful for us, the B.Cs. This procedure in Carlos code is very el
     * aborated and well done so lets try to do our own version here. */

    printf("----------------------------------------------------------\n");
    printf(" ++ Reading Boundary conditions...\n ");
    printf("----------------------------------------------------------\n\n");

    cg_nbocos( cg_file, baseid, zoneid, &nbocos);

    printf(" + Total number of B.Cs : %d\n", nbocos);


    /* Now lets loop through all boundary conditions man ! */

    int ib;

    for (ib=1; ib <= nbocos; ib++){

        cg_goto(cg_file, baseid,"Zone_t",1,"ZoneBC_t",1,"BC_t",ib,"end");
        cg_gridlocation_read(&igr);

        cg_boco_info( cg_file, baseid, zoneid, ib, boconame, &ibocotype,
                &iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);

        if (iptset != PointList){

            printf("\nError.  For this program, BCs must be set up as PointList type %s\n",PointSetTypeName[iptset]);
            exit(1);
        }


        printf("  + BC number: %i\n",ib);
        printf("  + name= %s\n",boconame);
        printf("  + type= %s\n",BCTypeName[ibocotype]);
        printf("  + no of elements= %i\n\n",(int)npts);

        cgsize_t ipnts[npts];

        cg_boco_read(cg_file, baseid, zoneid, ib, ipnts, &normallist);


        /* Get boundary condition type from the user */

        printf("\n===============================================\n");
        printf(" Please, select the type of B.C for: %s \n",boconame);
        printf("   - Supersonic Inlet  : 1\n");
        printf("   - Supersonic Outlet : 2\n");
        printf("   - Euler Wall        : 3\n");
        printf("   - Internal FLUID    : 4\n");
        printf("   - Symmetry          : 5\n");
        printf("===============================================\n");
        printf("Please, input boundary condition type: ");

        scanf("%d",&spaTyp);


        /* I dont want the fluid elemnets in my file ! */

        if (spaTyp != 4){

            int i;

            for ( i = 0; i < (int)npts-1; i++) {
                fprintf(f_boundc,"%d %d %d\n", spaTyp,ipnts[i],ipnts[i+1]);
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
