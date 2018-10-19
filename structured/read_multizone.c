#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cgnslib.h"


int main(void){

    printf("Read and feed the data structures of a multizone grid.\n");

    /* Open cgns file and return the id */

    int cg_file;

    /* Open cgns file and return the index */

    if (cg_open("m_zone_square.cgns",CG_MODE_READ,&cg_file)) cg_error_exit();

    /* Set the indexes as our mesh is very simple. */

    int index_base = 1; int nzone = 0;

    /* Read the number of zones. */

    cg_nzones(cg_file,index_base,&nzone);

    /* Set the sizes. */

    int size_imax_z[nzone];
    int size_jmax_z[nzone];

    /* Dump the number of zones in the mesh. */

    printf("The number of zones in the mesh is: %d\n",nzone);

    /* Loop through zones and get the size data we need. */

    char zonename[33]; cgsize_t isize[3][3],irmin[3],irmax[3];

    for (int nz = 1; nz <= nzone; nz++){

        /* Read the zone data. */

        cg_zone_read(cg_file,index_base,nz,zonename,isize[0]);

        /* Get the lower index. */

        irmin[0]=1;
        irmin[1]=1;
        irmin[2]=1;

        /* Get the upper index. */

        irmax[0]=isize[0][0];
        irmax[1]=isize[0][1];
        irmax[2]=isize[0][2];

        /* Store the sizes. */

        size_imax_z[nz-1] = irmax[0];
        size_jmax_z[nz-1] = irmax[1];

    }

    /* Get the maximun number of points. */

    int max_imax = 0; int max_jmax = 0;

    for (int nz = 0; nz < nzone; nz++) max_imax += size_imax_z[nz];
    for (int nz = 0; nz < nzone; nz++) max_jmax += size_jmax_z[nz];

    /* Dump for debug. */

    for (int nz = 0; nz < nzone; nz++)
        printf("Block # %d has %d nodes in i and %d nodes in j.\n",nz,size_imax_z[nz],size_jmax_z[nz]);

    /* Export to see what happens. */

    FILE * f_out = fopen("solution.dat", "w");

    /* Check the condition of the pointer. */

    if (f_out == NULL){
        printf("\nERROR: Output file cannot be opened !\n"); exit(1);
    }

    /* Print the header of the tecplot file. */

    fprintf(f_out,"TITLE = \" multi-zone 01 \"\n");
    fprintf(f_out,"VARIABLES = \"X\" , \"Y\" , \"Z\"\n");
    fprintf(f_out,"ZONE T =\"Zone-one\", I= %d ,J= %d , K=1, DATAPACKING=POINT\n",max_imax, max_jmax);

    /* Ree-do the loop and get the coordinates. */

    for (int nz = 1; nz <= nzone; nz++){

        /* Read the zone data. */

        cg_zone_read(cg_file,index_base,nz,zonename,isize[0]);

        /* Get the lower index. */

        irmin[0]=1;
        irmin[1]=1;
        irmin[2]=1;

        /* Get the upper index. */

        irmax[0]=isize[0][0];
        irmax[1]=isize[0][1];
        irmax[2]=isize[0][2];

        /* Set the aux vectors. */

        double x_aux[(int)irmax[1]][(int)irmax[0]];
        double y_aux[(int)irmax[1]][(int)irmax[0]];

        /* Get the coordinates. */

        cg_coord_read(cg_file,index_base,nz,"CoordinateX",CGNS_ENUMV(RealDouble),irmin,irmax,x_aux[0]);
        cg_coord_read(cg_file,index_base,nz,"CoordinateY",CGNS_ENUMV(RealDouble),irmin,irmax,y_aux[0]);

        printf("Reading coordinates at zone # %d\n",nz-1);

        /* Export the file. */

        for (int i = 0; i<(int)isize[0][0]; i++){
            for (int j = 0; j<(int)isize[0][1]; j++){

                printf(" Zone: %d | x coordinate: %lf | y coordinate %lf\n",nz,x_aux[i][j],y_aux[i][j]);

                fprintf(f_out,"%lf %lf %lf\n",x_aux[j][i],y_aux[j][i],0.0);
            }
        }
    }

    /* Close out file. */
    
    fclose(f_out);

    /* Close the cgns file. */

    cg_close(cg_file);

    return 0;
} 
