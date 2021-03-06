/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#include "output.h"
#include "particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void output_csv(system_t*system, int current_step){
    char filename[256];
    node_t* n;
    particle_t* p;
    sprintf(filename,"output_%u.csv",current_step);
    FILE*fp = fopen(filename,"w+");
    fprintf(fp, "id, x, y, z, vx, vy, vz, type, mass, rho, bvf_phi\n");
    for(n=system->particle_list->head; n!=NULL; n=n->next){
        p = n->data;
        fprintf(fp,"%u, %lf, %lf, %lf, %lf, %lf, %lf, %i, %lf, %lf, %lf\n",
            p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2],
            p->type, p->mass, p->rho, p->bvf_phi);
    }
    fclose(fp);
}


particle_t*output_buffer;
int output_buffer_size = 0;
int output_buffer_current_step;
int output_buffer_current_num_particles;
unsigned int* output_buffer_xx;
int output_buffer_xx_size = 0;
double* output_buffer_chem;
int output_buffer_chem_size = 0;

void output_vtk__sync_step(system_t*system, int current_step){
    output_buffer_current_step = current_step;
    if(output_buffer_size == 0){
        output_buffer_size = system->particle_list->count;
        output_buffer = (particle_t*) malloc(sizeof(particle_t)*output_buffer_size);
    }else if(output_buffer_size < system->particle_list->count){
        output_buffer = realloc(output_buffer, sizeof(particle_t)*output_buffer_size);
    }
    if(system->num_chem_species > 0){
        if(output_buffer_chem_size==0){
            output_buffer_chem_size = system->particle_list->count * system->num_chem_species;
            output_buffer_chem = (double*) malloc(sizeof(double)*output_buffer_chem_size);
        }else if(output_buffer_chem_size < system->particle_list->count * system->num_chem_species){
            output_buffer_chem_size = system->particle_list->count * system->num_chem_species;
            output_buffer_chem = realloc(output_buffer_chem, sizeof(double)*output_buffer_chem_size);
        }
    }
    int ncnt=0;
    node_t*n;
    for(n=system->particle_list->head; n!=NULL; n=n->next){
        memcpy( (void *) &output_buffer[ncnt++], (void *) n->data, sizeof(particle_t) );
        if(system->num_chem_species > 0){
            memcpy( (void *) &output_buffer_chem[n->data->id*system->num_chem_species], (void *) n->data->C, sizeof(double)*system->num_chem_species );
        }
    }
    output_buffer_current_num_particles = ncnt;
    if(system->rdme != NULL){
        // make a copy of the RDME state vector xx
        if(output_buffer_xx_size==0){
            output_buffer_xx_size = output_buffer_current_num_particles*system->num_stoch_species;
            output_buffer_xx = (unsigned int*) malloc(sizeof(unsigned int)*output_buffer_xx_size);
        }else if(output_buffer_xx_size < output_buffer_current_num_particles*system->num_stoch_species){
            output_buffer_xx_size = output_buffer_current_num_particles*system->num_stoch_species;
            output_buffer_xx = realloc(output_buffer_xx, sizeof(unsigned int)*output_buffer_xx_size);
        }
        /*
        printf("system->num_stoch_species*output_buffer_current_num_particles = %i\n",system->num_stoch_species*output_buffer_current_num_particles);
        printf("output_buffer_current_num_particles = %i\n",output_buffer_current_num_particles);
        printf("system->num_stoch_species = %i\n",system->num_stoch_species);
        int i;
        printf("xx = [ ");
        for(i=0;i<system->num_stoch_species*output_buffer_current_num_particles;i++){
            printf("%i ",system->rdme->xx[i]);
        }
        printf("]\n");
        printf("output_buffer_xx = [ ");
        for(i=0;i<system->num_stoch_species*output_buffer_current_num_particles;i++){
            printf("%i ",output_buffer_xx[i]);
        }
        printf("]\n");
        */
        //memcpy( (void*) &output_buffer_xx, (void*) system->rdme->xx, 
        //    sizeof(unsigned int)*system->num_stoch_species*output_buffer_current_num_particles);
        //for(int i=0;i<system->num_stoch_species*output_buffer_current_num_particles;i++){
        ncnt=0;
        for(n=system->particle_list->head; n!=NULL; n=n->next){
            memcpy( (void *) &output_buffer_xx[ncnt], (void *) n->data->xx, sizeof(unsigned int)*system->num_stoch_species );
            ncnt += system->num_stoch_species;
        }
    }
}
void output_vtk__async_step(system_t*system){
    FILE*fp;
    int i;
    char filename[256];
    int np = output_buffer_current_num_particles;
    if(output_buffer_current_step == 0){
        sprintf(filename, "output0_boundingBox.vtk");
        if(debug_flag){printf("Writing file '%s'\n", filename);}
        if((fp = fopen(filename,"w+"))==NULL){ 
            perror("Can't write 'output0_boundingBox.vtk'");exit(1);
        }
        fprintf(fp, "# vtk DataFile Version 4.1\n");
        fprintf(fp, "Generated by ssa_sdpd\n");
        fprintf(fp, "ASCII\n");
        fprintf(fp, "DATASET RECTILINEAR_GRID\n");
        fprintf(fp, "DIMENSIONS 2 2 2\n");
        fprintf(fp, "X_COORDINATES 2 double\n");
        fprintf(fp, "%lf %lf\n", system->xlo, system->xhi);
        fprintf(fp, "Y_COORDINATES 2 double\n");
        fprintf(fp, "%lf %lf\n", system->ylo, system->yhi);
        fprintf(fp, "Z_COORDINATES 2 double\n");
        fprintf(fp, "%lf %lf\n", system->zlo, system->zhi);
        fclose(fp);
    }
    sprintf(filename,"output%u.vtk",output_buffer_current_step);
    if(debug_flag){printf("Writing file '%s'\n", filename);}
    if((fp = fopen(filename,"w+"))==NULL){ 
        perror("Can't write output vtk file");exit(1);
    }
    fprintf(fp, "# vtk DataFile Version 4.1\n");
    fprintf(fp, "Generated by ssa_sdpd\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");
    fprintf(fp, "POINTS %i float\n",np);
    for(i=0;i<np;i++){
        //fprintf(fp, "%lf %lf %lf ",output_buffer[i].x[0],output_buffer[i].x[1],output_buffer[i].x[2]);
        fprintf(fp, "%.10e %.10e %.10e ",output_buffer[i].x[0],output_buffer[i].x[1],output_buffer[i].x[2]);
        if((i+1)%3==0){ fprintf(fp,"\n"); }
    }
    fprintf(fp,"\n");
    fprintf(fp, "VERTICES %i %i\n",np,2*np);
    for(i=0;i<np;i++){
        fprintf(fp,"1 %i\n",i);
    }
    fprintf(fp,"\n");
    fprintf(fp,"POINT_DATA %i\n", np);
    int num_fields = 7;
    if(system->rdme != NULL){
        num_fields += system->num_stoch_species;
    }
    if(system->num_chem_species > 0){
        num_fields += system->num_chem_species;
    }
    fprintf(fp,"FIELD FieldData %i\n",num_fields);//
    fprintf(fp,"id 1 %i int\n", np);
    for(i=0;i<np;i++){
        fprintf(fp, "%u ",output_buffer[i].id);
        if((i+1)%9==0){ fprintf(fp,"\n"); }
    }
    fprintf(fp,"\n");
    fprintf(fp,"type 1 %i int\n", np);
    for(i=0;i<np;i++){
        fprintf(fp, "%u ",output_buffer[i].type);
        if((i+1)%9==0){ fprintf(fp,"\n"); }
    }
    fprintf(fp,"\n");
    fprintf(fp,"v 3 %i double\n",np);
    for(i=0;i<np;i++){
        fprintf(fp, "%lf %lf %lf ",output_buffer[i].v[0],output_buffer[i].v[1],output_buffer[i].v[2]);
        if((i+1)%3==0){ fprintf(fp,"\n"); }
    }
    fprintf(fp,"\n");
    fprintf(fp,"rho 1 %i double\n", np);
    for(i=0;i<np;i++){
        fprintf(fp, "%lf ",output_buffer[i].rho);
        if((i+1)%9==0){ fprintf(fp,"\n"); }
    }
    fprintf(fp,"\n");
    fprintf(fp,"mass 1 %i double\n", np);
    for(i=0;i<np;i++){
        fprintf(fp, "%lf ",output_buffer[i].mass);
        if((i+1)%9==0){ fprintf(fp,"\n"); }
    }
    fprintf(fp,"\n");
    fprintf(fp,"bvf_phi 1 %i double\n", np);
    for(i=0;i<np;i++){
        fprintf(fp, "%lf ",output_buffer[i].bvf_phi);
        if((i+1)%9==0){ fprintf(fp,"\n"); }
    } 
    fprintf(fp,"\n");
    fprintf(fp,"nu 1 %i double\n", np);
    for(i=0;i<np;i++){
        fprintf(fp, "%lf ",output_buffer[i].nu);
        if((i+1)%9==0){ fprintf(fp,"\n"); }
    } 
    fprintf(fp,"\n");
    // loop here to check for continous species
    // c - concentration or continous? clarify at meeting
    if(system->num_chem_species > 0){
        int s;
        for(s=0;s<system->num_chem_species;s++){
            fprintf(fp,"C[%s] 1 %i double\n", system->species_names[s], np);
            for(i=0;i<np;i++){
                fprintf(fp, "%lf ",output_buffer_chem[i*system->num_chem_species+s] );
                if((i+1)%9==0){ fprintf(fp,"\n"); }
            }
            fprintf(fp,"\n");
        }
    }
    // d - discrete
    if(system->rdme != NULL){
        int s;
        for(s=0;s<system->num_stoch_species;s++){
            fprintf(fp,"D[%s] 1 %i int\n", system->species_names[s], np);
            for(i=0;i<np;i++){
                fprintf(fp, "%u ",output_buffer_xx[i*system->num_stoch_species+s] );
                if((i+1)%9==0){ fprintf(fp,"\n"); }
            }
            fprintf(fp,"\n");
        }
    }

    fclose(fp);

}

