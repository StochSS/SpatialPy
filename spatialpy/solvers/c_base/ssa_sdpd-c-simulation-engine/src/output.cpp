/**
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2019 - 2022 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "output.h"
#include "particle_system.hpp"

namespace Spatialpy{
    void output_csv(ParticleSystem*system, int current_step){
        char filename[256];
        Particle* p;
        sprintf(filename,"output_%u.csv",current_step);
        FILE*fp = fopen(filename,"w+");
        fprintf(fp, "id, x, y, z, vx, vy, vz, type, mass, rho, bvf_phi, C\n");
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p = &system->particles[i];
            fprintf(fp,"%u, %lf, %lf, %lf, %lf, %lf, %lf, %i, %lf, %lf, %lf, %f\n",
                p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2],
                p->type, p->mass, p->rho, p->bvf_phi, p->C[current_step]);
        }
        fclose(fp);
    }


    Particle *output_buffer;
    long unsigned int output_buffer_size = 0;
    int output_buffer_current_step;
    int output_buffer_current_num_particles;
    unsigned int* output_buffer_xx;
    long unsigned int output_buffer_xx_size = 0;
    double* output_buffer_chem;
    long unsigned int output_buffer_chem_size = 0;

    void output_vtk__sync_step(ParticleSystem*system, int current_step){
        output_buffer_current_step = current_step;
        if(output_buffer_size == 0){
            output_buffer_size = system->particles.size();
            output_buffer = (Particle*) malloc(sizeof(Particle)*output_buffer_size);
        }else if(output_buffer_size < system->particles.size()){
            output_buffer = (Particle*) realloc(output_buffer, sizeof(Particle)*output_buffer_size);
        }
        if(system->num_chem_species > 0){
            if(output_buffer_chem_size==0){
                output_buffer_chem_size = system->particles.size() * system->num_chem_species;
                output_buffer_chem = (double*) malloc(sizeof(double)*output_buffer_chem_size);
            }else if(output_buffer_chem_size < system->particles.size() * system->num_chem_species){
                output_buffer_chem_size = system->particles.size() * system->num_chem_species;
                output_buffer_chem = (double*) realloc(output_buffer_chem, sizeof(double)*output_buffer_chem_size) ;
            }
        }
        int ncnt=0;
        Particle* p;
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p = &system->particles[i] ;
            memcpy( (void *) &output_buffer[ncnt++], (void *) p, sizeof(Particle) );
            if(system->num_chem_species > 0){
                memcpy( (void *) &output_buffer_chem[p->id*system->num_chem_species], (void*) (p->C), sizeof(double)*system->num_chem_species );
            }
        }
        output_buffer_current_num_particles = ncnt;
        if(system->num_stoch_species > 0){
            // make a copy of the RDME state vector xx
            if(output_buffer_xx_size==0){
                output_buffer_xx_size = output_buffer_current_num_particles*system->num_stoch_species;
                output_buffer_xx = (unsigned int*) malloc(sizeof(unsigned int)*output_buffer_xx_size);
            }else if(output_buffer_xx_size < output_buffer_current_num_particles*system->num_stoch_species){
                output_buffer_xx_size = output_buffer_current_num_particles*system->num_stoch_species;
                output_buffer_xx = (unsigned int*) realloc(output_buffer_xx, sizeof(unsigned int)*output_buffer_xx_size);

            }
            ncnt=0;
            for(long unsigned int i = 0; i < system->particles.size(); i++){
                p = &system->particles[i] ;
                memcpy( (void *) &output_buffer_xx[ncnt], (void *) p->xx, sizeof(unsigned int)*system->num_stoch_species );
                ncnt += system->num_stoch_species;
            }
        }
    }
    void output_vtk__async_step(ParticleSystem*system){
        FILE*fp;
        int i;
        char filename[256];
        int np = output_buffer_current_num_particles;
        static unsigned int output_index = 0;
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
        sprintf(filename,"output%u.vtk", output_index++);
        if(debug_flag){printf("Writing file '%s'\n", filename);}
        if((fp = fopen(filename,"w+"))==NULL){
            perror("Can't write output vtk file");exit(1);
        }
        fprintf(fp, "# vtk DataFile Version 4.1\n");
        fprintf(fp, "Generated by SpatialPy\n");
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
        if(system->initialized){
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
            long unsigned int s;
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
        if(system->num_stoch_species > 0){
            long unsigned int s;
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
}
