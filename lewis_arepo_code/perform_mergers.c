#include "../allvars.h"
#include "../proto.h"

void perform_mergers(void)
  {
  if(NSinksAllTasks>1)
  {
    mpi_printf("SINK_MERGERS initialising \n");
    int merged[NSinksAllTasks];  /* keep track of deleted sinks, -1 means not merged */
    int j;
    int i;
    int z;
    double dx;
    double dy;
    double dz;
    double dist;
    double rad_sink_closest;
    double dvx;
    double dvy;
    double dvz;
    double vrad;
    double dax;
    double day;
    double daz;
    double arad;
    double dv;
    double sink_mass_i;
    double sink_mass_j;
    double ekin;
    double egrav;
    double SinkAccretionRadiusSquared = SinkAccretionRadius * SinkAccretionRadius;
    int index_merge;
    int index_survive;

    for(j = 0; j < NSinksAllTasks; j++)
      {
        merged[j]=-1;
      }

    /*for each sink, loop through all other sinks to check if they lie within thier accretion radius, are moving towards eachother, and are bound to eachother*/
    /*if all 3 criteria are met then the smaller sink will be logged to be merged with the more massive sink*/
    mpi_printf("SINK_MERGERS searching for candidates \n");
    for(j = 0; j < NSinksAllTasks; j++)
      {
        mpi_printf("SINK_MERGERS j = %d \n",j);
        if (merged[j]==-1)  /* if it no longer exists, do nothing */
          {
            for(i = 0; i < NSinksAllTasks; i++)
              {
                mpi_printf("SINK_MERGERS i = %d \n",i);
                if (merged[i]==-1)  /* can’t merge with a sink that doesn’t exist anymore */
                  {
                    if (SinkP[i].ID !=  SinkP[j].ID)  /* can’t merge with itself */
                      { 
                        if (merged[j]==-1) /*the current j sink could have been merged in the last i loop if mass[i]>mass[j]*/
                          {
                            /* Check separation*/
                            dx   = GRAVITY_NEAREST_X(SinkP[i].Pos[0] - SinkP[j].Pos[0]);
                            dy   = GRAVITY_NEAREST_Y(SinkP[i].Pos[1] - SinkP[j].Pos[1]);
                            dz   = GRAVITY_NEAREST_Z(SinkP[i].Pos[2] - SinkP[j].Pos[2]);
                            dist = dx * dx + dy * dy + dz * dz;
#ifdef SINK_MERGERS_DEBUG
                            mpi_printf("SINK_MERGERS distance from sink %d to %d = %g, r_accrete = %g \n",j,i,dist,SinkAccretionRadiusSquared);
#endif
                            if((i == 1) && (j==0))
                              rad_sink_closest = dist;
                            else
                              {
                                if(dist < rad_sink_closest)
                                   rad_sink_closest = dist;
                              }

#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
                              SinkAccretionRadiusSquared = SinkP[i].AccretionRadius * SinkP[i].AccretionRadius;
#endif

                              if(dist > SinkAccretionRadiusSquared)  /*close enough to merge*/
                               
#ifdef SINK_MERGERS_DEBUG
                                mpi_printf("SINK_MERGERS not close enough to merge \n");
#else
                                ;
#endif                                
                              else
                                {
                                  
                                  mpi_printf("SINK_MERGERS  pos (%g,%g,%g) and (%g,%g,%g) \n",SinkP[j].Pos[0],SinkP[j].Pos[1],SinkP[j].Pos[2],SinkP[i].Pos[0],SinkP[i].Pos[1],SinkP[i].Pos[2]);
                                  /* Check divergence*/
                                  dist = sqrt(dist);

                                  dvx  = SinkP[i].Vel[0] - SinkP[j].Vel[0];
                                  dvy  = SinkP[i].Vel[1] - SinkP[j].Vel[1];
                                  dvz  = SinkP[i].Vel[2] - SinkP[j].Vel[2];
                                  vrad = (dvx * dx + dvy * dy + dvz * dz) / dist;

                                  dax  = SinkP[i].Accel[0] - SinkP[j].Accel[0];
                                  day  = SinkP[i].Accel[1] - SinkP[j].Accel[1];
                                  daz  = SinkP[i].Accel[2] - SinkP[j].Accel[2];
                                  arad = (dax * dx + day * dy + daz * dz) / dist;

                                  if((vrad > 0) || (arad > 0)) /*moving towards each other*/
#ifdef SINK_MERGERS_DEBUG
                                    mpi_printf("SINK_MERGERS not oving towards eachother \n");
#else
                                    ;
#endif
                                  else
                                    {
                                      

                                      /* Check if bound*/
                                      dv   = dvx * dvx + dvy * dvy + dvz * dvz;
                                      if(SinkP[i].FormationTime == All.Time)   /*I think newly formed sinks don’t have mass yet?*/
                                        sink_mass_i = SinkP[i].FormationMass;
                                      else
                                        sink_mass_i = SinkP[i].Mass;

                                      if(SinkP[j].FormationTime == All.Time)  
                                        sink_mass_j = SinkP[j].FormationMass;
                                      else
                                        sink_mass_j = SinkP[j].Mass;

                                      ekin = 0.5 * sink_mass_i * dv; /*energies*/
                                      egrav  = All.G * sink_mass_i * sink_mass_j / dist;          
                                      if(All.ComovingIntegrationOn)            /*Some cosmo stuff idk*/
                                        {
                                          /*converting energies to physical units*/
                                          ekin /= (All.Time * All.Time);
                                          egrav /= All.Time;
                                        }

                                      int e_total= ekin- egrav; 
                                      if(e_total>0) /*bound*/
#ifdef SINK_MERGERS_DEBUG
                                        mpi_printf("SINK_MERGERS not bound \n");
#else
                                        ;
#endif
                                      else
                                        {
#ifdef SINK_MERGERS_DEBUG
                                          mpi_printf("SINK_MERGERS a sink has passed the test \n");
                                          mpi_printf("SINK_MERGERS closest sinks at d = %g while r_accrete = %g \n",rad_sink_closest,SinkAccretionRadiusSquared);	
#endif
                                          /*passed all merger tests*/
                                          /*going to merge massive sink with smaller sink*/
                                          if (sink_mass_i>=sink_mass_j)
                                            {
#ifdef SINK_MERGERS_DEBUG
                                              mpi_printf("SINK_MERGERS sink larger than current sink - destroying current sink \n");
#endif
                                              int keep=i;
                                              int destroy=j;
                                              merged[destroy]=keep; /*merged array is -1 for surviving sink, or the argument of the sink it will be eaten by*/
                                              for(z = 0; z < NSinksAllTasks; z++) /*i inherits all of the existing mergers to j  (otherwise they'll merge with nothing*/
                                               {
                                                 if (merged[z]==j);
#ifdef SINK_MERGERS_DEBUG
                                                 mpi_printf("SINK_MERGERS sink %d mergers transferred to sink  %d \n",j,i);
#endif
                                                 merged[z]=i;
                                               }
                                            }
                                          else
                                            {
#ifdef SINK_MERGERS_DEBUG
                                               mpi_printf("SINK_MERGERS sink smaller than current sink - keeping current sink \n");
#endif
                                               int keep=j;
                                               int destroy=i;
                                               merged[destroy]=keep; /*merged array is -1 for surviving sink, or the argument of the sink it will be eaten by*/
#ifdef SINK_MERGERS_DEBUG
                                           mpi_printf("SINK_MERGERS merger array ");

                                           for(z = 0; z < NSinksAllTasks; z++)
                                             {
                                               mpi_printf(" %d ",merged[z]);
                                             }
                                           mpi_printf("\n");      
#endif    
                                            }
                                        }
                                    }
                                }                                                         
                            
                          }
                      }
                  }
              }
          }
      } 
          
      
  
     int Nmerge = 0;
     for (int i = 0; i < NSinksAllTasks; i++)
       {
         if(merged[i]>-1)
           Nmerge+=1;
       }
     mpi_printf("SINK_MERGERS performing %d mergers \n",Nmerge);



     if (Nmerge>0)
       {
         /*Perform merge by transferring mass and linear momentum, moving to center of mass*/
         for (int i = 0; i < NSinksAllTasks; i++)          
           {
             if (merged[i]==-1) /*surviving sink*/
               {
                 if (SinkP[i].HomeTask == ThisTask) /*on this task*/
                   {
                     for (int j = 0; j < NSinksAllTasks; j++) /*find which sinks to eat*/
                       {  
                         if (merged[j]==i)
                           {
#ifdef SINK_MERGERS_DEBUG
                             mpi_printf("SINK_MERGERS transferring properties on task %d \n",ThisTask);
#endif
                             index_merge=SinkP[j].Index;
                             index_survive=SinkP[i].Index;
                             int iloc = SinkP[j].Index;  
#ifdef SINK_MERGERS_DEBUG
                             mpi_printf("SINK_MERGERS does index %d match %d \n",P[iloc].ID,SinkP[j].ID);
                             mpi_printf("SINK_MERGERS masses %g and %g \n",P[index_survive].Mass,P[index_merge].Mass);
                             mpi_printf("SINK_MERGERS old pos (%g,%g,%g) and (%g,%g,%g) \n",P[index_survive].Pos[0],P[index_survive].Pos[1],P[index_survive].Pos[2],P[index_merge].Pos[0],P[index_merge].Pos[1],P[index_merge].Pos[2]);
                             mpi_printf("SINK_MERGERS old vels (%g,%g,%g) and (%g,%g,%g) \n",P[index_survive].Vel[0],P[index_survive].Vel[1],P[index_survive].Vel[2],P[index_merge].Vel[0],P[index_merge].Vel[1],P[index_merge].Vel[2]);
#endif
                             /*move to centre of mass */
                             P[index_survive].Pos[0]=(P[index_survive].Mass*P[index_survive].Pos[0]+P[index_merge].Mass*P[index_merge].Pos[0])/(P[index_survive].Mass+P[index_merge].Mass);
                             P[index_survive].Pos[1]=(P[index_survive].Mass*P[index_survive].Pos[1]+P[index_merge].Mass*P[index_merge].Pos[1])/(P[index_survive].Mass+P[index_merge].Mass);
                             P[index_survive].Pos[2]=(P[index_survive].Mass*P[index_survive].Pos[2]+P[index_merge].Mass*P[index_merge].Pos[2])/(P[index_survive].Mass+P[index_merge].Mass);

                           
                           /*transfer linear momentum*/
                             P[index_survive].Vel[0]=(P[index_survive].Mass*P[index_survive].Vel[0]+P[index_merge].Mass*P[index_merge].Vel[0])/(P[index_survive].Mass+P[index_merge].Mass);
                             P[index_survive].Vel[1]=(P[index_survive].Mass*P[index_survive].Vel[1]+P[index_merge].Mass*P[index_merge].Vel[1])/(P[index_survive].Mass+P[index_merge].Mass);
                             P[index_survive].Vel[0]=(P[index_survive].Mass*P[index_survive].Vel[2]+P[index_merge].Mass*P[index_merge].Vel[2])/(P[index_survive].Mass+P[index_merge].Mass);
                             /*transfer mass*/
                             P[index_survive].Mass+=P[index_merge].Mass;

#ifdef SINK_MERGERS_DEBUG
                             mpi_printf("SINK_MERGERS sink %d has eaten sink %d \n" ,i,j);
                             mpi_printf("SINK_MERGERS new pos=(%g,%g,%g), vel=(%g,%g,%g), mass=%g \n",P[index_survive].Pos[0],P[index_survive].Pos[1],P[index_survive].Pos[2],P[index_survive].Vel[0],P[index_survive].Vel[1],P[index_survive].Vel[2],P[index_survive].Mass);
#endif
                           }
                       }
                   }
               }
           }       
       }


#ifdef SINK_MERGERS_DEBUG
     mpi_printf("SINK_MERGERS merger array ");

     for(z = 0; z < NSinksAllTasks; z++)
       {
         mpi_printf(" %d ",merged[z]);
       }
     mpi_printf("\n");
#endif


    if (Nmerge>0)
      {
        /*remove the eaten sinks from the type 5 array and set mass,vel to 0 */     
        for (int i = 0; i < NSinksAllTasks; i++)
          {
            if (merged[i]>-1)
              {
                mpi_printf("found a sink to delte, %d \n",i);
                if (SinkP[i].HomeTask == ThisTask)
                  {
                    mpi_printf("on task %d \n",ThisTask);
#ifdef SINK_MERGERS_DEBUG
                    mpi_printf("SINK_MERGERS deleting sink %d, making into type 3, on task %d \n",i,ThisTask);
#endif
                    mpi_printf("deleting");
                    index_merge=SinkP[i].Index;
                    P[index_merge].Type=3;
                    P[index_merge].Mass=0;
                    P[index_merge].Vel[0]=0;
                    P[index_merge].Vel[1]=0;
                    P[index_merge].Vel[2]=0;
                  }
              }  
          }
       }
    }  
    return;
  }

