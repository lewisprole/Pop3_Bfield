int merged[NSinksAllTasks] = { 0 };  /* keep track of deleted sinks, -1 means not merged */
for(j = 0; j < NSinksAllTasks; j++)
  {
    merged[j]=-1;
  }


/*for each sink, loop through all other sinks to check if they lie within thier accretion radius, are moving towards eachother, and are bound to eachother*/
/*if all 3 criteria are met then the smaller sink will be logged to be merged with the more massive sink*/

for(j = 0; j < NSinksAllTasks; j++)
  {
    if (merged[j]>-1)  /* if it no longer exists, do nothing */
      break;
    for(i = 0; i < NSinksAllTasks; i++)
      {
        if (merged[i]>-1)
          break;  /* can’t merge with a sink that doesn’t exist anymore */

        if (SinkP[i].ID ==  SinkP[j].ID)
          break;  /* can’t merge with itself */
        
        if (merged[j]>-1) /*the current j sink could have been merged in the last i loop if mass[i]>mass[j]*/
          break;
        /* Check separation*/
        dx   = GRAVITY_NEAREST_X(SinkP[i].Pos[0] - SinkP[j].Pos[0]);
        dy   = GRAVITY_NEAREST_Y(SinkP[i].Pos[1] - SinkP[j].Pos[1]);
        dz   = GRAVITY_NEAREST_Z(SinkP[i].Pos[2] - SinkP[j].Pos[2]);
        dist = dx * dx + dy * dy + dz * dz;


#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
        SinkAccretionRadiusSquared = SinkP[i].AccretionRadius * SinkP[i].AccretionRadius;
#endif

        if(dist > SinkAccretionRadiusSquare)  /*not close enough to merge*/
           break;

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

        if((vrad > 0) || (arad > 0)) /*not moving towards each other*/
          break;

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
        if(e_total>0) /*not bound*/
          break;
        /*passed all merger tests*/
        /*going to merge massive sink with smaller sink*/
        if (sink_mass_i=>sink_mass_j)
          {
            int keep=i;
            int destroy=j;
            for(z = 0; z < NSinksAllTasks; z++) /*i inherits all of the existing mergers to j  (otherwise they'll merge with nothing*/
              {
                if merged[z]==j
                  merged[z]==i;
              }
          }
        else
          {
            int keep=j;
            int destroy=i;
          }

        merged[destroy]=keep /*merged array is -1 for surviving sink, or the argument of the sink it will be eaten by*/

      } /*end of [i] loop*/
  }   
  


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
                   int index_merge=SinkP[j].Index;
                   int index_survive=SinkP[i].Index;
                   if(P[index_merge].ID != SinkP[j].ID || P[index_survive].ID != SinkP[i].ID)
                     terminate("Something is wrong, the sink is indexing the wrong particle, check maybe domain_rearrange.c where we eliminate sinks if we messed up the indexing there\n");

                   /*move to centre of mass */
                   P[index_survive].Pos[0]=(P[index_survive].Mass*P[index_survive].Pos[0]+P[index_merge].Mass*P[index_merge].Pos[0])/(P[index_survive].Mass+P[index_merge].Mass);
                   P[index_survive].Pos[1]=(P[index_survive].Mass*P[index_survive].Pos[1]+P[index_merge].Mass*P[index_merge].Pos[1])/(P[index_survive].Mass+P[index_merge].Mass);
                   P[index_survive].Pos[2]=(P[index_survive].Mass*P[index_survive].Pos[2]+P[index_merge].Mass*P[index_merge].Pos[2])/(P[index_survive].Mass+P[index_merge].Mass);

                   /*transfer linear momentum*/
                   P[index_survive].Vel[0]=(P[index_survive].Mass*P[index_survive].Vel[0]+P[index_merge].Mass*P[index_merge].Vel[0])/(P[index_survive].Mass+P[index_merge].Mass);
                   P[index_survive].Vel[1]=(P[index_survive].Mass*P[index_survive].Vel[1]+P[index_merge].Mass*P[index_merge].Vel[1])/(P[index_survive].Mass+P[index_merge].Mass);
                   P[index_survive].Vel[0]=(P[index_survive].Mass*P[index_survive].Vel[2]+P[index_merge].Mass*P[index_merge].Vel[2])/(P[index_survive].Mass+P[index_merge].Mass);

                   /*transfer mass*/
                   [index_survive].Mass+=P[index_merge].Mass;
                 }
             }
         }       
     }



/*remove the eaten sinks from the type 5 array and set mass,vel,accel to 0 */     
for (int i = 0; i < NSinksAllTasks; i++)
  {
    if (merged[i]>-1)
      {
        if (SinkP[i].HomeTask == ThisTask)
          {
            int index=SinkP[i].Index;
            P[index].Type==3
            P[index].Mass=0
            P[index].Vel[0]=0
            P[index].Vel[1]=0
            P[index].Vel[2]=0
            P[index].Accel[0]=0
            P[index].Accel[1]=0
            P[index].Accel[2]=0
          }
      }  
  }  


